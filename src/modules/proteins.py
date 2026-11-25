"""Protein haplotypes of a transcript (IMPLEMENTATION_NOTES section 3.5): a port of the reference's
process_transcript_haplotypes (legacy/src/modules/process_haplotypes.py), step for step, with these changes only:

- sequences are plain str instead of Bio.Seq (slicing, concatenation and comparison behave the same);
- per change string, the parse, get_rna_position, the reverse complement and the liftover check are computed once
  per transcript; the reference's affected codons are cached;
- when the changes come in increasing, non-overlapping reference order and every REF allele matches the reference
  cDNA, the haplotype's cDNA is built with one join (every check of the reference's loop then reads untouched
  reference sequence); otherwise the reference's loop runs as it is;
- translate() uses a codon table for sequences of A, C, G and T only, and Biopython for anything else;
- sequences are deduplicated in a dict by hash() and sorted by hash once, the order of the reference's sorted list;
- a merged row (the collapsing fix) collects its haplotypes and adds up counts and frequencies in the reference's
  order; population and sample strings are built only for rows that pass the threshold.

With prohap.py -fix_bugs 1 (Settings.fixes, IMPLEMENTATION_NOTES section 9) the fixes of the reference's bugs apply
on top; with none, the output is the reference's.

The reference's merge keys are kept: str(hash(IDs)) with a known reading frame (merge), hash(IDs) without (a later
haplotype replaces an earlier one, IMPLEMENTATION_NOTES section 7.4). Haplotype IDs are the rows' local indices; the
caller turns them into `<prefix>_<hex>`.
"""
import bisect
import re

import numpy as np
from Bio.Data import CodonTable
from Bio.Seq import Seq

from modules.fallback import Unsupported
from modules.haplotypes import conflicting
from modules.legacy_import import legacy

_coordinates = legacy.coordinates_toolbox
_check_protein_allele = legacy.process_haplotypes.check_protein_allele


def add_population_freqs(strings):
    """The reference's add_population_freqs applied left to right: add_population_freqs(add_population_freqs(s0,
    s1), s2)... Each step adds the frequencies of codes in either string and writes the sums with five decimals, which
    the next step reads back; codes keep the order of first appearance. Codes contain no ':' or ';' (Cohort checks)."""
    total = {}
    for elem in strings[0].split(';'):
        code, value = elem.split(':')
        total[code] = float(value)
    for string in strings[1:]:
        right = {}
        for elem in string.split(';'):
            code, value = elem.split(':')
            right[code] = float(value)
        for code in right:
            if code not in total:
                total[code] = 0
        total = {code: float('{:.5f}'.format(0 + total[code] + right[code]) if code in right
                             else '{:.5f}'.format(0 + total[code])) for code in total}
    if len(strings) == 1:
        return strings[0]
    return ';'.join([code + ':{:.5f}'.format(value) for code, value in total.items()])

# ---------------------------------------------------------------------------------------------------------------
# translation

_TABLE = CodonTable.unambiguous_dna_by_id[1]
_CODONS = {a + b + c: (_TABLE.forward_table.get(a + b + c) or '*') for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'}
assert sorted(c for c, aa in _CODONS.items() if aa == '*') == sorted(_TABLE.stop_codons)
_BASE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate('ACGT'):
    _BASE[ord(_b)] = _i
_LUT = np.frombuffer(''.join(_CODONS[a + b + c] for a in 'ACGT' for b in 'ACGT' for c in 'ACGT').encode('ascii'),
                     dtype=np.uint8)
_RC = str.maketrans('ACGT', 'TGCA')


_SHORT = {}                     # translations of short sequences (the codons around a change), by sequence


def translate(seq):
    """str(Seq(seq).transcribe().translate()): the standard code, '*' for stops, a trailing partial codon ignored."""
    n = len(seq) // 3
    if n <= 30:
        protein = _SHORT.get(seq)
        if protein is None:
            try:
                protein = ''.join([_CODONS[seq[i:i + 3]] for i in range(0, 3 * n, 3)]) \
                    if seq[3 * n:].strip('ACGT') == '' else _biopython(seq)
            except KeyError:
                protein = _biopython(seq)
            if len(_SHORT) >= 1 << 18:
                _SHORT.clear()
            _SHORT[seq] = protein
        return protein
    try:
        b = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    except UnicodeEncodeError:
        return _biopython(seq)
    codes = _BASE[b]
    if (codes == 255).any():
        return _biopython(seq)
    codes = codes[:3 * n].astype(np.intp)
    return _LUT[codes[0::3] * 16 + codes[1::3] * 4 + codes[2::3]].tobytes().decode('ascii')


def _biopython(seq):
    return str(Seq(seq).transcribe().translate())


def reverse_complement(seq):
    """str(Seq(seq).reverse_complement())"""
    if seq.strip('ACGT') == '':
        return seq.translate(_RC)[::-1]
    return str(Seq(seq).reverse_complement())


def get_affected_codons(cdna, allele_loc, allele_len, reading_frame, protein_start):
    """legacy coordinates_toolbox.get_affected_codons on str, with integer arithmetic (floor((a - b) / 3) is
    (a - b) // 3 for the integers involved)."""
    if reading_frame == -1:
        protein_location = [(allele_loc - rf) // 3 for rf in range(3)]
    else:
        protein_location = [(allele_loc - reading_frame) // 3 - protein_start]
    rf0 = max(reading_frame, 0)
    bp_from = (allele_loc - rf0) // 3 * 3 + rf0
    bp_from = max(max(bp_from, 0), reading_frame)
    bp_to = -((rf0 - allele_loc - allele_len) // 3) * 3 + rf0
    alleles_protein = [translate(cdna[bp_from:bp_to])] if bp_to - bp_from > 2 else ['-']
    if reading_frame == -1:
        for rf in (1, 2):
            bp_from = (allele_loc - rf) // 3 * 3 + rf
            bp_from = max(max(bp_from, 0), rf)
            bp_to = -((rf - allele_loc - allele_len) // 3) * 3 + rf
            alleles_protein.append(translate(cdna[bp_from:bp_to]) if bp_to - bp_from > 2 else '-')
    return alleles_protein, protein_location


# ---------------------------------------------------------------------------------------------------------------
# one transcript

class Transcript:
    """What process_transcript_haplotypes reads from the annotation and the cDNA FASTA for one transcript. exons:
    objects with start and end, ordered by start; start_codon / stop_codon: the start position of the first one, or
    None; cdna: the reference cDNA, or None when the transcript is not in the cDNA FASTA."""
    __slots__ = ('id', 'start', 'end', 'strand', 'exons', 'start_codon', 'stop_codon', 'biotype', 'cdna',
                 'exon_starts', 'exon_max_ends')

    def __init__(self, id, start, end, strand, exons, start_codon, stop_codon, biotype, cdna):
        self.id = id
        self.start = start
        self.end = end
        self.strand = strand
        self.exons = exons
        self.start_codon = start_codon
        self.stop_codon = stop_codon
        self.biotype = biotype
        self.cdna = cdna


class Settings:
    """prohap.py's options. fixes: the names of the fixes to apply (modules/fixes.py); empty, the reference's
    behaviour."""
    def __init__(self, chromosome, force_rf, min_freq, min_count, ignore_utr, skip_start_loss, output_cdnas,
                 fixes=frozenset()):
        self.chromosome = chromosome
        self.force_rf = force_rf
        self.min_freq = min_freq
        self.min_count = min_count
        self.ignore_utr = ignore_utr
        self.skip_start_loss = skip_start_loss
        self.output_cdnas = output_cdnas
        self.fixes = frozenset(fixes)


class _Merged:
    """A row of the result: the reference's 21 columns, with count and frequency added up as haplotypes merge, and
    the haplotypes behind its population and sample strings."""
    __slots__ = ('fields', 'parts')

    def __init__(self, fields, part):
        self.fields = fields
        self.parts = [part]


def check_protein_allele_fixed(change, start, stop):
    """legacy check_protein_allele with the fix 'stop' (IMPLEMENTATION_NOTES section 9.4): whether a protein change
    `loc:REF>loc:ALT` is in the coding sequence, which runs from the start (residue `start`) to the first stop codon
    (residue `stop`, the stop included). The reference requires every ALT residue to lie before `stop`, so the change
    that creates the first stop, whose ALT ends with it, was left out as a UTR change; the length it measures also
    counts the '(+fs)' or '(fs)' after a frameshift. Here a change is in when it starts at or after the start and its
    first ALT residue that differs from REF lies at or before the stop: a change that leaves the residues up to the
    stop as they are (an insertion right after the stop codon, whose window starts at the stop) is still 3' UTR, and
    one whose ALT and REF agree throughout (synonymous) is in exactly when the reference has it in."""
    fields = change.split(':')
    ref = fields[1].split('>')[0]
    loc = int(fields[1].split('>')[1]) + start
    alt = fields[2]
    for mark in ('(+fs)', '(fs)'):
        if alt.endswith(mark):
            alt = alt[:-len(mark)]
            break
    same = 0
    while same < min(len(ref), len(alt)) and ref[same] == alt[same]:
        same += 1
    return start <= loc and loc + same <= stop


def _trim(rna_location, ref_allele, alt_allele, reverse_strand, cdna_length):
    """A change as get_rna_position gives it (genome orientation, after truncation at splice sites) without the bases
    REF and ALT share: the shared prefix, then the shared suffix, trimmed in genome orientation as for the conflict
    check, then turned to the cDNA's orientation. Returns (cDNA location, REF, ALT) of the bases really changed; REF
    is empty for an insertion before the base at that location (fix 'conflicts', IMPLEMENTATION_NOTES section 9.2)."""
    n = min(len(ref_allele), len(alt_allele))
    k = 0
    while k < n and ref_allele[k] == alt_allele[k]:
        k += 1
    m = 0
    while m < n - k and ref_allele[len(ref_allele) - 1 - m] == alt_allele[len(alt_allele) - 1 - m]:
        m += 1
    ref_t = ref_allele[k:len(ref_allele) - m]
    alt_t = alt_allele[k:len(alt_allele) - m]
    loc = rna_location + k
    if reverse_strand:
        return cdna_length - loc - len(ref_t), reverse_complement(ref_t), reverse_complement(alt_t)
    return loc, ref_t, alt_t


def _apply_trimmed(cdna, anchored, trimmed):
    """The haplotype's cDNA from trimmed changes (fix 'conflicts'), or None with the reason when the changes cannot be
    placed: an anchored REF (what the VCF says is there, anchor base included) that is not in the reference cDNA, or
    trimmed changes that conflict in the cDNA. anchored: (location, REF) per change; trimmed: (location, REF, ALT)."""
    for loc, ref_allele in anchored:
        if loc < 0 or cdna[loc:loc + len(ref_allele)] != ref_allele:
            return None, 'Ref allele not matching the cDNA sequence, skipping haplotype!'
    changes = [c for c in trimmed if c[1] or c[2]]
    if conflicting([(loc, loc + len(ref_t)) for loc, ref_t, _ in changes]):
        return None, 'Changes overlap in the cDNA, skipping haplotype!'
    pieces = []
    previous_end = 0
    for loc, ref_t, alt_t in sorted(changes, key=lambda c: (c[0], c[0] + len(c[1]))):
        pieces.append(cdna[previous_end:loc])
        pieces.append(alt_t)
        previous_end = loc + len(ref_t)
    pieces.append(cdna[previous_end:])
    return ''.join(pieces), None


def _offsets(rna_locations, trimmed):
    """For each change, how much the changes before its (anchored) location lengthen the cDNA: the sum of
    len(ALT) - len(REF) of the other trimmed changes that end at or before it (an insertion at a location comes
    before the base there). Where the change's anchored alleles start in the haplotype's cDNA, for the annotation."""
    ends = sorted((loc + len(ref_t), len(alt_t) - len(ref_t)) for loc, ref_t, alt_t in trimmed)
    keys = [e for e, _ in ends]
    prefix = [0]
    for _, diff in ends:
        prefix.append(prefix[-1] + diff)
    result = []
    for x, (loc, ref_t, alt_t) in zip(rna_locations, trimmed):
        total = prefix[bisect.bisect_right(keys, x)]
        if loc + len(ref_t) <= x:
            total -= len(alt_t) - len(ref_t)
        result.append(total)
    return result


def start_change(start, changes, force_rf):
    """Fix 'start' (IMPLEMENTATION_NOTES section 9.5): what the changes do to the start codon at `start` of the
    reference cDNA. changes: trimmed changes (location, REF, ALT) in reference cDNA coordinates. Returns (lost, the
    indexes of the changes that touch the start codon, the start's location in the haplotype's cDNA). The codon is
    lost when a change overlaps it or inserts bases strictly inside it, or (without force_rf) a frameshift lies before
    it; changes before it move it. The reference checks each anchored change against a start it has already moved,
    mixing the coordinates of the reference and the haplotype, and counts an indel whose anchor base is in the codon."""
    touching = []
    shift = 0
    for k, (loc, ref_t, alt_t) in enumerate(changes):
        if not ref_t and not alt_t:
            continue
        end = loc + len(ref_t)
        if (ref_t and loc < start + 3 and start < end) or (not ref_t and start < loc < start + 3):
            touching.append(k)
        elif end <= start:
            shift += len(alt_t) - len(ref_t)
    lost = bool(touching) or (shift % 3 != 0 and not force_rf)
    return lost, touching, start + shift


def rescan_start(mutated_cdna, start, changes):
    """Fix 'start': where translation starts on a haplotype whose start codon is lost, in the scanning model: the
    first ATG of the haplotype's cDNA, in any frame, from where the lost start codon began (the start of the first
    change that touches it, if that comes earlier). Returns its location in the haplotype's cDNA, or -1."""
    _, touching, _ = start_change(start, changes, True)
    origin = min([start] + [changes[k][0] for k in touching])
    shift = sum(len(alt_t) - len(ref_t) for loc, ref_t, alt_t in changes
                if (ref_t or alt_t) and loc + len(ref_t) <= origin)
    return mutated_cdna.find('ATG', origin + shift)


def _reading_frame(t, cdna_length, force_rf):
    """(reading_frame, start_loc, protein_start) before any change, as the reference computes them per haplotype."""
    reverse_strand = t.strand == '-'
    try:
        if t.start_codon is not None:
            start_loc = _coordinates.get_rna_position_simple(t.id, t.start_codon, t.exons)
            if reverse_strand:
                start_loc = cdna_length - start_loc - 3
            reading_frame = start_loc % 3
            return reading_frame, start_loc, int((start_loc - reading_frame) / 3)
        if (t.stop_codon is not None) and force_rf:
            stop_loc = _coordinates.get_rna_position_simple(t.id, t.stop_codon, t.exons)
            if reverse_strand:
                stop_loc = cdna_length - stop_loc - 3
            return stop_loc % 3, 0, 0
    except Exception as e:      # the reference stops here: a start or stop codon outside the exons
        raise Unsupported('%s: %s' % (t.id, e))
    return -1, 0, 0


def process_transcript(t, rows, settings, cohort):
    """process_transcript_haplotypes for one transcript. rows: (local index, haplotypes.Row) of the transcript's table
    in order, the reference haplotype included. Returns (result rows, protein sequences, cDNA sequences): rows are
    the reference's 21 columns with HaplotypeID the local index; sequences are dicts with 'haplotypes' (local
    indices), 'sequence', 'start' and, for proteins, 'rfs'."""
    transcript_id = t.id
    if t.cdna is None:
        print('Transcript', transcript_id, 'not in cDNA database, skipping!')
        return [], [], []

    exons = t.exons
    cdna_sequence = t.cdna
    reverse_strand = t.strand == '-'
    force_rf = settings.force_rf
    fix_conflicts = 'conflicts' in settings.fixes
    fix_liftover = 'liftover' in settings.fixes
    check_protein_allele = check_protein_allele_fixed if 'stop' in settings.fixes else _check_protein_allele
    fix_start = 'start' in settings.fixes
    fix_frameshift = 'frameshift' in settings.fixes
    trimmed_cache = {}      # change string -> _trim (fixes 'conflicts' and 'liftover')
    skip_start_loss = settings.skip_start_loss
    output_cdnas = settings.output_cdnas
    initial_frame = None

    change_cache = {}       # change string -> values from the reference's first loop
    in_cds_cache = {}       # (protein change, start, stop) -> check_protein_allele
    ref_codons = {}         # (loc, length, reading frame, protein start) -> get_affected_codons on the reference

    local_result_data = {}
    protein_sequences = {}  # hash -> entry
    cdna_sequences = {}

    for index, row in rows:
        changes = row.changes
        if changes == 'REF':
            continue
        if initial_frame is None:
            initial_frame = _reading_frame(t, len(cdna_sequence), force_rf)
        reading_frame, start_loc, protein_start = initial_frame
        reading_frame_ref = reading_frame
        protein_start_ref = protein_start
        sequence_length_diff = 0

        cDNA_changes = []
        all_protein_changes = []
        protein_changes = []
        spl_junctions_affected = []
        frameshifts = []
        dna_var_types = []
        prot_var_types = []

        all_changes = changes.split(';')
        all_AFs = row.afs.split(';')
        all_vcf_IDs = row.ids.split(';')
        if reverse_strand:
            all_AFs.reverse()
            all_changes.reverse()
            all_vcf_IDs.reverse()

        ref_alleles = []
        alt_alleles = []
        rna_locations = []
        start_lost = False
        validity_check = True
        tmp_changes = []
        tmp_AFs = []
        tmp_VCF_IDs = []

        for ch_idx, change in enumerate(all_changes):
            cached = change_cache.get(change)
            if cached is None:
                ref_allele = re.split(r'\d+', change)[1].split('>')[0][1:]
                alt_allele = re.split(r'\d+', change)[1].split('>')[1]
                if ref_allele == '-':
                    ref_allele = ''
                if alt_allele == '-':
                    alt_allele = ''
                dna_location = int(re.split(r'[\-a-zA-Z\*]+', change)[0][:-1])
                rna_location, ref_allele, ref_len, alt_allele, alt_len, mutation_intersects_intron = \
                    _coordinates.get_rna_position(transcript_id, dna_location, ref_allele, alt_allele, exons)
                if fix_conflicts or fix_liftover or fix_start:
                    trimmed_cache[change] = _trim(rna_location, ref_allele, alt_allele, reverse_strand,
                                                  len(cdna_sequence))
                if reverse_strand:
                    ref_allele = reverse_complement(ref_allele)
                    alt_allele = reverse_complement(alt_allele)
                    rna_location = len(cdna_sequence) - rna_location - ref_len
                if fix_liftover:
                    # fix 'liftover' (IMPLEMENTATION_NOTES section 9.3): left out only when the reference cDNA has
                    # the ALT allele and not the REF allele, both without the bases they share (an indel's anchor)
                    loc_t, ref_t, alt_t = trimmed_cache[change]
                    keep = not (cdna_sequence[loc_t:loc_t + len(alt_t)] == alt_t and
                                cdna_sequence[loc_t:loc_t + len(ref_t)] != ref_t)
                else:
                    keep = alt_allele != cdna_sequence[rna_location:rna_location + alt_len]
                cached = change_cache[change] = (keep, rna_location, ref_allele, ref_len, alt_allele, alt_len,
                                                 mutation_intersects_intron)
            keep, rna_location, ref_allele, ref_len, alt_allele, alt_len, mutation_intersects_intron = cached
            # check if the alternative allele isn't the reference (error in liftover from older build)
            if not keep:
                continue
            tmp_changes.append(all_changes[ch_idx])
            tmp_AFs.append(all_AFs[ch_idx])
            tmp_VCF_IDs.append(all_vcf_IDs[ch_idx])

            if (mutation_intersects_intron is not None) and (mutation_intersects_intron not in spl_junctions_affected):
                spl_junctions_affected.append(int(mutation_intersects_intron))
            frameshifts.append((abs(ref_len - alt_len) % 3) != 0)
            if (mutation_intersects_intron is not None):
                dna_var_types.append('splice')
            elif (ref_len == alt_len):
                dna_var_types.append('SNP')
            else:
                dna_var_types.append('indel')

            if (t.start_codon is not None) and (reading_frame > -1) and not fix_start:
                start_loc, reading_frame = _coordinates.check_start_change(start_loc, reading_frame, rna_location,
                                                                           ref_len, alt_len, force_rf)
                if (start_loc == -1):
                    if (skip_start_loss):
                        validity_check = False
                        break
                    dna_var_types[-1] = 'start_lost'
                    start_loc = 0
                    protein_start = 0
                    protein_start_ref = 0
                    reading_frame_ref = -1
                    start_lost = True
                else:
                    protein_start = int((start_loc - reading_frame) / 3)

            rna_locations.append(rna_location)
            ref_alleles.append(ref_allele)
            alt_alleles.append(alt_allele)

        all_changes = tmp_changes
        all_AFs = tmp_AFs
        all_vcf_IDs = tmp_VCF_IDs
        if ((not validity_check) or (len(all_changes) == 0)):
            continue

        rescan = False          # fix 'start' with skip_start_lost: find the start again once the cDNA is built
        if fix_start and (t.start_codon is not None) and (reading_frame > -1):
            changes_t = [trimmed_cache[change] for change in all_changes]
            lost, touching, moved = start_change(start_loc, changes_t, force_rf)
            if not lost:
                start_loc = moved
                reading_frame = start_loc % 3
                protein_start = int((start_loc - reading_frame) / 3)
            else:
                for k in touching:
                    dna_var_types[k] = 'start_lost'
                start_lost = True
                if skip_start_loss:
                    rescan = True
                else:
                    # as the reference without skip_start_lost: translated in all three frames
                    start_loc = 0
                    protein_start = 0
                    protein_start_ref = 0
                    reading_frame_ref = -1
                    reading_frame = -1

        # the haplotype's cDNA
        mutated_cdna = None
        previous_end = 0
        ordered = True
        trimmed = None          # fix 'conflicts': the trimmed changes, when the anchored ones overlap
        for k in range(len(all_changes)):
            loc = rna_locations[k]
            ref_allele = ref_alleles[k]
            if loc < previous_end or loc + len(ref_allele) > len(cdna_sequence) or \
                    cdna_sequence[loc:loc + len(ref_allele)] != ref_allele:
                ordered = False
                break
            previous_end = loc + len(ref_allele)
        if ordered:
            pieces = []
            previous_end = 0
            for k in range(len(all_changes)):
                loc = rna_locations[k]
                pieces.append(cdna_sequence[previous_end:loc])
                pieces.append(alt_alleles[k])
                previous_end = loc + len(ref_alleles[k])
                cDNA_changes.append(str(loc) + ':' + ref_alleles[k] + '>' + alt_alleles[k])
            pieces.append(cdna_sequence[previous_end:])
            mutated_cdna = ''.join(pieces)
        elif fix_conflicts:
            # the anchored forms overlap where only an anchor base is shared (cases 1-5 of section 9.2): build the
            # cDNA from the bases each change really changes
            trimmed = [trimmed_cache[change] for change in all_changes]
            mutated_cdna, problem = _apply_trimmed(cdna_sequence, list(zip(rna_locations, ref_alleles)), trimmed)
            if mutated_cdna is None:
                print(problem)
                print(transcript_id + ' strand ' + t.strand + ' ' + ';'.join(all_changes))
                continue
            for ch_idx in range(len(all_changes)):
                cDNA_changes.append(str(rna_locations[ch_idx]) + ':' + ref_alleles[ch_idx] + '>' +
                                    alt_alleles[ch_idx])
        else:
            mutated_cdna = cdna_sequence
            for ch_idx, change in enumerate(all_changes):
                ref_allele = ref_alleles[ch_idx]
                alt_allele = alt_alleles[ch_idx]
                ref_len = len(ref_allele)
                alt_len = len(alt_allele)
                rna_location = rna_locations[ch_idx]
                cDNA_changes.append(str(rna_location) + ':' + ref_allele + '>' + alt_allele)
                rna_location += sequence_length_diff
                if (ref_allele != mutated_cdna[rna_location:rna_location + ref_len]):
                    print('Ref allele not matching the cDNA sequence, skipping haplotype!')
                    print(transcript_id + ' strand ' + t.strand + ' ' + change + ' cDNA: ' +
                          mutated_cdna[rna_location - 10:rna_location] + ' ' +
                          mutated_cdna[rna_location:rna_location + ref_len] + ' ' +
                          mutated_cdna[rna_location + ref_len:rna_location + ref_len + 10])
                    cDNA_changes = cDNA_changes[:-1]
                    validity_check = False
                    break
                mutated_cdna = mutated_cdna[:rna_location] + alt_allele + mutated_cdna[rna_location + ref_len:]
                sequence_length_diff += alt_len - ref_len
            if (not validity_check):
                continue

        exempt = None           # fix 'start', re-scanned: changes the UTR filter keeps (index -> True)
        if rescan:
            new_start = rescan_start(mutated_cdna, start_loc, changes_t)
            if new_start == -1:
                continue        # no ATG downstream: skipped, as the reference skips every lost start
            reading_frame = new_start % 3
            protein_start = new_start // 3
            # the changes from the lost start up to the new one decide where translation starts: the merge key keeps
            # them, with the changes that touch the lost start codon
            origin = min([start_loc] + [changes_t[k][0] for k in touching])
            exempt = [False] * len(all_changes)
            for k, (loc, ref_t, alt_t) in enumerate(changes_t):
                moved = loc + sum(len(a) - len(r) for j, (l2, r, a) in enumerate(changes_t)
                                  if j != k and (r or a) and l2 + len(r) <= loc)
                if k in touching or (loc >= origin and moved < new_start + 3):
                    exempt[k] = True

        # the changes on the protein level
        has_frameshift = False
        sequence_length_diff = 0
        non_synonymous = []     # per change, as all_protein_changes (for the fix 'start')
        if trimmed is None:
            order = range(len(all_changes))
            offsets = None
        else:
            # in the order the trimmed changes apply, each anchored allele where it starts in the haplotype's cDNA
            order = sorted(range(len(all_changes)),
                           key=lambda k: (trimmed[k][0], trimmed[k][0] + len(trimmed[k][1]), k))
            offsets = _offsets(rna_locations, trimmed)
            per_change = [None] * len(all_changes)
        for ch_idx in order:
            ref_len = len(ref_alleles[ch_idx])
            alt_len = len(alt_alleles[ch_idx])
            rna_location = rna_locations[ch_idx]
            rf_changes = []
            rf_conseq = []
            is_synonymous = []
            key = (rna_location, ref_len, reading_frame_ref, protein_start_ref)
            ref_side = ref_codons.get(key)
            if ref_side is None:
                ref_side = ref_codons[key] = get_affected_codons(cdna_sequence, rna_location, ref_len,
                                                                 reading_frame_ref, protein_start_ref)
            ref_alleles_protein, protein_location_ref = ref_side
            rna_location += sequence_length_diff if offsets is None else offsets[ch_idx]
            alt_alleles_protein, protein_location_alt = get_affected_codons(mutated_cdna, rna_location, alt_len,
                                                                            reading_frame, protein_start)
            for i, ref_allele_protein in enumerate(ref_alleles_protein):
                alt_allele_protein = alt_alleles_protein[i]
                loc_ref = protein_location_ref[i]
                loc_alt = protein_location_alt[i]
                if (dna_var_types[ch_idx] == 'splice'):
                    rf_conseq.append('splice_variant')
                elif (dna_var_types[ch_idx] == 'start_lost'):
                    rf_conseq.append('start_lost')
                elif (ref_allele_protein == alt_allele_protein) and not (fix_frameshift and frameshifts[ch_idx]):
                    rf_conseq.append('synonymous')
                elif ('*' in ref_allele_protein) and not ('*' in alt_allele_protein):
                    rf_conseq.append('stop_lost')
                elif not ('*' in ref_allele_protein) and ('*' in alt_allele_protein):
                    rf_conseq.append('stop_gained')
                elif (dna_var_types[ch_idx] == 'SNP'):
                    rf_conseq.append('SAV')
                elif (dna_var_types[ch_idx] == 'indel'):
                    if (frameshifts[ch_idx]):
                        rf_conseq.append('frameshift')
                    else:
                        rf_conseq.append('inframe_indel')
                is_synonymous.append(ref_allele_protein == alt_allele_protein and
                                     not (rescan and dna_var_types[ch_idx] == 'start_lost') and
                                     not (fix_frameshift and frameshifts[ch_idx]))
                protein_change = str(loc_ref) + ':' + ref_allele_protein + '>' + str(loc_alt) + ':' + \
                    alt_allele_protein
                if (frameshifts[ch_idx]):
                    protein_change += "(+fs)"
                elif (has_frameshift):
                    protein_change += "(fs)"
                    rf_conseq[-1] += '_after_fs'
                rf_changes.append(protein_change)
                has_frameshift = has_frameshift or (frameshifts[ch_idx] and (loc_ref >= 0))
            if offsets is not None:
                per_change[ch_idx] = (not all(is_synonymous), "|".join(rf_changes), '|'.join(rf_conseq))
                continue
            if not all(is_synonymous):
                protein_changes.append("|".join(rf_changes))
            all_protein_changes.append("|".join(rf_changes))
            non_synonymous.append(not all(is_synonymous))
            prot_var_types.append('|'.join(rf_conseq))
            sequence_length_diff += alt_len - ref_len
        if offsets is not None:
            for non_synonymous_change, protein_change, conseq in per_change:
                if non_synonymous_change:
                    protein_changes.append(protein_change)
                all_protein_changes.append(protein_change)
                non_synonymous.append(non_synonymous_change)
                prot_var_types.append(conseq)

        spl_junctions_affected_str = ';'.join([str(x) for x in spl_junctions_affected])
        if (len(spl_junctions_affected_str) == 0):
            spl_junctions_affected_str = '-'
        haplotype_id = index

        if (reading_frame > -1):
            protein_seq = translate(mutated_cdna[reading_frame:])
            cdna_start = reading_frame + protein_start * 3
            if (settings.ignore_utr):
                first_stop = protein_seq.find('*', protein_start)
                cdna_stop = reading_frame + first_stop * 3
                if (first_stop == -1):
                    first_stop = len(protein_seq)
                def in_cds(change):
                    key = (change, protein_start, first_stop)
                    inside = in_cds_cache.get(key)
                    if inside is None:
                        inside = in_cds_cache[key] = check_protein_allele(change, protein_start, first_stop)
                    return inside
                variant_filter = [in_cds(change) for change in all_protein_changes]
                if exempt is not None:
                    variant_filter = [exempt[idx] or variant_filter[idx] for idx in range(len(variant_filter))]
                    protein_changes = [ch for idx, ch in enumerate(all_protein_changes)
                                       if variant_filter[idx] and non_synonymous[idx]]
                all_vcf_IDs = [vcf_id for idx, vcf_id in enumerate(all_vcf_IDs) if variant_filter[idx]]
                all_changes = [ch for idx, ch in enumerate(all_changes) if variant_filter[idx]]
                all_AFs = [af for idx, af in enumerate(all_AFs) if variant_filter[idx]]
                cDNA_changes = [ch for idx, ch in enumerate(cDNA_changes) if variant_filter[idx]]
                all_protein_changes = [ch for idx, ch in enumerate(all_protein_changes) if variant_filter[idx]]
                prot_var_types = [v for idx, v in enumerate(prot_var_types) if variant_filter[idx]]
                if exempt is None:
                    protein_changes = [ch for ch in protein_changes if in_cds(ch)]
                if len(protein_changes) == 0:
                    continue
                protein_seq = protein_seq[protein_start:first_stop]
                mutated_cdna = mutated_cdna[cdna_start:cdna_stop]
                protein_start = 0
                cdna_start = 0

            haplo_hash = str(hash(';'.join(all_vcf_IDs)))
            frequency = row.frequency(cohort)
            if haplo_hash in local_result_data:
                merged = local_result_data[haplo_hash]
                merged.fields[16] += row.count
                merged.fields[17] += frequency
                merged.parts.append(row)
            else:
                local_result_data[haplo_hash] = _Merged([
                    transcript_id, settings.chromosome, t.biotype, haplotype_id, ';'.join(all_vcf_IDs),
                    ";".join(all_changes), ";".join(all_AFs), ';'.join(cDNA_changes), ';'.join(all_protein_changes),
                    ';'.join(prot_var_types), ';'.join(protein_changes), reading_frame, protein_start,
                    t.start_codon is None, start_lost, spl_junctions_affected_str, row.count, frequency, None, None,
                    None], row)
                seq_hash = hash(protein_seq)
                entry = protein_sequences.get(seq_hash)
                if entry is not None:
                    entry['haplotypes'].append(haplotype_id)
                    entry['rfs'].append(str(reading_frame))
                else:
                    protein_sequences[seq_hash] = {'haplotypes': [haplotype_id], 'sequence': protein_seq,
                                                   'start': protein_start, 'rfs': [str(reading_frame)]}
                if (output_cdnas):
                    cdna_hash = hash(mutated_cdna)
                    entry = cdna_sequences.get(cdna_hash)
                    if entry is not None:
                        entry['haplotypes'].append(haplotype_id)
                    else:
                        cdna_sequences[cdna_hash] = {'haplotypes': [haplotype_id], 'sequence': mutated_cdna,
                                                     'start': cdna_start}

        elif (len(protein_changes) > 0):
            for rf in range(0, 3):
                protein_seq = translate(mutated_cdna[rf:])
                seq_hash = hash(protein_seq)
                entry = protein_sequences.get(seq_hash)
                if entry is not None:
                    entry['haplotypes'].append(haplotype_id)
                    entry['rfs'].append(str(rf))
                else:
                    protein_sequences[seq_hash] = {'haplotypes': [haplotype_id], 'sequence': protein_seq,
                                                   'start': protein_start, 'rfs': [str(rf)]}
            if (output_cdnas):
                cdna_hash = hash(mutated_cdna)
                entry = cdna_sequences.get(cdna_hash)
                if entry is not None:
                    entry['haplotypes'].append(haplotype_id)
                else:
                    cdna_sequences[cdna_hash] = {'haplotypes': [haplotype_id], 'sequence': mutated_cdna, 'start': -1}
            haplo_hash = hash(';'.join(all_vcf_IDs))
            # not merged: a later haplotype with the same IDs replaces the row (the reference's behaviour)
            local_result_data[haplo_hash] = _Merged([
                transcript_id, settings.chromosome, t.biotype, haplotype_id, ';'.join(all_vcf_IDs),
                ";".join(all_changes), ";".join(all_AFs), ';'.join(cDNA_changes), ';'.join(all_protein_changes),
                ';'.join(prot_var_types), ';'.join(protein_changes), reading_frame, protein_start,
                t.start_codon is None, start_lost, spl_junctions_affected_str, row.count, row.frequency(cohort),
                None, None, None], row)

    # filter by frequency or count, after merging
    if (settings.min_freq != -1):
        kept = [m for m in local_result_data.values() if (m.fields[17] >= settings.min_freq)]
    else:
        kept = [m for m in local_result_data.values() if (m.fields[16] >= settings.min_count)]
    result_table = []
    for m in kept:
        pops = [row.frequency_strings(cohort) for row in m.parts]
        m.fields[18] = add_population_freqs([p[0] for p in pops])
        m.fields[19] = add_population_freqs([p[1] for p in pops])
        m.fields[20] = ';'.join([row.sample_string(cohort) for row in m.parts])
        result_table.append(m.fields)

    included = {fields[3] for fields in result_table}
    proteins = []
    for seq_hash in sorted(protein_sequences):
        seq = protein_sequences[seq_hash]
        keep = [(hap in included) for hap in seq['haplotypes']]
        seq['haplotypes'] = [hap for i, hap in enumerate(seq['haplotypes']) if keep[i]]
        seq['rfs'] = [rf for i, rf in enumerate(seq['rfs']) if keep[i]]
        if seq['haplotypes']:
            proteins.append(seq)
    cdnas = []
    if (output_cdnas):
        for cdna_hash in sorted(cdna_sequences):
            seq = cdna_sequences[cdna_hash]
            seq['haplotypes'] = [hap for hap in seq['haplotypes'] if (hap in included)]
            if seq['haplotypes']:
                cdnas.append(seq)
    return result_table, proteins, cdnas

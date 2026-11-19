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

The reference's merge keys are kept: str(hash(IDs)) with a known reading frame (merge), hash(IDs) without (a later
haplotype replaces an earlier one, IMPLEMENTATION_NOTES section 7.4). Haplotype IDs are the rows' local indices; the
caller turns them into `<prefix>_<hex>`.
"""
import re

import numpy as np
from Bio.Data import CodonTable
from Bio.Seq import Seq

from modules.fallback import Unsupported
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
    def __init__(self, chromosome, force_rf, min_freq, min_count, ignore_utr, skip_start_loss, output_cdnas):
        self.chromosome = chromosome
        self.force_rf = force_rf
        self.min_freq = min_freq
        self.min_count = min_count
        self.ignore_utr = ignore_utr
        self.skip_start_loss = skip_start_loss
        self.output_cdnas = output_cdnas


class _Merged:
    """A row of the result: the reference's 21 columns, with count and frequency added up as haplotypes merge, and
    the haplotypes behind its population and sample strings."""
    __slots__ = ('fields', 'parts')

    def __init__(self, fields, part):
        self.fields = fields
        self.parts = [part]


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
                if reverse_strand:
                    ref_allele = reverse_complement(ref_allele)
                    alt_allele = reverse_complement(alt_allele)
                    rna_location = len(cdna_sequence) - rna_location - ref_len
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

            if (t.start_codon is not None) and (reading_frame > -1):
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

        # the haplotype's cDNA
        mutated_cdna = None
        previous_end = 0
        ordered = True
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

        # the changes on the protein level
        has_frameshift = False
        sequence_length_diff = 0
        for ch_idx, change in enumerate(all_changes):
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
            rna_location += sequence_length_diff
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
                elif (ref_allele_protein == alt_allele_protein):
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
                is_synonymous.append(ref_allele_protein == alt_allele_protein)
                protein_change = str(loc_ref) + ':' + ref_allele_protein + '>' + str(loc_alt) + ':' + \
                    alt_allele_protein
                if (frameshifts[ch_idx]):
                    protein_change += "(+fs)"
                elif (has_frameshift):
                    protein_change += "(fs)"
                    rf_conseq[-1] += '_after_fs'
                rf_changes.append(protein_change)
                has_frameshift = has_frameshift or (frameshifts[ch_idx] and (loc_ref >= 0))
            if not all(is_synonymous):
                protein_changes.append("|".join(rf_changes))
            all_protein_changes.append("|".join(rf_changes))
            prot_var_types.append('|'.join(rf_conseq))
            sequence_length_diff += alt_len - ref_len

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
                        inside = in_cds_cache[key] = _check_protein_allele(change, protein_start, first_stop)
                    return inside
                variant_filter = [in_cds(change) for change in all_protein_changes]
                all_vcf_IDs = [vcf_id for idx, vcf_id in enumerate(all_vcf_IDs) if variant_filter[idx]]
                all_changes = [ch for idx, ch in enumerate(all_changes) if variant_filter[idx]]
                all_AFs = [af for idx, af in enumerate(all_AFs) if variant_filter[idx]]
                cDNA_changes = [ch for idx, ch in enumerate(cDNA_changes) if variant_filter[idx]]
                all_protein_changes = [ch for idx, ch in enumerate(all_protein_changes) if variant_filter[idx]]
                prot_var_types = [v for idx, v in enumerate(prot_var_types) if variant_filter[idx]]
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

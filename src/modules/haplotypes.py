"""The haplotypes of a transcript (IMPLEMENTATION_NOTES section 3.4): what the reference's get_haplotypes computes for
one transcript, from the copies that carry each record instead of a table with one column per sample.

The reference builds, for every haplotype copy of every included sample, the comma-separated list of the table rows
whose ALT allele the copy carries ('REF' when none), and groups the copies by that string, kept sorted as strings.
Here each copy's rows are hashed (two random 64-bit weights per row, summed), copies are grouped by (number of rows,
hashes), every member of a group is then checked row by row against the group's first member (a mismatch, i.e. a hash
collision, raises Unsupported), and the groups are sorted by the reference's string. Members stay in ascending copy
order, which is the reference's sample order.

Per group the reference's own steps follow: the changes `POS:REF>ALT`, allele frequencies from INFO and VCF IDs; a
group with conflicting variants (legacy cluster_conflicting_mutations) gives no row and its samples go to the removed
list; the changes and frequencies are sorted by POS (the IDs are not). Counts, frequencies and the population and
sample strings are computed only for the rows that need them (Row methods): at 500,000 samples the reference
haplotype alone has a million members.
"""
import numpy as np

from modules.fallback import Unsupported
from modules.legacy_import import legacy

_cluster = legacy.get_haplotypes.cluster_conflicting_mutations


class Cohort:
    """The included samples (the VCF's samples that are in the metadata, in VCF order) and their metadata, as
    get_gene_haplotypes uses them. rank r has copies 2r (':1') and 2r + 1 (':2')."""

    def __init__(self, names, male, populations, superpopulations, is_x, par1_to, par2_from):
        """names: included sample names; male: bool per sample; populations / superpopulations: (code per sample,
        {code: number of samples}, {code: number of male samples}) as the reference counts them."""
        self.names = list(names)
        self.n = len(self.names)
        self.male = np.asarray(male, dtype=bool)
        self.n_male = int(self.male.sum())
        self.is_x = is_x
        self.par1_to = par1_to
        self.par2_from = par2_from
        self.labels = [name + suffix for name in self.names for suffix in (':1', ':2')]
        self.pops = [self._groups(*populations), self._groups(*superpopulations)]

    @classmethod
    def from_metadata(cls, samples_df, sample_ids, is_x, par1_to, par2_from):
        """The cohort as get_gene_haplotypes derives it from the metadata table (pandas, as read by prohap.py) and the
        included sample names. Raises Unsupported where the reference would stop or the fast path would differ: no
        sample, a sample with two metadata rows, a population code that is missing, not text, or contains ':' or ';'
        (the reference splits frequency strings on them)."""
        if not sample_ids:
            raise Unsupported('no sample of the VCF is in the metadata')
        try:
            sample_info = samples_df[samples_df['Sample name'].isin(sample_ids)]
            counts = []
            for column in ('Population code', 'Superpopulation code'):
                count = sample_info[[column, 'Sample name']].groupby(column).count()
                male_count = sample_info[sample_info['Sex'] == 'male'][[column, 'Sample name']].groupby(column).count()
                counts.append(({k: int(v) for k, v in count['Sample name'].items()},
                               {k: int(v) for k, v in male_count['Sample name'].items()}))
            sample_info = sample_info.set_index('Sample name')
            if not sample_info.index.is_unique:
                raise Unsupported('a sample has more than one row in the metadata')
            sex = sample_info['Sex'].to_dict()
            population = sample_info['Population code'].to_dict()
            superpopulation = sample_info['Superpopulation code'].to_dict()
        except KeyError as e:
            raise Unsupported('the metadata lacks a column: %s' % e)
        pops = [population[s] for s in sample_ids]
        superpops = [superpopulation[s] for s in sample_ids]
        if not all(isinstance(c, str) for c in pops + superpops):
            raise Unsupported('a population or superpopulation code that is missing or not text')
        if any((':' in c) or (';' in c) for c in set(pops + superpops)):
            raise Unsupported("a population or superpopulation code contains ':' or ';'")
        return cls(sample_ids, [sex[s] == 'male' for s in sample_ids], (pops,) + counts[0], (superpops,) + counts[1],
                   is_x, par1_to, par2_from)

    @staticmethod
    def _groups(code_of_sample, count, male_count):
        codes = list(dict.fromkeys(code_of_sample))
        index = {code: i for i, code in enumerate(codes)}
        of_rank = [index[c] for c in code_of_sample]
        autosomal = [count[c] * 2 for c in codes]
        x = [male_count.get(c, 0) + (count[c] - male_count.get(c, 0)) * 2 for c in codes]
        return codes, of_rank, np.array(of_rank, dtype=np.int32), autosomal, x

    def is_autosomal(self, start, end):
        """The reference's test: not chromosome X, or the transcript lies in a pseudo-autosomal region."""
        return (not self.is_x) or ((start < self.par1_to) and (end <= self.par1_to)) or \
            ((start >= self.par2_from) and (end > self.par2_from))

    def total(self, autosomal):
        """count_freq's denominator."""
        if autosomal:
            return self.n * 2
        return self.n_male + (self.n - self.n_male) * 2

    def frequency_strings(self, members, autosomal):
        """count_freq_pop and count_freq_superpop for a haplotype whose copies are `members` (ascending)."""
        result = []
        small = len(members) <= 64
        if small:
            ranks = (members >> 1).tolist()
        for codes, of_rank, of_rank_array, totals_a, totals_x in self.pops:
            totals = totals_a if autosomal else totals_x
            if small:
                counts = {}         # in order of first appearance, as the reference's dict
                for r in ranks:
                    p = of_rank[r]
                    counts[p] = counts.get(p, 0) + 1
                pairs = counts.items()
            else:
                unique, first, n = np.unique(of_rank_array[members >> 1], return_index=True, return_counts=True)
                order = np.argsort(first, kind='stable')
                pairs = zip(unique[order].tolist(), n[order].tolist())
            parts = []
            for p, n_p in pairs:
                total = totals[p]
                if total > 0:
                    parts.append(codes[p] + ':{:.5f}'.format(n_p / total))
            result.append(';'.join(parts))
        return result

    def sample_string(self, members):
        labels = self.labels
        return ';'.join([labels[c] for c in members.tolist()])


class Row:
    """One row of the reference's haplotype table for a transcript: the records its haplotype carries (None for the
    reference haplotype) and the copies that have it. The changes, allele frequencies and VCF IDs are joined when
    asked for, as get_haplotypes joins them: changes and frequencies sorted by POS (stable), IDs in table order."""
    __slots__ = ('rows', 'columns', 'members', 'count', 'autosomal')

    def __init__(self, rows, columns, members, count, autosomal):
        self.rows = rows
        self.columns = columns
        self.members = members
        self.count = count
        self.autosomal = autosomal

    def _by_pos(self):
        return sorted(self.rows.tolist(), key=self.columns.positions.__getitem__)

    @property
    def changes(self):
        if self.rows is None:
            return 'REF'
        change = self.columns.changes
        return ';'.join([change[i] for i in self._by_pos()])

    @property
    def afs(self):
        if self.rows is None:
            return ''
        af = self.columns.afs
        return ';'.join([af[i] for i in self._by_pos()])

    @property
    def ids(self):
        if self.rows is None:
            return ''
        vcf_id = self.columns.ids
        return ';'.join([vcf_id[i] for i in self.rows.tolist()])

    def frequency(self, cohort):
        total = cohort.total(self.autosomal)
        if total == 0:
            return 0
        return self.count / total

    def frequency_strings(self, cohort):
        return cohort.frequency_strings(self.members, self.autosomal)

    def sample_string(self, cohort):
        return cohort.sample_string(self.members)


class _Columns:
    """Per record of a transcript: POS, `POS:REF>ALT`, the AF text and str() of the ID."""
    __slots__ = ('positions', 'changes', 'afs', 'ids')

    def __init__(self, positions, changes, afs, ids):
        self.positions = positions
        self.changes = changes
        self.afs = afs
        self.ids = ids


def allele_frequency(info):
    """get_haplotypes' AF column: the AF= value from INFO, or '-1' when INFO does not contain 'AF'."""
    if 'AF' in info:
        try:
            return info.split('AF=')[1].split(';')[0].split(maxsplit=1)[0]
        except IndexError:
            raise Unsupported('INFO contains AF but no value the reference can read: ' + info[:200])
    return '-1'


def _involved(positions, lengths):
    """Records whose span [POS, POS + len(REF)] meets another record's: the only ones that can conflict (the
    reference's clusters join spans that overlap or touch)."""
    n = len(positions)
    if n < 2:
        return np.zeros(n, dtype=bool)
    start = np.asarray(positions, dtype=np.int64)
    end = start + np.asarray(lengths, dtype=np.int64)
    order = np.argsort(start, kind='stable')
    s = start[order]
    e = end[order]
    previous_max = np.maximum.accumulate(np.concatenate(([np.iinfo(np.int64).min], e[:-1])))
    next_start = np.concatenate((s[1:], [np.iinfo(np.int64).max]))
    involved = np.empty(n, dtype=bool)
    involved[order] = (previous_max >= s) | (next_start <= e)
    return involved


_SCRATCH = {}                   # per process: arrays of one entry per haplotype copy, kept zero between calls


def group_copies(copies, n_copies):
    """Copies grouped by the set of records they carry: (record arrays of the groups, member arrays), in no
    particular order. copies: one array of carrier copies (0 .. n_copies - 1, no repeats) per record.

    Each copy's records are hashed into two sums of random 64-bit weights, and copies are grouped by (number of
    records, sums). The grouping is then checked exactly: for every record, each group must have all or none of its
    members among the record's carriers; a hash collision would break this and raises Unsupported."""
    if _SCRATCH.get('n') != n_copies:
        _SCRATCH.clear()
        _SCRATCH.update(n=n_copies, counts=np.zeros(n_copies, dtype=np.int32),
                        h1=np.zeros(n_copies, dtype=np.uint64), h2=np.zeros(n_copies, dtype=np.uint64),
                        group=np.zeros(n_copies, dtype=np.int64))
    counts, h1, h2, group = _SCRATCH['counts'], _SCRATCH['h1'], _SCRATCH['h2'], _SCRATCH['group']
    rng = np.random.default_rng(0x5eed)
    weights = rng.integers(0, np.iinfo(np.uint64).max, size=(2, len(copies)), dtype=np.uint64, endpoint=True)
    for r, c in enumerate(copies):
        if len(c):
            counts[c] += 1
            h1[c] += weights[0, r]
            h2[c] += weights[1, r]
    carriers = np.flatnonzero(counts)
    try:
        if len(carriers) == 0:
            return [], []
        n, a, b = counts[carriers], h1[carriers], h2[carriers]
        by_key = carriers[np.lexsort((b, a, n))]   # stable: members of a group stay in ascending copy order
        n, a, b = counts[by_key], h1[by_key], h2[by_key]
        first = np.empty(len(by_key), dtype=bool)
        first[0] = True
        first[1:] = (n[1:] != n[:-1]) | (a[1:] != a[:-1]) | (b[1:] != b[:-1])
        del n, a, b
        group_starts = np.flatnonzero(first)
        size = np.diff(np.append(group_starts, len(by_key)))
        group[by_key] = np.cumsum(first) - 1
        pair_group = []
        pair_record = []
        for r, c in enumerate(copies):
            if len(c) == 0:
                continue
            groups, members = np.unique(group[c], return_counts=True)
            if not np.array_equal(members, size[groups]):
                raise Unsupported('two different haplotypes with the same hash')
            pair_group.append(groups)
            pair_record.append(np.full(len(groups), r, dtype=np.int32))
        pair_group = np.concatenate(pair_group)
        order = np.argsort(pair_group, kind='stable')        # records stay in table order within a group
        pair_record = np.concatenate(pair_record)[order]
        bounds = np.searchsorted(pair_group[order], np.arange(len(group_starts) + 1))
        rows = [pair_record[bounds[g]:bounds[g + 1]] for g in range(len(group_starts))]
        members = [by_key[group_starts[g]:group_starts[g] + size[g]] for g in range(len(group_starts))]
        return rows, members
    finally:
        counts[carriers] = 0
        h1[carriers] = 0
        h2[carriers] = 0


def transcript_rows(records, ids, autosomal, cohort):
    """get_haplotypes for one transcript with at least one record. records: sweep records in VCF order; ids: str() of
    their IDs as the reference reads them. Returns (rows in the reference's order, copies removed for conflicting
    variants in the reference's order)."""
    copies = [r[6] for r in records]
    if not autosomal:
        # males have one copy outside the pseudo-autosomal regions: their second copy is not counted at all
        male = cohort.male
        copies = [c[~((c & 1).astype(bool) & male[c >> 1])] for c in copies]
    group_rows, group_members = group_copies(copies, 2 * cohort.n)
    keys = [','.join(map(str, rows.tolist())) for rows in group_rows]
    n_included = cohort.n * 2 if autosomal else cohort.n + (cohort.n - cohort.n_male)
    n_ref = n_included - sum(len(m) for m in group_members)
    order = sorted(range(len(keys)), key=keys.__getitem__)
    entries = [(keys[g], g) for g in order]
    del keys
    if n_ref > 0:
        entries.append(('REF', None))
        entries.sort(key=lambda e: e[0])

    positions = [r[1] for r in records]
    # the reference reads every carried record's AF for every haplotype; any failure would stop it
    carried = [len(c) > 0 for c in copies]
    columns = _Columns(positions, [r[0] + ':' + r[3] + '>' + r[4] for r in records],
                       [allele_frequency(r[5]) if carried[i] else None for i, r in enumerate(records)], ids)
    involved = _involved(positions, [len(r[3]) for r in records])
    result = []
    removed = []
    for key, g in entries:
        if g is None:
            result.append(Row(None, columns, None, n_ref, autosomal))
            continue
        rows = group_rows[g]
        if involved[rows].sum() >= 2:
            changelist = [[positions[i], records[i][3], records[i][4]] for i in rows.tolist()]
            clusters = _cluster([{'change': ch, 'id': i} for i, ch in enumerate(changelist)])
            if any(len(cluster) > 1 for cluster in clusters):
                removed.append(group_members[g])
                continue
        result.append(Row(rows, columns, group_members[g], len(group_members[g]), autosomal))
    removed = np.concatenate(removed) if removed else np.zeros(0, dtype=np.int64)
    return result, removed


def reference_only_row(cohort):
    """get_haplotypes' row for a transcript without records."""
    return Row(None, None, None, cohort.n * 2, True)

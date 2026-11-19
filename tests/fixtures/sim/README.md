# Simulated fixtures

Test data for comparing ProHap-fast (`src/`) with the reference code (`legacy/src`) beyond what the 1kG sample data
covers: larger cohorts, IDs that pandas retypes, chromosome X with haploid or mirrored males, very long transcripts
and big tables.

The data is simulated by [`bench/simulate.py`](../../../bench/simulate.py) (VCFs, sample metadata, transcript list)
and [`bench/simulate_fasta.py`](../../../bench/simulate_fasta.py) (protein FASTA) from fixed seeds, so it is not in
the repository: `make_data.py DATASET` writes it into `cache/` (git-ignored) the first time it is needed and again
whenever the definition in `make_data.py`, the simulator or the Ensembl release changes. Simulating a dataset takes
seconds; the 400,000-record FASTA about 10 s and 140 MB.

How the simulator works (in short; the details are in its docstring): it takes a block of consecutive MANE Select
transcripts with a start codon, puts variants on their exons with REF bases from the Ensembl cDNA, and gives them to
the samples' haplotype copies. Common variants (allele frequency 1-50%) are carried by 40 founder haplotypes per
50 kb, so many copies share a haplotype as in a real cohort; rare ones (below 1%) by a few random copies. The number
of variants grows with the cohort as in real data. 6% are deletions, 6% insertions, 6% have a second ALT allele, and
some deletions have an SNV inside them on the same copies, which ProHap reports as conflicting (the removed-samples
log).

## Fixtures

Each `NAME.yaml` is a workflow config; its data is the dataset named before the first `-`
(`mix-template.yaml` uses `mix`). `tests/legacy_vs_fast.sh` runs them as `sim/NAME`, in tiers.

| Fixture | Data | What it tests |
| --- | --- | --- |
| **small tier** (default, with the two 1kG configs) | | |
| `mix-threshold0` | 2,500 samples; about 30 kb of exons on each of chromosomes 20 (pandas-ambiguous IDs), 21 (rs IDs), 22 (`.` IDs, two VCF chunks) and X (5 PAR1 transcripts, one across the PAR1 boundary, 3 beyond it, haploid males); 5 VCF samples not in the metadata and 5 metadata rows not in the VCF | every ID mode; IDs such as `2e5`, `0012`, `NA`, `True` that pandas reads as numbers, missing values or booleans; merges keyed on equal IDs; haploid genotypes; PAR handling; conflicting variants. Thresholds 0, cDNA FASTA, plain outputs, 1 thread |
| `mix-template` | as above | the thresholds of `config_tre_template.yaml` (variants with AF ≥ 0.01, haplotypes seen ≥ 10 times), gzipped outputs, 4 threads |
| **medium tier** | | |
| `n2500`, `n5000`, `n10000` | the benchmark block (60 transcripts of chromosome 22 from ENST00000402142, 217 kb of exons) at 2,500, 5,000 and 10,000 samples; `.`, `.` and pandas-ambiguous IDs; `n5000` adds chromosome X with rs IDs and mirrored males (`1\|1`) | cohort sizes at which pandas reads the per-transcript tables in chunks of 256, 128 and 64 rows; 1, 2 and 4 threads; plain, gzipped cDNA and gzipped outputs |
| `whole22` | every MANE transcript of chromosome 22 (478) and 20 of chromosome 21, 2,500 samples | a haplotype table of more than 32,768 rows, which pandas reads in chunks in `merge_tables.py` and `haplo_extract_sample_names.py` |
| `dupseq` (not a workflow) | 400,000 protein records: the MANE reference proteins, the contaminants and haplotype proteins with duplicates | `merge_duplicate_seq.py` and `remove_UTR_seq.py` at a size where the reference's sorted-list insert is slow |
| **large tier** | | |
| `n26000` | the benchmark block at 26,000 samples | 32-row pandas chunks, the cohort size of the TRE runs |
| `n100000` | 100,000 samples, 3 short transcripts on each of chromosomes 21 (pandas-ambiguous IDs) and 22 (`.` IDs) | 8-row pandas chunks, where the type of an ID depends on its neighbours |
| `ttn` | TTN (ENST00000589042: 363 exons, 109 kb, − strand) at 2,500 samples | one transcript with 3,400 variants: long tables, many haplotypes, a long cDNA |

## Running them

```bash
conda activate prohap-fast
export PROHAP_REFS=~/prohap-refs
tests/legacy_vs_fast.sh                   # small tier
tests/legacy_vs_fast.sh medium            # or large, all, or fixtures by name: sim/n5000 dupseq
python3 tests/fixtures/sim/make_data.py n5000   # only generate (prints the directory)
```

Measured run times are in `tests/legacy_vs_fast.sh`'s header. The reference code's cost grows faster than the
cohort, so the large tier is for changes to the code that the cohort size reaches (pandas emulation, haplotype
grouping) and for releases. Run the larger tiers when a step's "done when" asks for them.

`pytest tests/test_bench.py` tests the simulator and checks that every config here names a dataset of
`make_data.py` with the same chromosomes, and that every tier of `tests/legacy_vs_fast.sh` names a fixture.

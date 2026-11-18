# 1kG fixture

The test data for comparing ProHap-fast (`src/`) with the reference code (`legacy/src`):

- **VCFs**: the phased 1000 Genomes data shipped with ProHap, `sample_data/sample_1kGP_common_global.tar.gz`
  (common variants, AF ≥ 0.01, 2,548 samples, every record with ID `.`). The fixture uses chromosomes 1 (two chunks:
  462 and 309 records), 21 (102), 22 (187) and X (7, all in PAR1).
- **Sample metadata**: `igsr_samples.tsv` in the repository root.
- **Transcripts**: `transcripts.csv`, 109 MANE Select transcripts of Ensembl release 113 on those chromosomes. On
  each VCF chunk it keeps 20 transcripts that overlap a variant and 5 that do not, spread along the chromosome; on X,
  the 4 that overlap a variant and 5 that do not. `make_transcripts.py` wrote it and says why: the reference code
  spends about 0.1 s on every transcript, variants or not, and only one MANE transcript in five overlaps a variant
  here, so the whole chromosomes would take minutes per run and add little.
- **References**: the pre-built references of release 113 (`$PROHAP_REFS/ensembl.v113`).

Two configs:

| Config | Settings |
| --- | --- |
| `threshold0.yaml` | every haplotype kept (`haplo_min_freq: 0`); cDNA FASTA on; plain-text outputs; 1 thread |
| `template.yaml` | the thresholds of `config_tre_template.yaml` (`haplo_min_freq: -1`, `haplo_min_count: 10`); gzipped table and FASTA; 4 threads |

## The check

```bash
conda activate prohap-fast            # or any environment with the ProHap packages and snakemake
PROHAP_REFS=~/prohap-refs tests/legacy_vs_fast.sh
```

For each config, the script runs the workflow with the reference code, then with ProHap-fast, on the same inputs,
and compares every output file after decompression. It prints each file as IDENTICAL or DIFFERENT (with the first
differing lines) and exits with status 0 only when everything is identical. About a minute and a half on 8 cores:
roughly 20 s per run, most of it the steps every run has (cleaning the reference proteome, merging with it).

`tests/test_identity.py` (pytest) runs `verify`, the user-facing version of the same check, on these configs, plus
a sample and transcript subset and a negative test (a copy of `src/` that changes one character must be reported
DIFFERENT). It takes about 3 minutes on 4 cores.

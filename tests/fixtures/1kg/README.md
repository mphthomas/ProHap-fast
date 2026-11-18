# 1kG fixture

Phased 1000 Genomes data shipped with ProHap: `sample_data/sample_1kGP_common_global.tar.gz` (common variants,
AF ≥ 0.01, 2,548 samples, every record with ID `.`) and the sample metadata `igsr_samples.tsv`. The fixture uses
chromosomes 2 (521 records), 21 (102), 22 (187) and X (7, all in PAR1), with the MANE Select transcripts of
Ensembl release 113 from the pre-built references.

| Config | Settings |
| --- | --- |
| `threshold0.yaml` | every haplotype kept (`haplo_min_freq: 0`); cDNA FASTA on; plain-text outputs; 1 thread |
| `template.yaml` | the thresholds of `config_tre_template.yaml` (`haplo_min_freq: -1`, `haplo_min_count: 10`); gzipped table and FASTA; 4 threads |

`tests/test_identity.py` sets up a work directory for each (the VCFs, metadata, `crap.fasta` and the config as
`config.yaml`) and runs `verify` there, which must report every output IDENTICAL. It also runs a subset of
`threshold0.yaml` on chromosomes 1 and 22 (chromosome 1 comes in two chunks) with `--samples 300 --transcripts 30`,
and a negative test: a copy of `src/` that changes one character of the removed-samples log must make `verify`
report that file DIFFERENT and exit with status 1.

To run a fixture by hand:

```bash
mkdir -p /tmp/1kg/sample_data && cd /tmp/1kg
tar -xzf ~/GitHub/ProHap-fast/sample_data/sample_1kGP_common_global.tar.gz -C sample_data
ln -s ~/GitHub/ProHap-fast/igsr_samples.tsv ~/GitHub/ProHap-fast/crap.fasta .
cp ~/GitHub/ProHap-fast/tests/fixtures/1kg/threshold0.yaml config.yaml
PROHAP_REFS=~/prohap-refs python3 ~/GitHub/ProHap-fast/verify/verify.py --cores 4
```

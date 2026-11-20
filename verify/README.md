# verify: check ProHap-fast against the reference on your own data

ProHap-fast promises output that is byte-identical to its reference implementation, the frozen copy of ProHap in
`legacy/`. `verify` lets you check that promise on your own cohort before you rely on it: it runs the workflow twice
on the same inputs, once with the reference and once with ProHap-fast, and compares every output.

## Running it

Run it from the directory that holds your `config.yaml`, the same directory a normal run uses. In the release image
it is the `verify` command, with the same mounts as a normal run (see the [README](../README.md#running-with-docker)):

```bash
cd /path/to/workdir
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" -v /data:/data:ro \
    mphthomas/prohap-fast:0.3.0-ensembl.v113 verify --cores 8 --samples 2000 --transcripts 50
singularity run --bind /data prohap-fast.sif verify --cores 8 --samples 2000 --transcripts 50
```

Without a container, with the ProHap environment active:

```bash
cd /path/to/workdir
python3 /path/to/ProHap-fast/verify/verify.py --cores 8 --samples 2000 --transcripts 50
```

Your `config.yaml` stays as it is. `verify` writes only under `verify_runs/<timestamp>/` in the current directory:

- `legacy/` and `fast/`: one work directory per run (an "arm"), each with its own copy of `config.yaml` that differs
  only in `implementation`. Relative input paths in your config are resolved against the current directory, and
  every output path moves into the arm's directory, so neither run touches your normal outputs.
- `input/`: the cut-down VCFs, metadata and transcript list, when you ask for a subset.
- `legacy.log`, `fast.log`: Snakemake's output for each run.
- `report.tsv`: the comparison, with full sha256 sums.

When every output is identical, `verify` deletes everything except `report.tsv`; `--keep` keeps it all. When an
output differs or a run fails, it keeps everything.

## Options

| Option | Meaning |
| --- | --- |
| `--cores N` | cores for each run (default 1). The runs happen one after the other. |
| `--chromosomes 21,22` | only these chromosomes of the config's `chromosomes`. |
| `--samples N` | N samples drawn at random from the VCF header (seeded, so a second run picks the same ones). |
| `--seed S` | the seed for `--samples` (default 0). |
| `--sample-list FILE` | the samples named in `FILE`, one per line, instead of a random draw. |
| `--transcripts N` | the first N transcripts of each chromosome in the transcript list. |
| `--keep` | keep the run directory when every output is identical. |
| `--dry-run` | prepare the run directory, check the input and dry-run both arms, without running them. |

Each subset is applied to both arms alike. `--samples` and `--sample-list` cut every VCF chunk to the chosen
samples and `sample_metadata_file` to the matching rows. The cut chunks keep their names and differ from yours only
by the other samples' columns: every header line, record, INFO value and line end stays as it is.

## Choosing a subset

The reference is slow at cohort scale: about 2,000 CPU seconds per 60 transcripts at 26,000 samples, and it grows
faster than the number of samples. Without a subset option, `verify` warns before it starts.

Use both kinds of subset, because they check different things:

- `--transcripts 20` with every sample checks ProHap-fast at your full sample count. This matters: the reference
  reads its tables in pandas chunks whose size depends on the number of samples (32 rows at 26,000 samples, 1 row at
  500,000), and some variant IDs come out differently depending on the chunk.
- `--samples 2000` checks every transcript of the chosen chromosomes in reasonable time.

A good first check is `--chromosomes 21,22 --samples 2000 --transcripts 50`, then `--transcripts 20` alone.

## What is compared

Every file that the scripts switched by `implementation` write, in workflow order: the cleaned reference proteome
(and the transcript list, when it is built in the work directory); the filtered VCF, haplotype table, protein FASTA
and removed-samples log (`log/prohap_chr<c>.log`) of each chromosome; the merged table; `haplo_table_file` and its sample-ID table; `haplo_fasta_file`; `haplo_cdna_file` if set; the cleaned and
de-duplicated intermediate FASTAs; `final_fasta_file`; and the simplified FASTA and header table. With ProVar on,
its table and FASTA too. Files are compared after decompression, because gzip headers carry a timestamp. A gzip file
never equals a plain one.

The table shows the first 12 characters of each sha256 (`report.tsv` has them in full), IDENTICAL or DIFFERENT, and
the first differing line; the lines themselves follow the table. MISSING means a file exists in only one arm. Then
come the wall and CPU time of each run.

The exit status is 0 only when every output is identical, 1 when any output differs or the ProHap-fast run fails,
and 2 when `verify` itself cannot run or the reference run fails.

ProHap-fast never guesses at input it does not handle: the script hands such input to the reference script of the
same name, which writes that step's outputs (for example a carriage return or a non-ASCII byte in the VCF, unphased
genotypes reaching `prohap.py`, or sample metadata the reference would stop on). `verify` lists these steps, with the
reason, after the ProHap-fast run. Their outputs are compared as usual, but for them the comparison checks the
reference against itself, and the step runs at the reference's speed.

## Checks before the runs

`verify` reads the header and the first 1,000 records (at most 256 MB) of every VCF chunk and warns about input that ProHap mishandles
or that ProHap-fast hands to the reference code: unphased genotypes, FORMAT fields other than GT, alleles written
with more than one character, sample names that are numbers (pandas reads them as numbers in the metadata, and the
samples are dropped), samples missing from `sample_metadata_file`, chunks with different sample lists, and chunks
out of order or overlapping (checked in full only for chunks read to the end or indexed with tabix). The warnings do
not stop the run.

## Things to know

- **Each implementation needs its own work directory.** Changing `implementation` in a directory that already holds
  outputs reruns nothing: Snakemake sees the same rules with the same code and considers the outputs up to date.
  `verify` gives each arm its own directory for this reason; do the same if you compare by hand.
- **Records with ID `.`** get their ID from their line number in the filtered VCF, header lines included. `verify`
  keeps every line, so the IDs match those of a full run. Tools that drop or add lines (`bcftools view` with `-c`,
  `-i` or a region, and even `bcftools view -s`, which adds missing header definitions) change the IDs; set IDs
  first with `bcftools annotate --set-id` if you need them stable.
- **Compare like with like.** The macOS build of pandas parses some 17-digit numbers differently from the Linux
  build, so output written on a Mac can differ from output written on Linux for large cohorts. `verify` runs both
  arms on the same machine, so its comparison is always valid; just do not compare a Mac run with a Linux run.

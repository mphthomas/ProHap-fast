# ProHap-fast

ProHap-fast builds protein haplotype databases for proteogenomics from phased genotypes, like
[ProHap](https://github.com/ProGenNo/ProHap), of which it is a fork. It is meant for biobank cohorts: hundreds of
thousands of phased whole-genome samples in one run, inside a Trusted Research Environment (TRE) with no network.

- **Same output as ProHap.** ProHap-fast's output is byte-identical to that of its reference implementation, a
  frozen copy of ProHap kept in [`legacy/`](legacy/README.md). You can check this on your own data with one command,
  [`verify`](#check-it-on-your-data-verify). The reference differs from upstream ProHap only where the
  [CHANGELOG](CHANGELOG.md) says so; chiefly, when two haplotypes merge into one protein haplotype, the merged row
  lists the samples of both, where upstream kept only the first list.
- **No downloads at run time.** Each image carries one Ensembl release with every reference ProHap derives from it
  already built: the gffutils databases, the MANE transcript list, the merged cDNA, the tagged reference proteome
  and the tagged contaminants.
- **Docker or Singularity/Apptainer**, as any user.

**Status: 1.0.0.** `src/`, the fast implementation, rewrites every step whose cost grows with the cohort and keeps the
output byte-identical. The VCF filter and `prohap.py`, now one step, grow linearly with the number of samples: on 60
simulated transcripts (a 315th of the exome) they take about 10 CPU seconds at 26,000 samples, where ProHap takes
2,644, and about 230 CPU seconds with at most 4 GB of memory at 500,000 samples, which ProHap cannot run. Every
genome-wide step after them streams its input, in about 3 GB at most at 500,000 samples. Whole-genome VCFs can be
read through their index, exons only (`restrict_to_exons`). Every release is checked against the reference on every
fixture, and CI checks every commit inside the image under Docker and Apptainer; see
[docs/IMPLEMENTATION_NOTES.md](docs/IMPLEMENTATION_NOTES.md) for the plan and the measurements, and
[sizing a large cohort](#sizing-a-large-cohort). ProVar (individual variants) runs unchanged from upstream.

When you use ProHap-fast, please cite the ProHap paper ([below](#citation-and-licence)).

## Contents

- [Quick start](#quick-start)
- [Images and Ensembl releases](#images-and-ensembl-releases)
- [Running with Docker](#running-with-docker) and [with Singularity or Apptainer](#running-with-singularity-or-apptainer)
- [Check it on your data: verify](#check-it-on-your-data-verify)
- [Input requirements](#input-requirements)
- [Configuration](#configuration)
- [Sizing a large cohort](#sizing-a-large-cohort)
- [Building the image](#building-the-image), [running without a container](#running-without-a-container)
- [Output](#output), [development](#development), [citation and licence](#citation-and-licence)

## Quick start

The 1000 Genomes example of this repository (2,548 samples, common variants, 23 chromosomes), with the image of
Ensembl release 113 from Docker Hub (about 9 GB; or [build it](#building-the-image)):

```bash
git clone https://github.com/mphthomas/ProHap-fast.git && cd ProHap-fast
IMAGE=mphthomas/prohap-fast:1.0.0-ensembl.v113
docker pull $IMAGE

tar -xzf sample_data/sample_1kGP_common_global.tar.gz -C sample_data
cp config_example1.yaml config.yaml
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" $IMAGE              # dry run: what would run
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" $IMAGE --cores 4 -p  # the run: about 20 minutes
```

The database is `results/example1_concatenated_database_simplified.fasta`, with its headers in
`results/example1_concatenated_database_header.tsv` and the haplotypes in `results/example1_haplotypes.tsv`.

## Images and Ensembl releases

An image holds one Ensembl release and is tagged with it: `prohap-fast:<version>-ensembl.v<release>`, for example
`prohap-fast:1.0.0-ensembl.v113` is ProHap-fast 1.0.0 with the references of Ensembl release 113 (GRCh38). The
release is also in the image's environment (`ENSEMBL_RELEASE`) and labels (`io.prohap.ensembl_release`, with the
version, commit and build date in the `org.opencontainers.image.*` labels).

- `ensembl_release` in `config.yaml` must equal the image's release; a config for another release stops before
  anything runs. To use another release, use (or [build](#building-the-image)) the image of that release. Outputs of
  different releases differ by design.
- Images are built for `linux/amd64` (Linux servers, HPC, TREs) and `linux/arm64` (Apple Silicon). They are about
  9 GB, of which the references take 8.2 GB.
- Images are published on Docker Hub as `mphthomas/prohap-fast:<version>-ensembl.v<release>`, from 0.1.0: one tag
  for both platforms, and the same tag with `-amd64` or `-arm64` for one. The examples below use the published
  name; an image you [build](#building-the-image) yourself is `prohap-fast:<version>-ensembl.v<release>`.
- For a TRE without registry access, transfer the image as a file. The smallest is a SIF built outside (about
  1.3 GB): `singularity build prohap-fast.sif docker-daemon:IMAGE`, or from a registry with `docker://`. A Docker
  archive also works: `docker save -o prohap-fast.tar IMAGE` (9 GB; gzip it for the transfer, to about 1.6 GB),
  then, uncompressed, `singularity build prohap-fast.sif docker-archive:prohap-fast.tar` inside. Apptainer 1.5
  does not read a gzipped archive.

Inside the image, `/opt/prohap` holds the workflow (`Snakefile`, `src/`, `legacy/src`), `verify/`, the references
in `refs/ensembl.v<release>/`, `crap.fasta` (the default contaminants) and `config_tre_template.yaml`. Everything
is read in place; nothing is copied into your directories.

## Running with Docker

Run from the work directory, the directory that holds `config.yaml`. It also receives the intermediates
(`data/`, `results/`, `tmp/`), the logs and Snakemake's state, and is mounted at its own path:

```bash
cd /path/to/workdir                       # contains config.yaml
docker run --rm -u "$(id -u):$(id -g)" \
    -v "$PWD":"$PWD" -w "$PWD" \
    -v /data/phased_vcfs:/data/phased_vcfs:ro \
    -v /results/prohap:/results/prohap \
    mphthomas/prohap-fast:1.0.0-ensembl.v113 --cores 16 -p
```

- Mount every path `config.yaml` refers to (VCFs, sample metadata, outputs) at the same path inside the container,
  so that the config resolves unchanged. Relative paths in the config are relative to the work directory.
- `-u "$(id -u):$(id -g)"` runs the container as you, so every file it writes is yours.
- After the image name comes a command:

| Command | What it does |
| --- | --- |
| (none) | a dry run: lists the jobs the workflow would run |
| `[snakemake] ARGS` | runs the workflow, `snakemake ARGS`, e.g. `--cores 16 -p` (`--rerun-triggers mtime` is added) |
| `verify [ARGS]` | runs the reference and ProHap-fast on the same input and compares every output ([below](#check-it-on-your-data-verify)) |
| `shell [ARGS]` | `bash` inside the image |
| `help` | the list of commands |

A run that stops can be resumed with the same command: Snakemake reruns only what is missing. Add `--keep-going`
to finish the other chromosomes when one fails.

## Running with Singularity or Apptainer

The same image runs under Singularity and Apptainer, which mount the current directory and run as you without
being asked:

```bash
singularity build prohap-fast.sif docker://mphthomas/prohap-fast:1.0.0-ensembl.v113   # or docker-archive:FILE.tar
cd /path/to/workdir                       # contains config.yaml
singularity run --bind /data/phased_vcfs,/results/prohap prohap-fast.sif --cores 16 -p
singularity run --bind /data/phased_vcfs prohap-fast.sif verify --cores 8 --chromosomes 21,22 --transcripts 20
```

Bind every path `config.yaml` refers to that is not under the current directory or your home directory. The
commands are those of the Docker table above. The image's references and environment are fixed: the host's
`PROHAP_REFS`, `PYTHONHASHSEED` and user site-packages (`~/.local`) do not reach the workflow. Snakemake's cache
goes to `.cache/` in the work directory, and its home to `.home/` there when the home directory is not writable
(Docker with `-u`, Singularity with `--no-home`).

## Check it on your data: verify

`verify` runs the workflow twice on your inputs, once with the reference (`implementation: legacy`) and once with
ProHap-fast (`implementation: fast`), each in its own directory under `verify_runs/<timestamp>/`, and compares
every output byte for byte after decompression. Your `config.yaml` and your outputs are not touched. It prints a
table (file, both sha256 sums, IDENTICAL or DIFFERENT, the first differing line) and both run times, and exits with
status 0 only when every output is identical.

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" -v /data:/data:ro \
    mphthomas/prohap-fast:1.0.0-ensembl.v113 verify --cores 8 --chromosomes 21,22 --samples 2000 --transcripts 50
```

The reference is slow at cohort scale (about 2,000 CPU seconds per 60 transcripts at 26,000 samples), so restrict
both runs: `--chromosomes 21,22`, `--samples N` (a seeded random draw) or `--sample-list FILE`, and
`--transcripts N` (the first N transcripts of each chromosome). `--transcripts 20` alone checks ProHap-fast at your
full sample count, which matters because the reference reads its tables in chunks whose size depends on the number
of samples. Details, and what each check covers: [verify/README.md](verify/README.md).

## Input requirements

- **Phased, multi-sample VCF**, one chromosome per file or several chunks per chromosome (a glob in
  `phased_vcf_file_name`). Text VCF, plain or gzip/bgzip; convert BCF with `bcftools view -Ov`, or index it and set
  `restrict_to_exons` (below). Records sorted by position.
- **Chunks** are read in alphabetical order and only the first header is kept, so their names must sort in genomic
  order, they must not overlap, and every chunk must list the same samples in the same order.
- **Genotypes**: `FORMAT` is `GT` alone, phased, one character per allele (`0|1`, `.|0`, haploid `1` for chrX in
  males). This is what the fast implementation is built for; input outside it (unphased or multi-character
  genotypes, other FORMAT fields, quotes, CR or non-ASCII bytes, duplicate sample or transcript IDs) is processed by
  the reference code, at the reference's speed, never guessed at.
- **Sample IDs must not be purely numeric**: pandas reads them as numbers from the metadata, and every sample is
  dropped without an error.
- **`sample_metadata_file`** lists every sample of the VCF, tab-separated, with the columns `Sample name`, `Sex`
  (`male` for males), `Population code` and `Superpopulation code`. Samples missing from it are dropped from the
  analysis and from every frequency.
- **Records with ID `.`** are named after their line number in the filtered VCF. Anything that adds or removes lines
  before them (another region, other filters, a different header) renames them; set IDs first with
  `bcftools annotate --set-id` if you need stable names.

`verify` checks the header and the first records of every chunk for most of these and warns before it runs.

## Configuration

Start from [`config_tre_template.yaml`](config_tre_template.yaml) (also in the image as
`/opt/prohap/config_tre_template.yaml`), save it as `config.yaml` in the work directory and edit the lines marked
`<<< EDIT >>>`. [`config_file_example`](config_file_example) documents every setting. The ones that matter most:

- `ensembl_release` must be the image's release. Nothing else names the release: the annotation and the reference
  transcript list follow from it.
- `implementation: fast` (default) or `legacy` picks the implementation; both write the same output.
- `phased_min_af` sets which variants are kept, and so the size of almost everything downstream.
- `haplo_min_count` applies only when `haplo_min_freq` is `-1`; any other `haplo_min_freq` thresholds by frequency
  instead. `haplo_min_freq: 0` keeps every private haplotype, which at cohort scale gives a database too large to
  search.
- `final_fasta_file` must not end in `.gz`; the other output names may.
- `restrict_to_exons: True` (default `False`) reads only the VCF records that overlap an exon of a listed
  transcript, through the index, before the filter; the rest of the workflow sees the same records in the same
  order. Worth it for whole-genome VCFs, where about 1 record in 50 is exonic. Every chunk must be indexed
  (bgzip + `tabix -p vcf`, or BCF + `bcftools index`; BCF is read only this way). Records with ID `.` get other
  names, since they are named by their line in the filtered VCF: set IDs first (`bcftools annotate --set-id`) if
  results must match a run without the option. `verify` applies the option to both arms alike.

Settings of upstream ProHap that downloaded files (`Ensembl_FTP_URL`, `annotationFilename`, `phased_FTP_URL`) are
ignored with a warning while they agree with `ensembl_release`, and stop the run when they do not.

## Sizing a large cohort

Snakemake runs one `compute_haplotypes` job per chromosome: the filter and `prohap.py`, which reads the VCF once and
hands each transcript to `max_cores` worker processes (default 4; more do not help, because the reader is the
bottleneck). The steps after it run once, on the whole genome. So the run takes about as long as its slowest
chromosome, chromosome 1 or 2, plus the genome-wide steps, provided every chromosome runs at once.

Each `compute_haplotypes` job tells Snakemake it needs `mem_mb_per_chromosome` (default 10000, 10 GB), so
`--resources mem_mb=N` caps the chromosomes in flight to what fits in N MB. Measured on simulated data (60
transcripts of chromosome 22, a 315th of the exome; IMPLEMENTATION_NOTES sections 3.8 and 8.2): at 500,000 samples
the job takes about 230 CPU seconds per block and peaks at about 4 GB with four workers (6.4 GB for TTN, the longest
transcript, alone); at 26,000 samples about 10 CPU seconds and under 1 GB.

For 500,000 samples, every chromosome at once:

```bash
docker run ... mphthomas/prohap-fast:1.0.0-ensembl.v113 --cores 96 --resources mem_mb=250000 -p
```

| Resource | 500,000 samples | Notes |
| --- | --- | --- |
| Cores | 96 | 23 chromosome jobs × 4 threads; each job also runs its reader and the filter |
| Memory | about 250 GB | 23 × 10 GB as scheduled (about 4–6.5 GB used); the genome-wide steps need about 3 GB |
| CPU | about 20–30 CPU-hours | the benchmark block × 315, on an Apple Silicon Mac; slower cores take longer (a TRE ran 4.4× slower) |
| Wall | chromosome 1 plus the genome-wide steps | about 1.5–2.5 h for chromosome 1 on a fast machine |
| Scratch | a few times the outputs | the per-chromosome tables and FASTAs (gzipped) until the merged ones are written; `merge_duplicate_seq.py` sorts through temporary files the size of its input |

With fewer cores or less memory, lower `--cores` or `--resources mem_mb`: Snakemake then runs fewer chromosomes at
a time and the run takes longer, not more memory. These figures are scaled from the benchmark block, not measured on
a whole genome: run one chromosome first (`verify --chromosomes 22 --transcripts 20` with every sample, then a run
with `chromosomes: ["22"]`) and measure before a genome-wide run. For whole-genome VCFs, `restrict_to_exons: True`
saves reading the records outside exons (see [Configuration](#configuration)).

The steps after the chromosome jobs run on the whole genome, one at a time. Merging the per-chromosome tables
(`merge_tables.py`, `haplo_extract_sample_names.py`) and the protein FASTAs (`merge_duplicate_seq.py`,
`remove_UTR_seq.py`) each stay at about 1 GB or less, and the last, `fasta_simplify_headers.py`, at about 2–3 GB
(0.3 GB plus about 80 MB per million haplotype IDs).

The reference code (`implementation: legacy`) is far slower and larger: about 2,600 CPU seconds and 12 GB per block
at 26,000 samples, and out of reach at 500,000.

## Building the image

On any machine with Docker and network access:

```bash
container/build_image.sh                                   # release 113, linux/amd64 and linux/arm64
container/build_image.sh --release 114 --platform linux/amd64 --no-cache
```

The Ensembl files (GTF, cDNA, ncRNA, protein FASTA) are downloaded while the image is built, and
`container/prebuild_refs.py` builds the references from them with the reference scripts, natively on the build
machine (about 4 minutes); the image of the other platform reuses them. A platform other than the build machine's
needs emulation: Docker Desktop on Apple Silicon builds `linux/amd64` with Rosetta. The script tags
`prohap-fast:<version>-ensembl.v<release>-<arch>` for each platform it builds, and
`prohap-fast:<version>-ensembl.v<release>` for the build machine's own. `<version>` comes from the `VERSION` file.
A first build takes about 10 minutes; after a change to the workflow's code, a few seconds.

## Running without a container

The workflow runs in a conda environment on Linux or macOS, given references built once with
`container/prebuild_refs.py`:

```bash
micromamba create -n ProHap -f envs/prohap.linux-64.lock          # or .linux-aarch64 / .osx-arm64
micromamba run -n ProHap python -m pip install --no-deps --only-binary=:all: --require-hashes -r envs/prohap.pip.lock
# the four Ensembl downloads of the release, as the Dockerfile fetches them, in ~/ensembl-113/
micromamba run -n ProHap python container/prebuild_refs.py --release 113 --downloads ~/ensembl-113 --out ~/prohap-refs
cd /path/to/workdir && PROHAP_REFS=~/prohap-refs micromamba run -n ProHap \
    snakemake --snakefile /path/to/ProHap-fast/Snakefile --cores 8 -p
```

`verify` runs the same way: `python /path/to/ProHap-fast/verify/verify.py --cores 8 ...`. On `linux-aarch64`, also
set `LD_LIBRARY_PATH` to the environment's `lib/` (the image does; see IMPLEMENTATION_NOTES, Pinned versions). On
macOS the workflow starts ProHap's worker processes with `fork`, as on Linux. The macOS build of pandas parses some
17-digit numbers differently from the Linux builds, so at cohort scale a Mac writes some `frequency` values
differently: compare Mac output only with Mac output. `verify` runs both implementations on the same machine, so its
comparison holds everywhere.

## Output

The workflow writes four kinds of output; the [ProHap wiki](https://github.com/ProGenNo/ProHap/wiki/Output-files)
describes them in detail.

1. *Concatenated FASTA* (`final_fasta_file`): the protein haplotypes (and ProVar variants, if enabled), the Ensembl
   reference proteome and the contaminants, ready for any search engine. With `fasta_simplify_headers`, the headers
   are moved into a separate table and a gene name is added to each protein.
2. *Haplotype table* (`haplo_table_file`): genomic coordinates, consequences, frequencies and counts of each
   haplotype, with the sample IDs in a table next to it (`<name>_sampleIDs.tsv`).
3. *Haplotype translations* (`haplo_fasta_file`): the translated haplotype sequences before they are merged with the
   reference proteins and contaminants. When variation in UTRs is ignored (the default), UTR sequences are left out.
4. *Haplotype cDNA* (`haplo_cdna_file`, optional): the cDNA sequences before translation.

To map peptide-spectrum matches back to haplotypes, variants and their genetic origin, use
[ProHap PeptideAnnotator](https://github.com/ProGenNo/ProHap_PeptideAnnotator).

Ready-made ProHap databases from the 1000 Genomes Project, the Haplotype Reference Consortium and the Human
Pangenome Reference Consortium are on Zenodo:
[10.5281/zenodo.10149277](https://doi.org/10.5281/zenodo.10149277),
[10.5281/zenodo.12671301](https://doi.org/10.5281/zenodo.12671301),
[10.5281/zenodo.12686818](https://doi.org/10.5281/zenodo.12686818).

## Development

- `src/` is the fast implementation, `legacy/src` the frozen reference ([never edited](legacy/README.md)); the
  config's `implementation` picks one, and every rule, path and reference is shared.
- `tests/legacy_vs_fast.sh` is the check every change to `src/` must pass. It runs the workflow with the reference
  code and then with `src/` on test data and compares every output: by default the 1kG fixture
  ([tests/fixtures/1kg](tests/fixtures/1kg/README.md)) and a small simulated one, in about three minutes;
  `tests/legacy_vs_fast.sh medium` and `large` add larger simulated cohorts, chromosome X with haploid or mirrored
  males, a table of more than 32,768 rows, 400,000 proteins for `merge_duplicate_seq.py`, and TTN
  ([tests/fixtures/sim](tests/fixtures/sim/README.md)). It runs in the ProHap environment with pre-built references
  (`PROHAP_REFS`), without a container.
- `bench/simulate.py` simulates a phased cohort of any size over a block of MANE transcripts, and `bench/bench.py`
  times the filter and `prohap.py` of each implementation on it (CPU, memory, phases).
- `pytest tests/` runs the unit tests in under a minute (pytest: `pip install --no-deps --require-hashes -r
  envs/test.pip.lock`). With pre-built references, `tests/test_identity.py` also runs `verify` on the fixture.
- CI (GitHub Actions): `unit` runs pytest on every push and pull request; `identity` builds the image from the
  commit and runs `verify` and `tests/legacy_vs_fast.sh` inside it, with the reference's outputs checked against
  golden checksums (`tests/ci_identity.sh`): under Docker on the small tier for every push and pull request, under
  Docker and Apptainer on the small and medium tiers nightly.
- Design, measurements and the rebuild recipes for each speed-up:
  [docs/IMPLEMENTATION_NOTES.md](docs/IMPLEMENTATION_NOTES.md).

Bug reports and contributions are welcome as [issues](https://github.com/mphthomas/ProHap-fast/issues) and pull
requests. A problem you can reproduce with `implementation: legacy` is a ProHap behaviour; report it upstream at
[ProGenNo/ProHap](https://github.com/ProGenNo/ProHap/issues) too.

## Citation and licence

When using ProHap-fast, or databases generated with it, please cite the ProHap paper: Vašíček, J., Kuznetsova,
K.G., Skiadopoulou, D. et al. ProHap enables human proteomic database generation accounting for population
diversity. _Nature Methods_ (2024). [https://doi.org/10.1038/s41592-024-02506-0](https://doi.org/10.1038/s41592-024-02506-0)

MIT licence, as ProHap ([LICENSE](LICENSE)). We follow the
[Contributor Covenant](https://www.contributor-covenant.org/) [Code of Conduct](docs/CODE_OF_CONDUCT.md).

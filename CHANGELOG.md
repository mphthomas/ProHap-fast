# Changelog

All notable changes to ProHap-fast are documented in this file. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The **Output changes** sections together list every way the output differs from upstream
[ProGenNo/ProHap](https://github.com/ProGenNo/ProHap) at commit `cafcc51` for the same inputs and settings. Each
change is also made in `legacy/src`, the reference implementation that ProHap-fast is checked against.

## [Unreleased]

### Added

- **`restrict_to_exons`** (opt-in, default `False`): read only the VCF records that overlap an exon of a listed
  transcript, through the index of each chunk (bgzip + tabix, or BCF + `bcftools index`), before the filter. The
  workflow then sees the same records in the same order, so the output is the same except for the names of records
  whose ID is `.`, which are numbered by their line in the filtered VCF. For whole-genome VCFs, where most records
  are outside exons. Indexed BCF can be read this way.

## [0.3.0] - 2026-11-19

### Changed

- **The per-chromosome haplotype tables are merged without holding them in memory; the output is the same.**
  `merge_tables.py` and `haplo_extract_sample_names.py` read the whole genome's haplotype table into pandas and wrote
  it back, which at 500,000 samples would need well over 100 GB. They now read it in pandas' own chunks of 32,768
  rows and let pandas re-type each chunk, so every value pandas rewrites (for example a frequency
  `0.17837519623233908` written back as `0.178375196232339`) comes out as before; where pandas would join two chunks
  by rewriting values (a column of whole numbers in one chunk and decimals in another), the step is handed to the
  reference code. On simulated tables of 500,000 samples (1.1 GB of text) each step takes 10–16 CPU seconds instead
  of 161–173, and stays at about 0.4 GB of memory on Linux however large the table.

## [0.2.0] - 2026-11-19

### Changed

- **`prohap.py` grows linearly with the cohort and runs at 500,000 samples; the output is the same.** The reference
  wrote, and read back, a table with one column per sample for every transcript, so its cost grew faster than the
  cohort: at 100,000 samples a chromosome with three short transcripts already needed 5–7 GB, and 500,000 samples
  were out of reach. `prohap.py` now reads the VCF once and keeps, for each
  variant that falls in an exon, only the haplotype copies that carry it; it groups the copies into haplotypes by
  hashing, checks every group exactly, and computes each haplotype's protein with the reference's own steps, in
  parallel workers, while it reads. The values the reference took from pandas, such as IDs that pandas retypes
  depending on the size of the cohort (`2e5` becoming `200000.0`), are reproduced. On the benchmark block (60
  simulated transcripts of chromosome 22, one thread), the filter and `prohap.py` together take 10.5 CPU-s at
  26,000 samples, where 0.1.0 took 499 and the reference 2,644; at 500,000 samples they take 296 CPU-s, and
  `prohap.py` peaks at 4.0 GB with four workers. TTN, the longest transcript, takes 22 CPU-s at 26,000 samples and
  517 at 500,000. Genome-wide, that is about 26 CPU-hours for 500,000 samples on an Apple Silicon Mac.
- **`vcf_filter_fix.py` is about 27 times faster.** It handles the genotypes as bytes with numpy instead of one
  Python string per sample, and writes its output at gzip level 1; 2.8 CPU-s instead of 76.5 for the block at 26,000
  samples.
- **Input the fast code does not handle is handed to the reference code, never guessed at.** A step that meets it
  (for example a carriage return or a non-ASCII byte in the VCF, unphased or multi-character genotypes reaching
  `prohap.py`, a FORMAT other than `GT` there, or sample metadata on which the reference stops) runs the reference
  script of the same name with the same arguments, which writes that step's outputs. Lines of the filter outside its
  fast path go through the reference's code for that line. `verify` lists the steps that were handed over.

## [0.1.0] - 2026-11-19

### Added

- **`verify` checks ProHap-fast against the reference on your own data.** Run from the directory that holds
  `config.yaml`, `verify/verify.py` runs the workflow twice, with `implementation: legacy` and with
  `implementation: fast`, each in its own work directory under `verify_runs/`, and compares every output byte for
  byte after decompression. It can restrict both runs to some chromosomes, samples or transcripts, checks the VCFs
  before it starts, and exits with status 0 only when every output is identical. See `verify/README.md`.
- **One image per Ensembl release, with every reference pre-built.** `container/build_image.sh` builds
  `prohap-fast:<version>-ensembl.v<release>` for `linux/amd64` and `linux/arm64`: the workflow, both
  implementations, `verify`, the locked environment and `ensembl.v<release>/` from `container/prebuild_refs.py`
  under `/opt/prohap/refs`. The Ensembl files are downloaded while the image is built, never when it runs. The
  version is in the new `VERSION` file and, with the release, the commit and the build date, in the image's labels.
  The image is about 9 GB.
- **The image runs the workflow, `verify` or a shell, in the current directory.** `IMAGE` alone shows a dry run,
  `IMAGE [snakemake] ARGS` runs the workflow, `IMAGE verify ARGS` runs `verify`, `IMAGE shell` opens bash and
  `IMAGE help` lists these. A config whose `ensembl_release` is not the image's release stops before Snakemake
  starts, naming both releases. The image's references, hash seed and packages are used whatever the host
  environment holds (`PROHAP_REFS`, `PYTHONHASHSEED`, a `~/.local` site-packages under Singularity).

### Changed

- **`prohap.py` is about 6× faster at 10,000–26,000 samples, uses less than 40% of the memory at 26,000, and writes
  the same output.** Four small changes to the reference code, each checked byte for byte against it on every test
  fixture: the variant columns of a transcript's table are read once as lists, where the reference built a row of
  every sample per haplotype and variant (the largest gain); sample sex and populations are looked up in dicts;
  haplotypes are grouped by transcript once; and each transcript's VCF lines are collected in a list and joined
  once, where the reference grew one string per transcript (the memory saving). On the benchmark block (60
  transcripts of chromosome 22, simulated, one thread) `prohap.py` takes 94.1 CPU-s instead of 530.6 at 10,000
  samples, and 421.9 instead of 2,567.2 at 26,000, with a peak of 4.7 GB instead of 12.3. It still grows
  faster than the cohort, because its per-transcript tables have one column per sample; the rewrite for 100,000
  samples and more comes in 0.2.0.
- **`merge_duplicate_seq.py` no longer slows down quadratically with the number of sequences.** It kept the
  sequences in a list sorted by hash and inserted each new one; it now collects them in a dict and sorts once, in
  the same order. 400,000 protein records take 9.0 CPU-s instead of 18.7, and the gain grows with the
  number of sequences: the design work estimated 4 minutes instead of 8.5 hours for 25 million.
- **Ensembl references are pre-built and read in place; nothing is downloaded.** Everything the workflow derives
  from the Ensembl release (whole-genome and per-chromosome gffutils DBs, the MANE transcript list, the merged cDNA,
  the tagged reference proteome and the tagged contaminants) is read from `<refs_dir>/ensembl.v<ensembl_release>/`,
  built once by `container/prebuild_refs.py` and shipped in the release image. `refs_dir` defaults to
  `$PROHAP_REFS`, then `/opt/prohap/refs`. The rules that downloaded or built these files in the work directory are
  gone. The output does not change.
- **`ensembl_release` is the only release setting.** `annotationFilename` follows from it. `Ensembl_FTP_URL`,
  `annotationFilename` and `phased_FTP_URL` are ignored with a warning when they agree with `ensembl_release`, and
  stop the run when they name another release or annotation. A missing reference directory, a `release.txt` that
  disagrees, or a missing chromosome DB also stops the run before any job starts, naming the config key.
- **`custom_transcript_list` can be left out.** It then means the reference list, `transcripts_reference_<release>.csv`
  (pre-built when `only_MANE_select` is set); the upstream value `data/transcripts_reference_<release>.csv` means the
  same list. An empty value still means every transcript of `included_transcript_biotypes`. A `contaminants_fasta`
  other than ProHap's `crap.fasta` is still tagged in the work directory.

### Output changes

- **Merged haplotypes list all their samples.**
  - ProHap merges haplotypes of a transcript into one row when the variants it keeps have the same VCF IDs. It
    drops variants whose ALT allele matches the reference cDNA, which includes many indels, and, when
    `haplo_ignore_UTR` is 1, variants outside the coding sequence. Different variants whose IDs read the same also
    merge: the alleles of a multi-allelic record with an rsID (the filter gives every allele the record's ID), or
    IDs that pandas reads as the same value, such as `NA` and `nan`. A merged row can therefore list samples whose
    haplotype also carried a dropped variant or a different allele, which its `DNA_changes` do not show.
  - The merged row keeps the first haplotype's ID and variant columns and adds up `occurrence_count`, `frequency`,
    `frequency_population` and `frequency_superpopulation`. Upstream kept only the first haplotype's `samples`.
    The row now lists the samples of every merged haplotype, separated by `;`: upstream's list first, then the
    others in the order ProHap processed them.
  - The number of entries in `samples` (one `SAMPLE:1` or `SAMPLE:2` per haplotype copy) now always equals
    `occurrence_count`. In upstream output, the affected rows are those where `occurrence_count` is larger.
  - Haplotypes with an unknown reading frame are not merged: as upstream, a later haplotype with the same IDs
    replaces the earlier one. These rows do not change.
  - Only `samples` changes: in the sample-ID table written next to `haplo_table_file`
    (`<name>_sampleIDs.tsv[.gz]` for `<name>.tsv[.gz]`), and in the intermediate tables `haplo_chr<chr>.tsv.gz`
    and `haplo_all.tsv.gz`. `haplo_table_file` has no `samples` column. Nothing else in the output changes.
- **A chromosome with no variants after filtering gets valid gzip files.** When a chromosome's filtered VCF has no
  records (for example, no variant passes `phased_min_af`), ProHap writes an empty haplotype table and protein
  FASTA for it. Upstream wrote them as plain text even when the file name ends in `.gz`, which the per-chromosome
  table's name always does, so the merge step stopped with `BadGzipFile`. They are now gzip-compressed when the
  name ends in `.gz`: the table holds only the header line and the FASTA is empty. Names without `.gz` are written
  as before. As upstream, no cDNA FASTA is written for such a chromosome.

[Unreleased]: https://github.com/mphthomas/ProHap-fast/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/mphthomas/ProHap-fast/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/mphthomas/ProHap-fast/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/mphthomas/ProHap-fast/compare/cafcc51...v0.1.0

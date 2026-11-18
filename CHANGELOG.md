# Changelog

All notable changes to ProHap-fast are documented in this file. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The **Output changes** sections together list every way the output differs from upstream
[ProGenNo/ProHap](https://github.com/ProGenNo/ProHap) at commit `cafcc51` for the same inputs and settings. Each
change is also made in `legacy/src`, the reference implementation that ProHap-fast is checked against.

## [Unreleased]

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

### Changed

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

[Unreleased]: https://github.com/mphthomas/ProHap-fast/compare/cafcc51...HEAD

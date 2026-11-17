# Changelog

All notable changes to ProHap-fast are documented in this file. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The **Output changes** sections together list every way the output differs from upstream
[ProGenNo/ProHap](https://github.com/ProGenNo/ProHap) at commit `cafcc51` for the same inputs and settings. Each
change is also made in `legacy/src`, the reference implementation that ProHap-fast is checked against.

## [Unreleased]

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

[Unreleased]: https://github.com/mphthomas/ProHap-fast/compare/cafcc51...HEAD

# legacy/: the frozen reference

`legacy/src` is the reference implementation that ProHap-fast's output is checked against: `verify` and the CI
identity tests run both and require byte-identical output. It is upstream
[ProGenNo/ProHap](https://github.com/ProGenNo/ProHap) `src/` at commit `cafcc51`
(`cafcc515ba3ba6b1b6121db6ea3ac612ca053124`) with `legacy-fixes.patch` applied, and nothing else.

**Do not edit anything in this directory.** No optimisation, refactoring or reformatting; speed work goes in
`src/`. `tests/test_legacy_frozen.py` rebuilds `legacy/src` from git history and fails if any file differs by a
single byte, or is added, removed, made executable or turned into a symlink. It ignores untracked `__pycache__/`
directories, which running the legacy arm creates, but fails if git tracks one.

## What `legacy-fixes.patch` changes

| Change | Files | Effect on output |
| --- | --- | --- |
| Merged haplotypes list the samples of every merged haplotype | `modules/process_haplotypes.py` | `samples` of merged rows ([CHANGELOG](../CHANGELOG.md)) |
| `empty_output` writes gzip when the file name ends in `.gz` | `modules/process_haplotypes.py` | files of a chromosome with no variants ([CHANGELOG](../CHANGELOG.md)) |
| Sets instead of lists for two membership tests | `prohap.py`, `modules/process_haplotypes.py` | none (speed only) |

## How it was built

```bash
mkdir legacy
git archive cafcc51 src | tar -x -C legacy/
git apply --directory=legacy legacy/legacy-fixes.patch
```

## Changing the reference

Only a deliberate output change touches `legacy/src`. It is made to `legacy/src` and `src/` together in its own
commit, with a CHANGELOG entry and updated test fixtures, and the same commit regenerates the patch so the test keeps
passing:

```bash
python tests/test_legacy_frozen.py --write-patch
```

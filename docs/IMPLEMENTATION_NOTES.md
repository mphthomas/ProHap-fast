# ProHap-fast implementation notes

Rebuild recipes for every planned change, recovered from the design work (Oct 2026). The prototypes were verified
byte-identical against the reference on all fixtures before they were lost; these notes record what they did and
why, including the traps that broke identity along the way.

Reference = `legacy/src` = upstream `cafcc51` + `legacy-fixes.patch`.

---

## Pinned versions

The reference's output, and the emulation in §3.3 and §4, depend on pandas internals: the row-chunk heuristic,
per-chunk typing, `concat_compat` and the float parser. The runtime is therefore pinned to what the TRE image
`mphthomas/prohap:v0.4-ens113` ran: Python 3.10.20, **pandas 2.3.3**, numpy 2.2.6, pyarrow 25.0.1, biopython 1.88,
gffutils 0.14, pyvcf3 1.0.4, pyyaml 6.0.3 and snakemake-minimal 7.32.4, plus bcftools and htslib 1.25 for `verify` and
region streaming.

- `envs/prohap.yaml` lists the pins. The conda locks `envs/prohap.{linux-64,linux-aarch64,osx-arm64}.lock` and
  `envs/prohap.pip.lock` fix every package build; the images install from them, and `envs/update_locks.sh`
  regenerates them after an edit to `prohap.yaml`.
- On linux-64 the lock has the TRE image's exact builds of python, numpy, biopython, gffutils, pyvcf3, pyyaml,
  snakemake-minimal, pyfaidx and simplejson, and the pandas and pyarrow wheels are the TRE image's own. Other
  packages are newer: system libraries (openssl, libsqlite, libgcc, icu, ...) and some of Snakemake's pure-Python
  dependencies. bcftools, htslib and their dependencies are additions. The ProHap scripts import none of the newer
  Python packages.
- pandas and pyarrow are pip wheels, as in the TRE image. Their C++ extensions must load the environment's
  libstdc++, so the image puts the environment's `lib/` on `LD_LIBRARY_PATH`; without it, importing gffutils after
  pandas fails on linux-aarch64, because the environment's ICU needs the newer C++ ABI. (On linux-64 numpy loads
  the environment's libstdc++ first; the same failure needs pyarrow imported before sqlite3, which no script does.)
- Checked 2026-11-17 on 1kG chr22, threshold 0, every output from the filtered VCF to the simplified FASTA,
  decompressed: the linux/amd64 and linux/arm64 images built from the locks (both arms), the linux/amd64 image from
  before the pins, the TRE image itself, and native macOS runs with and without the osx-arm64 lock are
  byte-identical. The gffutils DBs, MANE list, merged cDNA and tagged proteome were shared pre-built inputs; a
  separate check rebuilt the chr22 DB and the MANE list in the TRE image and both locked images, with identical
  contents (the DB files differ only in the SQLite version stamped in their header).
- **macOS is not Linux for floats.** The macOS pandas wheel parses some 17-digit floats differently from the Linux
  wheels: of 50,000 realistic rare-haplotype frequencies, 2,631 differed (`1.8032476490158774e-05` stays as is on
  the Mac and becomes `1.803247649015877e-05` on Linux and in the TRE image). The legacy merge therefore writes
  different `frequency` text on a Mac once frequencies are that small, i.e. in cohorts of tens of thousands; the
  1kG fixture is too small to show it. Compare Mac output only with Mac output (fast vs legacy); produce golden
  checksums and TRE comparisons in the Linux images.
- Changing any version, pandas above all, needs the full identity suite to pass before the commit lands.

---

## Pre-built references

`container/prebuild_refs.py` turns the four Ensembl downloads of a release (GTF `chr_patch_hapl_scaff`, cDNA, ncRNA,
pep) and `crap.fasta` into `<refs>/ensembl.v<rel>/`: the whole-genome gffutils DB, one DB per chromosome (1–22, X,
Y, MT), the MANE list `transcripts_reference_<rel>.csv`, `total_cdnas_<rel>.fa`, the tagged proteome and the tagged
contaminants, plus `release.txt`, `SOURCES.sha256` (the inputs) and `MANIFEST.sha256` (`sha256sum -c`). It runs the
`legacy/src` scripts with the arguments of the upstream rules it replaces (`parse_gtf_whole`, `split_gtf` +
`parse_gtf_chromosome`, `merge_cdnas_fasta`, `reference_filter_format`, `default_transcript_list` with -MANE 1,
`contaminants_fix_headers`), in parallel, and renames `ensembl.v<rel>.partial/` into place only when every step
has passed.

Measured 2026-11-18 on the Mac, v113, 8 jobs: 156 s wall (whole-genome DB 141 s, then the MANE list 11.8 s; the 25
chromosome DBs, 0.1–12 s each, run beside it), 0.37 GB peak RSS, 8.17 GB on disk. Every file is byte-identical to
what the M1 workflow built in its work directory from the same downloads with the same environment, the DBs
included; DBs built with another SQLite differ only in the version field of the header (offset 96).

### In the release image

The Dockerfile downloads the four Ensembl files and runs `prebuild_refs.py` in stages pinned to the build platform
(`FROM --platform=$BUILDPLATFORM`), so the references are built natively once and the image of the other platform
copies the same layer from the build cache. The prebuild stage's environment is the runtime environment built from
the same locks by the same command (`container/install_env.sh` picks the lock from `uname -m`), so the two are one
set of layers when build and target platform agree. Trap: a stage `FROM <another stage>` takes its files from that
stage but its platform from its own FROM line, which defaults to the target platform; without an explicit
`--platform=$BUILDPLATFORM` there, each target platform reran the prebuild as its own cache entry.

Measured 2026-11-18 on the Mac (Docker Desktop, 16 CPUs, 7.7 GB VM; `linux/amd64` through Rosetta): the first
build of `linux/arm64`, with nothing cached but the base images, took 522 s, of which the prebuild took 336 s
(8 jobs; the whole-genome DB 182 s). Image size 9.00 GB (`linux/arm64`) and 9.05 GB (`linux/amd64`), of which the
references are 8.17 GB; the M1 image `prohap-fast:m1.8` was 2.31 / 2.47 GB. The references built in the image (Linux, aarch64) are
byte-identical to those built natively on the Mac (`MANIFEST.sha256`), except `crap_tagged.fa.gz`, whose gzip header
carries the time it was written; decompressed it is identical too.

---

## Checking `src/` against `legacy/src`

`tests/legacy_vs_fast.sh` is the check every change to `src/` must pass. It runs on the host, in the ProHap
environment, with no container. For each fixture config it runs the workflow twice on the same inputs, once
with `--config implementation=legacy` and once with `--config implementation=fast`, each in its own work directory
(switching `implementation` in place reruns nothing). Then it compares every file under the two directories, except
`.snakemake/`, after decompression. A file that is gzip in one run and plain in the other is DIFFERENT.

The fixtures come in tiers. `small` (the default) is the two 1kG configs below and two configs of a small simulated
dataset; `medium` and `large` are simulated cohorts of 2,500 to 100,000 samples, a whole chromosome, TTN and
400,000 proteins for `merge_duplicate_seq.py` (`tests/fixtures/sim/README.md`, data from section 5's simulator,
generated on first use into a git-ignored cache). The tiers' run times are in the script's header.

The fixture keeps only 109 transcripts (`tests/fixtures/1kg/make_transcripts.py`), because the reference costs
about 0.1 s of CPU per transcript at 2,548 samples whether the transcript has variants or not. For every transcript,
`parse_vcf` writes a table with one column per sample, and `get_gene_haplotypes` reads it back. Measured
2026-11-18 on the Mac, legacy arm, threshold 0, 4 threads, 4 cores:

| Data | Wall time per run |
| --- | --- |
| chromosomes 1, 21, 22, X, all 3,647 MANE transcripts | 362 s (chr1 194 s, X 79 s with its 7 records) |
| chromosome 22 alone, 495 transcripts, 1 thread | 75 s (prohap 68 s) |
| 109 transcripts (a draft of the fixture's list), 4 cores | 42 s; the 12 s of cleaning the reference proteome held back the 4-thread prohap jobs |

With `-j 8` both fixture configs take 1 min 35 s for the four runs: 19–20 s per threshold0 run and 24–25 s per
template run. The floor is the steps every run repeats: `reference_remove_stop` (12 s), then mixing with the
reference proteome, `merge_duplicate_seq` and `fasta_simplify_headers` (about 5 s).

---

## CI and the golden checksums

Two GitHub Actions workflows, on GitHub-hosted `ubuntu-24.04` runners (free for this public repository):

- `unit.yml`, every push to `main` and every pull request: `pytest tests/` in the locked environment, created by
  `container/install_env.sh` with micromamba 1.5.10 (the image's), plus `envs/test.pip.lock` for pytest. No image,
  no references; the tests that need references skip themselves. `fetch-depth: 0` for `test_legacy_frozen.py`.
- `identity.yml`: builds the `linux/amd64` image from the commit with `build_image.sh`, its layers cached in GHCR
  (`ghcr.io/<owner>/prohap-fast-buildcache`, `mode=max`, written by pushes, nightly and manual runs, read by pull
  requests), then runs `tests/ci_identity.sh` twice inside it: `verify` through the entrypoint on a subset of the
  1kG fixture (chromosomes 1 and 22, `--samples 300 --transcripts 30`), as the runner's user with the work directory
  mounted at its own path; and `tests/legacy_vs_fast.sh -g` with `PROHAP_CODE=/opt/prohap`, the image's own code,
  environment and references (fixtures and the script from the checkout). Docker with the small tier on every push
  and pull request; Docker and Apptainer 1.5.5 (SIF from `docker-daemon:`) with the small and medium tiers nightly
  and on `workflow_dispatch`, which also takes the tiers to run. The large tier stays a local release check: its
  reference run on 26,000 samples needs about 12 GB of the runner's 16 (Morgan, 2026-11-20). Release images are not
  built in CI; they are built on a workstation and pushed to Docker Hub (`mphthomas/prohap-fast`).

`legacy_vs_fast.sh -g` checks the reference's outputs (the legacy arm, decompressed) against
`tests/fixtures/expected.v113.sha256`, entries `<fixture>/<file>`; `-G` writes the entries of the fixtures that ran
and keeps the others. Fast vs legacy cannot see a change that hits both arms (the Snakefile, the references, the
environment); the golden checksums can. Both refuse to run off Linux: the macOS pandas wheel parses some floats
differently (Pinned versions). The golden file was written with `-G small medium` in the `linux/amd64` 0.3.0 image
under Rosetta on the Mac, with `PROHAP_CODE` pointing at this commit's checkout (the image's environment and
references are those of this commit: `envs/`, `prebuild_refs.py` and `legacy/` did not change). 2026-11-20: 159 outputs
over the 9 fixtures of the small and medium tiers, 33 min with `-j 4` under Rosetta; the fast arm of this commit
(F2, F3, D4) was IDENTICAL on all of them on Linux.

Running the image under Apptainer locally: Apptainer runs in a `--privileged` Ubuntu 24.04 container on Docker
Desktop (PPA `apptainer/ppa`, arm64). As root it runs the SIF (kernel squashfs mount); as another user, the SIF's
FUSE mount fails there (`exec ... failed: invalid argument`), so convert it to a sandbox first
(`apptainer build --sandbox`). `APPTAINER_TMPDIR` must be on the container's own file system, not a bind mount from
the Mac, and `docker-archive:` needs an uncompressed `docker save` tar (Apptainer 1.5.5 rejects a gzipped one with
`gzip: invalid header`). (Found in M4, 2026-11-18, with the first version of this suite.)

---

## 1. Quick wins on the legacy code (P1–P4, D3): exact, small

Apply to a copy of `legacy/src` to make the first `src/`. All five were verified byte-identical. Rebuilt in M6
(0.1.0); `tests/test_quick_wins.py` also compares each changed function with the reference on small crafted
inputs, for the unit CI.

**P1, `modules/get_haplotypes.py` (biggest win, 12 lines).** Before the loop over `haplo_combinations`:

```python
col_POS, col_REF, col_ALT, col_ID, col_INFO = (vcf_df[c].tolist() for c in ('POS', 'REF', 'ALT', 'ID', 'INFO'))
```

and replace `row = vcf_df.iloc[idx]` + `row['X']` with `col_X[idx]`. Exact because `str()` of a list element equals
`str()` of the numpy scalar `iloc` returned (`200000.0`, `nan`, `123`, `True` all match).

**P2, same file.** After `sample_info.set_index(...)`: `sex_of`, `pop_of`, `superpop_of` from `.to_dict()`;
`male_samples` from `sex_of`; `male_set = set(male_samples)` and use it for `indiv not in male_set`; in
`count_freq_pop`/`count_freq_superpop` use `pop_of[...]`/`superpop_of[...]`; make `x_autosomal_transcripts` a set
before the `apply` calls. Keep `len(male_samples)` from the list (a sample the VCF lists twice counts twice).
A sample name listed twice in the metadata makes `.loc` return a frame, and the reference stops on the
`male_samples` line (ambiguous truth value); the rebuild keeps that lookup when the index is not unique, so
`src/` stops there too.

**P3, `modules/process_haplotypes.py`.**
`haplo_groups = {tid: g for tid, g in genes_haplo_df.groupby('TranscriptID', sort=False)}` before the closure; in it
`transcript_haplotypes = haplo_groups.get(transcript_id, genes_haplo_df.iloc[0:0])`. Replace
`result_data = result_data + list(...)` (and the two sequence lists) with `.extend(...)`.

**P4, `modules/vcf_reader.py`.** `'file_content': []`, `.append(vcf_file_line)`, and
`StringIO(VCF_header + ''.join(...))` at write time.

**D3, `merge_duplicate_seq.py`.** Replace the `bisect` + `KeyWrapper` sorted list with `dict[hash] -> entry`
(append to the entry's lists on a hit), then `result = [d[h] for h in sorted(d)]`. Same order: the original list
was kept sorted by hash with unique keys.

**Rejected, P6: `low_memory=False`** on the two wide `read_csv` calls. It saves another 20% but turns off pandas'
chunk-by-chunk typing, so `.`-derived IDs such as `1e2` stay `1e2` instead of the reference's `100.0`
(measured mismatch at 100k samples).

Measured, `prohap.py` CPU-s on the benchmark block (10k / 26k samples): legacy 589 / 2,008; P2 alone 512 / 1,898;
P1 alone 128 / 492; P1–P4 101 / 406; P1–P4+P6 89 / 315 (not exact).

Rebuilt in M6 (0.1.0) and measured 2026-11-19 in the tree of each commit, one tree at a time, on the quiet Mac of
§5's legacy baseline: `python3 bench/bench.py --samples 10000,26000 --arms fast`, which runs `prohap.py` on the
benchmark block (§5) with one thread. CPU seconds of `prohap.py`, the peak RSS of its process tree, and its phases at
26k:

| `src/` | 10k CPU-s | 26k CPU-s | 26k peak RSS | 26k phases (s): sweep / haplotypes / proteins |
| --- | ---: | ---: | ---: | --- |
| legacy (§5) | 530.6 | 2,567.2 | 12.3 GB | 243.3 / 2,302.5 / 41.1 |
| P1 | 124.4 | 509.6 | 12.3 GB | 250.0 / 213.1 / 42.5 |
| P1–P2 | 99.5 | 446.0 | 12.3 GB | 251.7 / 147.5 / 43.1 |
| P1–P3 | 99.9 | 451.1 | 12.3 GB | 256.8 / 147.4 / 42.7 |
| P1–P4 | 94.1 | 421.9 | 4.7 GB | 228.4 / 147.0 / 42.3 |

P1 does most of it: the haplotype phase falls from 2,302 to 213 s at 26k. P2 takes another third off that phase. P3
changes nothing measurable on 60 transcripts: the filter it replaces is one pass over the whole haplotype table per
transcript, which shows on whole chromosomes. P4 shortens the sweep by a tenth and cuts the peak memory at 26k from
12.3 to 4.7 GB (3.4 to 2.4 GB at 10k): the reference grew each transcript's VCF text as one string, copied whenever
a line was added, and transcripts that share lines each held their own copy; the list holds each line once. Against
the design work (P1 128 / 492, P1–P4 101 / 406), 10k is 3–7% faster and 26k 4% slower. Filter and `prohap.py`
together on the block: 117.6 / 499.4 CPU-s against 555.0 / 2,643.7 (4.7× / 5.3×); genome-wide at 26k by the scale
factor, 44 CPU-h on this Mac instead of 231.

D3 on the `dupseq` fixture (400,000 records, 341,225 distinct sequences), `merge_duplicate_seq.py` alone, quiet Mac:
18.7 → 9.0 CPU-s, 1.1 GB either way, output identical. The old sorted-list insert is quadratic in the number of
distinct sequences, so the gain grows with the cohort (the design work: about 8.5 h → 4 min at 25M sequences).

---

## 2. `filter_fast.py` (F1): drop-in for `vcf_filter_fix.py`

Same arguments, same stdout lines, same decompressed output. Output gzip level can drop to 1 (content identical).

- Read input as bytes (`gzip -dc` subprocess for `.gz`, a separate process doing the decompression). Any line with
  `\r` or non-ASCII bytes → re-exec the legacy script (text mode uses universal newlines and UTF-8, bytes mode does not).
- Split off the first 9 fields at the 9th tab. Fast path only if those 9 fields contain no other whitespace
  (`len(head.split()) == 9 and head.count('\t') == 8`) and FORMAT is exactly `GT`.
- Genotype section → `normalise_gts`:
  - Fixed width (`len % 4 == 3`): append a virtual tab, view as `uint32`, compare with the word for `"0|0\t"`; only
    non-matching tokens are checked: separator `|` or `/`, both alleles in `[0-9.]`, 4th byte a tab.
  - Mixed haploid/diploid (chrX males): find tabs, token lengths must be 1 or 3; 1-char tokens become `x|0`
    (the legacy rule: add `|0` when the token has neither `|` nor `/`).
  - Anything else → that line goes through the legacy per-sample code (copied verbatim into the script).
- Biallelic: output = first 8 fields + `\tGT\t` + normalised tokens + `\n`.
- Multi-allelic, allele `j` (1-based): output fields 0–3, the allele, fields 5–6, `\tMAF=<str(allele_maf)>\tGT\t`,
  then tokens remapped by a 256-entry table: digits `1`–`9` → `0`, digit `j` → `1` (if `j <= 9`), `0` and `.`
  unchanged, applied **only to tokens whose separator is `|`** (the legacy `str.replace('g|', ...)` chain never
  touches `/` tokens). `allele_maf` is int `0` when `-af 0` (prints `MAF=0`).
- AF handling (`get_MAF`) copied verbatim. Counters and the two final `print`s identical.

Measured: 26k block 65 → 2.6 s; 500k block 185 s (decompression 21 s; most of the rest is recompressing output,
which F3, fusing the filter into prohap, removes).

**As built in 0.2.0**, as `src/vcf_filter_fix.py` itself (same name, arguments and printed lines), with these
choices:

- Input through a `gzip -dc` process; output compressed by a thread with zlib at level 1 (the conda environment's
  zlib compresses this text twice as fast as macOS's `/usr/bin/gzip -1`).
- A biallelic line with FORMAT `GT` whose tokens are all `a|b` or `a/b` is written as it was read: the reference's
  output for it is the line itself.
- FORMAT with `GT` first and more fields (`GT:DS`, ...; added after 1.0.0): the reference writes
  `token.split(':')[0]` of each token, so the fast path cuts each token at its first `:` (`first_subfield_genotypes`:
  the GT subfield is one byte long when the byte after its first is `:`, tab or newline, three when the fourth is)
  and normalises as for `GT`. Any byte below `!` other than the tabs and the final newline sends the line to
  `reference_line`, since the reference splits on any whitespace, also inside later fields. Measured on 1,000
  records × 26,000 samples of `GT:DS` (260 MB, plain text): 0.52 CPU-s, against 2.05 through `reference_line`, 7.6
  for the reference script and 0.17 for the same records as `GT` alone.
- Up to nine ALT alleles on the fast path; lines outside it (FORMAT not starting with `GT`, other tokens, other whitespace,
  ten or more alleles) go through the reference's loop body for one line, copied into the script as
  `reference_line`. A carriage return or non-ASCII byte anywhere hands the whole run to the reference.
- `<file> is empty!` is printed at the end, before the counts, so stdout is the reference's even when a later file
  hands the run over.
- `tests/test_filter_fast.py` runs both scripts on a crafted VCF per path and per trigger.

Measured 2026-11-19 with `bench.py` on the quiet Mac, one thread (CPU-s of the filter step, wall seconds in
brackets): 0.9 (0.6) at 10k and 2.8 (1.6) at 26k on the benchmark block, against the reference's 24.4 and 76.5
(§5): 27× both; 174.1 (99.9) for the block at 500k, 0.5 GB; 86.2 (49.2) for TTN at 500k.

---

## 3. `prohap_fast.py` (R1, R2): drop-in for `prohap.py`

Same arguments. Imports from `legacy/src/modules`: `check_open_file`, `read_fasta`, `check_vcf_line_validity`,
`cluster_conflicting_mutations`, `result_columns`, `empty_output`, `check_protein_allele`, `add_population_freqs`,
`get_rna_position`, `get_rna_position_simple`, `check_start_change`. Reusing them is what keeps the biology identical.

### 3.1 Setup (identical to legacy)

Transcript list via pandas, `require_start` filter, `all_transcripts.sort(key=start)` (stable). Metadata via
`pd.read_csv(sep='\t')`. `sample_ids = [c for c in vcf_colnames if c in set(metadata names)]`.
Fallback to legacy if: duplicate transcript IDs, duplicate VCF column names, header not starting with the 9 standard
columns, a sample named like a fixed column, zero samples (legacy crashes there), duplicate metadata names among
included samples, non-string population codes, no transcripts.

Population tables exactly as legacy: `groupby('Population code').count()` on the included samples, male counts with
`KeyError → 0`. Precompute `copy_labels[2*i + h] = name + ':1' / ':2'`, `colrank[vcf column] = rank in sample_ids or -1`,
`male_rank`. Read the cDNA FASTA once and keep only this chromosome's transcripts. Put all of this in a module-level
dict before creating the `Pool` so workers inherit it by `fork`; each worker opens its own `gffutils.FeatureDB`.

### 3.2 Sweep (same logic as `parse_vcf` / `add_variants_to_transcripts`)

Same queue (insert by end with `bisect_left`), same validity check, same exon overlap test
(`exon.start < pos + len(REF)` and `exon.end >= pos`, with the same `break`s), same `vcf_linecount` and `.` ID
replacement (`'\t'.join(split(maxsplit=2)[:2]) + '\t' + hex(count)[2:] + '\t' + split(maxsplit=3)[3]`, applied to the
head only). Read bytes, decode only the head. Parse the genotypes of a record once, the first time it overlaps an
exon, and share the parsed record across transcripts.

`parse_genotypes`: body length must be `4 * n_samples - 1`; view as `uint32` against `"0|0\t"`; for non-zero tokens
require tab, `|` or `/`, alleles in `[0-9.]`. `hap1 = token[0] == '1' and sep == '|'`, `hap2 = token[2] == '1' and
sep == '|'` (legacy: `startswith('1|')`, `endswith('|1')`). Keep `int32` arrays of carrier column indices only.
Anything else → fallback.

### 3.3 Column values the reference sees (pandas emulation)

Legacy pushes each transcript's lines through `read_csv` → `to_csv` → `read_csv`, with chunk size
`B` = the largest power of two strictly below `2**20 // (n_samples + 9)`, minimum 1 (pandas: `b = 1; while b * 2 <
heuristic: b *= 2`), i.e. 512 at 2k samples, 32 at 26k, 8 at 100k, 1 at 500k. Each chunk is typed separately
and chunks are combined with `concat_compat`. So:

```python
STR_NA = {'-1.#IND','1.#QNAN','1.#IND','-1.#QNAN','#N/A N/A','#N/A','N/A','n/a','NA','<NA>','#NA',
          'NULL','null','NaN','-NaN','nan','-nan','None',''}
_numlike = re.compile(r'^\s*[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d*)?\s*$|^\s*[+-]?(inf|infinity|nan)\s*$', re.I)
def is_safe_str(t):   # pandas certainly keeps t as this exact string, whatever else is in its chunk
    return not (t in STR_NA or t.lower() in ('true','false') or _numlike.match(t) or '"' in t or t != t.strip())
```

- ID column: if every token is safe, use tokens verbatim. Otherwise emulate: split into chunks of `B`; a chunk of only
  safe tokens is an object array of the tokens, any other chunk goes through
  `pd.read_csv(StringIO('ID\tX\n' + ''.join(t + '\tx\n' ...)), sep='\t')['ID'].to_numpy()` (the dummy column stops
  blank lines being skipped); `concat_compat` the chunks; format with `DataFrame({'ID': arr, 'X': 'x'}).to_csv(...)`;
  repeat the chunked read on that text; `str()` each element.
- INFO must be safe, POS must be ASCII digits ≤18 chars, REF/ALT are already validated: otherwise fallback.
- Genotype tokens always contain `|` or `/`, so they are always strings.

### 3.4 Haplotype grouping (replaces the per-person loop)

For transcript variants `r = 0..V-1`: copies `2*rank` (hap1) and `2*rank+1` (hap2), dropping hap2 of males when the
transcript is not in a PAR on chrX (same PAR test as legacy). Build `(row, copy)` pairs, `lexsort` by copy then row;
`unique` copies with start/count. Hash each copy's rows with two random 64-bit weights per row summed with
`np.add.reduceat` (wrapping), group by `(h1, h2, count)`, then **verify every member's rows equal its group's first
member's** and fall back on any mismatch. Key string = `','.join(rows)`; REF members = all included copies that carry
nothing. Sort groups by key **as strings** (legacy uses `bisect` on strings: `'12,30' < '3' < 'REF'`). Members stay
in ascending copy order (= legacy's sample order).

Then per group, legacy logic verbatim: changes `str(POS)+':'+REF+'>'+ALT`, AF from INFO (`'AF' in INFO` →
`split('AF=')[1].split(';')[0].split(maxsplit=1)[0]`, else `'-1'`), `cluster_conflicting_mutations`, conflicting
haplotypes → their sample labels go to the removed list and the row is skipped, stable sort by `int(POS)`.
Transcript with no variants → one REF row (`count = 2N`, samples `'all'`). Row index = position among this
transcript's rows including REF rows.

Count, frequency, population strings and the sample string are **lazy**: only computed for rows that reach the output.
Frequency = `count / total` (chrX non-PAR total = males + 2 × females; `0` when total is 0). Population strings: first
appearance order, `'{:.5f}'.format(count / pop_total)`, skip pops with total 0.

### 3.5 Protein haplotypes (port of `process_transcript_haplotypes`)

Copy the legacy function and change only:

- Plain `str` instead of `Seq` (slicing, concatenation and comparison behave the same).
- Cache per change string: `get_rna_position` result, reverse complement, `keep = alt != cdna[loc:loc+alt_len]`.
  Cache `get_affected_codons` on the reference by `(loc, ref_len, rf_ref, protein_start_ref)`.
- Mutated cDNA: if changes are in increasing, non-overlapping reference order (`loc_k >= loc_{k-1} + len(ref_{k-1})`,
  first `>= 0`), every ref check reads the untouched reference, so build the haplotype with one `''.join(pieces)`;
  otherwise run the legacy loop.
- `translate(seq)`: ≤30 codons via a dict of the 64 ACGT codons; longer via numpy LUT
  (`A,C,G,T → 0..3`, index `16a+4b+c`), trailing partial codon ignored; any non-ACGT character → Biopython
  `str(Seq(seq).transcribe().translate())`. The table comes from `CodonTable.unambiguous_dna_by_id[1]`.
- Sequence dedup: `dict[hash(seq)]`, then `sorted` by hash (same as the legacy sorted-list insert).
- Merged rows (collapsing fix): collect sample strings in a list and `';'.join` once at the end. Note the legacy keys:
  `str(hash(...))` in the known-reading-frame branch, `hash(...)` (int) and **overwrite, not merge** in the unknown
  branch. Keep both.
- Haplotype IDs are local row indices in the worker; the parent turns them into `prefix + '_' + hex(offset + i)`.

### 3.6 Streaming and output

Transcripts leave the sweep queue in end order, but output must follow `all_transcripts` (start) order: give each
queue entry its index, submit to `Pool.apply_async` when it leaves the queue, keep at most `4 × threads` in flight,
park finished results in `done[index]` and write while `next_index in done`. The parent adds the row offset (all
haplotype rows, including REF and later-dropped ones), keeps the FASTA accession counter, and writes:

- TSV: header once, then `pd.DataFrame(rows, columns=result_columns).to_csv(fh, header=False, index=False, sep='\t')`
  every ~20k rows. Per-batch dtype inference equals whole-table inference because every column is homogeneous.
- FASTA / cDNA FASTA records exactly as legacy; removed-samples log line per transcript only if it produced ≥1 row.
- Empty VCF → `empty_output` (no log, no cDNA file), as legacy.

### 3.7 Fallback

`os.execv` the legacy script with the same arguments (it overwrites partial outputs). Triggers: unphased or
multi-character genotypes, extra FORMAT fields, quotes, CR or non-ASCII bytes, duplicate sample or transcript IDs,
hash collisions (§3.4), and the setup and record checks in §3.1–§3.3.

### 3.8 As built in 0.2.0

`src/prohap.py` keeps the reference's arguments and printed lines (plus one line saying how many workers it uses)
and is built from these modules, each tested against the reference functions it replaces:

| Module | Does | Tested by |
| --- | --- | --- |
| `fallback.py` | `to_legacy(__file__, reason)`: `os.execv` of the `legacy/src` script with the same arguments, through `tools/run_fork.py` where fork is not the default; the reason goes to stderr after `ProHap-fast: falling back to the reference code`. `PROHAP_FALLBACK=error` stops instead (exit 3). | `test_fast_helpers.py` |
| `legacy_import.py` | `legacy.<module>`: the `legacy/src/modules` files, imported once with `modules` pointing at the reference directory | `test_fast_helpers.py` |
| `pandas_emulation.py` | §3.3 | `test_pandas_emulation.py` |
| `sweep.py` | §3.2: `VCFReader` (bytes through `gzip -dc`), `read_header`, `RecordParser`, `assign_records` | `test_sweep.py` |
| `haplotypes.py` | §3.4: `Cohort` (metadata as `get_gene_haplotypes` reads it), `group_copies`, `transcript_rows`, lazy `Row` strings | `test_haplotypes.py` |
| `proteins.py` | §3.5: `process_transcript`, `translate`, `get_affected_codons`, `add_population_freqs` | `test_proteins.py` |
| `cdna.py` | the transcripts' cDNA as `read_fasta(..., True)` gives it, without a dict of the whole file | `test_cdna.py` |

`test_prohap_fast.py` runs both scripts on a crafted chromosome (a gffutils DB built from a GTF): the fast path with
1–3 workers and the threshold and flag settings, and one case per fallback trigger. `tests/legacy_vs_fast.sh` now
fails when a fast run fell back (it would compare the reference with itself), and `verify` lists such steps.

Choices made while building it, none of which changes output:

- **The setup reads the cDNA of the listed transcripts only** (`cdna.py`): the whole FASTA in one piece, records
  found by `\n>`, every header parsed with `read_fasta`'s own rules, later records with the same ID replacing
  earlier ones. 0.6 s instead of 2.2 s for the v113 cDNA, and the workers do not inherit a dict of 250,000
  sequences.
- **Transcript annotations are read in the parent** (exons, first start and stop codon, biotype) before the `Pool`
  forks, instead of a `FeatureDB` per worker.
- **Grouping** adds two random 64-bit weights per record into arrays of one entry per haplotype copy (kept zeroed
  between transcripts), groups the carrying copies by (number of records, two sums), and checks the grouping record
  by record: each group must have all or none of its members among a record's carriers, or the run falls back. One
  stable sort then gathers each group's records. Memory is one entry per copy plus one per (group, record) pair,
  never one per (copy, record) pair (§6).
- **Conflicts**: `cluster_conflicting_mutations` joins records whose spans `[POS, POS + len(REF)]` overlap or touch,
  so it is called only for groups with at least two records that meet another record; the others cannot conflict.
- **Population strings** are counted in Python for up to 64 copies and with `np.unique` above. Merged rows add them
  with an exact re-implementation of `add_population_freqs` applied left to right (each step reads back the
  five-decimal sums the previous one wrote); population codes with `:` or `;`, which the reference splits on, fall
  back.
- **`translate`** caches the translations of short sequences (the codons around a change); the reference's
  `check_protein_allele` results are cached per transcript.
- **The workers return finished rows**; the parent names the haplotypes (`<prefix>_<hex(offset + i)>`), numbers the
  accessions and writes the table through `DataFrame.to_csv` every 20,000 rows. Output tables and FASTAs ending in
  `.gz` are written at gzip level 1.
- **Fallback triggers** (§3.7), as built: in the VCF, a carriage return or non-ASCII byte anywhere, a header not
  starting with the nine VCF columns, no samples, an empty, repeated or quoted sample name, a data line with fewer
  than ten fields, whitespace other than single tabs, quotes or control characters in the first nine fields, POS
  other than plain digits, FORMAT other than `GT`, and, in records that fall in an exon, INFO that pandas might not
  read as text or genotypes other than `a|b` with one character per allele (unphased `a/b` included, as CLAUDE.md
  lists it); in the setup, `-af` above 0, a transcript listed twice, no transcript while the VCF has records (the
  reference stops), a metadata name like a VCF column, a sample name with `:` or `;`, a sample with two metadata
  rows, a population code that is missing, not text or contains `:` or `;`, an output table compressed other than by
  gzip; while running, a start or stop codon outside the exons (the reference stops) and a hash collision. The
  reference then stops where it stops (`test_prohap_fast.py`, cases `duplicate_transcript`,
  `duplicate_metadata_row`, `no_sample_in_metadata`, `no_transcript_for_the_records`).

Measured 2026-11-19 with `bench.py` on the quiet Mac (one thread unless stated; CPU seconds of the filter and of
`prohap.py`, the wall seconds and peak memory of `prohap.py`'s process tree; 0.1.0 from §1, the reference from §5):

| Data | Samples | filter | prohap | together | prohap wall | prohap peak | 0.1.0 | reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| block | 10,000 | 0.9 | 4.3 | 5.2 | 4.0 | 0.8 GB | 117.6 | 555.0 |
| block | 26,000 | 2.8 | 7.7 | 10.5 | 6.2 | 0.8 GB | 499.4 | 2,643.7 |
| block | 500,000 | 174.1 | 121.6 | 295.7 | 76.1 | 1.9 GB | | |
| block, 4 threads | 500,000 | 171.3 | 122.2 | 293.5 | 76.7 | 4.0 GB | | |
| TTN | 26,000 | 1.4 | 20.3 | 21.7 | 20.4 | 0.8 GB | | |
| TTN | 500,000 | 86.2 | 431.2 | 517.4 | 421.6 | 6.4 GB | | |
| `n100000` chr21 / chr22 | 100,000 | 0.3 / 0.3 | 1.3 / 1.3 | 1.6 / 1.6 | 1.3 / 1.4 | 0.9 GB | | 116.7 / 138.9 |

- Against the design work's targets: 10.5 CPU-s at 26k (13.5), 296 CPU-s and 4.0 GB with four workers at 500k
  (351, 4.3 GB), TTN 21.7 CPU-s at 26k (about 30) and 517 CPU-s, 6.4 GB at 500k (546, 6.7 GB). At 26k the block
  costs a 252nd of the reference's CPU time.
- Below 100,000 samples most of the memory is the setup (pandas, gffutils, reading the cDNA file once: 0.8 GB). At
  500,000 the reader is the bottleneck: `prohap.py` takes 76 s with one worker and with four.
- TTN at 500,000: 991,022 of the 1,000,000 haplotype copies have a haplotype of their own, each with a
  35,000-residue protein; the protein port takes about 290 s of the 431. Its largest process peaks at 4.8 GB (the
  worker: the transcript's haplotypes, its distinct proteins, and the result on its way to the parent). A
  chromosome-2 job with four workers adds three other workers' transcripts to that.
- Genome-wide by the scale factor (§5, block × 315): 26 CPU-h at 500,000 samples on this Mac (the design work: 31;
  the reference cannot run at that size), 0.9 CPU-h at 26,000 (the reference: 231, 0.1.0: 44).

---

## 4. Streaming `merge_tables.py` and `haplo_extract_sample_names.py` (D1, D2)

The recipe from the design work. Stream lines. Keep each line's fields; every `chunk` lines
(`pandas_chunk_rows(21) == 32768`), **and at the end of every input file** for the merge, build the text of all columns except `samples`, run
`pd.read_csv(..., sep='\t')` then `to_csv(header=False, index=False, sep='\t')`, and splice the re-formatted values
back. `samples` is copied verbatim (must contain `:` and not be an NA token). At the end, fall back to legacy if any
column saw both int and float chunks without an object chunk (that combination re-types ints as `x.0`). Quote
characters or `\r` → fallback. Extract writes the table without `samples` and `HaplotypeID\tsamples` from the same
re-formatted fields. `usecols` does not change pandas' chunking (verified), which is why this matches.

Measured there on a 135k-row table: merge 40 → 9 s, extract 42 → 10 s; memory one chunk instead of the whole table.

### 4.1 As built in 0.3.0: `src/modules/tsv_stream.py`

- **Chunks as pandas reads them.** `Table` reads a table's lines as bytes (`gzip -dc` in its own process) in blocks of
  `pandas_chunk_rows(width)` rows (32,768 for the 21 columns of a haplotype table), starting again with each input
  file, since the reference reads each file with its own `read_csv`. Each block goes through `pd.read_csv` and
  `to_csv` without its last column. A table one column narrower has chunks at least as long, so pandas types the
  block as one chunk, which is the reference's chunk.
- **`samples` is copied, not re-parsed,** when every value of the block contains `:` (`S1:1;S2:2`): pandas keeps such
  a word as text whatever else is in its chunk, and writes it back as read. Otherwise the block goes through pandas
  whole, last column included (no fixture takes this path; `test_tsv_stream.py` does).
- **Joining chunks and tables** (`combined_dtype`). A column keeps every value's text when all its chunks have one
  dtype, or when one of them is object: an object column holds Python ints, floats and booleans, whose `str()` is
  the text pandas writes for them in their own dtype (for floats: numpy's shortest repr, checked on 60,000 floats of
  every magnitude). Tables with no rows count for nothing. Checked against pandas' `concat_compat` (the chunks of one
  file) and `pd.concat` (the tables) for every mix of int64, uint64, float64, all-missing float64, bool, object and
  an empty table, in two and three parts (`test_combined_dtype_follows_pandas`). Every other mix rewrites values and
  hands the run to the reference: int64 next to float64 (`5` becomes `5.0`; also an int64 table next to one whose
  column is all missing), bool next to int64 (`True` becomes `1`; across tables pd.concat makes it int64, within a
  file concat_compat makes it an object of ones and zeros), uint64 next to int64 (float64). The check runs at the
  end of each table and of the merge, after the output is written; the reference then writes it again.
- **Fallback triggers**: a byte other than tab, newline or printable ASCII, or a quote (which pandas would unquote);
  a line whose number of fields is not the header's (pandas skips blank lines, pads short ones and stops on long
  ones); a header with duplicate or empty names, fewer than three columns, or a last column other than `samples`; an
  input or output whose name implies compression other than gzip.
- **Memory**: each line is held once (its last field as a memoryview into it, the other fields as one short copy),
  and a block is handed to a callback rather than yielded, so nothing of the previous block is alive while the next
  is read. `GzipWriter` (moved from `vcf_filter_fix.py` to `modules/gzip_writer.py`) hands writes of more than 4 MB
  to its thread in 4 MB slices, so its queue holds at most 32 MB beyond the block being written, not eight whole
  blocks. Outputs ending in `.gz` are written at gzip level 1 (the reference: pandas' default, 9).

### 4.2 `merge_tables.py` (D1)

Writes the first table's header, then the chunks of each table named by `-i`, in order. The headers must be equal
(pd.concat would align different columns by name); each column's dtypes are joined per table, then across tables.
`tests/test_tsv_stream.py` runs both scripts on crafted tables: four chromosomes with a header-only table among
them, plain files, chunk edges at 512 rows (1,030 columns) and at 32,768 (the real 21 columns) with values whose text
depends on their chunk's type, and one case per fallback trigger and per rewriting mix (each must fall back, and the
output is then the reference's). Shortening the chunk by one row, or accepting every dtype mix, fails cases.

Measured 2026-11-19 with `bench/simulate_tables.py`: tables like `prohap.py`'s at 500,000 samples (117,000 rows and
273 MB of text per chromosome, 89% of it `samples`; the benchmark block's output, IMPLEMENTATION_NOTES §5, is 117k
rows and 293 MB). CPU seconds and peak memory; macOS natively (busy machine), Linux in the arm64 image under Docker
Desktop (its wall times are dominated by the shared folder, so CPU only):

| Tables | Text | reference, macOS | fast, macOS | reference, Linux | fast, Linux |
| ---: | ---: | --- | --- | --- | --- |
| 4 | 1.1 GB | 161 CPU-s, 1.9 GB | 10.5 CPU-s, 1.1 GB | 173 CPU-s, 1.5 GB | 15.8 CPU-s, 0.38 GB |
| 16 | 4.4 GB | | 42 CPU-s, 1.3 GB | | 59 CPU-s, 0.40 GB |

- 15× less CPU, and memory no longer grows with the table. Python's own allocations peak at 0.3 GB (tracemalloc);
  macOS's allocator keeps freed blocks of a few hundred KB, hence 1.1–1.3 GB there. Most of the fast time is gzip:
  `gzip -dc` and the level-1 compression thread, which the main thread waits for about a third of the time; pandas
  takes about 2 s per GB.
- Genome-wide at 500,000 samples (315 blocks, about 90 GB of text): the reference would need 1.4–1.7 times the text
  in memory (125–155 GB) and about 3.7 CPU-h; the fast merge 15–20 CPU-min in 0.4 GB (Linux).
- Fallback seen on real-looking data: a table whose `chromosome` column is empty (all NaN, float64) next to tables of
  numbered chromosomes (int64) hands the merge to the reference, which writes `1.0`, `2.0`, ... for them.

### 4.3 `haplo_extract_sample_names.py` (D2)

Reads the merged table in its chunks (one file, so one chunking from start to end, as the reference's
`read_table`) and writes both outputs from the same re-formatted fields: the table without its last column, and
`HaplotypeID` with `samples`. A table without `HaplotypeID` goes to the reference, which stops. Tested by the same
cases as the merge, each run through both scripts in turn.

Measured as §4.2, on the merged tables:

| Tables | Text | reference, macOS | fast, macOS | fast, Linux |
| ---: | ---: | --- | --- | --- |
| 4 | 1.1 GB | 167 CPU-s, 1.9 GB | 10.3 CPU-s, 1.2 GB | 15.9 CPU-s, 0.40 GB |
| 16 | 4.4 GB | | 41 CPU-s, 1.4 GB | 61 CPU-s, 0.42 GB |

### 4.4 `fasta_simplify_headers.py` at about 25M IDs

The workflow's last step reads the final FASTA into a dict (`read_fasta`), the haplotype table without samples into
pandas (`haplo_df.set_index('HaplotypeID')`), and for every haplotype ID of every record looks up its transcript
(`haplo_df.loc[ID]`) and the transcript's gene in the whole-genome gffutils DB. Measured 2026-11-19 with
`bench/bench_simplify.py` (MANE reference proteins plus haplotype records with one or two IDs each; the script is
the same in both trees), macOS:

| Haplotype IDs | Records | CPU-s | Peak |
| ---: | ---: | ---: | ---: |
| 100,000 | 90,415 | 4.7 | 0.4 GB |
| 1,000,000 | 690,415 | 38.9 | 2.0 GB |
| 4,000,000 | 2,690,415 | 158.4 | 7.3 GB |

Linear: 39 µs and 1.75 GB per million IDs. At 25M IDs (500,000 samples genome-wide: 37M haplotype rows, 25M FASTA
records) that is about 1,000 CPU-s, which does not dominate (the fast merge and extract take 15–20 CPU-min each),
but **about 44 GB in one process**, four times the 10 GB per job of the project's goals and the largest of the
workflow at that size (`prohap.py` peaks at 6.4 GB). About 60% is the table (20 columns of Python objects for every
row, used only for HaplotypeID → TranscriptID) and most of the rest the FASTA dict. Follow-up: stream the FASTA and
keep only that mapping, typed as pandas types the two columns (`perf(S1)`, TODO M9).

### 4.5 As built: streaming `fasta_simplify_headers.py` (S1)

- **The FASTA is read record by record** by a copy of `read_fasta`'s loop (the same tag, accession, description and
  sequence, including its quirks: a last line without a newline loses its last character, blank lines add nothing).
  The reference's dict keeps one record per accession, in the first one's place with the last one's content; the
  script keeps a 64-bit hash of every accession (8 bytes per record) and hands the run to the reference at the end
  if two are equal (a hash collision among different accessions falls back too, harmlessly).
- **The haplotype table becomes HaplotypeID -> TranscriptID only**: the IDs as one sorted numpy array of fixed-width
  byte strings, and a 32-bit code per row for its transcript (the distinct transcript texts are kept once). Lookups
  are batched per block of 16,384 records (`np.searchsorted`). No chunk typing is needed, because a lookup is only
  served when pandas certainly keeps both words as their text whatever else is in their chunk (`is_safe_str`): every
  HaplotypeID of the table must be such a word and unique (otherwise `.loc` sees NaN, numbers or a non-unique index),
  and so must each TranscriptID that is looked up. Real tables (`haplo_chr<c>_<hex>`, `ENST...`) always are. A table
  that breaks the rule is fine as long as no record looks anything up, as for the reference. The lines are checked
  as tsv_stream checks them (field counts, quotes, CR, non-ASCII), since pandas would otherwise read other rows.
- **Gene names are cached per transcript** (one gffutils query each; about 24,000 in the benchmark).
- **The header table is written with `to_csv` per block of 16,384 rows.** Every value is a Python str (no NaN, no
  numbers), so pandas formats and quotes each row on its own (a tab or `"` in an accession is quoted the same way);
  the first block carries the header, an empty FASTA gives the header alone. The FASTA output is written one record
  at a time: joining a block into one string made macOS's allocator keep the pages of every such string (RSS grew by
  50–100 MB per block of 65,536 records while Python's own memory stayed at 0.1 GB).
- **Fallback triggers**: a carriage return, control or non-ASCII byte in the FASTA; the table's triggers above; a
  header-table name implying compression other than gzip (the FASTA output follows the reference: gzip when the name
  ends in `.gz`, else plain); an input that is also an output (the reference reads everything before writing); and
  everything on which the reference stops (a header without the four fields, a haplotype ID missing from the table or
  no table, a transcript missing from the DB, a variant ID without `ENST`), which then stops in the reference.
- Tested by `tests/test_simplify_fast.py`: 33 script cases against the reference (blocks of exactly and around
  16,384 records, a table of more than 2^20 rows with the longest ID after the first array, quotes and tabs in
  accessions, custom prefixes, every fallback trigger and every way the reference stops). Dropping the accession
  check, the uniqueness check, the HaplotypeID check or the per-block header each fails a case.

Measured 2026-11-20 with `bench/bench_simplify.py`, macOS, the reference and the fast script on the same inputs
(outputs byte-identical):

| Haplotype IDs | Records | reference | fast |
| ---: | ---: | --- | --- |
| 100,000 | 90,415 | 4.6 CPU-s, 0.39 GB | 1.4 CPU-s, 0.40 GB |
| 1,000,000 | 690,415 | 40.6 CPU-s, 1.9 GB | 6.2 CPU-s, 0.38 GB |
| 4,000,000 | 2,690,415 | 167 CPU-s, 6.3 GB | 21.2 CPU-s, 0.61 GB |

About 5 CPU-s and 80 MB per million IDs beyond a 0.3 GB base (imports, the gffutils DB, one block), so about 2 CPU-min
and 2.3 GB at 25M IDs, and about 3.2 GB at 37M table rows, against about 17 CPU-min and 44 GB for the reference.

---

## 5. Benchmark harness and simulator (`bench/`)

### Simulator (`bench/simulate.py`)

```bash
python3 bench/simulate.py --chrom 22 --first ENST00000402142 --transcripts 60 --samples 26000 --out DIR
```

writes `chr22.vcf.gz` (with `--chunks K`: `chr22.1.vcf.gz` ... `chr22.K.vcf.gz`, so the workflow glob is
`chr{chr}.*vcf.gz`), `samples.tsv`, `transcripts.csv` and `chr22.json` (arguments and counts), reading the references
directory (`--refs`, default the one `ensembl.v*` in `$PROHAP_REFS`). Several runs into one directory, one per
chromosome, make one dataset: `samples.tsv` must come out the same, and `transcripts.csv` keeps the other
chromosomes' rows. The model, as rebuilt from the design work:

- Block: `--transcripts` consecutive MANE transcripts with a start codon, in order of start (the order prohap.py
  sorts them in), from `--first` or after `--skip`; 0 means to the end of the chromosome. REF bases come from the
  Ensembl cDNA mapped back through the exons, complemented on the − strand; a transcript whose cDNA is not its exons
  spliced together gives no bases (none in the MANE blocks used here).
- Sites per exonic bp `0.3 × (N/500000)**0.43`, drawn without replacement from the block's exonic bases. Common
  sites, one per 300 exonic bp: AF log-uniform 0.01–0.5, carried by founders. Each 50 kb window has 40 founder
  haplotypes with Dirichlet(0.5) weights, every haplotype copy descends from one founder per window, and a site is
  carried by founders taken in random order until their weights reach its AF. Rare sites: allele count from a
  discrete power law (exponent 1.9) up to 1% of 2N, on random copies.
- 6% deletions and 6% insertions, 1–10 bp (geometric, p = 0.5); 6% with a second ALT (an SNV on other copies, with a
  rare count). A deletion may run past the end of its exon; its intronic REF bases are random, since nothing reads
  intronic bases.
- Added in the rebuild: 10% of deletions get an SNV inside the deleted bases (or on the anchor base, as a second
  record at the same POS) on half of the deletion's carriers plus 0–2 others, so prohap.py finds conflicting
  variants and logs removed samples (about 360 copies on the benchmark block at 2.5k samples; without it a
  conflict needs a rare and a common variant to overlap on one copy, and almost never happens). A common site whose
  founders no copy descends from (only in tiny cohorts) gets rare carriers instead, so no record has AC 0.
- IDs (`--ids`): `dot`; `rs` (unique random rs numbers); `ambiguous`: 30% tokens pandas may retype, of the kinds
  `2e5 123 7 0012 NA nan True False 1.50`, the rest `.`. The rest is `.` rather than rs because `.` becomes a hex
  line number, which pandas also often reads as a number, so 8-row chunks of only numeric-looking IDs occur at 100k
  samples; with rs IDs nearly every chunk would contain a string and stay text.
- chrX outside the PARs (`x_par1_to` 2,781,479, `x_par2_from` 155,701,383): males have one copy, written haploid
  (`1`, `--males haploid`) or mirrored (`1|1`); AC and AN count what is written.
- Samples `S000001`...: 26 1kG populations, sex 50/50. `--unlisted K` leaves K VCF samples out of `samples.tsv`
  and adds K metadata rows (`M000001`...) not in the VCF.
- FORMAT `GT` only; INFO `AC=..;AN=..;AF=..` with AF `%.6g`; QUAL `.`, FILTER `PASS`. gzip level 1 with mtime 0, so
  the same arguments write the same bytes. Random streams: the samples from `(seed)`, sites and genotypes from
  `(seed, chromosome)`, IDs from a third stream, so `--ids`, `--males` and `--chunks` change only what they name.

The benchmark block is 60 transcripts of chromosome 22 from ENST00000402142 (39.57–41.9 Mb), 217,387 exonic bp.
The design work recorded only "60 consecutive MANE transcripts on chr22, 217 kb"; this block has that size. It
has 6,717 / 9,059 / 12,198 / 18,414 records at 2.5k / 5k / 10k / 26k samples (20 MB gzipped at 26k), simulated in
1–3 s.

Scale factor: genome = block × 315 (autosomal MANE exon bp 68.2 Mb / block 217 kb). chr1 is 10.8% of it.

### Harness (`bench/bench.py`)

```bash
python3 bench/bench.py --samples 2500,5000,10000,26000 --arms legacy,fast [--threads 1]
python3 bench/bench.py --data DIR --arms legacy,/path/to/other/src
```

For each dataset (`--data`, or the benchmark block at each `--samples`, simulated once into `bench/data/`),
chromosome and arm (`legacy` = `legacy/src`, `fast` = `src/`, or a directory), it runs `vcf_filter_fix.py` and then
`prohap.py` with the Snakefile's arguments (`-af 0` and thresholds 0 by default, so every variant and haplotype is
kept), `PYTHONHASHSEED=0`, through `tools/run_fork.py` on macOS. Per step:

- wall time; CPU time (user + system) from the rusage `os.wait4()` returns for the step's process, which includes the
  Pool workers it waited for. `RUSAGE_CHILDREN` of the harness would also count the `ps` calls of the sampler.
- peak RSS of the process tree, summed over the processes `ps` lists every 0.25 s (forked workers share pages, so it
  overstates, but it is what a scheduler's memory limit sees), and the largest single-process RSS (`ru_maxrss`).
- prohap.py's phases, from the time each of its progress lines appears on stdout (`PYTHONUNBUFFERED=1`): setup,
  sweep (assigning variants to transcripts: the per-transcript tables), haplotypes (co-occurrence and frequencies),
  cdna (reading the cDNA FASTA), proteins, write.

It prints a table, appends to `bench/results.tsv` (git-ignored) and, when two arms ran, compares their outputs after
decompression (a convenience; `tests/legacy_vs_fast.sh` is the identity check).

Design-work measurements (prototypes, lost), CPU-s per block, legacy filter+prohap: 63 (2.5k), 202 (5k), 609 (10k),
2,072 (26k). Fast: 4.6, 6.1, 8.5, 13.5, then 23 (50k), 43 (100k), 97 (200k), 351 (500k). Fast prohap with 4 workers at
500k: 82 s wall (reader-bound). TTN alone at 500k: 546 CPU-s, 6.7 GB. 500k block output: 117k rows, 293 MB TSV text
(83% `samples`), 79k FASTA records.

### Legacy baseline

`legacy/src` on the rebuilt simulator, one thread, `-af 0` and thresholds 0, Apple Silicon Mac with nothing else
running, 2026-11-18:

```bash
python3 bench/bench.py --samples 2500,5000,10000,26000 --arms legacy --threads 1
python3 bench/bench.py --data tests/fixtures/sim/cache/ttn --data tests/fixtures/sim/cache/n100000 \
    --arms legacy --threads 1
```

| Data | Samples | Records | filter CPU-s | prohap CPU-s | Total CPU-s | Design work | Peak RSS | prohap phases (s): sweep / haplotypes / proteins |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| block | 2,500 | 6,717 | 3.8 | 61.9 | 65.7 | 63 | 1.6 GB | 3.8 / 45.1 / 10.4 |
| block | 5,000 | 9,059 | 9.1 | 168.3 | 177.4 | 202 | 2.5 GB | 11.1 / 139.1 / 15.3 |
| block | 10,000 | 12,198 | 24.4 | 530.6 | 555.0 | 609 | 3.4 GB | 37.5 / 467.6 / 23.1 |
| block | 26,000 | 18,414 | 76.5 | 2,567.2 | 2,643.7 | 2,072 | 12.3 GB | 243.3 / 2,302.5 / 41.1 |
| `ttn` (TTN) | 2,500 | 3,388 | 1.7 | 142.2 | 143.9 | | 3.5 GB | 4.6 / 107.3 / 30.3 |
| `n100000` chr21 (3 transcripts) | 100,000 | 431 | 5.4 | 111.3 | 116.7 | | 5.5 GB | 49.3 / 61.5 / 0.1 |
| `n100000` chr22 (3 transcripts) | 100,000 | 451 | 5.6 | 133.3 | 138.9 | | 7.3 GB | 57.1 / 76.5 / 0.2 |

Peak RSS is the process tree's (prohap's own `ru_maxrss` is 0–12% lower: forked workers share pages); filter stays
under 100 MB. Setup, cDNA and write take 0.2–2.2 s each. Wall time is CPU time + 0–4%.

- The haplotype phase (the per-person loop over the wide per-transcript tables) is 73–90% of prohap on the block and
  TTN and grows as N^1.7 (2.5k → 26k: 51×); the sweep, which builds those tables, as N^1.8 (64×); filter as N^1.3
  (20×). At 100k on 3 short transcripts the sweep is 44% of prohap.
- 2.5k–10k are within 12% of the design work. 26k costs 28% more CPU and 48% more memory (8.3 GB there); the design
  work did not record its record count or how it measured memory, so the cause is unknown. These numbers are the
  baseline from now on.
- Genome-wide at 26k by the scale factor: 2,644 × 315 = 833k CPU-s, 231 CPU-h on this Mac.
- At 100k samples a chromosome with 3 transcripts (3 kb of exons, about 440 records) already needs 5–7 GB: the
  per-transcript tables are N columns wide however few records they have.

---

## 6. Traps that broke identity during development

- Writing results in queue-pop (end) order instead of `all_transcripts` (start) order shifts every haplotype ID.
- Testing pandas emulation only at small N: chunk size is 512+ rows at 2k samples, so the width-dependent typing never
  shows. Include 100k-sample fixtures (8-row chunks) and a `--transcripts` mode in `verify`.
- `merge_tables.py` is not a neutral concatenation (float re-parse); a plain `cat` differs on `frequency`.
- `low_memory=False` looks harmless and is not.
- Comparing `.gz` files instead of decompressed bytes reports false differences.
- Comparing macOS output with Linux golden checksums: the macOS pandas wheel parses some floats differently (see
  Pinned versions).
- On macOS the default `spawn` start method breaks the legacy code's `Pool` of closure functions; use `fork`.
- prohap's stdout is not reproducible: Pool workers print `not in cDNA database, skipping!` and `Ref allele not
  matching` lines into block-buffered stdout, so their position changes between two runs of the same code (seen on
  1kG chr22 with all 8,237 transcripts, 2026-11-17). Compare stdout as a set of lines, never byte for byte. The
  removed-samples log is written by the parent and is stable.

Found while building 0.2.0:

- pandas reads booleans without regard to case (`tRuE` is `True`), and its float parser is not correctly rounded
  beyond 10^22 (`3e23` becomes `2.9999999999999997e+23`). `is_safe_str` is conservative and `_kind` types only plain
  integers, decimals of at most 15 digits and exponents up to 22 itself (`test_pandas_emulation.py`).
- `concat_compat` of a boolean chunk and an integer chunk casts the booleans to 1 and 0, then the column to object,
  so `True` reads back as `1`. Calling pandas' own `concat_compat` keeps this.
- A one-column `to_csv` writes an empty value as `""` (the csv module marks a row with one empty field), so the
  emulation writes a second, dummy column.
- Forked workers inherit stdout's unwritten buffer, and a worker that prints would write it out again; `prohap.py`
  flushes stdout before it creates the `Pool`.
- The fast scripts find `legacy/src` and `tools/` next to their own directory after resolving links, so `verify
  --fast-src` copies the directory instead of linking it.
- At 500,000 samples nearly every copy of TTN is its own haplotype (991,022 groups, 56M (copy, record) pairs): a
  grouping that sorts all pairs held 3.4 GB; scatter-adding hashes into arrays of one entry per copy and checking
  the groups record by record holds 1.5 GB.

---

## 7. Upstream behaviours the reference keeps

Found while checking the collapsing fix (2026-11-17). All five are in upstream `cafcc51`, `legacy/src`
reproduces them, and the fast implementation must too. The maintainer decided to keep them and report them
upstream; fixing one is a deliberate output change made in both trees. Line numbers are for
`legacy/src/modules/process_haplotypes.py` with `legacy-fixes.patch` applied.

1. **The liftover check drops most indels** (line 212). A change is dropped when its ALT allele equals the
   reference cDNA at its position, a check meant for alleles that became the reference in a liftover.
   `get_rna_position` (`coordinates_toolbox.py:9`) keeps the VCF anchor base, so on a + strand transcript every
   deletion matches (its ALT is the anchor base), as does an insertion whose inserted bases repeat the reference
   bases after it. On the − strand the comparison uses reverse complements, so a deletion is dropped only when its
   last REF base equals the anchor base. A haplotype whose only change is a dropped indel gets no row; one with
   other changes is keyed without the indel and merges with the haplotype that lacks it. 1kG, MANE, threshold 0:
   the haplotype tables of chr1, 2, 3, 21 and 22 contain no + strand deletion, against 10 on the − strand. Indels
   dropped on the + strand cause every merge on chr21 (6), chr3 (2) and chr22 (1), 4 of 8 on chr2 and 3 of 10 on
   chr1; UTR SNVs cause the rest.
2. **Stop-gained variants are filtered as UTR** when `haplo_ignore_UTR` is 1 (`check_protein_allele`, line 49,
   applied at 387 and 395). The coding sequence is cut at the first stop codon of the mutated protein (line 380),
   so the variant that creates that stop falls outside it. A haplotype whose only protein change is the stop-gain
   gets no row; with other changes it is keyed without the stop-gain and usually merges with the full-length
   haplotype, and the truncated protein never reaches the FASTA. Crafted test (Q45L; Q45L + stop; stop alone):
   three rows with `haplo_ignore_UTR` 0, one Q45L row with count 2 with 1.
3. **Merges are keyed on VCF ID strings, not alleles** (line 408). `vcf_filter_fix.py:111` gives every allele of a
   split multi-allelic record the record's ID, and pandas reads NA-like IDs (`NA`, `nan`, ...) as the same `nan`
   (`get_haplotypes.py:163`). Different variants with the same ID string therefore merge: count and frequencies
   add up, the row shows only the first haplotype's changes, and the second protein is not written to the FASTA
   (sequences are stored only when a row is created, from line 444). Since the collapsing fix (CHANGELOG: "Merged
   haplotypes list all their samples"), the row also lists the second haplotype's samples.
4. **The unknown-reading-frame branch overwrites instead of merging** (from line 463; int key at 486, assignment at
   488). A later haplotype with the same key replaces the earlier one's row, and the earlier one disappears from
   the table and, after the header filter at 522-531, from the FASTA. Happens only with `haplo_require_start: 0`
   (transcripts with neither start nor stop codon) or `haplo_skip_start_lost: 0` (start lost); the example and TRE
   configs use 1 for both.
5. **No cDNA FASTA for a chromosome with no variants.** `empty_output` (line 36) writes only the table and the
   protein FASTA, so with `haplo_cdna_file` set the `merge_fasta` rule fails: `cat: ...haplo_cdna_chr22.fa: No such
   file or directory`, under bash strict mode. Checked on 1kG chr22 with `phased_min_af: 1.1`.

---

## 8. Input preparation and defaults (M9)

### 8.1 `restrict_to_exons` (F2)

`prohap.py` gives a transcript a record when `exon.start < POS + len(REF)` and `exon.end >= POS` for one of its
exons (`add_variants_to_transcripts`), so a record that meets no exon of a listed transcript never reaches any
output, except through the line numbers that name `.` IDs. With whole-genome input those records are most of the
file, and both the filter and `prohap.py` read every one of them.

`tools/restrict_to_exons.py`, run by the `restrict_to_exons` rule for both arms when the config sets it:

- Regions: the exons of every transcript of the list on the chromosome (as `prohap.py` selects them, before its
  start-codon check, so a superset), from the chromosome's gffutils DB, sorted and merged.
- VCF chunks (bgzip + .tbi/.csi): `tabix -h -u -R regions.tsv` (1-based, inclusive regions file). tabix uses
  `[POS, POS + len(REF) - 1]` (or INFO END), the rule above. `-u` writes a record that meets several regions once;
  without it, a long deletion across a gap between two merged regions comes out twice, out of order
  (`test_crafted_edges`). Only the first chunk's header, as the filter keeps it.
- BCF chunks (.csi): `bcftools view -R --regions-overlap record -Ov`. bcftools rewrites the text (header lines,
  float formatting), which does not matter for BCF, which the reference cannot read at all. Not used for VCF text
  for that reason (decision of M3, `verify/vcfcut.py`).
- Output: bgzip level 1, `data/vcf/phased/chr<c>_exons.vcf.gz` (temp); the filter reads it instead of the glob.
- The contig is the chromosome's name with or without `chr` (`MT`/`chrM` alike); records of other contigs are
  dropped with a message, an unindexed chunk is an error naming the option.

Both arms see the same records, so `verify` holds; against a run without the option only `.` IDs change (and
the stdout counts of the filter). `tests/test_restrict_to_exons.py` checks the selection against the rule on
crafted and random chunks, and runs the 1kG fixture with IDs set: fast with the option = fast without (every output
but the filtered VCF), legacy with it = fast with it.

### 8.2 The filter fused into `prohap.py` (F3)

The legacy arm writes the filtered VCF gzipped (`filter_phased_vcf`), and `prohap.py` decompresses and reads it
again. In the fast arm `compute_haplotypes` runs `prohap.py -filter_glob GLOB -filter_af AF -filter_af_field FIELD`
instead of `-i FILE`, and `filter_phased_vcf` does not run (the Snakefile branches on `implementation`):

- `prohap.py` starts `vcf_filter_fix.py ... -o -` (`src/modules/fused.py`, `FilterReader`), which writes the
  filtered lines uncompressed to a pipe and its messages to stderr. The lines are the file's lines, so the line
  numbers that name `.` IDs are the same. The filter runs in its own process, beside `prohap.py`'s reader, which is
  the bottleneck at cohort scale; nothing is compressed or decompressed in between.
- The pipe can be read once, so the header is read once and the stream stays open for the sweep (the unfused path
  opens the file twice).
- `-o -` cannot hand input to the reference (the reference would write a file named `-`): the filter stops with exit
  status 3, which `FilterReader.finish` turns into `Unsupported`. On `Unsupported` from either side,
  `to_legacy_fused` runs `vcf_filter_fix.py` (which falls back itself if it must) into `-tmp_dir`, then the
  reference `prohap.py` on that file, removes the file and exits with the reference's status. `PROHAP_FALLBACK=error`
  stops with status 3 as before.
- The filtered VCF is no longer an output of the fast arm. `verify` and `tests/legacy_vs_fast.sh` no longer compare
  it (decision 2026-11-20); the fast filter is still checked against the reference filter by
  `tests/test_filter_fast.py`, including `-o -` against the file it writes otherwise.

Tests: `tests/test_prohap_fast.py` (`test_fused*`) runs the fused script against the reference's two steps on the
crafted cases, including a carriage return (the filter stops, the fallback runs both steps), unphased genotypes
(`prohap.py` stops), `-af` and an AF threshold in the filter.

Measured on the benchmark block (60 transcripts of chromosome 22, `-threads 4`, filter at AF 0), the filter and
`prohap.py` as two steps against the fused `prohap.py`, Mac with a load of about 5 from other jobs (2026-11-20):

| Samples | two steps: wall / CPU / peak | fused: wall / CPU / peak |
| --- | --- | --- |
| 26,000 | 4.8 s / 11.0 s / 0.81 GB | 4.6 s / 10.1 s / 0.84 GB |
| 500,000 | 187 s / 309 s / 3.8 GB | 150 s / 229 s / 4.1 GB |

At 500,000 samples the filtered text is 139 GB (69,551 lines of 2 MB) and 0.83 GB gzipped: the fused run saves its
compression and decompression (a quarter of the CPU) and overlaps the filter with the sweep, so the step takes about
as long as `prohap.py` alone did (1.25× wall here, less than the design's 2×, which assumed the two steps took
about as long). The filter's process adds about 0.3 GB.

### 8.3 `merge_duplicate_seq.py` and `remove_UTR_seq.py` record by record (D4)

Found while sizing 1.0.0: both scripts read the whole FASTA into a dict (`read_fasta`). On `bench/simulate_fasta.py
--records 3000000` (1.06 GB gzipped), D3's `merge_duplicate_seq.py` peaked at 7.9 GB (2.6 GB per million records,
so about 65 GB at the 25M records of 500,000 samples genome-wide), and the reference `remove_UTR_seq.py` at 3.5 GB on
its output.

- `src/modules/fasta_stream.py`: `records()` is `read_fasta`'s loop yielding (element ID, record) in file order;
  `Accessions` keeps a 64-bit `hash()` of every element ID and says at the end whether one repeated. `read_fasta`'s
  dict would have kept the last of such records in the place of the first, so a stream is exact only without
  repeats (equal hashes of different IDs count as a repeat: safe, slower).
- `merge_duplicate_seq.py`: the reference writes the merged entries in the order of `hash(sequence)`, numbered
  `prot_<hex>` in that order. Each record, reduced to (hash, tag, accession, sequence and the four parsed fields), is
  pickled (batches of 500, gzip level 1) into one of 128 temporary files in the output's directory by
  `(hash + 2**63) * 128 >> 64`, so the files cover consecutive hash ranges and each keeps input order. Then file by
  file: merge in a dict, sort by hash, write with the running entry number. A repeated accession: the run is redone
  the D3 way, in memory.
- `remove_UTR_seq.py`: each record is written as it is read; a repeated accession hands the run to the reference
  (its input comes from `merge_duplicate_seq.py`, whose accessions are unique).

At 3M records, identical output: `merge_duplicate_seq.py` 0.40 GB and 72 s (D3: 7.9 GB, 70 s), `remove_UTR_seq.py`
0.11 GB and 16 s (reference: 3.5 GB, 17 s); at 25M records about 1 GB. Tests: `tests/test_fasta_stream.py`
(`records` against `read_fasta` on odd files, both scripts against the reference with every header kind, a repeated
accession, an empty file); the `dupseq` fixture runs both at 400,000 records.


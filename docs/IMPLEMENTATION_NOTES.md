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
environment, with no container. For each 1kG fixture config it runs the workflow twice on the same inputs, once
with `--config implementation=legacy` and once with `--config implementation=fast`, each in its own work directory
(switching `implementation` in place reruns nothing). Then it compares every file under the two directories, except
`.snakemake/`, after decompression. A file that is gzip in one run and plain in the other is DIFFERENT.

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

## 1. Quick wins on the legacy code (P1–P4, D3): exact, small

Apply to a copy of `legacy/src` to make the first `src/`. All five were verified byte-identical.

**P1, `modules/get_haplotypes.py` (biggest win, 12 lines).** Before the loop over `haplo_combinations`:

```python
col_POS, col_REF, col_ALT, col_ID, col_INFO = (vcf_df[c].tolist() for c in ('POS', 'REF', 'ALT', 'ID', 'INFO'))
```

and replace `row = vcf_df.iloc[idx]` + `row['X']` with `col_X[idx]`. Exact because `str()` of a list element equals
`str()` of the numpy scalar `iloc` returned (`200000.0`, `nan`, `123`, `True` all match).

**P2, same file.** After `sample_info.set_index(...)`: `sex_of`, `pop_of`, `superpop_of` from `.to_dict()`;
`male_samples` from `sex_of`; `male_set = set(male_samples)` and use it for `indiv not in male_set`; in
`count_freq_pop`/`count_freq_superpop` use `pop_of[...]`/`superpop_of[...]`; make `x_autosomal_transcripts` a set
before the `apply` calls. Keep `len(male_samples)` from the list.

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

---

## 4. `merge_fast.py` / `extract_fast.py` (D1, D2)

Stream lines. Keep each line's fields; every `chunk` lines (`pandas_chunk_rows(21) == 32768`), **and at the end of
every input file** for the merge, build the text of all columns except `samples`, run
`pd.read_csv(..., sep='\t')` then `to_csv(header=False, index=False, sep='\t')`, and splice the re-formatted values
back. `samples` is copied verbatim (must contain `:` and not be an NA token). At the end, fall back to legacy if any
column saw both int and float chunks without an object chunk (that combination re-types ints as `x.0`). Quote
characters or `\r` → fallback. Extract writes the table without `samples` and `HaplotypeID\tsamples` from the same
re-formatted fields. `usecols` does not change pandas' chunking (verified), which is why this matches.

Measured on a 135k-row table: merge 40 → 9 s, extract 42 → 10 s; memory one chunk instead of the whole table.

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

### Harness

- Harness: run each arm with `PYTHONHASHSEED=0`; sample summed RSS of the process tree every 0.25 s; CPU from
  `RUSAGE_CHILDREN`; timestamp stdout lines to split legacy phases (sweep / haplotypes+freqs / proteins / write).

Measured CPU-s per block, legacy filter+prohap: 63 (2.5k), 202 (5k), 609 (10k), 2,072 (26k). Fast: 4.6, 6.1, 8.5,
13.5, then 23 (50k), 43 (100k), 97 (200k), 351 (500k). Fast prohap with 4 workers at 500k: 82 s wall (reader-bound).
TTN alone at 500k: 546 CPU-s, 6.7 GB. 500k block output: 117k rows, 293 MB TSV text (83% `samples`), 79k FASTA records.

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

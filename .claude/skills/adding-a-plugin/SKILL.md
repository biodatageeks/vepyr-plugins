---
name: adding-a-plugin
description: >-
  Add a new VEP-style plugin cache to vepyr-plugins (manifest + build + validate +
  upload), given a source URL or local file. Use whenever the user wants to add,
  integrate, or wire up a new annotation source (e.g. "add gnomAD as a plugin",
  "integrate REVEL scores") into the vepyr plugin_cache pipeline — even if they
  don't use the word "plugin" explicitly. Covers manifest authoring, source-format
  handling (CSV/TSV/VCF/Parquet), memory-safe per-chromosome building, validation,
  and Drive upload.
---

# Adding a new plugin to vepyr

A "plugin" here is: a `<name>.source.toml` manifest in `vepyr-plugins/plugins/<name>/`
that tells the `datafusion-bio-function-vep` builder how to turn a raw external
annotation source (position/allele → score(s)) into a per-chromosome parquet
shard vepyr can look up at annotation time. No Rust changes needed for a normal
plugin — only the manifest.

Reference implementations already in this repo: `plugins/cadd/cadd.source.toml`
(native TSV, two combined sources, `assume_unique`), `plugins/spliceai/spliceai.source.toml`
(VCF flattened to TSV, per-transcript `match_column`), `plugins/clinvar`,
`plugins/alphamissense`, `plugins/dbnsfp`. Read whichever is closest to your new
source before writing the manifest — copy the shape, don't invent one.

## 0. Before touching anything: source triage

Given a URL or file, answer these first — they decide almost everything else:

1. **What's the raw format?** CSV/TSV (native columns, e.g. CADD/dbNSFP/AlphaMissense)
   or VCF (INFO-field-packed, e.g. ClinVar/SpliceAI/gnomAD)?
   - `bcftools view -h <file> | head -50` or `zcat <file> | head -5` tells you fast.
2. **Is it tabix-indexed, and is a `.tbi`/`.csi` alongside it?** If yes, you can
   query it **remotely over HTTP** without downloading the whole file —
   `tabix <https-url> <chrom>` works directly against a public bucket/FTP as long
   as the matching index is fetchable at the same URL + `.tbi`. Proven this
   session against `krishna.gs.washington.edu` (CADD) and works the same way for
   any bgzip+tabix source (gnomAD's GCS bucket, Ensembl FTP, etc.) — **always try
   this before downloading a multi-GB/TB source locally.**
3. **How big is it, per chromosome?** Check the file listing (for GCS-hosted data:
   `curl -s "https://storage.googleapis.com/storage/v1/b/<bucket>/o?prefix=<prefix>" | python3 -m json.tool`
   lists objects+sizes with no auth for public buckets). This tells you whether a
   chromosome is even feasible on this machine — see the size table below.
4. **Does one row = one variant, or can the same runtime probe key --
   `(start, allele_string, <match columns, if any>)`, NOT just
   `(start, allele_string)` -- repeat with different values?** Get the key
   wrong and you'll silently keep the wrong row for a duplicate key —
   reason about the FULL key explicitly, don't guess.
   - CADD (no match column, key is just `(start, allele_string)`):
     structurally cannot repeat — all-possible-SNVs plus gnomAD-normalized
     indels, both unique by construction. Safe for `assume_unique = true`.
   - SpliceAI (match column `symbol`): the bare `(start, allele_string)`
     key DOES repeat at overlapping-gene loci -- confirmed empirically,
     101,919 duplicate `(pos,ref,alt)` keys in a 9.1M-row sample, e.g. the
     same variant scored once per overlapping gene. It's only safe for
     `assume_unique = true` because the FULL key includes `symbol`, which
     never repeats per variant. Checking `(start, allele_string)` alone
     here would have set the wrong flag.
   - AlphaMissense (match column `protein_variant`): overlapping UniProt
     entries CAN repeat even the full key with different scores — needs
     dedup, `assume_unique` must stay unset.

### Known source scale (build feasibility on a 16GB machine, learned the hard way)

| Source | Format | Biggest chrom (compressed) | Feasible? |
|---|---|---|---|
| CADD SNV | TSV | chr1 ≈ 14GB flat | yes (with streaming write, ~2-4.5h) |
| SpliceAI (masked) | VCF→TSV | chr1 ≈ 2.3GB flat | yes (~4.5h worst case, but usually <1h) |
| gnomAD v4.1 **exomes** | VCF | chr1 ≈ 19GB compressed | borderline — expect *worse* than CADD chr1 (VCF INFO parsing is heavier per row than flat TSV columns) |
| gnomAD v4.1 **genomes** | VCF | chr1 ≈ 44GB compressed | **no** — will not fit locally even with remote streaming |

If a candidate source's biggest chromosome is meaningfully bigger than CADD's
~14GB benchmark, say so explicitly and propose starting with the *smallest*
chromosome (usually chr21 or chrY) as a timing test before committing to the rest.

## 1. Source format handling

- **Native CSV/TSV** (delimited columns already, e.g. `chrom\tpos\tref\talt\tscore`):
  use `provider = "csv"` directly in `[[source]]`, pointing at the file (or the
  per-chrom flattened file — see step 2). No preprocessing needed.
- **VCF** (INFO-packed, e.g. `AC=...;AN=...;AF=...`): two options —
  - **Native VCF provider** (`provider = "vcf"`) — `ProviderKind::Vcf` is wired
    (`datafusion-bio-functions/.../plugin_cache/provider.rs`, landed with the
    streaming-write/VCF-provider PR): point `[[source]]` straight at the raw
    VCF and use `ingest_sql` to pull the INFO subfields you need as named
    columns (quote them, e.g. `"CLNSIG" AS clnsig` — VCF INFO tags are
    case-sensitive column names in the generated schema, and an unquoted
    identifier gets lowercase-folded by the SQL parser and won't resolve).
    It also doesn't split multiallelic records itself — it pipe-joins ALT
    alleles into one string (e.g. VCF's `A,C` becomes `"A|C"`), the exact
    same trap the flatten path's `bcftools norm -m -` step guards against;
    `UNNEST(string_to_array(alt, '|'))` in `ingest_sql` replicates that split.
    **Known gap, confirmed empirically (2026-07-29) on ClinVar**: the VCF
    provider's scan can run multi-partition, and the tier-inheritance LEFT
    JOIN (`plugin_cache/join.rs`) doesn't guarantee it preserves the probe
    side's row order the way a single-partition CSV/TSV scan implicitly does
    — the builder's `assert_start_monotonic` guard (by design) then refuses
    to write, with "tier shard write is not position-ascending". This needs
    a Rust-side fix (constrain the VCF-sourced scan to a single partition, or
    add an explicit re-sort after the tier join) before `provider = "vcf"` is
    safe to rely on for a real build — **use the flatten fallback below until
    that lands**, even though the native provider parses correctly today.
  - **Flatten fallback** (always works, proven on ClinVar + SpliceAI, and the
    currently-recommended default given the VCF-provider gap above): run
    `bcftools norm -m -` FIRST to split any multiallelic record into one
    biallelic record per ALT, THEN
    `bcftools query -r <chrom> -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/<FIELD1>\t...\n'`
    to explode the specific INFO subfields into a headerless TSV per chromosome,
    then feed that to a `provider = "csv"` source exactly like a native-TSV plugin.
    Skipping `bcftools norm -m -` is a silent-miss trap even if the source you're
    testing against happens to have zero multiallelic records today: a raw
    multiallelic `%ALT` comes back comma-joined (e.g. "A,C"), the `allele_string`
    built from it never matches a single-allele runtime probe, and that record
    quietly never annotates anything. This is also required if the INFO field is
    itself pipe/comma-packed (SpliceAI's masked `SpliceAI` tag needed a second
    `awk` pass after `bcftools query` to split 9 sub-values out of one INFO tag —
    look at `spliceai.source.toml`'s header comment for the exact one-liner).
- **Parquet**: `provider = "parquet"` already works, straightforward.

## 2. Per-chromosome flatten + sort (only if source isn't a single native file)

If you flattened from VCF, or the source has multiple files that combine per
chromosome (CADD's SNV + indel are two separate tabix'd files), **you must
globally position-sort the combined flat file before building**:

```bash
gsort -t $'\t' -k2,2n -S 1G --parallel=4 raw.tsv > sorted.tsv   # NOT coreutils `sort` — 10x+ slower on multi-GB files
```

Why: the builder's streaming write (see step 3) assumes the input arrives in
position-ascending order and skips an explicit in-memory sort for that reason —
if two source files are each internally sorted but concatenated (SNV block then
indel block), the *combined* file is NOT globally sorted and the on-disk shard
will silently violate its `(tier, start)`-sorted contract unless you `gsort` first.
(Root-caused this session on CADD; validated by diffing old vs. new builds.)

A single native VCF/TSV queried by `bcftools`/`tabix` for one chromosome IS
already globally position-sorted — no `gsort` needed in that case.

## 3. Build one chromosome

```python
import vepyr
result = vepyr.build_plugin_cache(
    '<plugin_name>', '<version>',          # version = the vepyr-plugins git ref/tag to resolve
                                            # <plugin_name>.source.toml from (NOT a
                                            # datafusion-bio-functions ref — that's fixed by
                                            # whatever datafusion-bio-function-vep build vepyr
                                            # itself is linked against).
    source_path='<path to flattened/sorted TSV, or raw file if no flattening needed>',
    cache_dir='<core Ensembl cache dir, e.g. .../116_GRCh38_merged>',
    plugin_cache_root='<plugin cache output root>',
    chroms=['<chrom>'],                     # one chromosome per call — never pass the whole genome at once
    plugins_repo='<path to this vepyr-plugins checkout>',
)
```

**`plugins_repo` + `version` resolve via `git worktree add <version>`, not your
working tree.** Editing a `.source.toml` locally and re-running the build
without committing first silently rebuilds from the *old, committed* manifest
— you'll see stale-schema errors that look unrelated to your edit (e.g. a CSV
field-count mismatch after switching `provider = "tsv"` to `"vcf"`). Commit
the manifest change in `plugins_repo` (a local commit is enough, no push
needed) before iterating.

Memory is bounded (~3-4GB RSS regardless of chromosome size, confirmed on
chromosomes up to 482M rows) as of the streaming-write fix in
`datafusion-bio-function-vep`'s `build.rs` — if you see RSS climbing past ~5-6GB
and not plateauing, something is wrong (report it, don't just add more RAM/timeout).
**Time, not memory, is the binding constraint** for big chromosomes — give `gtimeout`
a generous budget (8h, not 4h) for anything near CADD-chr1 scale; a timeout
(`rc=124`) is not a failure, just under-budgeted, and safely retryable.

`uv run` sometimes silently reuses a stale build even after a manifest/Rust
change — if the log doesn't show `Building vepyr @ file://...`, force it with
`uv sync --reinstall-package vepyr` once, then subsequent `uv run` calls pick it
up correctly.

## 4. Validate before trusting the result

Never assume a build is correct just because `rc=0`. If you have (or can quickly
build) an existing good shard for the same chromosome, diff them — but **sort
both by the full logical key first** (`tier, start, allele_string, <match_cols>`)
before comparing, since physical row order among ties is allowed to differ
(the runtime lookup is a pure hash-key probe, order-independent) and a naive
row-by-row diff will show false mismatches:

```python
import pyarrow.parquet as pq
old = pq.read_table('old.parquet')
new = pq.read_table('new.parquet')
# `match_cols` = this plugin's manifest [[match_column]] names, e.g.
# ['symbol'] for SpliceAI, ['protein_variant'] for AlphaMissense, [] for a
# per-variant plugin. Sorting by tier/start/allele_string alone leaves rows
# tied on those three in arbitrary relative order for a match-column plugin
# -- .equals() would then report a false mismatch on a logically-identical
# shard just because tied rows landed in a different physical order.
match_cols = []  # fill in from the manifest under test
keys = [(k, 'ascending') for k in ('tier', 'start', 'allele_string', *match_cols)]
assert old.sort_by(keys).combine_chunks().equals(new.sort_by(keys).combine_chunks())
```

If there's no prior-good shard to diff against (first time adding this plugin),
at minimum: spot-check row counts against `wc -l` on the flattened source, and
manually probe a handful of known variants against the source file.

## 5. Upload

```bash
rclone copy <plugin_cache_root>/plugin/<name>/chr<N>.parquet gdrive:plugin_cache/plugin/<name> --drive-root-folder-id <id> --low-level-retries 50
rclone copy <plugin_cache_root>/plugin/<name>/manifest.json  gdrive:plugin_cache/plugin/<name> --drive-root-folder-id <id> --low-level-retries 50
```

**`rclone copy` can silently "complete" without uploading on this network** (seen
repeatedly this session — process exits, no error, but the file isn't on Drive).
Always verify with `rclone lsf gdrive:... | grep '^chr<N>\.parquet$'` after every
upload before deleting the local copy or moving to the next chromosome.

Once a chromosome is confirmed uploaded, its local copy is safe to delete to
free disk for the next one — this machine's disk is the tighter constraint
than Drive storage.

## Common failure signatures (don't misdiagnose these)

| Symptom | Real cause | Fix |
|---|---|---|
| `rc=137` | genuine OOM (memory truly exceeded) | real capability failure — usually not retryable as-is |
| `rc=124` | `gtimeout` expired | not a crash — retry with a bigger timeout, same command |
| `rc=1` + "No space left on device" | disk exhaustion (often from macOS swapfile growth, not your data) | free disk (delete already-uploaded local shards), retry |
| Build process shows 0% CPU, no child `rustc`, for 10-30+ min | **often NOT stuck** — this project's release compiles have long legitimate gaps between compilation units; on this machine we killed a live build twice this session because of this, wasting hours | wait longer before concluding it's hung; only kill if truly no CPU/child activity for 30-45+ min AND `lsof` on the process shows nothing informative |
| `uv run` hangs holding `.venv/.lock` / `.cache/uv/*.lock`, 0% CPU, small RSS | genuine `uv` lock deadlock (different from the above) | `lsof -p <pid> \| grep lock` to confirm, then `pkill -9` the process tree and relaunch |

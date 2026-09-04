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
(native VCF, a packed `Number=.` INFO tag split with `array_element`+`split_part`,
per-transcript `match_column`), `plugins/clinvar/clinvar.source.toml` (native
VCF, `array_to_string` on `Number=.` fields, correlated multiallelic
expansion),
`plugins/alphamissense`, `plugins/dbnsfp`. Read whichever is closest to your
new source before writing the manifest — copy the shape, don't invent one.

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
5. **Where does it come from, and what are its exact bytes?** Every
   `[[source]]` MUST declare `url` and `md5` (the validator rejects a manifest
   without them):
   - `url` is the **canonical upstream** download — the publisher's FTP /
     bucket / release page, never a mirror, a Drive share or a local path. If
     the top-level file is a moving target (ClinVar's weekly `clinvar.vcf.gz`),
     pin the dated release you actually built from (`archive_2.0/<year>/
     clinvar_<date>.vcf.gz`); match its `.md5` against your copy to find which
     date that is. If the build input is a local re-compression of the
     upstream file (AlphaMissense's BGZF+tabix rebuild of a plain gzip), `url`
     and `md5` still name the upstream file — do NOT add a second digest key
     to the manifest. Write the preprocessing steps (commands, tabix columns,
     and the digest the derived file is expected to hash to, as information)
     in a `README.md` beside the manifest, following
     `plugins/alphamissense/README.md`. `build_plugin_cache` hashes
     `source_path` against `md5` before ingesting anything, so a derived
     input must be built with `verify_source=False` (or `"warn"`); say so in
     that README.
   - `md5` comes from the **publisher first**: CADD ships `MD5SUMs`, ClinVar a
     `.md5` beside every VCF, GCS exposes `md5Hash` (base64 → hex) via
     `https://storage.googleapis.com/storage/v1/b/<bucket>/o/<object>`, dbNSFP
     ships a `.md5` with its VEP-ready file. Check the source's `# Source:`
     URL and the plugin's header in `Ensembl/VEP_plugins` for the download
     location. Only when upstream publishes nothing (Ensembl's
     `variation_plugins/` FTP dir for SpliceAI) compute it on your downloaded
     copy, confirm the byte size matches the server's `Content-Length`, and
     say so in a comment.
   - A copy already on Drive doesn't need re-downloading to hash: Drive stores
     the MD5, and `rclone lsjson --hash --hash-type md5 -R gdrive-mw: --drive-root-folder-id <id>`
     returns it without transferring the file. The Drive MCP `get_file_metadata`
     does NOT expose it.

### Known source scale (build feasibility on a 16GB machine, learned the hard way)

| Source | Format | Biggest chrom (compressed) | Feasible? |
|---|---|---|---|
| CADD SNV | TSV | chr1 ≈ 14GB flat | yes (with streaming write, ~2-4.5h) |
| SpliceAI (masked) | VCF (native) | chr1 ≈ 2.3GB flat-equivalent | yes (~4.5h worst case, but usually <1h) |
| gnomAD v4.1 **exomes** | VCF | chr1 ≈ 19GB compressed | borderline — expect *worse* than CADD chr1 (VCF INFO parsing is heavier per row than flat TSV columns) |
| gnomAD v4.1 **genomes** | VCF | chr1 ≈ 44GB compressed | **no** — will not fit locally even with remote streaming |

If a candidate source's biggest chromosome is meaningfully bigger than CADD's
~14GB benchmark, say so explicitly and propose starting with the *smallest*
chromosome (usually chr21 or chrY) as a timing test before committing to the rest.

## 1. Source format handling

- **Native CSV/TSV** (delimited columns already, e.g. `chrom\tpos\tref\talt\tscore`):
  use `provider = "csv"` directly in `[[source]]`, pointing at the raw file or
  an unmodified per-chromosome slice. No preprocessing needed.
- **VCF** (INFO-packed, e.g. `AC=...;AN=...;AF=...`): use the native VCF
  provider (`provider = "vcf"`). `ProviderKind::Vcf` is wired
  (`datafusion-bio-functions/.../plugin_cache/provider.rs`), self-contained,
  no bcftools/awk. Point `[[source]]` straight at the raw VCF and use
  `ingest_sql` to pull the INFO subfields you need:
    - Quote the column, e.g. `"CLNSIG" AS clnsig` — VCF INFO tags are
      case-sensitive column names in the generated schema, and an unquoted
      identifier gets lowercase-folded by the SQL parser and won't resolve.
    - Wrap `Number=.` (multi-value) INFO fields with `array_to_string(..., ',')`
      — the provider exposes them as `List<Utf8>` even when a given record only
      has one element, so a bare reference renders as `[Uncertain_significance]`
      instead of `Uncertain_significance`. A `Number=1` (scalar) field must
      NOT be wrapped — `array_to_string` on a non-list column is a type error;
      check the field's `Number=` in the VCF header before deciding.
    - It doesn't split multiallelic records itself — it pipe-joins ALT alleles
      into one string (e.g. VCF's `A,C` becomes `"A|C"`). Split in SQL, but do
      **not** unnest ALT alone when a selected INFO field is allele-indexed.
      `Number=A` has one value per ALT; `Number=R` starts with the REF value and
      then has one value per ALT. Build the ALT list, remove the first element
      from each `Number=R` list, `arrays_zip` ALT with every allele-indexed INFO
      list, and unnest the zipped structs so the indexes remain correlated:

      ```sql
      WITH aligned AS (
          SELECT ...,
                 string_to_array(alt, '|') AS alts,
                 "AF" AS af_values, -- Number=A
                 array_slice("R_TAG", 2, array_length("R_TAG")) AS r_values
          FROM plugin_example_src
      ),
      exploded AS (
          SELECT ...,
                 unnest(arrays_zip(alts, af_values, r_values)) AS allele
          FROM aligned
      )
      SELECT ...,
             allele['c0'] AS alt,
             allele['c1'] AS af,
             allele['c2'] AS r_value
      FROM exploded
      ```

      Before building, run a cardinality query and require zero records where
      `array_length(alts) != array_length(Number=A)` or
      `array_length(alts) + 1 != array_length(Number=R)`; `arrays_zip` pads a
      shorter list with null and must not be used to conceal malformed input.
      A `Number=1` scalar is copied to every exploded allele. `Number=.` is
      source-defined: inspect its specification instead of assuming its values
      correspond to ALT indexes. ALT-only `unnest` is safe only when no emitted
      value is allele-indexed (as in the currently verified ClinVar release,
      which has zero multiallelic records).
    - `ALT = "."` (VCF's "no called alternate allele" marker) renders as an
      empty string here, not the literal `.` character — the provider parses
      VCF semantics, it doesn't echo raw text. A separately generated TSV or
      another tool that prints raw `%ALT` may keep the literal `.` instead. If
      cross-checking a manifest
      against a bcftools-derived reference, expect `allele_string` to differ
      on exactly these records — that's a source-parsing difference, not a
      bug in either path. See `clinvar.source.toml` for a full manifest using
      this provider.
    - A single INFO tag packing several sub-values into one delimited string
      (e.g. SpliceAI's `ALLELE|SYMBOL|DS_AG|...`) splits the same way:
      apply `split_part(..., '|', N)` to each packed entry. A `Number=.` header
      only says the list has variable length; it does **not** promise one entry.
      Use `array_element("TAG", 1)` only after a cardinality query proves
      `array_length("TAG") <= 1` for every record (as verified for the current
      SpliceAI source). If any record has several entries, preserve them with
      `unnest("TAG")`, or zip them with ALT when the source specification says
      the entries are allele-indexed, before splitting the packed fields. Never
      silently discard elements 2..N. See `spliceai.source.toml`.
    - Input and join order are irrelevant to shard correctness. The final
      DataFusion query always applies `ORDER BY tier, start`; its external
      sorter can spill within the bounded build pool, and the writer asserts
      monotonicity before publishing the shard.
- **BED** (`chrom, start, end` + optional extra columns): `provider = "bed"` —
  `ProviderKind::Bed` is wired via `datafusion-bio-format-bed`'s
  `BedTableProvider`. Its schema is only ever `chrom, start, end, name`
  regardless of the file's actual column count (BED4/5/6 select how many raw
  columns the reader parses per line, not how many get exposed) — a source
  needing more than one extra field packs it into `name` (e.g. `id|score`)
  and splits it back out in `ingest_sql` with `split_part`, the same trick
  used for a packed VCF INFO tag (see the VCF section above).
- **Parquet**: `provider = "parquet"` already works, straightforward.

## 2. Keep transformation and ordering in DataFusion

Do not flatten or externally sort a source to satisfy the cache builder.
Express allele normalization, packed-field expansion, filtering, and source
combination in `ingest_sql`. If a source requires shell reshaping to build, that
is an engine/provider gap to fix, not part of the manifest contract.

### Combine several source files in SQL

`[[source]]` is a list. Give every file a unique `part`; the builder registers
it as `plugin_<name>_src_<part>`, and `ingest_sql` can combine the tables:

```toml
[[source]]
provider = "csv"
part     = "snv"
path     = "whole_genome_SNVs.tsv.gz"
url      = "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
md5      = "88577a55f1cd519d44e0f415ba248eb9"   # from upstream MD5SUMs
  [source.csv]
  # ...

[[source]]
provider = "csv"
part     = "indel"
path     = "gnomad.genomes.r4.0.indel.tsv.gz"
url      = "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz"
md5      = "4b9c685c96d396af4d001c2f7dd9d8f9"   # from upstream MD5SUMs
  [source.csv]
  # ...
```

```sql
WITH combined AS (
    SELECT * FROM plugin_example_src_snv
    UNION ALL
    SELECT * FROM plugin_example_src_indel
)
SELECT ... FROM combined
```

Pass one real path per part when building (this requires the vepyr multi-source
API introduced with `source_path: str | dict[str, str]`):

```python
vepyr.build_plugin_cache(
    "example", "<ref>",
    source_path={"snv": ".../whole_genome_SNVs.tsv.gz",
                 "indel": ".../gnomad.genomes.r4.0.indel.tsv.gz"},
    ...
)
```

A single-source manifest still takes a plain string. A mapping must cover every
declared part and contain no unknown parts; otherwise a placeholder path could
be read accidentally, so vepyr rejects it. There is no combined temporary file
and no external sort.

### Row order is the engine's responsibility

The builder's final tier query owns physical ordering with an explicit
`ORDER BY tier, start`, independent of source order, execution partitions, or
the selected join algorithm. Its DataFusion sorter can spill within the bounded
build pool. A storage-boundary monotonicity assertion prevents publishing a
corrupt shard if that contract is ever violated. The engine regression
`sparse_plugin_with_disordered_source_still_builds_sorted` explicitly builds
from descending input and verifies the published Parquet shard is ascending.

An external position sort would not be sufficient anyway: the required `start`
is computed *after* `ingest_sql`. For example, CADD-style anchor trimming can
shift an indel from raw `pos = 100` to `start = 101` while SNVs at raw position
100 remain at `start = 100`. Only the DataFusion query knows the transformed
sort key.

A raw per-chromosome `tabix` slice is still useful to avoid repeatedly scanning
a whole-genome multi-gigabyte source; that is I/O pruning, not transformation.
For a multi-source manifest, slice each part independently and pass the paths in
the mapping above. Do not concatenate or sort the slices.

If a needed reshape cannot be expressed by the existing SQL functions or
native providers, extend the engine/provider surface rather than creating a
second preprocessing implementation whose semantics can drift.

## 3. Build one chromosome

```python
import vepyr
result = vepyr.build_plugin_cache(
    '<plugin_name>', '<version>',          # version = the vepyr-plugins git ref/tag to resolve
                                            # <plugin_name>.source.toml from (NOT a
                                            # datafusion-bio-functions ref — that's fixed by
                                            # whatever datafusion-bio-function-vep build vepyr
                                            # itself is linked against).
    # plain path for one [[source]]; {part: path} for several [[source]] entries
    source_path='<raw source or unmodified per-chrom slice>',
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
at minimum: reconcile row counts with an independent query of the raw source,
and manually probe a handful of known variants against that source.

When comparing `CSQ` output against Ensembl VEP, remember that the position of
a plugin's block is set per run, not by the manifest: pass the same order to
`annotate(plugins=[...])` that the VEP run gave its `--plugin` flags (the golden
116 runs use `spliceai, cadd, alphamissense, dbnsfp, clinvar`). `plugins=None`
falls back to alphabetical plugin-name order. Do not add a `csq_rank` key —
the validator rejects it.

## 5. Upload

```bash
rclone copy <plugin_cache_root>/plugin/<name>/chr<N>.parquet gdrive:plugin_cache/plugin/<name> --drive-root-folder-id <id> --low-level-retries 50
rclone copy <plugin_cache_root>/plugin/<name>/manifest.json  gdrive:plugin_cache/plugin/<name> --drive-root-folder-id <id> --low-level-retries 50
rclone check <plugin_cache_root>/plugin/<name> gdrive:plugin_cache/plugin/<name> \
  --drive-root-folder-id <id> --one-way --download \
  --include 'chr<N>.parquet' --include 'manifest.json'
```

**`rclone copy` can silently "complete" without uploading on this network** (seen
repeatedly this session — process exits, no error, but the file isn't on Drive).
The filename alone proves nothing when an older remote object already exists.
Require the filtered `rclone check --download` above to exit zero with no
differences or errors; it compares the fresh local bytes for both the chromosome
shard and `manifest.json` with their remote objects.

Only after that content check passes is the chromosome's local copy safe to
delete to free disk for the next one — this machine's disk is the tighter
constraint than Drive storage.

## Common failure signatures (don't misdiagnose these)

| Symptom | Real cause | Fix |
|---|---|---|
| `rc=137` | genuine OOM (memory truly exceeded) | real capability failure — usually not retryable as-is |
| `rc=124` | `gtimeout` expired | not a crash — retry with a bigger timeout, same command |
| `rc=1` + "No space left on device" | disk exhaustion (often from macOS swapfile growth, not your data) | free disk (delete already-uploaded local shards), retry |
| Build process shows 0% CPU, no child `rustc`, for 10-30+ min | **often NOT stuck** — this project's release compiles have long legitimate gaps between compilation units; on this machine we killed a live build twice this session because of this, wasting hours | wait longer before concluding it's hung; only kill if truly no CPU/child activity for 30-45+ min AND `lsof` on the process shows nothing informative |
| `uv run` hangs holding `.venv/.lock` / `.cache/uv/*.lock`, 0% CPU, small RSS | genuine `uv` lock deadlock (different from the above) | `lsof -p <pid> \| grep lock` to confirm, then `pkill -9` the process tree and relaunch |

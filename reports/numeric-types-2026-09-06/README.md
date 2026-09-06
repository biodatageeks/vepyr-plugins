Numeric plugin cache types — 2026-09-06

**Numeric storage is feasible, but most score columns need an explicit CSQ
rendering policy.** The ClinVar source manifest now stores `id` as Int32
instead of Utf8 while preserving the VCF body exactly. Six dbNSFP scalar columns were
also rebuilt as Float32 on chromosomes 22, 1 and 7: their native VCF output
differs, but a formatter using only the numeric output restores the exact VEP
body MD5 on all three chromosomes. CADD raw scores require Float64; Float32
loses information in a small number of real cache rows.

The implemented schema change is ClinVar `id`; float migrations require
additional engine support. The scripts and evidence here document both the
ClinVar validation and the float experiments. CSQ remains a textual VCF field;
the types discussed below are the physical plugin Parquet column types.

**Scope and method.** Inputs were the existing caches at
`/Users/mwiewior/workspace/data_vepyr/plugin_cache_116/plugin` and the installed
vepyr checkout at `/Users/mwiewior/research/git/vepyr`. All relevant score
columns were scanned in full on chromosomes 22, 1 and 7, with distinct strings
weighted by their actual occurrence counts. A standalone Rust probe uses the
same Float32 `Display` operation as the engine. Nulls are counted separately;
multi-value strings are recorded as invalid scalar conversions, not discarded.

| Plugin | Rows scanned across the three chromosomes |
|---|---:|
| CADD | 1,301,851,859 |
| SpliceAI | 554,824,622 |
| dbNSFP | 14,630,253 |
| ClinVar | 713,196 |

Additionally, all 25 existing ClinVar shards were checked: **4,439,569 IDs,
zero nulls, zero invalid Int32 casts, zero string round-trip changes**, range
2–4,857,410. AlphaMissense's score and SpliceAI's four position columns are
already numeric; unchanged AlphaMissense and SpliceAI caches participated in
every five-plugin VCF comparison.

**Decisions by field.** “With formatting” below means a future engine output
policy or the experimental postprocessor, not an option currently supported
by the source manifests.

| Plugin / cache columns | Audited cache type | Numeric candidate | Exact-output requirement and evidence |
|---|---|---|---|
| ClinVar `id` | Utf8 | Int32 | Direct conversion works now; full-cache ID check and three native body-MD5 passes. Applies to the pinned source; recheck IDs when changing the source release. |
| AlphaMissense `am_pathogenicity` | Float32 | Keep Float32 | Already numeric. |
| SpliceAI `dp_ag`, `dp_al`, `dp_dg`, `dp_dl` | Int32 | Keep Int32 | Already numeric. |
| SpliceAI `ds_ag`, `ds_al`, `ds_dg`, `ds_dl` | Utf8 | Float32 with formatting | Exactly two decimals, preserving negative zero. Zero formatted round-trip differences across all four complete columns on all three chromosomes. |
| CADD `cadd_phred`; dbNSFP `cadd_phred` | Utf8 | Float32 with formatting | Three decimals for values below 10, two below 20, one below 30, zero otherwise. Zero differences with this rule in both plugins. |
| CADD `cadd_raw`; dbNSFP `cadd_raw` | Utf8 | Float64 with formatting | Six decimals, preserving negative zero. Float32 fails even with six decimals; Float64 has zero differences on all scanned rows. Float64 currently requires engine support. |
| dbNSFP `metasvm_score`, `metalr_score` | Utf8 | Float32 with formatting | Exactly four decimals; preserve negative zero. Complete-column round trips and three rebuilt-cache/formatted-body MD5 passes. |
| dbNSFP `phylop100way_vertebrate`, `phastcons100way_vertebrate` | Utf8 | Float32 with formatting | Exactly six decimals; preserve negative zero. Complete-column round trips and three rebuilt-cache/formatted-body MD5 passes. |
| dbNSFP `gerp_rs` | Utf8 | Float32 with formatting | Source-style spelling: retain `.0` for integer values; nonzero magnitudes below 0.001 use uppercase `E`, an unpadded exponent and at least one mantissa decimal. Zero differences across 14,207,783 non-null cells; included in the rebuilt-cache experiment. |
| dbNSFP `sift4g_score`, `polyphen2_hdiv_score`, `polyphen2_hvar_score`, `mutationtaster_score`, `provean_score`, `vest4_score`, `revel_score` | Utf8 | Keep Utf8 under the current contract | These are frequently lists, including partial/all-missing lists such as `.,.`. Each column has 11,927,266 non-scalar cells across these chromosomes. A scalar float or `TRY_CAST` would lose content; typed lists would need a separate schema and serialization design. |
| Class, prediction, gene-symbol and other ClinVar fields | Utf8 | Keep Utf8 | Categorical/textual values. |

**Why precision alone does not fix MD5.** The engine's `ValueType`
contract supports only Utf8, Float32 and Int32. The authoritative engine source
for this run is the pinned `b145bc6` checkout:

- [source_manifest.rs](https://github.com/biodatageeks/datafusion-bio-functions/blob/b145bc6e3f4b4e757fe2d1e8fd8a5bd5e8e944ff/datafusion/bio-function-vep/src/plugin_cache/source_manifest.rs#L33): the value types and `ValueColumn` contract; there is no per-value output formatter.
- [write.rs](https://github.com/biodatageeks/datafusion-bio-functions/blob/b145bc6e3f4b4e757fe2d1e8fd8a5bd5e8e944ff/datafusion/bio-function-vep/src/plugin_cache/write.rs#L18): maps the declared type to the physical Parquet schema.
- [lookup.rs](https://github.com/biodatageeks/datafusion-bio-functions/blob/b145bc6e3f4b4e757fe2d1e8fd8a5bd5e8e944ff/datafusion/bio-function-vep/src/plugin_cache/lookup.rs#L312): decodes Float32, Int32 and text, with no Float64/Decimal path.
- [csq.rs](https://github.com/biodatageeks/datafusion-bio-functions/blob/b145bc6e3f4b4e757fe2d1e8fd8a5bd5e8e944ff/datafusion/bio-function-vep/src/plugin_cache/csq.rs#L32): renders floats using shortest round-trip `format!("{v}")`.

Consequently `1.420` becomes `1.42`, `0.00` becomes `0`, and `-0.00` becomes
`-0`. Widening to Float64 would not restore trailing zeros. Fixed precision
does not suit every field either: CADD PHRED and GERP have the rules above.

Float32 also has an actual precision limit. Across the three chromosomes,
six-decimal rendering changes **14 CADD raw rows** and **16 dbNSFP CADD raw
rows**; dbNSFP can contain the same genomic score under multiple amino-acid
contexts. Examples are `16.042138 → 16.042137` and
`-21.700286 → -21.700287`. Float64 plus six-decimal rendering has zero such
differences. These complete-cache checks matter: an HG002 body-MD5 pass alone
would not exercise every cached score.

**Actual rebuild and body-MD5 evidence.** Rebuilds used
`vepyr._core.build_plugin_cache`, SQL casts, the native tiering stage and native
Parquet writer. Their sources were the existing per-chromosome Parquet files,
not freshly downloaded upstream TSVs/VCFs. Keys and text fields were projected
unchanged; the builder regenerated the tiers. The rebuilt shards have the
expected numeric physical types and original row/warm/cold counts.

The edited production ClinVar manifest was also exercised directly with a
small tabix-indexed synthetic VCF. The native builder produced Int32 IDs and
preserved SNP/deletion/insertion normalization and explicit-dot versus absent
INFO values; see [clinvar_manifest_ingest.json](clinvar_manifest_ingest.json).

The annotation runs used the complete existing HG002 chromosome inputs,
`everything=True`, the merged release-116 core cache and the plugin order
`spliceai, cadd, alphamissense, dbnsfp, clinvar`. MD5 hashes contain every
unmodified decompressed body-line byte, including line endings. Headers and
compressed-file/Parquet bytes are excluded. Both the VEP reference and the
previous vepyr output were checked.

| Chromosome | Variants | VEP body MD5 | ClinVar Int32, native output | Six dbNSFP Float32 fields + experimental formatter |
|---|---:|---|---|---|
| 22 | 50,861 | `8c6ae1e17452776316ab3ee48def0ebf` | Match | Match |
| 1 | 323,430 | `2086a63b8676b82d3b52e97ca23533e3` | Match | Match |
| 7 | 234,522 | `f887ef32e2f2e05f1a72f675f1185dbf` | Match | Match |

The six dbNSFP columns were MetaSVM, MetaLR, GERP, phyloP, phastCons and CADD
PHRED; CADD raw stayed Utf8. Their native, unformatted output fails MD5 on
all three chromosomes. The experimental formatter changes 261, 1,134 and 496
VCF records respectively and then recovers the reference hashes. It reads
only that native numeric output and the CSQ header; it never consults original
plugin cells or reference values. This proves the feasibility of a rendering
hook, not that stock vepyr already provides one.

A separate chr22 negative control converted all seven dbNSFP scalar columns,
including CADD raw, to Float32. Its native body digest is
`73d28b785866ef549c8740a62fe5e502`, instead of the reference hash above.

Numeric types are not automatically smaller in this layout: the three ClinVar
Int32 shards are approximately 2.8–2.9% smaller; the six-Float32 dbNSFP shards
are approximately 1% larger than their original shards. These are observed
file sizes, not a performance benchmark.

**SQL casts and Decimal.** Numeric-to-Decimal-to-VARCHAR rendering works with
the policies above, but the final expression is Utf8. Putting it in
`ingest_sql` with an Utf8 value declaration preserves a string cache. Declaring
the value Float32 casts the rendered string back to a number and loses the
format again. Deployment therefore needs output rendering after numeric
lookup, plus Float64 support for CADD raw.

For a numeric Float32 SpliceAI score `v`, a tested two-decimal rendering is:

```sql
CASE
  WHEN CAST(v AS DECIMAL(18,2)) = 0
       AND starts_with(CAST(v AS VARCHAR), '-')
  THEN concat('-', CAST(CAST(v AS DECIMAL(18,2)) AS VARCHAR))
  ELSE CAST(CAST(v AS DECIMAL(18,2)) AS VARCHAR)
END
```

The sign guard is essential. Decimal storage alone cannot distinguish
`-0.00` from `0.00`; both occur in the caches. The guard uses decimal equality
because the tested DataFusion floating equality distinguishes negative and
positive zero. The same expression with scale six and `v` as DOUBLE preserves
CADD raw values. The full SQL probes have zero mismatches with the recommended
type, scale policy and sign guard; Float32 CADD raw retains the 14/16 precision
failures. GERP's custom policy was checked with Rust rather than a SQL rule.

Full SQL probes used the installed Python DataFusion 50.1.0. Representative
signed-zero, fixed-scale and CADD precision cases were also run through
vepyr's own DataFusion 53 native builder, confirming the same outcomes. Those
version differences and native results are recorded in the evidence.

**Implementation and remaining work.** The ClinVar source change uses
`CAST(id AS INT) AS id` in the final ingest projection and `type = "Int32"`
for its value column. For the float migrations, first add a per-value
serialization policy to the engine/source contract and propagate it through
plugin lookup to CSQ emission. Defaults should preserve current behavior.
Add Float64 to manifest parsing, schema writing and scalar decoding/rendering
before migrating either CADD raw column. Preserve nulls and signed zero, and
retain the multi-valued dbNSFP fields as text under the current contract.

**Reproduction and artifacts.** The complete generated caches, experimental
manifests, VCFs and logs are under
`/Users/mwiewior/workspace/data_vepyr/numeric_type_audit_20260906`.
The JSON files beside this report include complete-column profiles,
build schemas/counts, comparison hashes, SQL expressions/results and package
and checkout provenance. The environment's installed vepyr version was 0.4.0,
with bio-functions pinned to `b145bc6e3f4b4e757fe2d1e8fd8a5bd5e8e944ff`.

Run from the vepyr-plugins repository, using the existing vepyr environment:

```sh
audit_python=/Users/mwiewior/research/git/vepyr/.venv/bin/python
audit_data=/Users/mwiewior/workspace/data_vepyr
audit_output="$audit_data/numeric_type_audit_reproduction"

POLARS_MAX_THREADS=4 "$audit_python" scripts/audit_numeric_types.py \
  --cache-root "$audit_data/plugin_cache_116" --out "$audit_output/full"

POLARS_MAX_THREADS=4 "$audit_python" scripts/probe_numeric_sql.py \
  --cache-root "$audit_data/plugin_cache_116" --out "$audit_output/sql" \
  --native-cache-dir "$audit_data/cache/116_GRCh38_merged"

for audit_mode in clinvar_int32 dbnsfp_formatted; do
  "$audit_python" scripts/verify_numeric_cache.py \
    --data-root "$audit_data" --vepyr-repo /Users/mwiewior/research/git/vepyr \
    --out "$audit_output/$audit_mode" --mode "$audit_mode"
done

"$audit_python" scripts/verify_numeric_cache.py \
  --data-root "$audit_data" --vepyr-repo /Users/mwiewior/research/git/vepyr \
  --out "$audit_output/dbnsfp_float32" --mode dbnsfp_float32 --chroms 22
```

Cache-to-cache rebuilding validates the type migration, tiering and annotation
path. Full upstream raw-source parsing was not repeated; the production VCF
ingest check above covers a synthetic fixture. Float rendering results
cover every row on the three scanned chromosomes; only ClinVar ID conversion
was checked on all chromosomes. Production float rollout should retain these
checks and extend the complete-cache scan to the other contigs.

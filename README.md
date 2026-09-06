# vepyr-plugins

Plugin **source manifests** for [vepyr](https://github.com/biodatageeks/vepyr) —
a Rust reimplementation of Ensembl's Variant Effect Predictor, exposed as a
Python library.

A manifest is a declarative TOML file that tells vepyr's cache builder how to
turn a raw external annotation source (position + allele → score(s)) into a
frequency-tiered, per-chromosome Parquet cache, and which `CSQ` output fields
its values are emitted as. Adding a plugin needs **no Rust and no Python** —
only a `plugins/<name>/<name>.source.toml` file in this repository.

📚 **Full plugin documentation: <https://biodatageeks.org/vepyr/plugins/>**

| | |
|---|---|
| [Plugins](https://biodatageeks.org/vepyr/plugins/) | What a plugin cache is, building it, annotating with it |
| [Manifest structure](https://biodatageeks.org/vepyr/plugins/#manifest-structure) | Per-block reference for every key used in this repo |
| [Allele matching](https://biodatageeks.org/vepyr/plugins/#allele-matching-exact-vs-minimised) | `exact` vs `minimised`, and why it is a statement about upstream |
| [Supported plugins](https://biodatageeks.org/vepyr/plugins/#supported-plugins) | CSQ fields each plugin emits, in emitted order |
| [Download caches](https://biodatageeks.org/vepyr/downloads/#plugin-caches) | Prebuilt release-116 plugin caches on Hugging Face |
| [`build_plugin_cache()`](https://biodatageeks.org/vepyr/api/#vepyr.build_plugin_cache) | API reference |

## Manifests in this repository

| Plugin | Manifest | Raw source | Match | CSQ fields | Prebuilt cache |
|---|---|---|---|--:|---|
| **CADD** v1.7 | [`plugins/cadd`](plugins/cadd/cadd.source.toml) | 2 tabix TSVs (`snv` + `indel`) | per variant | 2 | ✅ |
| **SpliceAI** | [`plugins/spliceai`](plugins/spliceai/spliceai.source.toml) | VCF (packed `SpliceAI` INFO tag) | `{SYMBOL}` | 9 | ✅ |
| **AlphaMissense** | [`plugins/alphamissense`](plugins/alphamissense/alphamissense.source.toml) | tabix TSV | `{ref_aa}{Protein_position}{alt_aa}` | 2 | ✅ |
| **ClinVar** | [`plugins/clinvar`](plugins/clinvar/clinvar.source.toml) | VCF (`--custom`-style) | per variant | 6 | ✅ |
| **dbNSFP** | [`plugins/dbnsfp`](plugins/dbnsfp/dbnsfp.source.toml) | tabix TSV (505 columns) | `{ref_aa}/{alt_aa}` | 19 | ❌ licence |

All five are validated against golden Ensembl VEP 116 output. Four have a
prebuilt cache published on Hugging Face — see
[Plugin caches](https://biodatageeks.org/vepyr/downloads/#plugin-caches).
dbNSFP's licence forbids redistributing a converted cache, so that one is built
from your own registered download; the manifest itself is public.

## Using a manifest

vepyr resolves `plugins/<plugin>/<plugin>.source.toml` from this repository at
the git ref you pass as `version`, so different plugins can be pinned to
different revisions:

```python
import vepyr

vepyr.build_plugin_cache(
    plugin="alphamissense",                   # directory under plugins/
    version="v0.2.0",                         # git ref of THIS repo
    source_path="AlphaMissense_hg38.tsv.gz",  # the raw DATA (not stored here)
    cache_dir="/data/116_GRCh38_merged",      # Ensembl cache; supplies tiering
    plugin_cache_root="/data/plugin_cache",   # output: plugin/<name>/chr*.parquet
    # plugins_repo="/path/to/vepyr-plugins",  # local clone, for offline builds
)
```

The repository is cloned on demand; pass `plugins_repo` with a local clone to
build fully offline. A manifest with several `[[source]]` parts (CADD) takes a
`{part: path}` mapping as `source_path`. Then point `annotate()` at the same
`plugin_cache_root` and the plugin's fields appear in `CSQ` — see
[Annotating with plugins](https://biodatageeks.org/vepyr/plugins/#annotating-with-plugins).

### CSQ block order

Where a plugin's block sits in the `CSQ` string is a property of the
annotation run, not of the manifest: `annotate(plugins=[...])` emits the
blocks in the order the names are listed, and `plugins=None` (every plugin
under `plugin_cache_root`) falls back to alphabetical plugin-name order,
mirroring how Ensembl VEP orders `--plugin` flags. Manifests therefore carry
no rank; `field_order` only orders the fields *within* one plugin's block.

The golden Ensembl VEP 116 comparison the manifests are validated against
uses this order — pass it explicitly when reproducing the parity gate:

```python
plugins=["spliceai", "cadd", "alphamissense", "dbnsfp", "clinvar"]
```

## String scores and numeric casts

Some score columns deliberately use `Utf8` (Polars `String`) in the cache.
VEP emits the source's spelling: `1.420`, `0.00` and `-0.00` carry formatting
that a generic float-to-string conversion changes. Preserving those strings
lets vepyr match the **VCF body MD5**, including trailing zeros and negative
zero. Other strings contain several scores, such as dbNSFP's `0.01;.;0.03`,
and cannot be represented by one scalar float.

For filtering, aggregation and plotting, create numeric columns in your
DataFrame. A cast suitable for numerical analysis does **not** promise an
identical VCF string when cast back: binary floats approximate decimal values,
and formatting must be handled separately. Keep the original text columns if
you need to reproduce the source output. CADD raw scores need `Float64` to
retain their published six-decimal precision: a real value such as
`16.042138` becomes `16.042137` when stored as Float32 and printed to six
decimals. Decimal types can retain fixed scale, but cannot preserve the sign
of zero on their own.

CADD and ClinVar have one scalar value per variant. SpliceAI, AlphaMissense
and dbNSFP have one list element per consequence, aligned with `Consequence`.
The cache manifest determines each scalar or list element's initial type;
`lf.collect_schema()` shows the types in your build and cache.

The table uses **DataFrame/CSQ field names**, which are case-sensitive:
`CADD_RAW` belongs to CADD, while `CADD_raw` belongs to dbNSFP. Cache Parquet
files instead use each manifest's `column` names, such as `cadd_raw`.

| Plugin / DataFrame fields | Cast for numerical analysis | Why / limits |
|---|---|---|
| CADD `CADD_RAW` | `pl.Float64` | Retains the source's six-decimal precision; avoid Float32. |
| CADD `CADD_PHRED` | `pl.Float32` or `pl.Float64` | Float32 is sufficient for the source precision. VCF formatting uses three decimals below 10, two below 20, one below 30, then integers. |
| SpliceAI `SpliceAI_pred_DS_AG`, `DS_AL`, `DS_DG`, `DS_DL` (same prefix) | `pl.List(pl.Float32)` or `pl.List(pl.Float64)` | Two-decimal scores; the VCF spelling, including `-0.00`, is separate from the numeric value. |
| SpliceAI `SpliceAI_pred_DP_AG`, `DP_AL`, `DP_DG`, `DP_DL` (same prefix) | Already `pl.List(pl.Int32)` | Signed integer positions; no float cast is needed. |
| AlphaMissense `am_pathogenicity` | Already `pl.List(pl.Float32)` | Numeric score; `am_class` remains text. |
| dbNSFP `CADD_raw` | `pl.List(pl.Float64)` | Same precision requirement as CADD `CADD_RAW`. |
| dbNSFP `MetaSVM_score`, `MetaLR_score`, `GERP++_RS`, `phyloP100way_vertebrate`, `phastCons100way_vertebrate`, `CADD_phred` | `pl.List(pl.Float32)` or `pl.List(pl.Float64)` | One scalar score per consequence element. Float32 retained the source precision in the audited data, but decimal/scientific output formatting varies by field. |
| dbNSFP `SIFT4G_score`, `Polyphen2_HDIV_score`, `Polyphen2_HVAR_score`, `MutationTaster_score`, `PROVEAN_score`, `VEST4_score`, `REVEL_score` | Split each element into a numeric list; see below | An element may itself contain multiple scores and missing-value markers. Casting the whole string to one float loses that structure. |
| ClinVar `ClinVar` | `pl.Int32` for numeric IDs | All 4,439,569 IDs in the pinned source round-trip exactly through Int32. The initial dtype depends on the cache version. This is the plugin ID, not the input VCF's `id` or the core `clinvar_ids` column. |
| Class, prediction, gene-symbol and other ClinVar fields | Keep text | These are labels or identifiers, not scalar scores. |

The precision recommendations were checked on the complete relevant cache
columns for chromosomes 22, 1 and 7; ClinVar IDs were checked on all contigs.
They apply to the source versions pinned in these manifests. `Float64` is
also a useful default for downstream calculations on scalar scores; widening
an existing Float32 value cannot recover digits already lost to rounding.

### Scalar and per-consequence columns

Select the required columns, then add numeric companions. This works with
the `LazyFrame` returned by vepyr, or with a collected `DataFrame` using the
same `with_columns()` expressions:

```python
import polars as pl
import vepyr

lf = vepyr.annotate(
    "input.vcf.gz",
    "/data/116_GRCh38_merged",
    reference_fasta="GRCh38.fa",
    plugin_cache_root="/data/plugin_cache",
    plugins=["cadd", "spliceai", "dbnsfp", "clinvar"],
).select(
    "chrom", "start", "Consequence", "CADD_RAW", "CADD_PHRED",
    "SpliceAI_pred_DS_AG", "MetaSVM_score", "CADD_raw", "SIFT4G_score", "ClinVar",
)

numeric = lf.with_columns(
    pl.col("CADD_RAW").cast(pl.Float64).alias("CADD_RAW_num"),
    pl.col("CADD_PHRED").cast(pl.Float32).alias("CADD_PHRED_num"),
    pl.col("SpliceAI_pred_DS_AG")
      .cast(pl.List(pl.Float32)).alias("SpliceAI_pred_DS_AG_num"),
    pl.col("MetaSVM_score")
      .cast(pl.List(pl.Float32)).alias("MetaSVM_score_num"),
    pl.col("CADD_raw").cast(pl.List(pl.Float64)).alias("CADD_raw_num"),
    pl.col("ClinVar").cast(pl.Int32).alias("ClinVar_id"),
)
df = numeric.collect()
```

The list casts preserve consequence order and null elements; they do not
reduce multiple consequences to one score. Polars
[`cast()`](https://docs.pola.rs/api/python/stable/reference/expressions/api/polars.Expr.cast.html)
is strict by default, so unexpected non-numeric values raise an error.
`strict=False` would turn those values into nulls, including a whole
multi-score string such as `0.01&0.03`. Use it only when that loss is intended.
For independently imported text with explicit missing markers, replace `""`
and `"."` with null before casting; existing nulls already remain null.

### Multiple scores inside a dbNSFP element

There are two list levels: the outer list follows vepyr consequences; the
inner list holds the source scores carried by that consequence. In named
DataFrame columns, CSQ escaping represents the inner separator as
`&`, for example `"0.01&.&0.03"`. In a plugin cache Parquet column the same
value is `"0.01,.,0.03"`; the raw dbNSFP source uses `;`.

Use nested
[`list.eval()`](https://docs.pola.rs/api/python/stable/reference/expressions/api/polars.Expr.list.eval.html)
to preserve both levels and each missing score's position:

```python
with_sift4g = numeric.with_columns(
    pl.col("SIFT4G_score")
      .list.eval(
          pl.element().str.split("&").list.eval(
              pl.element().replace(["", "."], None).cast(pl.Float64)
          )
      )
      .alias("SIFT4G_score_values")
)
# List(List(Float64)): ["0.01&.&0.03", None] -> [[0.01, None, 0.03], None]
```

If reading the cache Parquet directly, there is no outer consequence list:
apply `pl.col("sift4g_score").str.split(",").list.eval(...)` with the same
inner missing-value replacement and cast. Avoid automatically taking the
first score or dropping nulls: source score lists are not a mapping to the
outer consequence list. Choose an aggregation only when it suits your
analysis. Changing these analytical columns does not change the cache or
vepyr's separate `output_vcf` serialization path.

## Layout

```
plugins/<name>/<name>.source.toml   one manifest per plugin
scripts/validate_manifests.py       static contract check (runs in CI)
scripts/next_version.sh             next semver tag (used by the release workflow)
scripts/release_notes.sh            manifest-change summary for release notes
.github/workflows/release.yml       manual tag + GitHub release (patch/minor/major)
.claude/skills/adding-a-plugin/     end-to-end guide for authoring a new plugin
```

## Adding a plugin

1. Read the manifest closest in shape to your source (native TSV → `cadd`,
   INFO-packed VCF → `spliceai`/`clinvar`, per-transcript amino-acid match →
   `alphamissense`/`dbnsfp`) and copy its shape.
2. Write `plugins/<name>/<name>.source.toml` against the
   [manifest reference](https://biodatageeks.org/vepyr/plugins/#manifest-structure).
   Note the TOML ordering rule: top-level scalars (`plugin_name`,
   `coordinate_system`, `ingest_sql`) must precede any table header.
3. Record provenance: every `[[source]]` carries the upstream `url` the raw
   file was downloaded from and the `md5` of that file (the publisher's
   checksum where one exists), so a cache can be traced back to exact input
   bytes. Never point `url` at a mirror or a Drive share. `md5` always
   describes the file at `url`. If the build input is a derived artifact of
   that file (AlphaMissense's BGZF+tabix re-compression of the upstream plain
   gzip), describe the preprocessing in a `README.md` next to the manifest
   (see [`plugins/alphamissense`](plugins/alphamissense/README.md)) rather
   than adding a second digest to the manifest, and build with
   `verify_source` disabled since the derived file cannot match `md5`.
4. Validate: `python scripts/validate_manifests.py` — checks plugin/filename
   agreement, providers, coordinate system, tabix/compression pairing, source
   `url`/`md5` presence and shape, value and match column uniqueness, CSQ
   field names, `allele_match` and `field_order`, and rejects the retired
   `csq_rank` key. CI runs it on every pull request.
5. Build one chromosome with `build_plugin_cache` and compare the resulting CSQ
   fields against an Ensembl VEP 116 run before opening a PR.

The [`adding-a-plugin`](.claude/skills/adding-a-plugin/SKILL.md) skill walks
through source triage, memory-safe per-chromosome building and the parity gate
in detail.

## Releasing

Releases are git tags of the form `vMAJOR.MINOR.PATCH`; users pin one as the
`version` of `build_plugin_cache()`. To cut a release run the **Release**
workflow from the Actions tab (or `gh workflow run release.yml -f bump=minor`)
on `master` and pick the bump:

| Bump | When |
|---|---|
| `patch` | A manifest fix that keeps the emitted CSQ fields and values as they were (URL/md5 update, typo, validator-only change). |
| `minor` | A new plugin, or new CSQ fields added to an existing one. |
| `major` | A manifest change that alters or removes existing CSQ fields, or a change in what `version` resolves to. |

The workflow re-runs the manifest validator, bumps the latest tag with
`scripts/next_version.sh`, then creates the tag and the GitHub release in one
step; the notes list the manifests changed since the previous tag followed by
the merged pull requests. Tick **dry_run** to preview the version and notes in the job
summary without tagging. `scripts/next_version.sh minor` prints the next tag
locally.

## Licence

Manifests: [Apache-2.0](LICENSE). The annotation sources they describe carry
their **own** terms — CADD, SpliceAI and AlphaMissense restrict use to academic
/ non-profit research, dbNSFP forbids redistribution of derived copies, ClinVar
is public domain. Check the upstream terms before building or sharing a cache.

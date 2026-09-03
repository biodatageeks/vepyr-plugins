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
   bytes. Never point `url` at a mirror or a Drive share. `build_plugin_cache`
   hashes the file `source_path` resolves to and refuses a mismatch, so when
   the build input is a derived artifact of `url` (AlphaMissense's BGZF+tabix
   re-compression of the upstream plain gzip) also declare `path_md5`, the
   digest of that artifact.
4. Validate: `python scripts/validate_manifests.py` — checks plugin/filename
   agreement, providers, coordinate system, tabix/compression pairing, source
   `url`/`md5` presence and shape (and `path_md5` shape), value and match column uniqueness, CSQ
   field names, `allele_match` and `field_order`.
   CI runs it on every pull request.
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

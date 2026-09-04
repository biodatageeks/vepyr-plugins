# AlphaMissense — preparing the build input

The manifest's `url` / `md5` describe the file DeepMind publishes:

| | |
|---|---|
| URL | <https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz> |
| MD5 | `9fd167735f16a1b87da6eb3e4c25fcb5` (from the GCS object metadata) |

That file is a **plain gzip** stream. It is not seekable, so vepyr's cache
builder, which reads one chromosome at a time through a tabix index, cannot
consume it directly. The build input the manifest points at
(`AlphaMissense_hg38.bgz.tsv.gz`) is the same data re-compressed as BGZF with
a tabix index alongside. It has to be produced locally, once, before
`build_plugin_cache()`.

## Steps

Requires `bgzip` and `tabix` from [htslib](https://www.htslib.org/).

```bash
# 1. Download the upstream file and confirm it is the published bytes.
curl -fL -o AlphaMissense_hg38.tsv.gz \
  https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
echo "9fd167735f16a1b87da6eb3e4c25fcb5  AlphaMissense_hg38.tsv.gz" | md5sum -c

# 2. Re-compress as BGZF. The decompressed bytes are unchanged; only the
#    container differs.
gzip -dc AlphaMissense_hg38.tsv.gz | bgzip -c > AlphaMissense_hg38.bgz.tsv.gz

# 3. Index. Column 1 is the contig (chr-prefixed), column 2 the 1-based
#    position, used for both begin and end; header/licence lines start
#    with '#'.
tabix -s 1 -b 2 -e 2 -c '#' AlphaMissense_hg38.bgz.tsv.gz
```

The result is `AlphaMissense_hg38.bgz.tsv.gz` plus
`AlphaMissense_hg38.bgz.tsv.gz.tbi`. Pass the `.bgz.tsv.gz` path as
`source_path`; the builder finds the `.tbi` beside it.

For reference, the file produced this way with htslib's default compression
level hashed to `46d0028375cf95088bd014ff6855cffd`. BGZF output is not
guaranteed byte-identical across htslib versions or thread counts, so treat
that as a sanity check, not a requirement.

## Building

`build_plugin_cache()` hashes the file `source_path` resolves to and compares
it with the manifest's `md5`. The BGZF file can never match the upstream
digest, so disable that check for this plugin:

```python
vepyr.build_plugin_cache(
    plugin="alphamissense",
    version="<vepyr-plugins tag>",
    source_path="AlphaMissense_hg38.bgz.tsv.gz",
    cache_dir="/data/116_GRCh38_merged",
    plugin_cache_root="/data/plugin_cache",
    verify_source=False,   # build input is a BGZF re-compression of `url`
)
```

Step 1 above is where the upstream bytes are verified instead.

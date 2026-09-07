#!/usr/bin/env python3
"""Profile numeric plugin text using full Parquet columns and Rust formatting.

Requires polars and rustc. Writes only to --out; existing caches are read-only.
Example: python scripts/audit_numeric_types.py --cache-root /data/cache --out /data/audit
Use --sample-rows 100000 for an explicitly sampled preliminary assessment.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import time

import polars as pl


SCALES = {
    "cadd": {"cadd_raw": 6, "cadd_phred": 3},
    "spliceai": {f"ds_{x}": 2 for x in ("ag", "al", "dg", "dl")},
    "clinvar": {"id": 0},
    "dbnsfp": {
        "sift4g_score": 3,
        "polyphen2_hdiv_score": 3,
        "polyphen2_hvar_score": 3,
        "mutationtaster_score": 6,
        "provean_score": 3,
        "vest4_score": 3,
        "metasvm_score": 4,
        "metalr_score": 4,
        "revel_score": 3,
        "gerp_rs": 3,
        "phylop100way_vertebrate": 6,
        "phastcons100way_vertebrate": 6,
        "cadd_raw": 6,
        "cadd_phred": 3,
    },
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--chroms", nargs="+", default=["22", "1", "7"])
    parser.add_argument("--plugins", nargs="+", default=list(SCALES))
    parser.add_argument("--sample-rows", type=int)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    binary = args.out / "numeric_format_probe"
    subprocess.run(
        [
            "rustc",
            "--edition=2021",
            "-O",
            str(Path(__file__).with_name("numeric_format_probe.rs")),
            "-o",
            str(binary),
        ],
        check=True,
    )
    for plugin in args.plugins:
        for chrom in args.chroms:
            dest = args.out / f"{plugin}_chr{chrom}.json"
            path = args.cache_root / "plugin" / plugin / f"chr{chrom}.parquet"
            scan = pl.scan_parquet(path)
            if args.sample_rows:
                scan = scan.head(args.sample_rows)
            result = {
                "source": str(path),
                "sample_rows": args.sample_rows,
                "columns": {},
            }
            for column, scale in SCALES[plugin].items():
                start = time.monotonic()
                values = (
                    scan.group_by(column)
                    .len(name="count")
                    .collect(engine="streaming")
                    .sort(column)
                )
                nulls = values.filter(pl.col(column).is_null())["count"].sum()
                values = values.filter(pl.col(column).is_not_null())
                data = "".join(f"{s}\t{n}\n" for s, n in values.iter_rows())
                probe = subprocess.run(
                    [str(binary), str(scale)],
                    input=data,
                    text=True,
                    capture_output=True,
                    check=True,
                )
                lines = probe.stdout.splitlines()
                keys = [
                    "non_null",
                    "invalid_numeric",
                    "f32_display_mismatch",
                    "f64_display_mismatch",
                    "f32_fixed_mismatch",
                    "f64_fixed_mismatch",
                ]
                stats = dict(zip(keys, map(int, lines[0].split("\t"))))
                stats.update(nulls=nulls, distinct=values.height, fixed_scale=scale)
                stats["examples"] = {k: [] for k in keys[1:]}
                for line in lines[1:]:
                    i, value = line.split("\t", 1)
                    stats["examples"][keys[int(i) + 1]].append(value)
                stats["decimal_places"] = (
                    values.with_columns(
                        pl.col(column)
                        .str.extract(r"^-?\d+\.(\d+)$", 1)
                        .str.len_chars()
                        .alias("scale")
                    )
                    .group_by("scale")
                    .agg(pl.col("count").sum())
                    .sort("scale")
                    .to_dicts()
                )
                numeric = values.with_columns(
                    pl.col(column).cast(pl.Float64, strict=False).alias("number")
                ).filter(pl.col("number").is_not_null())
                stats["min"] = numeric["number"].min()
                stats["max"] = numeric["number"].max()
                collisions = (
                    numeric.group_by("number")
                    .agg(pl.col(column).alias("spellings"), pl.col("count").sum())
                    .filter(pl.col("spellings").list.len() > 1)
                    .sort("number")
                )
                stats["numeric_spelling_collisions"] = collisions.height
                stats["collision_examples"] = collisions.head(5).to_dicts()
                result["columns"][column] = stats
                print(
                    f'{plugin}/chr{chrom}/{column}: rows={stats["non_null"] + nulls:,} invalid={stats["invalid_numeric"]:,} f32_diff={stats["f32_display_mismatch"]:,} f64_fixed_diff={stats["f64_fixed_mismatch"]:,} collisions={collisions.height} seconds={time.monotonic()-start:.1f}',
                    flush=True,
                )
            dest.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()

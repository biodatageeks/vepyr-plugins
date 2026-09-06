#!/usr/bin/env python3
"""Test numeric storage followed by explicit SQL rendering on real cache values.

The rendered expression is Utf8. Deploying this after a numeric plugin lookup
requires an engine rendering hook; ingest_sql alone cannot provide that hook.
"""

import argparse
import gc
import json
from pathlib import Path
import subprocess

from datafusion import SessionContext
import polars as pl

from audit_numeric_types import SCALES


def fixed(value, scale, preserve_sign=True):
    decimal = f"CAST(CAST({value} AS DECIMAL(18,{scale})) AS VARCHAR)"
    if not preserve_sign:
        return decimal
    # Arrow/DataFusion floating equality distinguishes -0.0 and +0.0 here.
    # Decimal equality detects the zero after removing its sign; take the sign
    # from the numeric value's string rendering, not from the source column.
    return f"CASE WHEN CAST({value} AS DECIMAL(18,{scale})) = 0 AND starts_with(CAST({value} AS VARCHAR), '-') THEN concat('-', {decimal}) ELSE {decimal} END"


def native_probe(out, cache_dir):
    """Run representative casts through vepyr's own (currently DF53) builder."""
    import vepyr

    cases = [
        (
            "raw",
            6,
            [
                "-0.000000",
                "0.000000",
                "-0.000001",
                "0.163600",
                "16.042138",
                "-21.700286",
            ],
        ),
        ("ds", 2, ["-0.00", "0.00", "0.10", "0.29", "1.00"]),
        ("meta", 4, ["-0.0000", "0.0000", "0.0010", "-0.1234"]),
    ]
    results = []
    for name, scale, values in cases:
        for dtype in ["FLOAT", "DOUBLE"]:
            directory = out / f"{name}_{dtype}"
            directory.mkdir(parents=True, exist_ok=True)
            source = directory / "input.parquet"
            pl.DataFrame(
                {
                    "chrom": ["22"] * len(values),
                    "start": list(range(1, len(values) + 1)),
                    "end": list(range(1, len(values) + 1)),
                    "allele_string": ["A/T"] * len(values),
                    "original": values,
                }
            ).write_parquet(source)
            sql = f'SELECT chrom, start, "end", allele_string, {fixed("v", scale)} AS score FROM (SELECT *, CAST(original AS {dtype}) AS v FROM plugin_numeric_sql_probe_src)'
            manifest = directory / "probe.source.toml"
            manifest.write_text(
                'plugin_name = "numeric_sql_probe"\ncoordinate_system = "1-based"\ningest_sql = '
                + json.dumps(sql)
                + '\n[[source]]\nprovider = "parquet"\npath = "input.parquet"\n[[value_columns]]\ncolumn = "score"\ncsq_field = "score"\ntype = "Utf8"\n'
            )
            vepyr._core.build_plugin_cache(
                str(manifest),
                str(source),
                str(cache_dir),
                str(directory / "cache"),
                ["22"],
                True,
                "skip",
                "numeric-sql-probe",
            )
            actual = (
                pl.read_parquet(
                    directory / "cache/plugin/numeric_sql_probe/chr22.parquet"
                )
                .sort("start")["score"]
                .to_list()
            )
            results.append(
                {
                    "case": name,
                    "numeric_type": dtype,
                    "source": values,
                    "rendered": actual,
                    "mismatches": sum(a != b for a, b in zip(values, actual)),
                    "sql": sql,
                }
            )
    (out / "native_sql.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--chroms", nargs="+", default=["22", "1", "7"])
    p.add_argument("--native-cache-dir", type=Path)
    p.add_argument("--native-only", action="store_true")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.native_only:
        if not args.native_cache_dir:
            p.error("--native-only requires --native-cache-dir")
        native_probe(args.out, args.native_cache_dir)
        return
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
    reports = []
    for plugin, columns in SCALES.items():
        if plugin == "clinvar":
            continue
        for column, scale in columns.items():
            if plugin == "dbnsfp" and column not in [
                "metasvm_score",
                "metalr_score",
                "gerp_rs",
                "phylop100way_vertebrate",
                "phastcons100way_vertebrate",
                "cadd_raw",
                "cadd_phred",
            ]:
                continue
            print(f"START {plugin}/{column}", flush=True)
            paths = [
                args.cache_root / "plugin" / plugin / f"chr{chrom}.parquet"
                for chrom in args.chroms
            ]
            values = (
                pl.scan_parquet(paths)
                .group_by(column)
                .agg(pl.len().cast(pl.UInt64).alias("count"))
                .collect(engine="streaming")
                .rename({column: "original"})
                .filter(pl.col("original").is_not_null())
                .sort("original")
            )
            policy = (
                100 if column == "cadd_phred" else 101 if column == "gerp_rs" else scale
            )
            data = "".join(f"{s}\t{n}\n" for s, n in values.iter_rows())
            run = subprocess.run(
                [str(binary), str(policy)],
                input=data,
                text=True,
                capture_output=True,
                check=True,
            )
            lines = run.stdout.splitlines()
            counts = list(map(int, lines[0].split("\t")))
            report = {
                "plugin": plugin,
                "column": column,
                "rows": counts[0],
                "policy": policy,
                "f32_render_mismatch": counts[4],
                "f64_render_mismatch": counts[5],
                "sql": [],
            }
            if column != "gerp_rs":
                context = SessionContext()
                context.register_record_batches(
                    "source", [values.to_arrow().to_batches()]
                )
                # Float32 for existing-supported scores, Float64 also for CADD raw.
                for dtype in ["FLOAT", "DOUBLE"] if column == "cadd_raw" else ["FLOAT"]:
                    for preserve_sign in (
                        [False, True] if column != "cadd_phred" else [True]
                    ):
                        rendered = fixed("v", scale, preserve_sign)
                        if column == "cadd_phred":
                            rendered = (
                                "CASE "
                                + " ".join(
                                    f'WHEN v < {bound} THEN {fixed("v", sc)}'
                                    for bound, sc in [(10, 3), (20, 2), (30, 1)]
                                )
                                + f' ELSE {fixed("v", 0)} END'
                            )
                        expression = f"SELECT original, count, {rendered} AS rendered FROM (SELECT *, CAST(original AS {dtype}) AS v FROM source)"
                        query = f"SELECT SUM(CASE WHEN original IS DISTINCT FROM rendered THEN CAST(count AS BIGINT) ELSE CAST(0 AS BIGINT) END) AS mismatches FROM ({expression})"
                        result = context.sql(query).to_pydict()["mismatches"][0]
                        report["sql"].append(
                            {
                                "numeric_type": dtype,
                                "preserve_negative_zero": preserve_sign,
                                "mismatches": result,
                                "expression": expression,
                            }
                        )
            reports.append(report)
            (args.out / "sql_rendering.json").write_text(
                json.dumps(reports, indent=2) + "\n"
            )
            print(
                f'{plugin}/{column}: rows={counts[0]:,} custom_f32_diff={counts[4]:,} custom_f64_diff={counts[5]:,} SQL={[(x["numeric_type"], x["preserve_negative_zero"], x["mismatches"]) for x in report["sql"]]}',
                flush=True,
            )
            if column != "gerp_rs":
                del context
            del values, data
            gc.collect()
    if args.native_cache_dir:
        native_probe(args.out, args.native_cache_dir)


if __name__ == "__main__":
    main()

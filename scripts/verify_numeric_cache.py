#!/usr/bin/env python3
"""Rebuild experimental caches from existing Parquet and compare strict VCF bodies.

Uses the installed vepyr native builder, including its SQL, tiering and writer.
The --data-root caches and --vepyr-repo references are read-only. All generated
manifests, caches, VCFs and results go under --out. This is cache re-ingestion,
not a fresh download/parse of the original upstream VCF/TSV sources.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import tomllib

import pyarrow.parquet as pq
import vepyr

PLUGINS = ["spliceai", "cadd", "alphamissense", "dbnsfp", "clinvar"]
SCALAR_DBNSFP = [
    "metasvm_score",
    "metalr_score",
    "gerp_rs",
    "phylop100way_vertebrate",
    "phastcons100way_vertebrate",
    "cadd_raw",
    "cadd_phred",
]


def body_digest(path):
    digest = hashlib.md5()
    count = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as handle:
        for line in handle:
            if not line.startswith(b"#"):
                digest.update(line)
                count += 1
    return {"records": count, "body_md5": digest.hexdigest()}


def make_manifest(plugin, changes, out):
    repo = Path(__file__).resolve().parents[1]
    original = tomllib.loads(
        (repo / "plugins" / plugin / f"{plugin}.source.toml").read_text()
    )
    keys = ["chrom", "start", "end", "allele_string"]
    keys += [m["column"] for m in original.get("match_column", [])]
    keys += [v["column"] for v in original["value_columns"]]
    select = [
        f'CAST("{k}" AS {changes[k][1]}) AS "{k}"' if k in changes else f'"{k}"'
        for k in keys
    ]
    lines = [
        f'plugin_name = "{plugin}"',
        'coordinate_system = "1-based"',
        f'field_order = "{original.get("field_order", "declared")}"',
        f'allele_match = "{original.get("allele_match", "exact")}"',
        'ingest_sql = """SELECT ' + ", ".join(select) + f' FROM plugin_{plugin}_src"""',
        "[[source]]",
        'provider = "parquet"',
        'path = "input.parquet"',
    ]
    for match in original.get("match_column", []):
        lines += ["[[match_column]]"] + [
            f"{k} = {json.dumps(v)}" for k, v in match.items()
        ]
    for value in original["value_columns"]:
        value = dict(value)
        if value["column"] in changes:
            value["type"] = changes[value["column"]][0]
        lines += ["[[value_columns]]"] + [
            f"{k} = {json.dumps(v)}" for k, v in value.items()
        ]
    path = out / f"{plugin}.source.toml"
    path.write_text("\n".join(lines) + "\n")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--vepyr-repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--chroms", nargs="+", default=["22", "1", "7"])
    parser.add_argument(
        "--mode",
        choices=["baseline", "clinvar_int32", "dbnsfp_float32", "dbnsfp_formatted"],
        required=True,
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    existing = args.data_root / "plugin_cache_116"
    core_cache = args.data_root / "cache" / "116_GRCh38_merged"
    cache_root = existing if args.mode == "baseline" else args.out / "cache"
    changed_plugin = {
        "clinvar_int32": "clinvar",
        "dbnsfp_float32": "dbnsfp",
        "dbnsfp_formatted": "dbnsfp",
    }.get(args.mode)
    builds = []
    if changed_plugin:
        (cache_root / "plugin").mkdir(parents=True, exist_ok=True)
        for plugin in PLUGINS:
            link = cache_root / "plugin" / plugin
            if plugin != changed_plugin and not link.exists():
                link.symlink_to(existing / "plugin" / plugin, target_is_directory=True)
        changes = (
            {"id": ("Int32", "INT")}
            if changed_plugin == "clinvar"
            else {c: ("Float32", "FLOAT") for c in SCALAR_DBNSFP}
        )
        if args.mode == "dbnsfp_formatted":
            # Full-cache scanning demonstrates that CADD raw needs Float64.
            # Float64 is unsupported by the current builder/reader, so retain it.
            changes.pop("cadd_raw")
        manifest = make_manifest(changed_plugin, changes, args.out)
        for chrom in args.chroms:
            start = time.monotonic()
            result, provenance = vepyr._core.build_plugin_cache(
                str(manifest),
                str(existing / "plugin" / changed_plugin / f"chr{chrom}.parquet"),
                str(core_cache),
                str(cache_root),
                [chrom],
                True,
                "skip",
                f"numeric-type-audit-20260906-{args.mode}",
            )
            shard = cache_root / "plugin" / changed_plugin / f"chr{chrom}.parquet"
            schema = pq.read_schema(shard)
            for column, (dtype, _) in changes.items():
                expected = "int32" if dtype == "Int32" else "float"
                assert str(schema.field(column).type) == expected, schema
            item = {
                "chrom": chrom,
                "result": result,
                "seconds": time.monotonic() - start,
                "bytes": shard.stat().st_size,
                "schema": str(schema),
                "sources": json.loads(provenance),
            }
            builds.append(item)
            (args.out / "builds.json").write_text(json.dumps(builds, indent=2) + "\n")
            print(
                f'BUILT {args.mode} chr{chrom} {result} seconds={item["seconds"]:.1f}',
                flush=True,
            )
    comparisons = []
    for chrom in args.chroms:
        run = args.vepyr_repo / "e2e-testing" / "results" / "116" / f"fast_chr{chrom}"
        output = args.out / f"vepyr_chr{chrom}.vcf.gz"
        start = time.monotonic()
        vepyr.annotate(
            str(run / f"input_chr{chrom}.vcf.gz"),
            str(core_cache),
            everything=True,
            reference_fasta=str(
                args.data_root / "input" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
            ),
            output_vcf=str(output),
            workers=4,
            plugin_cache_root=str(cache_root),
            plugins=PLUGINS,
            show_progress=False,
        )
        print(
            f"ANNOTATED {args.mode} chr{chrom} seconds={time.monotonic()-start:.1f}",
            flush=True,
        )
        reference = run / f"vep_chr{chrom}_merged_plugins.vcf"
        original = run / f"vepyr_parquet_chr{chrom}_merged_plugins.vcf.gz"
        item = {
            "chrom": chrom,
            "output": str(output),
            "reference": str(reference),
            "actual": body_digest(output),
            "vep": body_digest(reference),
            "previous_vepyr": body_digest(original),
        }
        item["matches_vep"] = item["actual"] == item["vep"]
        item["matches_previous_vepyr"] = item["actual"] == item["previous_vepyr"]
        if args.mode == "dbnsfp_formatted":
            from render_numeric_csq import render_vcf

            rendered_path = args.out / f"vepyr_chr{chrom}_formatted.vcf.gz"
            rendered, changed = render_vcf(output, rendered_path)
            item.update(
                experimental_formatted=rendered,
                formatted_output=str(rendered_path),
                formatted_changed_records=changed,
                formatted_matches_vep=rendered == item["vep"],
            )
        comparisons.append(item)
        (args.out / "comparisons.json").write_text(
            json.dumps(comparisons, indent=2) + "\n"
        )
        print(json.dumps(item), flush=True)
        if args.mode in ["baseline", "clinvar_int32"] and not item["matches_vep"]:
            raise RuntimeError(
                f"{args.mode} chr{chrom}: unexpected native body mismatch"
            )
        if args.mode == "dbnsfp_formatted" and not item["formatted_matches_vep"]:
            raise RuntimeError(
                f"chr{chrom}: experimental formatter did not restore the VEP body"
            )


if __name__ == "__main__":
    main()

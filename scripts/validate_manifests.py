#!/usr/bin/env python3
"""Validate the static contract of every plugin source manifest."""

from __future__ import annotations

import re
import sys
import tomllib
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_GLOB = "plugins/*/*.source.toml"
COORDINATE_SYSTEMS = {"1-based", "0-based-half-open"}
PROVIDERS = {"vcf", "csv", "tsv", "parquet", "bed"}
VALUE_TYPES = {"Utf8", "Float32", "Int32"}
ALLELE_MATCHES = {"exact", "minimised"}
FIELD_ORDERS = {"declared", "alphabetical"}
CSQ_FIELD = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.+-]*$")
RESERVED_VCF_KEYS = {
    "fileformat",
    "fileDate",
    "source",
    "reference",
    "phasing",
    "assembly",
    "pedigreeDB",
    "INFO",
    "FILTER",
    "FORMAT",
    "ALT",
    "contig",
    "SAMPLE",
    "PEDIGREE",
    "META",
    "VEP",
    "VEP-command-line",
    "datafusion-bio-function-vep",
    "datafusion-bio-function-vep-command-line",
}


def require_string(
    value: Any, path: Path, field: str, errors: list[str]
) -> str | None:
    if not isinstance(value, str) or not value:
        errors.append(f"{path}: {field} must be a non-empty string")
        return None
    return value


def duplicate_values(values: list[str]) -> set[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        else:
            seen.add(value)
    return duplicates


def validate_manifest(path: Path, errors: list[str]) -> None:
    try:
        manifest = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as exc:
        errors.append(f"{path}: cannot parse TOML: {exc}")
        return

    plugin_name = require_string(manifest.get("plugin_name"), path, "plugin_name", errors)
    expected_name = path.parent.name
    if plugin_name is not None and plugin_name != expected_name:
        errors.append(
            f"{path}: plugin_name {plugin_name!r} must match directory {expected_name!r}"
        )
    if path.name != f"{expected_name}.source.toml":
        errors.append(
            f"{path}: filename must be {expected_name}.source.toml"
        )

    coordinate_system = manifest.get("coordinate_system")
    if coordinate_system not in COORDINATE_SYSTEMS:
        errors.append(
            f"{path}: coordinate_system must be one of {sorted(COORDINATE_SYSTEMS)}"
        )
    require_string(manifest.get("ingest_sql"), path, "ingest_sql", errors)

    allele_match = manifest.get("allele_match", "exact")
    if allele_match not in ALLELE_MATCHES:
        errors.append(f"{path}: allele_match must be one of {sorted(ALLELE_MATCHES)}")
    field_order = manifest.get("field_order", "declared")
    if field_order not in FIELD_ORDERS:
        errors.append(f"{path}: field_order must be one of {sorted(FIELD_ORDERS)}")
    csq_rank = manifest.get("csq_rank")
    if csq_rank is not None and (
        isinstance(csq_rank, bool) or not isinstance(csq_rank, int) or csq_rank < 0
    ):
        errors.append(f"{path}: csq_rank must be a non-negative integer")
    if not isinstance(manifest.get("assume_unique", False), bool):
        errors.append(f"{path}: assume_unique must be a boolean")

    sources = manifest.get("source")
    if not isinstance(sources, list) or not sources:
        errors.append(f"{path}: source must contain at least one [[source]] table")
    else:
        table_parts: list[str] = []
        for index, source in enumerate(sources):
            label = f"source[{index}]"
            if not isinstance(source, dict):
                errors.append(f"{path}: {label} must be a table")
                continue
            provider = source.get("provider")
            if provider not in PROVIDERS:
                errors.append(f"{path}: {label}.provider must be one of {sorted(PROVIDERS)}")
            require_string(source.get("path"), path, f"{label}.path", errors)
            part = source.get("part", "")
            if not isinstance(part, str):
                errors.append(f"{path}: {label}.part must be a string")
            else:
                table_parts.append(part)
            record_layout = source.get("record_layout", False)
            if not isinstance(record_layout, bool):
                errors.append(f"{path}: {label}.record_layout must be a boolean")
            elif record_layout and provider != "vcf":
                errors.append(f"{path}: {label}.record_layout is supported only for VCF")
            csv = source.get("csv")
            if provider in {"csv", "tsv"} and not isinstance(csv, dict):
                errors.append(f"{path}: {label}.csv must be a table for {provider}")
            elif isinstance(csv, dict):
                text_index = csv.get("index")
                if text_index not in {None, "tabix"}:
                    errors.append(f"{path}: {label}.csv.index must be 'tabix'")
                if text_index == "tabix":
                    if provider not in {"csv", "tsv"}:
                        errors.append(
                            f"{path}: {label}.csv.index is supported only for csv/tsv"
                        )
                    if csv.get("compression") != "gzip":
                        errors.append(
                            f"{path}: {label}.csv.index='tabix' requires "
                            "compression='gzip' (BGZF)"
                        )
        for part in sorted(duplicate_values(table_parts)):
            errors.append(f"{path}: duplicate source part {part!r} creates a table collision")

    value_column_names: list[str] = []
    value_columns = manifest.get("value_columns")
    if not isinstance(value_columns, list) or not value_columns:
        errors.append(f"{path}: value_columns must be a non-empty array")
    else:
        csq_fields: list[str] = []
        for index, value in enumerate(value_columns):
            label = f"value_columns[{index}]"
            if not isinstance(value, dict):
                errors.append(f"{path}: {label} must be a table")
                continue
            column = require_string(value.get("column"), path, f"{label}.column", errors)
            csq_field = require_string(
                value.get("csq_field"), path, f"{label}.csq_field", errors
            )
            if column is not None:
                value_column_names.append(column)
            if csq_field is not None:
                csq_fields.append(csq_field)
                if not CSQ_FIELD.fullmatch(csq_field):
                    errors.append(f"{path}: invalid CSQ field name {csq_field!r}")
                if csq_field in RESERVED_VCF_KEYS:
                    errors.append(f"{path}: reserved VCF key used as CSQ field {csq_field!r}")
            if value.get("type") not in VALUE_TYPES:
                errors.append(f"{path}: {label}.type must be one of {sorted(VALUE_TYPES)}")
        for column in sorted(duplicate_values(value_column_names)):
            errors.append(f"{path}: duplicate value column {column!r}")
        for field in sorted(duplicate_values(csq_fields)):
            errors.append(f"{path}: duplicate CSQ field {field!r}")

    match_column_names: list[str] = []
    match_columns = manifest.get("match_column", [])
    if not isinstance(match_columns, list):
        errors.append(f"{path}: match_column must be an array of tables")
    else:
        for index, match in enumerate(match_columns):
            label = f"match_column[{index}]"
            if not isinstance(match, dict):
                errors.append(f"{path}: {label} must be a table")
                continue
            column = require_string(match.get("column"), path, f"{label}.column", errors)
            require_string(match.get("template"), path, f"{label}.template", errors)
            if column is not None:
                match_column_names.append(column)
        for column in sorted(duplicate_values(match_column_names)):
            errors.append(f"{path}: duplicate match column {column!r}")

    for column in sorted(set(value_column_names) & set(match_column_names)):
        errors.append(f"{path}: column {column!r} cannot be both a value and match column")


def main() -> int:
    paths = sorted(ROOT.glob(MANIFEST_GLOB))
    if not paths:
        print(f"error: no manifests found at {MANIFEST_GLOB}", file=sys.stderr)
        return 1

    errors: list[str] = []
    for path in paths:
        validate_manifest(path, errors)

    if errors:
        print("Manifest validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print(f"Validated {len(paths)} plugin source manifests")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

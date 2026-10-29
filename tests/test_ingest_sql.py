"""Run each manifest's ingest_sql on small fixture tables with DataFusion.

The fixtures mimic the source tables the cache builder registers (CSV columns from
the manifest schema; VCF tables with `alt` joined by '|', list-typed multi-value INFO
fields and `_vcf_info_keys`). Each test checks the normalized key (start,
allele_string) and the value columns for rows that must be kept as they are, rows
that must be rewritten and rows that must yield no value.
"""

import tomllib
from pathlib import Path

import pyarrow as pa
import pytest

datafusion = pytest.importorskip("datafusion")

ROOT = Path(__file__).resolve().parents[1]


def manifest(name: str) -> dict:
    return tomllib.loads((ROOT / "plugins" / name / f"{name}.source.toml").read_text(encoding="utf-8"))


def run(sql: str, tables: dict[str, pa.Table]) -> list[dict]:
    ctx = datafusion.SessionContext()
    for name, table in tables.items():
        batches = table.to_batches() or [pa.RecordBatch.from_pylist([], schema=table.schema)]
        ctx.register_record_batches(name, [batches])
    rows = ctx.sql(sql).to_arrow_table().to_pylist()
    return sorted(rows, key=lambda r: (r["start"], r["allele_string"]))


def utf8_table(columns: list[str], rows: list[dict]) -> pa.Table:
    return pa.table({c: pa.array([r.get(c) for r in rows], pa.string()) for c in columns})


# --- CADD -------------------------------------------------------------------

CADD_COLS = ["chrom", "pos", "ref", "alt", "raw_score", "phred"]


def cadd_rows(snv: list[tuple], indel: list[tuple]) -> list[dict]:
    m = manifest("cadd")
    t = lambda rows: utf8_table(CADD_COLS, [dict(zip(CADD_COLS, r)) for r in rows])  # noqa: E731
    return run(m["ingest_sql"], {"plugin_cadd_src_snv": t(snv), "plugin_cadd_src_indel": t(indel)})


def test_cadd_snv_passes_through_unshifted():
    rows = cadd_rows([("22", "100", "A", "G", "0.512", "10.31")], [])
    assert rows == [{"chrom": "22", "start": 100, "end": 100, "allele_string": "A/G",
                     "cadd_raw": "0.512", "cadd_phred": "10.31"}]


@pytest.mark.parametrize("ref, alt, start, allele", [
    ("AG", "A", 101, "G/-"),        # deletion: anchor base dropped, start + 1
    ("A", "AG", 101, "-/G"),        # insertion
    ("AGT", "AC", 101, "GT/C"),     # delins sharing the anchor base
])
def test_cadd_indel_anchor_trim(ref, alt, start, allele):
    rows = cadd_rows([], [("22", "100", ref, alt, "1.0", "15.0")])
    assert [(r["start"], r["end"], r["allele_string"]) for r in rows] == [(start, start, allele)]


def test_cadd_indel_without_shared_anchor_is_not_shifted():
    rows = cadd_rows([], [("22", "100", "AG", "TC", "1.0", "15.0")])
    assert [(r["start"], r["allele_string"]) for r in rows] == [(100, "AG/TC")]


def test_cadd_snv_and_indel_at_same_position_do_not_collide():
    rows = cadd_rows([("22", "100", "A", "G", "0.1", "1.0")], [("22", "100", "AG", "A", "2.0", "20.0")])
    keys = [(r["start"], r["allele_string"]) for r in rows]
    assert keys == [(100, "A/G"), (101, "G/-")]
    assert len(set(keys)) == len(keys)


# --- ClinVar ----------------------------------------------------------------

def clinvar_table(rows: list[dict]) -> pa.Table:
    lst = pa.list_(pa.string())
    return pa.table({
        "chrom": pa.array([r["chrom"] for r in rows], pa.string()),
        "start": pa.array([r["start"] for r in rows], pa.int64()),
        "ref": pa.array([r["ref"] for r in rows], pa.string()),
        "alt": pa.array([r["alt"] for r in rows], pa.string()),
        "id": pa.array([r["id"] for r in rows], pa.string()),
        "CLNSIG": pa.array([r.get("CLNSIG") for r in rows], lst),
        "CLNREVSTAT": pa.array([r.get("CLNREVSTAT") for r in rows], lst),
        "CLNDN": pa.array([r.get("CLNDN") for r in rows], lst),
        "CLNVC": pa.array([r.get("CLNVC") for r in rows], pa.string()),
        "CLNVI": pa.array([r.get("CLNVI") for r in rows], lst),
        "_vcf_info_keys": pa.array([r.get("keys", "") for r in rows], pa.string()),
    })


def clinvar_rows(rows: list[dict]) -> list[dict]:
    return run(manifest("clinvar")["ingest_sql"], {"plugin_clinvar_src": clinvar_table(rows)})


def test_clinvar_values_and_list_join():
    rows = clinvar_rows([{"chrom": "22", "start": 100, "ref": "A", "alt": "G", "id": "12345",
                          "CLNSIG": ["Pathogenic", "Likely_pathogenic"], "CLNREVSTAT": ["criteria_provided"],
                          "CLNDN": ["Disease_A"], "CLNVC": "single_nucleotide_variant", "CLNVI": None,
                          "keys": "CLNSIG;CLNREVSTAT;CLNDN;CLNVC"}])
    (r,) = rows
    assert (r["start"], r["allele_string"], r["id"]) == (100, "A/G", "12345")
    assert r["clnsig"] == "Pathogenic,Likely_pathogenic"
    assert r["clnvc"] == "single_nucleotide_variant"
    assert r["clnvi"] is None  # key absent from the record: no value


def test_clinvar_key_present_without_value_gives_dot():
    (r,) = clinvar_rows([{"chrom": "22", "start": 100, "ref": "A", "alt": "G", "id": "1",
                          "keys": "CLNSIG;CLNVI"}])
    assert r["clnsig"] == "." and r["clnvi"] == "."
    assert r["clnrevstat"] is None


def test_clinvar_multi_alt_is_exploded_and_only_anchor_trimmed():
    # Only the shared anchor base is dropped; AG>AGG stays G/GG and is matched
    # through allele_match = "minimised" at lookup time.
    rows = clinvar_rows([{"chrom": "22", "start": 100, "ref": "AG", "alt": "A|AGG", "id": "7", "keys": ""}])
    assert [(r["start"], r["allele_string"], r["id"]) for r in rows] == [(101, "G/-", "7"), (101, "G/GG", "7")]


# --- SpliceAI ---------------------------------------------------------------

def spliceai_rows(rows: list[tuple]) -> list[dict]:
    t = pa.table({
        "chrom": pa.array([r[0] for r in rows], pa.string()),
        "start": pa.array([r[1] for r in rows], pa.int64()),
        "ref": pa.array([r[2] for r in rows], pa.string()),
        "alt": pa.array([r[3] for r in rows], pa.string()),
        "SpliceAI": pa.array([[r[4]] for r in rows], pa.list_(pa.string())),
    })
    return run(manifest("spliceai")["ingest_sql"], {"plugin_spliceai_src": t})


def test_spliceai_unpacks_first_prediction():
    (r,) = spliceai_rows([("22", 100, "A", "G", "G|GENE1|0.01|0.00|0.85|0.02|-3|12|1|-40")])
    assert (r["start"], r["allele_string"], r["symbol"], r["symbol_out"]) == (100, "A/G", "GENE1", "GENE1")
    assert (r["ds_ag"], r["ds_al"], r["ds_dg"], r["ds_dl"]) == ("0.01", "0.00", "0.85", "0.02")
    assert (r["dp_ag"], r["dp_al"], r["dp_dg"], r["dp_dl"]) == (-3, 12, 1, -40)


def test_spliceai_same_variant_two_genes_keeps_both_rows():
    rows = spliceai_rows([("22", 100, "A", "G", "G|GENE1|0.1|0|0|0|1|2|3|4"),
                          ("22", 100, "A", "G", "G|GENE2|0.2|0|0|0|1|2|3|4")])
    assert sorted(r["symbol"] for r in rows) == ["GENE1", "GENE2"]


# --- dbNSFP -----------------------------------------------------------------

def dbnsfp_rows(rows: list[dict]) -> list[dict]:
    m = manifest("dbnsfp")
    cols = [c["name"] for c in m["source"][0]["csv"]["schema"]]
    return run(m["ingest_sql"], {"plugin_dbnsfp_src": utf8_table(cols, rows)})


def test_dbnsfp_stop_codon_and_missing_values():
    (r,) = dbnsfp_rows([{"chr": "22", "pos(1-based)": "100", "ref": "C", "alt": "T", "aaref": "R",
                         "aaalt": "X", "SIFT4G_score": ".", "REVEL_score": "0.7;0.8",
                         "MetaLR_pred": "D|T", "CADD_phred": "35"}])
    assert (r["start"], r["allele_string"], r["aa_change"]) == (100, "C/T", "R/*")
    assert r["sift4g_score"] is None      # '.' means no value
    assert r["revel_score"] == "0.7,0.8"  # per-transcript ';' becomes ','
    assert r["metalr_pred"] == "D&T"      # '|' would break CSQ, becomes '&'
    assert r["cadd_phred"] == "35"

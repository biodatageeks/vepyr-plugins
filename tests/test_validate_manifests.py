import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import validate_manifests as vm  # noqa: E402

GFF = '''\
plugin_name = "demo"
coordinate_system = "1-based"
lookup = "interval"
field_order = "alphabetical"
ingest_sql = "SELECT chrom, start, \\"end\\", gene_id, \\"Rat_gene_id\\" AS rat FROM plugin_demo_src"

[[source]]
provider = "gff"
path = "demo.gff3.gz"
url = "https://example.org/demo.gff3.gz"
md5 = "20e5401a198d7d3db66a982c037d3ad4"
index = "tabix"
  [source.gff]
  attributes = ["gene_id", "Rat_gene_id"]

[[match_column]]
column = "gene_id"
template = "{Gene}"

[[value_columns]]
column = "rat"
csq_field = "Demo_Rat"
type = "Utf8"
'''


def _errors(tmp_path: Path, body: str) -> list[str]:
    d = tmp_path / "plugins" / "demo"
    d.mkdir(parents=True, exist_ok=True)
    p = d / "demo.source.toml"
    p.write_text(body, encoding="utf-8")
    errors: list[str] = []
    vm.validate_manifest(p, errors)
    return errors


def test_gff_interval_manifest_is_valid(tmp_path):
    assert _errors(tmp_path, GFF) == []


def test_gff_requires_attributes(tmp_path):
    body = GFF.replace('  [source.gff]\n  attributes = ["gene_id", "Rat_gene_id"]\n', "")
    assert any("source[0].gff" in e and "required" in e for e in _errors(tmp_path, body))
    body = GFF.replace('attributes = ["gene_id", "Rat_gene_id"]', "attributes = []")
    assert any("non-empty" in e for e in _errors(tmp_path, body))


def test_gff_table_rejected_for_other_providers(tmp_path):
    body = GFF.replace('provider = "gff"', 'provider = "bed"').replace('index = "tabix"\n', "")
    assert any("source[0].gff is only valid for provider gff" in e for e in _errors(tmp_path, body))


def test_tabix_allowed_for_gff(tmp_path):
    assert not any("index='tabix'" in e for e in _errors(tmp_path, GFF))


def test_lookup_values(tmp_path):
    body = GFF.replace('lookup = "interval"', 'lookup = "span"')
    assert any("lookup must be one of" in e for e in _errors(tmp_path, body))


def test_interval_rejects_allele_match(tmp_path):
    body = GFF.replace('lookup = "interval"', 'lookup = "interval"\nallele_match = "minimised"')
    assert any("allele_match" in e and "interval" in e for e in _errors(tmp_path, body))


def test_point_default_still_accepts_existing_manifests():
    errors: list[str] = []
    for path in sorted(vm.ROOT.glob(vm.MANIFEST_GLOB)):
        vm.validate_manifest(path, errors)
    assert errors == []

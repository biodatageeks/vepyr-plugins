"""Tests for the parity gate.

The synthetic golden/vepyr pairs here exist to pin the two behaviours that
separate this harness from a naive `diff`:

* a plugin-field mismatch on a key where the **core agrees** with VEP is the
  manifest's fault and must FAIL;
* the *same* mismatch on a key where the **core also diverges** cannot be
  attributed to the manifest, so it is EXCLUDED — reported loudly, but not a
  failure (unless ``--strict``).

Getting that distinction wrong in either direction is fatal: blame the plugin
for a core bug and the gate is a liar; excuse a plugin bug as core drift and the
gate is worse than nothing.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from parity import (
    CORE_FIELDS,
    BuildHealthError,
    ParityConfig,
    Status,
    assess_build_health,
    evaluate,
    resolve_source,
)

# The CSQ layout of the synthetic files: VEP's core fields, then the plugin's,
# exactly as a real `--plugin` run appends them (last field last — the position
# the trailing-newline bug used to corrupt).
CSQ_LAYOUT: tuple[str, ...] = (
    "Allele",
    "Consequence",
    "Feature",
    "Amino_acids",
    "Protein_position",
    "am_class",
    "am_pathogenicity",
)
PLUGIN_FIELDS: tuple[str, ...] = ("am_class", "am_pathogenicity")


def _entry(**overrides: str) -> dict[str, str]:
    """One CSQ entry: a missense hit on ENST1 that AlphaMissense scored."""
    base = {
        "Allele": "G",
        "Consequence": "missense_variant",
        "Feature": "ENST1",
        "Amino_acids": "C/W",
        "Protein_position": "17",
        "am_class": "ambiguous",
        "am_pathogenicity": "0.4833",
    }
    return base | overrides


def _write_vcf(
    path: Path,
    variants: list[tuple[int, str, str, list[dict[str, str]]]],
    *,
    layout: tuple[str, ...] = CSQ_LAYOUT,
    trailing_newline: bool = True,
) -> Path:
    """Write a sites-only (8-column) VCF whose INFO carries a CSQ.

    Sites-only is the shape a region VEP run naturally produces, and the shape
    in which INFO is the final column — so the record terminator sits directly
    after the last CSQ field, which is where plugin fields live.
    """
    lines = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence '
        f'annotations. Format: {"|".join(layout)}">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for pos, ref, alt, entries in variants:
        if entries:
            csq = ",".join("|".join(e.get(f, "") for f in layout) for e in entries)
            info = f"CSQ={csq}"
        else:
            info = "."
        lines.append(f"chr22\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}")
    text = "\n".join(lines) + ("\n" if trailing_newline else "")
    path.write_text(text)
    return path


@pytest.fixture
def golden(tmp_path: Path) -> Path:
    """VEP's output: two variants, each one missense entry that AM scored."""
    return _write_vcf(
        tmp_path / "golden.vcf",
        [
            (22893742, "C", "G", [_entry()]),
            (22900000, "A", "T", [_entry(Feature="ENST2", Amino_acids="R/S",
                                         Protein_position="42",
                                         am_class="likely_benign",
                                         am_pathogenicity="0.1043")]),
        ],
    )


def test_identical_pair_is_clean(golden: Path, tmp_path: Path) -> None:
    """The degenerate case: a file against itself is a pass with nothing excluded."""
    report = evaluate(golden, golden, PLUGIN_FIELDS)

    assert report.status is Status.PASS
    assert report.is_clean
    assert report.excluded_keys == frozenset()
    assert report.surviving_mismatch_count == 0


def test_plugin_mismatch_on_core_agreeing_key_fails(golden: Path, tmp_path: Path) -> None:
    """Core agrees on transcript and amino acid, but the plugin's value differs.

    Nothing upstream can explain this away: it is the manifest's or the plugin
    cache's fault, and it must fail.
    """
    ours = _write_vcf(
        tmp_path / "ours.vcf",
        [
            # Same Feature / Amino_acids / Protein_position — only am_* differ.
            (22893742, "C", "G", [_entry(am_class="likely_pathogenic",
                                         am_pathogenicity="0.9999")]),
            (22900000, "A", "T", [_entry(Feature="ENST2", Amino_acids="R/S",
                                         Protein_position="42",
                                         am_class="likely_benign",
                                         am_pathogenicity="0.1043")]),
        ],
    )
    report = evaluate(golden, ours, PLUGIN_FIELDS)

    assert report.status is Status.FAIL
    assert report.excluded_keys == frozenset(), "core agreed; nothing may be excused"
    assert report.surviving_mismatch_count == 2  # am_class and am_pathogenicity
    assert "chr22\t22893742\tC\tG" in report.surviving_mismatch_keys["am_class"]


def test_same_mismatch_on_core_diverging_key_is_excluded(golden: Path, tmp_path: Path) -> None:
    """The identical plugin-field mismatch, but now the core disagrees too.

    vepyr called a different amino-acid change, so the plugin's discriminator
    (`{ref_aa}{Protein_position}{alt_aa}`) necessarily selected a different row.
    Blaming the manifest for that would be blaming it for a core bug.
    """
    ours = _write_vcf(
        tmp_path / "ours.vcf",
        [
            # Core diverges: Amino_acids C/W -> C/Y. The am_* values differ as a
            # direct consequence.
            (22893742, "C", "G", [_entry(Amino_acids="C/Y",
                                         am_class="likely_pathogenic",
                                         am_pathogenicity="0.9999")]),
            (22900000, "A", "T", [_entry(Feature="ENST2", Amino_acids="R/S",
                                         Protein_position="42",
                                         am_class="likely_benign",
                                         am_pathogenicity="0.1043")]),
        ],
    )
    report = evaluate(golden, ours, PLUGIN_FIELDS)

    assert report.status is Status.PASS, "core drift must not fail the plugin gate"
    assert report.excluded_keys == frozenset({"chr22\t22893742\tC\tG"})
    assert report.surviving_mismatch_count == 0
    # Loudly and separately: never folded into zero.
    assert "core drift" in report.render().lower()
    assert "1" in report.render()


def test_core_drift_fails_under_strict(golden: Path, tmp_path: Path) -> None:
    """`--strict` is the mode for PRs against the engine: drift fails too."""
    ours = _write_vcf(
        tmp_path / "ours.vcf",
        [
            (22893742, "C", "G", [_entry(Amino_acids="C/Y",
                                         am_class="likely_pathogenic",
                                         am_pathogenicity="0.9999")]),
            (22900000, "A", "T", [_entry(Feature="ENST2", Amino_acids="R/S",
                                         Protein_position="42",
                                         am_class="likely_benign",
                                         am_pathogenicity="0.1043")]),
        ],
    )
    report = evaluate(golden, ours, PLUGIN_FIELDS)

    assert report.status is Status.PASS
    assert report.status_under(strict=True) is Status.FAIL


def test_over_emission_fails(golden: Path, tmp_path: Path) -> None:
    """vepyr fills a field VEP left empty — one of the two real bugs found by hand.

    VEP emits AlphaMissense only for missense variants. A value on a synonymous
    entry means the plugin cache leaked, and no amount of core drift excuses it:
    the core agrees on every field here.
    """
    golden_syn = _write_vcf(
        tmp_path / "golden_syn.vcf",
        [
            (22893742, "C", "G", [_entry(Consequence="synonymous_variant",
                                         Amino_acids="", Protein_position="",
                                         am_class="", am_pathogenicity="")]),
        ],
    )
    ours = _write_vcf(
        tmp_path / "ours_syn.vcf",
        [
            # Same core fields; vepyr invented am_* where VEP emitted none.
            (22893742, "C", "G", [_entry(Consequence="synonymous_variant",
                                         Amino_acids="", Protein_position="",
                                         am_class="likely_benign",
                                         am_pathogenicity="0.1")]),
        ],
    )
    report = evaluate(golden_syn, ours, PLUGIN_FIELDS)

    assert report.status is Status.FAIL
    assert report.surviving_over_emission_count == 2
    assert "over-emission" in report.render().lower()


def test_csq_missing_in_test_fails(golden: Path, tmp_path: Path) -> None:
    """A total annotation dropout: vepyr emitted no CSQ at all for a variant."""
    ours = _write_vcf(
        tmp_path / "ours.vcf",
        [
            (22893742, "C", "G", []),  # no CSQ whatsoever
            (22900000, "A", "T", [_entry(Feature="ENST2", Amino_acids="R/S",
                                         Protein_position="42",
                                         am_class="likely_benign",
                                         am_pathogenicity="0.1043")]),
        ],
    )
    report = evaluate(golden, ours, PLUGIN_FIELDS)

    assert report.status is Status.FAIL
    assert report.gate.csq_missing_in_test == 1


def test_golden_without_plugin_field_fails(golden: Path, tmp_path: Path) -> None:
    """A golden that lacks the plugin's fields was produced without `--plugin`.

    That is a broken golden, not a pass. It must never read as "0 mismatches".
    """
    core_only = tuple(f for f in CSQ_LAYOUT if not f.startswith("am_"))
    no_plugin = _write_vcf(
        tmp_path / "no_plugin.vcf",
        [
            (22893742, "C", "G", [_entry()]),
            (22900000, "A", "T", [_entry(Feature="ENST2", Amino_acids="R/S",
                                         Protein_position="42")]),
        ],
        layout=core_only,
    )
    report = evaluate(no_plugin, golden, PLUGIN_FIELDS)

    assert report.status is Status.FAIL
    assert report.gate.fields_missing_from_truth


def test_no_trailing_newline_does_not_manufacture_a_mismatch(
    golden: Path, tmp_path: Path
) -> None:
    """The last CSQ field is where plugin fields live — guard the regression."""
    same = _write_vcf(
        tmp_path / "same_no_nl.vcf",
        [
            (22893742, "C", "G", [_entry()]),
            (22900000, "A", "T", [_entry(Feature="ENST2", Amino_acids="R/S",
                                         Protein_position="42",
                                         am_class="likely_benign",
                                         am_pathogenicity="0.1043")]),
        ],
        trailing_newline=False,
    )
    report = evaluate(golden, same, PLUGIN_FIELDS)

    assert report.status is Status.PASS, "a phantom mismatch on the final field"


# --------------------------------------------------------------------------
# Build health: caught before any comparison, because the diff would otherwise
# be a wall of empties instead of one clear sentence.
# --------------------------------------------------------------------------


def test_warm_zero_with_rows_is_a_build_failure() -> None:
    """Not one source row joined the variation cache: the plugin cache is dead."""
    with pytest.raises(BuildHealthError, match="warm"):
        assess_build_health([("chr22", 54624, 0, 54624)])


def test_warm_zero_with_no_rows_is_fine() -> None:
    """No rows in the region is not a build failure — there was nothing to join."""
    assess_build_health([("chr22", 0, 0, 0)])  # must not raise


def test_healthy_build_passes() -> None:
    """The real AlphaMissense mini-cache build shape."""
    assess_build_health([("chr22", 54624, 147, 54477)])  # must not raise


# --------------------------------------------------------------------------
# Licence-gated plugins: skipped, but never silently.
# --------------------------------------------------------------------------


def _write_toml(path: Path, *, redistributable: bool) -> Path:
    path.write_text(
        f"""
plugin          = "alphamissense"
csq_fields      = ["am_class", "am_pathogenicity"]
region          = "chr22:22000000-23500000"
redistributable = {str(redistributable).lower()}

[vep]
release     = "115"
plugin_args = "AlphaMissense,file={{source}}"
"""
    )
    return path


def test_non_redistributable_without_fixture_skips_loudly(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """CC-BY-NC data cannot be committed, so CI cannot build it — say so out loud.

    A licence-gated plugin that quietly goes green is the worst of both worlds.
    """
    cfg = ParityConfig.load(_write_toml(tmp_path / "parity.toml", redistributable=False))
    source = resolve_source(cfg, override=None)

    assert source is None  # -> SKIP
    out = capsys.readouterr().out
    assert "SKIP" in out
    assert "redistributable" in out.lower()
    assert cfg.plugin in out


def test_redistributable_without_fixture_is_an_error(tmp_path: Path) -> None:
    """If the data is redistributable, a missing fixture is a broken port, not a skip."""
    cfg = ParityConfig.load(_write_toml(tmp_path / "parity.toml", redistributable=True))

    with pytest.raises(FileNotFoundError, match="fixture"):
        resolve_source(cfg, override=None)


def test_config_round_trip(tmp_path: Path) -> None:
    """The schema the plan specifies, parsed."""
    cfg = ParityConfig.load(_write_toml(tmp_path / "parity.toml", redistributable=False))

    assert cfg.plugin == "alphamissense"
    assert cfg.csq_fields == ("am_class", "am_pathogenicity")
    assert cfg.chrom == "chr22"
    assert cfg.span == (22_000_000, 23_500_000)
    assert cfg.vep_release == "115"
    assert not cfg.redistributable
    assert CORE_FIELDS == ("Feature", "Consequence", "Amino_acids", "Protein_position")

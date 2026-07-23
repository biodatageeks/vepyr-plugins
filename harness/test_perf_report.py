"""Tests for the cost report: perf JSON -> side-by-side table -> markdown / PDF."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from perf_report import CostRow, build_cost_table, render_markdown, render_pdf

_HAVE_MPL = importlib.util.find_spec("matplotlib") is not None


def _write_perf_json(
    dir_: Path,
    plugin: str,
    param: str,
    *,
    base_s: float,
    plugin_s: float,
    cost_rss_bytes: int,
    disk_bytes: int,
) -> Path:
    """Write a perf_<plugin>_<param>.json in the exact shape the driver emits."""
    doc = {
        "plugin": plugin,
        "param": param,
        "build": {"elapsed_s": 0.1, "peak_rss_bytes": 1},
        "baseline": {"elapsed_s": base_s, "peak_rss_bytes": 500_000},
        "withplugin": {"elapsed_s": plugin_s, "peak_rss_bytes": 500_000 + cost_rss_bytes},
        "cache_bytes": disk_bytes,
        "baseline_out_bytes": 10,
        "withplugin_out_bytes": 12,
        "cost": {
            "elapsed_s": round(plugin_s - base_s, 6),
            "peak_rss_bytes": cost_rss_bytes,
            "disk_bytes": disk_bytes,
        },
    }
    path = dir_ / f"perf_{plugin}_{param}.json"
    path.write_text(json.dumps(doc, indent=2))
    return path


@pytest.fixture
def perf_jsons(tmp_path: Path) -> list[Path]:
    """Two plugins with different cost so sorting and % are both exercised."""
    return [
        _write_perf_json(
            tmp_path, "alphamissense", "hgvs",
            base_s=3.0, plugin_s=4.5,  # cost 1.5s -> 50%
            cost_rss_bytes=2 * 1024 * 1024, disk_bytes=4 * 1024 * 1024,
        ),
        _write_perf_json(
            tmp_path, "spliceai", "hgvs",
            base_s=2.0, plugin_s=6.0,  # cost 4.0s -> 200%
            cost_rss_bytes=8 * 1024 * 1024, disk_bytes=16 * 1024 * 1024,
        ),
    ]


def test_row_fields_and_derived_values(perf_jsons: list[Path]) -> None:
    rows = build_cost_table(perf_jsons)
    assert all(isinstance(r, CostRow) for r in rows)

    by_plugin = {r.plugin: r for r in rows}
    am = by_plugin["alphamissense"]
    assert am.param == "hgvs"
    assert am.base_s == 3.0
    assert am.plugin_s == 4.5
    assert am.cost_s == pytest.approx(1.5)
    assert am.cost_rss_mb == pytest.approx(2.0)   # 2 MiB
    assert am.cost_disk_mb == pytest.approx(4.0)  # 4 MiB
    assert am.cost_pct == pytest.approx(50.0)     # 1.5 / 3.0 * 100


def test_cost_pct_is_cost_over_base_times_100(perf_jsons: list[Path]) -> None:
    rows = build_cost_table(perf_jsons)
    for r in rows:
        assert r.cost_pct == pytest.approx(r.cost_s / r.base_s * 100.0)


def test_rows_sorted_by_cost_s_desc(perf_jsons: list[Path]) -> None:
    rows = build_cost_table(perf_jsons)
    costs = [r.cost_s for r in rows]
    assert costs == sorted(costs, reverse=True)
    assert rows[0].plugin == "spliceai"  # 4.0s cost outranks alphamissense's 1.5s


def test_zero_baseline_does_not_divide_by_zero(tmp_path: Path) -> None:
    """A degenerate zero-second baseline must not crash the table."""
    p = _write_perf_json(
        tmp_path, "edge", "everything",
        base_s=0.0, plugin_s=0.0, cost_rss_bytes=0, disk_bytes=0,
    )
    rows = build_cost_table([p])
    assert rows[0].cost_pct == 0.0


def test_render_markdown_is_a_table(perf_jsons: list[Path]) -> None:
    rows = build_cost_table(perf_jsons)
    md = render_markdown(rows)

    assert "| plugin" in md or "plugin" in md.splitlines()[0]
    assert "spliceai" in md and "alphamissense" in md
    # A GitHub-flavoured table has a header separator row of dashes.
    assert "---" in md
    # spliceai (higher cost) is listed before alphamissense.
    assert md.index("spliceai") < md.index("alphamissense")


def test_build_cost_table_accepts_paths_or_dicts(tmp_path: Path) -> None:
    """A caller may pass loaded dicts as well as file paths."""
    p = _write_perf_json(
        tmp_path, "am", "hgvs",
        base_s=1.0, plugin_s=2.0, cost_rss_bytes=1024 * 1024, disk_bytes=1024 * 1024,
    )
    doc = json.loads(p.read_text())
    rows = build_cost_table([doc])
    assert rows[0].plugin == "am"
    assert rows[0].cost_s == pytest.approx(1.0)


@pytest.mark.skipif(not _HAVE_MPL, reason="matplotlib not installed")
def test_render_pdf_writes_a_file(perf_jsons: list[Path], tmp_path: Path) -> None:
    rows = build_cost_table(perf_jsons)
    out = tmp_path / "perf.pdf"
    returned = render_pdf(rows, out)

    assert returned == out
    assert out.is_file()
    assert out.stat().st_size > 0
    assert out.read_bytes()[:4] == b"%PDF"

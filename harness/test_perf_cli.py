"""Tests for the sbatch-facing CLI entrypoints.

The gate sbatch shells out to ``python -m perf_plugins`` (once per plugin, looping
its ``--params``) and ``python -m perf_report`` (globs the driver's JSON, writes a
markdown + PDF). These tests pin the arg-parsing and dispatch without touching the
heavy ``vepyr``/matplotlib paths: :func:`run_perf` and :func:`render_pdf` are
monkeypatched, so the CLIs are exercised on a bare machine.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import perf_plugins
import perf_report


# ---------------------------------------------------------------------------
# perf_plugins CLI: loop --params, call run_perf per (plugin, param).
# ---------------------------------------------------------------------------


def test_perf_plugins_cli_dispatches_run_perf_per_param(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    calls: list[dict[str, object]] = []

    def fake_run_perf(
        plugin, param, region_vcf, mini_cache, plugins_repo, source_path, out_dir,
        *, version="HEAD", overwrite=True,
    ):
        calls.append(
            {
                "plugin": plugin,
                "param": param,
                "region_vcf": Path(region_vcf),
                "mini_cache": Path(mini_cache),
                "plugins_repo": Path(plugins_repo),
                "source_path": Path(source_path),
                "out_dir": Path(out_dir),
                "version": version,
                "overwrite": overwrite,
            }
        )
        # Emulate the real driver: create out_dir, then drop its JSON.
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        (Path(out_dir) / f"perf_{plugin}_{param}.json").write_text("{}")
        return object()

    monkeypatch.setattr(perf_plugins, "run_perf", fake_run_perf)

    out = tmp_path / "out"
    rc = perf_plugins.main(
        [
            "--plugin", "alphamissense",
            "--params", "everything,hgvs",
            "--region-vcf", str(tmp_path / "r.vcf"),
            "--mini-cache", str(tmp_path / "mc"),
            "--plugins-repo", str(tmp_path / "repo"),
            "--source-path", str(tmp_path / "am.tsv.gz"),
            "--out", str(out),
        ]
    )

    assert rc == 0
    # One call per param, in order, all for the one plugin.
    assert [c["param"] for c in calls] == ["everything", "hgvs"]
    assert all(c["plugin"] == "alphamissense" for c in calls)
    # Paths are forwarded verbatim.
    assert calls[0]["region_vcf"] == tmp_path / "r.vcf"
    assert calls[0]["mini_cache"] == tmp_path / "mc"
    assert calls[0]["plugins_repo"] == tmp_path / "repo"
    assert calls[0]["source_path"] == tmp_path / "am.tsv.gz"
    assert calls[0]["out_dir"] == out
    # Defaults: HEAD manifest, overwrite on.
    assert calls[0]["version"] == "HEAD"
    assert calls[0]["overwrite"] is True


def test_perf_plugins_cli_splits_and_trims_params(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Whitespace + trailing commas in --params are tolerated."""
    seen: list[str] = []
    monkeypatch.setattr(
        perf_plugins, "run_perf",
        lambda plugin, param, *a, **k: seen.append(param) or object(),
    )
    rc = perf_plugins.main(
        [
            "--plugin", "am",
            "--params", " everything , hgvs ,",
            "--region-vcf", str(tmp_path / "r.vcf"),
            "--mini-cache", str(tmp_path),
            "--plugins-repo", str(tmp_path),
            "--source-path", str(tmp_path / "am.tsv"),
            "--out", str(tmp_path / "o"),
        ]
    )
    assert rc == 0
    assert seen == ["everything", "hgvs"]


def test_perf_plugins_cli_requires_plugin(tmp_path: Path) -> None:
    """A missing required arg is an argparse usage error (exit 2)."""
    with pytest.raises(SystemExit) as ei:
        perf_plugins.main(["--params", "hgvs", "--out", str(tmp_path)])
    assert ei.value.code == 2


# ---------------------------------------------------------------------------
# perf_report CLI: glob perf_*.json, write markdown + PDF.
# ---------------------------------------------------------------------------


def _write_perf_json(dir_: Path, plugin: str, param: str, base_s: float, plugin_s: float) -> None:
    doc = {
        "plugin": plugin,
        "param": param,
        "build": {"elapsed_s": 0.1, "peak_rss_bytes": 1},
        "baseline": {"elapsed_s": base_s, "peak_rss_bytes": 1},
        "withplugin": {"elapsed_s": plugin_s, "peak_rss_bytes": 2},
        "cache_bytes": 1024,
        "baseline_out_bytes": 1,
        "withplugin_out_bytes": 1,
        "cost": {
            "elapsed_s": round(plugin_s - base_s, 6),
            "peak_rss_bytes": 1,
            "disk_bytes": 1024,
        },
    }
    (dir_ / f"perf_{plugin}_{param}.json").write_text(json.dumps(doc))


def test_perf_report_cli_globs_and_writes_md_and_pdf(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    indir = tmp_path / "in"
    indir.mkdir()
    _write_perf_json(indir, "alphamissense", "hgvs", 2.0, 3.0)
    _write_perf_json(indir, "spliceai", "hgvs", 2.0, 6.0)
    # A stray non-perf JSON must be ignored by the glob.
    (indir / "verdict.json").write_text("{}")

    rendered: dict[str, object] = {}
    monkeypatch.setattr(
        perf_report,
        "render_pdf",
        lambda rows, out: rendered.setdefault("pdf", (list(rows), Path(out))) or Path(out),
    )

    md = tmp_path / "perf.md"
    pdf = tmp_path / "perf.pdf"
    rc = perf_report.main(["--in", str(indir), "--md", str(md), "--pdf", str(pdf)])

    assert rc == 0
    assert md.is_file()
    text = md.read_text()
    assert "plugin" in text
    assert "spliceai" in text and "alphamissense" in text
    # PDF render was invoked with the same rows and the requested path.
    assert "pdf" in rendered
    rows, pdf_out = rendered["pdf"]  # type: ignore[misc]
    assert pdf_out == pdf
    assert len(rows) == 2  # only the two perf_*.json, not verdict.json


def test_perf_report_cli_errors_on_empty_dir(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        perf_report.main(
            ["--in", str(tmp_path), "--md", str(tmp_path / "m.md"), "--pdf", str(tmp_path / "p.pdf")]
        )

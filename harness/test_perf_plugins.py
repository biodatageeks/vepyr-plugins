"""Tests for the per-plugin perf driver.

Two layers, mirroring how ``test_parity_harness.py`` separates pure logic from
the heavy real-``vepyr`` path:

* **pure unit tests** — the cost arithmetic, the ``PARAM_FLAGS`` → ``annotate()``
  kwarg mapping, and the ``perf_<plugin>_<param>.json`` round-trip. No ``vepyr``,
  no cache; always run.
* **one env-gated integration test** — the real ``run_perf`` against the local
  region mini-cache + AlphaMissense slice. Skipped (never failed) when ``vepyr``
  or the fixtures are absent, exactly like the harness gates its real-VEP tests.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from perf import PerfSample
from perf_plugins import PARAM_FLAGS, CostDelta, PerfResult, run_perf

# --------------------------------------------------------------------------
# Local fixtures for the opt-in integration test (this Mac).
# --------------------------------------------------------------------------
_MINI_CACHE = Path(os.environ.get("VEPYR_MINI_CACHE", "/tmp/mini_cache_region"))
_REGION_VCF = Path(__file__).resolve().parent / "regions" / "chr22-22.0-23.5Mb.vcf"
_REPO_ROOT = Path(__file__).resolve().parent.parent
_AM_SLICE = Path(
    os.environ.get(
        "VEPYR_AM_SLICE",
        "/Users/wojtek/Documents/vepyr/_parity_data/alphamissense/"
        "AlphaMissense_hg38.chr22_region.tsv.gz",
    )
)


def _vepyr_available() -> bool:
    try:
        import vepyr  # noqa: F401
    except Exception:
        return False
    return True


_INTEGRATION_READY = (
    _vepyr_available()
    and _MINI_CACHE.is_dir()
    and _REGION_VCF.is_file()
    and _AM_SLICE.is_file()
)


# --------------------------------------------------------------------------
# Pure: the param → annotate() kwargs mapping.
# --------------------------------------------------------------------------


def test_param_flags_are_the_two_profiles() -> None:
    """Exactly the {everything, hgvs} profiles the benchmark reports."""
    assert set(PARAM_FLAGS) == {"everything", "hgvs"}


def test_param_flags_map_to_real_annotate_kwargs() -> None:
    """Each profile is expressible as ``vepyr.annotate`` boolean kwargs.

    Derived from ``annotate()``'s signature (named ``everything`` / ``hgvs``
    flags), not a VEP-flag string — the toolkit exposes them directly.
    """
    assert PARAM_FLAGS["everything"] == {"everything": True}
    # VEP --hgvs = hgvsc + hgvsp; the toolkit's hgvs=True implies both.
    assert PARAM_FLAGS["hgvs"] == {"hgvs": True}


# --------------------------------------------------------------------------
# Pure: cost = with-plugin − baseline, per metric.
# --------------------------------------------------------------------------


def _result(**over: object) -> PerfResult:
    base = dict(
        plugin="demo",
        param="hgvs",
        build=PerfSample(elapsed_s=3.0, peak_rss_bytes=100),
        baseline=PerfSample(elapsed_s=2.0, peak_rss_bytes=500_000),
        withplugin=PerfSample(elapsed_s=5.5, peak_rss_bytes=900_000),
        cache_bytes=4_000_000,
        baseline_out_bytes=1_000,
        withplugin_out_bytes=1_500,
    )
    base.update(over)
    return PerfResult.from_samples(**base)


def test_cost_is_withplugin_minus_baseline_per_metric() -> None:
    r = _result()
    assert r.cost.elapsed_s == pytest.approx(3.5)  # 5.5 − 2.0
    assert r.cost.peak_rss_bytes == 400_000  # 900_000 − 500_000
    assert r.cost.disk_bytes == 4_000_000  # the plugin cache is the marginal disk


def test_negative_time_cost_is_kept_not_clamped() -> None:
    """Noise can make +plugin faster than baseline; report the honest delta."""
    r = _result(withplugin=PerfSample(elapsed_s=1.9, peak_rss_bytes=900_000))
    assert r.cost.elapsed_s == pytest.approx(-0.1)


def test_unknown_param_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="param"):
        run_perf(
            "alphamissense",
            "nonsense",
            _REGION_VCF,
            _MINI_CACHE,
            _REPO_ROOT,
            _AM_SLICE,
            tmp_path,
        )


# --------------------------------------------------------------------------
# Pure: the JSON the report consumes.
# --------------------------------------------------------------------------


def test_write_json_round_trips(tmp_path: Path) -> None:
    r = _result(plugin="alphamissense", param="hgvs")
    path = r.write_json(tmp_path)

    assert path == tmp_path / "perf_alphamissense_hgvs.json"
    loaded = json.loads(path.read_text())
    assert loaded["plugin"] == "alphamissense"
    assert loaded["param"] == "hgvs"
    assert loaded["baseline"]["elapsed_s"] == 2.0
    assert loaded["withplugin"]["peak_rss_bytes"] == 900_000
    assert loaded["cost"]["elapsed_s"] == pytest.approx(3.5)
    assert loaded["cost"]["disk_bytes"] == 4_000_000
    assert loaded["cache_bytes"] == 4_000_000


# --------------------------------------------------------------------------
# Integration: the real driver, opt-in and env-gated.
# --------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.skipif(
    not _INTEGRATION_READY,
    reason="needs vepyr + region mini-cache + AlphaMissense slice",
)
def test_run_perf_real_alphamissense(tmp_path: Path) -> None:
    """Full baseline-vs-plugin run against the local mini-cache."""
    result = run_perf(
        "alphamissense",
        "hgvs",
        _REGION_VCF,
        _MINI_CACHE,
        _REPO_ROOT,
        _AM_SLICE,
        tmp_path,
    )

    # The cost identity the whole benchmark rests on.
    assert result.cost.elapsed_s == pytest.approx(
        result.withplugin.elapsed_s - result.baseline.elapsed_s, abs=1e-6
    )
    # A plugin cache was built and has real bytes on disk.
    assert result.cache_bytes > 0
    assert result.cost.disk_bytes == result.cache_bytes
    # Both runs produced output; the JSON landed where the report reads it.
    assert result.baseline.elapsed_s > 0
    assert result.withplugin.elapsed_s > 0
    assert (tmp_path / "perf_alphamissense_hgvs.json").is_file()

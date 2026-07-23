"""Tests for the perf-capture primitives.

These are pure: no ``vepyr``, no cache, no subprocess. They pin the two things
every downstream perf number is built from — an elapsed/peak-RSS sample from a
context manager, and a byte-accounting helper for files and directories.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from perf import PerfSample, dir_bytes, file_bytes, measure, measure_call


def _alloc_mb(mb: int) -> int:
    """Allocate and touch ~``mb`` MiB so it counts against RSS; return its length.

    Module-level (so ``spawn`` can pickle it) and page-touching (a zero-filled
    ``bytearray`` is resident) — a deterministic, dependency-free RSS workload.
    """
    buf = bytearray(mb * 1024 * 1024)
    buf[::4096] = b"\x01" * len(buf[::4096])  # fault every page in
    return len(buf)


def test_measure_captures_elapsed_and_rss() -> None:
    """``with measure() as s`` exposes elapsed + peak RSS after the block exits."""
    with measure() as s:
        _ = [0] * 2_000_000  # allocate ~16 MB
        time.sleep(0.05)

    assert isinstance(s, PerfSample)
    assert s.elapsed_s >= 0.05
    # ru_maxrss is the process's monotonic peak, so an absolute reading is always
    # positive; the point of the sample is that the value is populated on exit.
    assert s.peak_rss_bytes > 0


def test_sample_as_dict_round_trips() -> None:
    """The JSON the driver dumps is exactly the sample's two fields."""
    with measure() as s:
        pass
    d = s.as_dict()
    assert set(d) == {"elapsed_s", "peak_rss_bytes"}
    assert d["elapsed_s"] == s.elapsed_s
    assert d["peak_rss_bytes"] == s.peak_rss_bytes


def test_file_and_dir_bytes(tmp_path: Path) -> None:
    """``file_bytes`` is one file; ``dir_bytes`` is a recursive sum."""
    p = tmp_path / "a.bin"
    p.write_bytes(b"x" * 1234)
    assert file_bytes(p) == 1234

    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "b").write_bytes(b"y" * 10)
    assert dir_bytes(tmp_path) == 1244


def test_dir_bytes_ignores_directories_themselves(tmp_path: Path) -> None:
    """Only regular files count; empty subdirs contribute nothing."""
    (tmp_path / "empty").mkdir()
    assert dir_bytes(tmp_path) == 0


# --------------------------------------------------------------------------
# Subprocess isolation: each measured call gets its own process, so ru_maxrss
# is that call's own peak and never the high-water mark of a prior run.
# --------------------------------------------------------------------------


def test_measure_call_returns_result_and_sample() -> None:
    """The child's return value comes back, with a populated sample."""
    result, s = measure_call(_alloc_mb, 8)
    assert result == 8 * 1024 * 1024
    assert isinstance(s, PerfSample)
    assert s.elapsed_s >= 0.0
    assert s.peak_rss_bytes > 0


def test_measure_call_peak_is_isolated_per_call() -> None:
    """A big run must not inflate a later small run's peak.

    ``ru_maxrss`` is a per-process monotonic high-water mark: two sequential
    workloads in one process would make the second report at least the first's
    peak. Running each in a fresh process is exactly what breaks that coupling,
    so the 5 MiB run lands far below the 200 MiB run despite running second.
    """
    _, big = measure_call(_alloc_mb, 200)
    _, small = measure_call(_alloc_mb, 5)
    assert small.peak_rss_bytes < big.peak_rss_bytes
    # And the gap is real, not a rounding wobble: at least ~100 MiB apart.
    assert big.peak_rss_bytes - small.peak_rss_bytes > 100 * 1024 * 1024


def test_measure_call_surfaces_child_exceptions() -> None:
    """A failure in the child is re-raised in the parent, not swallowed."""
    with pytest.raises(RuntimeError, match="isolated call"):
        measure_call(_boom)


def _boom() -> None:
    raise ValueError("kaboom in the child")

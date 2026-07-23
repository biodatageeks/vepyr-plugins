"""Tests for the perf-capture primitives.

These are pure: no ``vepyr``, no cache, no subprocess. They pin the two things
every downstream perf number is built from — an elapsed/peak-RSS sample from a
context manager, and a byte-accounting helper for files and directories.
"""

from __future__ import annotations

import time
from pathlib import Path

from perf import PerfSample, dir_bytes, file_bytes, measure


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

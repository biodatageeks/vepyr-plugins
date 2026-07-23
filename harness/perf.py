"""Perf-capture primitives for the per-plugin benchmark.

Three things, kept deliberately small and dependency-free so they are trivially
testable and safe to run inside a spawned child process:

* :class:`PerfSample` — one measurement: wall-clock elapsed + peak RSS.
* :func:`measure` — a context manager that fills a sample for the block it wraps.
* :func:`file_bytes` / :func:`dir_bytes` — on-disk cost accounting.

Peak-RSS semantics
------------------
``resource.getrusage(RUSAGE_SELF).ru_maxrss`` is the process's *monotonic* peak,
so :func:`measure` reports the process peak at block exit — an absolute reading,
not a delta. That is honest only when one workload owns the process: a second
workload in the same process inherits the first's high-water mark. The perf
driver therefore runs each ``annotate()`` in its own fresh process via
:func:`measure_call`, where the peak is that run's alone.
"""

from __future__ import annotations

import os
import resource
import sys
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterator

# ru_maxrss is KiB on Linux, bytes on macOS/BSD — normalise to bytes.
_MAXRSS_UNIT: int = 1 if sys.platform == "darwin" else 1024


def _peak_rss_bytes() -> int:
    """This process's peak resident set size, in bytes (platform-normalised)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * _MAXRSS_UNIT


@dataclass
class PerfSample:
    """A single (elapsed, peak-RSS) measurement.

    Mutable by design: :func:`measure` yields a placeholder instance and fills it
    in the block's ``finally``, so ``with measure() as s: ...; s.elapsed_s`` reads
    the finished value without any proxy indirection.
    """

    elapsed_s: float = 0.0
    peak_rss_bytes: int = 0

    def as_dict(self) -> dict[str, float | int]:
        """The two fields, ready to serialise into the driver's JSON."""
        return asdict(self)


@contextmanager
def measure() -> Iterator[PerfSample]:
    """Time + peak RSS of the enclosed block.

    Peak RSS is the process-wide peak (``RUSAGE_SELF``), so the reading is only
    attributable when a single workload runs per process — see the module note
    and :func:`measure_call`.
    """
    sample = PerfSample()
    t0 = time.perf_counter()
    try:
        yield sample
    finally:
        sample.elapsed_s = round(time.perf_counter() - t0, 6)
        sample.peak_rss_bytes = _peak_rss_bytes()


def file_bytes(p: Path | str) -> int:
    """Size of a single file, in bytes."""
    return os.path.getsize(p)


def dir_bytes(p: Path | str) -> int:
    """Recursive on-disk size of every regular file under ``p``, in bytes."""
    return sum(f.stat().st_size for f in Path(p).rglob("*") if f.is_file())

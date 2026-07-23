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

import multiprocessing as mp
import os
import resource
import sys
import time
import traceback
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


# ---------------------------------------------------------------------------
# Subprocess isolation: honest per-run peak RSS.
# ---------------------------------------------------------------------------

# "spawn" (not "fork") gives a pristine interpreter with its own address space,
# so ru_maxrss starts from the interpreter baseline rather than inheriting the
# parent's high-water mark — the whole point of isolating each annotate().
_SPAWN = mp.get_context("spawn")


def _child_target(conn, fn, args, kwargs) -> None:  # pragma: no cover - child process
    """Run ``fn(*args, **kwargs)`` in the child; ship (result, sample) back.

    The child is fresh, so its ``ru_maxrss`` is this call's own peak. On failure
    the formatted traceback is sent instead, so the error surfaces in the parent
    rather than vanishing with the process.
    """
    try:
        with measure() as sample:
            result = fn(*args, **kwargs)
        conn.send(("ok", result, sample.elapsed_s, sample.peak_rss_bytes))
    except BaseException:  # noqa: BLE001 - relayed and re-raised in the parent
        conn.send(("err", traceback.format_exc(), 0.0, 0))
    finally:
        conn.close()


def measure_call[T](fn, /, *args, **kwargs) -> tuple[T, PerfSample]:
    """Run ``fn(*args, **kwargs)`` in a fresh spawned process; measure it there.

    Returns ``(result, PerfSample)`` where ``peak_rss_bytes`` is the child's own
    peak — uncontaminated by anything the parent, or a prior call, allocated.
    ``fn`` and its arguments must be picklable (``spawn`` re-imports the child):
    the module-level ``vepyr.annotate`` and its str/bool arguments are.

    Raises:
        RuntimeError: the child raised (its traceback is embedded) or died
            without sending a result (e.g. an OOM kill).
    """
    parent_conn, child_conn = _SPAWN.Pipe(duplex=False)
    proc = _SPAWN.Process(target=_child_target, args=(child_conn, fn, args, kwargs))
    proc.start()
    child_conn.close()  # parent holds only the read end, so recv unblocks on exit
    name = getattr(fn, "__name__", repr(fn))
    try:
        payload = parent_conn.recv()
    except EOFError as exc:
        proc.join()
        raise RuntimeError(
            f"isolated call to {name!r} produced no result "
            f"(child exit code {proc.exitcode})"
        ) from exc
    finally:
        parent_conn.close()
    proc.join()

    status, result, elapsed_s, peak_rss_bytes = payload
    if status == "err":
        raise RuntimeError(f"isolated call to {name!r} failed:\n{result}")
    return result, PerfSample(elapsed_s=elapsed_s, peak_rss_bytes=peak_rss_bytes)

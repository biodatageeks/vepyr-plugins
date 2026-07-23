#!/usr/bin/env python3
"""The cost report: per-plugin perf JSON -> one side-by-side table.

Reads the ``perf_<plugin>_<param>.json`` files the driver (:mod:`perf_plugins`)
emits and folds them into cost rows — baseline vs +plugin elapsed, plus the
marginal RSS/disk/percent overhead the plugin adds — rendered as GitHub-flavoured
markdown (for the PR comment / job summary) and a one-page matplotlib PDF (the
gate artifact). Rows are sorted by time cost, worst offender first.

``cost_pct = cost_s / base_s * 100`` — the plugin's slowdown as a fraction of the
engine's own baseline, which is the number a reviewer actually reasons about.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

_MIB: float = 1024.0 * 1024.0


@dataclass(frozen=True, slots=True)
class CostRow:
    """One benchmarked ``(plugin, param)`` reduced to the numbers a table shows."""

    plugin: str
    param: str
    base_s: float
    plugin_s: float
    cost_s: float
    cost_rss_mb: float
    cost_disk_mb: float
    cost_pct: float

    @classmethod
    def from_perf(cls, doc: dict[str, object]) -> "CostRow":
        """Build a row from one perf JSON document."""
        baseline = doc["baseline"]  # type: ignore[index]
        withplugin = doc["withplugin"]  # type: ignore[index]
        cost = doc["cost"]  # type: ignore[index]
        base_s = float(baseline["elapsed_s"])  # type: ignore[index]
        cost_s = float(cost["elapsed_s"])  # type: ignore[index]
        return cls(
            plugin=str(doc["plugin"]),  # type: ignore[index]
            param=str(doc["param"]),  # type: ignore[index]
            base_s=base_s,
            plugin_s=float(withplugin["elapsed_s"]),  # type: ignore[index]
            cost_s=cost_s,
            cost_rss_mb=float(cost["peak_rss_bytes"]) / _MIB,  # type: ignore[index]
            cost_disk_mb=float(cost["disk_bytes"]) / _MIB,  # type: ignore[index]
            # Guard the degenerate zero-second baseline rather than raise: a
            # sub-millisecond region can round to 0.0 and must not crash a report.
            cost_pct=(cost_s / base_s * 100.0) if base_s else 0.0,
        )

    def as_dict(self) -> dict[str, str | float]:
        return {
            "plugin": self.plugin,
            "param": self.param,
            "base_s": self.base_s,
            "plugin_s": self.plugin_s,
            "cost_s": self.cost_s,
            "cost_rss_mb": self.cost_rss_mb,
            "cost_disk_mb": self.cost_disk_mb,
            "cost_pct": self.cost_pct,
        }


def _load(item: dict[str, object] | Path | str) -> dict[str, object]:
    """Accept an already-loaded perf dict, or a path to a perf JSON file."""
    if isinstance(item, dict):
        return item
    return json.loads(Path(item).read_text())


def build_cost_table(
    perf_jsons: Iterable[dict[str, object] | Path | str],
) -> list[CostRow]:
    """Fold perf JSON (paths or dicts) into cost rows, worst time cost first."""
    rows = [CostRow.from_perf(_load(item)) for item in perf_jsons]
    rows.sort(key=lambda r: r.cost_s, reverse=True)
    return rows


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

_HEADERS: tuple[str, ...] = (
    "plugin",
    "param",
    "base_s",
    "plugin_s",
    "cost_s",
    "cost_rss_mb",
    "cost_disk_mb",
    "cost_pct",
)


def _cells(row: CostRow) -> tuple[str, ...]:
    """One row's display strings, numbers fixed to a sane precision."""
    return (
        row.plugin,
        row.param,
        f"{row.base_s:.3f}",
        f"{row.plugin_s:.3f}",
        f"{row.cost_s:+.3f}",
        f"{row.cost_rss_mb:+.1f}",
        f"{row.cost_disk_mb:.1f}",
        f"{row.cost_pct:+.1f}%",
    )


def render_markdown(rows: Sequence[CostRow]) -> str:
    """A GitHub-flavoured markdown table of the cost rows."""
    lines = [
        "| " + " | ".join(_HEADERS) + " |",
        "| " + " | ".join("---" for _ in _HEADERS) + " |",
    ]
    lines += ["| " + " | ".join(_cells(r)) + " |" for r in rows]
    return "\n".join(lines) + "\n"


def render_pdf(rows: Sequence[CostRow], out: Path | str) -> Path:
    """Render the cost table to a one-page PDF at ``out``; return its path.

    Imported lazily so the module (and its markdown path) load without
    matplotlib on a bare machine.
    """
    import matplotlib

    matplotlib.use("Agg")  # headless: no display needed for a file render
    import matplotlib.pyplot as plt

    out = Path(out)
    cells = [list(_cells(r)) for r in rows] or [["—"] * len(_HEADERS)]

    fig_h = 1.2 + 0.4 * max(len(cells), 1)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.axis("off")
    ax.set_title("vepyr per-plugin performance cost (with-plugin − baseline)", pad=12)

    table = ax.table(
        cellText=cells,
        colLabels=list(_HEADERS),
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.4)
    for col in range(len(_HEADERS)):
        table[0, col].set_facecolor("#2b2b2b")
        table[0, col].set_text_props(color="white", fontweight="bold")

    fig.tight_layout()
    fig.savefig(out, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# CLI: the gate sbatch runs ``python -m perf_report`` after the driver, globbing
# every ``perf_*.json`` in a directory into one markdown + one PDF.
# ---------------------------------------------------------------------------


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="perf_report",
        description="Fold every perf_*.json in a directory into one markdown + PDF cost table.",
    )
    parser.add_argument(
        "--in", dest="in_dir", required=True, type=Path, help="Directory of perf_*.json files."
    )
    parser.add_argument("--md", required=True, type=Path, help="Markdown output path.")
    parser.add_argument("--pdf", required=True, type=Path, help="PDF output path.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Build the cost table from ``--in`` and write both ``--md`` and ``--pdf``."""
    args = _parse_args(argv)
    paths = sorted(args.in_dir.glob("perf_*.json"))
    if not paths:
        raise SystemExit(f"no perf_*.json under {args.in_dir}")
    rows = build_cost_table(paths)
    args.md.parent.mkdir(parents=True, exist_ok=True)
    args.md.write_text(render_markdown(rows))
    args.pdf.parent.mkdir(parents=True, exist_ok=True)
    render_pdf(rows, args.pdf)
    print(f"[perf] report: {len(rows)} rows -> {args.md}, {args.pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

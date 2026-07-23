#!/usr/bin/env python3
"""The per-plugin perf driver: what does adding a vepyr plugin cost?

For one ``(plugin, param)`` pair it builds the plugin cache from a PR checkout,
then annotates the same region twice — ``baseline`` (``plugin_cache_root=None``)
and ``+plugin`` (``plugin_cache_root=<built cache>``) — measuring elapsed time,
peak RSS and disk for each. The reported *cost* is ``with-plugin − baseline`` per
metric, so the plugin's own overhead is isolated from the engine's baseline.

Modelled on ``harness/parity.py::run_check``: same manifest resolution
(``plugins_repo`` = the PR checkout, at git ``version``), same mini-cache used as
both the tiering source for ``build_plugin_cache`` and the annotation
``cache_dir``, same ``reference_fasta`` discovery.

``vepyr`` is imported lazily inside :func:`run_perf` so this module (and its pure
unit tests) load on a bare machine without the native toolkit.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Final

from perf import PerfSample, dir_bytes, file_bytes, measure, measure_call

# The two annotation profiles the benchmark reports, expressed as the exact
# keyword arguments ``vepyr.annotate`` accepts. Derived from ``annotate()``'s
# signature — the toolkit exposes ``everything`` and ``hgvs`` as named booleans,
# so no VEP-flag string parsing is needed. ``hgvs=True`` implies hgvsc + hgvsp
# (VEP ``--hgvs``); ``everything=True`` is the full 80-field CSQ. Both require a
# ``reference_fasta``, which :func:`run_perf` always supplies from the mini-cache.
PARAM_FLAGS: Final[dict[str, dict[str, bool]]] = {
    "everything": {"everything": True},
    "hgvs": {"hgvs": True},
}

_FASTA_NAME: Final[str] = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"


# ---------------------------------------------------------------------------
# Result model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CostDelta:
    """The marginal cost of the plugin: with-plugin minus baseline, per metric."""

    elapsed_s: float
    peak_rss_bytes: int
    disk_bytes: int

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class PerfResult:
    """One ``(plugin, param)`` benchmark: the three runs and their delta."""

    plugin: str
    param: str
    build: PerfSample
    baseline: PerfSample
    withplugin: PerfSample
    cache_bytes: int
    baseline_out_bytes: int
    withplugin_out_bytes: int
    cost: CostDelta

    @classmethod
    def from_samples(
        cls,
        plugin: str,
        param: str,
        build: PerfSample,
        baseline: PerfSample,
        withplugin: PerfSample,
        cache_bytes: int,
        baseline_out_bytes: int,
        withplugin_out_bytes: int,
    ) -> "PerfResult":
        """Assemble a result, computing ``cost`` = with-plugin − baseline.

        The time delta is kept signed: measurement noise can make the plugin run
        marginally faster, and clamping that to zero would quietly bias the
        report. Disk cost is the plugin cache itself — the artefact only the
        ``+plugin`` run reads, so it is the plugin's marginal on-disk footprint.
        """
        cost = CostDelta(
            elapsed_s=round(withplugin.elapsed_s - baseline.elapsed_s, 6),
            peak_rss_bytes=withplugin.peak_rss_bytes - baseline.peak_rss_bytes,
            disk_bytes=cache_bytes,
        )
        return cls(
            plugin=plugin,
            param=param,
            build=build,
            baseline=baseline,
            withplugin=withplugin,
            cache_bytes=cache_bytes,
            baseline_out_bytes=baseline_out_bytes,
            withplugin_out_bytes=withplugin_out_bytes,
            cost=cost,
        )

    def to_json_dict(self) -> dict[str, object]:
        """The on-disk shape the report (:mod:`perf_report`) consumes."""
        return {
            "plugin": self.plugin,
            "param": self.param,
            "build": self.build.as_dict(),
            "baseline": self.baseline.as_dict(),
            "withplugin": self.withplugin.as_dict(),
            "cache_bytes": self.cache_bytes,
            "baseline_out_bytes": self.baseline_out_bytes,
            "withplugin_out_bytes": self.withplugin_out_bytes,
            "cost": self.cost.as_dict(),
        }

    def write_json(self, out_dir: Path | str) -> Path:
        """Write ``perf_<plugin>_<param>.json`` into ``out_dir``; return its path."""
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"perf_{self.plugin}_{self.param}.json"
        path.write_text(json.dumps(self.to_json_dict(), indent=2) + "\n")
        return path


# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------


def _chrom_of(region_vcf: Path) -> str:
    """The contig of the region's first data record (for ``chroms=[...]``)."""
    with region_vcf.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                return line.split("\t", 1)[0]
    raise ValueError(f"{region_vcf}: no data records to infer a chromosome from")


def run_perf(
    plugin: str,
    param: str,
    region_vcf: Path | str,
    mini_cache: Path | str,
    plugins_repo: Path | str,
    source_path: Path | str,
    out_dir: Path | str,
    *,
    version: str = "HEAD",
    overwrite: bool = True,
) -> PerfResult:
    """Benchmark one plugin under one param profile; write + return the result.

    Args:
        plugin: Plugin name, e.g. ``"alphamissense"``.
        param: ``"everything"`` or ``"hgvs"`` (see :data:`PARAM_FLAGS`).
        region_vcf: The biallelic region to annotate.
        mini_cache: The variation cache — tiering source for the plugin build and
            the annotation ``cache_dir`` both, mirroring ``parity.run_check``.
        plugins_repo: The vepyr-plugins checkout the manifest is resolved from at
            git ``version`` (the PR checkout, so a PR benchmarks what it ships).
        source_path: The plugin's source data slice (e.g. the AlphaMissense TSV).
        out_dir: Where the plugin cache, both output VCFs and the JSON are written.
        version: Git revision to materialise the manifest at (default ``HEAD``).
        overwrite: Rebuild an existing plugin cache (default ``True``).

    Returns:
        The :class:`PerfResult`, also serialised to
        ``out_dir/perf_<plugin>_<param>.json``.

    Raises:
        ValueError: ``param`` is not one of :data:`PARAM_FLAGS`.
    """
    if param not in PARAM_FLAGS:
        raise ValueError(
            f"unknown param {param!r}; expected one of {sorted(PARAM_FLAGS)}"
        )

    import vepyr  # noqa: PLC0415 - heavy native import; only the real run needs it

    region_vcf = Path(region_vcf).resolve()
    mini_cache = Path(mini_cache).resolve()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    chrom = _chrom_of(region_vcf)
    fasta = mini_cache / _FASTA_NAME
    flags = PARAM_FLAGS[param]

    # 1. Build the plugin cache — measured (this is the plugin's build/disk cost).
    plugin_cache = out_dir / f"plugin_cache_{plugin}"
    with measure() as build:
        vepyr.build_plugin_cache(
            plugin,
            version,
            source_path=str(source_path),
            cache_dir=str(mini_cache),
            plugin_cache_root=str(plugin_cache),
            chroms=[chrom],
            plugins_repo=str(plugins_repo),
            overwrite=overwrite,
        )

    # 2. Baseline: no plugin cache attached. Run in its OWN process so the peak
    #    RSS is this run's alone — ru_maxrss is a per-process monotonic
    #    high-water mark, so sharing a process with the build or the +plugin run
    #    would silently fold their peaks into the baseline's.
    base_out = out_dir / f"{plugin}_{param}_baseline.vcf"
    _, baseline = measure_call(
        vepyr.annotate,
        str(region_vcf),
        str(mini_cache),
        reference_fasta=str(fasta),
        plugin_cache_root=None,
        output_vcf=str(base_out),
        show_progress=False,
        **flags,
    )

    # 3. +plugin: the same run with the built plugin cache attached, likewise in
    #    a fresh process so baseline and +plugin peaks are independently honest.
    plug_out = out_dir / f"{plugin}_{param}_withplugin.vcf"
    _, withplugin = measure_call(
        vepyr.annotate,
        str(region_vcf),
        str(mini_cache),
        reference_fasta=str(fasta),
        plugin_cache_root=str(plugin_cache),
        output_vcf=str(plug_out),
        show_progress=False,
        **flags,
    )

    result = PerfResult.from_samples(
        plugin,
        param,
        build,
        baseline,
        withplugin,
        cache_bytes=dir_bytes(plugin_cache),
        baseline_out_bytes=file_bytes(base_out),
        withplugin_out_bytes=file_bytes(plug_out),
    )
    result.write_json(out_dir)
    return result

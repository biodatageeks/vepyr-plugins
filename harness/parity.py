#!/usr/bin/env python3
"""The parity gate: does a vepyr plugin manifest reproduce Ensembl VEP exactly?

Porting a VEP plugin to vepyr is writing a TOML manifest. The manifest is only
worth anything if it produces the *same annotations as the Perl plugin*, so this
harness is the gate that decides whether a port is real:

    build the plugin cache -> annotate a region -> diff the plugin's CSQ fields
    against real Ensembl VEP output -> require 100%.

Two modes:

``--check`` (hermetic; what CI runs)
    Builds, annotates and compares against a committed golden. No Perl, no Rust,
    no 34 GB cache.

``--refresh-golden`` (local; the only place Perl runs)
    Runs real Ensembl VEP with ``--plugin <X>`` and commits the result.

The blame-attribution rule
--------------------------
Plugin CSQ fields are *derived from core engine attributes* — AlphaMissense's
row discriminator is ``{ref_aa}{Protein_position}{alt_aa}``. If vepyr's core
disagrees with VEP about the transcript or the amino-acid change, the plugin
field comes out wrong through no fault of the manifest, and a naive diff blames
the plugin. The core *does* have known divergences from VEP, so this is not
hypothetical; the predictable outcome of a naive gate is someone "fixing" parity
by loosening it.

So the gate is computed from *two* comparisons over the same pair of files:

1. **Core agreement** over :data:`CORE_FIELDS` — the attributes the discriminator
   is built from. This is not the gate. It only determines the set of variant
   keys on which vepyr and VEP already agree, using the comparator's complete,
   uncapped :attr:`~vepyr.parity.ComparisonResult.mismatch_keys`.
2. **The port gate** over the plugin's own CSQ fields, restricted to that agreed
   subset. A mismatch that survives the subtraction is unambiguously the
   manifest's or the plugin cache's fault.

Keys excluded for core drift are reported loudly and separately — never folded
into zero. A rising exclusion rate is a signal about the *core*, and it must stay
visible. ``--strict`` makes core drift fail the build too: the mode for PRs
against the engine.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import tempfile
import tomllib
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import Final, Self

from vepyr.parity import ComparisonResult, compare_csq_fields

REPO_ROOT: Final[Path] = Path(__file__).resolve().parent.parent

CORE_FIELDS: Final[tuple[str, ...]] = (
    "Feature",
    "Consequence",
    "Amino_acids",
    "Protein_position",
)
"""The CSQ fields a plugin's row discriminator is derived from.

Deliberately *not* ``ref``/``alt``: those are not CSQ fields in this comparator
at all but part of the variant key, so a disagreement there cannot pair and
surfaces as ``keys_only_in_*``. Asking for them as ``fields=`` would fail
unclean with a spurious ``fields_missing_from_*``.
"""

_UNCAPPED: Final[int] = 1 << 30
"""``max_examples`` large enough that the examples are exhaustive.

:attr:`~vepyr.parity.ComparisonResult.mismatch_keys` is complete by construction,
but ``over_emissions`` is only a *count*; the keys behind it live in the capped
examples. Uncapping is how over-emission is made key-attributable, so that it can
be subjected to the same blame rule as any other mismatch.
"""

_REGION_VCF: Final[Path] = REPO_ROOT / "harness" / "regions" / "chr22-22.0-23.5Mb.vcf"


class HarnessError(RuntimeError):
    """The harness cannot render a verdict — a misconfiguration, not a result."""


class BuildHealthError(HarnessError):
    """The plugin cache built, but is dead on arrival."""


class Status(StrEnum):
    """The three things a parity run can conclude."""

    PASS = "PASS"
    FAIL = "FAIL"
    SKIP = "SKIP"


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ParityConfig:
    """A plugin's ``parity.toml``: what to build, what to compare, what to trust."""

    path: Path
    plugin: str
    csq_fields: tuple[str, ...]
    region: str
    redistributable: bool
    version: str
    vep_release: str
    vep_plugin_args: str
    source_url: str | None
    source_sha256: str | None

    @classmethod
    def load(cls, path: str | Path) -> Self:
        """Parse a ``parity.toml``, failing loudly on a missing required key."""
        path = Path(path).resolve()
        with path.open("rb") as fh:
            raw = tomllib.load(fh)

        vep = raw.get("vep", {})
        source = raw.get("source", {})
        try:
            return cls(
                path=path,
                plugin=raw["plugin"],
                csq_fields=tuple(raw["csq_fields"]),
                region=raw["region"],
                redistributable=bool(raw["redistributable"]),
                version=raw.get("version", "HEAD"),
                vep_release=str(vep.get("release", "")),
                vep_plugin_args=vep.get("plugin_args", ""),
                source_url=source.get("url"),
                source_sha256=source.get("sha256"),
            )
        except KeyError as exc:
            raise HarnessError(f"{path}: missing required key {exc}") from exc

    @property
    def plugin_dir(self) -> Path:
        """The plugin's directory — the manifest, fixtures and golden live here."""
        return self.path.parent

    @property
    def golden(self) -> Path:
        """Where ``--refresh-golden`` writes and ``--check`` reads."""
        return self.plugin_dir / "golden" / f"{self.plugin}.vcf"

    @property
    def chrom(self) -> str:
        """The region's contig, e.g. ``chr22``."""
        return self.region.split(":", 1)[0]

    @property
    def span(self) -> tuple[int, int]:
        """The region's 1-based inclusive start and end."""
        start, end = self.region.split(":", 1)[1].split("-", 1)
        return int(start), int(end)


def resolve_source(cfg: ParityConfig, *, override: str | Path | None) -> Path | None:
    """Locate the plugin's source data, or decide that this run must be skipped.

    Returns ``None`` — and *says so on stdout* — when the data is licence-gated
    and absent. A licence-gated plugin that quietly goes green is the worst of
    both worlds, so the skip is never silent.

    Raises:
        FileNotFoundError: The data is redistributable, so a missing fixture is a
            broken port rather than a licensing fact.
    """
    if override is not None:
        path = Path(override).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"--source-path does not exist: {path}")
        return path

    fixtures = sorted((cfg.plugin_dir / "fixtures").glob(f"{cfg.plugin}.*"))
    if fixtures:
        return fixtures[0]

    if cfg.redistributable:
        raise FileNotFoundError(
            f"{cfg.plugin}: redistributable = true, but no fixture at "
            f"{cfg.plugin_dir / 'fixtures'}/{cfg.plugin}.* and no --source-path. "
            "A redistributable plugin must ship the region slice it is tested on."
        )

    print(
        f"\n{'=' * 78}\n"
        f"  SKIP  {cfg.plugin}: source data is not redistributable\n"
        f"{'=' * 78}\n"
        f"  parity.toml sets redistributable = false, so the region slice cannot\n"
        f"  be committed to a public repo and CI has nothing to build from.\n"
        f"  This plugin is NOT covered by this run. To gate it, re-run on a machine\n"
        f"  that holds the data:\n\n"
        f"      harness/parity.py --check {cfg.plugin} --source-path <the data>\n"
        f"{'=' * 78}\n"
    )
    return None


# ---------------------------------------------------------------------------
# Build health
# ---------------------------------------------------------------------------


def assess_build_health(rows: Sequence[tuple[str, int, int, int]]) -> None:
    """Fail fast if the plugin cache built but joined nothing.

    ``warm == 0`` with ``rows > 0`` means not one source row matched the variation
    cache. The cache is dead; annotating with it would emit empty plugin fields
    everywhere and the diff would be a confusing wall of empties rather than one
    clear sentence. Say the sentence.

    Raises:
        BuildHealthError: Some chromosome ingested rows but warmed none of them.
    """
    for chrom, n_rows, warm, _cold in rows:
        if n_rows > 0 and warm == 0:
            raise BuildHealthError(
                f"{chrom}: {n_rows:,} source rows ingested but warm = 0 — not one "
                "row joined the variation cache. The plugin cache is dead; the "
                "manifest's key (chrom/start/allele_string) almost certainly does "
                "not match the cache's."
            )


# ---------------------------------------------------------------------------
# The verdict
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class BlameReport:
    """The gate's verdict, with core drift held apart from plugin failure."""

    plugin: str
    core: ComparisonResult
    gate: ComparisonResult
    excluded_keys: frozenset[str]
    surviving_mismatch_keys: dict[str, frozenset[str]]
    surviving_over_emission_keys: dict[str, frozenset[str]]
    excluded_over_emission_keys: dict[str, frozenset[str]]

    @property
    def surviving_mismatch_count(self) -> int:
        """Field-level mismatches that no core divergence can explain away."""
        return sum(len(keys) for keys in self.surviving_mismatch_keys.values())

    @property
    def surviving_over_emission_count(self) -> int:
        """vepyr populated a field VEP left empty, on a key where the core agreed."""
        return sum(len(keys) for keys in self.surviving_over_emission_keys.values())

    @property
    def hard_failures(self) -> tuple[str, ...]:
        """Everything that fails the gate regardless of ``--strict``."""
        problems: list[str] = []
        if self.gate.fields_missing_from_truth:
            problems.append(
                "the golden does not carry "
                f"{', '.join(self.gate.fields_missing_from_truth)} — it was produced "
                "WITHOUT --plugin, so it cannot gate anything"
            )
        if self.gate.fields_missing_from_test:
            problems.append(
                "vepyr's output does not carry "
                f"{', '.join(self.gate.fields_missing_from_test)} — the plugin cache "
                "was not attached, or the manifest emits different csq_field names"
            )
        if self.surviving_over_emission_count:
            problems.append(
                f"{self.surviving_over_emission_count} over-emission(s): vepyr filled "
                "a field VEP left empty, on keys where the core agreed"
            )
        residual = self.surviving_mismatch_count - self.surviving_over_emission_count
        if residual:
            problems.append(f"{residual} value mismatch(es) on core-agreeing keys")
        if self.gate.csq_missing_in_test:
            problems.append(
                f"{self.gate.csq_missing_in_test} variant(s) VEP annotated but vepyr "
                "left with no CSQ at all"
            )
        if self.gate.csq_missing_in_truth:
            problems.append(
                f"{self.gate.csq_missing_in_truth} variant(s) vepyr annotated but VEP "
                "left with no CSQ at all"
            )
        return tuple(problems)

    @property
    def core_drift(self) -> tuple[str, ...]:
        """Divergences that are the *engine's*, not the manifest's."""
        signals: list[str] = []
        if self.excluded_keys:
            signals.append(f"{len(self.excluded_keys)} key(s) excluded: core drift")
        if excluded_oe := sum(len(v) for v in self.excluded_over_emission_keys.values()):
            signals.append(
                f"{excluded_oe} of the excluded mismatches were over-emissions "
                "(a core divergence made vepyr match a row VEP could not)"
            )
        if self.core.entry_count_mismatch:
            signals.append(
                f"{self.core.entry_count_mismatch} variant(s) where vepyr and VEP "
                "emitted a different NUMBER of CSQ entries"
            )
        if self.core.entry_order_mismatch:
            signals.append(
                f"{self.core.entry_order_mismatch} variant(s) with a CSQ entry-order "
                "difference"
            )
        if self.core.keys_only_in_truth:
            signals.append(
                f"{self.core.keys_only_in_truth} variant(s) only VEP emitted"
            )
        if self.core.keys_only_in_test:
            signals.append(
                f"{self.core.keys_only_in_test} variant(s) only vepyr emitted"
            )
        return tuple(signals)

    @property
    def is_clean(self) -> bool:
        """True when nothing survives the blame rule."""
        return not self.hard_failures

    @property
    def status(self) -> Status:
        """The default (non-strict) verdict."""
        return self.status_under(strict=False)

    def status_under(self, *, strict: bool) -> Status:
        """The verdict, optionally holding the core to the same standard."""
        if self.hard_failures:
            return Status.FAIL
        if strict and self.core_drift:
            return Status.FAIL
        return Status.PASS

    def render(self, *, strict: bool = False, examples: int = 5) -> str:
        """The human-readable verdict. Core drift is never folded into zero."""
        out: list[str] = []
        add = out.append
        status = self.status_under(strict=strict)

        add("")
        add("=" * 78)
        add(f"  PARITY GATE  {self.plugin}")
        add("=" * 78)

        add("")
        add(f"  variants compared      : {self.gate.keys_compared:,}")
        add(f"  CSQ entries compared   : {self.gate.field_totals.get(next(iter(self.gate.fields_compared), ''), 0):,}")
        add(f"  plugin fields gated    : {', '.join(self.gate.fields_compared) or '(none)'}")

        add("")
        add("  --- the gate: plugin fields on keys where the core AGREES with VEP")
        for f in self.gate.fields_compared:
            total = self.gate.field_totals.get(f, 0)
            surviving = len(self.surviving_mismatch_keys.get(f, frozenset()))
            over = len(self.surviving_over_emission_keys.get(f, frozenset()))
            verdict = "ok" if surviving == 0 else f"FAIL ({surviving} bad keys)"
            note = f", {over} over-emission" if over else ""
            add(f"      {f:<24} {total - surviving:,}/{total:,} entries agree  [{verdict}{note}]")

        for f, keys in self.surviving_mismatch_keys.items():
            if not keys:
                continue
            add("")
            add(f"      {f}: {len(keys)} unexplained mismatch(es)")
            for ex in self.gate.field_mismatch_examples.get(f, [])[:examples]:
                if ex.key in keys:
                    tag = "  <- OVER-EMISSION" if ex.is_over_emission else ""
                    add(
                        f"        {ex.key.replace(chr(9), ' ')}: "
                        f"vep={ex.truth!r} vepyr={ex.test!r}{tag}"
                    )

        add("")
        add("  --- excluded: core drift (reported, NOT folded into the gate)")
        if not self.core_drift:
            add("      none — vepyr's core agreed with VEP on every compared key")
        for signal in self.core_drift:
            add(f"      {signal}")
        if self.excluded_keys:
            add("")
            add("      what the core disagreed about:")
            for f in self.core.fields_compared:
                exs = [
                    ex
                    for ex in self.core.field_mismatch_examples.get(f, [])
                    if ex.key in self.excluded_keys
                ]
                for ex in exs[:examples]:
                    add(
                        f"        {f} @ {ex.key.replace(chr(9), ' ')}: "
                        f"vep={ex.truth!r} vepyr={ex.test!r}"
                    )

        for problem in self.hard_failures:
            add("")
            add(f"  !! {problem}")

        add("")
        add("-" * 78)
        if status is Status.PASS and self.core_drift and not strict:
            add(f"  {status}  (plugin parity clean; core drift present and reported above)")
        else:
            add(f"  {status}")
        add("-" * 78)
        add("")
        return "\n".join(out)


def _over_emission_keys(result: ComparisonResult) -> dict[str, frozenset[str]]:
    """Per field, the keys where vepyr emitted a value and VEP emitted none.

    Only sound because the comparison was run uncapped (:data:`_UNCAPPED`): the
    over-emission *count* is public, but its keys live in the examples.
    """
    keyed: dict[str, frozenset[str]] = {}
    for f in result.fields_compared:
        keyed[f] = frozenset(
            ex.key
            for ex in result.field_mismatch_examples.get(f, [])
            if ex.is_over_emission
        )
        if len(keyed[f]) != result.over_emissions.get(f, 0):
            raise HarnessError(
                f"{f}: {result.over_emissions.get(f, 0)} over-emissions counted but "
                f"{len(keyed[f])} keyed — the examples were capped, so the blame rule "
                "cannot be applied soundly."
            )
    return keyed


def evaluate(
    golden_vcf: str | Path,
    test_vcf: str | Path,
    csq_fields: Sequence[str],
    *,
    plugin: str = "",
    ignore_entry_order: bool = False,
) -> BlameReport:
    """Run both comparisons and apply the blame-attribution rule.

    Args:
        golden_vcf: Real Ensembl VEP output — the truth side.
        test_vcf: vepyr's output with the plugin cache attached — the test side.
        csq_fields: The plugin's own CSQ fields. These, and only these, are gated.
        plugin: Name, for the report header.
        ignore_entry_order: Forwarded to the comparator.

    Returns:
        The verdict, with the core-drift exclusion set held separately.

    Raises:
        HarnessError: The golden lacks a core field, which makes the exclusion
            set uncomputable — the blame rule would silently degrade to a naive
            diff, so refuse rather than mislead.
    """
    core = compare_csq_fields(
        golden_vcf,
        test_vcf,
        CORE_FIELDS,
        ignore_entry_order=ignore_entry_order,
        max_examples=_UNCAPPED,
    )
    if core.fields_missing_from_truth or core.fields_missing_from_test:
        missing = set(core.fields_missing_from_truth) | set(core.fields_missing_from_test)
        raise HarnessError(
            f"core fields absent from the compared VCFs: {', '.join(sorted(missing))}. "
            "Without them the core-agreement set cannot be computed and the gate "
            "would degrade to a naive diff that blames the plugin for core bugs."
        )

    # THE exclusion set: every key on which vepyr's core already disagrees with
    # VEP about the transcript or the amino-acid change the plugin keys off.
    excluded: frozenset[str] = frozenset().union(*core.mismatch_keys.values()) if core.mismatch_keys else frozenset()

    gate = compare_csq_fields(
        golden_vcf,
        test_vcf,
        csq_fields,
        ignore_entry_order=ignore_entry_order,
        max_examples=_UNCAPPED,
    )
    gate_over = _over_emission_keys(gate)

    surviving = {f: frozenset(keys) - excluded for f, keys in gate.mismatch_keys.items()}
    surviving_over = {f: keys - excluded for f, keys in gate_over.items()}
    excluded_over = {f: keys & excluded for f, keys in gate_over.items()}

    return BlameReport(
        plugin=plugin,
        core=core,
        gate=gate,
        excluded_keys=excluded,
        surviving_mismatch_keys=surviving,
        surviving_over_emission_keys=surviving_over,
        excluded_over_emission_keys=excluded_over,
    )


# ---------------------------------------------------------------------------
# --check : build, annotate, compare. No Perl.
# ---------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    """Content digest, for provenance."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def run_check(
    cfg: ParityConfig,
    *,
    mini_cache: Path,
    source_override: str | Path | None,
    region_vcf: Path,
    strict: bool,
    workdir: Path,
) -> Status:
    """Build the plugin cache, annotate the region, and gate the result."""
    import vepyr  # noqa: PLC0415 — heavy native import; only --check needs it

    source = resolve_source(cfg, override=source_override)
    if source is None:
        return Status.SKIP

    if not cfg.golden.exists():
        raise HarnessError(
            f"{cfg.plugin}: no golden at {cfg.golden}. A gate without a golden is "
            "not a gate. Generate one on a machine with Ensembl VEP:\n"
            f"    harness/parity.py --refresh-golden {cfg.plugin} --source-path <data>"
        )

    print(f"[{cfg.plugin}] building plugin cache from {source.name} ...")
    plugin_cache = workdir / "plugin_cache"
    rows = vepyr.build_plugin_cache(
        cfg.plugin,
        cfg.version,
        source_path=str(source),
        cache_dir=str(mini_cache),
        plugin_cache_root=str(plugin_cache),
        chroms=[cfg.chrom],
        # The manifest under test is the one in THIS working tree — never one
        # fetched from a tag, or the PR would gate a different file than it ships.
        plugins_repo=str(REPO_ROOT),
        overwrite=True,
    )
    for chrom, n_rows, warm, cold in rows:
        print(f"[{cfg.plugin}]   {chrom}: rows={n_rows:,} warm={warm:,} cold={cold:,}")
    assess_build_health(rows)

    print(f"[{cfg.plugin}] annotating {region_vcf.name} with the plugin cache ...")
    ours = workdir / f"{cfg.plugin}.vepyr.vcf"
    vepyr.annotate(
        str(region_vcf),
        str(mini_cache),
        reference_fasta=str(mini_cache / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        plugin_cache_root=str(plugin_cache),
        output_vcf=str(ours),
        show_progress=False,
    )

    report = evaluate(cfg.golden, ours, cfg.csq_fields, plugin=cfg.plugin)
    print(report.render(strict=strict))
    return report.status_under(strict=strict)


# ---------------------------------------------------------------------------
# --refresh-golden : the only place Perl runs.
# ---------------------------------------------------------------------------


def _vep_plugin_version(dir_plugins: Path, plugin_pm: str) -> str:
    """Identify the Perl plugin that produced a golden.

    A golden whose provenance is unknown is not a golden, it is a rumour.
    """
    try:
        rev = subprocess.run(
            ["git", "-C", str(dir_plugins), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        branch = subprocess.run(
            ["git", "-C", str(dir_plugins), "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        return f"VEP_plugins@{branch}:{rev}:{plugin_pm}"
    except (subprocess.CalledProcessError, FileNotFoundError):
        return f"VEP_plugins@unknown:{plugin_pm}"


def run_refresh_golden(
    cfg: ParityConfig,
    *,
    vep: Path,
    vep_cache: Path,
    fasta: Path,
    dir_plugins: Path,
    source_override: str | Path | None,
    region_vcf: Path,
    workdir: Path,
) -> Status:
    """Run real Ensembl VEP with ``--plugin`` and write the golden."""
    source = resolve_source(cfg, override=source_override)
    if source is None:
        raise HarnessError(
            f"{cfg.plugin}: cannot refresh a golden without the source data. "
            "Pass --source-path."
        )

    for needed, what in ((vep, "the vep script"), (vep_cache, "the VEP cache"),
                         (dir_plugins, "the VEP_plugins checkout")):
        if not needed.exists():
            raise HarnessError(
                f"cannot run Ensembl VEP: {what} is missing at {needed}. "
                "Refusing to fabricate a golden — a fabricated golden would "
                "silently validate a broken port forever."
            )

    # `{source}` in plugin_args is substituted with the resolved data path, so the
    # committed parity.toml carries no machine-specific path.
    plugin_args = cfg.vep_plugin_args.replace("{source}", str(source))
    raw = workdir / "vep_raw.vcf"
    cmd = [
        str(vep),
        "--input_file", str(region_vcf),
        "--output_file", str(raw),
        "--vcf",
        "--cache",
        "--offline",
        "--dir_cache", str(vep_cache),
        "--assembly", "GRCh38",
        "--fasta", str(fasta),
        "--no_stats",
        "--force_overwrite",
        "--dir_plugins", str(dir_plugins),
        "--plugin", plugin_args,
    ]
    print(f"[{cfg.plugin}] $ {' '.join(cmd)}\n")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise HarnessError(
            "Ensembl VEP failed — NOT writing a golden.\n"
            f"  exit {proc.returncode}\n"
            f"  stdout: {proc.stdout[-2000:]}\n"
            f"  stderr: {proc.stderr[-2000:]}"
        )
    if proc.stderr.strip():
        print(f"[{cfg.plugin}] vep stderr:\n{proc.stderr.strip()}\n")

    plugin_pm = cfg.vep_plugin_args.split(",", 1)[0]
    provenance = [
        f"##vepyr_parity_generated={datetime.now(UTC).isoformat(timespec='seconds')}",
        f"##vepyr_parity_plugin={cfg.plugin}",
        f"##vepyr_parity_region={cfg.region}",
        f"##vepyr_parity_vep_release={cfg.vep_release}",
        f"##vepyr_parity_vep_plugin={_vep_plugin_version(dir_plugins, plugin_pm)}",
        f"##vepyr_parity_plugin_args={cfg.vep_plugin_args}",
        f"##vepyr_parity_source_file={Path(source).name}",
        f"##vepyr_parity_source_sha256={_sha256(Path(source))}",
        f"##vepyr_parity_region_vcf_sha256={_sha256(region_vcf)}",
    ]
    if cfg.source_url:
        provenance.append(f"##vepyr_parity_source_url={cfg.source_url}")

    lines = raw.read_text().splitlines()
    # Provenance goes directly after ##fileformat, where a reader looks first.
    head = 1 if lines and lines[0].startswith("##fileformat") else 0
    stamped = lines[:head] + provenance + lines[head:]

    cfg.golden.parent.mkdir(parents=True, exist_ok=True)
    cfg.golden.write_text("\n".join(stamped) + "\n")

    body = sum(1 for line in stamped if not line.startswith("#"))
    print(f"[{cfg.plugin}] wrote {cfg.golden} ({body:,} variants)")
    for p in provenance:
        print(f"[{cfg.plugin}]   {p}")
    return Status.PASS


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _discover(plugin: str | None) -> list[ParityConfig]:
    """Every ``plugins/*/parity.toml``, or just the named one."""
    if plugin:
        path = REPO_ROOT / "plugins" / plugin / "parity.toml"
        if not path.exists():
            raise HarnessError(f"no parity.toml for plugin {plugin!r} at {path}")
        return [ParityConfig.load(path)]
    return [ParityConfig.load(p) for p in sorted((REPO_ROOT / "plugins").glob("*/parity.toml"))]


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point. Exit code is the verdict: 0 pass or skip, 1 fail."""
    ap = argparse.ArgumentParser(
        prog="harness/parity.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="the hermetic gate (CI)")
    mode.add_argument(
        "--refresh-golden",
        action="store_true",
        help="run real Ensembl VEP to (re)generate the golden (local only)",
    )
    ap.add_argument("plugin", nargs="?", help="plugin name; default: every plugin")
    ap.add_argument(
        "--strict",
        action="store_true",
        help="core drift fails the build too (the mode for PRs against the engine)",
    )
    ap.add_argument(
        "--source-path",
        help="the plugin's source data; required when it is not redistributable",
    )
    ap.add_argument(
        "--mini-cache",
        default=os.environ.get("VEPYR_MINI_CACHE", "/tmp/mini_cache_region"),
        help="the region mini-cache [env VEPYR_MINI_CACHE]",
    )
    ap.add_argument("--region-vcf", default=str(_REGION_VCF))
    ap.add_argument(
        "--vep",
        default=os.environ.get("VEP", "/Users/wojtek/Documents/vepyr/ensembl-vep/vep"),
    )
    ap.add_argument(
        "--vep-cache",
        default=os.environ.get("VEP_CACHE", "/Users/wojtek/Documents/vepyr/_cache_v115"),
    )
    ap.add_argument(
        "--fasta",
        default=os.environ.get(
            "VEP_FASTA",
            "/Users/wojtek/Documents/vepyr/_cache_v115/"
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        ),
    )
    ap.add_argument(
        "--dir-plugins",
        default=os.environ.get("VEP_PLUGINS", "/Users/wojtek/Documents/vepyr/VEP_plugins"),
    )
    ap.add_argument("--keep-workdir", action="store_true")
    args = ap.parse_args(argv)

    try:
        configs = _discover(args.plugin)
    except HarnessError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    region_vcf = Path(args.region_vcf).resolve()
    statuses: dict[str, Status] = {}

    for cfg in configs:
        workdir = Path(tempfile.mkdtemp(prefix=f"parity_{cfg.plugin}_"))
        try:
            if args.check:
                statuses[cfg.plugin] = run_check(
                    cfg,
                    mini_cache=Path(args.mini_cache),
                    source_override=args.source_path,
                    region_vcf=region_vcf,
                    strict=args.strict,
                    workdir=workdir,
                )
            else:
                statuses[cfg.plugin] = run_refresh_golden(
                    cfg,
                    vep=Path(args.vep),
                    vep_cache=Path(args.vep_cache),
                    fasta=Path(args.fasta),
                    dir_plugins=Path(args.dir_plugins),
                    source_override=args.source_path,
                    region_vcf=region_vcf,
                    workdir=workdir,
                )
        except (HarnessError, FileNotFoundError) as exc:
            print(f"\n  FAIL  {cfg.plugin}: {exc}\n", file=sys.stderr)
            statuses[cfg.plugin] = Status.FAIL
        finally:
            if args.keep_workdir:
                print(f"[{cfg.plugin}] workdir kept at {workdir}")
            else:
                subprocess.run(["rm", "-rf", str(workdir)], check=False)

    if len(statuses) > 1:
        print("=" * 78)
        for name, status in statuses.items():
            print(f"  {status:<5} {name}")
        print("=" * 78)

    return 1 if Status.FAIL in statuses.values() else 0


if __name__ == "__main__":
    sys.exit(main())

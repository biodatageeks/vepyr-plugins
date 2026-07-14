# Hand-off — VEP plugins → vepyr: the plugin port factory

**Date:** 2026-07-14
**Owner:** sitekwb@gmail.com
**State:** Plan A (engine) and Plan B (toolkit) **merged**. Plan C (catalogue + parity gate) **in flight** — the gate is built, green on its first client, and has already caught a real engine bug.

---

## 0. Read these first (they hold the detail; this doc does not repeat it)

| Artifact | Where |
|---|---|
| **Design spec** (the factory: parity harness, mini-cache, blame-attribution rule) | `datafusion-bio-functions:master-sitekwb` → `docs/superpowers/specs/2026-07-13-vep-plugin-port-factory-design.md` |
| **Plan A** — engine (VCF provider + manifest hardening) | same repo → `docs/superpowers/plans/2026-07-13-plugin-factory-engine.md` |
| **Plan B** — toolkit (`vepyr.parity`, mini-cache) | `vepyr:master-sitekwb` → `docs/superpowers/plans/2026-07-14-plugin-factory-toolkit.md` |
| **Plan C** — catalogue (gate + 3 clients) | `vepyr-plugins:docs/plan-c` → `docs/plans/2026-07-14-plugin-factory-catalogue.md` |
| **Feasibility analysis** (all 78 plugins bucketed) | `datafusion-bio-functions:origin/dev-test-plugins-feasibility` → `docs/superpowers/specs/2026-07-10-vep-plugins-systematic-port-feasibility.md` |

**DO NOT trust the 2026-07-12 `HANDOFF veppluginsport.md`.** Four of its premises are false against the code (it claims a Python `build_plugin_cache` "already shipped" that didn't exist; its seed REVEL/PrimateAI manifests are wrong; it calls Bucket A "manifest only" while the VCF provider returned `NotImplemented`). Corrections are recorded in the spec's §1.

---

## 1. Branch policy (hard rule — the user set this mid-session)

**`master-sitekwb`** exists in all three repos, cut from `master`, and is treated as `main`.
Feature branches off it; PRs into it. **NEVER commit to `master`/`main`.**

Why not `dev-test`: it diverged from `master` at `v0.10.0` and both moved. `master` took #181 (rewrote `annotate_provider.rs`, deleted `variant_lookup_exec.rs`); `dev-test` has 45 commits including an engine feature that exists nowhere else. Merging them is a **feature port into a rewritten hot path**, not a sync — see §4.

---

## 2. What is merged

**Plan A** — `datafusion-bio-functions` PRs #194 (spec+plan), #195 (engine) → `master-sitekwb`.
`provider = "vcf"` wired to `VcfTableProvider` (it was `NotImplemented` while the feasibility analysis called Bucket A "manifest only"). Plus: unknown manifest keys rejected, `--overwrite` honoured *and parsed*, `--source-path`/`--chrom` fixed, `SourceManifest::validate()`, a build-time guard rejecting pipe-joined `allele_string`, a strict value-column cast, and a `warm == 0` warning.

**Plan B** — `vepyr` PR #33 → `master-sitekwb`.
`vepyr.parity` (the CSQ comparator, extracted from `e2e-testing/scripts/run_annotation_fast.py`, **proven behaviour-preserving** byte-for-byte), `scripts/build_mini_cache.py`, and the engine pin moved to the Plan A merge.

**Plan C** — `vepyr-plugins` PR **#1** (`feat/parity-harness` → `master-sitekwb`) — **open, needs review/merge.**
The harness (`--check` / `--refresh-golden`), `parity.toml`, and client 1 (AlphaMissense) green against a **real** Ensembl VEP golden.

---

## 3. Facts about the environment you will otherwise rediscover painfully

- **The mini-cache must be BUILT, not sliced.** The local variation cache (`_cache_v115/parquet/`) has **no `tier` column** (and no `chrom_manifest.json`, so the current runtime cannot even open it). `plugin_cache::join` selects `tier`. `vepyr/scripts/build_mini_cache.py` rebuilds the region with the current builder.
- **Region: `chr22:22,000,000–23,500,000`**, +5 kb flank. The flank is load-bearing: VEP annotates transcripts up to `--distance` (default 5000 bp) away, and a bare overlap filter dropped one 4,723 bp outside — which would have manufactured phantom "core drift" and poisoned the blame rule. The fixture test (full-cache vs mini-cache annotation must be **body-identical**) caught it. It passes now.
- Mini-cache: **14 MB** tarball, `sha256 e4617ac28750902e0e1ad56f01708e46673cdc75ed4d3f37203e37b9cb46db45`. **Not published** — awaiting user approval. Plan C's CI cannot fetch it until it is.
- **Real Ensembl VEP is now installed on this machine** (115.2 / API v115 / GENCODE 49, with BioPerl + `Bio::DB::HTS`; the install needed two fixes). It was not there before.
- The **VCF reader has no `pos` column** — POS is `start` (1-based, UInt32); `end` = `POS + len(REF) - 1` for indels. INFO keys are **bare and case-sensitive** (`` `AF` ``, not `info_af` — the crate's own docs lie). Multi-allelic ALTs arrive **pipe-joined** (`G|T`), so `ingest_sql` must `unnest(string_to_array(alt, '|'))`.
- **obelix** (HPC `$HOME`): 663 GB / 1 TB after this session's cleanup (57 GB reclaimed; summaries archived to `bvp/results/_archive/*.tar.gz`). `vepyr/data/vep_native_cache_116` (107 GB) and `ground_truth_vep_116` (41 GB) are the reference-VEP + golden data — **do not delete, Plan C stands on them.**

---

## 4. The engine bug the gate found — and why it changes `dev-test`'s status

**vepyr corrupts the CSQ record on multi-allelic variants** (on `master`/`master-sitekwb`):

1. It emits **one entry per transcript, not per (allele × transcript)** — VEP emitted exactly 2× vepyr's count on every multi-allelic site.
2. It writes the ALTs into `Allele` **joined by `|`** — the CSQ separator. Entries carry **77 tokens against their own 76-field header**, so any conformant parser reads `Consequence` as an allele and `Feature` as `"Transcript"`. Bi-allelic entries: exactly 76. Multi-allelic: *every* one is 77.

It did not affect the AlphaMissense gate (none of those variants are missense, and the blame rule correctly quarantined them), but **on real WGS it hits every multi-allelic site.**

**The fix already exists on `dev-test`**: `PerAltCtx` + `vcf_to_vep_allele_multi` — "multi-ALT CSQ per-allele expansion", which exists nowhere else. Decoupling it from the factory was right (it blocked nothing), but this **reclassifies it from "someone's parallel work" to "the fix for a live data-corruption bug on master"**. Port it, with a parity test on multi-allelic sites — the gate can now verify it.

⚠️ Someone was mid-merge on `sync/master-into-dev-test-2026-07-07` in the main `datafusion-bio-functions` worktree (uncommitted changes to `annotate_provider.rs`, `transcript_consequence.rs`). **Do not clobber it.**

---

## 5. Open decisions — blocked on the user

1. **Publish the mini-cache** as a release asset (14 MB). Without it Plan C's CI has nothing to fetch. *Outward-facing → needs explicit approval.*
2. **Data for the remaining two clients.** REVEL needs registration (non-commercial); the VCF client (Mastermind or gnomADMt) has its own terms. Which VCF plugin, and may data be downloaded?
3. **Enrich the region VCF with missense-dense variants.** The gate currently rests on **3 missense variants / 14 CSQ entries** — the other ~10,950 are a true-negative check. A negative control (corrupted source → gate failed with correct blame) proves it has teeth, but 3 positives is thin, and REVEL is also a missense scorer, so it inherits the weakness.

---

## 6. Known follow-ups (recorded, not filed — the user did not authorise issue creation)

**Engine (`datafusion-bio-functions`):**
- A wrong `[[match_column]]` **format** (e.g. `p.His101Tyr` vs `H101Y`) builds `warm=2, cold=0` — perfectly healthy — and annotates nothing. The variation join is on `(chrom, start, allele_string)` only and never touches match columns, so `warm == 0` is structurally blind to it. **The last silent-failure path.** Needs design, not a patch.
- `schema_matches` compares only value/match columns, so a filtered rebuild preserves chroms built from a **different** `ingest_sql`/`coordinate_system`. Hash the manifest into `CacheManifest`.
- `reject_pipe_joined_alleles` handles Utf8/Utf8View/LargeUtf8 and errors otherwise — but the `env_logger` init in `examples/build_plugin.rs` is still not pinned by an automated test.
- A corrupt `manifest.json` is silently treated as "no prior build".
- Hoist the bio-formats git tag into `[workspace.dependencies]` (the `-core`/`-vcf` tag pin is enforced only by a comment).

**Toolkit (`vepyr`):**
- `vepyr.annotate()` has **no off switch for regulatory/motif**. Goldens must use VEP's `--regulatory`, else the exclusion set inflates 6 → 102 — and excluded keys are exempt from the gate, so that is a place a real plugin bug could hide.

---

## 7. Suggested skills for the next session

- `superpowers:subagent-driven-development` — how Plans A/B/C were executed (fresh subagent per task, two-stage review). It works; keep it.
- `superpowers:test-driven-development`, `superpowers:verification-before-completion` — the load-bearing ones here.
- `superpowers:requesting-code-review` / fresh, context-free reviewers. **This is the single highest-value habit from this session:** independent reviewers who verified *by execution* rather than by reading caught three false claims of mine, a `--overwrite` fix that never reached the CLI, and a guard that failed open on `Utf8View`.
- `iipw-infrastructure` — for anything touching obelix/HPC.

## 8. The lesson worth carrying over

**The same defect recurred three times, each time one layer up, with green tests underneath.** `--overwrite` was fixed in the builder but never parsed by the CLI. The `warm == 0` warning was added but no logger was initialised. Then that same working warning was **muted through the Python API** because `env_logger` defaults to `error`. Each layer looked correct in isolation. When you add a guard, **prove it reaches a human on the path a user actually takes** — not that the function computing it returns the right string.

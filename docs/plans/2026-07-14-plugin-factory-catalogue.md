# Plugin Factory — Plan C: The Catalogue (the parity gate, and the three clients that prove it)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Make a PR that adds a plugin manifest **testable by itself** — CI builds the plugin's cache, annotates a region, and diffs against real Ensembl VEP output. 100% on the plugin's CSQ fields, or the PR is red.

**Architecture:** `vepyr-plugins` becomes data + configuration + goldens + one harness. It compiles no Rust. Its CI needs exactly three things: `pip install vepyr`, the mini-cache release asset, and the committed goldens. Real Ensembl VEP (Perl) runs **only** locally, only when a golden is (re)generated.

**Tech Stack:** Python (harness), GitHub Actions, TOML.

**Depends on (both merged/underway in Plan B):**
- `vepyr.build_plugin_cache(manifest_path, source_path, variation_cache_dir, plugin_cache_root, chroms=None, overwrite=False) -> list[(chrom, rows, warm, cold)]`
- `vepyr.parity.compare_csq_fields(truth_vcf, test_vcf, fields=None, *, ignore_entry_order=False) -> ComparisonResult`
- the region mini-cache release asset (`chr22:22,000,000–23,500,000`, **built** with the current builder — the old cache has no `tier` column)

**Source spec:** `datafusion-bio-functions/docs/superpowers/specs/2026-07-13-vep-plugin-port-factory-design.md` (§6, §10).

**Branch policy:** cut `master-sitekwb` from `master`, treat it as `main`, PR into it. Never commit to `master`.

---

## Prerequisite — one small addition to `vepyr.parity`

The blame-attribution rule (below) needs **the full set of variant keys where the core annotation disagrees with VEP** — not a sample of them. `ComparisonResult` currently exposes per-field *examples*, which are capped. Add to `vepyr.parity` (Plan B repo, its own PR):

```python
mismatch_keys: dict[str, set[str]]   # field -> every key that mismatched, uncapped
```

Without it the harness can only approximate the exclusion set, and an approximate blame rule is worse than none — it would quietly attribute core bugs to plugins on the long tail it cannot see.

---

## The rule this whole plan exists to enforce

Plugin CSQ fields are **derived from core engine attributes**: AlphaMissense's discriminator is
`{ref_aa}{Protein_position}{alt_aa}`. So if vepyr's core disagrees with VEP about the transcript or
the amino-acid change, the plugin field comes out wrong — and a naive diff blames the plugin. The
core *does* have known divergences (`vepyr/e2e-testing/reports/issue88_remaining_unresolved.md`), so
this is not hypothetical. The predictable outcome of a naive gate is someone "fixing" parity by
loosening it.

From the same pair of files, with no extra runs, compute **two verdicts**:

1. **Core agreement** — compare the CSQ fields the discriminator depends on: **`Feature`,
   `Consequence`, `Amino_acids`, `Protein_position`** — and NOT `ref`/`alt`, which are not CSQ
   fields at all in this comparator but part of the variant key, so a disagreement there cannot
   pair and surfaces as `keys_only_in_*`. Passing them as `fields=` would fail unclean with a
   spurious `fields_missing_from_*`. This is **not** the gate. It only
   determines the set of variant keys where vepyr and VEP already agree.
2. **The port gate** — compare the plugin's CSQ fields, **restricted to that agreed subset**. A
   mismatch there is unambiguously the manifest's or `plugin_cache`'s fault.

Lines excluded for core drift are reported **loudly and separately** (`excluded: core drift, N keys`),
never folded into zero. A rising exclusion rate is a signal about the core, and it must be visible.

`--strict` makes core drift fail the build too — the mode for PRs against the engine.

---

## File Structure

```
plugins/<name>/
  <name>.source.toml     # the manifest (exists for alphamissense)
  parity.toml            # NEW: csq_fields, region, VEP + plugin versions, redistributable
  fixtures/<name>.<ext>  # source-data slice for the region (only if redistributable)
  golden/<name>.vcf      # real VEP --plugin <X> output (written by --refresh-golden)
harness/parity.py        # NEW: the two-mode harness
harness/regions/chr22-22.0-23.5Mb.vcf   # NEW: the input VCF, committed (tens of KB)
.github/workflows/parity.yml            # NEW
```

### `parity.toml` schema

```toml
plugin          = "alphamissense"
csq_fields      = ["am_class", "am_pathogenicity"]   # THE gate's fields
region          = "chr22:22000000-23500000"
redistributable = false          # AlphaMissense is CC-BY-NC -> no fixture in the public repo

[vep]
release      = "115"
plugin_args  = "AlphaMissense,file=/data/AlphaMissense_hg38.tsv.gz"

[source]
url    = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
sha256 = "<of the full release, for provenance>"
```

`redistributable = false` ⇒ the plugin's CI job is **skipped with a loud, visible message** and runs
only in the nightly job on a machine that holds the data. Never a silent skip: a licence-gated plugin
that quietly goes green is the worst of both worlds.

---

### Task 1: The harness — `--check` (CI, hermetic)

**Files:** `harness/parity.py`, `harness/test_parity_harness.py`

- [ ] **Step 1: Write the failing tests first.** Over synthetic golden/vepyr VCF pairs in `tmp_path`, cover the semantics that actually matter:
  - a plugin-field mismatch on a **core-agreeing** key ⇒ FAIL;
  - the same mismatch on a key where the **core also disagrees** ⇒ EXCLUDED, reported, and the run still passes (unless `--strict`);
  - **over-emission** (vepyr fills a field VEP left empty) ⇒ FAIL — this was one of two real bugs the manual AlphaMissense run caught;
  - **`warm == 0` with `rows > 0`** from `build_plugin_cache` ⇒ FAIL before comparing anything. Not one row joined the variation cache; the cache is dead and the diff would be a confusing wall of empties rather than one clear sentence;
  - `redistributable = false` with no fixture ⇒ SKIP, and the skip is *printed*, not silent.

- [ ] **Step 2: Implement.** Flow:
  1. read `parity.toml`;
  2. `vepyr.build_plugin_cache(...)` against the mini-cache, `chroms=["chr22"]`;
  3. **assert `warm > 0`** — the free health signal (see Step 1);
  4. annotate the region VCF with `plugin_cache_root` set;
  5. `compare_csq_fields(golden, ours, fields=CORE_FIELDS)` → the excluded key set;
  6. `compare_csq_fields(golden, ours, fields=csq_fields)` → the gate, minus the excluded keys;
  7. exit non-zero on any surviving mismatch, over-emission, or `csq_missing_*`.

- [ ] **Step 3: Commit.**

---

### Task 2: The harness — `--refresh-golden` (local, the only place Perl lives)

- [ ] **Step 1:** Run real Ensembl VEP (`/Users/wojtek/Documents/vepyr/ensembl-vep/vep`, cache v115 on disk) with `--plugin <X>` over the region, writing `plugins/<name>/golden/<name>.vcf`.
- [ ] **Step 2:** Record in the golden's header **which VEP release and which plugin version produced it** — a golden whose provenance is unknown is not a golden, it is a rumour.
- [ ] **Step 3:** Emit the VCF **sites-only** (8 columns) if that is what the region run naturally produces — and note that `vepyr.parity` was fixed for exactly this shape (a trailing-newline bug that manufactured phantom mismatches on the **last** CSQ field, which is precisely where plugin fields are appended). Verify no phantom mismatch appears.
- [ ] **Step 4:** Commit the golden.

---

### Task 3: Client 1 — AlphaMissense (proves the HARNESS is correct)

The only manifest that exists, and the only one with an independently confirmed parity result. Used
**strictly as a fixture** — its manifest is owned elsewhere and is **not touched**.

- [ ] Generate the golden with `--refresh-golden`; run `--check`; require a clean gate.
- [ ] **Be honest in the report:** the manual run's "1,912/1,912" figure was full-chr22 and **cannot** be reproduced here. What is reproduced is *parity on the mini-cache region*. Say so; do not quietly imply otherwise.

---

### Task 4: Client 2 — REVEL (proves "a port is just a manifest" for someone who didn't write the engine)

- [ ] Author `plugins/REVEL/REVEL.source.toml` **from the AlphaMissense pattern** — not from the 2026-07-12 handoff's seed, which is wrong (`has_header = true` with no `schema` block; `provider.rs::csq_schema` builds the Arrow schema *solely* from that block, so an absent block yields a zero-field schema).
- [ ] Open decision to settle **against real data and the golden, not by preference**: ship REVEL as a pure per-variant plugin, or per-transcript with a `[[match_column]]` aa-gate? Run both, diff against VEP, keep the one that matches.

---

### Task 5: Client 3 — Mastermind or gnomADMt (proves the new VCF provider end-to-end)

A provider with no client is an unproven provider. This is the one that exercises Plan A's VCF wiring
for real.

- [ ] **The manifest MUST split multi-allelic ALTs.** The reader pipe-joins them (`alt = "G|T"`), so
  the natural `concat(ref, '/', alt)` produces `A/G|T`, which can never match the probe key. The engine
  now **rejects** this at build time rather than silently writing dead rows — so a naive manifest fails
  loudly. Use `unnest(string_to_array(alt, '|'))` so each ALT is its own row.
- [ ] INFO columns are **bare, case-sensitive** keys (`` `AF` ``, not `info_af`), and the reader has
  **no `pos` column** — VCF POS is `start` (1-based), with `end`. Both traps are documented on
  `VcfParams`; both will bite anyone who writes the manifest from memory.

---

### Task 6: CI

**Files:** `.github/workflows/parity.yml`

- [ ] Matrix over `plugins/*/parity.toml`. Per job: `pip install vepyr`; restore the mini-cache release asset via `actions/cache` keyed on its checksum; run `harness/parity.py --check <plugin>`.
- [ ] **No Rust toolchain. No Perl. No 34 GB.** If a step needs any of those, the design has slipped — stop and say so.
- [ ] A `redistributable = false` plugin prints a loud SKIP and exits 0.

---

## Definition of done

`parity --check` is green for three clients, each proving a different axis: **AlphaMissense** (the
harness is correct), **REVEL** (a newcomer can port a TSV plugin from a manifest alone), **a VCF
plugin** (the new provider works end-to-end). Only then does Wave 1 — the ~20 Bucket A/B scoring
plugins — become the pure manifest work the feasibility analysis promised.

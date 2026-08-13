# LCMS Consensus Feature Grouping — Design

**Date:** 2026-08-06 (updated 2026-08-13)  
**Status:** Living design for MR !248 / issue #272 — **Stage 2 adducts WIP with known correctness issues**  
**Branch:** `272_lcms_consensus_feature_grouping`  
**Approach:** A′ — RT ∩ chemical Δm edges → Pearson apex-height gate → roll-up / merge labels  

This document tracks **what we are building and review stages**. Implementation may land in **one MR**; review proceeds by stage.  
**Not committed to the library tree by default** (local planning artifact under `docs/superpowers/`).

**WIP caveat (2026-08-13):** Expanded adduct list, polarity filter, ordered ion types, and neutral-mass multi-charge edges are in code, but **multi-form labeling is not production-ready**. See **§12 Known challenges**. Do not treat complex multi-adduct groups on real panels as validated.

---

## 1. Goal

After consensus mass features are formed and gap-filled, and **before** molecular annotation (formula / MS2 search), group consensus features into families that represent the same underlying analyte under different **natural-abundance isotopic** and **adduct** forms.

### In scope (all stages of this work)

1. Natural-abundance isotope families (e.g. ¹²C / ¹³C, other rare forms above a natural-abundance floor).
2. Adduct families via pairwise mass offsets from `ion_type_dict` (e.g. `[M+H]+` ↔ `[M+NH4]+`).
3. Later: in-source fragments (MS2-informed), then multicore scale-out.

### Explicit non-goals

- Tracer / enriched / labeled isotope experiments (e.g. ¹³C metabolic labeling).
- Changing how consensus features are built (grouping **consumes** consensus clusters).
- Removing `LCMSBase.find_c13_mass_features` (kept optional; see collision policy).
- Exporting correlation scores as product columns.
- Inventing monoisotopic parents when geometry does not support them.

### Hard constraints

1. **Do not invent parents.** Only assign isotope roles when a chemical mono is present in the group via interpretable geometry.
2. Unexplained spacing → no isotope/adduct labels for those pairs.
3. Isotope mass deltas from **`Atoms`** only (never hard-coded `1.003355`).
4. Adduct mass deltas from **`ion_type_dict`** + `Atoms` (add/subtract atom counts).
5. No masscube / MSAC dependencies.
6. No mono-vs-family **height prior** (removed: incorrect for large lipids / multi-charge where M+1 can dominate).

---

## 2. When it runs

| Order | Step |
|-------|------|
| 1 | Per-sample feature finding / prep |
| 2 | `align_lcms_objects()` |
| 3 | `add_consensus_mass_features()` |
| 4 | Gap-fill via `process_consensus_features` (or equivalent) |
| 5 | **`group_consensus_features()`** |
| 6 | Annotation (MS1 formula, MS2 search) — **after** grouping |

**Pipeline hook:** `process_consensus_features(..., group_features=False)`.

- Runtime kwarg only (default **off**).
- When `True`, run after gap-fill and before formula/MS2 if those are also requested.
- Thin wrapper around `group_consensus_features()`.

---

## 3. Architecture

### 3.1 Public API

| Piece | Location |
|-------|----------|
| Entry | `LCMSCollection.group_consensus_features()` |
| Pipeline flag | `group_features` on `process_consensus_features` |
| Algorithm | `corems/mass_spectra/calc/feature_grouping.py` |
| Settings | `LCMSCollectionSettings` (`feature_group_*`) |
| Isotope masses | `Atoms.atomic_masses` |
| Adduct definitions | `ion_type_dict` in `corems.mass_spectra.output.export` |

### 3.2 Operating unit: consensus cluster

Grouping is on **consensus features** (one row per `cluster`), not every sample-level feature id.

| Field | Source |
|-------|--------|
| `cluster` | Cluster id |
| `mz` | `mz_median` (or equivalent) |
| `rt` | `scan_time_aligned_median` |
| Height vector length \(S\) | Per-sample apex **intensity** for that cluster |

Height matrix \(H \in \mathbb{R}^{N \times S}\). Missing heights → **0.0** after gap-fill.

### 3.3 Output columns (no correlation score)

| Column | Type | Meaning |
|--------|------|---------|
| `feature_group_id` | int or None | Shared family id; `None` if unlabeled |
| `ion_role` | str or None | `"mono"` \| `"isotope"` \| `None` — **chemical** role only |
| `ion_type` | str or None | `ion_type_dict` key, e.g. `"[M+H]+"`, `"[M+NH4]+"` |
| `isotope_state` | str or None | e.g. `"M+0"`, `"13C1"`, `"13C2"` |
| `mono_cluster_id` | int or None | Chemical monoisotope cluster for this row's `ion_type` (self if mono; that form's mono if isotope) |

**Invariants**

1. If `ion_role == "isotope"`, `mono_cluster_id` points at a mono in the same `feature_group_id` (same ion form).
2. Never emit isotope labels without a chemical mono in the same group/form.
3. A group is created only when ≥1 non-mono member is accepted (isotope and/or linked adduct form). Singletons stay unlabeled.
4. **No** assumption that mono is the tallest peak in the envelope.
5. Edge tables may still use `parent_idx`/`child_idx` for directed geometry; that is not the label column `mono_cluster_id`.

### 3.4 Charge

```text
unit_spacing = (m_rare - m_mono) / abs(z)     # isotopes
adduct_spacing = (offset(type_b) - offset(type_a)) / abs(z)
```

- `z` from `feature_group_min_charge` … `feature_group_max_charge` (default **1…1**).

### 3.5 Collision with per-file C13

| Path | Role |
|------|------|
| `LCMSBase.find_c13_mass_features` | Optional single-file / legacy |
| `group_consensus_features` | Preferred for **collection** families |

Do not treat both as authoritative in one workflow. Collection labels are a separate layer from per-file `monoisotopic_mf_id` / `isotopologue_type`. `drop_isotopologues` continues to use **per-file** mono labels only.

---

## 4. Algorithm (Approach A′)

Target scale (eventual): ~20k consensus features × ~200 samples. Multicore is **last** (Stage 4).

### 4.1 Geometry candidates

Sparse RT ∩ m/z pair search (KDTree / sparse distance matrices):

- RT window: `alignment_rt_tol`
- m/z window sized to largest relevant unit spacing + ppm pad
- m/z tol: `max(mz_i, mz_j) * alignment_mz_tol_ppm * 1e-6`

### 4.2 Isotope edges (natural abundance)

For each mono element in `feature_group_isotope_atoms` (default `("C",)`):

1. Rare forms from `Atoms` with natural abundance ≥ `feature_group_min_isotope_abundance`.
2. **Unit** spacing only: `(m_rare − m_mono) / |z|` (not direct M+n edges).
3. Higher-order isotopologues (¹³C₂, …) via **roll-up** along successive unit edges (no ¹³C₂ without ¹³C₁). Cap depth with `feature_group_max_isotope_offset` (default 4).
4. Parent = chemical mono **side** of the unit step from Atoms (lighter for ¹³C; can be heavier for ⁵⁴Fe).

### 4.3 Adduct edges (pairwise ion types)

For every unordered pair `(type_a, type_b)` in `feature_group_ion_types`:

1. `offset(type) = mass(atoms_to_add) − mass(atoms_to_subtract)` from `ion_type_dict`.
2. Spacing = `|offset(a) − offset(b)| / |z|`.
3. Lighter peak → lower-offset `ion_type`; heavier → higher-offset `ion_type`.
4. **No designated base form**; intensity does not choose roles.

Default types: `("[M+H]+", "[M+NH4]+")`. Empty or single type → isotopes only.

### 4.4 Quant gate (Pearson on apex intensity)

Fixed policy:

- **Pearson only** (no Spearman/cosine switch).
- **Apex `intensity` only** (not area).
- Pairwise-complete: only samples where **both** heights &gt; 0.
- `n_shared ≥ ceil(feature_group_min_shared_sample_fraction * n_samples)` (default **0.15**).
- `r ≥ feature_group_corr_threshold` (default **0.80**).
- Do **not** store `r` on output tables.

Applied to isotope edges and adduct edges.

### 4.5 Labeling

1. **Isotope roll-up** on correlated unit edges → preliminary families; `ion_role` mono/isotope; `ion_type` set if only one configured type, else deferred.
2. **Merge** across correlated adduct edges: union `feature_group_id`; assign `ion_type` per side; keep chemical mono/isotope roles per form; `mono_cluster_id` = mono of that form.
3. No mono-vs-family height filter.

### 4.6 Idempotence

Re-running `group_consensus_features()` clears prior group columns, then recomputes.

### 4.7 Multicore (Stage 4 — last)

Not required for Stage 1–3 correctness. When needed:

- RT partitions + halo ≥ `alignment_rt_tol`
- `feature_group_partition_size` (default 5000), `cores` from collection settings
- Commit labels only for partition core RT; renumber `feature_group_id` globally

Until Stage 4, implementation may remain single-process.

---

## 5. Settings (`LCMSCollectionSettings`)

| Setting | Default | Purpose |
|---------|---------|---------|
| `alignment_rt_tol` | existing | RT window (reused) |
| `alignment_mz_tol_ppm` | existing | ppm window (reused) |
| `feature_group_min_charge` | `1` | Min \|z\| |
| `feature_group_max_charge` | `1` | Max \|z\| |
| `feature_group_isotope_atoms` | `("C",)` | Mono elements for natural-abundance isotope edges |
| `feature_group_min_isotope_abundance` | `0.01` | Min natural abundance for rare forms |
| `feature_group_max_isotope_offset` | `4` | Max roll-up depth |
| `feature_group_corr_threshold` | `0.80` | Pearson gate |
| `feature_group_min_shared_sample_fraction` | `0.15` | Min shared non-zero sample fraction |
| `feature_group_ion_types` | `("[M+H]+", "[M+NH4]+")` | Pairwise adduct forms (`ion_type_dict` keys) |
| `feature_group_partition_size` | `5000` | Multicore (Stage 4) |

No separate `feature_group_rt_tol` / `feature_group_mz_tol_ppm`.  
Removed: `feature_group_mono_height_fraction`, `feature_group_base_ion_type`, `feature_group_adduct_ion_types`.

---

## 6. Error handling

| Condition | Behavior |
|-----------|----------|
| No consensus / empty cluster summary | `ValueError` |
| \(N < 2\) | No-op; all unlabeled |
| Invalid settings / unknown atom or `ion_type` | `ValueError` at start |
| Parallel worker failure (Stage 4) | Surface exception; no silent partial commit |

No `print()` in importable library paths; use existing verbose/logging patterns.

---

## 7. Testing and local validation

### 7.1 Unit tests (`tests/test_lcms_feature_grouping.py`)

Examples (not exhaustive):

| Case | Expectation |
|------|-------------|
| Mono + ¹³C₁, high corr | One group; roles/states/parent correct |
| + NH₄ + ¹³C–NH₄ | Same `feature_group_id`; each form has mono/isotope; `ion_type` set |
| RT outside tol | No group |
| Pearson below threshold | No group |
| Mono shorter than ¹³C₁ | Still labels (no height prior) |
| ¹³C₁ + ¹³C₂ chain | Roll-up; no direct mono→M+2 edge required |
| No ¹³C₂ without ¹³C₁ | No skip-level label |
| Fe ⁵⁴ / ⁵⁶ | Mono on most-abundant side |
| Charge 2 spacing | Pluggable charge works |
| Re-run | Overwrites cleanly |
| `Atoms` spacing | No hard-coded 1.003355 |

### 7.2 Debugger / eval scripts (local, not CI)

| Path | Role |
|------|------|
| `tmp_data/debug_consensus_feature_grouping.py` | End-to-end collection → grouping; **keep** |
| `tmp_data/feature_grouping_eval/` | Multi-panel (RP/HILIC) preprocess + corr sweeps |

These stay developer tools under `tmp_data/` (not required in the library package).

---

## 8. Export

- Group columns should appear where cluster summary / mass feature tables are exported.
- No standalone report product required for early stages.
- Do not map collection labels onto per-file `monoisotopic_mf_id` in early stages.

### 8.1 Persistence today vs later

| Layer | Today (Stages 1–2) | Future (see Stage 3 polish / §9.1) |
|-------|--------------------|-------------------------------------|
| In-memory | `feature_group_dataframe` + merge into `cluster_summary` / mass-feature tables after `group_consensus_features()` | Unchanged |
| Tabular dump | Debug CSVs under `tmp_data/` (e.g. `feature_group_labels.csv`) | Optional library-side table export |
| Collection HDF5 | **`LCMSCollectionExport` / `ReadSavedLCMSCollection` do not yet save or restore** `feature_group_*` columns | **Required** for reload without re-running grouping |
| Per-sample HDF5 | Cluster assignments / induced features only | Group labels stay collection-level (not remapped to per-file C13 fields) |

Until §9.1 is implemented, reloading a gap-filled collection and inspecting groups means either re-running `group_consensus_features()` (cheap vs gap-fill) or re-attaching labels from a side CSV. That is acceptable for debug, not a product contract.

---

## 9. Review stages (one MR, staged review)

Work ships in **one merge request** (!248). Review and implementation proceed in stages:

| Stage | Scope | Status (as of 2026-08-13) |
|-------|--------|---------------------------|
| **1 — Isotopes** | Natural-abundance unit edges, Pearson, roll-up, settings, unit tests, API hook | **Core complete** (reviewable) |
| **2 — Adducts** | Pairwise `feature_group_ion_types`, merge labeling, tests; polarity filter; ordered ion types; multi-charge | **WIP** — geometry + merge bugs (§12); harden before shipping multi-adduct |
| **3 — ISF / polish** | In-source fragments (MS2-informed); **save/reload feature grouping**; export/docs polish; collision notes | Not started |
| **4 — Multicore** | RT partition + halo; scale to large N | **Last**; not blocking 1–3 |

### Stage 1 acceptance (review checklist)

- [x] `group_consensus_features()` + `group_features=` wired (default off)
- [x] Settings for isotopes + quant gates
- [x] Unit tests green for isotope cases
- [x] No mono-height prior
- [x] Natural-abundance scope documented in code
- [ ] MR description reflects Stage 1 + Stage 2 WIP (this update)
- [ ] Design doc matches code (this update)

### Stage 2 acceptance (target)

- [x] Pairwise ion-type edges from `ion_type_dict` (no base form)
- [x] Mono + ¹³C + NH₄ + ¹³C–NH₄ unit test
- [x] Polarity filter on `feature_group_ion_types` before adduct search (drop opposite-sign keys such as `[M+HCOO]-` / `[M+CH3COO]-` on positive collections; mixed-polarity collection raises)
- [x] Expanded `ion_type_dict` common adduct set (incl. multi-charge `[M+2H]2+`, etc.); adduct edges use per-type `|z|` + neutral-mass match; debug list + unit tests
- [x] Literature-ordered `DEFAULT_ION_TYPES` + rank-sum tie-break on exact residual ties
- [x] Polarity filter; stage timings for debug
- [ ] **Blocker:** merge / multi-charge correctness (§12) — not optional polish
- [ ] Real-panel eval with adducts enabled (reuse `tmp_data/feature_grouping_eval`) after §12 fixes
- [ ] Harden merge edge cases + regression tests from rp_pos g7 / g16 / g24 / g38

### Stage 3–4 acceptance

- Documented when those stages start; Stage 3 includes **§9.1 save/reload**.

### 9.0 Follow-up (out of MR !248) — file issue for feature-group visualization

**Do not implement visualization in this MR.** A prototype extension of `plot_cluster` (MS1 stems/labels for co-group isotopes/adducts on feature group 0 / rp_pos) was tried for feedback and **reverted**: full-spectrum x-range plus dense isotopologue labels were hard to read; product design needs a dedicated issue.

**Follow-up action**

1. **File a GitLab issue** (after or alongside !248) for **visualizing consensus feature groups**.
2. Suggested problem statement: when inspecting a consensus mass feature that belongs to a `feature_group_id`, highlight and label associated natural-abundance isotopes and adducts (and later ISF) on existing collection plot APIs (prefer extending `plot_cluster` / MS1—not a new plot type unless needed).
3. Capture prototype lessons in the issue:
   - Auto-zoom (or dual view) to the group m/z envelope rather than the full MS1 scan range
   - Avoid overlapping rotated labels on isotope ladders; consider stacked annotations, a side legend table, or peak markers with a key
   - Color/style by `ion_type` and mono vs isotope; mark the focus cluster clearly
   - Depend on grouping labels (and ideally §9.1 save/reload for iterative plot work without re-gap-filling)
4. Out of scope for the issue until product owners decide: interactive GUI, export-only SVG reports, replacing per-file C13 plots.

**Acceptance for this design note:** issue filed with link recorded here when created; no library plot change required in Stages 1–2.

### 9.1 Future step — save / reload feature grouping (Stage 3 polish)

**Motivation.** Gap-fill is the expensive step in debug/eval loops. Collection save/load already restores alignment, consensus clusters, and induced features via `LCMSCollectionExport.export_to_hdf5` + `ReadSavedLCMSCollection`. Feature-group labels are **not** part of that round-trip yet, so every reload either re-runs grouping or depends on ad-hoc CSVs.

**Goal.** Persist and restore collection-level grouping so a saved, gap-filled collection can be reopened with the same `feature_group_id` / `ion_role` / `ion_type` / `isotope_state` / `mono_cluster_id` without re-executing the algorithm (unless the user opts to recompute).

**Proposed scope (not implemented yet)**

1. **Write path (`LCMSCollectionExport`)**  
   - If `feature_group_dataframe` is present and non-empty, write a dedicated group in the collection HDF5 (e.g. `feature_group_labels`) keyed by `cluster`, with the five `GROUP_COLUMNS` (nullable dtypes preserved).  
   - Optionally mirror the same columns into any collection-level cluster-summary / mass-feature tables that the exporter already writes.  
   - Do **not** invent a separate product file format; stay on the existing collection HDF5 (+ parameters TOML/JSON).

2. **Read path (`ReadSavedLCMSCollection.get_lcms_collection`)**  
   - After cluster assignments (and induced features) are restored, load `feature_group_labels` when present.  
   - Set `lcms_collection.feature_group_dataframe`.  
   - Re-merge columns into `mass_features_dataframe` / `induced_mass_features_dataframe` / `cluster_summary_dataframe` the same way `group_consensus_features()` does today (shared helper preferred for DRY).  
   - Missing group dataset → leave unlabeled (backward compatible with older collection files).

3. **Recompute vs restore**  
   - Default reload: **restore** labels when present.  
   - Explicit recompute: caller runs `group_consensus_features()` again (already idempotent: clears then rewrites).  
   - Debug script may grow flags such as `use_saved_collection` / `rerun_grouping` later; that is a `tmp_data/` convenience, not a library API requirement for this step.

4. **Tests**  
   - Round-trip: group → export → load → assert label tables equal (dtypes and NA).  
   - Older HDF5 without the group still loads.  
   - Restored labels available for any later visualization work (§9.0) without re-grouping.

5. **Out of scope for this step**  
   - Persisting intermediate edge/correlation tables.  
   - Mapping collection labels onto per-file `monoisotopic_mf_id` / `isotopologue_type`.  
   - Skipping gap-fill in the algorithm itself (gap-fill remains a separate pipeline stage; save/reload only avoids **re-running** it).

**Acceptance when Stage 3 polish lands**

- [ ] Collection HDF5 includes feature-group labels when grouping was run  
- [ ] `ReadSavedLCMSCollection` restores `feature_group_dataframe` and merged summary columns  
- [ ] Unit/integration test for export → import label fidelity  
- [ ] Design §8 / decisions log updated to “implemented”

---

## 10. Design decisions log

| Topic | Decision |
|-------|----------|
| Algorithm family | A′: RT ∩ Δm → Pearson → roll-up/merge |
| Isotope edges | **Unit** only; higher n via roll-up |
| Mono identity | Geometry / Atoms side of unit step — **not** “tallest peak” |
| Height prior | **Removed** (large lipids / multi-charge) |
| Quant metric | Pearson of apex intensity; not exported |
| Missing heights | 0-fill |
| Charge | \|z\|=1 default; pluggable range |
| Adduct model | Pairwise among `feature_group_ion_types`; no base form |
| Adduct polarity | Filter ion types by collection polarity before edge search (trailing `+`/`-`); mixed-polarity collection raises |
| Ion type order | Most → least common (literature); preserved after polarity filter; used to break exact Δm / neutral-mass ties (e.g. prefer `[M+H-H2O]+` over `[M+H+H2O]+` vs `[M+H]+`) |
| `ion_role` | Chemical mono/isotope only; form in `ion_type` |
| Pipeline default | `group_features=False` |
| Per-file C13 | Keep optional; document collision |
| Multicore | **Stage 4 last** |
| Delivery | One MR; staged review |
| Local validation | Debug script + panel eval under `tmp_data/` |
| Feature-group persistence | **Future Stage 3** (§9.1): extend `LCMSCollectionExport` / `ReadSavedLCMSCollection`; not in Stages 1–2 |
| Feature-group visualization | **Out of !248** (§9.0): file a follow-up issue; prototype `plot_cluster` MS1 annotate was tried and **reverted** |

---

## 11. Current code map (implementation)

| Piece | Path |
|-------|------|
| Algorithm | `corems/mass_spectra/calc/feature_grouping.py` |
| Collection API | `LCMSCollection.group_consensus_features` in `lc_calc.py` |
| Settings | `LCMSCollectionSettings` in `processingSetting.py` |
| Tests | `tests/test_lcms_feature_grouping.py` |
| Debug | `tmp_data/debug_consensus_feature_grouping.py` |
| Panel eval | `tmp_data/feature_grouping_eval/` |
| Collection save/load (gap-fill / clusters today; **groups later §9.1**) | `LCMSCollectionExport`, `ReadSavedLCMSCollection` |
| Feature-group visualization | **Not in library**; follow-up issue per §9.0 |
| Known multi-adduct bugs | §12; WIP commit note on branch |

---

## 12. Known challenges (Stage 2 adducts — WIP)

Diagnosed on **rp_pos** debug outputs (`tmp_data/feature_grouping_debug/grouping_outputs_rp_pos_adducts/`) after expanded common adducts + ordered ion types. Grouping itself is fast (~1 s / 800 clusters); **correctness** is the issue.

### 12.1 Invariant that is violated

Within one `feature_group_id`, each `(ion_type, isotope_state)` should correspond to a unique cluster m/z, and all forms should imply the **same neutral mass**

\[
M = |z| \cdot m/z - \mathrm{offset}(\mathrm{ion\_type})
\]

within ppm. Real panels currently produce groups that break this (multiple `[M+H]+` monos at different \(M\), cross-m/z “families”).

### 12.2 Failure modes (with example groups)

| ID | Failure | Example |
|----|---------|---------|
| **M1** | Merge **repaints entire** current group sides with one `ion_type` | Later adduct edge overwrites forms already labeled |
| **M2** | **Transitive union** without shared-\(M\) check | A–B and B–C merge even if \(M_A \neq M_C\) |
| **M3** | **Multi-charge cross-links** too weak | Same peak interpreted as mono-z1 of one \(M\) and high-\(z\) of another; bridges unrelated m/z (**g7**: 391 `[M+H]+` glued to 798 via false 2+) |
| **M4** | **First edge wins**; later edges do not re-type | Rank-sum preference lost when a worse exact edge merges first (**g24**: labels show NH₃ / water-adduct instead of preferred water-loss / NH₄) |
| **M5** | Multiple **monos of same ion_type** at different m/z | **g16**, **g38** — soup of coeluting lipids stamped `[M+H]+` |
| **M6** | Exact **Δm collisions** among ion-type pairs | e.g. \|Δm\| = 18.0106 for water loss vs water adduct vs second −H₂O; \|Δm\| = 17.0265 for NH₄ vs −NH₃; order helps pairwise ties only, not M1–M5 |

### 12.3 Concrete panel examples (rp_pos)

- **g7** (`[M+H]+` / `[M+Na]+` at \(M \approx 390\)): correct core 64/233/26; **false** multi-charge edges pull in 2341/2322 at \(M \approx 797\)–802.  
- **g16**: good 144/179 NH₃-loss / `[M+H]+` pair; extras 30/406 labeled same form with incompatible \(M\).  
- **g24**: geometry admits preferred water-loss/NH₄ and alternate NH₃/water-adduct; **labels follow first merge**, not global preference.  
- **g38**: many `[M+H]+` M+0 rows at different m/z under one `feature_group_id` / shared `mono_cluster_id`.

### 12.4 Multi-charge “disappeared?”

Not removed. After ordered types, **`[M+2H]2+` / `[M+2H-NH3]2+` still appear** (few rows; e.g. multi-form g68) but:

- Rare in intensity-sorted “top multi-form” lists.  
- Used as **bridges** in bad groups (g7) as much as true mono↔2+ pairs.  
- On pure Δm ties with mono-charge pairs, **earlier mono-charge ranks win** by design.

### 12.5 Performance note

On rp_pos (803 clusters, ~18 pos ion types after polarity filter): gap-fill ~200 s; **grouping ~1 s** (adduct_edges ~0.7 s). Scaling concern for Stage 4 later; **not** the current blocker.

### 12.6 Fix direction (scope in follow-on design / plan)

Priority P0 (correctness):

1. Paint **form subtrees only** on adduct merge (never re-type other forms).  
2. **One interpretation per cluster**; reject conflicting type assignment.  
3. **Shared-\(M\) check** on merge (and/or post-pass split).  
4. Restrict multi-charge: whitelist same-series mono↔multi (e.g. `[M+H]+`↔`[M+2H]2+`) or same-\(|z|\) general pairs + separate strict multi pass.  
5. Global edge selection (score all edges; greedy under constraints) so first-edge-wins ends.  
6. Validation invariant: unique `(ion_type, isotope_state)` per group; form monos share \(M\).

Regression fixtures: synthetic + member sets from **g7 / g16 / g24 / g38**.

Detailed approaches and test matrix: brainstorming / follow-on design under `docs/superpowers/specs/` (Stage 2 correctness pass).

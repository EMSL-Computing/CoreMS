# Feature-Group Adduct Correctness Pass — Design

**Date:** 2026-08-13  
**Status:** Draft for review (approved direction; not yet implemented)  
**Branch:** `272_lcms_consensus_feature_grouping`  
**Parent design:** `2026-08-06-lcms-consensus-feature-grouping-design.md` (§12 known challenges)  
**Related:** MR !248 / issue #272, Stage 2 adducts  

---

## 1. Goal

Make multi-form consensus feature grouping **chemically consistent** on real panels so that:

1. Within one `feature_group_id`, every `(ion_type, isotope_state)` maps to a **unique** cluster m/z.
2. All ion forms in the group imply the **same neutral mass**  
   \(M = |z| \cdot m/z - \mathrm{offset}(\mathrm{ion\_type})\) within ppm.
3. Multi-charge forms (e.g. `[M+2H]2+`) still link when real, without gluing unrelated high-m/z clusters into low-m/z families (rp_pos **g7**-style failures).
4. Ion-type preference (literature order) is applied **globally**, not only on pairwise ties (**g24**).

This is a **correctness pass** on the existing Approach A′ pipeline (geometry → Pearson → isotope roll-up → adduct merge). It is not ISF, save/reload, visualization, or multicore.

---

## 2. Background (what is broken)

Diagnosed on rp_pos debug outputs after expanded adducts + ordered ion types. Grouping is fast (~1 s / 800 clusters); labels are not trustworthy for complex multi-form groups.

| ID | Failure | Symptom |
|----|---------|---------|
| **M1** | Merge repaints **entire** current group sides with one `ion_type` | Later edge overwrites other forms |
| **M2** | Transitive `feature_group_id` union without shared-\(M\) check | A–B and B–C merge though \(M_A \neq M_C\) |
| **M3** | Free multi-charge ↔ mono-charge neutral-mass links | Same peak as mono of \(M_1\) and high-\(z\) of \(M_2\); false bridges (**g7**) |
| **M4** | First edge wins; later edges do not re-type | Preferred ion pair loses to earlier worse edge (**g24**) |
| **M5** | Multiple monos of the same `ion_type` at different m/z | **g16**, **g38** |
| **M6** | Exact Δm collisions among type pairs | Order helps pure ties only; does not fix M1–M5 |

**Invariant (target):** one shared \(M\) per group; one chemical mono per ion form; isotopes only of that form’s mono.

---

## 3. Approved approach

**Surgical constraints + scored adduct edges** (not a form-centric full rewrite).

1. Keep isotope stage as-is (unit Atoms edges, roll-up, multi-\(|z|\) for ¹³C of multi-charge forms when charge range allows).
2. **Pass 1 — same-\(|z|\) adducts:** only pairs where `ion_type_charge(a) == ion_type_charge(b)`.
3. **Pass 2 — mono↔multi series:** only pairs listed in an explicit **series map** (related mono/multi of the same chemistry, including multiple multi-charge partners per mono).
4. **Merge:** paint **form isotope subtree only**; shared-\(M\) check; ≤1 `ion_type` per cluster.
5. **Edge order:** sort all accepted edges by residual ↑, then rank-sum ↑, then process once under constraints.

---

## 4. Architecture

```text
mz, rt, heights (apex)
        │
        ▼
  find_isotope_edges  →  Pearson  →  assign_isotope_labels
        │                      (forms: mono + isotopes, mono_cluster_id)
        ▼
  Pass 1: same-|z| adduct edges
  (ion_type_dict offsets / same-z Δm or neutral-mass with za=zb)
        │
        ▼
  Pass 2: mono↔multi edges only if (type_a, type_b) ∈ SERIES_PAIRS
  (neutral mass M = z·mz − offset)
        │
        ▼
  Sort edges: residual, rank_sum(type_a)+rank_sum(type_b)
        │
        ▼
  merge under constraints (form paint, shared M, one type/cluster)
        │
        ▼
  validate_group_invariants (optional split / unlabel violators)
```

### 4.1 Components

| Piece | Responsibility |
|-------|----------------|
| `ion_type_charge` / polarity filter | Unchanged |
| `DEFAULT_ION_TYPES` order | Rank for tie-break / edge scoring |
| `find_adduct_edges` (or split helpers) | Pass 1 same-\(|z|\); Pass 2 series-only multi |
| `SERIES_PAIRS` | Frozen set/map of allowed mono↔multi type pairs |
| `merge_adduct_edges_into_labels` | Form-local paint; shared \(M\); no full-group retype |
| Edge sort | Residual then rank-sum before merge |
| Validation | Group-level invariant check |

### 4.2 Series map (Pass 2)

Declarative pairs of `ion_type_dict` keys. Direction does not matter. Minimum initial set (extend as needed):

| Mono / lower form | Multi / related |
|-------------------|-----------------|
| `[M+H]+` | `[M+2H]2+`, `[M+3H]3+` |
| `[M+Na]+` | `[M+2Na]2+`, `[M+H+Na]2+` |
| `[M+K]+` | `[M+H+K]2+` |
| `[M+NH4]+` | (none required in v1 unless a 2+ NH₄ form is added to `ion_type_dict`) |
| `[M-H]-` | `[M-2H]2-` |

**Policy:** `[M+Na]+` may link to **both** `[M+2Na]2+` and `[M+H+Na]2+` when neutral mass matches. Do **not** allow arbitrary cross-type multi-charge (e.g. light peak as `[M+2H-NH3]2+` of a heavy `[M+H]+` stranger).

Implementation: undirected set of frozensets `{type_a, type_b}`, or ordered pairs with normalization. Unknown multi-charge combinations → **no edge** in Pass 2.

### 4.3 Same-\(|z|\) Pass 1

For every unordered pair of configured ion types with equal `ion_type_charge`:

- Use existing neutral-mass residual  
  \(|z \cdot mz_i - off_i - (z \cdot mz_j - off_j)| \le\) ppm of \(M\)  
  (equivalent to constant Δm when \(z\) matches).
- RT coelution gate unchanged (`alignment_rt_tol` / params.rt_tol).
- Pearson height gate unchanged.

No multi-charge in Pass 1.

### 4.4 Merge rules (replaces current paint-all-members)

When processing edge \((p, c)\) with types \(t_p, t_c\):

1. **Form of \(p\):** all clusters currently sharing `mono_cluster_id` with \(p\)’s chemical mono **and** the same prior `ion_type` (or unlabeled isotope family still tied only by isotope stage to that mono).  
   Practical rule v1:  
   - If \(p\) has `ion_role` mono/isotope and `mono_cluster_id` set: form members = all clusters with that `mono_cluster_id` **whose `ion_type` is null or equal to \(t_p\)**.  
   - If \(p\) unlabeled singleton: form = \(\{p\}\).

2. **Assign** \(t_p\) only to form(\(p\)); \(t_c\) only to form(\(c\)). Never write \(t_p\) onto members that already have a **different** non-null `ion_type`.

3. **Shared \(M\):** compute \(M\) from each form’s mono cluster (or endpoint if mono) using assigned types. Reject edge if \(|M_p - M_c| >\) tol.

4. **One type per cluster:** if endpoint already has a different `ion_type`, reject edge (or keep existing if better score—v1: **reject** to avoid thrashing).

5. **Union** `feature_group_id` only if edge accepted. Update `mono_cluster_id` within each form only (mono of that form).

6. **Isotope roles** stay mono/isotope from isotope stage; do not demote form mono of the heavier adduct to “isotope.”

### 4.5 Edge ordering

Before merge:

1. Collect Pass 1 + Pass 2 edges that pass geometry + Pearson.  
2. Sort by:  
   - `residual` ascending  
   - `rank(t_p) + rank(t_c)` ascending (`rank` = index in `params.ion_types`)  
   - stable secondary key (e.g. parent/child cluster id)  
3. Process in that order under §4.4.

This replaces pure discovery-order first-edge-wins for multi-edge coelution webs.

### 4.6 Validation pass (after merge)

For each `feature_group_id`:

- For each `ion_type` present: exactly one `ion_role == "mono"` (or zero if only isotopes—should not happen).  
- All isotopes of that form share that form’s `mono_cluster_id`.  
- Implied \(M\) from each form mono agrees within ppm.  
- Unique `(ion_type, isotope_state)` among labeled members.

**On violation (v1):** remove the worst member(s) from the group (clear group columns for those clusters) rather than inventing parents. Log counts in timings if useful. Prefer deterministic order (e.g. highest residual from form mono).

---

## 5. Data flow / API

| Surface | Change |
|---------|--------|
| `group_features_arrays` | Orchestrate pass1 → pass2 → sort → merge → validate; keep `timings_out` keys for each stage |
| `FeatureGroupParams` | No new user knobs required for v1 (series map is code constant; optional later) |
| `feature_group_ion_types` | Still ordered most→least common; polarity filter unchanged |
| Public columns | Unchanged: `feature_group_id`, `ion_role`, `ion_type`, `isotope_state`, `mono_cluster_id` |
| `group_consensus_features` | Unchanged call signature |

No change to gap-fill pipeline order.

---

## 6. Error handling

| Condition | Behavior |
|-----------|----------|
| Empty series map / no multi types configured | Pass 2 no-op |
| Conflicting type on cluster | Reject edge |
| Shared \(M\) fail | Reject edge |
| Validation failure | Unlabel violator clusters; keep consistent core |
| Unknown ion type in series map | `ValueError` at validate_params / import time |

No `print` in library paths; timings via `timings_out` only.

---

## 7. Testing

### 7.1 Unit tests (required)

| Case | Expectation |
|------|-------------|
| Same-\(z\) NH₄ / H | Still one group; roles mono/isotope per form |
| Water loss vs water adduct tie | Prefer `[M+H-H2O]+` / `[M+H]+` via order + sort |
| **g7-like false multi-charge** | 391 H/Na family **not** unioned with 798 via 2+ bridge |
| True mono↔`[M+2H]2+` | Same `feature_group_id`; types correct; shared \(M\) |
| Na series | `[M+Na]+` can link to `[M+2Na]2+` **and** separately to `[M+H+Na]2+` when \(M\) matches (two scenarios) |
| Form paint | Second edge does not re-type first form’s ¹³C members to the new form’s type |
| **g24-like dual exact Δm** | Preferred pair wins when both residual 0 (sorted edges) |
| **g16/g38-like multi mono** | No two `[M+H]+` M+0 at incompatible \(M\) in one group after validation |
| Polarity filter | Unchanged |
| Pass 1 rejects different \(z\) | No edge for arbitrary z-mismatch type pair outside series map |

### 7.2 Regression data (optional but preferred)

Member m/z/RT tables extracted from rp_pos **g7, g16, g24, g38** (fixed arrays in test file; no raw files required). Assert post-fix membership and types.

### 7.3 Manual debug

Re-run `tmp_data/debug_consensus_feature_grouping.py` on rp_pos; inspect groups 7, 16, 24, 38 and multi-charge counts; confirm timing still ~1 s class for grouping.

---

## 8. Implementation sketch (for plan skill later)

1. Add `SERIES_PAIRS` + `is_allowed_multi_charge_pair(a, b)`.  
2. Refactor `find_adduct_edges` into same-\(z\) finder + multi series finder (or filter after).  
3. Sort combined edges; rewrite merge (form members helper, shared \(M\), reject conflicts).  
4. Add `validate_feature_group_labels` post-pass.  
5. Tests from §7.  
6. Re-run rp_pos debug; update parent design §12 status when fixed.

**Estimated risk:** medium (merge rewrite); isotope path untouched if roll-up output is treated as form seeds only.

---

## 9. Out of scope

- In-source fragment linking (Stage 3).  
- Collection HDF5 save/reload of labels (§9.1 parent design).  
- Feature-group visualization (§9.0).  
- Multicore partition (Stage 4).  
- Changing Pearson / apex-only quant gate.  
- Expanding chemistry of `ion_type_dict` beyond what correctness needs.  
- Intensity-based adduct priors beyond existing ion-type order.

---

## 10. Success criteria

- [ ] g7-like fixture: only chemically consistent H/Na family; no 798 bridge.  
- [ ] True `[M+H]+` ↔ `[M+2H]2+` still groups.  
- [ ] Na mono can group with either allowed multi-charge partner.  
- [ ] No group with two `[M+H]+` monos at different \(M\).  
- [ ] g24-like preference stable under dual exact geometry.  
- [ ] Existing isotope + NH₄ unit tests still green.  
- [ ] Parent design §12 updated to “addressed” with residual caveats if any.  
- [ ] Grouping wall time remains acceptable on rp_pos (order of ~seconds, not minutes).

---

## 11. Decisions log

| Topic | Decision |
|-------|----------|
| Approach | Surgical + scored edges (not full form-centric rewrite) |
| Multi-charge | Pass 1 same \|z\| only; Pass 2 series map only |
| Na multi | Allow both `[M+2Na]2+` and `[M+H+Na]2+` from `[M+Na]+` |
| Merge paint | Form subtree only |
| Shared \(M\) | Required on merge |
| Edge order | Residual then rank-sum |
| Conflict | Reject edge (v1) |
| Validation | Unlabel violators |

---

## 12. Spec self-review

- No TBD placeholders left for v1 scope.  
- Consistent with parent §12 and user multi-charge policy.  
- Single implementation plan size (one correctness pass).  
- Explicit: Pass 1 vs Pass 2, series map, merge rules, tests.  

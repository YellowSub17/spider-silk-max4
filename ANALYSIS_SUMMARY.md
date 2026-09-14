# Spider silk SAXS/WAXS analysis — summary for handoff

**Revision history**: v1 used `compare_groups.GROUPS` hand-curated from free-text log comments (mislabeled several conditions, missed the capillary screen). v2 rebuilt `GROUPS` from `data/scans_clean.csv`'s authoritative `Sample` column. v3 added a within-run self-consistency check that dropped two internally-inconsistent conditions, a separate capillary pipeline, and a drift-vs-effect-size check. **This version (v4)** discovers that the "position N" label used throughout v1-v3 is not a reliable physical-position identifier, replaces it with real stage coordinates, finds that the mid-run intensity cycle is genuinely position-locked (not just periodic), adds position-matched (rather than pooled) protein-vs-blank subtraction, and **retracts the v1-v3 "no mixing-time effect" conclusion** — see Overall state below. This is the current state of the analysis.

## Data & codebase

- `data/scans_clean.csv` — **authoritative** per-ScanID metadata: `scanID`, `Sample` (format `[CHANNEL1]/[CHANNEL2]`, e.g. `Y2F+kpi/aa` = Y2F protein in KPi buffer, mixed in-channel with acetic acid; `kpi/aa` is the matched protein-free blank), `description` ("position N" label — **not physically reliable across samples, see stage_positions.py below**), `sams4_x`/`sams4_y` (the actual stage coordinates -- now the authoritative position source), `start time`.
- `data/datasetsTable2026-6-1.csv` — superseded, not used for grouping.
- `data/raw/scan-<id>*.h5`, `data/process/azint/scan-<id>_*_integrated.h5` — raw and azimuthally-integrated detector data (Eiger=SAXS, Pilatus=WAXS).
- `src/` (pre-existing, extended over the course of this analysis):
  - `scan.py` / `combo_scan.py` — `Scan`/`ComboScan`: load and merge SAXS+WAXS into `self.qs`/`self.Is`; `Scan` computes a log-log linear-fit `r_squared` streak-quality metric over `det.LINEAR_SAXS_RANGE` (q=0.005–0.015 Å⁻¹), used everywhere as the frame-quality filter (`RSQ_LIM=0.999`).
  - `scale_fns.py` — `norm_i0()` (normalize each frame by incident flux `i0`), alongside pre-existing `norm_qrange()`/`norm_max()`.
  - `saxs_metrics.py` — `guinier_fit` (Rg/I0/r; `q_min`/`q_max`/`refine`/`refine_min_points` params — had a real instability bug, see Key Bugs), `porod_exponent`, `peak_analysis` (WAXS peak height/area/position via baseline subtraction, `q_base` configurable — also had a real bug, see Key Bugs), `kratky`, `dimensionless_kratky`.

## `compare_groups.py` — the foundation everything else imports

Builds `GROUPS` from `scans_clean.csv`'s `Sample` column (`MIN_SCANS_PER_GROUP=3`). Provides:
- `SAMPLE_INFO[name]` — parsed dict: `protein`, `channel1`, `channel2`, `is_capillary`, `has_acid`.
- `BLANK_OF[name]` — matching protein-free blank group name, or `None`.
- `RAW_NAME[name]` — original Sample string.
- `NORM_METHODS` (`linear_saxs`/`neutral`/`i0`), `load_group`/`load_scan_ids`, `kept_frames`/`mean_curve`/**`kept_mask`** (new: exposes the streak-filter boolean mask directly, so callers can line up kept frames with their originating scan IDs — needed by the position-based scripts below), `coeff_of_variation`.
- `drift_check(name, norm, q_base)` — splits a group's own scans into chronological early/late halves and runs the scale-corrected WAXS-peak test *between them* (no condition change, just time). Used by `protein_vs_buffer.py` to sanity-check results against the blank's own internal drift.
- **`SKIP_FIRST_CYCLE = True`** (new, default) — every condition's first full "position cycle" (N chronologically-earliest scans, N = that condition's number of distinct physical stage spots) is now dropped by default before any analysis, via `first_cycle_ids()`/`load_group(..., skip_first_cycle=...)`. Motivated by `drift_within_runs.py` finding the first cycle of every run is consistently a settling transient (e.g. YR2A_h2o_h2o: total intensity 12.2, 12.9, 18.6, 18.0, then ~14 for the rest of the run). `report_first_cycle_skips()` prints how many scans this drops per condition (typically N=5, i.e. 5 scans; up to 20 for the messiest condition, `Y2F_kpi_aa`).
- `DROPPED_CONDITIONS` — conditions removed from `GROUPS` entirely because `drift_within_runs.py` found them internally inconsistent: `Y2F_kpi_kpi` (108.8% early/late RMSD, a real mid-run event) and `resin` (46.6% RMSD, a step-change mid-run).

**Current `GROUPS` (12 conditions):**

| Group name | Raw Sample | n scans | Protein | Acid? |
|---|---|---|---|---|
| `Y2F_kpi_aa` | Y2F+kpi/aa | 251 | Y2F | yes |
| `Y2F_h2o_kpi` | Y2F+h2o/kpi | 165 | Y2F | no |
| `kpi_h2o` | kpi/h2o | 145 | — | no |
| `Y2F_h2o_h2o` | Y2F+h2o/h2o | 120 | Y2F | no |
| `kpi_aa` | kpi/aa | 82 | — | yes |
| `YR2A_h2o_kpi` | YR2A+h2o/kpi | 45 | YR2A | no |
| `Y2F_kpi_h2o` | Y2F+kpi/h2o | 42 | Y2F | no |
| `YR2A_h2o_h2o` | YR2A+h2o/h2o | 38 | YR2A | no |
| `h2o_h2o` | h2o/h2o | 30 | — | no |
| `cap_YF_kpi_shear` | cap:YF/kpi/shear | 7 | Y2F (capillary) | no |
| `cap_empty` | cap:empty | 4 | — (capillary) | no |
| `cap_YF_kpi_aa` | cap:YF/kpi/aa | 4 | Y2F (capillary) | yes |

`BLANK_OF`: `Y2F_kpi_aa`→`kpi_aa`, `Y2F_h2o_h2o`→`h2o_h2o`, `Y2F_h2o_kpi`→`kpi_h2o`, `YR2A_h2o_kpi`→`kpi_h2o`, `YR2A_h2o_h2o`→`h2o_h2o`, `Y2F_kpi_h2o`→`kpi_h2o`.

A richer capillary screen (third protein "WT" = wild-type/natural spider silk, shear/no-shear/gel labels) exists below `MIN_SCANS_PER_GROUP` and is handled separately by `capillary_analysis.py`.

## `stage_positions.py` — authoritative physical position (new, foundational for this revision)

**Key discovery**: the `description` "position N" label is assigned independently *within each sample's own run* — it is NOT a shared physical coordinate system. Concretely verified in this data: `kpi/aa` has only 3 labelled positions at y = 26.999/27.999/32.500mm, while the 5-position flow groups have y = 26.975/27.476/27.976/29.977/32.477mm — `kpi/aa`'s label 2 is physically the same spot as the 5-position groups' label 3, and its label 3 is their label 5. `kpi/aa` also sits at a different x (~117.04mm) than every other flow group (~116.76mm). There's also a per-cycle anomaly: "position 1" in every group sits at a measurably different x (~19.6 µm lateral offset) than positions 2-5 and carries much larger y jitter — one odd stop per cycle, likely a settling/backlash artifact of the stage move sequence.

This module clusters every scan's `(sams4_x, sams4_y)` into global "spot" IDs (0.05mm tolerance) shared across all samples, independent of label:
- `SPOT_OF[scanID]` → global spot id; `SPOT_COORD[spot_id]` → `(x, y)`; `spots_of_group(name)` → `{spot_id: [scanIDs]}`.
- `UNRELIABLE_SAMPLES` — (sample, label) pairs whose scans don't cleanly cluster (label spans multiple physical spots, or >10x-tolerance coordinate spread). Confirmed `Y2F+kpi/aa`'s labels 1/2/3 all span 5-6 distinct physical spots each (label 2 alone spans 4.94mm of y) — this condition's position labels are essentially meaningless.
- `print_diagnostic_table()` (run directly) prints sample/label/spot/coordinate/flag for every condition used downstream.

## `position_locking.py` — is the intensity cycle really position-locked? (new)

For every condition: groups i0-normalized per-frame total intensity and WAXS peak height by global spot id, runs one-way ANOVA + Kruskal-Wallis across spots, and reports η² (variance explained by spot), with and without the first cycle. **Figures**: `figures/position/<condition>_spot_effect.png`, `figures/position/spot_effect_summary.png`.

**Finding: yes, decisively.** 9/9 testable conditions are Kruskal-Wallis-significant (p<0.05) for total intensity, both with and without the first cycle — and η² is often large: `h2o_h2o` 82.3%, `YR2A_h2o_h2o` 81.0% (excl. 1st cycle), `Y2F_kpi_aa` 52.5%, `Y2F_kpi_h2o` 47.5%. Excluding the settling transient generally *increases* η² rather than explaining it away — the position-lock is a real, separate effect, not just the transient. This directly resolves the "cyclic, not monotonic" structure visible in `drift_within_runs.py`'s scan-order plots: it was a spatial cycle (period = number of distinct spots) aliased into an apparent temporal one.

## `position_series.py` — channel-position (mixing-time) effect, rewritten

Previously regressed each metric against the "position N" label's integer index (1..5) — wrong on two counts, both fixed: (1) index treats the 2.0mm position-3→4 step identically to the 0.5mm position-1→2 step, so a distance-linear effect looks nonlinear and underpowered against index; (2) the label itself isn't reliable (see `stage_positions.py`). Now groups frames by physical spot, regresses against real y-offset in mm (`stage_positions.SPOT_COORD`), keeps the old index-based fit printed alongside for comparison, adds a Kruskal-Wallis test across spots (catches non-monotonic dependence a linear fit would miss), and adds total integrated intensity as a 7th metric. **Figures**: `<group>_position_curves.png`, `<group>_position_metrics.png` (now a 2×4 grid, x-axis in mm) for all 9 flow groups.

**Finding: the effect is real but non-monotonic** — both the index-based and the new distance-based linear fits still find almost nothing (0-1 significant out of 9 conditions per metric), but Kruskal-Wallis finds a significant spot-dependence in the *majority* of condition/metric pairs: Porod exponent 8/9, peak height 7/9, peak area 6/9, total intensity 6/9, Rg 5/9, Guinier I0 5/9, peak position 3/9. Visually (e.g. `figures/Y2F_h2o_kpi_position_metrics.png`) the pattern is a dip-then-rise or similar non-monotonic shape, not a trend — exactly why linear regression against either index or distance missed it. **This directly overturns the v1-v3 conclusion of "no mixing-time effect" (which was based solely on the index-based linear fit, 2/54 significant).**

## `protein_vs_buffer.py` — protein vs matched blank, now with position-matched subtraction

`PAIRS` built from `BLANK_OF` (6 pairs). Diagnostics: CV ratio, scaling-factor fit, peak significance (raw + scale-corrected, `DIFF_Q_BASE`-anchored), replicate envelope, relative difference, drift check (all as in v3) — **plus new diagnostic #8, position-matched subtraction** (`position_matched_diff`): instead of pooling all frames, averages protein and blank *separately within each shared physical spot*, combines per-spot differences with equal weight (removing any position-mix bias by construction), and reports each shared spot's own peak significance individually. Explicitly skips (never force-matches to a nearest spot) when coverage is too poor — verified this correctly triggers for `Y2F_kpi_aa` vs `kpi_aa`, which share **zero** physical spots (their x coordinates differ by ~0.28mm, entirely non-overlapping clusters).

**Figures**: `<pair>__<norm>.png` (pooled, as before) plus new `<pair>__<norm>__matched.png` (position-matched).

**Finding — the key result**: `Y2F_kpi_h2o` vs `kpi_h2o` survives position-matched subtraction and gets *stronger* under the best-behaved normalization: pooled scale-corrected 8.1σ → matched 10.5σ (i0 norm), with full 5/5-spot coverage (100% of both groups' frames used — the position-mix skew does not explain this result away). The per-spot breakdown adds an important nuance: the signal is concentrated at the two furthest-downstream spots (y=3.0mm: 4.7σ, y=5.5mm: 6.3σ) and essentially absent at the first three (y=0/0.5/1.0mm: no peak, 0.6σ, 0.4σ) — a pattern that *increases with distance along the channel*, consistent with a real mixing-time-dependent effect rather than a single-spot geometry artifact (which would show one outlier spot, not a monotonic-with-distance pattern). This remains the single most credible positive finding in the analysis, now on firmer footing, though still one pair pending replication. All other pairs' matched results are consistent with their pooled results (no significant peak either way), except `Y2F_kpi_aa` vs `kpi_aa`, correctly skipped for zero coverage.

## `kratky_plot.py` — unchanged this revision

Same 6 `PAIRS`. Guinier Rg (q<0.02) ~505–528 Å for every curve regardless of pair; q<0.1 gives a common ~36–39 Å instead (r²≈0.65–0.70). Neither window separates protein from blank.

## `map_llps_by_rg.py` — systematic Rg-artifact map, unchanged this revision

All 12 `GROUPS` conditions show Rg≈505–528 Å, r²≈0.98–0.99, including `cap_empty`. Not acid- or protein-specific — a shared instrumental/geometric artifact (now additionally explained: the fit always lands on the same ~10 points nearest the detector's lowest measurable q, see Key Bugs).

## `capillary_analysis.py` — separate pipeline for the static capillary screen, unchanged this revision

All 14 capillary conditions, including third protein construct WT (wild-type/natural silk). The ~500 Å Rg artifact also appears here, including in a literally empty capillary. Y2F vs YR2A vs WT indistinguishable under matched sheared-KPi conditions; shear vs no-shear shows no visible effect for Y2F; Y2F in water shows a visibly higher WAXS bump than in KPi. Qualitative only, n=2 per condition.

## `drift_within_runs.py` — within-run self-consistency check, unchanged this revision (now also benefits from SKIP_FIRST_CYCLE)

Computes early-vs-late relative RMSD per condition. **2/14 original conditions were chaotic** (`Y2F_kpi_kpi` 108.8%, `resin` 46.6%) — both dropped (see `DROPPED_CONDITIONS`). Of the remaining 12, now with `SKIP_FIRST_CYCLE` applied by default, **0/9 testable conditions are chaotic**, and RMSDs improved further versus v3 (e.g. `Y2F_kpi_h2o` 5.6%→1.9%, `YR2A_h2o_h2o` 3.4%→2.3%) — consistent with the first cycle being a real, now-removed settling transient rather than the source of the position-locking found above.

## Key bugs found and fixed (all sessions, kept together for reference)

1. **Baseline-anchoring bug in `peak_analysis`**: default baseline windows (q=0.03–0.08, 1.9–2.1) are fine for a raw curve but land inside the low-q Guinier-mismatch region on a *difference* curve, inflating significance by orders of magnitude (a spurious 56σ collapsed to <1σ once fixed via `DIFF_Q_BASE=((0.5,0.7),(1.5,1.7))`). A second copy of the same bug in the replicate-level peak calculations was found and fixed separately.
2. **Guinier refinement instability** (fixed via `refine_min_points` guard + exposing `q_fit_max`/`n_points`): the two-stage fit could collapse onto as few as 5 points hugging the detector's lowest measurable q. In this dataset it always lands on the same ~10 points regardless of sample — why every condition reports the same ~500 Å Rg.
3. **(v4) The "position N" label is not a physical position** — see `stage_positions.py` above. This invalidated the v1-v3 `position_series.py` analysis (which used the label's index as if it were both physically real and evenly spaced) and is the reason the "no mixing-time effect" conclusion has been retracted.

## Overall state / open questions (v4)

1. **Position-locking is real and should inform every future analysis on this dataset.** The mixing-channel intensity (and several derived metrics) depends significantly on physical stage spot in most conditions, non-monotonically, with "position 1" behaving anomalously in essentially every run (see `stage_positions.py`, `position_locking.py`). Any future pooled comparison across conditions should check spot-mix balance first (`stage_positions.spots_of_group`) or use the position-matched approach in `protein_vs_buffer.py`.
2. **"No mixing-time effect" is retracted.** Kruskal-Wallis across physical spots is significant for most metrics in most conditions (`position_series.py`); the earlier null result was an artifact of fitting a non-monotonic effect with a linear model against the wrong (label-index) x-axis. What remains true: no *simple linear* mixing-time trend was found, with either the wrong (index) or correct (distance) x-axis.
3. **`Y2F_kpi_h2o` vs `kpi_h2o` is the strongest surviving lead**, now confirmed by position-matched subtraction (10.5σ, i0 norm, full spot coverage) with a physically sensible per-spot pattern (signal grows with distance downstream, absent at the earliest spots) that argues against it being a geometry artifact. Still just one pair — worth targeted replication (more `Y2F_kpi_h2o`/`kpi_h2o` repeats, ideally position-balanced and temporally interleaved) before treating it as confirmed protein signal.
4. **No LLPS (acid-specific large-Rg) signature** — the ~500 Å feature is universal, confirmed in flow-channel, static-capillary, and zero-sample (`cap_empty`, `resin`-before-drop) geometries. Unaffected by this revision's findings.
5. **`Y2F_kpi_aa` vs `kpi_aa`** (the acid-exposed pair) **has zero shared physical spots** and cannot be position-matched at all — its pooled result (no significant peak, high drift-check failure) should be treated as unreliable/uninformative rather than a negative result; a real test of this condition needs new data collected at matching stage coordinates.
6. **Capillary screen** (`capillary_analysis.py`) still has cleaner data and a third protein (WT) but only 2 repeats per condition — no statistics attempted.
7. Two conditions (`Y2F_kpi_kpi`, `resin`) remain excluded from all analysis as internally unreliable (`compare_groups.DROPPED_CONDITIONS`).

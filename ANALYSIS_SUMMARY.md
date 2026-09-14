# Spider silk SAXS/WAXS analysis — summary for handoff

**Revision history**: v1 used `compare_groups.GROUPS` hand-curated from free-text log comments (mislabeled several conditions, missed the capillary screen). v2 rebuilt `GROUPS` from `data/scans_clean.csv`'s authoritative `Sample` column. **This version (v3)** adds a within-run self-consistency check that found two conditions internally too inconsistent to trust, drops them, adds a separate capillary-screen pipeline, and adds an automated drift-vs-effect-size check to the protein/blank comparison. This is the current state of the analysis.

## Data & codebase

- `data/scans_clean.csv` — **authoritative** per-ScanID metadata: `scanID`, `Sample` (format `[CHANNEL1]/[CHANNEL2]`, e.g. `Y2F+kpi/aa` = Y2F protein in KPi buffer, mixed in-channel with acetic acid; `kpi/aa` is the matched protein-free blank), `description` ("position N" — position along the flow-channel / mixing time), `start time`.
- `data/datasetsTable2026-6-1.csv` — superseded, not used for grouping; its free-text `Comment` is still a readable human log if useful for context.
- `data/raw/scan-<id>*.h5`, `data/process/azint/scan-<id>_*_integrated.h5` — raw and azimuthally-integrated detector data (Eiger=SAXS, Pilatus=WAXS).
- `src/` (pre-existing, extended this session):
  - `scan.py` / `combo_scan.py` — `Scan`/`ComboScan`: load and merge SAXS+WAXS into `self.qs`/`self.Is`; `Scan` computes a log-log linear-fit `r_squared` streak-quality metric over `det.LINEAR_SAXS_RANGE` (q=0.005–0.015 Å⁻¹), used everywhere as the frame-quality filter (`RSQ_LIM=0.999`).
  - `scale_fns.py` — added `norm_i0()` (normalize each frame by incident flux `i0`), alongside pre-existing `norm_qrange()`/`norm_max()`.
  - `saxs_metrics.py` — added: `guinier_fit` (Rg/I0/r; `q_min`/`q_max`/`refine`/`refine_min_points` params — see Key Findings, this had a real instability bug), `porod_exponent`, `peak_analysis` (WAXS peak height/area/position via baseline subtraction, `q_base` configurable — also had a real bug, see Key Findings), `kratky`, `dimensionless_kratky`.

## `compare_groups.py` — the foundation everything else imports

Builds `GROUPS` from `scans_clean.csv`'s `Sample` column (`MIN_SCANS_PER_GROUP=3`). Provides:
- `SAMPLE_INFO[name]` — parsed dict: `protein`, `channel1`, `channel2`, `is_capillary`, `has_acid`.
- `BLANK_OF[name]` — matching protein-free blank group name (matches `{channel1}/{channel2}` or its reverse), or `None`.
- `RAW_NAME[name]` — original Sample string.
- `NORM_METHODS` (`linear_saxs`/`neutral`/`i0`), `load_group`/`load_scan_ids`, `kept_frames`/`mean_curve`, `coeff_of_variation` — the shared loading/normalization/streak-filter backbone every other script uses.
- **`drift_check(name, norm, q_base)`** (new) — splits a group's own scans into chronological early/late halves and runs the scale-corrected WAXS-peak test *between them*, i.e. with no condition change at all, just time. Used by `protein_vs_buffer.py` to sanity-check every result against the blank's own internal drift.
- **`DROPPED_CONDITIONS`** (new) — conditions removed from `GROUPS` entirely because `drift_within_runs.py` found them internally inconsistent:
  - `Y2F_kpi_kpi`: 108.8% early/late RMSD, total intensity swung +387% within 55 minutes — a real mid-run event (clogging/aggregation/bubble), not steady drift.
  - `resin`: 46.6% RMSD over just 5 min/10 scans — a clean step-change between the first 3 scans and the rest, most likely a beam/position shift mid-run.

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

`COMPARISONS`: buffer effect on Y2F, Y2F vs YR2A in water/KPi, buffer-only controls, and **"Y2F in KPi with vs without acid"** (`Y2F_kpi_aa` vs `kpi_aa` only — the closer no-acid-same-buffer baseline `Y2F_kpi_kpi` was dropped, see above).

A richer capillary screen exists in `scans_clean.csv` with a **third protein construct "WT" (wild-type/natural spider silk)** and shear/no-shear/gel-state labels, mostly below `MIN_SCANS_PER_GROUP` — handled separately by `capillary_analysis.py` (below), not folded into `GROUPS`.

## Scripts and their figures

### `compare_groups.py`
Figures: `buffer_effect_on_y2f.png`, `y2f_vs_yr2a_in_water.png`, `y2f_vs_yr2a_in_kpi.png`, `buffer-only_controls.png`, `y2f_in_kpi_with_vs_without_acid.png`.
**Finding**: run-to-run CV is large across the board (12–42%) — the first-pass signal that any condition difference needs checking against noise.

### `position_series.py` — channel-position (mixing-time) effect
`CHANNEL_GROUPS` = every non-capillary group in `GROUPS` (9, built generically via `SAMPLE_INFO`). Position labels from `scans_clean.csv`'s `description` column.
**Figures**: `<group>_position_curves.png`, `<group>_position_metrics.png` for all 9 groups.
**Finding**: **no robust channel-position trend** — 2/54 metric×group tests hit p<0.05 (expected false-positive rate).

### `protein_vs_buffer.py` — protein vs matched blank, with reliability diagnostics
`PAIRS` built programmatically from `BLANK_OF` (6 pairs): `Y2F_kpi_aa`, `Y2F_h2o_h2o`, `Y2F_h2o_kpi`, `YR2A_h2o_kpi`, `YR2A_h2o_h2o`, `Y2F_kpi_h2o`, each vs its matched blank. Diagnostics: CV ratio, scaling-factor fit, peak significance (raw + scale-corrected, baseline anchored near the peak via `DIFF_Q_BASE` — imported from `compare_groups.DIFF_Q_BASE`), replicate envelope, relative difference, **and now a drift check** (#7, below) — each run under `linear_saxs`/`neutral`/`i0` normalization.
**Figures**: `<pair>__<norm>.png`, 18 total (6 pairs × 3 norms).
**Finding**: `Y2F_kpi_h2o` vs `kpi_h2o` is still the best-looking lead (up to 9.5σ scale-corrected under `i0`), but **the drift check now flags every pair anchored on `kpi_h2o` or `kpi_aa` as blank** — those two blanks' own early-vs-late split reproduces 59–441% of the claimed protein effect at the WAXS-peak level (see below), so none of those 4 pairs' results can be trusted as protein signal without correcting for the blank's internal drift first. Only the `h2o_h2o`-anchored pairs (`Y2F_h2o_h2o`, `YR2A_h2o_h2o`) pass the drift check cleanly, but neither shows a significant peak to begin with.

### `kratky_plot.py` — Kratky-shape diagnostic
Same 6 `PAIRS`. **Figures**: `<pair>__kratky.png`.
**Finding**: unchanged — Guinier Rg (q<0.02) ~505–528 Å for every curve regardless of pair; q<0.1 gives a common ~36–39 Å instead (r²≈0.65–0.70). Neither window separates protein from blank.

### `map_llps_by_rg.py` — systematic Rg-artifact map
Fits Guinier Rg for all 12 `GROUPS` conditions; `has_acid` from `SAMPLE_INFO`. `ACID_PROGRESSION` = `Y2F_h2o_kpi → Y2F_kpi_aa → kpi_aa` (updated after dropping `Y2F_kpi_kpi`; the first step is now a rougher, not-single-variable comparison since `Y2F_h2o_kpi`'s channel1 is water, not KPi).
**Figures**: `rg_by_condition.png`, `guinier_acid_progression.png`.
**Finding, unchanged and still the key result**: **all 12/12 conditions show Rg≈505–528 Å, r²≈0.98–0.99**, including `cap_empty` (a literally empty capillary). Not acid- or protein-specific — a shared instrumental/geometric artifact. Two real bugs in `src/saxs_metrics.guinier_fit` were found and fixed while establishing this (see Key Findings).

### `capillary_analysis.py` — separate pipeline for the static capillary screen (new)
Kept apart from the flow-channel scripts (different measurement modality: pre-made static samples, ~2 exposures each, no mixing-channel position; `compare_groups.load_scan_ids` was factored out of `load_group` specifically to support this). Covers all 14 capillary conditions (no `MIN_SCANS_PER_GROUP` filter), including the **third protein construct WT (wild-type/natural spider silk)** absent from the flow data, plus shear/no-shear/gel-state labels.
**Figures** (in `figures/capillary/`): `cap_rg_by_condition.png` (repeats the Rg-artifact check in this independent geometry), `cap_protein_kpi_noshear.png`, `cap_protein_kpi_shear.png`, `cap_protein_h2o.png` (Y2F vs YR2A vs WT), `cap_shear_yr2a/y2f/wt.png` (shear vs no-shear per protein), `cap_y2f_context.png` (Y2F across all its capillary conditions), `cap_buffer_acid.png` (KPi vs KPi+acid, no protein).
**Findings**: (1) the ~500 Å Rg artifact **also appears in every capillary condition**, including a literally empty capillary — stronger evidence it's beamline/detector-wide, not flow-cell-specific. (2) Curves are visually much cleaner here (tight 2-repeat overlap) than the flow data. Y2F vs YR2A vs WT look indistinguishable under matched sheared-KPi conditions; shear vs no-shear shows no visible effect for Y2F. Y2F in plain water shows a visibly higher WAXS bump than in KPi (sheared/unsheared/gelled), with KPi+acid intermediate — qualitative only, n=2 per condition, no statistics run.

### `drift_within_runs.py` — within-run self-consistency check (new, motivated the two drops above)
For every condition in `GROUPS` (run *before* the two drops, on all 14 original conditions): loads every kept frame, i0-normalized only (no q-region pinned, so drift can't be hidden or manufactured), plots I(q) colored by chronological scan order, tracks total intensity and WAXS peak height vs scan order, and computes the relative RMSD between the mean curve of the first half vs second half of a condition's own scans — a single "how much does this run's own data drift" number, independent of any cross-condition comparison.
**Figures** (in `figures/drift/`): `<condition>_drift.png` per condition, `drift_summary.png` (ranked bar chart).
**Finding**: **2/12 conditions were chaotic** (RMSD≥25%: `Y2F_kpi_kpi` 108.8%, `resin` 46.6%) — both subsequently dropped from `GROUPS` (see `compare_groups.DROPPED_CONDITIONS`). After the drop, **0/10 remaining conditions are chaotic**; most are well under 10% RMSD. Also resolved an apparent contradiction: `kpi_h2o` looks internally *stable* by this whole-curve metric (8.1% RMSD) even though its early/late split reproduces most of the `Y2F_kpi_h2o` "protein peak" — because the WAXS peak region is a tiny fraction of the whole curve, and for a protein-free buffer that peak metric is just noise oscillating around zero, so a large *relative* drift there is invisible in a *whole-curve* RMSD. Both checks are valid; they answer different questions (bulk curve reliability vs. one specific feature's reliability).

## Two real bugs found and fixed in `src/saxs_metrics`/`protein_vs_buffer.py` (matters for reuse)

1. **Baseline-anchoring bug in `peak_analysis`**: default baseline windows (q=0.03–0.08, 1.9–2.1) are fine for a raw curve but land inside the low-q Guinier-mismatch region on a *difference* curve, inflating significance by orders of magnitude (a spurious 56σ collapsed to <1σ once fixed via `DIFF_Q_BASE=((0.5,0.7),(1.5,1.7))`). This fix was initially applied only to the main peak test in `protein_vs_buffer.py`; a **second copy of the same bug** was later found in the replicate-level peak calculations (missing `q_base=DIFF_Q_BASE`) and fixed too — that's what made the printed "replicate variance" ranges look wildly inconsistent (`[-0.008, 0.47]`) before the fix vs. tightly one-signed after (`[0.006, 0.10]`).
2. **Guinier refinement instability** (fixed via `refine_min_points` guard, default 10, plus exposing `q_fit_max`/`n_points`): the two-stage fit could collapse onto as few as 5 points hugging the detector's lowest measurable q whenever the broad fit's Rg came out large. In this dataset it always lands on the same 10 points (q≈0.003–0.006 Å⁻¹) regardless of sample — which is *why* every condition reports the same ~500 Å Rg (the `map_llps_by_rg.py` finding).

## Overall state / open questions

1. **No LLPS (acid-specific large-Rg) signature** — the ~500 Å feature is universal, confirmed in both the flow-channel and the completely independent static-capillary geometry, and present with zero sample in the beam. To look for real LLPS you'd need to subtract/work around this shared feature first, or use a different technique (turbidity, imaging).
2. **No time-resolved (channel-position) effect** in any of the 9 flow-channel conditions.
3. **No protein-vs-blank WAXS peak survives both the normalization/scale-correction checks AND the drift check.** `Y2F_kpi_h2o` vs `kpi_h2o` remains the most interesting candidate (9.5σ scale-corrected, `i0` norm) but its blank (`kpi_h2o`) shows internal drift at the WAXS-peak level large enough (76–88% of the claimed effect) to not trust the result as-is. Next step if this is worth pursuing: drift-correct the blank (e.g. interpolate its early/late trend to the protein run's timestamp) before subtracting, rather than using its plain mean.
4. **Capillary screen** (`capillary_analysis.py`) has cleaner data and a third protein (WT) but only 2 repeats per condition — no statistics attempted; the qualitative water-vs-KPi WAXS bump difference for Y2F there might be worth a dedicated look if more capillary repeats can be collected.
5. Two conditions (`Y2F_kpi_kpi`, `resin`) are excluded from all analysis as internally unreliable — see `compare_groups.DROPPED_CONDITIONS` for why, and `drift_within_runs.py`'s figures for the visual evidence (a real mid-run event for the former, an apparent beam/position step-change for the latter).

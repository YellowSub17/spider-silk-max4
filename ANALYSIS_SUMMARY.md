# Spider silk SAXS/WAXS analysis — summary for handoff

**Note**: this replaces an earlier version of this summary. That version was built when `compare_groups.GROUPS` was hand-curated from free-text log comments in `data/datasetsTable2026-6-1.csv`, which turned out to have mislabeled several conditions (e.g. a group called "Y2F_buffer_pH8" that actually contained protein) and completely missed a rich capillary-screening dataset (including a third protein construct, "WT"). `GROUPS` has since been rebuilt from `data/scans_clean.csv`'s authoritative `Sample` column, and every downstream script rerun against the corrected groups. This version reflects that rebuild.

## Data & codebase

- `data/scans_clean.csv` — **authoritative** per-ScanID metadata. Key columns: `scanID`, `Sample` (ground-truth condition, format `[CHANNEL1]/[CHANNEL2]`, e.g. `Y2F+kpi/aa` = Y2F protein in KPi buffer, mixed in-channel with acetic acid; a protein-free run like `kpi/aa` is the matched blank), `description` (rotator "position N" label — position along the mixing channel / time since the two channels merged).
- `data/datasetsTable2026-6-1.csv` — the file used *before* this rebuild; no longer used for grouping. Its free-text `Comment` column is still a handy human-readable log if you want to cross-check something, but not authoritative for sample identity.
- `data/raw/scan-<id>*.h5`, `data/process/azint/scan-<id>_*_integrated.h5` — raw and azimuthally-integrated detector data (Eiger=SAXS, Pilatus=WAXS).
- `src/` (pre-existing, lightly extended this session):
  - `scan.py` / `combo_scan.py` — `Scan`/`ComboScan`: load and merge SAXS+WAXS into `self.qs`/`self.Is`; `Scan` computes a log-log linear-fit `r_squared` streak-quality metric over `det.LINEAR_SAXS_RANGE` (q=0.005–0.015 Å⁻¹).
  - `det.py`, `const.py` — detector geometry/q-ranges, paths.
  - `scale_fns.py` — added `norm_i0()` (normalize each frame by incident flux `i0`), alongside pre-existing `norm_qrange()`/`norm_max()`.
  - `saxs_metrics.py` — new: `guinier_fit` (Rg/I0/r, with `q_min`/`q_max`/`refine`/`refine_min_points` params — see Key Findings, this had a real bug), `porod_exponent`, `peak_analysis` (WAXS peak height/area/position via baseline subtraction, `q_base` configurable), `kratky`, `dimensionless_kratky`.

## `compare_groups.py` — the foundation everything else imports

Rebuilt to derive `GROUPS` directly from `scans_clean.csv`'s `Sample` column (`MIN_SCANS_PER_GROUP=3` filters out unusably small conditions). Also provides:
- `SAMPLE_INFO[name]` — parsed dict per group: `protein` (`Y2F`/`YR2A`/`WT`/`None`), `channel1`, `channel2`, `is_capillary`, `has_acid`.
- `BLANK_OF[name]` — the matching protein-free blank group name for a protein-containing group (matches on `{channel1}/{channel2}` or its reverse — the final mixed composition is the same regardless of which channel the protein started in), or `None` if no matching blank exists in the data.
- `RAW_NAME[name]` — the original (unsanitized) Sample string.
- `NORM_METHODS` (`linear_saxs`/`neutral`/`i0`), `load_group`, `kept_frames`/`mean_curve` (apply normalization + the `r_squared>0.999` streak filter), `coeff_of_variation` — all unchanged, all still generic over whatever's in `GROUPS`.

**Current `GROUPS` (14 conditions, ≥3 scans each):**

| Group name | Raw Sample | n scans | Protein | Acid? |
|---|---|---|---|---|
| `Y2F_kpi_aa` | Y2F+kpi/aa | 251 | Y2F | yes |
| `Y2F_h2o_kpi` | Y2F+h2o/kpi | 165 | Y2F | no |
| `kpi_h2o` | kpi/h2o | 145 | — | no |
| `Y2F_h2o_h2o` | Y2F+h2o/h2o | 120 | Y2F | no |
| `Y2F_kpi_kpi` | Y2F+kpi/kpi | 96 | Y2F | no |
| `kpi_aa` | kpi/aa | 82 | — | yes |
| `YR2A_h2o_kpi` | YR2A+h2o/kpi | 45 | YR2A | no |
| `Y2F_kpi_h2o` | Y2F+kpi/h2o | 42 | Y2F | no |
| `YR2A_h2o_h2o` | YR2A+h2o/h2o | 38 | YR2A | no |
| `h2o_h2o` | h2o/h2o | 30 | — | no |
| `resin` | resin | 15 | — | no |
| `cap_YF_kpi_shear` | cap:YF/kpi/shear | 7 | Y2F (capillary) | no |
| `cap_empty` | cap:empty | 4 | — (capillary) | no |
| `cap_YF_kpi_aa` | cap:YF/kpi/aa | 4 | Y2F (capillary) | yes |

`BLANK_OF`: `Y2F_kpi_aa`→`kpi_aa`, `Y2F_h2o_h2o`→`h2o_h2o`, `Y2F_h2o_kpi`→`kpi_h2o`, `YR2A_h2o_kpi`→`kpi_h2o`, `YR2A_h2o_h2o`→`h2o_h2o`, `Y2F_kpi_h2o`→`kpi_h2o`. **`Y2F_kpi_kpi` has no blank** (no `kpi/kpi` sample exists) — can't be background-subtracted.

There's also a much richer capillary screen in the raw data (2–7 scans each, below `MIN_SCANS_PER_GROUP` for most) including a **third protein construct "WT" (wild-type)** and shear/no-shear/gel-state labels (`cap:WT/kpi/shear`, `cap:YF/kpi/gel`, etc.) — present in `scans_clean.csv` but not big enough individually to analyze with this pipeline; worth a dedicated look if useful.

`COMPARISONS` (for `compare_groups.py`'s own plots) updated to the new names: buffer effect on Y2F (`Y2F_h2o_h2o`/`Y2F_h2o_kpi`/`Y2F_kpi_aa`), Y2F vs YR2A in water/KPi, buffer-only controls (`h2o_h2o`/`kpi_h2o`/`kpi_aa`), and the corrected acid progression (`Y2F_kpi_kpi`→`Y2F_kpi_aa`→`kpi_aa`).

## Scripts and their figures

### `compare_groups.py`
`plot_comparison()` per `COMPARISONS` entry → `buffer_effect_on_y2f.png`, `y2f_vs_yr2a_in_water.png`, `y2f_vs_yr2a_in_kpi.png`, `buffer-only_controls.png` (renamed from `..._over_time.png`), `acetic_acid_progression_on_y2f.png`.
**Finding**: run-to-run CV is large everywhere (12–68%; `Y2F_kpi_kpi` is worst at 68%) — the first-pass signal that any condition difference needs checking against noise.

### `position_series.py` — channel-position (mixing-time) effect
`CHANNEL_GROUPS` is now every non-capillary, non-`resin` group in `GROUPS` (9 of them) — built generically via `SAMPLE_INFO`, not a hardcoded list. Position labels also now sourced from `scans_clean.csv`'s `description` column (previously a different CSV), so sample-identity and position labels are guaranteed consistent.
**Figures**: `<group>_position_curves.png`, `<group>_position_metrics.png` for all 9 groups.
**Finding**: still **no robust channel-position trend** — only 2/54 metric×group tests hit p<0.05 (expected false-positive rate), same conclusion as before the rebuild.

### `protein_vs_buffer.py` — protein vs matched blank, with reliability diagnostics
`PAIRS` is now built programmatically from `BLANK_OF` (6 pairs, up from 4 hand-picked ones): `Y2F_kpi_aa`, `Y2F_h2o_h2o`, `Y2F_h2o_kpi`, `YR2A_h2o_kpi`, `YR2A_h2o_h2o`, `Y2F_kpi_h2o`, each vs its matched blank. Diagnostics unchanged: CV ratio, scaling-factor fit (neutral q=0.1–1.0 Å⁻¹ window), peak significance (raw + scale-corrected, baseline anchored near the peak — see Key Findings), replicate envelope, relative difference — each run under `linear_saxs`/`neutral`/`i0` normalization.
**Figures**: `<pair>__<norm>.png`, 18 total (6 pairs × 3 norms).
**Finding**: results are noisy and normalization-sensitive as before, with one new, more solid result: `Y2F_kpi_h2o` vs `kpi_h2o` shows a scale-corrected peak significance that's consistent and elevated under `i0` (9.5σ) and present (if smaller) under `neutral` (2.6σ) — worth a closer look. Nothing else clears the bar consistently across all three normalizations.

### `kratky_plot.py` — Kratky-shape diagnostic
Runs on the same 6 `PAIRS`. **Figures**: `<pair>__kratky.png`.
**Finding**: unchanged conclusion — Guinier Rg (q<0.02) is ~505–528 Å for every curve regardless of pair, protein, or blank; a wider q<0.1 fit gives a common ~36–39 Å instead, also shared between protein and blank, with poor linearity (r²≈0.65–0.70). Neither window separates protein from blank.

### `map_llps_by_rg.py` — systematic Rg-artifact map (the key result)
Fits Guinier Rg for every group in `GROUPS` (now 14, up from 20 — some previously-included ranges no longer exist as distinct groups under the corrected labeling); `has_acid` now read from `SAMPLE_INFO` instead of a substring guess. `ACID_PROGRESSION` corrected to `Y2F_kpi_kpi → Y2F_kpi_aa → kpi_aa` (same buffer throughout, acid added, then protein removed — a real single/double-variable progression, unlike the old mislabeled one).
**Figures**: `rg_by_condition.png`, `guinier_acid_progression.png`.
**Finding, reconfirmed and strengthened**: all 14/14 conditions show Rg≈505–528 Å with r²≈0.98–0.99, **including `cap_empty` — a literally empty capillary with no sample at all**, and `resin` (channel/resin material, no liquid). This is about as strong a smoking gun as it gets: the ~500 Å feature cannot be protein- or acid-related; it's a shared instrumental/geometric artifact common to every acquisition on this setup.

## Two real bugs found and fixed in `src/saxs_metrics.guinier_fit` (matters for anyone reusing it)

1. **Baseline-anchoring bug in `peak_analysis`** (fixed via `protein_vs_buffer.py`'s `DIFF_Q_BASE=((0.5,0.7),(1.5,1.7))` override): the default baseline windows (q=0.03–0.08, 1.9–2.1) are fine for a raw curve but land inside the huge low-q Guinier-mismatch region when applied to a *difference* curve, inflating peak significance by orders of magnitude (a spurious 56σ collapsed to <1σ once fixed).
2. **Guinier refinement instability** (fixed via `refine_min_points` guard, default 10, plus exposing `q_fit_max`/`n_points` in the returned dict): the two-stage fit (broad window → refine to q·Rg<1.3) could collapse onto as few as 5 points hugging the detector's lowest measurable q whenever the broad fit's Rg came out large. In this dataset it always lands on exactly the same 10 points near qmin (q≈0.003–0.006 Å⁻¹) regardless of sample — which turned out to be *why* every condition reports the same ~500 Å Rg (see `map_llps_by_rg.py` finding above). The guard doesn't change the conclusion (the 10-point window is genuinely the best-supported fit available at that Rg scale), but it does make the fit reproducible/inspectable rather than silently unstable, and `map_llps_by_rg.py` now reports `n_points` alongside `Rg`/`r²` so this is visible rather than hidden.

## Overall state / open questions

1. **No LLPS (acid-specific large-Rg) signature** — the ~500 Å feature is universal, present even with no sample in the beam path at all. If you want to look for real LLPS, you need to work around/subtract this shared feature first, or use a technique other than absolute low-q Guinier fitting on this setup (e.g. turbidity, imaging).
2. **No time-resolved (channel-position) effect** found in any of the 9 flow-channel conditions.
3. **One moderately promising, not-yet-fully-vetted lead**: `Y2F_kpi_h2o` (Y2F pre-dissolved in KPi, mixed with water in-channel) vs its `kpi_h2o` blank shows an elevated, cross-normalization-present peak significance (up to 9.5σ under `i0`, the most assumption-free normalization) — this is the best surviving candidate for a real protein WAXS signature in the whole dataset and could reward a closer, dedicated look (e.g. full nonlinear peak fit with an explicit background model instead of baseline subtraction).
4. **An unanalyzed capillary screen exists** in `scans_clean.csv` with a third protein construct (WT) and shear/no-shear/gel conditions — each individual condition has too few scans (2–7) for this pipeline's per-condition statistics, but might be worth pooling or examining qualitatively if that construct/those conditions matter to the underlying question.




## scripts

[compare_groups.py] Foundational module: builds GROUPS directly from data/scans_clean.csv's Sample column (format [CHANNEL1]/[CHANNEL2]), parsed into SAMPLE_INFO (protein/channel1/channel2/is_capillary/has_acid) and auto-matched protein-free blanks via BLANK_OF. Provides the shared loading/normalization backbone (load_group, load_scan_ids, kept_frames, mean_curve, three normalization schemes: linear_saxs/neutral/i0) that every other script imports. Its own comparison plots show run-to-run CV is large everywhere (12–68%), meaning any condition difference needs statistical scrutiny before trusting it — the throughline for everything that follows.

[position_series.py] Tests for a channel-position (mixing-time) effect: splits each of the 9 flow-channel conditions by rotator "position N" and fits Rg, Porod exponent, and WAXS peak metrics vs. position. Result: no significant position trend anywhere — only 2/54 metric×condition tests hit p<0.05, consistent with pure chance. No time-resolved mixing signature detected at this position spacing.

[protein_vs_buffer.py] Core reliability-diagnostic suite: for each of 6 protein-vs-matched-blank pairs, computes CV ratio, a scaling-factor fit (flags normalization mismatch), WAXS peak significance (both raw and scale-corrected, using a baseline anchor fixed near the peak after finding the metric's defaults were unreliable on difference curves), a per-replicate envelope, and a relative-difference curve — each repeated across all 3 normalization schemes. Result: nothing survives consistently. The one number that looked strong early on (Y2F/KPi, 11.5σ) collapsed to <1σ once normalization and baseline artifacts were controlled for. The most promising surviving lead is Y2F_kpi_h2o vs its kpi_h2o blank (up to 9.5σ under the most assumption-free normalization, i0) — worth a closer look but not yet confirmed.

[kratky_plot.py] Global-shape (Kratky) diagnostic, chosen specifically because it doesn't depend on picking one q-window to trust the way single-peak testing does. Result: uninformative here — Guinier Rg (q<0.02) comes out ~505–528 Å for every protein and blank alike (see map_llps_by_rg.py), which collapses the dimensionless Kratky plot into an unreachable corner; the raw Kratky plot is dominated by the un-decayed high-q edge, a known artifact when data extends into WAXS. A wider q<0.1 Guinier window gives a common ~36–39 Å for everything too, with poor linearity (r²≈0.65–0.70) — no separable protein signal at any window tried.

[map_llps_by_rg.py] The key result of the whole analysis. Systematically fits Guinier Rg across all 14 flow-channel + capillary-background conditions to test whether the ~500 Å apparent Rg is an LLPS (droplet) signature specific to acetic-acid conditions. Result: it is not. All 14/14 conditions show Rg≈505–528 Å with r²≈0.98–0.99, identically whether or not acid or protein is present — including resin (pure channel material, no liquid). Along the way this surfaced and fixed two real bugs in src/saxs_metrics.guinier_fit: a baseline-anchoring bug that was inflating peak-significance numbers by orders of magnitude on difference curves, and a refinement-window instability that was silently collapsing every fit onto the same 10 points nearest the detector's lowest measurable q — which turned out to be why every condition reports the same Rg.

[capillary_analysis.py] Separate pipeline for the static capillary screen (a different measurement modality — pre-made samples, ~2 exposures each, no flow/position), which also uniquely contains a third protein construct, WT (wild-type/natural spider silk). Repeats the Rg-artifact check in this independent geometry: the same ~500 Å artifact reappears in all 14 capillary conditions, including a literally empty capillary — stronger evidence the artifact is beamline/detector-wide, not specific to the flow cell. Direct comparisons here are visually much cleaner (tight repeat overlap) than the flow data: Y2F vs YR2A vs WT are indistinguishable under matched sheared-KPi conditions, and shear vs. no-shear shows no visible effect for Y2F — but Y2F in plain water shows a noticeably higher WAXS bump than in KPi (sheared/unsheared/gelled), with KPi+acetic acid intermediate. These are qualitative observations only (n=2 per condition, no statistics run).

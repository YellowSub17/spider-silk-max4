# Protein vs. Buffer Comparison Figures

This folder contains detailed diagnostic plots for each protein-vs-blank pair, showing the raw curves, difference curves, and reliability diagnostics.

## Figure Naming Convention

Files are named: `<protein>_vs_<blank>__<normalization>[__matched].png`

- **protein** = The protein-containing sample (e.g., `y2f_h2o_kpi`)
- **blank** = The matched protein-free buffer (e.g., `kpi_h2o`)
- **normalization** = How the curves are scaled
  - `i0` = Incident flux normalization
  - `linear_saxs` = Low-q region normalization
  - `neutral` = Mid-q region normalization
- **matched** suffix = Position-matched subtraction (only for pairs with sufficient spot overlap)

**Example**: `y2f_h2o_kpi_vs_kpi_h2o_blank__i0.png` shows Y2F in KPi buffer vs. KPi buffer alone, normalized by incident flux.

## The Three-Panel Diagnostic Plot

### Panel 1: Raw Curves (top, log-log scale)
Shows the mean intensities before subtraction.

- **Red line + shading** = Protein sample mean curve, ±1σ variability (std dev across frames)
- **Blue line + shading** = Blank sample mean curve, ±1σ variability
- **Gray vertical band** = Scaling-fit window (q = 0.1–1.0 Å⁻¹)
  - This region is used to fit an optimal scale factor `s` such that `s * blank ≈ protein` over this q-range
  - The assumption: this region should be identical between samples (solvent-dominated scattering)

**What to look for**:
- Do the curves overlap? If yes, no visible protein signal.
- Do they diverge uniformly (parallel on log-log)? Suggests a scaling difference, not a localized feature.
- Is there a localized bump (WAXS peak at q ≈ 1–1.5 Å⁻¹) in one curve but not the other? Potential protein signal.

### Panel 2: Difference & Diagnostics (middle, symlog scale)
The key panel—this is where you judge whether a "protein signal" is real.

**Three traces**:
1. **Green solid line** (`mean diff`) = `I_protein - I_blank`
   - The raw difference curve, what you'd get by simple subtraction
   - Green shading = ±1 SEM (standard error of the mean across frames)

2. **Gray filled region** (`replicate envelope`, 5th–95th percentile)
   - Shows what single-replicate differences look like
   - Computed as: each individual protein frame minus the blank mean, AND the protein mean minus each individual blank frame
   - **If this envelope encompasses the green line, the signal is within noise** — any single run alone could produce that peak

3. **Orange dashed line** (`scale-corrected diff`) = `I_protein - s * I_blank`
   - Subtracts the *scaled* blank instead of the raw blank
   - The scale factor `s` is printed (e.g., `s = 0.98`)
   - **This is crucial**: If a peak disappears after scale correction, it was an overall intensity difference (likely normalization), not a real protein feature
   - If the peak **persists and is similar**, the signal is localized and independent of overall scaling

**Scale correction explained**:
- Suppose `I_blank` is systematically 10% lower than `I_protein` across all q (e.g., due to concentration, geometry, or normalization artifacts)
- A naive difference `I_protein - I_blank` would have a broad, featureless ~10% offset everywhere
- Scale correction fits `s ≈ 0.91` and computes `I_protein - 0.91 * I_blank`, which brings the overall levels into agreement
- Now any remaining peak in the diff curve is a *localized* feature, not an offset

**Interpretation**:
- **Peak visible in green AND orange, similar position and height**: Real, localized protein signal
- **Peak visible in green but disappears in orange**: Artifact of overall scaling mismatch
- **Peak smaller than replicate envelope**: Within noise; not resolvable
- **Peak at the edge of the q-window (q ≈ 0.9 or 1.9 Å⁻¹)**: Likely baseline artifact, not real

### Panel 3: Relative Difference (bottom)
`(protein - blank) / blank` — the fractional difference.

- Less sensitive to overall scaling artifacts
- Easier to spot when one condition is systematically higher/lower across q
- Example: If protein is 5% higher than blank everywhere, this shows a flat ~0.05 line

---

## Interpreting the Diagnostics (Console Output)

When you run `protein_vs_buffer.py`, each pair prints 8 diagnostics:

1. **CV ratio** (Coefficient of Variation): protein_CV / blank_CV
   - Should be < 2.0
   - If > 2, protein is noisier than blank; "signal" may be detecting instability

2. **Scaling factor**: The `s` value printed on the orange dashed line
   - Should be close to 1.0 (±5%)
   - If s = 0.5 or 2.0, the groups differ drastically in overall intensity

3. **Peak significance (raw diff)**: Reported in sigma (σ)
   - `peak_z >= 3 sigma` = resolvable above noise
   - `peak_z < 3 sigma` = below threshold, don't trust it

4. **Replicate variance in diff**: The range of peak heights across single-replicate diffs
   - If the range is huge, the "peak" is noise

5. **Peak significance (scale-corrected)**: The peak after accounting for overall scaling
   - **This is the more credible number** if scale ≠ 1
   - Must be ≥ 3σ to believe the result

6. **Drift check**: Split the blank into early-vs-late halves, run the identical test
   - If drift alone produces a peak ≥ 50% the size of the claimed protein peak, the result is unreliable
   - Indicates the difference is session drift, not protein

7. **Position-matched subtraction** (if applicable):
   - Groups frames by physical stage spot, matches protein and blank separately per spot
   - Reports shared spots, coverage %, and per-spot peak significance
   - A peak present at all shared spots is real; one at only one spot is a geometry artifact

8. **Normalization comparison**: A table across `i0`, `linear_saxs`, and `neutral` norms
   - Does the signal persist across norms? If yes, more credible

---

## Case Studies

### Y2F_h2o_h2o vs h2o_h2o (No Signal)
- Curves overlay almost perfectly (panel 1)
- Difference is tiny (panel 2, green line near zero)
- Replicate envelope dwarfs the diff (panel 2, gray >> green)
- Conclusion: **No protein signal; difference is noise**

### Y2F_h2o_kpi vs kpi_h2o (Weak Signal)
- Protein curve slightly higher than blank (panel 1)
- Small high-q peak visible in green (panel 2)
- Orange line shows scale factor ≈ 0.98 (nearly 1.0, good sign)
- Peak significance ≈ 8σ raw, ≈ 8σ scale-corrected (marginal)
- Position-matched version with full spot coverage: **10.5σ** — credible

### Y2F_kpi_aa vs kpi_aa (No Comparison Possible)
- Protein and blank don't share any physical stage spots (different x coordinates by ~0.28 mm)
- Position-matched subtraction is **skipped** (not forced to a nearest-spot match)
- Pooled result is unreliable because the two groups sampled different parts of the mixing channel
- **Conclusion**: Need new data collected at matching coordinates

---

## Using These Figures for a Presentation

1. **Start with panel 1** (raw curves): "Do you see a visible difference?"
2. **Move to panel 2** (difference curve): "Is the peak real or noise? Is it a localized feature or just scaling?"
3. **Check the diagnostics**: "Peak significance, drift check, replicate envelope—is this above the bar?"
4. **Show the position-matched version** (if available): "The signal persists and grows with distance—not a geometry artifact"

This progression convinces a skeptical audience that you've done due diligence before claiming a protein signal.


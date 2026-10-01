# Buffer Subtraction & Normalization Methods

## Overview

The goal of protein-vs-buffer analysis is to isolate the X-ray scattering contribution from the protein itself by subtracting a matched protein-free blank (buffer only) from the protein sample. However, this is not a simple subtraction—careful normalization is required to account for differences in beam intensity, sample composition, and instrumental drift.

## The Challenge

Raw scattering intensities I(q) depend on:
- **Protein structure** (what we want)
- **Buffer/solvent scattering** (want to remove)
- **Incident beam intensity (i0)** (instrumental, not sample)
- **Sample positioning** (affects total collected intensity)
- **Detector geometry** (fixed, but affects low-q region)

Without proper normalization, the "protein signal" you extract is contaminated by these confounders, and false positives are common.

## Three Normalization Strategies

### 1. **linear_saxs** (Low-q normalization)
- **What it does**: Normalizes each frame by the sum of intensity in the low-q "linear SAXS" region (q = 0.005–0.015 Å⁻¹)
- **Assumption**: This q-region should be the same between protein and buffer; differences there reflect only normalization mismatch
- **Best for**: Stable regions where protein and buffer truly scatters identically
- **Weakness**: If the protein actually *changes* scattering in the low-q Guinier region (due to concentration, hydration, or aggregation), this method over-corrects

### 2. **neutral** (Mid-q normalization)
- **What it does**: Normalizes each frame by the sum of intensity in the mid-q "neutral" region (q = 0.1–1.0 Å⁻¹)
- **Assumption**: This region is dominated by solvent/background scattering, not protein structure
- **Best for**: Detecting local (e.g., WAXS) features without assuming the low-q region is identical
- **Weakness**: Assumes the buffer scattering profile is flat in this range, which may not hold

### 3. **i0** (Incident flux normalization)
- **What it does**: Normalizes each frame by its incident flux (beam intensity delivered to the sample)
- **Assumption**: None about which q-region should match; just corrects for beam fluctuations
- **Best for**: No assumptions about sample equivalence; purely a beam-intensity correction
- **Strength**: Most model-independent, best for comparing conditions with potentially different bulk properties
- **Weakness**: If sample loading/position varies, this won't correct for geometric effects

## Figure Collection: `01_normalization_overview/`

**File**: `sample_blank_curves__i0.png` (and variants for other norms)

This shows all 6 protein-vs-blank pairs side-by-side, normalized by i0 (or the specified method). Each panel plots:
- **Red line**: Protein sample (n frames kept after quality filter)
- **Blue line**: Matched blank (buffer only)
- **Gray shaded region**: ±1σ run-to-run variability

### What to look for:
1. **Curves track together**: If protein and blank curves overlap (especially in low-q), the protein may be largely invisible—or the normalization is over-correcting.
2. **Systematic separation**: If one curve is consistently higher or lower across q, check the normalization assumption is valid.
3. **High-q bump**: The ~1 Å⁻¹ peak is a WAXS feature; it's visible in all samples because the diffraction patterns are dominated by buffer structure, not protein.

---

## Protein-vs-Buffer Analysis: `02_protein_vs_blank_comparisons/`

Each pair has **up to 7 diagnostic figures** per normalization method:

### Figure Structure (Example: `y2f_h2o_kpi_vs_kpi_h2o_blank__i0.png`)

A three-panel plot:

#### Panel 1: Raw Curves (top)
- **Red + shaded region**: Mean protein curve ± 1σ (std dev across frames)
- **Blue + shaded region**: Mean blank curve ± 1σ
- **Gray shaded region**: "Scaling-fit window" (q = 0.1–1.0 Å⁻¹), the neutral region used to fit a scale factor

**Why this matters**: If the two curves are nearly identical, subtracting them will give near-zero result. If they diverge, you need to check whether it's real protein signal or a normalization artifact.

#### Panel 2: Difference Curve (middle)
Shows three traces:

1. **Green solid line** (`mean diff`): `I_protein - I_blank`
   - This is the claimed "protein signal"
   - Green shading: ±1 SEM (standard error of the mean)

2. **Orange dashed line** (`scale-corrected diff`): `I_protein - scale * I_blank`
   - Scale factor `s` is fitted over the neutral window (q=0.1–1.0) to minimize the overall difference
   - If s ≠ 1, the two groups differ in overall intensity; this corrects for that before looking for a *localized* peak
   - **This is the key**: A true protein feature should survive scale correction. If it disappears after scaling, it was just an overall intensity difference (likely normalization, not protein).

3. **Gray shaded region** (`replicate envelope`): The 5th–95th percentile of all single-replicate differences
   - Each individual protein frame compared to the mean blank (and vice versa)
   - **If this envelope is much larger than the mean diff, the "signal" is just noise** — run-to-run variability could alone produce that peak

#### Panel 3: Relative Difference (bottom)
- `(protein - blank) / blank`
- Less sensitive to an overall scaling mismatch than the absolute difference
- Easier to spot when one condition is systematically higher/lower

### Reading the Diagnostics

Each pair prints a text summary. Key flags:

- **CV ratio** (Coefficient of Variation protein / blank): Should be < 2.0. If protein CV >> blank CV, the "signal" may just be detecting instability, not structure.
- **Scaling factor** `s`: Should be near 1.0 (±5%). If s = 0.5 or 2.0, the two groups' overall intensities differ drastically—check normalization.
- **Peak significance**: Reported in sigma (σ). Need ≥ 3σ to be resolvable above noise (analogous to 3σ in particle physics).
- **Replicate envelope**: If the envelope encompasses the diff curve, the signal is within the noise of single-replicate variation.
- **Drift check**: The blank condition is split into early-vs-late halves, and the identical test is run with no protein. If that produces a peak ≥ 50% as large as the claimed protein peak, the result may be session drift, not protein.

---

## Why Most Pairs Show "No Protein Signal"

**Example: Y2F_h2o_h2o vs h2o_h2o**

Looking at the first figure above:
- Panel 1: Protein (red) and blank (blue) curves nearly overlay
- Panel 2: The difference is tiny (green line hugs zero)
- Panel 2: The replicate envelope (gray) is enormous — single-replicate diffs scatter widely
- **Conclusion**: No resolvable protein signal; the difference is within run-to-run noise

This is **expected**. Protein scattering at low q is dominated by Guinier (diffuse) scattering from the protein's overall size. A 100 kDa protein in a sample with similar-sized solvent molecules produces only ~1–2% relative intensity difference in most q regions. Measuring that 1–2% difference requires:
1. Very high signal-to-noise (good beam, stable sample)
2. Excellent run-to-run reproducibility (tight temperature, flow control)
3. Large number of replicates to beat down statistical noise

This dataset does not have all three; hence most pairs are null.

---

## The One Exception: Y2F_kpi_h2o vs kpi_h2o

**Example: Y2F_h2o_kpi vs kpi_h2o_blank (second figure above)**

- Panel 1: Protein (red) is *consistently higher* than blank (blue) across most of the curve
- Panel 2: Green line (mean diff) shows a clear high-q (WAXS) peak around q = 1–1.5 Å⁻¹
- Panel 2: Orange dashed line (scale-corrected) shows the same peak, slightly moved
  - Scale factor s = 0.98 (nearly 1.0, so no over-all scaling artifact)
  - This strengthens the case: the peak isn't just a scaling mismatch
- Panel 2: Replicate envelope is much smaller than the peak; single replicates don't alone produce this
- **Conclusion**: This pair shows a 10.5σ peak under position-matched analysis (see next section)

This is the **single most credible positive result** in the entire dataset.

---

## Position-Matched vs. Pooled Subtraction

### Pooled Subtraction (shown above)
- Compute mean I(q) of all protein frames: `<I_p>`
- Compute mean I(q) of all blank frames: `<I_b>`
- Difference: `<I_p> - <I_b>`
- **Assumption**: Both groups visited the same physical spots in the mixing channel with similar frequency

### Position-Matched Subtraction
- Group protein frames by physical stage spot (e.g., y = 0.0 mm, y = 0.5 mm, ...)
- Group blank frames by the same spots
- For each *shared* spot, compute: `<I_p(spot)> - <I_b(spot)>`
- **Equal-weight average** of all shared spots
- **Why**: If protein frames cluster at spots 1–3 and blank frames cluster at spots 3–5, pooled subtraction mixes the position-mismatch into the result. Position-matched removes this bias by construction.

**Finding**: For Y2F_kpi_h2o vs kpi_h2o, position-matched subtraction:
- Increases the peak from 8.1σ (pooled, i0 norm) → 10.5σ (matched, i0 norm)
- Shows the signal is concentrated at the two *furthest-downstream* spots (y ≈ 3–5.5 mm)
- Signal is absent at the earliest spots (y ≈ 0–1 mm)
- **Pattern**: Signal *increases with distance*, consistent with a real mixing-time effect (more time in the channel = more structural change), not a single-spot geometry artifact

---

## Interpretation Summary

| Sample | Result | Explanation |
|--------|--------|-------------|
| Y2F_h2o_h2o vs h2o_h2o | No signal | Protein/buffer curves overlap; difference within replicate noise |
| Y2F_h2o_kpi vs kpi_h2o | Weak signal | Tiny peak, borderline significance |
| **Y2F_kpi_h2o vs kpi_h2o** | **Strong, position-dependent** | **10.5σ; grows downstream; real mixing-time effect** |
| Y2F_kpi_aa vs kpi_aa | Skipped | No shared physical spots between protein and blank (different x coordinates) |
| YR2A pairs | No signal | Similar to Y2F pairs; structure dominates noise only in one case |

---

## Next Steps: Replication

The Y2F_kpi_h2o vs kpi_h2o signal is compelling but is only **one pair**. To confirm:
1. Collect additional Y2F_kpi_h2o and kpi_h2o repeats
2. Balance spatial coverage (ensure both groups visit all spots equally)
3. Interleave collection temporally (alternate sample/blank) to avoid session drift
4. Test whether the downstream-increasing pattern holds


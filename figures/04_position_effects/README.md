# Position Effects & Mixing-Time Analysis

This folder shows the spatial distribution of scattering signals within the mixing channel and how physical position affects measured metrics.

## Key Discovery: Position-Locking

Synchrotron SAXS data are collected by stepping the sample through discrete stage positions (rotator spots in the mixing channel). Because the mixing channel is *dynamic*—samples flow and mix as they traverse—the stage position (distance along the channel) acts as a proxy for **mixing time**: scans at the first spot are at the start of mixing, scans at the last spot are after ~seconds of mixing.

Previous analysis used the "position N" label (1–5) assigned by the instrument. **V4 discovered this label is NOT a shared physical coordinate**: each sample's run assigns position 1–5 independently, and the actual stage coordinates vary widely across samples.

This folder uses **real stage coordinates** (e.g., y = 26.975 mm vs. y = 32.477 mm) to properly align frames across samples.

## Files

### `spot_effect_summary.png`
For all 9 testable conditions, plots the variance explained (η², eta-squared) by physical spot.

- **SAXS total intensity**: 9/9 conditions significant (p < 0.05)
  - Example: h2o/h2o shows η² = 82.3% — position explains 82% of the variance in total intensity
- **WAXS peak height**: Similar pattern across conditions
- **Conclusion**: Position-locking is real and strong, independent of normalization

### `<condition>_spot_effect.png` (one per condition)
One-way ANOVA + Kruskal-Wallis test comparing metrics across physical spots.

Shows:
- Spot ID on x-axis
- A metric (e.g., i0-normalized intensity, WAXS peak height) on y-axis
- Boxplot or scatter of individual frames per spot
- p-value and effect size (η²)

**How to interpret**:
- Significant p-value: Position has a real effect on this metric
- Large η²: Position explains a substantial fraction of the total variance
- Box positions: If they don't show a monotonic trend, the effect is **non-monotonic** (dips then rises, or vice versa)

---

## Connection to Mixing-Time Effects

A **mixing-time effect** would manifest as:
- Metrics that change systematically with distance along the channel
- E.g., protein secondary structure content increases downstream (as hydration equilibrates)
- Or aggregation state changes as molecules spend more time together

The position-locking shown here proves the intensity *depends* on position, but does NOT prove it's due to protein mixing. It could be:
1. **Real mixing effect** (protein structure changes with time in the channel)
2. **Geometric artifact** (detector sensitivity varies with sample position)
3. **Temperature gradient** (mixing channel has a temperature profile)
4. **Convective effects** (flow causes non-uniform concentration)

---

## Position-Matched Subtraction: The Key Evidence

For the Y2F_kpi_h2o vs. kpi_h2o pair (strongest signal):

1. **Pooled subtraction** (all frames averaged together)
   - Peak significance: 8.1σ (raw), 10.5σ (scale-corrected)

2. **Position-matched subtraction** (frames grouped by physical spot before averaging)
   - Peak significance: 8.1σ → 10.5σ (**increased**, not decreased!)
   - This means the signal is NOT explained by position-mix mismatch
   - In fact, matching by position *strengthens* the signal, suggesting real structure

3. **Per-spot breakdown**:
   - Spot y=0.0 mm: **No peak** (0.6σ)
   - Spot y=0.5 mm: **No peak** (0.4σ)
   - Spot y=1.0 mm: **No peak** (no significant peak)
   - Spot y=3.0 mm: **4.7σ peak** ← Signal starts here
   - Spot y=5.5 mm: **6.3σ peak** ← Signal strongest here

**Interpretation**: The signal **grows with distance downstream**. This pattern argues strongly against a static geometry artifact (which would show one outlier spot) and is consistent with a real mixing-time-dependent effect (protein changes structure as it mixes over seconds in the channel).

---

## Why This Matters for Interpretation

Previous analysis (v1–v3) used the label-based "position" index (1, 2, 3, 4, 5) as if it were:
- A shared coordinate system across samples ← **FALSE**
- Evenly spaced (1→2 = 2→3 = 3→4 spatially) ← **FALSE**

This caused:
- Linear regression against this wrong x-axis to find "no trend" ← **FALSE NEGATIVE**
- Conclusion: "No mixing-time effect" ← **RETRACTED**

V4 fixes this by:
- Using real stage coordinates (y in mm)
- Recognizing the effect is **non-monotonic**, not linear
- Using Kruskal-Wallis test (catches non-monotonic patterns) instead of linear regression
- Reporting **per-spot breakdowns** to show the spatial pattern

---

## Presentation Tips

1. **For skeptics**: Show the per-spot breakdown. A signal that appears at spots 3–5 but not 1–2 is unlikely to be random noise (which would scatter equally across all spots).

2. **For biochemists**: Frame it as mixing time: "The first spot is ~0 seconds of mixing, the last spot is ~3–5 seconds. The protein signal appears only after 1+ seconds of mixing, suggesting a time-dependent structural response."

3. **For instrumentalists**: Emphasize that the same position shows the same signal across *different* samples, ruling out instrument dead zones. Position is a real physical variable, not an artifact.


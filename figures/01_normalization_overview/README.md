# Normalization Overview

This folder shows all protein-vs-buffer pairs under three different normalization schemes, displayed side-by-side for comparison.

## Understanding Unnormalized Data First

**→ START HERE**: Read [README_UNNORMALIZED.md](README_UNNORMALIZED.md) to understand why normalization is critical and what raw data looks like.

In brief: Raw X-ray scattering intensities are **not directly comparable** between samples because they're contaminated by:
- Beam intensity fluctuations (±10-30%)
- Sample positioning/absorption differences
- Geometric variations

Normalization removes these instrumental artifacts so the true protein signal becomes visible.

---

## The Three Normalization Methods

### Files

- **sample_blank_curves__i0.png** — Incident flux normalization
  - Each frame normalized by its beam intensity (i0)
  - Most model-independent; removes beam fluctuations only
  - **No assumptions** about which q-region should match
  
- **sample_blank_curves__linear_saxs.png** — Low-q normalization  
  - Normalizes to the low-q region (q = 0.005–0.015 Å⁻¹)
  - Assumes protein/buffer have identical low-q "Porod" scattering
  - **Risk**: Over-corrects if protein actually changes low-q structure
  
- **sample_blank_curves__neutral.png** — Mid-q normalization
  - Normalizes to the neutral mid-q region (q = 0.1–1.0 Å⁻¹)
  - Assumes this region is solvent-dominated (same in both samples)
  - **Safer**: Less likely to hide real low-q differences

## How to Read These Plots

Each panel shows a protein-vs-blank pair at the same scale (after normalization):
- **Red line** = Protein (with ±1σ shading)
- **Blue line** = Matched blank (with ±1σ shading)
- **Gray shaded region** = Scaling-fit window (q = 0.1–1.0 Å⁻¹)
  - Used to estimate the optimal scale factor between the two curves
  - If curves are perfectly scaled, scale = 1.0

**Key comparison**: Look across the three normalizations:
- Do curves remain nearly identical under all three norms? → Protein is truly invisible
- Does separation appear only in one norm? → Likely an artifact of that norm's assumption
- Is separation consistent across all three? → More credible; not a normalization artifact

## Looking for Real Protein Signal

A true protein signal should show:
1. **Visible separation** between protein and blank curves
2. **Persistence across normalizations** — same signal in i0, linear_saxs, AND neutral
3. **A localized WAXS-region peak** (high q, ~1–1.5 Å⁻¹) — not just a broad low-q shape change
4. **Reproducibility** — appears in position-matched analysis and across multiple replicates


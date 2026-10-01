# Normalization Overview

This folder shows all protein-vs-buffer pairs under three different normalization schemes, displayed side-by-side for comparison.

## Files

- **sample_blank_curves__i0.png** — Incident flux normalization
  - Each frame normalized by its beam intensity (i0)
  - Most model-independent; assumes only beam fluctuations
  
- **sample_blank_curves__linear_saxs.png** — Low-q normalization  
  - Normalizes to the low-q region (q = 0.005–0.015 Å⁻¹)
  - Assumes protein/buffer have identical low-q scattering
  
- **sample_blank_curves__neutral.png** — Mid-q normalization
  - Normalizes to the neutral mid-q region (q = 0.1–1.0 Å⁻¹)
  - Assumes this region is solvent-dominated, not protein-specific

## How to Read

Each panel shows a protein-vs-blank pair at the same scale:
- **Red line** = Protein (with ±1σ shading)
- **Blue line** = Matched blank (with ±1σ shading)
- **Gray shading** = Scaling-fit window (q = 0.1–1.0 Å⁻¹), used to estimate the relative scale factor between the two curves

**Key observation**: Compare across the three normalization methods. If the protein and blank curves remain nearly identical under all three norms, the protein is truly invisible (or the normalization differences don't matter). If one norm separates them while others don't, investigate which assumptions are being violated.

## Looking for Signal

A true protein signal should:
1. Be visible in the **raw curves** (panel 1)—the protein curve should be noticeably higher or lower than the blank
2. Persist across **multiple normalization methods**—if signal appears only under one norm, it's likely an artifact of that norm's assumption
3. Have a **clear WAXS-region peak** (high q, ~1–1.5 Å⁻¹), not just a low-q shape change


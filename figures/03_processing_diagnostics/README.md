# Processing Diagnostics

This folder contains figures that verify data quality and check for systematic issues like instrumental drift and invariant structural features.

## Files

### `rg_by_condition.png`
Shows the radius of gyration (Rg) fitted in two q-windows for every condition.

- **Blue points** (q < 0.02 Å⁻¹): Low-q Guinier fit
  - Expected Rg for proteins: ~20–50 Å
  - **Finding**: All conditions show Rg ≈ 505–528 Å — a universal artifact

- **Orange points** (q < 0.1 Å⁻¹): Extended low-q fit
  - Shows Rg ≈ 36–39 Å with r² ≈ 0.65–0.70
  - Also universal, independent of sample

**Interpretation**: The ~500 Å feature is a shared instrumental/geometric artifact (confirmed by its appearance even in empty capillaries and blank buffers). It is NOT protein-specific and should not be interpreted as LLPS or aggregation.

### `guinier_acid_progression.png`
Guinier analysis (low-q region, q < 0.02 Å⁻¹) comparing Y2F samples with vs. without acid.

- Overlays curves for Y2F in different acid conditions
- The ~500 Å Rg artifact appears in all, independent of acid presence
- **Conclusion**: No LLPS (liquid-liquid phase separation) signature detected; acid does not induce large-scale aggregation in this dataset

### `drift/` subdirectory
Within-run stability checks (drift_within_runs.py output).

- Shows normalized intensity vs. scan order (chronological time)
- Reveals early/late RMSD and whether a condition is internally consistent
- Most conditions are stable after skipping the first cycle (settling transient)
- Two dropped conditions (`Y2F_kpi_kpi`, `resin`) showed >40% mid-run drift and were excluded

---

## What These Diagnostics Tell You

1. **Rg artifact is universal** — not confounded by sample condition
2. **No acid-induced large-scale aggregation** — if LLPS were occurring, we'd see q-dependent Rg shifts or two-population Guinier fits
3. **Most runs are internally consistent** — after removing the first cycle, RMSDs are small (~2–5%)
4. **Frame quality is high** — thanks to the r² streak-quality filter (r² > 0.999)

---

## Use in Presentations

When defending against the claim "your protein signal is actually an artifact," point to:
1. **rg_by_condition.png**: "The ~500 Å Rg is universal, even in blanks. It's instrumental, not protein-specific."
2. **drift/ figures**: "Each condition is internally stable. No signs of systematic drift that could masquerade as protein signal."


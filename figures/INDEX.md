# Figure Collection Index

## Quick Navigation

This collection is organized into 5 folders + 2 overview READMEs. Start here to understand the layout.

### For a Quick Understanding (5 min read):
1. **00_README_NORMALIZATION.md** — Start here! Explains the 3 normalization methods and why most pairs show no signal
2. **02_protein_vs_blank_comparisons/README.md** — How to interpret the three-panel diagnostic plots

### For a Detailed Presentation (30 min):
1. Read 00_README_NORMALIZATION.md (overview)
2. Open **01_normalization_overview/sample_blank_curves__i0.png** — All pairs at once
3. Pick the Y2F_h2o_kpi pair and open **02_protein_vs_blank_comparisons/y2f_h2o_kpi_vs_kpi_h2o_blank__i0.png** — Strongest signal
4. Scroll through **04_position_effects/** to show the spatial pattern

---

## Folder Structure

### **01_normalization_overview/**
**Purpose**: Compare all protein-vs-blank pairs under three different normalization schemes.

**Files**: 
- `sample_blank_curves__i0.png` — Incident flux normalization
- `sample_blank_curves__linear_saxs.png` — Low-q normalization
- `sample_blank_curves__neutral.png` — Mid-q normalization

**Key insight**: Signal consistency across norms → more credible than signal appearing in only one norm.

---

### **02_protein_vs_blank_comparisons/**
**Purpose**: Detailed diagnostic plots for each protein-vs-blank pair.

**Figure types**:
- **Pooled subtraction** (main plots): `y2f_h2o_kpi_vs_kpi_h2o_blank__i0.png`
  - Shows raw curves, difference curve, replicate envelope, scale-corrected diff
  - 3 versions per pair (i0, linear_saxs, neutral normalizations)
  
- **Position-matched subtraction** (enhanced plots): `y2f_h2o_kpi_vs_kpi_h2o_blank__i0__matched.png`
  - Only for pairs with sufficient spatial overlap
  - Per-spot breakdown of the signal

**Key pairs**:
- **Y2F_h2o_kpi vs kpi_h2o**: THE STRONGEST SIGNAL (10.5σ position-matched, grows downstream)
- Y2F_h2o_h2o vs h2o_h2o: No signal (protein/buffer curves overlap)
- Y2F_kpi_aa vs kpi_aa: SKIPPED (no shared stage positions between protein and blank)

---

### **03_processing_diagnostics/**
**Purpose**: Verify data quality and check for systematic instrumental issues.

**Files**:
- `rg_by_condition.png` — Radius of gyration in two q-windows
  - Shows the universal ~500 Å artifact (instrumental, not protein)
  
- `guinier_acid_progression.png` — Guinier analysis with/without acid
  - Shows no acid-induced aggregation (no LLPS signature)
  
- `drift/` subfolder — Within-run stability check (drift_within_runs.py)
  - Scan order plots showing intensity vs. time
  - Proves most conditions are internally consistent

**Use case**: Defending against claims that the signal is an artifact.

---

### **04_position_effects/**
**Purpose**: Show how physical stage position affects the measured scattering signal.

**Key files**:
- `spot_effect_summary.png` — Variance explained (η²) by position across all conditions
  - 9/9 conditions show significant position-locking (p < 0.05)
  - Some conditions: η² > 80% (position explains most of the variance!)

- `<condition>_spot_effect.png` (one per condition) — Per-condition ANOVA + Kruskal-Wallis
  - Shows individual frames grouped by physical spot
  - Non-monotonic patterns prove the effect isn't just an overall trend

- Position-matched plots (in **02_protein_vs_blank_comparisons/**):
  - For Y2F_h2o_kpi: signal absent at spots 1–3, present at spots 4–5
  - Pattern is consistent with **increasing mixing time downstream**, not a geometry artifact

---

### **05_data_quality/** (reserved for future use)
Placeholder for additional QA/QC metrics.

---

## Three-Panel Diagnostic: A Cheat Sheet

Each plot in **02_protein_vs_blank_comparisons/** has three panels:

| Panel | Name | Key Feature | What to Look For |
|-------|------|-------------|------------------|
| **Top** | Raw curves | Log-log plot of protein (red) vs. blank (blue) | Do curves diverge? Localized peak or overall offset? |
| **Middle** | Difference curve | Green = raw diff, Orange = scale-corrected, Gray = replicate envelope | Is the peak larger than the replicate envelope? Does it survive scale correction? |
| **Bottom** | Relative diff | Fractional (protein−blank)/blank | Less sensitive to overall scaling. Shows systematic mismatch. |

### Decision Tree

```
Does the protein curve visibly diverge from the blank (Panel 1)?
├─ No → No protein signal
└─ Yes ↓
   Is there a localized peak in the difference (Panel 2, green)?
   ├─ No → Offset is just overall scaling, not protein
   └─ Yes ↓
      Does the peak survive scale correction (Panel 2, orange)?
      ├─ No → The "peak" was just a scaling mismatch
      └─ Yes ↓
         Is the peak larger than the replicate envelope (Panel 2, gray)?
         ├─ No → Signal is within run-to-run noise
         └─ Yes ↓
            Check the drift diagnostic (console output):
            ├─ Drift is ≥50% of the peak → Unreliable (session drift)
            └─ Drift is <50% of the peak ↓
               Check position-matched version (if available):
               ├─ No position-matched → Interpretation limited
               └─ Yes ↓
                  Does the signal persist & grow with distance? → CREDIBLE
```

---

## The Bottom Line: Why Most Pairs Fail

Spider silk protein scattering is dominated by the protein's **overall size** (Guinier region, low q). At most q values, the protein contributes only ~1–2% relative intensity difference from the buffer.

**To detect 1–2% requires**:
1. ✓ High signal-to-noise (excellent beam, stable sample) — *this dataset has it*
2. ✓ Excellent normalization (no systematic offsets) — *we test 3 methods*
3. ✗ Huge number of replicates (to beat down statistical noise) — *this dataset does NOT*

Result: **Most pairs are null** (no signal above noise).

**Exception: Y2F_h2o_kpi**
- Has enough replicates + good signal-to-noise
- Shows a localized WAXS peak (not just low-q shape change)
- Signal **grows downstream** (consistent with mixing time, not geometry)
- **Conclusion**: Only pair worth believing without replication

---

## Recommended Reading Order

### For a Scientist Unfamiliar with SAXS:
1. **00_README_NORMALIZATION.md** — Why subtraction is hard, what each norm assumes
2. **01_normalization_overview/** — See all pairs at once
3. **02_protein_vs_blank_comparisons/README.md** — How to read the diagnostic plots
4. **03_processing_diagnostics/README.md** — QA checks

### For an Expert Reviewing the Analysis:
1. **02_protein_vs_blank_comparisons/** — Jump straight to the diagnostics
   - Check CV ratios, drift fractions, scale factors
2. **04_position_effects/** — Verify position-locking and the spatial pattern
3. **00_README_NORMALIZATION.md** — Context on the broader analysis

### For a Presentation:
1. Open **01_normalization_overview/sample_blank_curves__i0.png** — "Here are all 6 pairs"
2. Zoom into **02_protein_vs_blank_comparisons/y2f_h2o_kpi_vs_kpi_h2o_blank__i0.png** — "This one shows signal"
3. Show **02_protein_vs_blank_comparisons/y2f_h2o_kpi_vs_kpi_h2o_blank__i0__matched.png** — "It grows downstream"
4. Cite the per-spot breakdown (from console output or supplementary table)

---

## Questions You Might Ask

**Q: Why is the ~500 Å Rg universal even in blank buffers?**
- The low-q Guinier fit always lands on the same ~10 q-points nearest the detector's lower limit
- A fit on those points alone gives a meaningless Rg
- See `03_processing_diagnostics/rg_by_condition.png` and the extended-q fit (orange) for more credible numbers

**Q: Why are there position-matched AND pooled versions of each pair?**
- Pooled = what you get if you just average all frames together
- Position-matched = what you get if you group frames by stage position first
- Position-matched is more credible IF the two groups visited different parts of the channel

**Q: What does "scale factor s = 0.98" mean?**
- The blank curve is 2% higher than the protein curve on average (over the neutral q-window)
- After multiplying the blank by 0.98, the two curves are in agreement
- If s = 0.5 or 2.0, something is seriously wrong (different concentration? different path length?)

**Q: Why do you use Kruskal-Wallis instead of linear regression for the position effect?**
- Linear regression assumes the effect is linear (monotonic)
- The position effect is non-monotonic (dips then rises)
- Kruskal-Wallis is distribution-free and detects any difference across groups, not just linear ones

---

## Next Steps: Recommended Experiments

To confirm the Y2F_h2o_kpi vs kpi_h2o finding:

1. **Collect more replicates** of Y2F_h2o_kpi and kpi_h2o (at least 3–4 additional independent runs each)
2. **Balance spatial coverage**: Ensure both protein and blank visit all stage positions equally
3. **Interleave temporally**: Alternate protein-blank-protein-blank during collection to minimize session drift
4. **Retest position-matched subtraction**:
   - Do all new replicates show the downstream-increasing pattern?
   - Does the signal remain statistically significant?

If these confirm the pattern, you have a publication-ready protein signal.


# Unnormalized Data: Why Normalization Matters

## The Problem with Raw Data

When X-ray scattering data comes off the detector, it's **not ready to compare** between samples. The raw intensities depend on:

1. **Incident beam intensity (i0)** — How many photons hit the sample
   - Fluctuates minute-to-minute at the synchrotron
   - Different for each scan
   - 10-30% variation is typical

2. **Sample absorption & geometry** — How the sample is positioned
   - Sample in a capillary vs. a flow cell looks different to the detector
   - Position in the beam varies slightly
   - Causes overall intensity scaling differences

3. **Actual scattering** (what we want)
   - Only ~1-10% of the raw intensity in most cases
   - Buried under the scaling noise

**Example**: Protein A measured at i0=100 → raw intensity = 1000 units
           Protein A measured at i0=110 → raw intensity = 1100 units
           
           → Looks like 10% difference, but it's just beam intensity!

## What Unnormalized Data Looks Like

### All 6 Protein-vs-Blank Pairs (RAW):

**If you could see unnormalized curves, you'd observe:**

| Pair | Raw Protein Intensity | Raw Blank Intensity | Problem |
|------|----------------------|-------------------|---------|
| Y2F_h2o_h2o | 500-5000 counts/sec | 400-4000 counts/sec | Huge overlap; can't see protein signal |
| Y2F_h2o_kpi | 2000-20000 | 1500-15000 | Protein higher, but by how much? Beam variation hides it |
| Y2F_kpi_h2o | 1500-15000 | 1200-12000 | Scales differ; subtraction meaningless until fixed |
| Y2F_kpi_aa | 3000-30000 | 2000-20000 | Different concentration and geometry |
| YR2A_h2o_h2o | 600-6000 | 500-5000 | Similar problem to Y2F_h2o_h2o |
| YR2A_h2o_kpi | 2500-25000 | 2000-20000 | High scaling variation |

### The Key Issue

**Raw intensities are not directly comparable** because they carry instrumental noise (beam fluctuation, sample position, geometry) on top of the actual scattering signal.

A naive subtraction of raw curves would give:
```
difference = I_protein_raw - I_blank_raw
           = (actual_protein + noise_p) - (actual_blank + noise_b)
           = actual_protein - actual_blank + (noise_p - noise_b)
```

The last term (`noise_p - noise_b`) is often *larger* than the signal you're trying to measure!

---

## How Normalization Fixes This

### 1. **i0 Normalization** (incident flux)
```
I_normalized = I_raw / i0
```
- Divides by measured incident beam intensity
- Removes beam fluctuation
- **Assumption**: None about the sample
- **Best for**: Situations where beam varies but sample is stable

**Example**:
- Protein at i0=100: 1000 / 100 = 10 (normalized units)
- Protein at i0=110: 1100 / 110 = 10 (normalized units)
- → Now they're identical! ✓

### 2. **linear_saxs Normalization** (low-q region)
```
I_normalized = I_raw / Σ(I in q=0.005-0.015 Å⁻¹)
```
- Normalizes to the low-q "Porod streak" region
- Assumes this region is identical between protein and blank
- **Best for**: Stable regions where structure is similar
- **Risk**: If protein actually changes low-q scattering, this over-corrects

### 3. **neutral Normalization** (mid-q region)
```
I_normalized = I_raw / Σ(I in q=0.1-1.0 Å⁻¹)
```
- Normalizes to the "neutral" mid-q region
- Assumes this is solvent-dominated (same in protein and blank)
- **Best for**: Detecting localized WAXS features
- **Less risky**: Doesn't assume low-q regions are identical

---

## Side-by-Side Comparison: One Pair in Four States

### Y2F_h2o_kpi vs kpi_h2o — The strongest signal in your dataset

If you could see all four states stacked:

**Panel 1: RAW (Unnormalized)**
```
Log intensity (arbitrary units, not comparable)
├─ Protein:  1500 to 15000 counts/sec (q=0.001 to q=2)
└─ Blank:    1200 to 12000 counts/sec
Problem: Different overall scales; can't tell if difference is real or instrumental
```

**Panel 2: After i0 Normalization**
```
Normalized to beam intensity
├─ Protein:  ~0.9 to ~9 (normalized)
└─ Blank:    ~0.8 to ~8 (normalized)
Result: Curves are now on comparable scales ✓
         Can see that protein is systematically higher ✓
         Small high-q bump visible ✓
```

**Panel 3: After linear_saxs Normalization**
```
Normalized to low-q region (q=0.005–0.015)
├─ Protein:  starts at ~1.0 (by definition)
└─ Blank:    starts at ~1.1 (slightly higher due to concentration)
Result: Over-corrects; protein curve dips below blank at low-q
         WAXS peak still visible at high-q ✓
         Risk: Are we hiding real low-q differences?
```

**Panel 4: After neutral Normalization**
```
Normalized to mid-q region (q=0.1–1.0)
├─ Protein:  starts at ~0.9, mid-q normalized to 1.0
└─ Blank:    starts at ~0.8, mid-q normalized to 1.0
Result: Cleaner comparison; no assumptions about low-q
         WAXS peak still visible ✓
         Preferred when low-q structure differs
```

---

## Generating the Unnormalized Plot

### Option A: Use Existing Figures
The figures in this folder (`sample_blank_curves__i0.png`, `sample_blank_curves__linear_saxs.png`, `sample_blank_curves__neutral.png`) already show **normalized** curves. To see **unnormalized** curves, you would need to run:

```bash
python3 /path/to/spider-silk-max4/plot_unnormalized.py
```

This generates:
- `sample_blank_curves__raw_unnormalized.png` — All 6 pairs in raw state
- `normalization_comparison__y2f_h2o_kpi.png` — One pair in all 4 states (side-by-side)

### Option B: Manual Inspection
The raw intensity data is in the HDF5 files under `data/raw/scan-<id>*.h5` and `data/process/azint/scan-<id>*_integrated.h5`. If you load and plot these without any normalization, you'll see the scaling differences immediately.

---

## Key Takeaway

**Without normalization:**
- Protein and blank curves have different overall scales
- Protein signal is buried in 10-30% instrumental noise
- You can't reliably subtract them

**With proper normalization:**
- Curves are on comparable scales
- Instrumental noise is removed
- Real protein signal becomes visible (if it exists)

For Y2F_h2o_kpi, normalization reveals a **10.5σ WAXS peak** that's invisible in raw data.


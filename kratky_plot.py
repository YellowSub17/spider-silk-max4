"""
Kratky-plot analysis of the protein-vs-blank pairs from protein_vs_buffer.py.

A Kratky plot (I(q)*q^2 vs q) is a standard SAXS shape diagnostic:
  - a compact/globular particle gives a bell-shaped peak that decays back
    toward the baseline at high q,
  - an extended, intrinsically disordered, or aggregated/fractal particle
    instead plateaus or keeps rising at high q (no decay).

This is attractive here because it's a *global shape* statistic rather than
a single-peak search, so it isn't vulnerable to the baseline-anchoring
artifact that inflated the WAXS peak-significance numbers in
protein_vs_buffer.py. It also doesn't require picking one specific q-window
to trust the way the earlier peak/scale diagnostics did.

For each protein/blank pair, using the "neutral" normalization (the one that
gave the best-behaved, near-1.0 scale factor in protein_vs_buffer.py), this
script plots:
  1. Raw Kratky (I*q^2) for protein and blank overlaid, plus their
     background-subtracted difference -- does the protein contribute a
     compact-particle peak that decays, or does the signal keep rising?
  2. Dimensionless ("Receveur-Carrey") Kratky plot using each curve's own
     Guinier Rg/I0, which lets you compare shape independent of overall
     scale, and check against the canonical globular peak at
     (q*Rg, y) = (sqrt(3), ~1.1).

Run directly: `python kratky_plot.py`
Figures are written to ./figures/.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

import compare_groups as cg
import protein_vs_buffer as pvb
import src

NORM = "neutral"  # best-behaved normalization found in protein_vs_buffer.py


def high_q_plateau_ratio(qs, kratky_y, q_peak_search=(0.01, 0.3), q_plateau=(1.5, 2.0)):
    """
    Crude globular-vs-disordered indicator: ratio of the Kratky value at high
    q (where a globular particle should have decayed back down) to its peak
    value in the low/mid-q region. Ratio near/above 1 => plateau/rise (no
    decay, non-globular); ratio << 1 => decays as expected for a compact
    particle.
    """
    peak_mask = (qs >= q_peak_search[0]) & (qs <= q_peak_search[1])
    plateau_mask = (qs >= q_plateau[0]) & (qs <= q_plateau[1])
    if not np.any(peak_mask) or not np.any(plateau_mask):
        return np.nan
    peak_val = np.nanmax(kratky_y[peak_mask])
    plateau_val = np.nanmean(kratky_y[plateau_mask])
    if peak_val == 0 or not np.isfinite(peak_val):
        return np.nan
    return plateau_val / peak_val


def analyze_kratky(title, protein_name, blank_name, norm=NORM, out_dir="figures"):
    os.makedirs(out_dir, exist_ok=True)
    print(f"\nPair: {title}  [normalization: {norm}]")

    protein = cg.load_group(protein_name)
    blank = cg.load_group(blank_name)
    qs, Is_p, nk_p, nt_p = cg.kept_frames(protein, norm=norm)
    qs_b, Is_b, nk_b, nt_b = cg.kept_frames(blank, norm=norm)
    assert np.allclose(qs, qs_b)

    I_p, std_p = Is_p.mean(axis=0), Is_p.std(axis=0)
    I_b, std_b = Is_b.mean(axis=0), Is_b.std(axis=0)
    diff = I_p - I_b

    kratky_p = src.saxs_metrics.kratky(qs, I_p)
    kratky_b = src.saxs_metrics.kratky(qs, I_b)
    kratky_diff = src.saxs_metrics.kratky(qs, diff)

    ratio_p = high_q_plateau_ratio(qs, kratky_p)
    ratio_b = high_q_plateau_ratio(qs, kratky_b)
    print(f"  Raw Kratky high-q/peak ratio: protein={ratio_p:.2f}, blank={ratio_b:.2f}"
          + ("  (>~0.5 => plateauing/rising, not globular)" if np.isfinite(ratio_p) else ""))

    # Dimensionless Kratky, from each curve's own Guinier fit. Two windows:
    # the default (q<0.02, auto-refined to q*Rg<1.3 -- picks up the huge,
    # shared low-q aggregate/background feature) vs a fixed q<0.1 window
    # (no refinement) that reaches further out, past that low-q feature, to
    # see whether a smaller, protein-scale Guinier region shows up instead.
    g_p = src.saxs_metrics.guinier_fit(qs, I_p)
    g_b = src.saxs_metrics.guinier_fit(qs, I_b)
    g_p_wide = src.saxs_metrics.guinier_fit(qs, I_p, q_max=0.1, refine=False)
    g_b_wide = src.saxs_metrics.guinier_fit(qs, I_b, q_max=0.1, refine=False)
    for label, g, g_wide in [("protein", g_p, g_p_wide), ("blank", g_b, g_b_wide)]:
        if g is None:
            print(f"  Guinier fit (default, q<0.02) failed for {label}")
        else:
            print(f"  Guinier fit (default, q<0.02, {label}): "
                  f"Rg={g['Rg']:.1f} A, I0={g['I0']:.3g}, r={g['r']:.4f}")
        if g_wide is None:
            print(f"  Guinier fit (q<0.1, {label}) failed -- not linear over this window")
        else:
            print(f"  Guinier fit (q<0.1, {label}): "
                  f"Rg={g_wide['Rg']:.1f} A, I0={g_wide['I0']:.3g}, r={g_wide['r']:.4f}")

    fig, (ax_raw, ax_dim) = plt.subplots(1, 2, figsize=(13, 5))

    ax_raw.plot(qs, kratky_p, color="tab:red", label=f"{protein_name}")
    ax_raw.fill_between(qs, src.saxs_metrics.kratky(qs, I_p - std_p),
                         src.saxs_metrics.kratky(qs, I_p + std_p), color="tab:red", alpha=0.2, linewidth=0)
    ax_raw.plot(qs, kratky_b, color="tab:blue", label=f"{blank_name}")
    ax_raw.fill_between(qs, src.saxs_metrics.kratky(qs, I_b - std_b),
                         src.saxs_metrics.kratky(qs, I_b + std_b), color="tab:blue", alpha=0.2, linewidth=0)
    ax_raw.plot(qs, kratky_diff, color="tab:green", label="protein - blank")
    ax_raw.axhline(0, color="k", linewidth=0.8)
    ax_raw.set_xscale("log")
    ax_raw.set_xlabel("q [1/A]")
    ax_raw.set_ylabel("I(q) * q^2")
    ax_raw.set_title("Raw Kratky plot")
    ax_raw.legend(fontsize=8)

    for label, g, I_curve, color in [
        ("protein, q<0.02 fit", g_p, I_p, "tab:red"),
        ("blank, q<0.02 fit", g_b, I_b, "tab:blue"),
    ]:
        if g is None:
            continue
        x, y = src.saxs_metrics.dimensionless_kratky(qs, I_curve, g["Rg"], g["I0"])
        ax_dim.plot(x, y, color=color, linestyle="-", label=f"{label} (Rg={g['Rg']:.0f} A)")
    for label, g, I_curve, color in [
        ("protein, q<0.1 fit", g_p_wide, I_p, "tab:red"),
        ("blank, q<0.1 fit", g_b_wide, I_b, "tab:blue"),
    ]:
        if g is None:
            continue
        x, y = src.saxs_metrics.dimensionless_kratky(qs, I_curve, g["Rg"], g["I0"])
        ax_dim.plot(x, y, color=color, linestyle="--", label=f"{label} (Rg={g['Rg']:.0f} A)")
    ax_dim.axvline(np.sqrt(3), color="gray", linestyle=":", label="globular peak (qRg=sqrt(3))")
    ax_dim.axhline(1.1, color="gray", linestyle=":")
    ax_dim.set_xlim(0, 5)
    ax_dim.set_xlabel("q * Rg")
    ax_dim.set_ylabel("(q*Rg)^2 * I/I0")
    ax_dim.set_title("Dimensionless (Receveur-Carrey) Kratky plot")
    ax_dim.legend(fontsize=8)

    fig.suptitle(title)
    fig.tight_layout()
    out_path = os.path.join(out_dir, title.lower().replace(" ", "_") + "__kratky.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    for title, protein_name, blank_name in pvb.PAIRS:
        analyze_kratky(title, protein_name, blank_name)
    plt.show()

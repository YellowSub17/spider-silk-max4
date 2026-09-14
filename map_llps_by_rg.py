"""
Map apparent Guinier Rg across every condition in compare_groups.GROUPS, to
test the LLPS (liquid-liquid phase separation) hypothesis: Y2F/YR2A mixed
with acetic acid are expected to form droplets, which would show up as a
very large apparent Rg (~500 Å) in the low-q Guinier region.

Earlier analysis (kratky_plot.py) found Rg~500-530 Å in *both* a Y2F/KPi run
and its protein-free KPi blank -- i.e. in a condition with no acid at all --
which would argue against a clean LLPS signature and for a shared background/
aggregate feature instead. This script checks that systematically across
every condition rather than the handful spot-checked before.

For every group:
  1. Load + combine (streak-filtered via r_squared > 0.999, see
     compare_groups.mean_curve).
  2. Normalize with the "neutral" method (mid-q pinned; gave the best-behaved
     scale factors in protein_vs_buffer.py).
  3. Fit the Guinier region (q < 0.02 1/A, auto-refined to q*Rg < 1.3) via
     src.saxs_metrics.guinier_fit, extracting Rg, I0, r^2.

Outputs:
  - figures/rg_by_condition.png: Rg vs condition, sorted, colored by whether
    the condition name indicates acetic acid, annotated with fit r^2, with a
    reference line at Rg=500 A.
  - figures/guinier_acid_progression.png: 3-panel ln(I) vs q^2 Guinier plots
    with fit overlay for Y2F_buffer_pH8 -> Y2F_1M_acetic -> buffer_acetic_final.
  - A printed summary table, and a verdict on whether large Rg is
    acid-specific or widespread.

Run directly: `python map_llps_by_rg.py`
"""

import os

import matplotlib.pyplot as plt
import numpy as np

import compare_groups as cg
import src

NORM = "neutral"
RG_LLPS_LINE = 500.0   # reference line: the apparent Rg seen in earlier spot-checks
RG_FLAG = 200.0        # flag as "large" above this
R2_FLAG = 0.95         # flag as "poor fit" below this

# Y2F+KPi with no acid (baseline) -> Y2F+KPi+acid (expected LLPS) -> KPi+acid
# alone, no protein (control for the acid/buffer background itself).
# Note: Y2F_kpi_kpi (the closest true "same buffer, no acid" baseline) was
# dropped -- see compare_groups.DROPPED_CONDITIONS, it was internally
# inconsistent (108.8% early/late RMSD). Y2F_h2o_kpi is used instead as the
# no-acid reference; its channel1 is water rather than KPi, so this is no
# longer a single-variable (acid only) progression -- treat the first step
# as a rougher comparison than the other two.
ACID_PROGRESSION = ["Y2F_h2o_kpi", "Y2F_kpi_aa", "kpi_aa"]


def has_acid(name):
    """Uses the parsed Sample identity (compare_groups.SAMPLE_INFO), not a
    substring guess -- exact, since Sample strings are structured data."""
    return cg.SAMPLE_INFO[name]["has_acid"]


def fit_group(name, norm=NORM):
    """Load a group, normalize, streak-filter, and Guinier-fit the mean curve."""
    combo = cg.load_group(name)
    qs, I_mean, I_std, n_keep, n_total = cg.mean_curve(combo, norm=norm)
    g = src.saxs_metrics.guinier_fit(qs, I_mean, q_max=0.02)
    row = {
        "condition": name, "has_acid": has_acid(name),
        "n_kept": n_keep, "n_total": n_total,
        "qs": qs, "I_mean": I_mean,
    }
    if g is None:
        row.update(Rg=np.nan, I0=np.nan, r2=np.nan, q_fit_max=np.nan, n_points=0)
    else:
        row.update(Rg=g["Rg"], I0=g["I0"], r2=g["r"] ** 2,
                   q_fit_max=g["q_fit_max"], n_points=g["n_points"])
    return row


def build_table(norm=NORM):
    print(f"Fitting Guinier region for all {len(cg.GROUPS)} conditions "
          f"[normalization: {norm}]...")
    rows = []
    for name in cg.GROUPS:
        print(f"\n{name}")
        rows.append(fit_group(name, norm=norm))
    return rows


def print_summary(rows):
    print("\n" + "=" * 70)
    print(f"{'Condition':<24}{'Rg (A)':>10}{'I0':>12}{'r^2':>8}{'n_pts':>7}{'Has_Acid?':>11}")
    print("-" * 70)
    for r in sorted(rows, key=lambda r: (-r["Rg"] if np.isfinite(r["Rg"]) else np.inf)):
        rg_str = f"{r['Rg']:.1f}" if np.isfinite(r["Rg"]) else "fit failed"
        i0_str = f"{r['I0']:.3g}" if np.isfinite(r["I0"]) else "-"
        r2_str = f"{r['r2']:.4f}" if np.isfinite(r["r2"]) else "-"
        flags = []
        if np.isfinite(r["Rg"]) and r["Rg"] > RG_FLAG:
            flags.append("LARGE Rg")
        if np.isfinite(r["r2"]) and r["r2"] < R2_FLAG:
            flags.append("POOR FIT")
        if r["n_points"] and r["n_points"] < 15:
            flags.append(f"only {r['n_points']} pts -- fit may be edge-dominated")
        flag_str = ("  <-- " + ", ".join(flags)) if flags else ""
        print(f"{r['condition']:<24}{rg_str:>10}{i0_str:>12}{r2_str:>8}{r['n_points']:>7}"
              f"{str(r['has_acid']):>11}{flag_str}")
    print("=" * 70)


def plot_rg_by_condition(rows, out_dir="figures"):
    os.makedirs(out_dir, exist_ok=True)
    valid = [r for r in rows if np.isfinite(r["Rg"])]
    valid.sort(key=lambda r: r["Rg"])

    names = [r["condition"] for r in valid]
    rgs = [r["Rg"] for r in valid]
    r2s = [r["r2"] for r in valid]
    colors = ["tab:red" if r["has_acid"] else "tab:blue" for r in valid]

    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(names)), 6))
    x = np.arange(len(names))
    ax.bar(x, rgs, color=colors)
    for xi, rg, r2 in zip(x, rgs, r2s):
        ax.annotate(f"r2={r2:.2f}", (xi, rg), textcoords="offset points",
                    xytext=(0, 3), ha="center", fontsize=7, rotation=90, va="bottom")

    line_handle = ax.axhline(RG_LLPS_LINE, color="k", linestyle="--", linewidth=1,
                              label=f"Rg={RG_LLPS_LINE:.0f} A (LLPS-like)")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Apparent Rg [A]")
    ax.set_title(f"Guinier Rg by condition (normalization: {NORM})")

    from matplotlib.patches import Patch
    handles = [
        Patch(color="tab:red", label="acetic-acid-containing"),
        Patch(color="tab:blue", label="no acid"),
        line_handle,
    ]
    ax.legend(handles=handles, fontsize=8)

    fig.tight_layout()
    out_path = os.path.join(out_dir, "rg_by_condition.png")
    fig.savefig(out_path, dpi=150)
    print(f"\nSaved {out_path}")


def plot_guinier_progression(rows_by_name, out_dir="figures"):
    os.makedirs(out_dir, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

    for ax, name in zip(axes, ACID_PROGRESSION):
        r = rows_by_name[name]
        qs, I_mean = r["qs"], r["I_mean"]

        # Context window for the scatter (a bit wider than the fit itself)
        ctx_mask = (qs > 0) & (qs <= 0.05) & (I_mean > 0)
        x_ctx = qs[ctx_mask] ** 2
        y_ctx = np.log(I_mean[ctx_mask])
        ax.plot(x_ctx, y_ctx, ".", color="gray", markersize=3, label="data (q<0.05)")

        if np.isfinite(r["Rg"]):
            rg, i0 = r["Rg"], r["I0"]
            # Use the exact window guinier_fit actually used (it's a 2-stage
            # fit: an initial broad-window fit decides whether the narrower
            # q*Rg<1.3 refinement engages, so the final Rg does not by
            # itself tell you which window produced it).
            fit_mask = (qs > 0) & (qs <= r["q_fit_max"]) & (I_mean > 0)
            x_fit = qs[fit_mask] ** 2
            y_fit = np.log(I_mean[fit_mask])
            ax.plot(x_fit, y_fit, "o", color="tab:red", markersize=4,
                    label=f"fit window (n={r['n_points']})")

            x_line = np.array([0, x_fit.max() if x_fit.size else 0.02 ** 2])
            y_line = np.log(i0) - (rg ** 2 / 3) * x_line
            ax.plot(x_line, y_line, "k--", linewidth=1.5,
                    label=f"fit: Rg={rg:.1f} A, r2={r['r2']:.3f}")

        ax.set_xlabel("q^2 [1/A^2]")
        ax.set_ylabel("ln(I)")
        ax.set_title(name)
        ax.legend(fontsize=8)

    fig.suptitle("Guinier fits along the acetic-acid progression")
    fig.tight_layout()
    out_path = os.path.join(out_dir, "guinier_acid_progression.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


def print_progression_summary(rows_by_name):
    print("\nAcetic-acid progression detail:")
    print(f"  {'Condition':<22}{'Rg (A)':>10}{'I0':>12}{'r^2':>8}")
    for name in ACID_PROGRESSION:
        r = rows_by_name[name]
        rg_str = f"{r['Rg']:.1f}" if np.isfinite(r["Rg"]) else "fit failed"
        i0_str = f"{r['I0']:.3g}" if np.isfinite(r["I0"]) else "-"
        r2_str = f"{r['r2']:.4f}" if np.isfinite(r["r2"]) else "-"
        print(f"  {name:<22}{rg_str:>10}{i0_str:>12}{r2_str:>8}")


def print_verdict(rows):
    valid = [r for r in rows if np.isfinite(r["Rg"])]
    large = [r for r in valid if r["Rg"] > RG_FLAG]
    large_acid = [r for r in large if r["has_acid"]]
    large_no_acid = [r for r in large if not r["has_acid"]]

    print("\n" + "-" * 70)
    print("VERDICT: is the large apparent Rg an LLPS signature or a shared artifact?")
    print(f"  Conditions with Rg > {RG_FLAG:.0f} A: {len(large)}/{len(valid)}")
    print(f"    - acid-containing: {len(large_acid)} ({[r['condition'] for r in large_acid]})")
    print(f"    - no acid:         {len(large_no_acid)} ({[r['condition'] for r in large_no_acid]})")
    if large_no_acid and large_acid:
        print("  -> Large Rg appears in BOTH acid and non-acid conditions: "
              "NOT an acid-specific (LLPS) signature. Most likely a shared "
              "measurement artifact or background feature common to all runs.")
    elif large_acid and not large_no_acid:
        print("  -> Large Rg appears ONLY in acid-containing conditions: "
              "consistent with an LLPS-specific signature.")
    elif large_no_acid and not large_acid:
        print("  -> Large Rg appears only in non-acid conditions: "
              "inconsistent with LLPS being the cause.")
    else:
        print("  -> No condition shows a large apparent Rg.")
    print("-" * 70)


if __name__ == "__main__":
    rows = build_table()
    rows_by_name = {r["condition"]: r for r in rows}

    print_summary(rows)
    plot_rg_by_condition(rows)
    plot_guinier_progression(rows_by_name)
    print_progression_summary(rows_by_name)
    print_verdict(rows)

    plt.show()

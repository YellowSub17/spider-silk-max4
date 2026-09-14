"""
Look for a channel-position (i.e. mixing-time) dependent effect within each
flow-channel condition.

Previously this fit each metric against the "position N" LABEL's integer
index (1..5). That is wrong on two counts, both fixed here using
stage_positions.py's physical stage coordinates:

  1. The physical spacing between positions is NOT uniform: y-offsets are
     0, 0.5, 1.0, 3.0, 5.5 mm (from stage_positions.py) -- index treats the
     position-3-to-4 step (2.0 mm) as identical to position-1-to-2 (0.5 mm).
     An effect that's actually linear in physical distance/mixing-time looks
     nonlinear (and has much less statistical power) when fit against index.
     Fixed: regress against each spot's real y-offset in mm instead. Both
     fits are printed side by side so you can see how much this matters.
  2. The "position N" LABEL is not reliable in general (see
     stage_positions.py -- it can span multiple physical spots for messy
     conditions like Y2F_kpi_aa). Fixed: group frames by global stage SPOT
     ID (stage_positions.SPOT_OF), not by the raw label.

Also added:
  - A Kruskal-Wallis test across spots per metric (non-parametric, so a
    non-monotonic spot dependence -- which a linear distance fit would
    miss entirely -- still shows up).
  - Total integrated intensity as a 7th metric (drift_within_runs.py found
    it's the metric with the clearest position-locked cyclic structure, but
    it was never tested here before).

For each condition (every non-capillary group in compare_groups.GROUPS):
  1. Groups kept frames by physical spot.
  2. Plots I(q) per spot (viridis: dark = smallest y-offset/most upstream,
     light = largest).
  3. For each of 7 scalar metrics: per-spot mean+SEM, a distance-based
     linear fit, an index-based linear fit (for comparison), and a
     Kruskal-Wallis test across spots.

Run directly: `python position_series.py`
Figures are written to ./figures/.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import compare_groups as cg
import stage_positions as sp
import src
from drift_within_runs import total_intensity
from position_locking import spot_groups

NORM = "linear_saxs"  # matches this script's original default

# Every flow-channel (non-capillary) condition -- these used the position
# rotator as a mixing channel, worth checking for a position-dependent
# trend. Capillary (`cap:`) samples are a different, static measurement.
# analyze_group() itself skips any group with fewer than 2 usable spots.
CHANNEL_GROUPS = [name for name in cg.GROUPS if not cg.SAMPLE_INFO[name]["is_capillary"]]


def frame_metrics(qs, I):
    """All scalar metrics for one frame, including total intensity (not part
    of src.saxs_metrics.all_metrics, added here per this script's brief)."""
    m = src.saxs_metrics.all_metrics(qs, I)
    return {
        "Rg [A]": m["guinier"]["Rg"] if m["guinier"] else None,
        "Guinier I0": m["guinier"]["I0"] if m["guinier"] else None,
        "Porod exponent n": m["porod"]["n"] if m["porod"] else None,
        "Peak height": m["peak"]["height"] if m["peak"] else None,
        "Peak area": m["peak"]["area"] if m["peak"] else None,
        "Peak position [1/A]": m["peak"]["q_peak"] if m["peak"] else None,
        "Total intensity": total_intensity(qs, I),
    }


METRIC_LABELS = ["Rg [A]", "Guinier I0", "Porod exponent n", "Peak height",
                  "Peak area", "Peak position [1/A]", "Total intensity"]


def analyze_group(group_name, out_dir="figures"):
    combo = cg.load_group(group_name)  # default skip_first_cycle=True
    cg.NORM_METHODS[NORM](combo)
    keep = cg.kept_mask(combo)
    qs = combo.qs
    Is = combo.Is[keep]
    n_total = len(keep)
    n_keep = int(keep.sum())

    by_spot_rows = spot_groups(combo, keep, list(combo.scan_ids))
    spots = sorted(by_spot_rows)
    if len(spots) < 2:
        print(f"{group_name}: fewer than 2 distinct spots found, skipping")
        return

    # Real physical y-offset (mm) from this condition's most-upstream spot,
    # and an index (1..N) ranking spots in that same order -- for the
    # side-by-side index-vs-distance comparison.
    y_of_spot = {s: sp.SPOT_COORD[s][1] for s in spots}
    y0 = min(y_of_spot.values())
    offset_of_spot = {s: y_of_spot[s] - y0 for s in spots}
    spots_by_y = sorted(spots, key=lambda s: y_of_spot[s])
    index_of_spot = {s: i + 1 for i, s in enumerate(spots_by_y)}

    print(f"\n{group_name} (n_kept={n_keep}/{n_total})")
    for s in spots_by_y:
        print(f"  spot {s}: index={index_of_spot[s]}, y-offset={offset_of_spot[s]:.3f} mm, "
              f"n={len(by_spot_rows[s])}")

    os.makedirs(out_dir, exist_ok=True)
    cmap = plt.colormaps.get_cmap("viridis")
    fig_curve, ax_curve = plt.subplots(figsize=(6, 5))
    metric_values_by_spot = {s: {} for s in spots}  # spot -> {label: array}
    for i, s in enumerate(spots_by_y):
        rows = by_spot_rows[s]
        I_mean = Is[rows].mean(axis=0)
        color = cmap(i / max(len(spots_by_y) - 1, 1))
        ax_curve.plot(qs, I_mean, color=color,
                      label=f"spot {s} (y+{offset_of_spot[s]:.1f}mm, n={len(rows)})")

        frame_dicts = [frame_metrics(qs, I) for I in Is[rows]]
        for label in METRIC_LABELS:
            vals = np.array([d[label] for d in frame_dicts if d[label] is not None], dtype=float)
            metric_values_by_spot[s][label] = vals

    ax_curve.set_xscale("log")
    ax_curve.set_yscale("log")
    ax_curve.set_xlim(src.det.FULL_QRANGE)
    ax_curve.set_xlabel("q [1/A]")
    ax_curve.set_ylabel("Inten. (norm.)")
    ax_curve.set_title(f"{group_name}: I(q) by physical spot")
    ax_curve.legend(fontsize=8)
    fig_curve.tight_layout()
    curve_path = os.path.join(out_dir, f"{group_name}_position_curves.png")
    fig_curve.savefig(curve_path, dpi=150)
    print(f"Saved {curve_path}")

    # --- Metric grid: one panel per scalar metric, vs y-offset (mm) --------
    fig, axes = plt.subplots(2, 4, figsize=(19, 8))
    for ax, label in zip(axes.flat, METRIC_LABELS):
        offsets = np.array([offset_of_spot[s] for s in spots_by_y])
        indices = np.array([index_of_spot[s] for s in spots_by_y])
        means = np.array([metric_values_by_spot[s][label].mean()
                          if metric_values_by_spot[s][label].size else np.nan for s in spots_by_y])
        sems = np.array([
            metric_values_by_spot[s][label].std(ddof=1) / np.sqrt(metric_values_by_spot[s][label].size)
            if metric_values_by_spot[s][label].size > 1 else 0.0
            for s in spots_by_y
        ])
        valid = np.isfinite(means)

        ax.errorbar(offsets[valid], means[valid], yerr=sems[valid], marker="o", color="k")

        dist_line = ""
        if valid.sum() >= 2:
            slope_d, intercept_d, r_d, p_d, se_d = stats.linregress(offsets[valid], means[valid])
            xs = np.array([offsets.min(), offsets.max()])
            ax.plot(xs, intercept_d + slope_d * xs, "r--", label=f"vs distance: r={r_d:.2f}, p={p_d:.3g}")
            dist_line = f"slope={slope_d:.4g}, r={r_d:.3f}, p={p_d:.4g}"

            slope_i, intercept_i, r_i, p_i, se_i = stats.linregress(indices[valid], means[valid])
            idx_line = f"slope={slope_i:.4g}, r={r_i:.3f}, p={p_i:.4g}"
        else:
            idx_line = "not enough valid spots"
            dist_line = "not enough valid spots"

        # Kruskal-Wallis across spots' raw per-frame values (not just means)
        groups = [metric_values_by_spot[s][label][np.isfinite(metric_values_by_spot[s][label])]
                  for s in spots_by_y]
        groups = [g for g in groups if len(g) >= 1]
        kw_p = np.nan
        if len(groups) >= 2:
            try:
                kw_p = stats.kruskal(*groups).pvalue
            except Exception:
                pass
        kw_str = f"{kw_p:.3g}" if np.isfinite(kw_p) else "n/a"

        d_sig = "SIGNIFICANT" if (valid.sum() >= 2 and p_d < 0.05) else "not significant"
        i_sig = "SIGNIFICANT" if (valid.sum() >= 2 and p_i < 0.05) else "not significant"
        kw_sig = "SIGNIFICANT" if np.isfinite(kw_p) and kw_p < 0.05 else "not significant"
        print(f"  {label}:")
        print(f"    vs distance (mm): {dist_line} -> {d_sig}")
        print(f"    vs index (1..N):  {idx_line} -> {i_sig}")
        print(f"    Kruskal-Wallis across spots: p={kw_str} -> {kw_sig}")

        ax.legend(fontsize=7)
        ax.set_xlabel("y-offset from upstream spot [mm]")
        ax.set_ylabel(label)

    for ax in axes.flat[len(METRIC_LABELS):]:
        ax.axis("off")

    fig.suptitle(f"{group_name}: SAXS/WAXS metrics vs physical position (mm)")
    fig.tight_layout()
    metrics_path = os.path.join(out_dir, f"{group_name}_position_metrics.png")
    fig.savefig(metrics_path, dpi=150)
    print(f"Saved {metrics_path}")


if __name__ == "__main__":
    for g in CHANNEL_GROUPS:
        analyze_group(g)
    plt.show()

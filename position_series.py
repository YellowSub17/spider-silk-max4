"""
Look for a channel-position (i.e. mixing-time) dependent effect within each
condition. The CSV's Description column records the rotator/channel position
("position 1".."position 5") each scan was taken at; for the flow-mixing runs
that position corresponds to distance along the mixing channel, i.e. time
since the protein and buffer met.

For each condition group, this script:
  1. Splits scans by position.
  2. Combines repeats at each position and computes a mean I(q) curve.
  3. Overlays the position curves (viridis: dark = position 1, light = last).
  4. Computes a simple beta-sheet peak / baseline intensity ratio per position
     and fits a line against position number to check for a monotonic trend.

Run directly: `python position_series.py`
Figures are written to ./figures/.
"""

import csv
import os
import re

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import compare_groups as cg
import src

# --- Map ScanID -> channel position, from the "position N" description label ---
# (sourced from scans_clean.csv, the same file compare_groups.py builds GROUPS
# from, so position and sample-identity labels are guaranteed consistent)
def load_scan_positions():
    positions = {}
    with open(cg.SCANS_CLEAN_CSV, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            sid = row["scanID"].strip()
            if not sid:
                continue
            m = re.fullmatch(r"position (\d+)", row["description"].strip())
            if m:
                positions[int(sid)] = int(m.group(1))
    return positions


POSITIONS = load_scan_positions()

# Every flow-channel (non-capillary) condition group -- these used the
# position rotator as a mixing channel, i.e. worth checking for a
# position-dependent trend. Capillary (`cap:`) samples are a different, static
# measurement and excluded. analyze_group() itself skips any group that
# doesn't actually have >=2 distinct positions recorded.
CHANNEL_GROUPS = [
    name for name in cg.GROUPS
    if not cg.SAMPLE_INFO[name]["is_capillary"] and name != "resin"
]

def split_by_position(group_name):
    """Return {position: [scan_ids]} for scans in this group carrying a position label."""
    by_pos = {}
    for sid in cg.GROUPS[group_name]:
        pos = POSITIONS.get(sid)
        if pos is not None:
            by_pos.setdefault(pos, []).append(sid)
    return by_pos


# Scalar metrics extracted per-frame via src.saxs_metrics, and how to pull
# each one out of the dict that all_metrics() returns.
METRIC_SPECS = [
    ("Rg [A]",            lambda m: m["guinier"]["Rg"] if m["guinier"] else None),
    ("Guinier I0",        lambda m: m["guinier"]["I0"] if m["guinier"] else None),
    ("Porod exponent n",  lambda m: m["porod"]["n"] if m["porod"] else None),
    ("Peak height",       lambda m: m["peak"]["height"] if m["peak"] else None),
    ("Peak area",         lambda m: m["peak"]["area"] if m["peak"] else None),
    ("Peak position [1/A]", lambda m: m["peak"]["q_peak"] if m["peak"] else None),
]


def load_position(scan_ids, load_n=1):
    """Load+combine scans at one position, skipping any that fail (see compare_groups)."""
    good_ids = []
    for sid in scan_ids:
        try:
            src.Scan(sid, load_imgs=False, load_n=load_n)
        except Exception:
            continue
        good_ids.append(sid)
    return src.ComboScan(good_ids, load_imgs=False, load_n=load_n)


def analyze_group(group_name, out_dir="figures"):
    by_pos = split_by_position(group_name)
    positions = sorted(by_pos)
    if len(positions) < 2:
        print(f"{group_name}: fewer than 2 distinct positions found, skipping")
        return

    os.makedirs(out_dir, exist_ok=True)
    cmap = plt.colormaps.get_cmap("viridis")

    # metric_values[label] = list (per position) of arrays of per-frame values
    metric_values = {label: [] for label, _ in METRIC_SPECS}

    fig_curve, ax_curve = plt.subplots(figsize=(6, 5))
    print(f"\n{group_name}")
    for i, pos in enumerate(positions):
        scan_ids = by_pos[pos]
        combo = load_position(scan_ids)
        qs, I_mean, I_std, n_keep, n_total = cg.mean_curve(combo)
        print(f"  position {pos}: {n_keep}/{n_total} scans kept")

        color = cmap(i / max(len(positions) - 1, 1))
        ax_curve.plot(qs, I_mean, color=color, label=f"pos {pos} (n={n_keep})")

        keep = combo.r_squared > cg.RSQ_LIM
        if not np.any(keep):
            keep = np.ones_like(combo.r_squared, dtype=bool)

        frame_metrics = [src.saxs_metrics.all_metrics(qs, I) for I in combo.Is[keep]]
        for label, extract in METRIC_SPECS:
            vals = np.array([v for m in frame_metrics if (v := extract(m)) is not None])
            metric_values[label].append(vals)

    ax_curve.set_xscale("log")
    ax_curve.set_yscale("log")
    ax_curve.set_xlim(src.det.FULL_QRANGE)
    ax_curve.set_xlabel("q [1/A]")
    ax_curve.set_ylabel("Inten. (norm.)")
    ax_curve.set_title(f"{group_name}: I(q) by channel position")
    ax_curve.legend()
    fig_curve.tight_layout()
    curve_path = os.path.join(out_dir, f"{group_name}_position_curves.png")
    fig_curve.savefig(curve_path, dpi=150)
    print(f"Saved {curve_path}")

    # --- Metric grid: one panel per scalar SAXS metric, vs channel position ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, (label, _) in zip(axes.flat, METRIC_SPECS):
        vals_per_pos = metric_values[label]
        means = np.array([v.mean() if len(v) else np.nan for v in vals_per_pos])
        sems = np.array([v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
                          for v in vals_per_pos])
        valid = np.isfinite(means)

        ax.errorbar(np.array(positions)[valid], means[valid], yerr=sems[valid],
                     marker="o", color="k")

        if valid.sum() >= 2:
            slope, intercept, r, p, se = stats.linregress(
                np.array(positions)[valid], means[valid])
            xs = np.array([min(positions), max(positions)])
            ax.plot(xs, intercept + slope * xs, "r--",
                     label=f"r={r:.2f}, p={p:.3g}")
            ax.legend(fontsize=8)
            verdict = "SIGNIFICANT (p<0.05)" if p < 0.05 else "not significant"
            print(f"  {label}: slope={slope:.4g}, r={r:.3f}, p={p:.4g} -> {verdict}")
        else:
            print(f"  {label}: not enough valid positions to fit a trend")

        ax.set_xlabel("channel position")
        ax.set_ylabel(label)

    fig.suptitle(f"{group_name}: SAXS/WAXS metrics vs channel position")
    fig.tight_layout()
    metrics_path = os.path.join(out_dir, f"{group_name}_position_metrics.png")
    fig.savefig(metrics_path, dpi=150)
    print(f"Saved {metrics_path}")


if __name__ == "__main__":
    for g in CHANNEL_GROUPS:
        analyze_group(g)
    plt.show()

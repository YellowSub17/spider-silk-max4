"""
How self-consistent is each condition's own run, with NO background
subtraction and regardless of whether it contains protein?

This is a prerequisite check for everything else in this analysis: if a
single condition's repeated exposures already disagree with each other more
than any protein-vs-blank difference we're looking for, there's no point
chasing that difference further -- it's below the noise floor set by the
measurement itself. (This was motivated directly by protein_vs_buffer.py's
drift_check: comparing kpi_h2o's own early half to its own late half
reproduced 60-440% of several claimed "protein signatures".)

For every condition in compare_groups.GROUPS, this script:
  1. Loads every kept (streak-filtered) frame, i0-normalized only (flux
     correction per frame, no q-region pinned to match anything -- the
     least assumption-laden normalization available, so it doesn't hide or
     manufacture drift).
  2. Plots I(q) for every individual frame, colored by chronological scan
     order (dark viridis = earliest, light = latest) -- drift is visible
     directly as a color gradient in the curves.
  3. Tracks two scalar diagnostics per frame vs scan order: total
     integrated intensity, and WAXS peak height (both on the RAW curve,
     not a difference -- so src.saxs_metrics.peak_analysis's default
     raw-curve baseline windows are the correct ones here, unlike in
     protein_vs_buffer.py).
  4. Computes one summary number per condition: the relative RMSD between
     the mean curve of the first half of its scans and the mean curve of
     the second half -- "how much does this condition's own curve change
     from start to end of its run", independent of any other condition.

Run directly: `python drift_within_runs.py`
Figures are written to ./figures/drift/.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np

import compare_groups as cg
import src

OUT_DIR = "figures/drift"
NORM = "i0"  # flux-only normalization -- doesn't pin any q-region, so real
             # drift (in shape, not just scale) isn't masked or manufactured
MIN_FOR_SPLIT = 6  # need at least this many kept frames for an early/late split

RMSD_GOOD = 0.10     # <10% relative RMSD: run looks internally consistent
RMSD_MARGINAL = 0.25  # 10-25%: marginal; >25%: chaotic, treat with real suspicion


def load_scan_times():
    times = {}
    with open(cg.SCANS_CLEAN_CSV, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            sid = row["scanID"].strip()
            if sid and row["start time"].strip():
                times[int(sid)] = row["start time"].strip()
    return times


SCAN_TIMES = load_scan_times()


def elapsed_minutes(name):
    ids = sorted(cg.GROUPS[name])
    t0, t1 = SCAN_TIMES.get(ids[0]), SCAN_TIMES.get(ids[-1])
    if not (t0 and t1):
        return np.nan
    import datetime
    fmt = "%Y-%m-%dT%H:%M:%S.%f"
    dt = datetime.datetime.strptime(t1, fmt) - datetime.datetime.strptime(t0, fmt)
    return dt.total_seconds() / 60


def total_intensity(qs, I):
    return np.trapezoid(I, qs)


def drift_metrics(name, norm=NORM):
    """Load every kept frame for a condition (in chronological order) and
    compute per-frame + summary drift diagnostics."""
    combo = cg.load_group(name)
    qs, Is, n_kept, n_total = cg.kept_frames(combo, norm=norm)
    if n_kept == 0:
        return None

    totals = np.array([total_intensity(qs, I) for I in Is])
    peaks = np.array([
        (pk["height"] if (pk := src.saxs_metrics.peak_analysis(qs, I)) else np.nan)
        for I in Is
    ])

    result = {
        "name": name, "qs": qs, "Is": Is, "n_kept": n_kept, "n_total": n_total,
        "totals": totals, "peaks": peaks,
        "elapsed_min": elapsed_minutes(name),
        "rel_rmsd": np.nan, "total_intensity_drift": np.nan,
    }

    if n_kept >= MIN_FOR_SPLIT:
        half = n_kept // 2
        I_early = Is[:half].mean(axis=0)
        I_late = Is[half:].mean(axis=0)
        mask = I_early > 0
        result["rel_rmsd"] = np.sqrt(np.mean(((I_late[mask] - I_early[mask]) / I_early[mask]) ** 2))
        result["total_intensity_drift"] = (totals[half:].mean() - totals[:half].mean()) / totals[:half].mean()

    return result


def plot_condition(result, out_dir=OUT_DIR):
    os.makedirs(out_dir, exist_ok=True)
    name, qs, Is, n = result["name"], result["qs"], result["Is"], result["n_kept"]

    fig, (ax_curves, ax_total, ax_peak) = plt.subplots(1, 3, figsize=(16, 5))
    cmap = plt.colormaps.get_cmap("viridis")
    for i, I in enumerate(Is):
        ax_curves.plot(qs, I, color=cmap(i / max(n - 1, 1)), alpha=0.6, linewidth=0.8)
    ax_curves.set_xscale("log")
    ax_curves.set_yscale("log")
    ax_curves.set_xlabel("q [1/A]")
    ax_curves.set_ylabel("Inten. (i0-norm.)")
    ax_curves.set_title("I(q), colored by scan order\n(dark=early, light=late)")

    idx = np.arange(n)
    ax_total.plot(idx, result["totals"], "o-", color="tab:blue", markersize=3)
    ax_total.set_xlabel("scan order (chronological)")
    ax_total.set_ylabel("total intensity (trapz over q)")
    ax_total.set_title("Total intensity vs time")

    ax_peak.plot(idx, result["peaks"], "o-", color="tab:red", markersize=3)
    ax_peak.set_xlabel("scan order (chronological)")
    ax_peak.set_ylabel("WAXS peak height")
    ax_peak.set_title("WAXS peak height vs time")

    rmsd_str = f"{result['rel_rmsd']:.1%}" if np.isfinite(result["rel_rmsd"]) else "n/a (too few scans to split)"
    span_str = f"{result['elapsed_min']:.0f} min" if np.isfinite(result["elapsed_min"]) else "unknown span"
    fig.suptitle(f"{name} (n={n} kept/{result['n_total']} total, {span_str}) "
                 f"-- early-vs-late curve RMSD: {rmsd_str}")
    fig.tight_layout()
    out_path = os.path.join(out_dir, f"{name}_drift.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


def summary_plot(results, out_dir=OUT_DIR):
    os.makedirs(out_dir, exist_ok=True)
    valid = [r for r in results if np.isfinite(r["rel_rmsd"])]
    valid.sort(key=lambda r: r["rel_rmsd"])
    names = [r["name"] for r in valid]
    vals = [r["rel_rmsd"] * 100 for r in valid]
    colors = ["tab:green" if v < RMSD_GOOD * 100 else
              ("tab:orange" if v < RMSD_MARGINAL * 100 else "tab:red") for v in vals]

    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(names)), 6))
    ax.bar(range(len(names)), vals, color=colors)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("early-vs-late relative RMSD [%]")
    ax.axhline(RMSD_GOOD * 100, color="gray", linestyle="--", linewidth=1,
               label=f"{RMSD_GOOD:.0%} (self-consistent)")
    ax.axhline(RMSD_MARGINAL * 100, color="gray", linestyle=":", linewidth=1,
               label=f"{RMSD_MARGINAL:.0%} (marginal / chaotic above this)")
    ax.set_title("Within-run drift by condition (green=stable, red=chaotic)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out_path = os.path.join(out_dir, "drift_summary.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


def print_summary(results):
    print(f"\n{'Condition':<24}{'n_kept':>7}{'span':>10}{'RMSD':>8}{'tot.drift':>10}{'CV':>8}")
    print("-" * 70)
    ordered = sorted(results, key=lambda r: r["rel_rmsd"] if np.isfinite(r["rel_rmsd"]) else -1,
                      reverse=True)
    for r in ordered:
        rmsd_str = f"{r['rel_rmsd']:.1%}" if np.isfinite(r["rel_rmsd"]) else "n/a"
        totd_str = f"{r['total_intensity_drift']:+.1%}" if np.isfinite(r["total_intensity_drift"]) else "n/a"
        span_str = f"{r['elapsed_min']:.0f}min" if np.isfinite(r["elapsed_min"]) else "?"
        cv = cg.coeff_of_variation(r["qs"], r["Is"].mean(axis=0), r["Is"].std(axis=0))
        flag = ""
        if np.isfinite(r["rel_rmsd"]):
            if r["rel_rmsd"] >= RMSD_MARGINAL:
                flag = "  <-- CHAOTIC"
            elif r["rel_rmsd"] >= RMSD_GOOD:
                flag = "  <-- marginal"
        print(f"{r['name']:<24}{r['n_kept']:>7}{span_str:>10}{rmsd_str:>8}{totd_str:>10}{cv:>7.1%}{flag}")

    n_chaotic = sum(1 for r in results if np.isfinite(r["rel_rmsd"]) and r["rel_rmsd"] >= RMSD_MARGINAL)
    n_split = sum(1 for r in results if np.isfinite(r["rel_rmsd"]))
    print(f"\n{n_chaotic}/{n_split} conditions with enough scans to check are CHAOTIC "
          f"(early-vs-late RMSD >= {RMSD_MARGINAL:.0%}).")


if __name__ == "__main__":
    results = []
    for name in cg.GROUPS:
        print(f"\n{name}")
        res = drift_metrics(name)
        if res is None:
            print("  no usable frames, skipping")
            continue
        results.append(res)
        plot_condition(res)

    summary_plot(results)
    print_summary(results)
    plt.show()

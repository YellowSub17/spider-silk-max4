"""
Quantify whether the cyclic structure seen in drift_within_runs.py's
"total intensity vs scan order" plots is truly POSITION-locked (tied to a
physical stage spot, i.e. spatial) rather than merely periodic with the
condition's cycle length (which would look identical in a scan-order plot
but wouldn't prove a spatial cause).

For every condition in compare_groups.GROUPS, i0-normalized (same
normalization and the same total_intensity/WAXS-peak-height metrics as
drift_within_runs.py, reused directly from there):
  1. Groups per-frame total intensity and WAXS peak height by global stage
     spot id (stage_positions.SPOT_OF).
  2. Reports per-spot n/mean/std, and peak-to-peak spread across spots as a
     percentage of the grand mean.
  3. Runs a one-way ANOVA and a Kruskal-Wallis test across spots (the latter
     doesn't assume normality/equal variance, appropriate given how skewed
     some of these per-frame metrics are), plus the fraction of total
     variance explained by spot (eta-squared: between-spot SS / total SS).
  4. Does all of the above twice: once on every frame, once excluding each
     condition's first full position cycle (compare_groups.SKIP_FIRST_CYCLE)
     -- so the settling transient's own contribution to any spot effect is
     visible by comparing the two.

Figures:
  figures/position/<condition>_spot_effect.png -- total intensity and WAXS
    peak height vs spot id (strip plot, colored by whether that frame was in
    the excluded first cycle), one per condition.
  figures/position/spot_effect_summary.png -- conditions ranked by variance
    explained by spot (with-first-cycle vs without, side by side).

Run directly: `python position_locking.py`
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import compare_groups as cg
import stage_positions as sp
import src
from drift_within_runs import total_intensity

OUT_DIR = "figures/position"
NORM = "i0"
MIN_FRAMES_PER_SPOT = 2   # spots with fewer kept frames than this are dropped
                          # from the stats (can't estimate their own spread)
MIN_SPOTS_FOR_TEST = 3    # need at least this many spots with enough frames


def per_frame_metrics(qs, Is):
    """Total integrated intensity and raw WAXS peak height per frame."""
    totals = np.array([total_intensity(qs, I) for I in Is])
    peaks = np.array([
        (pk["height"] if (pk := src.saxs_metrics.peak_analysis(qs, I)) else np.nan)
        for I in Is
    ])
    return totals, peaks


def spot_groups(combo, keep_mask, scan_ids):
    """{spot_id: row-indices into the kept-frame arrays} for the frames in
    this combo that passed keep_mask, using each frame's originating scan ID
    to look up its global stage spot (stage_positions.SPOT_OF)."""
    kept_ids = [sid for sid, k in zip(scan_ids, keep_mask) if k]
    by_spot = {}
    for row, sid in enumerate(kept_ids):
        spot = sp.SPOT_OF.get(sid)
        if spot is not None:
            by_spot.setdefault(spot, []).append(row)
    return by_spot


def variance_explained(values_by_group):
    """Eta-squared: between-group SS / total SS, for a dict of group ->
    array of values. Returns nan if fewer than 2 usable groups."""
    groups = [v for v in values_by_group.values() if len(v) >= 1]
    if len(groups) < 2:
        return np.nan
    all_vals = np.concatenate(groups)
    grand_mean = np.mean(all_vals)
    ss_total = np.sum((all_vals - grand_mean) ** 2)
    if ss_total == 0:
        return np.nan
    ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups)
    return ss_between / ss_total


def run_tests(values_by_spot):
    """One-way ANOVA + Kruskal-Wallis across spots with >=MIN_FRAMES_PER_SPOT
    values each. Returns dict with n_spots_used, anova_p, kw_p, eta2,
    pk_pk_pct (peak-to-peak spread across spot means, as % of grand mean)."""
    usable = {s: np.asarray(v)[np.isfinite(v)] for s, v in values_by_spot.items()}
    usable = {s: v for s, v in usable.items() if len(v) >= MIN_FRAMES_PER_SPOT}

    result = {"n_spots_used": len(usable), "anova_p": np.nan, "kw_p": np.nan,
              "eta2": np.nan, "pk_pk_pct": np.nan}
    if len(usable) < MIN_SPOTS_FOR_TEST:
        return result

    arrays = list(usable.values())
    means = [np.mean(a) for a in arrays]
    grand_mean = np.mean(np.concatenate(arrays))
    if grand_mean != 0:
        result["pk_pk_pct"] = (max(means) - min(means)) / abs(grand_mean) * 100

    try:
        result["anova_p"] = stats.f_oneway(*arrays).pvalue
    except Exception:
        pass
    try:
        result["kw_p"] = stats.kruskal(*arrays).pvalue
    except Exception:
        pass
    result["eta2"] = variance_explained(usable)
    return result


def analyze_condition(name):
    """Returns dict with per-spot metrics and test results, both including
    and excluding the first cycle, plus the raw per-frame arrays (for
    plotting)."""
    combo_all = cg.load_group(name, skip_first_cycle=False)
    scan_ids_all = list(combo_all.scan_ids)

    cg.NORM_METHODS[NORM](combo_all)
    keep_all = cg.kept_mask(combo_all)
    Is_all = combo_all.Is[keep_all]
    totals_all, peaks_all = per_frame_metrics(combo_all.qs, Is_all)
    by_spot_all_totals, by_spot_all_peaks = {}, {}
    spot_rows_all = spot_groups(combo_all, keep_all, scan_ids_all)
    for spot, rows in spot_rows_all.items():
        by_spot_all_totals[spot] = totals_all[rows]
        by_spot_all_peaks[spot] = peaks_all[rows]

    n_dropped = len(cg.first_cycle_ids(name))
    if n_dropped:
        # Excluding-first-cycle pass: same kept mask, but with the first
        # n_dropped chronological frames additionally excluded.
        keep_excl = cg.kept_mask(combo_all, skip_first_n=n_dropped)
        Is_excl = combo_all.Is[keep_excl]
        totals_excl, peaks_excl = per_frame_metrics(combo_all.qs, Is_excl)
        spot_rows_excl = spot_groups(combo_all, keep_excl, scan_ids_all)
        by_spot_excl_totals = {s: totals_excl[r] for s, r in spot_rows_excl.items()}
        by_spot_excl_peaks = {s: peaks_excl[r] for s, r in spot_rows_excl.items()}
    else:
        keep_excl = keep_all
        totals_excl, peaks_excl = totals_all, peaks_all
        by_spot_excl_totals, by_spot_excl_peaks = by_spot_all_totals, by_spot_all_peaks

    return {
        "name": name, "qs": combo_all.qs,
        "scan_ids": scan_ids_all, "keep_all": keep_all, "keep_excl": keep_excl,
        "totals_all": totals_all, "peaks_all": peaks_all,
        "totals_excl": totals_excl, "peaks_excl": peaks_excl,
        "by_spot_all_totals": by_spot_all_totals, "by_spot_all_peaks": by_spot_all_peaks,
        "by_spot_excl_totals": by_spot_excl_totals, "by_spot_excl_peaks": by_spot_excl_peaks,
        "n_dropped": n_dropped,
        "tests_all_totals": run_tests(by_spot_all_totals),
        "tests_all_peaks": run_tests(by_spot_all_peaks),
        "tests_excl_totals": run_tests(by_spot_excl_totals),
        "tests_excl_peaks": run_tests(by_spot_excl_peaks),
    }


def print_result(r):
    name = r["name"]
    print(f"\n{name} (n_kept={r['keep_all'].sum()}, first-cycle drop={r['n_dropped']})")
    print(f"  {'spot':>6}{'n':>5}{'mean total':>13}{'std total':>12}"
          f"{'mean peak':>12}{'std peak':>11}")
    for spot in sorted(r["by_spot_all_totals"]):
        t = r["by_spot_all_totals"][spot]
        p = r["by_spot_all_peaks"][spot]
        x, y = sp.SPOT_COORD.get(spot, (np.nan, np.nan))
        print(f"  {spot:>6}{len(t):>5}{np.mean(t):>13.4g}{np.std(t):>12.4g}"
              f"{np.nanmean(p):>12.4g}{np.nanstd(p):>11.4g}   (x={x:.3f}, y={y:.3f})")

    def fmt_test(label, t):
        anova = f"{t['anova_p']:.3g}" if np.isfinite(t["anova_p"]) else "n/a"
        kw = f"{t['kw_p']:.3g}" if np.isfinite(t["kw_p"]) else "n/a"
        eta2 = f"{t['eta2']:.1%}" if np.isfinite(t["eta2"]) else "n/a"
        pkpk = f"{t['pk_pk_pct']:.1f}%" if np.isfinite(t["pk_pk_pct"]) else "n/a"
        sig = ""
        if np.isfinite(t["kw_p"]):
            sig = "  <-- POSITION-LOCKED (p<0.05)" if t["kw_p"] < 0.05 else ""
        print(f"  {label:<32}n_spots={t['n_spots_used']:<4}ANOVA p={anova:<9}"
              f"KW p={kw:<9}eta^2={eta2:<8}pk-pk={pkpk}{sig}")

    fmt_test("total intensity, all frames:", r["tests_all_totals"])
    fmt_test("total intensity, excl. 1st cycle:", r["tests_excl_totals"])
    fmt_test("WAXS peak height, all frames:", r["tests_all_peaks"])
    fmt_test("WAXS peak height, excl. 1st cycle:", r["tests_excl_peaks"])


def plot_condition(r, out_dir=OUT_DIR):
    os.makedirs(out_dir, exist_ok=True)
    name = r["name"]
    fig, (ax_tot, ax_peak) = plt.subplots(1, 2, figsize=(12, 5))

    spots_sorted = sorted(r["by_spot_all_totals"])
    x_of_spot = {s: i for i, s in enumerate(spots_sorted)}
    kept_ids = [sid for sid, k in zip(r["scan_ids"], r["keep_all"]) if k]
    dropped_ids = set(cg.first_cycle_ids(name))
    is_dropped = np.array([sid in dropped_ids for sid in kept_ids])

    for ax, values, ylabel in [(ax_tot, r["totals_all"], "total intensity"),
                                (ax_peak, r["peaks_all"], "WAXS peak height")]:
        spot_ids = np.array([sp.SPOT_OF.get(sid, -1) for sid in kept_ids])
        xs = np.array([x_of_spot.get(s, -1) for s in spot_ids])
        valid = xs >= 0
        ax.scatter(xs[valid & ~is_dropped], values[valid & ~is_dropped],
                   color="tab:blue", alpha=0.6, s=20, label="kept")
        ax.scatter(xs[valid & is_dropped], values[valid & is_dropped],
                   color="tab:red", alpha=0.8, s=20, marker="x", label="1st-cycle (excluded)")
        ax.set_xticks(range(len(spots_sorted)))
        ax.set_xticklabels(spots_sorted, fontsize=8)
        ax.set_xlabel("global spot id")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8)

    kw_t = r["tests_excl_totals"]["kw_p"]
    kw_p = r["tests_excl_peaks"]["kw_p"]
    kw_str = (f"KW p (totals, excl. 1st cycle)={kw_t:.3g}" if np.isfinite(kw_t) else "") + \
             (f", (peak)={kw_p:.3g}" if np.isfinite(kw_p) else "")
    fig.suptitle(f"{name}: metric vs stage spot  [{kw_str}]")
    fig.tight_layout()
    out_path = os.path.join(out_dir, f"{name}_spot_effect.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


def plot_summary(results, out_dir=OUT_DIR):
    os.makedirs(out_dir, exist_ok=True)
    names = [r["name"] for r in results]
    eta2_all = [r["tests_all_totals"]["eta2"] for r in results]
    eta2_excl = [r["tests_excl_totals"]["eta2"] for r in results]

    order = np.argsort([-(v if np.isfinite(v) else -1) for v in eta2_excl])
    names = [names[i] for i in order]
    eta2_all = [eta2_all[i] for i in order]
    eta2_excl = [eta2_excl[i] for i in order]

    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(names)), 6))
    width = 0.35
    ax.bar(x - width / 2, [v * 100 if np.isfinite(v) else 0 for v in eta2_all],
           width, label="all frames", color="tab:gray")
    ax.bar(x + width / 2, [v * 100 if np.isfinite(v) else 0 for v in eta2_excl],
           width, label="excl. 1st cycle", color="tab:blue")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("variance in total intensity explained by spot [%]")
    ax.set_title("Position-locking strength by condition (total intensity, eta^2)")
    ax.legend()
    fig.tight_layout()
    out_path = os.path.join(out_dir, "spot_effect_summary.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    results = []
    for name in cg.GROUPS:
        r = analyze_condition(name)
        print_result(r)
        plot_condition(r)
        results.append(r)
    plot_summary(results)
    plt.show()

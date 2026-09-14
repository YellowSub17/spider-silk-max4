"""
Load groups of related scans and plot them against each other for comparison.

GROUPS is built directly from data/scans_clean.csv's `Sample` column -- the
authoritative per-scan condition label -- rather than from free-text log
comments. Sample values follow the convention `[CHANNEL1]/[CHANNEL2]`: the
two solutions being merged in the mixing channel, e.g. `Y2F+kpi/aa` is a
solution of Y2F protein in KPi buffer, mixed in-channel with acetic acid;
`Y2F+kpi/h2o` is the same protein solution mixed with water instead. A
protein-free run like `kpi/aa` is KPi buffer alone mixed with acetic acid --
the matched blank for `Y2F+kpi/aa`. `cap:`-prefixed samples are static
capillary measurements (some also carrying a shear/no-shear or gel-state
third field), a different measurement modality from the flow-channel mixing
runs and not part of the channel1/channel2 blank-matching logic below.

Repeats within a group are combined via ComboScan; COMPARISONS defines which
groups get overlaid on the same axes.

Run directly: `python compare_groups.py`
Figures are written to ./figures/.
"""

import csv
import os
import re

import matplotlib.pyplot as plt
import numpy as np

import src

SCANS_CLEAN_CSV = f"{src.const.DATA_PATH}/scans_clean.csv"

MIN_SCANS_PER_GROUP = 3  # drop conditions with too few scans to be usable


def _sanitize(name):
    """Turn a raw Sample string into a usable Python/filename-safe identifier."""
    return re.sub(r"[^0-9a-zA-Z]+", "_", name).strip("_")


def load_sample_column(path=SCANS_CLEAN_CSV):
    """Map ScanID -> raw Sample string, for every row where Sample is filled in."""
    samples = {}
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            sid = row["scanID"].strip()
            sample = row["Sample"].strip()
            if sid and sample:
                samples[int(sid)] = sample
    return samples


def parse_sample(raw_sample):
    """
    Break a raw Sample string into its meaning: protein (if any), the two
    mixing-channel contents, whether it's a static capillary measurement, and
    whether acetic acid is involved anywhere. Returns a dict; fields that
    don't apply are None.
    """
    is_capillary = raw_sample.startswith("cap:")
    body = raw_sample[len("cap:"):] if is_capillary else raw_sample
    parts = body.split("/")

    protein = None
    channel1 = parts[0] if parts else None
    if channel1 and "+" in channel1:
        protein, channel1 = channel1.split("+", 1)

    channel2 = parts[1] if len(parts) > 1 else None
    extra = parts[2] if len(parts) > 2 else None  # e.g. shear/no-shear/gel

    has_acid = any(tok in ("aa", "aca") for tok in (channel1, channel2, extra) if tok)

    return {
        "raw": raw_sample, "is_capillary": is_capillary, "protein": protein,
        "channel1": channel1, "channel2": channel2, "extra": extra,
        "has_acid": has_acid,
    }


def find_blank(raw_sample, all_raw_samples):
    """
    For a protein-containing, non-capillary Sample, find the matching
    protein-free blank: the same two channel contents (in either order --
    the final mixed composition is the same regardless of which channel the
    protein solution started in), just without a protein prefix. Returns the
    raw blank Sample string, or None if no matching blank exists in the data.
    """
    info = parse_sample(raw_sample)
    if info["protein"] is None or info["is_capillary"]:
        return None
    c1, c2 = info["channel1"], info["channel2"]
    for candidate in (f"{c1}/{c2}", f"{c2}/{c1}"):
        if candidate in all_raw_samples:
            return candidate
    return None


# --- Build GROUPS from the Sample column -----------------------------------
_sample_by_scan = load_sample_column()
_scans_by_raw_sample = {}
for _sid, _raw in _sample_by_scan.items():
    _scans_by_raw_sample.setdefault(_raw, []).append(_sid)

# name (sanitized) -> sorted scan-id list; SAMPLE_INFO/RAW_NAME let you get
# back the parsed meaning and original Sample string for a given group name.
GROUPS = {}
SAMPLE_INFO = {}
RAW_NAME = {}
for _raw, _sids in _scans_by_raw_sample.items():
    if len(_sids) < MIN_SCANS_PER_GROUP:
        continue
    _name = _sanitize(_raw)
    GROUPS[_name] = sorted(_sids)
    SAMPLE_INFO[_name] = parse_sample(_raw)
    RAW_NAME[_name] = _raw

# name -> matched protein-free blank's group name (or None), for every
# flow-channel (non-capillary) group in GROUPS.
BLANK_OF = {}
for _name, _raw in RAW_NAME.items():
    _blank_raw = find_blank(_raw, _scans_by_raw_sample)
    BLANK_OF[_name] = _sanitize(_blank_raw) if _blank_raw else None

# Conditions dropped after drift_within_runs.py found them internally
# self-inconsistent (early-vs-late relative RMSD, i0-normalized, no q-pinning):
#   Y2F_kpi_kpi: 108.8% RMSD, total intensity swung +387% within 55 minutes
#     (a real mid-run event -- clogging/aggregation/bubble -- not steady drift).
#   resin: 46.6% RMSD over just 5 min/10 scans, a clean step-change between
#     the first 3 scans and the rest (likely a beam/position shift mid-run).
# Neither was usable as a background/reference anyway; both are excluded from
# GROUPS (and everything downstream that iterates over it) entirely.
DROPPED_CONDITIONS = {
    "Y2F_kpi_kpi": "internally inconsistent (108.8% early/late RMSD, see drift_within_runs.py)",
    "resin": "internally inconsistent (46.6% early/late RMSD, step-change mid-run)",
}
for _dropped in DROPPED_CONDITIONS:
    GROUPS.pop(_dropped, None)
    SAMPLE_INFO.pop(_dropped, None)
    RAW_NAME.pop(_dropped, None)
    BLANK_OF.pop(_dropped, None)


# --- Groups worth comparing to each other ----------------------------------
COMPARISONS = [
    ("Buffer effect on Y2F", ["Y2F_h2o_h2o", "Y2F_h2o_kpi", "Y2F_kpi_aa"]),
    ("Y2F vs YR2A in water", ["Y2F_h2o_h2o", "YR2A_h2o_h2o"]),
    ("Y2F vs YR2A in KPi", ["Y2F_h2o_kpi", "YR2A_h2o_kpi"]),
    ("Buffer-only controls", ["h2o_h2o", "kpi_h2o", "kpi_aa"]),
    # Y2F_kpi_kpi (the "no acid" baseline) was dropped -- see DROPPED_CONDITIONS
    # -- so this is now just protein+KPi+acid vs the matched protein-free blank.
    ("Y2F in KPi with vs without acid", ["Y2F_kpi_aa", "kpi_aa"]),
]

RSQ_LIM = 0.999  # keep only frames with a clean log-log linear streak fit

# Neutral q-window used by the "neutral" normalization and by the scaling-factor
# diagnostic in protein_vs_buffer.py: mid-q, away from both the Guinier region
# and the WAXS peaks, where protein and blank curves should look similar.
Q_NEUTRAL = (0.1, 1.0)

# Normalization strategies. Each mutates combo.Is in place (same convention as
# ComboScan.norm_qrange/norm_max/norm_i0). r_squared (the streak-quality metric)
# is computed before any of these run and is invariant to per-frame scalar
# rescaling, so RSQ_LIM filtering is unaffected by which of these is used.
NORM_METHODS = {
    # Pin the low-q Porod/streak region -- the original approach used throughout
    # this analysis. Assumes that region is the same between groups.
    "linear_saxs": lambda combo: combo.norm_qrange(*src.det.LINEAR_SAXS_RANGE),
    # Pin the mid-q "neutral" region instead -- assumes *that* region is the
    # same between groups (reasonable if it's dominated by solvent/background
    # rather than the low-q Guinier signal or the high-q protein peaks).
    "neutral": lambda combo: combo.norm_qrange(*Q_NEUTRAL),
    # Normalize by incident flux (i0) only -- no assumption about which q
    # region should match between groups, just corrects for delivered beam
    # intensity per frame.
    "i0": lambda combo: combo.norm_i0(),
}


# Whether load_group()/drift_check() drop each condition's "first cycle" --
# the earliest N chronological scans, where N is the number of distinct
# physical stage spots (stage_positions.spots_of_group) that condition
# visits -- by default. Motivated by drift_within_runs.py: the first full
# position cycle of a run is consistently a settling transient (e.g.
# YR2A_h2o_h2o's total intensity: 12.2, 12.9, 18.6, 18.0, then ~14 for the
# rest of the run). Override per-call with the skip_first_cycle argument.
SKIP_FIRST_CYCLE = True


def first_cycle_ids(name):
    """
    The scan IDs comprising a condition's first full position cycle: the N
    chronologically-earliest scans, where N = the number of distinct
    physical stage spots (see stage_positions.py) that condition's scans
    fall into. Returns [] if that count is 0 or >= the condition's total
    scan count (nothing sensible to drop).

    Lazily imports stage_positions (which itself imports this module) to
    avoid a circular import at module load time.
    """
    import stage_positions  # noqa: local import, see docstring
    ids = sorted(GROUPS[name])
    n_spots = len(stage_positions.spots_of_group(name))
    if n_spots <= 0 or n_spots >= len(ids):
        return []
    return ids[:n_spots]


def report_first_cycle_skips():
    """Print, for every condition in GROUPS, how many scans first_cycle_ids
    would drop -- i.e. what SKIP_FIRST_CYCLE actually removes per condition."""
    print(f"{'condition':<20}{'n_total':>9}{'n_spots':>9}{'n_dropped':>11}")
    import stage_positions
    for name in GROUPS:
        n_total = len(GROUPS[name])
        n_spots = len(stage_positions.spots_of_group(name))
        dropped = first_cycle_ids(name)
        print(f"{name:<20}{n_total:>9}{n_spots:>9}{len(dropped):>11}")


def load_scan_ids(scan_ids, load_n=1, label=None):
    """Load and combine an arbitrary list of scan IDs, skipping any that fail
    to load (a handful of scans have corrupt/incomplete raw HDF5 groups)."""
    all_ids = list(scan_ids)
    good_ids = []
    skipped = []
    for sid in all_ids:
        try:
            src.Scan(sid, load_imgs=False, load_n=load_n)
        except Exception:
            skipped.append(sid)
        else:
            good_ids.append(sid)
    if skipped:
        print(f"  skipping failed/incomplete scans: {skipped}")
    tag = f"'{label}' " if label else ""
    print(f"Loading {tag}({len(good_ids)}/{len(all_ids)} scans)...")
    return src.ComboScan(good_ids, load_imgs=False, load_n=load_n)


def load_group(name, load_n=1, skip_first_cycle=None):
    """
    Load and combine all scans in a named GROUPS entry.

    skip_first_cycle (default: SKIP_FIRST_CYCLE) drops that condition's
    settling-transient first cycle (see first_cycle_ids) BEFORE loading --
    so the dropped scans are never even read, and every downstream
    kept_frames/mean_curve call on the resulting ComboScan automatically
    excludes them without needing its own skip logic.
    """
    if skip_first_cycle is None:
        skip_first_cycle = SKIP_FIRST_CYCLE
    ids = sorted(GROUPS[name])
    if skip_first_cycle:
        dropped = first_cycle_ids(name)
        if dropped:
            print(f"  skipping first-cycle settling transient for '{name}': "
                  f"{len(dropped)} scans {dropped}")
            ids = [i for i in ids if i not in dropped]
    return load_scan_ids(ids, load_n=load_n, label=name)


def kept_mask(combo, rsq_lim=RSQ_LIM, skip_first_n=0):
    """
    Boolean mask (over combo.Is / combo.scan_ids rows, in load order) of
    frames that pass the streak-quality cut, with skip_first_n additionally
    excluded. Factored out of kept_frames so callers that need to line up
    kept frames with their originating scan IDs (e.g. position_locking.py,
    matching per-frame metrics to stage spots) can get the same mask
    kept_frames/mean_curve would use, without duplicating the logic.
    """
    keep = combo.r_squared > rsq_lim
    if skip_first_n > 0:
        keep = keep.copy()
        keep[:skip_first_n] = False
    if not np.any(keep):
        print(f"  warning: no frames passed r_squared > {rsq_lim}, using all frames")
        keep = np.ones_like(combo.r_squared, dtype=bool)
    return keep


def kept_frames(combo, rsq_lim=RSQ_LIM, norm="linear_saxs", skip_first_n=0):
    """Normalize and return the individual per-frame I(q) curves that pass the
    streak-quality cut (r_squared > rsq_lim), i.e. the replicate curves that
    mean_curve() would otherwise collapse into a single mean+std.

    `norm` selects one of NORM_METHODS (default: the original low-q pinning).

    `skip_first_n` additionally excludes the first N frames of the (already
    chronologically-ordered) combo -- a defensive/manual equivalent of
    load_group's skip_first_cycle, for combos built without going through
    load_group (e.g. a manually-supplied scan-id list). When a combo came
    from load_group with its default skip_first_cycle=True, those frames
    are already absent and skip_first_n should stay 0 (the default).

    Returns qs, Is_kept (n_kept x n_q), n_kept, n_total.
    """
    NORM_METHODS[norm](combo)
    keep = kept_mask(combo, rsq_lim, skip_first_n)
    return combo.qs, combo.Is[keep], keep.sum(), len(keep)


def mean_curve(combo, rsq_lim=RSQ_LIM, norm="linear_saxs", skip_first_n=0):
    """Normalize and average I(q) over frames that pass the streak-quality cut.

    Returns qs, mean, std, n_kept, n_total. std is the frame-to-frame (run-to-run)
    standard deviation at each q -- i.e. how much the repeats within this group
    disagree with each other, not measurement error on a single frame.
    """
    qs, Is_kept, n_kept, n_total = kept_frames(combo, rsq_lim, norm=norm, skip_first_n=skip_first_n)
    return qs, Is_kept.mean(axis=0), Is_kept.std(axis=0), n_kept, n_total


def coeff_of_variation(qs, I_mean, I_std, qrange=None):
    """Median run-to-run coefficient of variation (std/mean) over qrange."""
    if qrange is None:
        qrange = src.det.FULL_QRANGE
    qloc = (qs >= qrange[0]) & (qs <= qrange[1]) & (I_mean > 0)
    return np.median(I_std[qloc] / I_mean[qloc])


# Baseline windows for peak_analysis when applied to a *difference* curve
# rather than a raw I(q) curve -- anchored close to the WAXS peak itself
# rather than peak_analysis's raw-curve defaults, which land inside the
# low-q Guinier-mismatch region on a difference curve and badly distort the
# baseline (see protein_vs_buffer.py's module docstring for the history).
DIFF_Q_BASE = ((0.5, 0.7), (1.5, 1.7))

MIN_SCANS_FOR_DRIFT_CHECK = 20  # need enough scans in each half to be meaningful


def drift_check(name, norm="i0", q_base=DIFF_Q_BASE, skip_first_cycle=None):
    """
    Split a named group's scans into chronological early/late halves (scan ID
    order is a good proxy for time order here) and run the same
    scale-corrected WAXS-region peak analysis between them that
    protein_vs_buffer.py uses for a protein-vs-blank comparison -- but with
    NO protein or condition difference at all, just time.

    Use this to sanity-check a protein-vs-blank result: if a condition's own
    early-vs-late split produces a "peak" of comparable size to the claimed
    protein signature, the run-to-run/session time gap between the protein
    and blank measurements could be producing (or contributing to) that
    signature by itself, rather than real protein structure.

    skip_first_cycle (default: SKIP_FIRST_CYCLE) drops the settling-transient
    first cycle from the full id list before splitting into early/late, same
    as load_group -- otherwise the transient itself (which is large and
    concentrated at the very start) could dominate the "early" half and
    inflate this check's own result.

    Returns dict(scale, peak, n_early, n_late) or None if there aren't enough
    scans in the group to split meaningfully (see MIN_SCANS_FOR_DRIFT_CHECK).
    peak is whatever src.saxs_metrics.peak_analysis returns (or None).
    """
    if skip_first_cycle is None:
        skip_first_cycle = SKIP_FIRST_CYCLE
    ids = sorted(GROUPS[name])
    if skip_first_cycle:
        dropped = first_cycle_ids(name)
        ids = [i for i in ids if i not in dropped]
    if len(ids) < MIN_SCANS_FOR_DRIFT_CHECK:
        return None
    half = len(ids) // 2
    early, late = ids[:half], ids[half:]

    combo_e = load_scan_ids(early, label=f"{name} (early half, drift check)")
    combo_l = load_scan_ids(late, label=f"{name} (late half, drift check)")
    qs_e, Is_e, n_e, _ = kept_frames(combo_e, norm=norm)
    qs_l, Is_l, n_l, _ = kept_frames(combo_l, norm=norm)
    if n_e == 0 or n_l == 0:
        return None
    I_e, I_l = Is_e.mean(axis=0), Is_l.mean(axis=0)

    mask = (qs_e >= Q_NEUTRAL[0]) & (qs_e <= Q_NEUTRAL[1])
    denom = np.sum(I_e[mask] ** 2)
    scale = np.sum(I_l[mask] * I_e[mask]) / denom if denom > 0 else np.nan

    diff = I_l - scale * I_e
    peak = src.saxs_metrics.peak_analysis(qs_e, diff, q_base=q_base)
    return {"scale": scale, "peak": peak, "n_early": n_e, "n_late": n_l}


def plot_comparison(title, group_names, out_dir="figures"):
    os.makedirs(out_dir, exist_ok=True)

    fig, (ax, ax_cv) = plt.subplots(2, 1, figsize=(7, 8), sharex=True,
                                     gridspec_kw={"height_ratios": [3, 1]})
    cmap = plt.colormaps.get_cmap("tab10")

    print(f"\n{title}")
    for i, name in enumerate(group_names):
        combo = load_group(name)
        qs, I_mean, I_std, n_keep, n_total = mean_curve(combo)
        cv = coeff_of_variation(qs, I_mean, I_std)
        print(f"  {name}: median run-to-run CV (std/mean) = {cv:.2%} over {n_keep} frames")

        color = cmap(i)
        label = f"{name} ({n_keep}/{n_total}, CV={cv:.0%})"
        ax.plot(qs, I_mean, color=color, label=label)
        ax.fill_between(qs, I_mean - I_std, I_mean + I_std, color=color, alpha=0.2, linewidth=0)

        ax_cv.plot(qs, I_std / np.where(I_mean > 0, I_mean, np.nan), color=color)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Inten. (norm.)")
    ax.set_xlim(src.det.FULL_QRANGE)
    ax.set_title(title)
    ax.legend()

    ax_cv.set_yscale("log")
    ax_cv.set_xlabel("q [1/A]")
    ax_cv.set_ylabel("run-to-run\nstd / mean")
    fig.tight_layout()

    out_path = os.path.join(out_dir, title.lower().replace(" ", "_") + ".png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    for title, group_names in COMPARISONS:
        plot_comparison(title, group_names)
    plt.show()

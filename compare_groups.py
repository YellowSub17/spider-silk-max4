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


# --- Groups worth comparing to each other ----------------------------------
COMPARISONS = [
    ("Buffer effect on Y2F", ["Y2F_h2o_h2o", "Y2F_h2o_kpi", "Y2F_kpi_aa"]),
    ("Y2F vs YR2A in water", ["Y2F_h2o_h2o", "YR2A_h2o_h2o"]),
    ("Y2F vs YR2A in KPi", ["Y2F_h2o_kpi", "YR2A_h2o_kpi"]),
    ("Buffer-only controls", ["h2o_h2o", "kpi_h2o", "kpi_aa"]),
    ("Acetic acid progression on Y2F", ["Y2F_kpi_kpi", "Y2F_kpi_aa", "kpi_aa"]),
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


def load_group(name, load_n=1):
    """Load and combine all scans in a named GROUPS entry."""
    return load_scan_ids(GROUPS[name], load_n=load_n, label=name)


def kept_frames(combo, rsq_lim=RSQ_LIM, norm="linear_saxs"):
    """Normalize and return the individual per-frame I(q) curves that pass the
    streak-quality cut (r_squared > rsq_lim), i.e. the replicate curves that
    mean_curve() would otherwise collapse into a single mean+std.

    `norm` selects one of NORM_METHODS (default: the original low-q pinning).

    Returns qs, Is_kept (n_kept x n_q), n_kept, n_total.
    """
    NORM_METHODS[norm](combo)
    keep = combo.r_squared > rsq_lim
    if not np.any(keep):
        print(f"  warning: no frames passed r_squared > {rsq_lim}, using all frames")
        keep = np.ones_like(combo.r_squared, dtype=bool)
    return combo.qs, combo.Is[keep], keep.sum(), len(keep)


def mean_curve(combo, rsq_lim=RSQ_LIM, norm="linear_saxs"):
    """Normalize and average I(q) over frames that pass the streak-quality cut.

    Returns qs, mean, std, n_kept, n_total. std is the frame-to-frame (run-to-run)
    standard deviation at each q -- i.e. how much the repeats within this group
    disagree with each other, not measurement error on a single frame.
    """
    qs, Is_kept, n_kept, n_total = kept_frames(combo, rsq_lim, norm=norm)
    return qs, Is_kept.mean(axis=0), Is_kept.std(axis=0), n_kept, n_total


def coeff_of_variation(qs, I_mean, I_std, qrange=None):
    """Median run-to-run coefficient of variation (std/mean) over qrange."""
    if qrange is None:
        qrange = src.det.FULL_QRANGE
    qloc = (qs >= qrange[0]) & (qs <= qrange[1]) & (I_mean > 0)
    return np.median(I_std[qloc] / I_mean[qloc])


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

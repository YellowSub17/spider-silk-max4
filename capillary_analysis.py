"""
Separate analysis pipeline for the static capillary screen (Sample values
starting "cap:" in data/scans_clean.csv), kept apart from the flow-channel
pipeline (compare_groups.py / protein_vs_buffer.py / position_series.py)
because it's a different measurement modality:
  - flow-channel runs: two solutions mixed in a channel, measured at several
    positions along the channel (a proxy for time-since-mixing), tens to
    hundreds of repeat exposures per condition.
  - capillary runs: individual, static, pre-made samples loaded one per
    rotator slot ("Sample N" = which physical capillary, not a mixing
    position), only ~2 repeat exposures per condition.

This screen also covers a THIRD protein construct not present in the flow
data: "WT" = wild-type, i.e. natural (non-recombinant) spider silk protein,
alongside the two recombinant constructs Y2F ("YF" in these labels) and
YR2A ("YRA"). It also records a shear vs no-shear mechanical treatment, and
one "gel" state label -- both absent from the flow-channel screen.

Sample naming: `cap:[PROTEIN+]CHANNEL1[/CHANNEL2][/EXTRA]`, e.g.:
  cap:YF/kpi/shear     -> Y2F in KPi buffer, sheared
  cap:WT/kpi/no-shear  -> wild-type silk in KPi buffer, not sheared
  cap:YF/kpi/aa        -> Y2F in KPi buffer, acetic-acid exposed
  cap:YF/kpi/gel       -> Y2F in KPi buffer, visibly gelled
  cap:kpi, cap:kpi/aca -> buffer alone (+/- acetic acid), no protein
  cap:empty            -> empty capillary (instrument background)

Because each condition only has ~2 exposures, this script does NOT try to
build the same run-to-run statistics (CV, replicate envelopes) as the flow
pipeline -- there isn't enough data for that. Instead it:
  1. Maps out every capillary condition (protein, buffer, shear state).
  2. Repeats the universal-Rg-artifact check from map_llps_by_rg.py in this
     completely different (static, no flow) measurement geometry -- if the
     same ~500 A apparent Rg shows up here too, that's strong independent
     evidence it's a beamline/detector artifact rather than anything to do
     with the flow channel; if it does NOT show up, that would instead point
     at the flow channel/cell itself as the source.
  3. Plots direct comparisons that this dataset uniquely enables: protein
     identity (Y2F vs YR2A vs WT) at matched buffer/shear conditions, shear
     vs no-shear per protein, and the fullest single-protein context (Y2F
     across buffer/shear/acid/gel states).

Run directly: `python capillary_analysis.py`
Figures are written to ./figures/capillary/.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np

import compare_groups as cg
import src

OUT_DIR = "figures/capillary"
NORM = "neutral"

# Y2F/YR2A were called "YF"/"YRA" in the capillary labels; unify for display.
PROTEIN_ALIAS = {"YF": "Y2F", "YRA": "YR2A", "WT": "WT"}


def load_capillary_samples():
    """All cap: Sample values -> sorted scan-id list, no minimum-count filter
    (unlike compare_groups.GROUPS, since these conditions only have ~2-7
    scans each)."""
    return {
        raw: sorted(ids)
        for raw, ids in cg._scans_by_raw_sample.items()
        if raw.startswith("cap:")
    }


CAP_SAMPLES = load_capillary_samples()          # raw Sample -> scan ids
CAP_INFO = {s: cg.parse_sample(s) for s in CAP_SAMPLES}  # raw Sample -> parsed dict


def display_protein(info):
    """
    compare_groups.parse_sample() only recognizes a protein when it's written
    as "PROTEIN+channel1" (the flow-data convention). Capillary labels instead
    put the protein directly as channel1 (e.g. "cap:YF/kpi/shear"), so check
    there too.
    """
    if info["protein"]:
        return PROTEIN_ALIAS.get(info["protein"], info["protein"])
    if info["channel1"] in PROTEIN_ALIAS:
        return PROTEIN_ALIAS[info["channel1"]]
    return None


def print_overview():
    print(f"{len(CAP_SAMPLES)} capillary conditions, "
          f"{sum(len(v) for v in CAP_SAMPLES.values())} scans total:\n")
    print(f"{'Sample':<24}{'n':>4}{'protein':>9}{'channel1':>10}{'channel2':>10}{'extra':>10}")
    for raw, ids in sorted(CAP_SAMPLES.items(), key=lambda kv: -len(kv[1])):
        info = CAP_INFO[raw]
        print(f"{raw:<24}{len(ids):>4}{str(display_protein(info)):>9}"
              f"{str(info['channel1']):>10}{str(info['channel2']):>10}{str(info['extra']):>10}")


def load_curves(raw_sample, norm=NORM):
    """Load+normalize+streak-filter all frames for one capillary condition.
    Returns qs, Is_kept (n_kept x n_q) -- individual curves, not just a mean,
    since n is small enough that overlaying them directly is more honest than
    collapsing to mean+-std."""
    combo = cg.load_scan_ids(CAP_SAMPLES[raw_sample], label=raw_sample)
    qs, Is_kept, n_kept, n_total = cg.kept_frames(combo, norm=norm)
    return qs, Is_kept, n_kept, n_total


# --- 1. Universal-Rg-artifact check, repeated in this different geometry ---

def rg_artifact_check(out_dir=OUT_DIR):
    os.makedirs(out_dir, exist_ok=True)
    print("\n--- Guinier Rg check across capillary conditions ---")
    rows = []
    for raw in CAP_SAMPLES:
        qs, Is_kept, n_kept, n_total = load_curves(raw)
        if n_kept == 0:
            continue
        I_mean = Is_kept.mean(axis=0)
        g = src.saxs_metrics.guinier_fit(qs, I_mean, q_max=0.02)
        rows.append({
            "sample": raw, "n_kept": n_kept,
            "Rg": g["Rg"] if g else np.nan,
            "r2": g["r"] ** 2 if g else np.nan,
            "n_points": g["n_points"] if g else 0,
        })

    print(f"{'Sample':<24}{'Rg (A)':>10}{'r^2':>8}{'n_pts':>7}{'n_kept':>8}")
    for r in sorted(rows, key=lambda r: (-r["Rg"] if np.isfinite(r["Rg"]) else np.inf)):
        rg_str = f"{r['Rg']:.1f}" if np.isfinite(r["Rg"]) else "fit failed"
        r2_str = f"{r['r2']:.4f}" if np.isfinite(r["r2"]) else "-"
        print(f"{r['sample']:<24}{rg_str:>10}{r2_str:>8}{r['n_points']:>7}{r['n_kept']:>8}")

    valid = [r for r in rows if np.isfinite(r["Rg"])]
    n_large = sum(1 for r in valid if r["Rg"] > 200)
    print(f"\n{n_large}/{len(valid)} capillary conditions show Rg > 200 A "
          f"(including cap:empty: "
          f"{'yes' if any(r['sample']=='cap:empty' and r['Rg']>200 for r in valid) else 'no'})")
    if n_large == len(valid) and len(valid) > 0:
        print("-> Same universal artifact seen in the flow-channel data ALSO appears in "
              "this static capillary geometry -- strong evidence it's a beamline/detector "
              "artifact (shared low-q feature), not specific to the flow cell.")
    elif n_large == 0:
        print("-> The artifact does NOT appear in capillary data -- points at the "
              "flow channel/cell itself as the source, not the beamline.")
    else:
        print("-> Mixed result -- inspect which specific conditions show it.")

    # Bar plot, same style as map_llps_by_rg.py's rg_by_condition.png
    valid.sort(key=lambda r: r["Rg"])
    names = [r["sample"] for r in valid]
    rgs = [r["Rg"] for r in valid]
    colors = ["tab:red" if CAP_INFO[n]["has_acid"] else
              ("tab:gray" if n in ("cap:empty",) else "tab:blue") for n in names]

    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(names)), 6))
    x = np.arange(len(names))
    ax.bar(x, rgs, color=colors)
    ax.axhline(500, color="k", linestyle="--", linewidth=1, label="Rg=500 A (flow-data artifact level)")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Apparent Rg [A]")
    ax.set_title("Guinier Rg by capillary condition")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out_path = os.path.join(out_dir, "cap_rg_by_condition.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


# --- 2. Direct comparisons this dataset uniquely enables --------------------

def plot_overlay(title, raw_samples, out_name, out_dir=OUT_DIR, labels=None):
    """Overlay individual (not just mean) curves for a list of raw Sample
    values -- useful here since n~2 per condition."""
    os.makedirs(out_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.colormaps.get_cmap("tab10")

    print(f"\n{title}")
    for i, raw in enumerate(raw_samples):
        if raw not in CAP_SAMPLES:
            print(f"  (skipping '{raw}': not present in this dataset)")
            continue
        qs, Is_kept, n_kept, n_total = load_curves(raw)
        color = cmap(i)
        label = labels[i] if labels else raw
        for j, I in enumerate(Is_kept):
            ax.plot(qs, I, color=color, alpha=0.5, linewidth=1,
                    label=f"{label} (n={n_kept})" if j == 0 else None)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("q [1/A]")
    ax.set_ylabel("Inten. (norm.)")
    ax.set_xlim(src.det.FULL_QRANGE)
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()

    out_path = os.path.join(out_dir, out_name)
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


PROTEIN_COMPARISONS = [
    ("Protein identity in KPi, no shear",
     ["cap:YRA/kpi/no-shear", "cap:YF/kpi/no-shear", "cap:WT/kpi/no-shear"],
     ["YR2A", "Y2F", "WT (natural silk)"], "cap_protein_kpi_noshear.png"),
    ("Protein identity in KPi, sheared",
     ["cap:YRA/kpi/shear", "cap:YF/kpi/shear", "cap:WT/kpi/shear"],
     ["YR2A", "Y2F", "WT (natural silk)"], "cap_protein_kpi_shear.png"),
    ("Protein identity in water",
     ["cap:YRA/h2o", "cap:YF/h2o", "cap:WT/h2o"],
     ["YR2A", "Y2F", "WT (natural silk)"], "cap_protein_h2o.png"),
]

SHEAR_COMPARISONS = [
    ("YR2A: shear vs no-shear (KPi)",
     ["cap:YRA/kpi/no-shear", "cap:YRA/kpi/shear"], ["no-shear", "shear"], "cap_shear_yr2a.png"),
    ("Y2F: shear vs no-shear (KPi)",
     ["cap:YF/kpi/no-shear", "cap:YF/kpi/shear"], ["no-shear", "shear"], "cap_shear_y2f.png"),
    ("WT: shear vs no-shear (KPi)",
     ["cap:WT/kpi/no-shear", "cap:WT/kpi/shear"], ["no-shear", "shear"], "cap_shear_wt.png"),
]

Y2F_CONTEXT = (
    "Y2F across capillary conditions",
    ["cap:YF/h2o", "cap:YF/kpi/no-shear", "cap:YF/kpi/shear", "cap:YF/kpi/gel", "cap:YF/kpi/aa"],
    ["water", "KPi, no shear", "KPi, sheared", "KPi, gelled", "KPi + acetic acid"],
    "cap_y2f_context.png",
)

BUFFER_ACID_CHECK = (
    "Buffer alone: KPi vs KPi+acetic acid (no protein)",
    ["cap:kpi", "cap:kpi/aca"], ["KPi", "KPi + acetic acid"], "cap_buffer_acid.png",
)


if __name__ == "__main__":
    print_overview()
    rg_artifact_check()

    for title, samples, labels, out_name in PROTEIN_COMPARISONS:
        plot_overlay(title, samples, out_name, labels=labels)
    for title, samples, labels, out_name in SHEAR_COMPARISONS:
        plot_overlay(title, samples, out_name, labels=labels)
    title, samples, labels, out_name = Y2F_CONTEXT
    plot_overlay(title, samples, out_name, labels=labels)
    title, samples, labels, out_name = BUFFER_ACID_CHECK
    plot_overlay(title, samples, out_name, labels=labels)

    plt.show()

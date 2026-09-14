"""
Isolate the protein's contribution to the scattering by comparing each
protein-containing run against its matching protein-free blank (buffer alone,
or water alone) run: same flow-cell setup, same buffer chemistry, minus the
protein.

Because a difference curve this small can easily be an artifact of scaling,
noise, or normalization rather than real protein structure, for each pair
this script runs a set of reliability diagnostics before trusting the
difference:

  1. Run-to-run stability (CV): coefficient of variation within each group.
     If the protein group is much noisier than its blank (ratio >> 1), any
     "protein signature" may just be picking up that instability.
  2. Scaling factor test: fits an optimal multiplicative scale s such that
     s * blank ~= protein over a "neutral" mid-q window (0.1-1.0 1/A, away
     from both the Guinier region and the WAXS peaks), where real protein
     and buffer structure should look similar. s far from 1.0 flags a
     normalization/scaling mismatch between the two groups.
  3. Peak significance: how many sigma the WAXS difference peak sits above
     zero, using the pooled SEM of both groups. <3 sigma isn't resolvable.
  4. Replicate-level subtraction: instead of only subtracting mean curves,
     each individual protein replicate is subtracted from the blank mean
     (and vice versa) to build an envelope of "what if we'd used a
     different single run" -- shown as a shaded band under the diff curve.
     A peak that's much smaller than this envelope isn't reliable.
  5. Relative difference curve (protein-blank)/blank, which is less
     sensitive to a shared scaling artifact than the absolute difference.
  6. Scale-corrected peak significance: repeats the peak test on
     (protein - scale*blank) instead of (protein - blank), using the fitted
     scale from #2. This removes any overall multiplicative offset between
     the two groups (real or a normalization artifact) before looking for a
     *localized* WAXS peak, so a broad Porod-region mismatch can't masquerade
     as -- or mask -- a real structural peak.
  7. Drift check (compare_groups.drift_check): the protein and blank runs are
     often collected tens of minutes to hours apart (sometimes with other
     conditions run in between), so an apparent "protein peak" could really
     be session drift. This splits the *blank* condition's own scans into
     chronological early/late halves and runs the identical scale-corrected
     peak test between them -- i.e. the same analysis, with no protein and no
     condition change at all, just time. If that alone produces a peak
     that's a large fraction of the claimed protein-vs-blank peak, treat the
     result with real suspicion (see the Y2F_kpi_h2o vs kpi_h2o case, where
     an early/late split of kpi_h2o alone reproduced ~75% of the claimed
     peak height).
  8. Position-matched subtraction (position_matched_diff): stage_positions.py
     found that the "position N" rotator label is NOT a shared physical
     location across samples/conditions, that mid-run intensity cycles are
     genuinely position-locked (position_locking.py), and that at least one
     flow group's frames are skewed across physical spots relative to its
     blank. The pooled subtraction above (mean of all protein frames minus
     mean of all blank frames) therefore silently assumes the two groups
     sampled the same mix of physical spots -- if they didn't, the "protein"
     difference partly (or entirely) reflects that position-mix mismatch,
     not protein structure. This alternative instead averages protein and
     blank separately WITHIN each shared physical spot, then combines the
     per-spot differences with equal weight -- removing the position-mix
     bias by construction -- and separately reports each shared spot's own
     peak significance (a peak present at only one spot is geometry; one
     consistent across all shared spots is real). Reported alongside, never
     replacing, the pooled result above; skipped outright (not forced to a
     nearest-spot match) when coverage between the two groups is too poor.

Run directly: `python protein_vs_buffer.py`
Figures are written to ./figures/.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

import compare_groups as cg
import src
import stage_positions as sp
from position_locking import spot_groups, MIN_FRAMES_PER_SPOT

# (title, protein_group, blank_group) -- blank is the matching protein-free run
# (same two mixing-channel contents, minus the protein) that isolates the
# protein's contribution when subtracted. Built from compare_groups.BLANK_OF,
# which matches purely on the parsed Sample string (see compare_groups.py's
# module docstring for the [CHANNEL1]/[CHANNEL2] naming convention) -- so this
# automatically covers every protein-containing flow-channel condition that
# has a matching protein-free run in the data, skipping ones that don't
# (e.g. Y2F_kpi_kpi has no "kpi/kpi" blank available).
PAIRS = [
    (f"{name} vs {cg.BLANK_OF[name]} blank", name, cg.BLANK_OF[name])
    for name in cg.GROUPS
    if cg.SAMPLE_INFO[name]["protein"] and cg.BLANK_OF[name] is not None
]

# Neutral window for the scaling-factor fit: mid-q, away from the low-q
# Guinier/Porod region and the high-q WAXS peaks, where protein and blank
# should scatter similarly if properly scaled.
Q_NEUTRAL = (0.1, 1.0)

CV_RATIO_WARN = 2.0     # protein CV / blank CV above this -> flag
SCALE_DEVIATION_WARN = 0.05  # |s - 1| above this -> flag
SIGMA_WARN = 3.0        # peak significance below this -> flag


def _flag(is_ok):
    return "✓" if is_ok else "⚠ WARNING"


def fit_scale_factor(qs, I_p, I_b, q_range=Q_NEUTRAL):
    """Least-squares scale s minimizing ||s*I_b - I_p||^2 over q_range."""
    mask = (qs >= q_range[0]) & (qs <= q_range[1])
    b = I_b[mask]
    p = I_p[mask]
    denom = np.sum(b ** 2)
    if denom <= 0:
        return np.nan
    return np.sum(p * b) / denom


# src.saxs_metrics.peak_analysis's default baseline anchors (q=0.03-0.08 and
# q=1.9-2.1) are tuned for a raw I(q) curve, where that low-q window is flat.
# On a *difference* curve that window instead sits inside the region most
# dominated by leftover low-q Guinier mismatch (order 1-10 vs the ~1e-4 WAXS
# signal), so the baseline interpolation across the peak window ends up wildly
# sloped -- inflating (or hiding) the apparent peak area/significance. Anchor
# to windows immediately flanking the WAXS peak(s) instead. (Shared with
# compare_groups.drift_check, which uses the exact same peak test.)
DIFF_Q_BASE = cg.DIFF_Q_BASE


def resolve_peak(qs, diff, diff_err, q_base=DIFF_Q_BASE):
    """
    Run src.saxs_metrics.peak_analysis on a difference curve and judge whether
    the result is a real, localized peak or just baseline-window noise (the
    argmax sitting right at the search-window edge).

    Returns dict(peak, peak_z, resolvable, ok).
    """
    peak = src.saxs_metrics.peak_analysis(qs, diff, q_base=q_base)
    if peak is None:
        return {"peak": None, "peak_z": np.nan, "resolvable": False, "ok": False}

    i_peak = np.argmin(np.abs(qs - peak["q_peak"]))
    at_edge = np.isclose(
        peak["q_peak"], qs[(qs >= 0.9) & (qs <= 1.9)][[0, -1]]
    ).any()
    resolvable = peak["area"] > 0 and not at_edge
    if not resolvable:
        return {"peak": peak, "peak_z": np.nan, "resolvable": False, "ok": False}

    peak_z = peak["height"] / diff_err[i_peak] if diff_err[i_peak] > 0 else np.nan
    ok = np.isfinite(peak_z) and peak_z >= SIGMA_WARN
    return {"peak": peak, "peak_z": peak_z, "resolvable": True, "ok": ok}


def replicate_envelope(qs, Is_p, Is_b, scale=1.0):
    """
    Per-replicate difference curves: each individual protein frame minus the
    (optionally scaled) blank mean, and the protein mean minus each individual
    scaled blank frame. This approximates "what the diff curve would look
    like had we only collected one replicate", to show how much run-to-run
    noise alone could produce.

    Returns diff_reps (n_reps x n_q).
    """
    I_p_mean = Is_p.mean(axis=0)
    I_b_mean = Is_b.mean(axis=0)
    reps_from_protein = Is_p - scale * I_b_mean       # each protein frame vs blank mean
    reps_from_blank = I_p_mean - scale * Is_b         # protein mean vs each blank frame
    return np.vstack([reps_from_protein, reps_from_blank])


MIN_SHARED_SPOTS = 2  # need at least this many shared, well-populated spots
                       # to trust a position-matched result at all


def position_matched_diff(protein_name, blank_name, norm="linear_saxs"):
    """
    Alternative to the pooled-mean subtraction in analyze_pair(): match
    protein and blank frames by PHYSICAL stage spot (stage_positions.SPOT_OF)
    and only compare frames taken at spots BOTH groups actually visited with
    enough repeats (>= MIN_FRAMES_PER_SPOT each), combining the per-spot
    differences with equal weight. This removes any bias from the two
    groups sampling the mixing-channel positions in different proportions
    (position_locking.py showed the intensity cycle is genuinely
    position-locked, so such a mismatch would leak straight into a pooled
    difference).

    Returns None (and prints exactly why) if coverage is too poor to trust
    -- fewer than MIN_SHARED_SPOTS shared, well-populated spots -- rather
    than silently falling back to a nearest-spot match. This is expected to
    happen for at least one real pair in this dataset (Y2F+kpi/aa vs kpi/aa
    sit at measurably different x, ~0.28 mm apart -- see stage_positions.py
    -- so they may share zero spots even though both use "position 1-3/5"
    labels).
    """
    protein = cg.load_group(protein_name)
    blank = cg.load_group(blank_name)

    cg.NORM_METHODS[norm](protein)
    cg.NORM_METHODS[norm](blank)
    keep_p = cg.kept_mask(protein)
    keep_b = cg.kept_mask(blank)

    spots_p = spot_groups(protein, keep_p, list(protein.scan_ids))
    spots_b = spot_groups(blank, keep_b, list(blank.scan_ids))

    qs = protein.qs
    Is_p_kept = protein.Is[keep_p]
    Is_b_kept = blank.Is[keep_b]

    all_spots_p = set(spots_p)
    all_spots_b = set(spots_b)
    shared = sorted(
        s for s in (all_spots_p & all_spots_b)
        if len(spots_p[s]) >= MIN_FRAMES_PER_SPOT and len(spots_b[s]) >= MIN_FRAMES_PER_SPOT
    )

    n_frames_p_total, n_frames_b_total = int(keep_p.sum()), int(keep_b.sum())
    n_frames_p_used = sum(len(spots_p[s]) for s in shared)
    n_frames_b_used = sum(len(spots_b[s]) for s in shared)
    discarded_p = sorted(all_spots_p - set(shared))
    discarded_b = sorted(all_spots_b - set(shared))

    print(f"\n  Position-matched coverage ({protein_name} vs {blank_name}):")
    print(f"    protein spots: {sorted(all_spots_p)} (kept frames: {n_frames_p_total})")
    print(f"    blank spots:   {sorted(all_spots_b)} (kept frames: {n_frames_b_total})")
    print(f"    shared, well-populated spots (n>={MIN_FRAMES_PER_SPOT} both sides): {shared}")
    if discarded_p:
        print(f"    discarded (protein-only or under-populated) spots: {discarded_p}")
    if discarded_b:
        print(f"    discarded (blank-only or under-populated) spots: {discarded_b}")

    if len(shared) < MIN_SHARED_SPOTS:
        print(f"    SKIPPED: only {len(shared)} shared spot(s) with enough frames "
              f"(need >= {MIN_SHARED_SPOTS}). Position-matched subtraction is not "
              f"trustworthy for this pair -- most likely a systematic x/y offset "
              f"between the two conditions' stage coordinates (see stage_positions.py). "
              f"NOT forcing a nearest-spot match.")
        return None

    frac_p = n_frames_p_used / n_frames_p_total if n_frames_p_total else 0.0
    frac_b = n_frames_b_used / n_frames_b_total if n_frames_b_total else 0.0
    print(f"    frames used: protein {n_frames_p_used}/{n_frames_p_total} ({frac_p:.0%}), "
          f"blank {n_frames_b_used}/{n_frames_b_total} ({frac_b:.0%})")

    per_spot = {}
    for s in shared:
        Ip, Ib = Is_p_kept[spots_p[s]], Is_b_kept[spots_b[s]]
        I_p_mean, I_p_std = Ip.mean(axis=0), Ip.std(axis=0)
        I_b_mean, I_b_std = Ib.mean(axis=0), Ib.std(axis=0)
        sem_p = I_p_std / np.sqrt(len(Ip))
        sem_b = I_b_std / np.sqrt(len(Ib))
        diff = I_p_mean - I_b_mean
        diff_err = np.sqrt(sem_p ** 2 + sem_b ** 2)
        per_spot[s] = {"I_p_mean": I_p_mean, "I_b_mean": I_b_mean,
                       "diff": diff, "diff_err": diff_err, "n_p": len(Ip), "n_b": len(Ib)}

    # Equal-weight combination across shared spots -- this is what removes
    # the position-mix bias: each spot counts once, regardless of how many
    # frames either group happened to collect there.
    S = len(shared)
    I_p_matched = np.mean([per_spot[s]["I_p_mean"] for s in shared], axis=0)
    I_b_matched = np.mean([per_spot[s]["I_b_mean"] for s in shared], axis=0)
    matched_diff = I_p_matched - I_b_matched
    matched_diff_err = np.sqrt(np.sum([per_spot[s]["diff_err"] ** 2 for s in shared], axis=0)) / S

    scale = fit_scale_factor(qs, I_p_matched, I_b_matched)
    matched_diff_scaled = I_p_matched - scale * I_b_matched

    res = resolve_peak(qs, matched_diff, matched_diff_err)
    res_sc = resolve_peak(qs, matched_diff_scaled, matched_diff_err)

    # Per-spot peak diagnostics -- the key check: does a "peak" appear at
    # every shared spot (consistent, likely real) or only one/a few
    # (geometry-dependent, likely a stage-position artifact)?
    per_spot_peaks = {s: resolve_peak(qs, per_spot[s]["diff"], per_spot[s]["diff_err"]) for s in shared}

    return {
        "qs": qs, "shared_spots": shared, "per_spot": per_spot, "per_spot_peaks": per_spot_peaks,
        "I_p_matched": I_p_matched, "I_b_matched": I_b_matched,
        "matched_diff": matched_diff, "matched_diff_err": matched_diff_err,
        "matched_diff_scaled": matched_diff_scaled, "scale": scale,
        "res": res, "res_sc": res_sc,
        "n_frames_p_used": n_frames_p_used, "n_frames_p_total": n_frames_p_total,
        "n_frames_b_used": n_frames_b_used, "n_frames_b_total": n_frames_b_total,
    }


def print_matched_summary(matched):
    if matched is None:
        return
    res, res_sc = matched["res"], matched["res_sc"]

    def fmt(label, r):
        if not r["resolvable"]:
            print(f"  Matched peak significance ({label}): no resolvable peak above baseline {_flag(False)}")
        else:
            print(f"  Matched peak significance ({label}): {r['peak_z']:.1f} sigma {_flag(r['ok'])}")

    fmt("raw", res)
    fmt(f"scale-corrected, s={matched['scale']:.3f}", res_sc)

    print("  Per-spot peak breakdown (does the peak show up everywhere, or just one spot?):")
    for s in matched["shared_spots"]:
        r = matched["per_spot_peaks"][s]
        x, y = sp.SPOT_COORD.get(s, (np.nan, np.nan))
        if r["resolvable"]:
            print(f"    spot {s} (x={x:.3f}, y={y:.3f}): {r['peak_z']:.1f} sigma {_flag(r['ok'])}")
        else:
            print(f"    spot {s} (x={x:.3f}, y={y:.3f}): no resolvable peak")


def plot_matched(title, protein_name, blank_name, norm, matched, out_dir="figures"):
    if matched is None:
        return
    os.makedirs(out_dir, exist_ok=True)
    qs = matched["qs"]

    fig, (ax_curve, ax_diff) = plt.subplots(2, 1, figsize=(7, 8), sharex=True,
                                             gridspec_kw={"height_ratios": [2, 1.3]})
    ax_curve.plot(qs, matched["I_p_matched"], color="tab:red",
                  label=f"{protein_name} (position-matched, S={len(matched['shared_spots'])})")
    ax_curve.plot(qs, matched["I_b_matched"], color="tab:blue",
                  label=f"{blank_name} (position-matched)")
    ax_curve.set_xscale("log")
    ax_curve.set_yscale("log")
    ax_curve.set_ylabel("Inten. (norm.)")
    ax_curve.set_xlim(src.det.FULL_QRANGE)
    ax_curve.set_title(f"{title} -- position-matched")
    ax_curve.legend(fontsize=8)

    ax_diff.axhline(0, color="k", linewidth=0.8)
    cmap = plt.colormaps.get_cmap("tab10")
    for i, s in enumerate(matched["shared_spots"]):
        ax_diff.plot(qs, matched["per_spot"][s]["diff"], color=cmap(i), alpha=0.5, linewidth=1,
                     label=f"spot {s} (n_p={matched['per_spot'][s]['n_p']}, "
                           f"n_b={matched['per_spot'][s]['n_b']})")
    ax_diff.plot(qs, matched["matched_diff"], color="k", linewidth=2, label="equal-weight matched diff")
    ax_diff.fill_between(qs, matched["matched_diff"] - matched["matched_diff_err"],
                          matched["matched_diff"] + matched["matched_diff_err"],
                          color="k", alpha=0.15, linewidth=0)
    ax_diff.set_yscale("symlog", linthresh=1e-5)
    ax_diff.set_xlabel("q [1/A]")
    ax_diff.set_ylabel("protein - blank\n(position-matched)")
    ax_diff.legend(fontsize=7, ncol=2)
    fig.tight_layout()

    out_path = os.path.join(out_dir, title.lower().replace(" ", "_") + f"__{norm}__matched.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


def analyze_pair(title, protein_name, blank_name, norm="linear_saxs", out_dir="figures"):
    os.makedirs(out_dir, exist_ok=True)

    print(f"\nPair: {title}  [normalization: {norm}]")
    protein = cg.load_group(protein_name)
    blank = cg.load_group(blank_name)

    qs_p, Is_p, nk_p, nt_p = cg.kept_frames(protein, norm=norm)
    qs_b, Is_b, nk_b, nt_b = cg.kept_frames(blank, norm=norm)
    assert np.allclose(qs_p, qs_b), "protein/blank groups are on different q grids"
    qs = qs_p

    I_p, std_p = Is_p.mean(axis=0), Is_p.std(axis=0)
    I_b, std_b = Is_b.mean(axis=0), Is_b.std(axis=0)

    diff = I_p - I_b
    sem_p = std_p / np.sqrt(max(nk_p, 1))
    sem_b = std_b / np.sqrt(max(nk_b, 1))
    diff_err = np.sqrt(sem_p ** 2 + sem_b ** 2)
    rel_diff = diff / np.where(I_b != 0, I_b, np.nan)

    # --- 1. Run-to-run stability -------------------------------------------------
    cv_p = cg.coeff_of_variation(qs, I_p, std_p)
    cv_b = cg.coeff_of_variation(qs, I_b, std_b)
    cv_ratio = cv_p / cv_b if cv_b > 0 else np.nan
    cv_ok = cv_ratio < CV_RATIO_WARN

    # --- 2. Scaling factor test ---------------------------------------------------
    scale = fit_scale_factor(qs, I_p, I_b)
    scale_ok = abs(scale - 1.0) <= SCALE_DEVIATION_WARN

    # --- 3 & 4. Peak significance + replicate envelope (raw diff) -----------------
    diff_reps = replicate_envelope(qs, Is_p, Is_b)
    rep_lo = np.percentile(diff_reps, 5, axis=0)
    rep_hi = np.percentile(diff_reps, 95, axis=0)

    res = resolve_peak(qs, diff, diff_err)
    peak, peak_z, resolvable, peak_ok = res["peak"], res["peak_z"], res["resolvable"], res["ok"]
    rep_heights = np.array([])
    if resolvable:
        rep_peaks = [src.saxs_metrics.peak_analysis(qs, row, q_base=DIFF_Q_BASE) for row in diff_reps]
        rep_heights = np.array([r["height"] for r in rep_peaks if r is not None])

    # --- 6. Scale-corrected peak significance --------------------------------------
    # Same test, but subtracting scale*blank instead of blank, so an overall
    # multiplicative offset (from #2) can't be mistaken for -- or hide -- a
    # localized structural peak.
    diff_scaled = I_p - scale * I_b
    diff_scaled_err = np.sqrt(sem_p ** 2 + (scale * sem_b) ** 2)
    diff_reps_scaled = replicate_envelope(qs, Is_p, Is_b, scale=scale)

    res_sc = resolve_peak(qs, diff_scaled, diff_scaled_err)
    peak_sc, peak_z_sc, resolvable_sc, peak_ok_sc = (
        res_sc["peak"], res_sc["peak_z"], res_sc["resolvable"], res_sc["ok"]
    )
    rep_heights_sc = np.array([])
    if resolvable_sc:
        rep_peaks_sc = [src.saxs_metrics.peak_analysis(qs, row, q_base=DIFF_Q_BASE) for row in diff_reps_scaled]
        rep_heights_sc = np.array([r["height"] for r in rep_peaks_sc if r is not None])

    # --- Plot ----------------------------------------------------------------
    fig, (ax_curve, ax_diff, ax_rel) = plt.subplots(
        3, 1, figsize=(7, 11), sharex=True,
        gridspec_kw={"height_ratios": [2, 1.3, 1]},
    )

    ax_curve.plot(qs, I_p, color="tab:red", label=f"{protein_name} (n={nk_p})")
    ax_curve.fill_between(qs, I_p - std_p, I_p + std_p, color="tab:red", alpha=0.2, linewidth=0)
    ax_curve.plot(qs, I_b, color="tab:blue", label=f"{blank_name} (n={nk_b})")
    ax_curve.fill_between(qs, I_b - std_b, I_b + std_b, color="tab:blue", alpha=0.2, linewidth=0)
    ax_curve.axvspan(*Q_NEUTRAL, color="gray", alpha=0.15, label="scaling-fit window")
    ax_curve.set_xscale("log")
    ax_curve.set_yscale("log")
    ax_curve.set_ylabel("Inten. (norm.)")
    ax_curve.set_xlim(src.det.FULL_QRANGE)
    ax_curve.set_title(title)
    ax_curve.legend(fontsize=8)

    ax_diff.axhline(0, color="k", linewidth=0.8)
    ax_diff.fill_between(qs, rep_lo, rep_hi, color="gray", alpha=0.3, linewidth=0,
                          label="replicate envelope (5-95%)")
    ax_diff.plot(qs, diff, color="tab:green", label="mean diff")
    ax_diff.fill_between(qs, diff - diff_err, diff + diff_err, color="tab:green", alpha=0.3, linewidth=0)
    ax_diff.plot(qs, diff_scaled, color="tab:orange", linestyle="--",
                 label=f"scale-corrected diff (s={scale:.2f})")
    # symlog: the low-q Guinier region dwarfs the high-q WAXS-peak difference on a
    # linear scale, so use a threshold near the WAXS-scale noise floor to keep
    # both regions legible.
    ax_diff.set_yscale("symlog", linthresh=1e-5)
    ax_diff.set_ylabel("protein - blank")
    ax_diff.legend(fontsize=8)

    ax_rel.axhline(0, color="k", linewidth=0.8)
    ax_rel.plot(qs, rel_diff, color="tab:purple")
    ax_rel.set_yscale("symlog", linthresh=0.01)
    ax_rel.set_xlabel("q [1/A]")
    ax_rel.set_ylabel("(protein-blank)\n/ blank")
    fig.tight_layout()

    out_path = os.path.join(out_dir, title.lower().replace(" ", "_") + f"__{norm}.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")

    # --- Diagnostic summary ----------------------------------------------------
    print(f"  Protein CV: {cv_p:.0%}, Buffer CV: {cv_b:.0%}, "
          f"Ratio: {cv_ratio:.2f} {_flag(cv_ok)}"
          + ("" if cv_ok else " (protein noisier than buffer -- diff may reflect instability, not structure)"))

    print(f"  Scaling factor (fit over q={Q_NEUTRAL[0]}-{Q_NEUTRAL[1]} 1/A): "
          f"{scale:.3f} {_flag(scale_ok)}"
          + ("" if scale_ok else " (deviates >5% from 1.0 -- check normalization/scaling)"))

    if not resolvable:
        edge_q = peak["q_peak"] if peak is not None else None
        print(f"  Peak significance (raw diff): no resolvable peak above baseline"
              + (f" (argmax at window edge q={edge_q:.3f}, likely noise)" if edge_q is not None else "")
              + f" {_flag(False)}")
    else:
        print(f"  Peak significance (raw diff): {peak_z:.1f} sigma {_flag(peak_ok)}"
              + ("" if peak_ok else f" (below {SIGMA_WARN:.0f} sigma threshold)"))
        if rep_heights.size:
            print(f"  Replicate variance in diff curve: peak height range "
                  f"[{rep_heights.min():.3g}, {rep_heights.max():.3g}] "
                  f"(mean={rep_heights.mean():.3g}, n={rep_heights.size} single-replicate diffs)")
        else:
            print("  Replicate variance in diff curve: no replicate produced a resolvable peak")

    if not resolvable_sc:
        edge_q_sc = peak_sc["q_peak"] if peak_sc is not None else None
        print(f"  Peak significance (scale-corrected, s={scale:.3f}): no resolvable peak above baseline"
              + (f" (argmax at window edge q={edge_q_sc:.3f}, likely noise)" if edge_q_sc is not None else "")
              + f" {_flag(False)}")
    else:
        print(f"  Peak significance (scale-corrected, s={scale:.3f}): {peak_z_sc:.1f} sigma {_flag(peak_ok_sc)}"
              + ("" if peak_ok_sc else f" (below {SIGMA_WARN:.0f} sigma threshold)"))
        if rep_heights_sc.size:
            print(f"  Replicate variance in scale-corrected diff: peak height range "
                  f"[{rep_heights_sc.min():.3g}, {rep_heights_sc.max():.3g}] "
                  f"(mean={rep_heights_sc.mean():.3g}, n={rep_heights_sc.size} single-replicate diffs)")

    # --- Drift check: is the blank_name condition's own early-vs-late split
    # producing a peak of comparable size, from nothing but time/session
    # drift? If so, the protein-vs-blank result above may be partly or
    # entirely a time-gap artifact rather than real protein signal.
    drift_frac = np.nan
    drift = cg.drift_check(blank_name, norm=norm)
    if drift is None:
        print(f"  Drift check: skipped ('{blank_name}' has too few scans to split, "
              f"need >= {cg.MIN_SCANS_FOR_DRIFT_CHECK})")
    elif drift["peak"] is None or drift["peak"]["area"] <= 0:
        print(f"  Drift check: no drift-only peak found in '{blank_name}' "
              f"(early n={drift['n_early']}, late n={drift['n_late']}) {_flag(True)}")
    else:
        drift_height = drift["peak"]["height"]
        # Compare against the scale-corrected peak, not the raw one -- the
        # drift check's own diff is itself scale-corrected (early vs
        # scale*late), so this is the apples-to-apples comparison.
        ref_height = peak_sc["height"] if resolvable_sc else np.nan
        drift_ok = True
        if resolvable_sc and np.isfinite(ref_height) and ref_height > 0:
            drift_frac = drift_height / ref_height
            drift_ok = drift_frac < 0.5
        print(f"  Drift check ('{blank_name}' early n={drift['n_early']} vs late "
              f"n={drift['n_late']}, scale={drift['scale']:.3f}): "
              f"drift-only peak height={drift_height:.3g}"
              + (f", = {drift_frac:.0%} of the protein-vs-blank peak height ({ref_height:.3g})"
                 if np.isfinite(drift_frac) else "")
              + f" {_flag(drift_ok)}"
              + ("" if drift_ok else
                 " (drift alone produces a peak >=50% the size of the claimed effect -- "
                 "treat the result above with caution)"))

    # --- 8. Position-matched subtraction ---------------------------------------
    matched = position_matched_diff(protein_name, blank_name, norm=norm)
    print_matched_summary(matched)
    plot_matched(title, protein_name, blank_name, norm, matched, out_dir=out_dir)

    matched_peak_z = matched["res"]["peak_z"] if matched and matched["res"]["resolvable"] else np.nan
    matched_peak_z_sc = matched["res_sc"]["peak_z"] if matched and matched["res_sc"]["resolvable"] else np.nan
    matched_n_shared = len(matched["shared_spots"]) if matched else 0

    return {
        "norm": norm, "cv_p": cv_p, "cv_b": cv_b, "cv_ratio": cv_ratio,
        "scale": scale, "peak_z": peak_z, "resolvable": resolvable,
        "peak_z_scaled": peak_z_sc, "resolvable_scaled": resolvable_sc,
        "drift_frac": drift_frac,
        "matched_peak_z": matched_peak_z, "matched_peak_z_sc": matched_peak_z_sc,
        "matched_n_shared": matched_n_shared,
    }


# Normalization strategies to compare -- see compare_groups.NORM_METHODS for
# what each one actually does to the data.
NORM_METHODS_TO_COMPARE = ["linear_saxs", "neutral", "i0"]


def compare_normalizations(title, protein_name, blank_name):
    """Run analyze_pair under each normalization method and print a comparison table."""
    results = [
        analyze_pair(title, protein_name, blank_name, norm=norm)
        for norm in NORM_METHODS_TO_COMPARE
    ]

    print(f"\n  == {title}: normalization comparison ==")
    header = (f"  {'method':<14}{'scale':>10}{'CV ratio':>12}"
              f"{'peak sig.':>14}{'scale-corr. peak':>20}{'drift %':>10}"
              f"{'matched peak':>15}{'(S spots)':>11}")
    print(header)
    for r in results:
        peak_str = f"{r['peak_z']:.1f} sigma" if r["resolvable"] and np.isfinite(r["peak_z"]) else "n/a"
        peak_sc_str = (f"{r['peak_z_scaled']:.1f} sigma"
                       if r["resolvable_scaled"] and np.isfinite(r["peak_z_scaled"]) else "n/a")
        drift_str = f"{r['drift_frac']:.0%}" if np.isfinite(r["drift_frac"]) else "n/a"
        matched_str = f"{r['matched_peak_z_sc']:.1f} sigma" if np.isfinite(r["matched_peak_z_sc"]) else "n/a"
        s_str = f"S={r['matched_n_shared']}" if r["matched_n_shared"] else "skipped"
        print(f"  {r['norm']:<14}{r['scale']:>10.3f}{r['cv_ratio']:>12.2f}"
              f"{peak_str:>14}{peak_sc_str:>20}{drift_str:>10}"
              f"{matched_str:>15}{s_str:>11}")
    return results


if __name__ == "__main__":
    for title, protein_name, blank_name in PAIRS:
        compare_normalizations(title, protein_name, blank_name)
    plt.show()

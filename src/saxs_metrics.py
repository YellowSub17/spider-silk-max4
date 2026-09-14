"""
Small set of standard SAXS/WAXS reduction metrics, applied to the merged
I(q) curves produced by Scan/ComboScan (src.det.FULL_QRANGE).

These curves are detector-merged, scaled and normalized (not absolute
intensity), so treat the outputs as *relative*/comparative quantities
(good for tracking change between conditions or positions), not as
calibrated physical values.
"""

import numpy as np
from scipy import stats


def guinier_fit(qs, I, q_min=0, q_max=0.02, min_points=5, refine=True,
                 refine_min_points=10):
    """
    Guinier analysis: ln(I) = ln(I0) - (Rg^2/3) q^2, fit over q_min < q < q_max,
    then (if refine) narrow once more to q < 1.3/Rg (the standard Guinier
    validity window) -- unless that window is *wider* than [q_min, q_max]
    (refine only tightens, never widens, an explicitly requested q-range),
    or unless the narrower window would leave fewer than `refine_min_points`
    points, in which case the initial (broader) fit is kept instead.

    That last guard matters: when the initial fit's Rg comes out large
    (e.g. from a shared low-q upturn/background rather than real particle
    size), 1.3/Rg can land just above the detector's minimum measurable q,
    so "refining" collapses onto only a handful of points hugging the qmin
    edge. A straight line through ~5-10 points near one edge of the data
    will almost always look deceptively linear (high r^2) regardless of
    whether it's a real Guinier region, and can swing Rg by 2x+ from the
    broader fit. Requiring more points before trusting the refinement
    avoids reporting that unstable, edge-dominated number.

    q_min lets you skip a low-q upturn (e.g. from aggregation or a shared
    background) and fit a higher-q "Guinier-like" region instead, at the
    cost of it no longer being a rigorous Guinier fit if that region isn't
    truly linear in ln(I) vs q^2.

    Returns dict(Rg, I0, r) or None if the fit isn't usable (not enough
    points, non-linear/positive slope, etc).
    """
    def _fit(mask):
        if mask.sum() < min_points:
            return None
        x = qs[mask] ** 2
        y = np.log(I[mask])
        if not np.all(np.isfinite(y)):
            return None
        slope, intercept, r, p, se = stats.linregress(x, y)
        if slope >= 0:  # unphysical: I(q) must decay across the fit window
            return None
        rg = np.sqrt(-3 * slope)
        i0 = np.exp(intercept)
        return rg, i0, r

    mask = (qs > q_min) & (qs <= q_max) & (I > 0)
    fit = _fit(mask)
    if fit is None:
        return None
    rg, i0, r = fit
    q_fit_max = q_max
    n_points = int(mask.sum())

    if refine and 1.3 / rg < q_max:
        mask2 = (qs > q_min) & (qs <= 1.3 / rg) & (I > 0)
        if mask2.sum() >= refine_min_points:
            fit2 = _fit(mask2)
            if fit2 is not None:
                rg, i0, r = fit2
                q_fit_max = 1.3 / fit[0]  # window was set from the *initial* rg
                n_points = int(mask2.sum())

    return {"Rg": rg, "I0": i0, "r": r, "q_fit_max": q_fit_max, "n_points": n_points}


def porod_exponent(qs, I, q_range=(0.05, 0.3)):
    """
    Local power-law exponent: I(q) ~ q^-n over q_range, fit as the slope of
    log(I) vs log(q). n=4 is ideal (sharp interface) Porod behavior; lower
    values indicate mass-fractal/aggregated structure, higher values
    (with a negative deviation) indicate diffuse interfaces.

    Returns dict(n, r) or None.
    """
    mask = (qs >= q_range[0]) & (qs <= q_range[1]) & (I > 0)
    if mask.sum() < 5:
        return None
    slope, intercept, r, p, se = stats.linregress(np.log(qs[mask]), np.log(I[mask]))
    return {"n": -slope, "r": r}


def peak_analysis(qs, I, q_window=(0.9, 1.9), q_base=((0.03, 0.08), (1.9, 2.1))):
    """
    Characterize the silk beta-sheet WAXS peak(s) in q_window: subtracts a
    linear baseline anchored at the mean intensity in the two q_base windows,
    then reports the max peak height/position and the integrated (trapezoid)
    area above baseline.

    Returns dict(q_peak, height, area) or None.
    """
    mask = (qs >= q_window[0]) & (qs <= q_window[1])
    if mask.sum() < 3:
        return None

    def _base_level(rng):
        m = (qs >= rng[0]) & (qs <= rng[1])
        return np.mean(I[m]) if np.any(m) else np.nan

    b_lo = _base_level(q_base[0])
    b_hi = _base_level(q_base[1])
    if not (np.isfinite(b_lo) and np.isfinite(b_hi)):
        return None

    q_win = qs[mask]
    I_win = I[mask]
    baseline = np.interp(q_win, [q_win[0], q_win[-1]], [b_lo, b_hi])
    I_sub = I_win - baseline

    i_max = np.argmax(I_sub)
    area = np.trapezoid(np.clip(I_sub, 0, None), q_win)

    return {"q_peak": q_win[i_max], "height": I_sub[i_max], "area": area}


def kratky(qs, I):
    """Kratky representation: I(q)*q^2 vs q. A compact/globular particle's
    Kratky plot rises then decays back toward zero at high q; an extended,
    disordered, or aggregated/fractal one plateaus or keeps rising."""
    return qs ** 2 * I


def dimensionless_kratky(qs, I, rg, i0):
    """
    Normalized ("Receveur-Carrey") Kratky plot: (q*Rg)^2 * I/I0 vs q*Rg.
    A well-folded globular particle peaks near (x, y) = (1.73, 1.1) (the
    q*Rg=sqrt(3) point); deviation from that peak position/height, or a
    plateau/rise at high q*Rg, indicates non-globular (extended/disordered/
    aggregated) structure. Requires a valid Guinier fit (Rg, I0 > 0).
    """
    x = qs * rg
    y = (x ** 2) * (I / i0)
    return x, y


def all_metrics(qs, I):
    """Compute every metric for one I(q) curve; missing/failed fits are None."""
    return {
        "guinier": guinier_fit(qs, I),
        "porod": porod_exponent(qs, I),
        "peak": peak_analysis(qs, I),
    }

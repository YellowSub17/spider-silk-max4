"""
Authoritative physical sample-stage position for every scan, built from
scans_clean.csv's `sams4_x`/`sams4_y` columns rather than the `description`
("position N") label.

Why this exists: the "position N" label is assigned independently *within*
each sample's own run -- it is NOT a shared physical coordinate system.
Confirmed directly in this data: kpi/aa has only 3 labelled positions
(1/2/3) at y = 26.999 / 27.999 / 32.500 mm, while the non-acid flow groups
have 5 labelled positions (1-5) at y = 26.975 / 27.476 / 27.976 / 29.977 /
32.477 mm. kpi/aa's label 2 (y=27.999) is physically the same spot as the
5-position groups' label 3 (y=27.976); kpi/aa's label 3 (y=32.500) is their
label 5 (y=32.477). kpi/aa also sits at a different x (~117.04 vs ~116.76
mm) than every other flow group. Treating "position 2" as the same physical
location across different samples is therefore wrong in general.

There's also a per-cycle anomaly: within every group that has one, "position
1" sits at a measurably different x than positions 2-5 (a consistent ~19.6
um lateral offset) and carries much larger y jitter -- i.e. it is one odd
stop in an otherwise tight, repeated cycle, most likely a settling/backlash
artifact of the stage's move sequence rather than a different intended spot.

This module clusters every scan's (x, y) coordinate into global "spot" IDs
shared across all samples (same physical location -> same ID, regardless of
what label it carried), and flags samples whose scans do not cleanly cluster
(i.e. a single label spans multiple physical spots, or vice versa) --
Y2F+kpi/aa is the clearest example: its label "position 2" alone spans 4.94
mm of y and 0.30 mm of x, clearly straddling more than one physical stop.

Public API:
  SPOT_OF[scanID]        -> global spot id (int)
  SPOT_COORD[spot_id]    -> (x, y) mean coordinate of that spot
  spots_of_group(name)   -> {spot_id: [scanIDs]} for a compare_groups.GROUPS name
  UNRELIABLE_SAMPLES      -> {raw Sample string: reason} for samples whose
                              labelled positions don't cleanly correspond to
                              tight physical clusters

Run directly (`python stage_positions.py`) to print the sample/label/spot
diagnostic table described in the project brief.
"""

import csv
from collections import defaultdict

import numpy as np

import compare_groups as cg

SCANS_CLEAN_CSV = cg.SCANS_CLEAN_CSV

# Coordinate tolerance for clustering into the same physical spot (mm on
# each axis). The real position-to-position spacing in this data is >=0.5
# mm, and repeat measurements at the same spot agree to sub-micron, so this
# has wide margin on both sides.
SPOT_TOLERANCE_MM = 0.05

# A (sample, label) group is flagged unreliable if its scans' coordinate
# spread exceeds this many multiples of SPOT_TOLERANCE_MM on either axis --
# i.e. that label does not correspond to one tight physical spot.
UNRELIABLE_SPREAD_FACTOR = 10.0


def load_scan_coords(path=SCANS_CLEAN_CSV):
    """scanID -> (x, y, raw Sample, raw label)."""
    coords = {}
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            sid = row["scanID"].strip()
            if not sid:
                continue
            coords[int(sid)] = (
                float(row["sams4_x"]), float(row["sams4_y"]),
                row["Sample"].strip(), row["description"].strip(),
            )
    return coords


SCAN_COORDS = load_scan_coords()


def _cluster(points, tol=SPOT_TOLERANCE_MM):
    """
    Greedy single-pass clustering: for each point (in input order), join the
    nearest existing cluster if within tol on both axes (Chebyshev distance
    to the running centroid), else start a new cluster. Good enough here
    because real spots are separated by >=0.5 mm (10x tol) -- there's no
    ambiguous middle ground to get wrong.

    points: list of (key, x, y). Returns {key: cluster_index}, and a list of
    (mean_x, mean_y) per cluster index.
    """
    centroids = []  # [x_sum, y_sum, n]
    assign = {}
    for key, x, y in points:
        best = None
        for i, (xs, ys, n) in enumerate(centroids):
            cx, cy = xs / n, ys / n
            if abs(x - cx) <= tol and abs(y - cy) <= tol:
                best = i
                break
        if best is None:
            centroids.append([x, y, 1])
            best = len(centroids) - 1
        else:
            centroids[best][0] += x
            centroids[best][1] += y
            centroids[best][2] += 1
        assign[key] = best
    coords = [(xs / n, ys / n) for xs, ys, n in centroids]
    return assign, coords


def _build_spots():
    points = [(sid, x, y) for sid, (x, y, _, _) in SCAN_COORDS.items()]
    raw_assign, raw_coords = _cluster(points)

    # Relabel spot IDs in a stable, readable order: sorted by y then x
    # (roughly upstream-to-downstream along the channel).
    order = sorted(range(len(raw_coords)), key=lambda i: (raw_coords[i][1], raw_coords[i][0]))
    relabel = {old: new for new, old in enumerate(order)}

    spot_of = {sid: relabel[c] for sid, c in raw_assign.items()}
    spot_coord = {relabel[c]: raw_coords[c] for c in range(len(raw_coords))}
    return spot_of, spot_coord


SPOT_OF, SPOT_COORD = _build_spots()


def spots_of_group(name):
    """{spot_id: [scanIDs]} for every scan in a compare_groups.GROUPS entry."""
    by_spot = defaultdict(list)
    for sid in cg.GROUPS[name]:
        if sid in SPOT_OF:
            by_spot[SPOT_OF[sid]].append(sid)
    return dict(by_spot)


# Only samples actually used downstream (compare_groups.GROUPS, via RAW_NAME)
# -- the raw data also contains blank-Sample exploratory/alignment scans from
# before any condition was defined, which would otherwise flood the
# diagnostic table and unreliable-sample list with noise irrelevant to the
# analysis pipeline.
RELEVANT_SAMPLES = set(cg.RAW_NAME.values())


def _find_unreliable_samples():
    """
    For every (raw Sample, label) combination among RELEVANT_SAMPLES, check
    whether its scans fall into a single tight physical cluster. Flags any
    combination whose coordinate spread (max-min on either axis) exceeds
    UNRELIABLE_SPREAD_FACTOR * SPOT_TOLERANCE_MM -- i.e. the label does not
    correspond to one physical spot -- and separately flags any (raw Sample,
    label) that maps to >1 distinct spot id (a related but not identical
    check: a label can span multiple spot ids while technically each spot's
    scans are tight, if the scans alternate between two nearby spots).
    """
    by_sample_label = defaultdict(list)
    for sid, (x, y, sample, label) in SCAN_COORDS.items():
        if sample in RELEVANT_SAMPLES:
            by_sample_label[(sample, label)].append(sid)

    unreliable = {}
    for (sample, label), sids in by_sample_label.items():
        xs = [SCAN_COORDS[s][0] for s in sids]
        ys = [SCAN_COORDS[s][1] for s in sids]
        x_spread = max(xs) - min(xs)
        y_spread = max(ys) - min(ys)
        spot_ids = {SPOT_OF[s] for s in sids}
        threshold = UNRELIABLE_SPREAD_FACTOR * SPOT_TOLERANCE_MM
        if x_spread > threshold or y_spread > threshold or len(spot_ids) > 1:
            reasons = []
            if x_spread > threshold:
                reasons.append(f"x spread {x_spread:.3f}mm")
            if y_spread > threshold:
                reasons.append(f"y spread {y_spread:.3f}mm")
            if len(spot_ids) > 1:
                reasons.append(f"spans {len(spot_ids)} distinct spots {sorted(spot_ids)}")
            unreliable[(sample, label)] = ", ".join(reasons)
    return unreliable


UNRELIABLE_SAMPLES = _find_unreliable_samples()


def print_diagnostic_table():
    """sample | label | spot id | n | x mean+-spread | y mean+-spread | y offset from upstream spot.

    Restricted to RELEVANT_SAMPLES (see its docstring) -- the raw scan log
    also contains blank-Sample exploratory/alignment scans that would
    otherwise dominate this table with noise unrelated to any condition
    used downstream.
    """
    by_key = defaultdict(list)  # (sample, label, spot_id) -> [scanIDs]
    for sid, (x, y, sample, label) in SCAN_COORDS.items():
        if sample in RELEVANT_SAMPLES:
            by_key[(sample, label, SPOT_OF[sid])].append(sid)

    # per-sample most-upstream spot (smallest mean y among that sample's spots)
    sample_spots = defaultdict(set)
    for (sample, label, spot_id) in by_key:
        sample_spots[sample].add(spot_id)
    upstream_y = {
        sample: min(SPOT_COORD[s][1] for s in spots)
        for sample, spots in sample_spots.items()
    }

    header = (f"{'sample':<22}{'label':<14}{'spot':>5}{'n':>5}"
              f"{'x mean':>11}{'x spread':>10}{'y mean':>11}{'y spread':>10}{'y offset':>10}  flags")
    print(header)
    print("-" * len(header))
    for (sample, label, spot_id), sids in sorted(by_key.items(), key=lambda kv: (kv[0][0], kv[0][2])):
        xs = [SCAN_COORDS[s][0] for s in sids]
        ys = [SCAN_COORDS[s][1] for s in sids]
        x_mean, y_mean = np.mean(xs), np.mean(ys)
        x_spread, y_spread = max(xs) - min(xs), max(ys) - min(ys)
        y_offset = y_mean - upstream_y[sample]
        flag = " <-- UNRELIABLE" if (sample, label) in UNRELIABLE_SAMPLES else ""
        print(f"{sample:<22}{label:<14}{spot_id:>5}{len(sids):>5}"
              f"{x_mean:>11.4f}{x_spread:>10.4f}{y_mean:>11.4f}{y_spread:>10.4f}{y_offset:>10.4f}{flag}")

    print(f"\n{len(UNRELIABLE_SAMPLES)} (sample, label) combinations flagged unreliable:")
    for (sample, label), reason in UNRELIABLE_SAMPLES.items():
        print(f"  {sample!r} label={label!r}: {reason}")


if __name__ == "__main__":
    print(f"{len(SPOT_COORD)} global spots found across {len(SCAN_COORDS)} scans "
          f"(tolerance={SPOT_TOLERANCE_MM} mm)\n")
    print_diagnostic_table()

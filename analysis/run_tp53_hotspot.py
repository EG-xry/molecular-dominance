#!/usr/bin/env python3
"""TP53 hotspot-anchored calibration sensitivity (revision of 8 November 2026).

    python3 run_tp53_hotspot.py

The main TP53 analysis (run_tp53.py) shows that synonymous/nonsense reference
classes do not supply a common dominant-negative endpoint for the two
wild-type-background screens: nonsense variants in the scanned region truncate
before the tetramerisation domain and score near neutral there.  This script
asks the complementary question.  If the most frequent TP53 missense hotspot
mutations, fixed in advance from cancer mutation frequencies rather than from
these screens' scores, are taken as the common dominant-negative endpoint, how
far do the two mixed-state slopes converge?

Each axis is anchored separately: z = (score - median_synonymous) /
(median_hotspot - median_synonymous), so synonymous variants sit at 0 and the
hotspot class at 1 on the alone-state axis and on each mixed-state axis.  The
through-origin slope is then computed on the same common set of missense
variants used in run_tp53.py.

Caveats carried into every output:
  * circularity: the dominant-negative behaviour of hotspot mutants was itself
    demonstrated partly with these screens (Giacomelli 2018; Boettcher 2019);
  * hotspot mutants differ in dominant-negative potency, so a class median is
    an anchor choice, not a measured maximal endpoint;
  * two-point anchoring does not make nonlinear response functions equivalent.

Writes tp53_hotspot.csv (per-variant anchor scores, leave-one-out slopes) and
tp53_hotspot.json (slopes, ratios, bootstrap intervals, settings).
"""
import json
from datetime import date

import numpy as np
import pandas as pd

from run_tp53 import ALONE_STATE_SCREENS, SCREENS, cohorts, load

# Pre-specified anchor set: the eight most frequent TP53 missense hotspot
# substitutions in human cancers (see HOTSPOT_SOURCE).  Exploratory choice used consistently when computing
# any anchored slope; not selected using the screens analysed here.
HOTSPOTS = ["R175H", "G245S", "R248Q", "R248W", "R249S", "R273C", "R273H", "R282W"]
HOTSPOT_SOURCE = ("most frequent TP53 missense hotspot substitutions in human cancers; "
                  "citation recorded in the manuscript and S1 Text")
NBOOT = 5000
SEED = 20261108
ALONE_FILE = "mavedb_00000068-b-1.csv"


def anchors(frame, hotspots):
    """Synonymous and hotspot medians of one deposited score axis."""
    syn = frame.loc[frame.syn, "score"]
    hot = frame.loc[frame.key.isin(hotspots), "score"]
    if len(syn) < 2 or len(hot) < 2:
        raise ValueError("insufficient anchors")
    return float(syn.median()), float(hot.median()), int(len(syn)), int(len(hot))


def anchored(values, syn_med, hot_med):
    return (np.asarray(values, dtype=float) - syn_med) / (hot_med - syn_med)


def slope(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float(np.sum(x * y) / np.sum(x * x))


def anchored_slopes(alone, mixed, common, hotspots):
    """Through-origin slopes after anchoring both axes, on the common set."""
    sx, hx, _, _ = anchors(alone, hotspots)
    out = {}
    for name, frame in mixed.items():
        sy, hy, _, _ = anchors(frame, hotspots)
        j = alone.merge(frame, on="key", suffixes=("_x", "_y"), validate="one_to_one")
        j = j[j.key.isin(common) & ~j.syn_x]
        out[name] = slope(anchored(j.score_x, sx, hx), anchored(j.score_y, sy, hy))
    return out


def main():
    alone = load(ALONE_FILE)
    mixed = {name: load(fn, col) for name, (fn, col) in SCREENS.items()}
    alone_state = {name: load(fn, col) for name, (fn, col) in ALONE_STATE_SCREENS.items()}
    common, _ = cohorts(alone, mixed, alone_state)

    names = list(mixed)
    present = {h: [h in set(alone.key)] + [h in set(mixed[n].key) for n in names]
               for h in HOTSPOTS}
    missing = [h for h, flags in present.items() if not all(flags)]
    if missing:
        raise ValueError(f"hotspots missing from an axis: {missing}")

    # Raw slopes on the common set, recomputed here for a self-contained ratio.
    raw = {}
    for name, frame in mixed.items():
        j = alone.merge(frame, on="key", suffixes=("_x", "_y"), validate="one_to_one")
        j = j[j.key.isin(common) & ~j.syn_x]
        raw[name] = slope(j.score_x, j.score_y)

    anchor_rows = []
    for label, frame in [("Giacomelli, p53-null + nutlin-3 (alone axis)", alone), *mixed.items()]:
        s, h, ns, nh = anchors(frame, HOTSPOTS)
        anchor_rows.append({"axis": label, "synonymous_median": s, "hotspot_median": h,
                            "n_synonymous": ns, "n_hotspot": nh, "hotspot_minus_synonymous": h - s})

    primary = anchored_slopes(alone, mixed, common, HOTSPOTS)
    loo = {}
    for h in HOTSPOTS:
        loo[h] = anchored_slopes(alone, mixed, common, [k for k in HOTSPOTS if k != h])

    # Response-only anchoring (alone axis left at its deposited origin and
    # scale) as a secondary check on the dependence on the predictor anchor.
    response_only = {}
    for name, frame in mixed.items():
        sy, hy, _, _ = anchors(frame, HOTSPOTS)
        j = alone.merge(frame, on="key", suffixes=("_x", "_y"), validate="one_to_one")
        j = j[j.key.isin(common) & ~j.syn_x]
        response_only[name] = slope(j.score_x, anchored(j.score_y, sy, hy))

    # Joint residue bootstrap: resample residues once per draw (multinomial
    # counts) and re-estimate every anchor and slope from the same resampled
    # residues.  Counts weight both the anchor medians and the slope sums.
    rng = np.random.default_rng(SEED)
    frames = {"alone": alone, **mixed}
    residues = np.array(sorted(set().union(*(set(f.pos) for f in frames.values()))))
    max_pos = int(residues.max())

    def axis_parts(f):
        syn = f[f.syn]
        hot = f[f.key.isin(HOTSPOTS)]
        return (syn.score.to_numpy(), syn.pos.to_numpy(), hot.score.to_numpy(), hot.pos.to_numpy())

    parts = {k: axis_parts(f) for k, f in frames.items()}
    paired = {}
    for n in names:
        j = alone.merge(mixed[n], on="key", suffixes=("_x", "_y"), validate="one_to_one")
        j = j[j.key.isin(common) & ~j.syn_x]
        paired[n] = (j.score_x.to_numpy(), j.score_y.to_numpy(), j.pos_x.to_numpy())

    def wmedian(values, weights):
        rep = np.repeat(values, weights)
        if rep.size < 2:
            raise ValueError("insufficient anchors")
        return float(np.median(rep))

    boot = {n: [] for n in names}
    boot_ratio = []
    invalid = 0
    for _ in range(NBOOT):
        cnt = np.zeros(max_pos + 1, dtype=int)
        np.add.at(cnt, residues[rng.integers(0, len(residues), len(residues))], 1)
        try:
            med = {}
            for k, (ss, sp, hs, hp) in parts.items():
                med[k] = (wmedian(ss, cnt[sp]), wmedian(hs, cnt[hp]))
            vals = {}
            for n in names:
                x, y, pos = paired[n]
                w = cnt[pos]
                xa = (x - med["alone"][0]) / (med["alone"][1] - med["alone"][0])
                ya = (y - med[n][0]) / (med[n][1] - med[n][0])
                vals[n] = float(np.sum(w * xa * ya) / np.sum(w * xa * xa))
        except ValueError:
            invalid += 1
            continue
        for n in names:
            boot[n].append(vals[n])
        boot_ratio.append(vals[names[1]] / vals[names[0]])

    def ci(v):
        return [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))]

    result = {
        "generated": date.today().isoformat(),
        "common_missense_n": len(common),
        "hotspot_set": HOTSPOTS,
        "hotspot_set_rationale": HOTSPOT_SOURCE,
        "anchor_definition": "z = (score - median_synonymous) / (median_hotspot - median_synonymous) on every axis",
        "anchors": anchor_rows,
        "raw_slopes": raw,
        "raw_ratio_second_over_first": raw[names[1]] / raw[names[0]],
        "anchored_slopes": primary,
        "anchored_ratio_second_over_first": primary[names[1]] / primary[names[0]],
        "anchored_bootstrap_ci": {n: ci(boot[n]) for n in names},
        "anchored_ratio_bootstrap_ci": ci(boot_ratio),
        "leave_one_out": loo,
        "leave_one_out_range": {n: [min(v[n] for v in loo.values()), max(v[n] for v in loo.values())]
                                for n in names},
        "response_only_anchored_slopes": response_only,
        "bootstrap": {"n_boot": NBOOT, "seed": SEED, "invalid_draws": invalid,
                      "scope": "residues resampled jointly across all three axes; anchors re-estimated in every draw; conditional on deposited scores"},
        "caveats": [
            "Hotspot dominant-negative behaviour was demonstrated partly with these screens, so the anchor is not independent of the data.",
            "Hotspot mutants differ in dominant-negative potency; the class median is an anchor choice, not a measured maximal endpoint.",
            "Two-point anchoring does not establish equivalent nonlinear response functions or molecular activity scales.",
        ],
    }
    with open("tp53_hotspot.json", "w") as fh:
        json.dump(result, fh, indent=2, allow_nan=False)

    rows = []
    for h in HOTSPOTS:
        row = {"hotspot": h, "alone_score": float(alone.loc[alone.key == h, "score"].iloc[0])}
        for n in names:
            row[f"score::{n}"] = float(mixed[n].loc[mixed[n].key == h, "score"].iloc[0])
            row[f"loo_slope::{n}"] = loo[h][n]
        rows.append(row)
    pd.DataFrame(rows).to_csv("tp53_hotspot.csv", index=False)

    print(f"common missense set: {len(common)}")
    print(pd.DataFrame(anchor_rows).to_string(index=False))
    for n in names:
        print(f"{n}: raw {raw[n]:.3f}; anchored {primary[n]:.3f} "
              f"(95% CI {ci(boot[n])[0]:.3f}-{ci(boot[n])[1]:.3f}); "
              f"LOO {result['leave_one_out_range'][n][0]:.3f}-{result['leave_one_out_range'][n][1]:.3f}; "
              f"response-only {response_only[n]:.3f}")
    print(f"ratio raw {result['raw_ratio_second_over_first']:.3f}; anchored "
          f"{result['anchored_ratio_second_over_first']:.3f} "
          f"(95% CI {ci(boot_ratio)[0]:.3f}-{ci(boot_ratio)[1]:.3f}); invalid draws {invalid}")
    print("wrote tp53_hotspot.csv/json")


if __name__ == "__main__":
    main()

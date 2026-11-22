#!/usr/bin/env python3
"""Alternative uncertainty intervals for the reported slopes (revision of 8 November 2026).

    python3 run_uncertainty.py

The main analyses report residue-cluster percentile intervals, which treat the
panel's residues as a sample from residues ascertained in the same way and
condition on the published variant means. This script reports, side by side
for every slope given with an interval:

  * residue: the same residue-cluster percentile bootstrap;
  * two-stage: residues are resampled, then each sampled variant's alone and
    mixed values are perturbed by Gaussian noise with its deposited SEs;
  * wild cluster: a Rademacher wild bootstrap of the through-origin regression
    with residue clusters, suited to panels with few clusters.

Two-stage intervals are unavailable where the deposits give no per-variant
SEs (the TP53 screens). Perturbing the alone-state values also shows how
measurement error attenuates a through-origin slope: the two-stage bootstrap
median can fall below the point estimate. None of these intervals covers
shared reference, normalization, or batch uncertainty. Writes uncertainty.json.
"""
import json
from datetime import date

import numpy as np
import pandas as pd

import dominance as D
from run_tp53 import ALONE_STATE_SCREENS, SCREENS, cohorts, load as load_tp53

NBOOT = 10000
SEED = 20261110


def kv7_panel(path):
    d = pd.read_csv(path)
    d = d[(d.sx > 0) & (d.sy > 0)].copy()
    d = d[[D.classify_variant(v)[3] == "missense" for v in d.variant]]
    return d.x.to_numpy(float), d.y.to_numpy(float), d.sx.to_numpy(float), d.sy.to_numpy(float), d.pos.to_numpy()


def kcnq2(atag, mtag, col, wt):
    """MaveDB 00000664 pairs with SEs and residues (panel or controls)."""
    from run_kv7_rank import load
    a, m = load(atag), load(mtag)
    key = "hgvs_nt" if a.hgvs_nt.notna().mean() > 0.9 else "hgvs_pro"
    other = "hgvs_pro" if key == "hgvs_nt" else "hgvs_nt"
    j = a.merge(m[[key, col, "SEM", "score"]], on=key, suffixes=("_a", "_m"), validate="one_to_one")
    j = j.dropna(subset=[f"{col}_a", f"{col}_m", "SEM_a", "SEM_m"])
    if wt is None:
        x, y = j[f"{col}_a"].to_numpy(float), j[f"{col}_m"].to_numpy(float)
        sx, sy = j.SEM_a.to_numpy(float), j.SEM_m.to_numpy(float)
    else:
        x, y = wt - j[f"{col}_a"].to_numpy(float), wt - j[f"{col}_m"].to_numpy(float)
        sx = np.abs((j.SEM_a * j[f"{col}_a"] / j.score_a).to_numpy(float))
        sy = np.abs((j.SEM_m * j[f"{col}_m"] / j.score_m).to_numpy(float))
    pro = j[other] if key == "hgvs_nt" else j[key]
    nt = j[key] if key == "hgvs_nt" else j[other]
    codon = np.ceil(nt.astype(str).str.extract(r"c\.(\d+)")[0].astype(float) / 3)
    res = pro.astype(str).str.extract(r"(\d+)")[0].astype(float).fillna(codon).to_numpy()
    syn = pro.astype(str).str.contains("=", regex=False).to_numpy()
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(sx) & np.isfinite(sy) & np.isfinite(res)
    return x[keep], y[keep], sx[keep], sy[keep], res[keep], syn[keep]


def intervals(x, y, clusters, sx=None, sy=None, rng=None):
    x, y = np.asarray(x, float), np.asarray(y, float)
    b = float(np.sum(x * y) / np.sum(x * x))
    labels, inv = np.unique(clusters, return_inverse=True)
    G = len(labels)
    sxy = np.bincount(inv, weights=x * y, minlength=G)
    sxx = np.bincount(inv, weights=x * x, minlength=G)
    counts = rng.multinomial(G, np.full(G, 1.0 / G), size=NBOOT)
    res = (counts @ sxy) / (counts @ sxx)
    out = {"b": b, "n": int(len(x)), "clusters": int(G),
           "residue_ci": [float(np.quantile(res, 0.025)), float(np.quantile(res, 0.975))]}
    if sx is not None:
        groups = [np.where(inv == g)[0] for g in range(G)]
        two = np.empty(NBOOT)
        for k in range(NBOOT):
            idx = np.concatenate([groups[g] for g in rng.integers(0, G, G)])
            xs = x[idx] + sx[idx] * rng.standard_normal(len(idx))
            ys = y[idx] + sy[idx] * rng.standard_normal(len(idx))
            two[k] = np.sum(xs * ys) / np.sum(xs * xs)
        out["two_stage_ci"] = [float(np.quantile(two, 0.025)), float(np.quantile(two, 0.975))]
        out["two_stage_median"] = float(np.median(two))
    u = y - b * x
    score = np.bincount(inv, weights=x * u, minlength=G)
    v = rng.choice([-1.0, 1.0], size=(NBOOT, G))
    wild = b + (v @ score) / np.sum(x * x)
    out["wild_cluster_ci"] = [float(np.quantile(wild, 0.025)), float(np.quantile(wild, 0.975))]
    return out


def main():
    rng = np.random.default_rng(SEED)
    rows = {}
    for name, path in (("KCNQ4 activation-voltage panel", "data/kv7/fit_KCNQ4_primary.csv"),
                       ("KCNQ2 activation-voltage panel", "data/kv7/fit_KCNQ2_primary.csv"),
                       ("KCNQ2/KCNQ3 activation-voltage panel", "data/kv7/fit_KCNQ2_KCNQ3_primary.csv")):
        x, y, sx, sy, pos = kv7_panel(path)
        rows[name] = intervals(x, y, pos, sx, sy, rng)
    endpoints = {"activation-voltage shift": (("a-2", "a-5"), ("a-14", "a-17"), "Shifted V1/2 (mV)", None),
                 "current at +40 mV": (("a-1", "a-4"), ("a-13", "a-16"), "Normalized current density at +40 mV", 1.0),
                 "current at -40 mV": (("a-3", "a-6"), ("a-15", "a-18"), "Normalized current density at -40 mV", 1.0)}
    for label, (ptags, ctags, col, wt) in endpoints.items():
        x, y, sx, sy, res, _ = kcnq2(*ptags, col, wt)
        rows[f"KCNQ2 disease panel, {label}"] = intervals(x, y, res, sx, sy, rng)
        x, y, sx, sy, res, syn = kcnq2(*ctags, col, wt)
        rows[f"KCNQ2 benign controls, {label}"] = intervals(x, y, res, sx, sy, rng)
        rows[f"KCNQ2 synonymous controls, {label}"] = intervals(x[syn], y[syn], res[syn], sx[syn], sy[syn], rng)
    peak = pd.read_csv("scn5a_peak.csv")
    ww = float(pd.read_csv("scn5a_controls.csv").set_index("variant").loc["WT+WT", "het"])
    mis = peak[peak.variant != "W822X"]
    scn5a_scale = ww / 100
    x = (ww - scn5a_scale * mis.homo).to_numpy(float)
    y = (ww - mis.het).to_numpy(float)
    pos = mis.variant.str.extract(r"(\d+)")[0].astype(int).to_numpy()
    scn5a_label = "SCN5A matched-reference scaling slope (s=2.184, W=218.4)"
    rows[scn5a_label] = intervals(x, y, pos, scn5a_scale * mis.homo_se.to_numpy(float),
                                 mis.het_se.to_numpy(float), rng)
    rows[scn5a_label]["scale_assumptions"] = {
        "reference": ww, "combined_scale": scn5a_scale, "p": 1-100/ww,
        "scope": "Fixed reference-conversion/dose factor paired with reference; reference and scale uncertainty omitted."}
    alone = load_tp53("mavedb_00000068-b-1.csv")
    mixed = {n: load_tp53(f, c) for n, (f, c) in SCREENS.items()}
    other = {n: load_tp53(f, c) for n, (f, c) in ALONE_STATE_SCREENS.items()}
    common, _ = cohorts(alone, mixed, other)
    for n, frame in {**mixed, **other}.items():
        j = alone.merge(frame, on="key", suffixes=("_x", "_y"), validate="one_to_one").dropna(subset=["score_x", "score_y"])
        j = j[~j.syn_x]
        if n in mixed:
            j = j[j.key.isin(common)]
        rows[f"TP53 {n}"] = intervals(j.score_x, j.score_y, j.pos_x.to_numpy(), rng=rng)
    out = {"generated": date.today().isoformat(), "n_boot": NBOOT, "seed": SEED,
           "superpopulation": ("residues of the same protein region, ascertained and selected as in each panel; "
                               "intervals condition on the published variant means except that the two-stage "
                               "version also propagates deposited per-variant SEs"),
           "wild_cluster": "Rademacher weights per residue cluster on through-origin residuals",
           "limits": ["No interval covers shared reference, normalization, or batch uncertainty.",
                      "TP53 deposits give no per-variant SEs for the mixed screens or the common alone axis."],
           "slopes": rows}
    with open("uncertainty.json", "w") as fh:
        json.dump(out, fh, indent=2)
    for k, r in rows.items():
        two = (f"two-stage {r['two_stage_ci'][0]:.3f}-{r['two_stage_ci'][1]:.3f} (median {r['two_stage_median']:.3f}); "
               if "two_stage_ci" in r else "two-stage n/a; ")
        print(f"{k}: b {r['b']:.3f}, n {r['n']}, clusters {r['clusters']}; residue {r['residue_ci'][0]:.3f}-"
              f"{r['residue_ci'][1]:.3f}; {two}wild {r['wild_cluster_ci'][0]:.3f}-{r['wild_cluster_ci'][1]:.3f}")
    print("wrote uncertainty.json")


if __name__ == "__main__":
    main()

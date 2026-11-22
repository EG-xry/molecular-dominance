

import os as _os


def _data(vendored, fallback):
    """Prefer the copy vendored under dominance/data/, so a checkout of this
    directory alone reproduces everything; fall back to the sibling tree for
    anyone working in the full repository."""
    here = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "data", vendored)
    return here if _os.path.exists(here) else fallback

#!/usr/bin/env python3
"""Measure the Kv7 rank-order claim instead of asserting it.

    python3 run_kv7_rank.py

Describe agreement between per-variant ratios on different KCNQ2 endpoints.
The filter requires |x| to exceed the control RMS on each endpoint. The same
filter leaves only 5-6 controls per endpoint pair, so a matched control test
is inconclusive. Unfiltered controls provide a descriptive sensitivity
analysis; their near-zero denominators make them unsuitable for establishing
that ratio ranks are immune to reference errors or score offsets.

Controls are the deposition's 57 OddsPath benign variants, including 43
synonymous variants. Clinical benignity does not establish zero assay effect.
Intervals resample residues, conditional on published scores and fixed RMS
thresholds; they do not resample shared wild-type references or batches.

Writes kv7_rank.csv.
"""
import json
from datetime import date

import numpy as np
import pandas as pd
from scipy import stats

CACHE = _data("kv7", "../heterogeneity/newdata")
rng = np.random.default_rng(20260928)
NBOOT = 20000

PANEL = {"dV1/2": ("a-2", "a-5", "Shifted V1/2 (mV)", None),
         "+40 mV": ("a-1", "a-4", "Normalized current density at +40 mV", 1.0),
         "-40 mV": ("a-3", "a-6", "Normalized current density at -40 mV", 1.0)}
CONTROL = {"dV1/2": ("a-14", "a-17", "Shifted V1/2 (mV)", None),
           "+40 mV": ("a-13", "a-16", "Normalized current density at +40 mV", 1.0),
           "-40 mV": ("a-15", "a-18", "Normalized current density at -40 mV", 1.0)}
# A transfer ratio y/x is meaningless where x is within the noise, and the
# axes have very different scales, so the threshold is set per axis from that
# axis's own benign controls: the alone-state effect must exceed the RMS
# alone-state value of the control panel. This is a screening scale, not an
# identified measurement-error variance.


def load(tag):
    d = pd.read_csv(f"{CACHE}/mavedb_00000664-{tag}.csv")
    d.columns = [c.strip() for c in d.columns]
    return d


def paired_measurements(atag, mtag, col, wt):
    a, m = load(atag), load(mtag)
    key = "hgvs_nt" if a.hgvs_nt.notna().mean() > 0.9 else "hgvs_pro"
    other = "hgvs_pro" if key == "hgvs_nt" else "hgvs_nt"
    if a[key].isna().any() or m[key].isna().any():
        raise ValueError(f"missing {key} identifier in {atag}/{mtag}")
    j = (a[[key, other, col]].merge(m[[key, col]], on=key, suffixes=("_x", "_y"),
                                   validate="one_to_one")
         .rename(columns={f"{col}_x": "x", f"{col}_y": "y"}).dropna(subset=["x", "y"]))
    if wt is not None:
        j["x"], j["y"] = wt - j.x, wt - j.y
    pro = j[other] if key == "hgvs_nt" else j[key]
    nt = j[key] if key == "hgvs_nt" else j[other]
    codon = (nt.astype(str).str.extract(r"c\.(\d+)")[0].astype(float) / 3).apply(np.ceil)
    j["res"] = pro.astype(str).str.extract(r"(\d+)")[0].astype(float).fillna(codon)
    # the identifier must be the UNIQUE key for this score set: hgvs_nt where
    # the deposition carries it, hgvs_pro otherwise. Preferring hgvs_pro would
    # give all 43 synonymous controls the same id ("p.=").
    j["id"] = j[key].astype(str)
    return j.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["res", "x", "y"])[["id", "res", "x", "y"]]


def ratios(atag, mtag, col, wt):
    j = paired_measurements(atag, mtag, col, wt)
    j["ratio"] = j.y / j.x.replace(0, np.nan)
    return j.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["res", "ratio"])[["id", "res", "ratio", "x"]]


# NO cross-convention normalisation. The comparison is always endpoint against
# endpoint WITHIN one set -- controls against controls, panel against panel --
# so each set keeps its own identifier. An earlier version normalised both to a
# protein-change string, which collapsed all 43 synonymous controls onto "p.="
# and formed a 43x43 cartesian product, inflating 5 real pairs into 129.


def spearman_ci(a, b, res):
    r = float(stats.spearmanr(a, b).statistic)
    groups = [np.array(g) for g in pd.Series(range(len(a))).groupby(res.values).apply(list)]
    out = []
    for _ in range(NBOOT):
        s = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
        if len(set(a.values[s])) < 3 or len(set(b.values[s])) < 3:
            continue
        out.append(stats.spearmanr(a.values[s], b.values[s]).statistic)
    out = np.array([v for v in out if np.isfinite(v)])
    return r, float(np.quantile(out, 0.025)), float(np.quantile(out, 0.975))


def main():
    NOISE = {k: float(np.sqrt(np.mean(paired_measurements(*v).x ** 2)))
             for k, v in CONTROL.items()}
    print("screening scale per axis, from benign controls (RMS alone-state):")
    for k, v in NOISE.items():
        print(f"  {k:8s} {v:.3f}")
    print()

    MIN_N = 10
    rows = []
    # The |x| >= noise filter is right for a panel, where a ratio with a near-zero
    # denominator is unstable. On a putatively neutral control panel it removes almost
    # everything by construction, so the control comparison is reported BOTH ways:
    # unfiltered, where all 57 controls pair, and filtered, where 5-6 remain.
    for label, src in (("disease panel", PANEL), ("benign controls", CONTROL),
                       ("benign controls, unfiltered", CONTROL)):
        R = {}
        for k, v in src.items():
            t = ratios(*v)
            R[k] = t if label.endswith("unfiltered") else t[t.x.abs() >= NOISE[k]]
        pairs = [("dV1/2", "+40 mV"), ("dV1/2", "-40 mV"), ("+40 mV", "-40 mV")]
        for p, q in pairs:
            j = R[p].merge(R[q], on="id", suffixes=("_p", "_q"), validate="one_to_one")
            if len(j) < MIN_N:
                rows.append({"set": label, "endpoints": f"{p} vs {q}", "n": len(j),
                             "n_residues": int(j.res_p.nunique()),
                             "spearman": None, "ci_lo": None, "ci_hi": None,
                             "excludes_zero": None,
                             "interpretation": "inconclusive: fewer than 10 matched variants"})
                continue
            r, lo, hi = spearman_ci(j.ratio_p, j.ratio_q, j.res_p)
            rows.append({"set": label, "endpoints": f"{p} vs {q}", "n": len(j),
                         "n_residues": int(j.res_p.nunique()),
                         "spearman": round(r, 4), "ci_lo": round(lo, 4),
                         "ci_hi": round(hi, 4),
                         "excludes_zero": not (lo <= 0 <= hi),
                         "interpretation": ("descriptive sensitivity; unstable small denominators"
                                            if label.endswith("unfiltered") else
                                            "within-panel descriptive association")})

    df = pd.DataFrame(rows)
    df.to_csv("kv7_rank.csv", index=False)
    json.dump({"generated": date.today().isoformat(), "urn": "urn:mavedb:00000664",
               "control_rms_per_axis": NOISE, "n_boot": NBOOT, "seed": 20260928,
               "control_basis": "OddsPath benign variants; zero-effect assumption unverified",
               "bootstrap_scope": "residues conditional on published scores and fixed RMS cutoffs",
               "rows": rows},
              open("kv7_rank.json", "w"), indent=2)
    print("Per-variant transfer-ratio agreement between endpoints, KCNQ2\n")
    print(df.to_string(index=False))
    d = df[df.set == "disease panel"].dropna(subset=["spearman"])
    print(f"\ndisease panel: rho {d.spearman.min():+.3f} to {d.spearman.max():+.3f}; "
          f"{int(d.excludes_zero.sum())}/{len(d)} intervals exclude zero")
    print("The matched filtered-control comparison is inconclusive (6/5/5 variants).")
    print("Unfiltered controls are a sensitivity analysis and do not establish rank immunity.")
    print("\nwrote kv7_rank.csv, kv7_rank.json")


if __name__ == "__main__":
    main()

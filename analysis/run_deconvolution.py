

import os as _os


def _data(vendored, fallback):
    """Prefer the copy vendored under dominance/data/, so a checkout of this
    directory alone reproduces everything; fall back to the sibling tree for
    anyone working in the full repository."""
    here = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "data", vendored)
    return here if _os.path.exists(here) else fallback

#!/usr/bin/env python3
"""Sensitivity of a conditional moment correction to the assay origin.

    python3 run_deconvolution.py

For through-origin slopes b = E[xy]/E[x^2], the candidate correction is

    d = (E[x_p y_p] - E[x_c y_c]) / (E[x_p^2] - E[x_c^2]).

These are raw second moments, not variances. Interpreting d as biological
transfer requires negligible true control effects, transferable nuisance
moments, no signal-nuisance cross terms, and a meaningful common effect scale.
Under additive common offsets, the published-origin version violates the
zero-cross-term condition unless signal means or offsets satisfy additional
restrictions. Control centring removes the offsets; the intercept version
targets a different, centred estimand. The 57 OddsPath benign controls do not
establish neutral effects or transferable nuisance moments. The control
slope is a descriptive statistic, not an identified shared-noise fraction.

Each specification is bootstrapped by residue within the disease and control
panels independently. These intervals condition on the published variant
means and cannot recover uncertainty in shared references or experimental
batches. Similar replicate counts do not establish equal nuisance variances.

Writes deconvolution.csv and deconvolution.json.
"""
import json
from datetime import date

import numpy as np
import pandas as pd

CACHE = _data("kv7", "../heterogeneity/newdata")
rng = np.random.default_rng(20260928)
NBOOT = 20000

# axis: (panel alone, panel mixed, control alone, control mixed, column, wt level)
AXES = {
    "dV1/2 (mV)":        ("a-2", "a-5", "a-14", "a-17", "Shifted V1/2 (mV)", None),
    "current at +40 mV": ("a-1", "a-4", "a-13", "a-16", "Normalized current density at +40 mV", 1.0),
    "current at -40 mV": ("a-3", "a-6", "a-15", "a-18", "Normalized current density at -40 mV", 1.0),
}


def load(tag):
    d = pd.read_csv(f"{CACHE}/mavedb_00000664-{tag}.csv")
    d.columns = [c.strip() for c in d.columns]
    return d


def pair(atag, mtag, col, wt, with_se=False):
    """Join a variant's alone and mixed measurements.

    The deposition keys the disease panels on hgvs_pro (their hgvs_nt is
    empty) and the control panels on hgvs_nt (their hgvs_pro is "p.=" for all
    43 synonymous controls), so the key is chosen per panel rather than fixed.
    """
    a, m = load(atag), load(mtag)
    key = "hgvs_nt" if a.hgvs_nt.notna().mean() > 0.9 else "hgvs_pro"
    other = "hgvs_pro" if key == "hgvs_nt" else "hgvs_nt"
    if a[key].isna().any() or m[key].isna().any():
        raise ValueError(f"missing {key} identifier in {atag}/{mtag}")
    extra = ["SEM", "score"] if with_se else []
    j = (a[[key, other, col] + extra]
         .merge(m[[key, col] + extra], on=key, suffixes=("_x", "_y"), validate="one_to_one")
         .rename(columns={f"{col}_x": "x", f"{col}_y": "y"})
         .dropna(subset=["x", "y"]))
    if with_se:
        j["sx"], j["sy"] = j.SEM_x, j.SEM_y
        if wt is not None:
            j["sx"] = (j.SEM_x * j.x / j.score_x).abs()
            j["sy"] = (j.SEM_y * j.y / j.score_y).abs()
    if wt is not None:
        j["x"], j["y"] = wt - j.x, wt - j.y
    nt = j[key] if key == "hgvs_nt" else j[other]
    pro = j[other] if key == "hgvs_nt" else j[key]
    codon = (nt.astype(str).str.extract(r"c\.(\d+)")[0]
             .astype(float) / 3).apply(np.ceil)
    j["res"] = (pro.astype(str).str.extract(r"(\d+)")[0]
                .astype(float).fillna(codon))
    return j.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["res", "x", "y"] + (["sx", "sy"] if with_se else []))[
            ["x", "y", "res"] + (["sx", "sy"] if with_se else [])].reset_index(drop=True)


def slope(x, y):
    return float(np.sum(x * y) / np.sum(x * x))


def correct(xp, yp, xc, yc, origin="as-published"):
    """Candidate moment correction under an explicitly specified origin.

      as-published      slopes and second moments about zero
      control-centred   both axes shifted by the control means
      intercept         slopes from OLS with a constant, central variances
    """
    if origin not in ("as-published", "control-centred", "intercept"):
        raise ValueError(f"unknown origin specification: {origin}")
    xp, yp, xc, yc = [np.asarray(v, dtype=float) for v in (xp, yp, xc, yc)]
    if origin == "control-centred":
        mx, my = xc.mean(), yc.mean()
        xp, yp, xc, yc = xp - mx, yp - my, xc - mx, yc - my
    if origin == "intercept":
        vp, vc = float(np.var(xp, ddof=1)), float(np.var(xc, ddof=1))
        cp = float(np.cov(xp, yp, ddof=1)[0, 1])
        cc = float(np.cov(xc, yc, ddof=1)[0, 1])
    else:
        vp, vc = float(np.mean(xp ** 2)), float(np.mean(xc ** 2))
        cp, cc = float(np.mean(xp * yp)), float(np.mean(xc * yc))
    bp, bc = cp / vp if vp > 0 else np.nan, cc / vc if vc > 0 else np.nan
    d = (cp - cc) / (vp - vc) if vp > vc else np.nan
    return d, bp, bc, vp / vc if vc > 0 else np.nan


def clusters(df):
    """Positional index arrays, one per residue."""
    pos = {k: i for i, k in enumerate(df.index)}
    return [np.array([pos[i] for i in g.index]) for _, g in df.groupby("res")]


def resample(groups):
    pick = rng.integers(0, len(groups), len(groups))
    return np.concatenate([groups[i] for i in pick])


def precision(tag):
    d = load(tag)
    return float(d.replicates.median())


ORIGINS = {"as-published": "as_published", "control-centred": "control_centred",
           "intercept": "intercept"}


def se_correct(x, y, sx):
    """EIV comparator assuming zero-mean errors and no x/y error covariance."""
    x, y, sx = [np.asarray(v, dtype=float) for v in (x, y, sx)]
    den = np.sum(x * x - sx * sx)
    return float(np.sum(x * y) / den) if den > 0 else np.nan


def percentile_summary(values):
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    return {"ci_lo": float(np.quantile(finite, 0.025)),
            "ci_hi": float(np.quantile(finite, 0.975)),
            "n_invalid": len(values) - len(finite)}


def main():
    rows, difference_rows, control_rows, se_rows = [], [], [], []
    for axis, (pa, pm, ca, cm, col, wt) in AXES.items():
        P, C = pair(pa, pm, col, wt, with_se=True), pair(ca, cm, col, wt)
        Px, Py, Cx, Cy = P.x.values, P.y.values, C.x.values, C.y.values
        d, bp, bc, mr = correct(Px, Py, Cx, Cy)
        by_origin = {o: correct(Px, Py, Cx, Cy, o)[0] for o in ORIGINS}
        mean_share = float(C.x.mean() * C.y.mean() / np.mean(Cx * Cy))
        t_x = float(C.x.mean() / (C.x.std(ddof=1) / np.sqrt(len(C))))
        t_y = float(C.y.mean() / (C.y.std(ddof=1) / np.sqrt(len(C))))
        Pg, Cg = clusters(P), clusters(C)
        draws = {o: np.empty(NBOOT) for o in ORIGINS}
        means = np.empty((NBOOT, 2))
        se_draws = np.empty((NBOOT, 2))
        for i in range(NBOOT):
            ip, ic = resample(Pg), resample(Cg)
            for origin in ORIGINS:
                draws[origin][i] = correct(Px[ip], Py[ip], Cx[ic], Cy[ic], origin)[0]
            means[i] = Cx[ic].mean(), Cy[ic].mean()
            se_draws[i, 0] = se_correct(Px[ip], Py[ip], P.sx.values[ip])
            se_draws[i, 1] = se_correct(Px[ip] - means[i, 0], Py[ip] - means[i, 1], P.sx.values[ip])
        row = {
            "axis": axis, "n_panel": len(P), "n_control": len(C),
            "b_panel_raw": round(bp, 4), "b_control": round(bc, 4),
            "second_moment_ratio_panel_over_control": round(mr, 2),
            "median_reps_panel": precision(pa), "median_reps_control": precision(ca),
            "d_origin_span": round(max(by_origin.values()) - min(by_origin.values()), 4),
            "control_mean_t_x": round(t_x, 2), "control_mean_t_y": round(t_y, 2),
            "mean_product_share_of_control_numerator": round(mean_share, 3),
        }
        for origin, key in ORIGINS.items():
            summary = percentile_summary(draws[origin])
            row[f"d_{key}"] = round(by_origin[origin], 4)
            row[f"ci_lo_{key}"] = round(summary["ci_lo"], 4)
            row[f"ci_hi_{key}"] = round(summary["ci_hi"], 4)
            row[f"n_invalid_boot_{key}"] = summary["n_invalid"]
        for first, second in (("control-centred", "as-published"),
                              ("intercept", "as-published"), ("control-centred", "intercept")):
            difference_rows.append({"axis": axis, "contrast": f"{first} minus {second}",
                "difference": by_origin[first] - by_origin[second],
                **percentile_summary(draws[first] - draws[second])})
        for arm, values, j in (("alone", Cx, 0), ("mixed", Cy, 1)):
            control_rows.append({"axis": axis, "arm": arm, "mean": float(values.mean()),
                                 **percentile_summary(means[:, j])})
        for j, origin in enumerate(("as-published", "control-centred")):
            mx, my = (0.0, 0.0) if j == 0 else (Cx.mean(), Cy.mean())
            se_rows.append({"axis": axis, "origin": origin,
                            "estimate": se_correct(Px - mx, Py - my, P.sx.values),
                            **percentile_summary(se_draws[:, j])})
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv("deconvolution.csv", index=False)
    pd.DataFrame(difference_rows).to_csv("deconvolution_origin_differences.csv", index=False)
    pd.DataFrame(control_rows).to_csv("deconvolution_control_means.csv", index=False)
    pd.DataFrame(se_rows).to_csv("deconvolution_se_comparator.csv", index=False)
    with open("deconvolution.json", "w") as f:
        json.dump({
            "generated": date.today().isoformat(), "urn": "urn:mavedb:00000664",
            "n_boot": NBOOT, "seed": 20260928,
            "control_basis": "OddsPath benign variants; zero-effect assumption unverified",
            "interpretation": "specification-dependent descriptive estimates; published-origin correction is inconsistent under general additive offsets",
            "bootstrap_scope": "independent residue resampling within panels, shared across origin specifications; control means reestimated; shared references fixed",
            "nuisance_model": {
                "equations": "x=S+c_x+e; y=T+c_y+u; controls have S=T=0; residual errors have zero means, no signal-error cross moments, and transferable second moments",
                "as_published_limit": "(E[ST]+c_y E[S]+c_x E[T])/(E[S^2]+2 c_x E[S])",
                "control_centred_target": "E[ST]/E[S^2]",
                "intercept_target": "Cov(S,T)/Var(S)",
                "caveat": "Targets require neutral controls and transferable residual nuisance moments; neither assumption is identified from these data.",
            },
            "origin_definitions": {
                "as-published": "raw cross moments about deposited origins; offset-biased under the additive-offset model",
                "control-centred": "subtract common offsets using control means; targets through-origin signal moment slope under stated assumptions",
                "intercept": "separately centred sample covariance and variance (ddof=1); targets centred signal slope",
            },
            "se_comparator_assumptions": "sum(x*y)/(sum(x*x)-sum(sx*sx)); assumes unbiased per-variant SEs, zero-mean errors and zero alone/mixed error covariance; omits shared-reference uncertainty; control-centred version reestimates control means",
            "rows": rows, "paired_origin_differences": difference_rows,
            "control_mean_intervals": control_rows, "se_comparators": se_rows}, f, indent=2)
        f.write("\n")

    print(df.to_string(index=False))
    print("\nJoint-bootstrap origin differences (conditional on shared references):")
    print(pd.DataFrame(difference_rows).to_string(index=False))
    print("\nControl means:")
    print(pd.DataFrame(control_rows).to_string(index=False))
    print("\nSE-based comparators:")
    print(pd.DataFrame(se_rows).to_string(index=False))
    print("\nThe published-origin version is not a valid biological correction under general additive offsets.")
    print("Control centring and the intercept target different signal slopes; all nuisance assumptions remain unverified.")
    print("Wrote deconvolution.csv/json and three diagnostic CSVs.")


if __name__ == "__main__":
    main()

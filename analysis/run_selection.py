#!/usr/bin/env python3
"""Selection on the alone-state axis (revision of 8 November 2026).

    python3 run_selection.py

Several paired panels were assembled by selecting variants on their measured
alone-state effect; Zheng et al. co-expressed with wild type only the KCNQ4
variants whose alone-state current or V1/2 departed from wild type.  The
KCNQ2 MaveDB panel was instead ascertained clinically, and its 301 missense
variants were measured in both arms, so it can show what such a selection rule
does to a through-origin slope.

Two calculations, on each KCNQ2 endpoint:
  1. empirical: the slope on all variants and on the subset whose alone-state
     effect differs from wild type by more than 1.96 deposited standard errors;
  2. simulation under constant transfer: observed alone-state effects are used
     as true effects, every true transfer ratio equals the full-panel slope,
     Gaussian errors with the deposited SEs are added to both arms, and the
     same selection rule is applied.  The simulated slopes separate
     selection-induced bias from changes in panel composition.

The simulation is illustrative.  Using observed effects as true effects
overstates their spread, and independent arm errors ignore shared plate or
reference effects.  Writes selection.json.
"""
import json
from datetime import date

import numpy as np
import pandas as pd

from run_deconvolution import CACHE

ENDPOINTS = {
    "current at +40 mV": ("a-1", "a-4", "Normalized current density at +40 mV", 1.0),
    "current at -40 mV": ("a-3", "a-6", "Normalized current density at -40 mV", 1.0),
    "dV1/2 (mV)": ("a-2", "a-5", "Shifted V1/2 (mV)", None),
}
Z = 1.959963984540054
NSIM = 2000
SEED = 20261108


def load(tag):
    d = pd.read_csv(f"{CACHE}/mavedb_00000664-{tag}.csv")
    d.columns = [c.strip() for c in d.columns]
    return d


def effects(atag, mtag, col, wt):
    """Paired effects and standard errors on the analysis scale.

    Normalized currents are score / wild-type mean, so their SE is the
    deposited SEM multiplied by normalized/score.  Effects are 1 - normalized
    current; voltage shifts are used as deposited."""
    a, m = load(atag), load(mtag)
    j = a.merge(m, on="hgvs_pro", suffixes=("_a", "_m"), validate="one_to_one")
    j = j.dropna(subset=[f"{col}_a", f"{col}_m", "SEM_a", "SEM_m"])
    if wt is None:
        x, y = j[f"{col}_a"].to_numpy(float), j[f"{col}_m"].to_numpy(float)
        sx, sy = j.SEM_a.to_numpy(float), j.SEM_m.to_numpy(float)
    else:
        x, y = wt - j[f"{col}_a"].to_numpy(float), wt - j[f"{col}_m"].to_numpy(float)
        sx = (j.SEM_a * j[f"{col}_a"] / j.score_a).to_numpy(float)
        sy = (j.SEM_m * j[f"{col}_m"] / j.score_m).to_numpy(float)
    keep = np.isfinite(x) & np.isfinite(y) & (sx > 0) & (sy > 0)
    return x[keep], y[keep], sx[keep], sy[keep]


def slope(x, y):
    return float(np.sum(x * y) / np.sum(x * x))


def main():
    rng = np.random.default_rng(SEED)
    out = {}
    for name, (atag, mtag, col, wt) in ENDPOINTS.items():
        x, y, sx, sy = effects(atag, mtag, col, wt)
        sel = np.abs(x) / sx > Z
        b_all, b_sel = slope(x, y), slope(x[sel], y[sel])
        d0 = b_all
        sim_all, sim_sel, n_sel = [], [], []
        for _ in range(NSIM):
            xo = x + rng.normal(0, sx)
            yo = d0 * x + rng.normal(0, sy)
            s = np.abs(xo) / sx > Z
            sim_all.append(slope(xo, yo))
            sim_sel.append(slope(xo[s], yo[s]))
            n_sel.append(int(s.sum()))
        sim_all, sim_sel = np.array(sim_all), np.array(sim_sel)
        out[name] = {
            "n_all": int(len(x)), "b_all": b_all,
            "n_selected": int(sel.sum()), "b_selected": b_sel,
            "selected_minus_all": b_sel - b_all,
            "simulation_true_ratio": d0,
            "simulation_mean_b_all": float(sim_all.mean()),
            "simulation_mean_b_selected": float(sim_sel.mean()),
            "simulation_relative_bias_all": float(sim_all.mean() / d0 - 1),
            "simulation_relative_bias_selected": float(sim_sel.mean() / d0 - 1),
            "simulation_mean_n_selected": float(np.mean(n_sel)),
            "median_abs_x_over_se": float(np.median(np.abs(x) / sx)),
        }
    # Signal-to-noise of the selected historical panels (descriptive).
    import dominance as D
    panels = {"KCNQ4 (Zheng; co-expressed after alone-state selection)": "data/kv7/fit_KCNQ4_primary.csv",
              "KCNQ2/KCNQ3 (Vanoye; heteromeric alone arm)": "data/kv7/fit_KCNQ2_KCNQ3_primary.csv"}
    snr = {}
    for label, path in panels.items():
        d = pd.read_csv(path)
        d = d[(d.sx > 0) & (d.sy > 0)].copy()
        d = d[[D.classify_variant(v)[3] == "missense" for v in d.variant]]
        r = np.abs(d.x.to_numpy(float)) / d.sx.to_numpy(float)
        x, y, sx, sy = (d[c].to_numpy(float) for c in ("x", "y", "sx", "sy"))
        d0 = slope(x, y)
        sims = []
        for _ in range(NSIM):
            xo = x + rng.normal(0, sx)
            yo = d0 * x + rng.normal(0, sy)
            sims.append(slope(xo, yo))
        snr[label] = {"n": int(len(d)), "median_abs_x_over_se": float(np.median(r)),
                      "share_abs_x_over_se_below_1_96": float(np.mean(r < Z)),
                      "simulation_relative_bias_no_further_selection": float(np.mean(sims) / d0 - 1)}
    result = {
        "generated": date.today().isoformat(),
        "historical_panels_signal_to_noise": snr,
        "panel": "MaveDB urn:mavedb:00000664 disease panels (clinically ascertained; both arms measured)",
        "selection_rule": "|alone-state effect| / deposited SE > 1.96, mimicking selection on a significant alone-state departure",
        "n_sim": NSIM, "seed": SEED,
        "simulation": "constant true transfer ratio equal to the full-panel slope; observed alone-state effects used as true effects; independent Gaussian errors with deposited SEs in both arms",
        "endpoints": out,
        "limits": ["Observed effects used as true effects overstate their spread and understate attenuation.",
                   "Independent arm errors omit shared plate, transfection, and reference effects."],
    }
    with open("selection.json", "w") as fh:
        json.dump(result, fh, indent=2)
    for name, r in out.items():
        print(f"{name}: all {r['n_all']} b={r['b_all']:.3f}; selected {r['n_selected']} "
              f"b={r['b_selected']:.3f}; simulated relative bias all {r['simulation_relative_bias_all']:+.3f}, "
              f"selected {r['simulation_relative_bias_selected']:+.3f}; median |x|/se {r['median_abs_x_over_se']:.2f}")
    for label, r in result["historical_panels_signal_to_noise"].items():
        print(f"{label}: n {r['n']}, median |x|/se {r['median_abs_x_over_se']:.2f}, "
              f"share below 1.96 {r['share_abs_x_over_se_below_1_96']:.2f}, "
              f"noise-only relative bias {r['simulation_relative_bias_no_further_selection']:+.3f}")
    print("wrote selection.json")


if __name__ == "__main__":
    main()

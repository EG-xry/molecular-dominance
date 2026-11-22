

import os as _os


def _data(vendored, fallback):
    """Prefer the copy vendored under dominance/data/, so a checkout of this
    directory alone reproduces everything; fall back to the sibling tree for
    anyone working in the full repository."""
    here = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "data", vendored)
    return here if _os.path.exists(here) else fallback

#!/usr/bin/env python3
"""Zero-effect control slopes for the Kv7 panels, with intervals.

    python3 run_kv7_controls.py

The KCNQ2 deposition (urn:mavedb:00000664) includes 57 control variants
measured on the same axes as the disease panel: a-14 alone, a-17 co-expressed
with wild type.  43 are synonymous.  If the estimator returned a biological
transfer ratio, a panel whose variants change nothing should return a slope
indistinguishable from zero.  It does not.

Synonymous variants are *putative* neutral controls, not variants demonstrated
to have exactly zero effect: synonymous changes can alter splicing, stability
or translation rate.  The argument here needs only that their true effects are
far smaller than the disease panel's, which the alone-state spread confirms.

Writes kv7_controls.csv and kv7_controls.json.
"""
import json
import sys
from datetime import date

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
import dominance as D

CACHE = _data("kv7", "../heterogeneity/newdata")
# The 57 control variants, each measured alone and co-expressed with wild type,
# on the three axes the disease panel uses.  Score-set pairs are (alone, mixed).
# Each axis must place wild type at the origin, or a through-origin slope is
# not a transfer ratio.  "Shifted V1/2" is already a shift from wild type.
# "Normalized current density" is on a wild-type = 1 scale, so the effect is
# the LOSS, 1 - value.  `wt` below is the wild-type value on the raw scale;
# None means the column is already a deviation.
AXES = {
    "dV1/2 (mV)":        ("a-14", "a-17", "Shifted V1/2 (mV)", None),
    "current at +40 mV": ("a-13", "a-16", "Normalized current density at +40 mV", 1.0),
    "current at -40 mV": ("a-15", "a-18", "Normalized current density at -40 mV", 1.0),
}
KEY = "hgvs_nt"     # unique; hgvs_pro is "p.=" for all 43 synonymous


def load(tag):
    d = pd.read_csv(f"{CACHE}/mavedb_00000664-{tag}.csv")
    d.columns = [c.strip() for c in d.columns]
    return d


rows = []
for axis, (atag, mtag, col, wt) in AXES.items():
    a, m = load(atag), load(mtag)
    j = (a[[KEY, "hgvs_pro", col, "SEM"]]
         .merge(m[[KEY, col, "SEM"]], on=KEY, suffixes=("_x", "_y"))
         .rename(columns={f"{col}_x": "x", f"{col}_y": "y",
                          "SEM_x": "sx", "SEM_y": "sy"})
         .dropna(subset=["x", "y", "sx", "sy"]))
    if wt is not None:                       # convert level to effect
        j["x"], j["y"] = wt - j.x, wt - j.y
    j["variant"] = j[KEY]
    # residue for the clustered bootstrap: from the protein change where there
    # is one, otherwise from the codon implied by the nucleotide position
    pos = j.hgvs_pro.str.extract(r"(\d+)")[0]
    pos = pos.fillna(
        (j[KEY].str.extract(r"c\.(\d+)")[0].astype(float) / 3).apply(np.ceil))
    j["pos"] = pos.astype(float)
    j = j[(j.sx > 0) & (j.sy > 0)]
    print(f"{axis}: {atag} alone, {mtag} mixed -> {len(j)} paired controls")
    syn = j[j.hgvs_pro.astype(str).str.contains("=", na=False)]
    for label, panel in (("synonymous", syn), ("all controls", j)):
        if len(panel) < 5:
            continue
        r = D.estimate(panel, units="physical", allow_uncorrelated=True,
                       min_pearson=-1.0)
        rows.append({"axis": axis, "panel": label, "n": r.n,
                     "residues": r.n_residues, "b": round(r.b, 4),
                     "ci_lo": round(r.ci[0], 4), "ci_hi": round(r.ci[1], 4),
                     "pearson_r": round(r.pearson_r, 4),
                     "axis_origin": ("already a deviation" if wt is None
                                     else f"effect = {wt:g} - value"),
                     "mean_abs_x": round(float(panel.x.abs().mean()), 3),
                     "excludes_zero": not (r.ci[0] <= 0 <= r.ci[1]),
                     "passes_gate": abs(r.b) < 0.25})

df = pd.DataFrame(rows)
df.to_csv("kv7_controls.csv", index=False)
json.dump({"generated": date.today().isoformat(),
           "source": {k: {"alone": v[0], "mixed": v[1], "column": v[2],
                          "wt_level": v[3]} for k, v in AXES.items()},
           "urn": "urn:mavedb:00000664",
           "n_boot": 20000, "seed": 20260926,
           "bootstrap": "residue-cluster percentile intervals from dominance.estimate defaults",
           "gate": "|b| < 0.25",
           "gate_note": ("analyst-chosen descriptive screen for a near-zero control slope; "
                         "not pre-registered, passed by no control slope, and not used in any inference"),
           "rows": rows},
          open("kv7_controls.json", "w"), indent=2)
print(df.to_string(index=False))
print("\nA control panel that passes requires |b| < 0.25 AND an interval covering zero.")
print("wrote kv7_controls.csv, kv7_controls.json")

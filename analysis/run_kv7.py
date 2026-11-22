"""Apply dominance.py to the three Kv7 panels.

Their axes are in mV, so a ratio of two effects is at least scale-free -- which
is necessary for an absolute reading but, as section 2 of the manuscript shows,
not sufficient: the zero-effect controls on these same axes return slopes of
+0.29 to +0.91. The panels are reported as descriptive summaries, and the slope
printed here is an effect-weighted panel summary, not a dominance coefficient.
"""


import os as _os


def _data(vendored, fallback):
    """Prefer the copy vendored under dominance/data/, so a checkout of this
    directory alone reproduces everything; fall back to the sibling tree for
    anyone working in the full repository."""
    here = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "data", vendored)
    return here if _os.path.exists(here) else fallback

import sys, numpy as np, pandas as pd
sys.path.insert(0, ".")
import dominance as D

PANELS = {
    "KCNQ4 (deafness), C+E":        _data("kv7/fit_KCNQ4_primary.csv", "../fix-run/codex/scratch/fit_KCNQ4_primary.csv"),
    "KCNQ2 (epilepsy), MaveDB":     _data("kv7/fit_KCNQ2_primary.csv", "../fix-run/codex/scratch/fit_KCNQ2_primary.csv"),
    "KCNQ2 in KCNQ3, Vanoye":       _data("kv7/fit_KCNQ2_KCNQ3_primary.csv", "../fix-run/codex/scratch/fit_KCNQ2_KCNQ3_primary.csv"),
}
print("=" * 78)
print("Effect-transfer slope, Kv7 activation-voltage shifts (both axes in mV)")
print("=" * 78)
res = {}
for name, path in PANELS.items():
    d = pd.read_csv(path)
    d = d[(d.sx > 0) & (d.sy > 0)].copy()
    cls = d.variant.map(D.classify_variant)
    d["klass"] = [c[3] for c in cls]
    gates = D.class_gates(d)
    mis = d[d.klass == "missense"]
    r = D.estimate(mis, units="physical")
    res[name] = r
    print(f"\n### {name}")
    print(f"  classes present: {dict(d.klass.value_counts())}")
    print("  " + str(r).replace("\n", "\n  "))

import json
from datetime import date
json.dump({"generated": date.today().isoformat(),
           "estimator": "through-origin slope; residue-cluster percentile bootstrap (dominance.estimate)",
           "n_boot": 20000, "seed": 20260926,
           "panels": {name: {"source": PANELS[name], "b": r.b, "ci": list(r.ci), "n": r.n,
                             "n_residues": r.n_residues, "pearson_r": r.pearson_r,
                             "units": r.units} for name, r in res.items()},
           "interpretation": "descriptive panel summaries on deposited activation-voltage axes; not dominance coefficients"},
          open("kv7_panels.json", "w"), indent=2)
print("\nwrote kv7_panels.json")
print("\nThese slopes are panel summaries on the deposited physical axes.")
print("They do not measure assembly-competent subunit fractions or identify")
print("a stoichiometric profile. See the manuscript for control limitations.")

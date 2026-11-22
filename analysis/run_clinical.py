

import os as _os


def _data(vendored, fallback):
    """Prefer the copy vendored under dominance/data/, so a checkout of this
    directory alone reproduces everything; fall back to the sibling tree for
    anyone working in the full repository."""
    here = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "data", vendored)
    return here if _os.path.exists(here) else fallback

#!/usr/bin/env python3
"""Do floor-breaching SCN5A variants carry a higher Brugada case burden?

    python3 run_clinical.py [path/to/workbook.xlsx]

A null we report rather than omit.  Cohort counts come from the same O'Neill
supplementary workbook as the electrophysiology, sheet "Population Data":
Walsh consortium Brugada cases and gnomAD population controls per variant.

Two definitions the reader needs.

BREACH.  Two definitions are run and both are reported, because the result
depends entirely on which is used.

  NUMERICAL   het < 54.6%, the simple comparison -- 21 of 51 variants.
              This is the definition behind the numbers reported in the text.
  SIGNIFICANT het significantly below 54.6% under specification (i) of
              section 4, Bonferroni-corrected across the 51-variant family --
              9 of 51.

The numerical definition is the permissive one and is used as primary only
because the significant set, and still more the four-variant intersection, is
too small for a cohort test. Under the significant definition the association
disappears, which is reported alongside rather than buried.

PARTIAL CORRELATION.  Spearman's rho between the predictor and the per-variant
case fraction, with alone-state current partialled out of both by rank
regression; p from the two-sided t approximation on n - 3 degrees of freedom.
Alone-state current is controlled because breaching variants are also more
damaging alone, which is the confound that makes the pooled statistic
misleading.

Writes clinical_null.csv.
"""
import json
import sys
from datetime import date

import numpy as np
import pandas as pd
from scipy import stats

WB = sys.argv[1] if len(sys.argv) > 1 else \
    _data("NIHMS1819453-supplement-ONeill_GiM_Supplemental_xlsx.xlsx", "../NIHMS1819453-supplement-ONeill_GiM_Supplemental_xlsx.xlsx")

ep = pd.read_csv("scn5a_specifications.csv")            # from run_scn5a.py
ctl = pd.read_csv("scn5a_controls.csv").set_index("variant")
FLOOR = float(ctl.loc["WT+WT", "het"]) / 4

pop = pd.ExcelFile(WB).parse("Population Data", header=0)
pop.columns = ["variant", "gnomad", "maf", "walsh"]
for c in ("gnomad", "walsh"):
    pop[c] = pd.to_numeric(pop[c], errors="coerce")
pop = pop[~pop.variant.isin(["WT", "WT+WT"])]

d = ep.merge(pop, on="variant").dropna(subset=["gnomad", "walsh"])
d["carriers"] = d.walsh + d.gnomad
d["case_fraction"] = d.walsh / d.carriers

DEFS = {"numerical (het < floor)": d.het < FLOOR,
        "significant (spec i, Bonferroni)": d.p_fixed < 0.05}
print(f"{len(d)} variants with both electrophysiology and cohort counts; "
      f"floor = {FLOOR:.1f}%")
for k, v_ in DEFS.items():
    print(f"  {int(v_.sum()):2d}/{len(d)} breach under the {k} definition")
print()

rows = []

# --- pooled 2x2: the statistic that inflates -------------------------------
def fisher(sub, label):
    t = [[float(sub.loc[sub.breach, "walsh"].sum()),
          float(sub.loc[sub.breach, "gnomad"].sum())],
         [float(sub.loc[~sub.breach, "walsh"].sum()),
          float(sub.loc[~sub.breach, "gnomad"].sum())]]
    orr, p = stats.fisher_exact(t)
    rows.append({"analysis": label, "unit": "pooled carriers",
                 "n_variants": len(sub), "statistic": "odds ratio",
                 "value": round(float(orr), 3), "p": round(float(p), 5),
                 "note": "treats each carrier as independent; the variant is "
                         "the unit of analysis, so this overstates"})
    return orr, p


# --- variant level: the unit that is correct -------------------------------
v_all = d[d.carriers > 0].copy()


def partial_spearman(a, b, given):
    R = np.array([pd.Series(x).rank().values for x in (a, b, given)])
    R = (R - R.mean(axis=1, keepdims=True)) / R.std(axis=1, keepdims=True)
    ra = R[0] - np.dot(R[0], R[2]) / len(R[2]) * R[2]
    rb = R[1] - np.dot(R[1], R[2]) / len(R[2]) * R[2]
    r = float(np.corrcoef(ra, rb)[0, 1])
    n = len(a)
    t = r * np.sqrt((n - 3) / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), n - 3))


for defname, mask in DEFS.items():
    d["breach"] = mask
    v = v_all.assign(breach=mask.loc[v_all.index])
    fisher(d, f"all variants, pooled [{defname}]")
    fisher(d[d.homo < 10], f"matched to LoF alone (<10% WT), pooled [{defname}]")
    r, p = partial_spearman(v.breach.astype(float), v.case_fraction, v.homo)
    rows.append({"analysis": f"breach indicator, controlling for alone-state "
                             f"current [{defname}]",
                 "unit": "variant", "n_variants": len(v),
                 "statistic": "partial Spearman rho",
                 "value": round(r, 4), "p": round(p, 4), "note": ""})

v = v_all.assign(breach=DEFS["numerical (het < floor)"].loc[v_all.index])
r, p = partial_spearman(v.het, v.case_fraction, v.homo)
rows.append({"analysis": "heterozygous current (continuous), controlling for "
                         "alone-state current",
             "unit": "variant", "n_variants": len(v),
             "statistic": "partial Spearman rho",
             "value": round(r, 4), "p": round(p, 4), "note": ""})

# --- does the dichotomy add anything over the continuous measurement? ------
only = v[v.gnomad == 0]
both = v[v.gnomad > 0]


def auc(score, label):
    s, l = np.asarray(score, float), np.asarray(label, bool)
    pos, neg = s[l], s[~l]
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


lab = (v.gnomad == 0).values
scores = [(f"breach indicator [{k}]", m.loc[v.index].astype(float))
          for k, m in DEFS.items()]
scores += [("heterozygous current", -v.het), ("alone-state current", -v.homo)]
for name, sc in scores:
    rows.append({"analysis": f"{name}: discriminate case-only variants",
                 "unit": "variant", "n_variants": len(v), "statistic": "AUC",
                 "value": round(auc(sc, lab), 4), "p": None,
                 "note": f"{int(lab.sum())} case-only vs {int((~lab).sum())} also in gnomAD"})

df = pd.DataFrame(rows)
df.to_csv("clinical_null.csv", index=False)
json.dump({"generated": date.today().isoformat(), "workbook": WB,
           "floor": FLOOR, "breach_specification": "fixed floor, two-sided "
           "normal, Bonferroni across 51", "rows": rows},
          open("clinical_null.json", "w"), indent=2)
print(df.to_string(index=False))
print("""
Two things make this a null rather than a finding.

1. The pooled odds ratio and the variant-level correlations disagree because
   they use different units of analysis. Pooling treats each carrier as an
   independent observation when the variant is the unit; the variant-level
   result is the correct one and does not reach significance.

2. The association depends entirely on the permissive breach definition.
   Requiring the breach to be statistically significant -- the definition used
   everywhere else in this paper -- removes it: the partial correlation falls
   from +0.28 to +0.02 and the AUC from 0.68 to 0.52. A signal that survives
   only under the looser of two definitions of its own predictor is not
   evidence of clinical utility.""")
print("\nwrote clinical_null.csv, clinical_null.json")

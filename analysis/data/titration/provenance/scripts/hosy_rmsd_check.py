"""Try to reproduce the RMSD values printed in the Hosy & Vivaudou 2014 Fig. 4
legend (0.327 n=4, 0.124 n=3, 0.119 n=2, 0.352 n=1, 0.031 n=2*, 0.071 linear)
from the digitized data, under alternative computational choices.

Usage: python3 -I hosy_rmsd_check.py <fig3I_json> <fig4B_json> <out_json>
"""
import itertools
import json
import sys

import numpy as np

F3, F4, OUT = sys.argv[1:4]
f3 = json.load(open(F3))
f4 = json.load(open(F4))
R = {b["wt_fraction"]: b["mean"] for b in f3["bars"]}
D4 = {r["p_nominal"]: r["value"] for r in f4["points"]}
P4x = {r["p_nominal"]: r["p"] for r in f4["points"]}   # x positions as plotted in Fig 4B
PAPER = {"n4": 0.327, "n3": 0.124, "n2": 0.119, "n1": 0.352, "n2s": 0.031, "lin": 0.071}


def model(name, p, grid=None):
    def f(p):
        q = 1 - p
        P4, P3, P2 = p ** 4, 4 * p ** 3 * q, 6 * p ** 2 * q ** 2
        return {"n4": P4, "n3": P4 + P3, "n2": P4 + P3 + P2, "n1": 1 - q ** 4,
                "n2s": P4 + P3 + 0.5 * P2, "n2adj46": P4 + P3 + (4 / 6) * P2, "lin": p}[name]
    if grid is None:
        return f(p)
    xs = np.arange(0, 1 + 1e-9, grid)
    return float(np.interp(p, xs, [f(x) for x in xs]))


norms = {
    "(R-R0)/(R1-R0)": lambda r, R: (r - R[0.0]) / (R[1.0] - R[0.0]),
    "R/R1": lambda r, R: r / R[1.0],
    "(R-R0)/(1-R0)": lambda r, R: (r - R[0.0]) / (1 - R[0.0]),
    "raw R": lambda r, R: r,
    # endpoints rounded as quoted in the Results text ("0.05 ... 1.1")
    "(R-0.05)/(1.1-0.05), endpoints fixed 1/0": lambda r, R: (1.0 if r == R[1.0] else 0.0 if r == R[0.0]
                                                             else (r - 0.05) / (1.1 - 0.05)),
}
point_sets = {"all 8": [1.0, 0.8, 0.6, 0.5, 0.4, 0.2, 0.16, 0.0],
              "6 interior": [0.8, 0.6, 0.5, 0.4, 0.2, 0.16],
              "7 (no p=1)": [0.8, 0.6, 0.5, 0.4, 0.2, 0.16, 0.0],
              "7 (no p=0)": [1.0, 0.8, 0.6, 0.5, 0.4, 0.2, 0.16]}
rows = []
for (nname, nf), (sname, ps), ddof, xmode, grid, src in itertools.product(
        norms.items(), point_sets.items(), (0, 1), ("nominal p", "Fig4B-plotted p"), (None, 0.05, 0.1),
        ("Fig3I", "Fig4B")):
    if src == "Fig4B" and nname != "(R-R0)/(R1-R0)":
        continue
    vals = {}
    for p in ps:
        d = nf(R[p], R) if src == "Fig3I" else D4[p]
        x = p if xmode == "nominal p" else P4x[p]
        vals[p] = (x, d)
    out = {}
    for m in ("n4", "n3", "n2", "n1", "n2s", "n2adj46", "lin"):
        sq = [(d - model(m, x, grid)) ** 2 for x, d in vals.values()]
        out[m] = float(np.sqrt(np.sum(sq) / (len(sq) - ddof)))
    score = float(np.sqrt(np.mean([(out[k] - PAPER[k]) ** 2 for k in PAPER])))
    rows.append({"normalization": nname, "points": sname, "ddof": ddof, "x": xmode, "model_grid": grid,
                 "data": src, "rmsd": out, "rms_mismatch_vs_paper": score})
rows.sort(key=lambda r: r["rms_mismatch_vs_paper"])
json.dump(rows, open(OUT, "w"), indent=1)
print("best 12 combinations (rms mismatch of the 6 RMSDs vs paper):")
for r in rows[:12]:
    o = r["rmsd"]
    print("%.4f | %-15s %-11s ddof=%d %-15s grid=%-4s %-6s | n4 %.3f n3 %.3f n2 %.3f n1 %.3f 2* %.3f (4/6: %.3f) lin %.3f" % (
        r["rms_mismatch_vs_paper"], r["normalization"], r["points"], r["ddof"], r["x"], r["model_grid"], r["data"],
        o["n4"], o["n3"], o["n2"], o["n1"], o["n2s"], o["n2adj46"], o["lin"]))
print("paper:                                                            | n4 0.327 n3 0.124 n2 0.119 n1 0.352 2* 0.031 lin 0.071")

# How large a change in the data would be needed to reproduce the paper's six
# RMSDs (standard normalization, 8 points)? Minimise ||delta||^2 subject to a
# good RMSD match (penalty), interior points only.
from scipy.optimize import minimize
ps = [0.8, 0.6, 0.5, 0.4, 0.2, 0.16]
base = np.array([(R[p] - R[0.0]) / (R[1.0] - R[0.0]) for p in ps])
def rm(d):
    o = {}
    for m in PAPER:
        sq = [(di - model(m, p)) ** 2 for di, p in zip(d, ps)]
        o[m] = np.sqrt(np.sum(sq) / 8.0)
    return o
def obj(delta, lam=1e3):
    o = rm(base + delta)
    return float(np.sum(delta ** 2) + lam * np.sum([(o[k] - PAPER[k]) ** 2 for k in PAPER]))
best = minimize(obj, np.zeros(len(ps)), method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-6, "fatol": 1e-10})
o = rm(base + best.x)
print("\nminimal data perturbation that best reproduces the paper's RMSDs (8-point denominator):")
for p, b0, dl in zip(ps, base, best.x):
    print("  p=%.2f  digitized %.4f  required %.4f  delta %+.4f" % (p, b0, b0 + dl, dl))
print("  resulting RMSDs:", {k: round(float(v), 3) for k, v in o.items()})
json.dump({"combos": rows[:50], "perturbation": {"p": ps, "digitized": base.tolist(), "delta": best.x.tolist(),
           "resulting_rmsd": {k: float(v) for k, v in o.items()}}}, open(OUT, "w"), indent=1)

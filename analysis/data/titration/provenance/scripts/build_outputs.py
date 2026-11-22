"""Assemble titration_digitized.csv and titration_model_predictions.csv from the
JSON outputs of the digitization scripts.

Usage: python3 -I build_outputs.py <work_dir> <out_dir>
"""
import csv
import json
import sys

WORK, OUTDIR = sys.argv[1:3]
f3 = json.load(open(f"{WORK}/hosy_fig3I.json"))
f4 = json.load(open(f"{WORK}/hosy_fig4B.json"))
xc = json.load(open(f"{WORK}/hosy_crosscheck.json"))
cv = json.load(open(f"{WORK}/clatot_vector_points.json"))

COLS = ["study", "figure_panel", "construct", "condition_label", "wt_fraction_input", "mutant_fraction_input",
        "measure", "units_or_normalization", "value", "err_low", "err_high", "err_type", "n",
        "resolution_data_units", "source", "notes"]
rows = []
fmt = lambda v: "NA" if v is None else f"{v:.4f}"

# ------------------------------- Hosy 2014 --------------------------------
HOSY = "Hosy & Vivaudou 2014 Front Physiol 5:11 (doi:10.3389/fphys.2014.00011)"
CONSTRUCT_H = "Kir6.2 + SUR2A WT/SUR2A-DNN (D832N+D1469N) cRNA mixture"
MEAS_H = ("I(ADP)/I(control): steady-state inside-out macropatch current in 300 uM MgADP (1 mM Mg2+, -50 mV) "
          "divided by the current measured just before in nucleotide-free solution")
UNIT_H = "dimensionless per-patch ratio (1 = no net effect of ADP); not normalized to WT"
CTRL_NA = {"1::0": "4.4 +- 0.6", "0.8::0.2": "1 +- 0.5", "0.6::0.4": "1.3 +- 0.4", "0.5::0.5": "5.3 +- 1.4",
           "0.4::0.6": "1.4 +- 0.4", "0.2::0.8": "1.6 +- 0.7", "0.16::0.84": "3.5 +- 1", "0::1": "1.5 +- 0.3"}
res_h = 1.0 / f3["calibration"]["px_per_unit"]
f2d = xc["fig2D"]["bars"]; f1a = xc["fig1A"]["bars"]
EXTRA_H = {
    "1::0": (f"Cross-checks (same 87 patches re-plotted): Fig 2D {f2d['1::0 (n=87)']['mean']:.3f} "
             f"(SEM {f2d['1::0 (n=87)']['sem']:.3f}); Fig 1A 300 uM WT {f1a['WT 300 uM (n=87)']['mean']:.3f} "
             f"(SEM {f1a['WT 300 uM (n=87)']['sem']:.3f}). Text: ratio 'ranged from 0.05 ... to 1.1 for pure WT'."),
    "0.5::0.5": (f"Cross-check (same 19 patches): Fig 2D {f2d['0.5::0.5 (n=19)']['mean']:.3f} "
                 f"(SEM {f2d['0.5::0.5 (n=19)']['sem']:.3f}). Text (Fig 1C, equimolar mix): 300 uM MgADP led to "
                 "a ~50% decrease in current."),
    "0::1": ("Cap overlaps the bar outline (seen only as a ~1.6 px bump above it); cap centre = bump upper edge + "
             "half cap stroke, so the SEM is resolution-limited (+-0.003). "
             f"Cross-checks (same 38 patches): Fig 1A 300 uM DNN {f1a['DNN 300 uM (n=38)']['mean']:.3f} "
             f"(SEM ~{f1a['DNN 300 uM (n=38)']['sem']:.3f}, cap also merged); Fig 2D {f2d['0::1 (n=38)']['mean']:.3f} "
             f"(SEM {f2d['0::1 (n=38)']['sem']:.3f}, cap separate). Text: 'ratio ... 0.05 for pure mutant'."),
    "0.2::0.8": "Axis label printed all in black (other labels show the DNN part in red); no other ambiguity.",
}
for b in f3["bars"]:
    lab = b["label"]
    wt = b["wt_fraction"]
    note = (f"PRIMARY. Bar chart: value = centre of the bar's top outline; err_high = centre of the cap. "
            f"Only the upper half of the SEM bar is drawn: err_low not drawn (symmetric lower end would be "
            f"{b['mean'] - b['sem']:.4f}; SEM = {b['sem']:.4f}). n = patches (number printed above bar). "
            f"Pre-ADP current in 0 ATP (legend): {CTRL_NA[lab]} nA (mean +- SEM). "
            f"cRNA: 2 ng Kir6.2 + 6 ng SUR2A (WT+DNN mixed at the stated proportion) per oocyte.")
    if lab in EXTRA_H:
        note += " " + EXTRA_H[lab]
    rows.append({"study": HOSY, "figure_panel": "Fig 3I", "construct": CONSTRUCT_H,
                 "condition_label": f"SUR2A WT::DNN {lab}", "wt_fraction_input": f"{wt:g}",
                 "mutant_fraction_input": f"{1 - wt:g}" if wt not in (0.16,) else "0.84",
                 "measure": MEAS_H, "units_or_normalization": UNIT_H, "value": fmt(b["mean"]),
                 "err_low": "NA", "err_high": fmt(b["err_high"]), "err_type": "SEM (upper half drawn only)",
                 "n": str(b["n"]), "resolution_data_units": f"{res_h:.4f}", "source": "digitized", "notes": note})

# Fig 4B (authors' normalization of Fig 3I)
MEAS_4B = ("normalized ADP-activated current = (R(p) - R(0)) / (R(1) - R(0)), R = Fig 3I ratio, "
           "p = WT fraction (authors' normalization; curves = fraction of channels with >= n WT subunits)")
res_4 = 1.0 / f4["calibration"]["px_per_unit_y"]
rad_4 = f4["marker_radius_px"] / f4["calibration"]["px_per_unit_y"]
lab_of = {1.0: "1::0", 0.8: "0.8::0.2", 0.6: "0.6::0.4", 0.5: "0.5::0.5", 0.4: "0.4::0.6", 0.2: "0.2::0.8",
          0.16: "0.16::0.84", 0.0: "0::1"}
n_of = {b["wt_fraction"]: b["n"] for b in f3["bars"]}
for r in sorted(f4["points"], key=lambda r: -r["p_nominal"]):
    p = r["p_nominal"]
    overlap = p in (0.2, 0.16)
    note = (f"DERIVED, not independent: the authors' normalized re-plot of Fig 3I. Marker = black disc "
            f"(radius {f4['marker_radius_px']:.1f} px = {rad_4:.3f} units), centre from distance-transform maximum. "
            f"Plotted at p = {r['p']:.3f}. Normalizing my Fig 3I readings gives {r['fig3I_normalized_mean']:.4f} "
            f"(diff {r['diff_vs_fig3I']:+.4f}); normalized Fig 3I SEM = {r['fig3I_normalized_sem']:.4f}.")
    if r["err_high"] is None and r["err_low"] is None:
        note += " Error bar hidden inside the marker (half-length < marker radius)."
    if p == 1.0:
        note += " Marker sits on the y-axis line; caps read left of the axis."
    if overlap:
        note += (" Markers for p = 0.2 and 0.16 overlap; separated cleanly by distance-transform peaks, but "
                 "resolution set to the marker radius as instructed.")
    if p == 0.16:
        note += " Note: plotted ~0.008 to the right of p = 0.16 (at p = 0.168)."
    if p == 0.2:
        note += (" Drawn error bar (+-~0.029) is shorter than the normalized Fig 3I SEM (0.040): "
                 "inter-figure inconsistency.")
    rows.append({"study": HOSY, "figure_panel": "Fig 4B", "construct": CONSTRUCT_H,
                 "condition_label": f"p = {p:g} (WT::DNN {lab_of[p]})", "wt_fraction_input": f"{p:g}",
                 "mutant_fraction_input": f"{round(1 - p, 2):g}", "measure": MEAS_4B,
                 "units_or_normalization": "dimensionless; 1 = pure WT, 0 = pure DNN by construction",
                 "value": fmt(r["value"]), "err_low": fmt(r["err_low"]), "err_high": fmt(r["err_high"]),
                 "err_type": "SEM (as plotted, presumably Fig 3I SEM / (R(1)-R(0)))", "n": str(n_of[p]),
                 "resolution_data_units": f"{(rad_4 if overlap else res_4):.4f}", "source": "digitized",
                 "notes": note})

# ------------------------------- Clatot 2017 -------------------------------
CLAT = "Clatot et al. 2017 Nat Commun 8:2077 (doi:10.1038/s41467-017-02262-0)"
NT1 = {  # Supplementary Table 1, n per mutant fraction
    "WT+L325R": {0: 62, 0.1: 16, 0.25: 33, 0.5: 51, 0.75: 49, 0.9: 30, 1: 12},
    "WT+R104W": {0: 55, 0.1: None, 0.25: 19, 0.5: 39, 0.75: 13, 0.9: 13, 1: 15},
    "WT+R878C": {0: 31, 0.1: None, 0.25: 17, 0.5: 31, 0.75: 19, 0.9: None, 1: 12},
    "S460A+L325R": {0: 45, 0.1: 30, 0.25: 16, 0.5: 12, 0.75: 39, 0.9: 16, 1: 12},
    "S460A+R104W": {0: 46, 0.1: 15, 0.25: 28, 0.5: 14, 0.75: 18, 0.9: 14, 1: 12},
    "WT+L325R+Difopein": {0: 22, 0.1: 27, 0.25: 32, 0.5: 30, 0.75: 30, 0.9: 12, 1: 12},
}
SER = {("fig1c", "circle"): ("WT+L325R", "Nav1.5 WT + L325R (DN)", "WT"),
       ("fig1c", "triangle"): ("WT+L325R+Difopein", "Nav1.5 WT + L325R (DN) + difopein (14-3-3 inhibitor)", "WT"),
       ("fig1c", "square"): ("S460A+L325R", "Nav1.5 S460A + L325R (DN); S460A is the 'wild-type' partner", "S460A"),
       ("figS3C", "circle"): ("WT+R104W", "Nav1.5 WT + R104W (DN)", "WT"),
       ("figS3C", "square"): ("S460A+R104W", "Nav1.5 S460A + R104W (DN); S460A is the 'wild-type' partner", "S460A"),
       ("figS3C", "triangle"): ("WT+R878C", "Nav1.5 WT + R878C (non-conducting, non-DN control)", "WT")}
PANEL = {"fig1c": "Fig 1c", "figS3C": "Supplementary Fig 3C"}
MEAS_C = ("peak whole-cell INa density at -20 mV (HEK293, I/V protocol from -120 mV holding), normalized to the "
          "mean of the 0%-mutant condition of the same series")
TEXT_RATIO = {0.1: "1:10", 0.25: "1:4", 0.5: "1:1", 0.75: "4:1", 0.9: "10:1"}
LIT = {0.1: 1 / 11, 0.25: 1 / 5, 0.5: 0.5, 0.75: 4 / 5, 0.9: 10 / 11}
res_c = {k: max(abs(v["calibration"]["x_max_tick_residual"]), abs(v["calibration"]["y_max_tick_residual"]), 0.0003)
         for k, v in cv.items()}
for key in ("fig1c", "figS3C"):
    cal = cv[key]["calibration"]
    for r in cv[key]["points"]:
        series, cons, wtname = SER[(key, r["kind"])]
        x = round(r["x_val"] * 20) / 20 if abs(r["x_val"] - 0.1) > 0.02 and abs(r["x_val"] - 0.9) > 0.02 \
            else round(r["x_val"], 1)
        x = round(x, 2)
        n = NT1[series].get(x if x not in (0.0, 1.0) else int(x))
        mut = series.split("+")[1]
        lab = (f"{mut}:{wtname} mutant fraction {x:g}" +
               (f" (text: {mut}:WT {TEXT_RATIO[x]}; SI: {x * 100:g}% + {100 - x * 100:g}%)" if x in TEXT_RATIO else
                (" (0% mutant; normalization reference)" if x == 0 else " (100% mutant)")))
        notes = [f"PRIMARY. Read from the vector drawing in the PDF (marker = {r['kind']}); "
                 f"x plotted at {round(r['x_val'], 4) + 0.0:.4f}."]
        if x in LIT and x != 0.5:
            notes.append(f"Plotted x matches the SI fractions; a literal {TEXT_RATIO[x]} cDNA ratio would be "
                         f"{LIT[x]:.3f}.")
        if r["top_val"] is None:
            notes.append("No error bar drawn (0 within resolution; text: no current detected for 100% mutant)."
                         if x == 1 else "No error bar drawn.")
        if x == 0:
            notes.append("Equals 1 by construction (each series normalized to its own 0% mean); "
                         "the bar shows the SEM of that normalization group.")
        if r["top_val"] is not None and not r["top_cap"]:
            notes.append(f"Upper bar clipped at the axis maximum ({cal['y_axis_top_val']:.4f}), no cap drawn; "
                         f"symmetric upper end would be {2 * r['mean_val'] - r['bottom_val']:.4f}.")
        if r["top_val"] is not None and r["top_cap"]:
            notes.append(f"SEM (half-length) = {r['mean_val'] - r['bottom_val']:.4f} (lower), "
                         f"{(r['top_val'] - r['mean_val']):.4f} (upper).")
        elif r["top_val"] is not None:
            notes.append(f"SEM (half-length) = {r['mean_val'] - r['bottom_val']:.4f} (lower half; upper half clipped).")
        notes.append("cDNA: 0.3 ug total per 35-mm dish; mutant fraction varied at constant total (Methods).")
        if n is None:
            notes.append("Supplementary Table 1 lists ND for this cell.")
        if series == "WT+L325R+Difopein":
            notes.append("Difopein plasmid amount not stated.")
        if series.startswith("S460A"):
            notes.append("S460A-Nav1.5 replaces WT as the partner; wt_fraction = S460A cDNA fraction.")
        if series == "WT+R878C":
            notes.append("Control: trafficking-competent, non-conducting, non-DN mutant (peer-review response: "
                         "'a reduction of only 50%' at 50:50).")
        if series in ("WT+L325R", "WT+R104W") and x == 0.5:
            notes.append("Peer-review response: DN mutants give 'a 75% decrease in INa' at 50:50.")
        rows.append({"study": CLAT, "figure_panel": PANEL[key], "construct": cons, "condition_label": lab,
                     "wt_fraction_input": f"{1 - x:g}", "mutant_fraction_input": f"{x:g}", "measure": MEAS_C,
                     "units_or_normalization": "fraction of the 0%-mutant mean current density (dimensionless)",
                     "value": fmt(r["mean_val"]),
                     "err_low": fmt(r["bottom_val"]) if r["top_val"] is not None else "NA",
                     "err_high": fmt(r["top_val"]) if (r["top_val"] is not None and r["top_cap"]) else "NA",
                     "err_type": "SEM" if r["top_val"] is not None else "SEM (no bar drawn)",
                     "n": "NA" if n is None else str(n),
                     "resolution_data_units": f"{res_c[key]:.4f}", "source": "digitized", "notes": " ".join(notes)})

with open(f"{OUTDIR}/titration_digitized.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS)
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f"wrote {len(rows)} rows")


# ---------------------------- model predictions ----------------------------
def hosy_models(p):
    q = 1 - p
    P4, P3, P2 = p ** 4, 4 * p ** 3 * q, 6 * p ** 2 * q ** 2
    return [("n>=4 WT (p^4)", P4, "plotted"), ("n>=3 WT", P4 + P3, "plotted"),
            ("n>=2 WT", P4 + P3 + P2, "plotted"), ("n>=1 WT (1-(1-p)^4)", 1 - q ** 4, "plotted"),
            ("n>=2* as plotted (p^4+4p^3q+3p^2q^2)", P4 + P3 + 0.5 * P2, "plotted; verified against curve pixels"),
            ("no co-assembly (linear, p)", p, "plotted (dotted green)"),
            ("2 adjacent WT, combinatorial weight 4/6 (p^4+4p^3q+4p^2q^2)", P4 + P3 + (4 / 6) * P2,
             "NOT in paper; reference only"),
            ("2 opposite WT, combinatorial weight 2/6 (p^4+4p^3q+2p^2q^2)", P4 + P3 + (2 / 6) * P2,
             "NOT in paper; reference only")]


mrows = []
for p in (1.0, 0.8, 0.6, 0.5, 0.4, 0.2, 0.16, 0.0):
    for name, v, st in hosy_models(p):
        mrows.append(["Hosy & Vivaudou 2014", "Fig 4A/B", name, f"{p:g}", f"{round(1 - p, 2):g}", f"{v:.4f}",
                      "fraction of channels activatable by ADP (compare with Fig 4B normalized data)", st])
for x in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0):
    for name, v, st in (("Expected monomer, 1-x", 1 - x, "plotted; matches vector paths (max dev <=0.0005 in both panels)"),
                        ("Expected dimer, (1-x)^2", (1 - x) ** 2, "plotted; matches vector paths (max dev 0.0007 Fig 1c, 0.0004 Supp Fig 3C)"),
                        ("Trimer per Methods, (1-x)^3", (1 - x) ** 3, "stated in Methods but NOT what is plotted"),
                        ("Curve labelled 'Expected trimer' as plotted, (1-x)^4", (1 - x) ** 4,
                         "plotted dash-dot curve; matches (1-x)^4 (max dev 0.0010 Fig 1c, 0.0013 Supp Fig 3C), i.e. a tetramer curve; deviates from (1-x)^3 by up to 0.106")):
        mrows.append(["Clatot et al. 2017", "Fig 1c and Supplementary Fig 3C", name, f"{1 - x:g}", f"{x:g}",
                      f"{v:.4f}", "normalized current density predicted for a DN mutant", st])
with open(f"{OUTDIR}/titration_model_predictions.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["study", "figure_panel", "model", "wt_fraction", "mutant_fraction", "predicted_value",
                "quantity", "status"])
    w.writerows(mrows)
print(f"wrote {len(mrows)} model rows")

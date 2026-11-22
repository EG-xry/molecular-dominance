"""Print the Markdown tables used in DIGITIZATION_REPORT.md from the JSON
outputs (so that no number in the report is hand-transcribed).

Usage: python3 -I report_tables.py <work_dir> > report_tables.md
"""
import json
import sys

W = sys.argv[1]
f3 = json.load(open(f"{W}/hosy_fig3I.json"))
f4 = json.load(open(f"{W}/hosy_fig4B.json"))
xc = json.load(open(f"{W}/hosy_crosscheck.json"))
rm = json.load(open(f"{W}/hosy_rmsd_check.json"))
cv = json.load(open(f"{W}/clatot_vector_points.json"))
px = json.load(open(f"{W}/clatot_pixel_calibration.json"))

print("### T1. Hosy Fig 3I y-axis calibration (image px; tick strokes in cols 197-208)\n")
print("| tick label | row (px) | stroke width (px) | fit residual (units) |\n|---|---|---|---|")
for t in f3["calibration"]["ticks"]:
    print(f"| {t['value']:.1f} | {t['row']:.3f} | {t['stroke_px']:.2f} | {t['fit_residual_units']:+.4f} |")
c = f3["calibration"]
print(f"\nvalue = {c['value_at_row0']:.6f} {c['value_per_row']:+.8f} x row; {c['px_per_unit']:.2f} px per unit "
      f"(1 px = {1 / c['px_per_unit']:.4f}); max |residual| {c['max_abs_residual_units']:.4f}.\n")
print("Checks: dashed y = 1 line reads " + ", ".join(f"{d['value']:.4f}" for d in f3["checks"]["dashed_line_y1"]) +
      "; bottom outlines of the 8 bars read " +
      ", ".join(sorted(set(f"{d['value']:+.4f}" for d in f3["checks"]["bar_bottoms"]))) + ".\n")

print("### T2. Hosy Fig 3I readings\n")
print("| WT::DNN | n | bar-top row | mean | cap row | mean+SEM | SEM | stem col |\n|---|---|---|---|---|---|---|---|")
for b in f3["bars"]:
    print(f"| {b['label']} | {b['n']} | {b['bar_top_row']:.2f} | {b['mean']:.4f} | {b['cap_row']:.2f} | "
          f"{b['err_high']:.4f} | {b['sem']:.4f} | {b['stem_col']} |")

print("\n### T3. Hosy Fig 4B calibration and readings\n")
c4 = f4["calibration"]
print("x ticks (p, col px): " + ", ".join(f"{p:g} @ {col:.2f}" for p, col in c4["x_ticks"]) +
      f"; {c4['px_per_unit_x']:.2f} px per unit; max |residual| {max(abs(r) for r in c4['x_residuals']):.4f}.\n")
print("y ticks (value, row px): " + ", ".join(f"{v:g} @ {r:.2f}" for v, r in c4["y_ticks"]) +
      f"; {c4['px_per_unit_y']:.2f} px per unit; max |residual| {max(abs(r) for r in c4['y_residuals']):.4f}.\n")
print("| p (nominal) | p as plotted | value | lower end | upper end | Fig 3I normalized | diff | Fig 3I SEM normalized |\n"
      "|---|---|---|---|---|---|---|---|")
for r in sorted(f4["points"], key=lambda r: -r["p_nominal"]):
    lo = "hidden" if r["err_low"] is None else f"{r['err_low']:.4f}"
    hi = "hidden" if r["err_high"] is None else f"{r['err_high']:.4f}"
    print(f"| {r['p_nominal']:g} | {r['p']:.3f} | {r['value']:.4f} | {lo} | {hi} | {r['fig3I_normalized_mean']:.4f} | "
          f"{r['diff_vs_fig3I']:+.4f} | {r['fig3I_normalized_sem']:.4f} |")

print("\n### T4. Same patches, three renderings (Hosy)\n")
b3 = {b["label"]: b for b in f3["bars"]}
f2 = xc["fig2D"]["bars"]; f1 = xc["fig1A"]["bars"]
print("| condition (n) | Fig 3I mean / SEM | Fig 2D mean / SEM | Fig 1A (300 uM) mean / SEM |\n|---|---|---|---|")
print(f"| 1::0 (87) | {b3['1::0']['mean']:.4f} / {b3['1::0']['sem']:.4f} | {f2['1::0 (n=87)']['mean']:.4f} / "
      f"{f2['1::0 (n=87)']['sem']:.4f} | {f1['WT 300 uM (n=87)']['mean']:.4f} / {f1['WT 300 uM (n=87)']['sem']:.4f} |")
print(f"| 0.5::0.5 (19) | {b3['0.5::0.5']['mean']:.4f} / {b3['0.5::0.5']['sem']:.4f} | "
      f"{f2['0.5::0.5 (n=19)']['mean']:.4f} / {f2['0.5::0.5 (n=19)']['sem']:.4f} | not shown |")
print(f"| 0::1 (38) | {b3['0::1']['mean']:.4f} / {b3['0::1']['sem']:.4f} | {f2['0::1 (n=38)']['mean']:.4f} / "
      f"{f2['0::1 (n=38)']['sem']:.4f} | {f1['DNN 300 uM (n=38)']['mean']:.4f} / {f1['DNN 300 uM (n=38)']['sem']:.4f} |")
fx = xc["fig2D"]
print(f"\nFig 2D ticks (col px): 0 @ {fx['ticks'][0][1]:.2f}, 0.5 @ {fx['ticks'][1][1]:.2f}, 1 @ {fx['ticks'][2][1]:.2f}; "
      f"spacings {fx['spacing_0_to_0.5_px']:.2f} and {fx['spacing_0.5_to_1_px']:.2f} px; a linear axis puts 0.5 at "
      f"{fx['where_0.5_tick_should_be_if_linear']:.2f}; bar base at {fx['bar_base_col']:.2f} "
      f"(= {fx['bar_base_value_0_1_ticks']:+.4f}). Fig 1A: {xc['fig1A']['px_per_unit']:.2f} px per unit, "
      f"dashed y = 1 line reads {xc['fig1A']['dashed_line_y1_value']:.4f}.\n")

print("### T5. RMSD of the normalized data from each model (Hosy Fig 4 legend vs recomputed)\n")
best_std = [r for r in rm["combos"] if r["normalization"] == "(R-R0)/(R1-R0)" and r["points"] == "all 8"
            and r["ddof"] == 0 and r["x"] == "nominal p" and r["model_grid"] is None and r["data"] == "Fig3I"]
best_rnd = [r for r in rm["combos"] if r["normalization"].startswith("(R-0.05)") and r["points"] == "all 8"
            and r["ddof"] == 0 and r["x"] == "nominal p" and r["model_grid"] is None and r["data"] == "Fig3I"]
b4 = [r for r in rm["combos"] if r["normalization"] == "(R-R0)/(R1-R0)" and r["points"] == "all 8"
      and r["ddof"] == 0 and r["x"] == "nominal p" and r["model_grid"] is None and r["data"] == "Fig4B"]
names = [("n4", "n >= 4"), ("n3", "n >= 3"), ("n2", "n >= 2"), ("n1", "n >= 1"), ("n2s", "n >= 2*"),
         ("lin", "no co-assembly (linear)"), ("n2adj46", "2 adjacent, weight 4/6 (not in paper)")]
paper = {"n4": "0.327", "n3": "0.124", "n2": "0.119", "n1": "0.352", "n2s": "0.031", "lin": "0.071", "n2adj46": "-"}
print("| model | legend | Fig 3I normalized with plotted endpoints | Fig 4B digitized | Fig 3I normalized with endpoints 1.1 / 0.05 |\n"
      "|---|---|---|---|---|")
for k, nm in names:
    print(f"| {nm} | {paper[k]} | {best_std[0]['rmsd'][k]:.3f} | {b4[0]['rmsd'][k]:.3f} | {best_rnd[0]['rmsd'][k]:.3f} |")
print(f"\nRMS mismatch of the six legend values: {best_std[0]['rms_mismatch_vs_paper']:.4f} (plotted endpoints), "
      f"{b4[0]['rms_mismatch_vs_paper']:.4f} (Fig 4B), {best_rnd[0]['rms_mismatch_vs_paper']:.4f} (rounded endpoints).\n")

for key, title in (("fig1c", "Clatot Fig 1c (main PDF p. 2)"), ("figS3C", "Clatot Supplementary Fig 3C (SI PDF p. 3)")):
    cal = cv[key]["calibration"]; pc = px[key]
    print(f"### Calibration, {title}\n")
    print("| axis | label | PDF pt | 600-dpi crop px |\n|---|---|---|---|")
    for v, p_, pt in pc["x_ticks_value_px_pt"]:
        print(f"| x | {v:g} | {pt:.4f} | {p_:.2f} |")
    for v, p_, pt in pc["y_ticks_value_px_pt"]:
        print(f"| y | {v:g} | {pt:.4f} | {p_:.2f} |")
    print(f"\n{cal['pt_per_x_unit']:.3f} pt per x unit, {cal['pt_per_y_unit']:.3f} pt per y unit; max tick-fit residual "
          f"{cal['x_max_tick_residual']:.5f} (x), {cal['y_max_tick_residual']:.5f} (y); y-axis line ends at "
          f"{cal['y_axis_top_val']:.4f}. Crop origin {pc['crop_origin_px']} px at 600 dpi.\n")
    print("| marker | x (mutant fraction) | mean | lower end | upper end | upper cap drawn |\n|---|---|---|---|---|---|")
    for r in cv[key]["points"]:
        lo = "none" if r["bottom_val"] is None else f"{r['bottom_val']:.4f}"
        hi = "none" if r["top_val"] is None else f"{r['top_val']:.4f}"
        print(f"| {r['kind']} | {r['x_val']:.4f} | {r['mean_val']:.4f} | {lo} | {hi} | {r['top_cap']} |")
    print()

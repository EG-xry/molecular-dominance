"""Digitize Hosy & Vivaudou 2014 Fig. 4B (normalized ADP-activated current vs
WT fraction p, with the authors' binomial model curves).

Usage: python3 -I hosy_fig4B_digitize.py <fig4_png> <fig3I_json> <out_json> <overlay_png>

Image: lossless 1450 x 1336 px Figure 4 embedded in the publisher PDF
(pdfimages -all, page 5, image 3). blackness = 255 - max(R, G, B).

  * x axis (p, reversed: 1 at left, 0 at right): tick centres from half-max
    edges; the p = 0 tick is hidden under the p = 0 marker, so it is fitted
    from the other five ticks.
  * y axis: ticks 0.25, 0.5, 0.75 and the x-axis line (y = 0); the y = 1 tick
    is hidden under the p = 1 marker.
  * markers: black discs (~26 px). Thin strokes (curves, error bars) are removed
    with a binary opening (disc radius 7) and disc centres are taken at
    distance-transform maxima (this separates the overlapping p = 0.2 and
    p = 0.16 markers).
  * error bars: vertical stem through the marker centre; the cap row is the
    half-max-edge midpoint of the cap stroke (stem columns excluded). Bars
    shorter than the marker radius are hidden and reported as such.
  * model curves are NOT digitized: they are recomputed from the formulas in
    the paper and spot-checked against the coloured curve pixels.
"""
import json
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi

SRC, F3I, OUT, OVL = sys.argv[1:5]
img = Image.open(SRC).convert("RGB")
A = np.asarray(img).astype(float)
K = 255.0 - A.max(axis=2)


def stroke_centre(prof, coords):
    prof = np.asarray(prof, float)
    i = int(np.argmax(prof))
    peak = prof[i]
    b0 = float(np.median(prof[:3])); b1 = float(np.median(prof[-3:]))
    h0 = (peak + b0) / 2.0; h1 = (peak + b1) / 2.0
    j = i
    while prof[j - 1] > h0:
        j -= 1
    e0 = coords[j - 1] + (h0 - prof[j - 1]) / (prof[j] - prof[j - 1])
    j = i
    while prof[j + 1] > h1:
        j += 1
    e1 = coords[j] + (prof[j] - h1) / (prof[j] - prof[j + 1])
    return float((e0 + e1) / 2.0), float(e1 - e0)


# ------------------------------ calibration ------------------------------
xt = []
for p, c in ((1.0, 173), (0.8, 300), (0.6, 426), (0.4, 554), (0.2, 680)):
    if p == 1.0:
        continue  # the p=1 tick coincides with the y-axis line (col 173); use as check only
    cc = np.arange(c - 7, c + 8, dtype=float)
    xt.append((p, stroke_centre(K[1232:1241, c - 7:c + 8].mean(axis=0), cc)[0]))
# y-axis line column (p = 1 position)
cc = np.arange(165, 182, dtype=float)
yaxis_col = stroke_centre(K[1000:1100, 165:182].mean(axis=0), cc)[0]
xt.append((1.0, yaxis_col))
xt.sort(key=lambda t: t[1])
bx, ax = np.polyfit([c for p, c in xt], [p for p, c in xt], 1)
x_res = [p - (ax + bx * c) for p, c in xt]
p_of = lambda col: float(ax + bx * col)
col_of = lambda p: float((p - ax) / bx)

yt = []
for v, r in ((0.25, 1109), (0.5, 991), (0.75, 872)):
    rr = np.arange(r - 7, r + 8, dtype=float)
    yt.append((v, stroke_centre(K[r - 7:r + 8, 160:171].mean(axis=1), rr)[0]))
rr = np.arange(1218, 1239, dtype=float)
xaxis_row = stroke_centre(K[1218:1239, 300:400].mean(axis=1), rr)[0]
yt.append((0.0, xaxis_row))
yt.sort(key=lambda t: t[1])
by, ay = np.polyfit([r for v, r in yt], [v for v, r in yt], 1)
y_res = [v - (ay + by * r) for v, r in yt]
y_of = lambda row: float(ay + by * row)
row_of = lambda v: float((v - ay) / by)

# -------------------------------- markers --------------------------------
r0, r1, c0, c1 = 700, 1260, 130, 840
mask = K[r0:r1, c0:c1] > 128
yy, xx = np.mgrid[-7:8, -7:8]
op = ndi.binary_opening(mask, structure=(xx ** 2 + yy ** 2) <= 49)
dt = ndi.distance_transform_edt(op)
peaks = np.argwhere((dt == ndi.maximum_filter(dt, size=9)) & (dt > 10))
# refine each centre: mean of the pixels within 2 px of the DT maximum value
centres = []
for y, x in peaks:
    win = dt[y - 4:y + 5, x - 4:x + 5]
    m = win >= win.max() - 0.6
    ys, xs = np.nonzero(m)
    centres.append((x + c0 - 4 + xs.mean(), y + r0 - 4 + ys.mean(), float(dt[y, x])))
centres.sort()
# merge duplicate maxima of the same disc (flat DT plateau): keep one per 6 px
merged = []
for c in centres:
    if merged and abs(c[0] - merged[-1][0]) < 6 and abs(c[1] - merged[-1][1]) < 6:
        a, b = merged[-1], c
        merged[-1] = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, max(a[2], b[2]))
    else:
        merged.append(c)
centres = merged
radius = float(np.median([c[2] for c in centres]))

pts = []
CAP_W = 3.0      # cap stroke width (px), measured on the clearly separated caps below
CAP_HALF = 6     # caps span about +-6.5 px around the stem (13 px wide)


def outer_edge(colset, cy, sgn, rad):
    """walk from the marker centre outward (sgn=-1 up, +1 down) through the
    connected dark region (disc + stem + cap) in the given columns; return the
    sub-pixel half-maximum position of its outer boundary."""
    prof_rows = np.arange(int(round(cy)), int(round(cy)) + sgn * 70, sgn)
    prof = K[prof_rows][:, colset].mean(axis=1)
    j = 0
    while prof[j + 1] > 127.5:
        j += 1
    # crossing between j (dark) and j+1 (light)
    frac = (prof[j] - 127.5) / (prof[j] - prof[j + 1])
    return float(prof_rows[j] + sgn * frac)


for cx, cy, rad in centres:
    rec = {"col": cx, "row": cy, "p": p_of(cx), "value": y_of(cy), "marker_radius_px": rad}
    ic = int(round(cx))
    on_axis = abs(cx - yaxis_col) < 3
    # stem columns: the 3 centre columns, except for the p = 1 marker whose stem
    # lies inside the y-axis line; there the caps are read in the 4 columns
    # left of the axis (inside the cap span, outside the axis stroke)
    cols = list(range(ic - 6, ic - 2)) if on_axis else [ic - 1, ic, ic + 1]
    # where the bare disc would end in these columns
    dx = max(abs(c - cx) for c in cols)
    disc_half = np.sqrt(max(rad ** 2 - dx ** 2, 0.0))
    for tag, sgn in (("top", -1), ("bottom", +1)):
        if on_axis:
            # here the cap is separated from the disc by a white gap: find the
            # first stroke beyond the disc edge and take its half-max midpoint
            rows_ = np.arange(int(round(cy)) + sgn * int(disc_half + 2), int(round(cy)) + sgn * 50, sgn)
            prof = K[rows_][:, cols].mean(axis=1)
            k0 = 0
            while prof[k0] > 60:
                k0 += 1
            hit = k0 + int(np.argmax(prof[k0:] > 100))
            rc = int(rows_[hit])
            lo, hi = rc - 4, rc + 5
            cap_c, cap_w = stroke_centre(K[lo:hi, cols].mean(axis=1), np.arange(lo, hi, dtype=float))
            rec[f"{tag}_cap_row"] = cap_c
            rec[f"{tag}_cap_stroke_px"] = cap_w
            rec[f"{tag}_beyond_disc_px"] = float(abs(cap_c - cy) - rad)
            continue
        e = outer_edge(cols, cy, sgn, rad)
        beyond = sgn * (e - cy) - (rad + 0.5 if not on_axis else disc_half + 0.5)
        rec[f"{tag}_outer_edge_row"] = e
        rec[f"{tag}_beyond_disc_px"] = float(beyond)
        if beyond > 1.0:
            rec[f"{tag}_cap_row"] = e - sgn * CAP_W / 2
        else:
            rec[f"{tag}_cap_row"] = None   # error bar hidden inside the marker disc
    if on_axis:
        rec["note"] = "marker on the y-axis line; caps read in the 4 columns left of the axis"
    rec["err_high"] = y_of(rec["top_cap_row"]) if rec["top_cap_row"] else None
    rec["err_low"] = y_of(rec["bottom_cap_row"]) if rec["bottom_cap_row"] else None
    pts.append(rec)

# ------------------------ comparison with Fig. 3I ------------------------
f3 = json.load(open(F3I))
bars = {b["wt_fraction"]: b for b in f3["bars"]}
R1, R0 = bars[1.0]["mean"], bars[0.0]["mean"]
for rec in pts:
    p_nom = min(bars, key=lambda w: abs(w - rec["p"]))
    b = bars[p_nom]
    rec["p_nominal"] = p_nom
    rec["fig3I_normalized_mean"] = (b["mean"] - R0) / (R1 - R0)
    rec["fig3I_normalized_sem"] = b["sem"] / (R1 - R0)
    rec["diff_vs_fig3I"] = rec["value"] - rec["fig3I_normalized_mean"]


# ------------------- model curves (recomputed) + checks ------------------
def models(p):
    q = 1 - p
    P4, P3, P2 = p ** 4, 4 * p ** 3 * q, 6 * p ** 2 * q ** 2
    return {"n>=4": P4, "n>=3": P4 + P3, "n>=2": P4 + P3 + P2, "n>=1": 1 - q ** 4,
            "n>=2* (paper: 3p^2q^2 for 2 adjacent WT)": P4 + P3 + 0.5 * P2,
            "n>=2 adjacent, 4 of 6 arrangements (4p^2q^2)": P4 + P3 + (4 / 6) * P2,
            "no co-assembly (linear)": p}


def rmsd(vals_by_p, key, ps):
    return float(np.sqrt(np.mean([(vals_by_p[p] - models(p)[key]) ** 2 for p in ps])))

ps_all = [1.0, 0.8, 0.6, 0.5, 0.4, 0.2, 0.16, 0.0]
ps_int = [0.8, 0.6, 0.5, 0.4, 0.2, 0.16]
norm3 = {p: (bars[p]["mean"] - R0) / (R1 - R0) for p in ps_all}
dig4 = {rec["p_nominal"]: rec["value"] for rec in pts}
paper = {"n>=4": 0.327, "n>=3": 0.124, "n>=2": 0.119, "n>=1": 0.352,
         "n>=2* (paper: 3p^2q^2 for 2 adjacent WT)": 0.031, "no co-assembly (linear)": 0.071}
rms_table = {}
for key in models(0.5):
    rms_table[key] = {"paper": paper.get(key),
                      "fig3I_norm_8pts": rmsd(norm3, key, ps_all), "fig3I_norm_6interior": rmsd(norm3, key, ps_int),
                      "fig4B_digitized_8pts": rmsd(dig4, key, ps_all), "fig4B_digitized_6interior": rmsd(dig4, key, ps_int)}

# spot-check the plotted 2* (solid black) and n>=2 (dashed black) curves: at a
# few p values well away from markers, find black curve pixels in the column
curve_checks = []
for p in (0.9, 0.7, 0.3, 0.1):
    col = int(round(col_of(p)))
    colk = K[740:1230, col]
    runs, inr = [], False
    for i, v in enumerate(colk):
        if v > 128 and not inr:
            s = i; inr = True
        if v <= 128 and inr:
            runs.append((740 + s, 740 + i - 1)); inr = False
    found = [y_of((a + b) / 2) for a, b in runs]
    curve_checks.append({"p": p_of(col), "black_ink_values_in_column": found,
                         "pred_2star": models(p_of(col))["n>=2* (paper: 3p^2q^2 for 2 adjacent WT)"],
                         "pred_n2": models(p_of(col))["n>=2"],
                         "pred_2adj_4of6": models(p_of(col))["n>=2 adjacent, 4 of 6 arrangements (4p^2q^2)"]})

out = {"source_image": SRC,
       "calibration": {"x_ticks": xt, "x_fit": [ax, bx], "x_residuals": x_res, "px_per_unit_x": float(-1 / bx),
                       "y_ticks": yt, "y_fit": [ay, by], "y_residuals": y_res, "px_per_unit_y": float(-1 / by),
                       "x_scale": "linear, reversed (p=1 left)", "y_scale": "linear"},
       "marker_radius_px": radius, "points": pts, "rmsd": rms_table, "curve_checks": curve_checks,
       "normalization_check": {"R(1)": R1, "R(0)": R0}}
json.dump(out, open(OUT, "w"), indent=1)

# -------------------------------- overlay --------------------------------
try:
    FONT = ImageFont.load_default(size=18)
except TypeError:
    FONT = ImageFont.load_default()
X0, Y0, X1, Y1 = 60, 690, 860, 1300
S = 2
cr = img.crop((X0, Y0, X1, Y1)).resize(((X1 - X0) * S, (Y1 - Y0) * S), Image.NEAREST)
d = ImageDraw.Draw(cr)
P = lambda x, y: ((x - X0 + 0.5) * S, (y - Y0 + 0.5) * S)
for p, c in xt:
    d.line([P(c, 1229), P(c, 1250)], fill=(0, 120, 255), width=2)
for v, r in yt:
    d.line([P(150, r), P(176, r)], fill=(0, 120, 255), width=2)
for rec in pts:
    x, y = P(rec["col"], rec["row"])
    d.ellipse([x - radius * S, y - radius * S, x + radius * S, y + radius * S], outline=(255, 0, 255), width=2)
    d.line([(x - 6, y), (x + 6, y)], fill=(255, 0, 255), width=2)
    for tag in ("top_cap_row", "bottom_cap_row"):
        if rec[tag]:
            d.line([P(rec["col"] - 12, rec[tag]), P(rec["col"] + 12, rec[tag])], fill=(0, 190, 0), width=2)
    d.text((x + 30, y - 10), f"{rec['value']:.3f}", fill=(200, 0, 160), font=FONT)
d.text((10, (Y1 - Y0) * S - 30), "magenta: located marker discs (centre = mean); green: error-bar caps; blue: ticks",
       fill=(0, 0, 0), font=FONT)
cr.save(OVL)

print("x ticks", [(p, round(c, 2)) for p, c in xt], "px/unit %.2f" % (-1 / bx), "resid", [round(r, 4) for r in x_res])
print("y ticks", [(v, round(r, 2)) for v, r in yt], "px/unit %.2f" % (-1 / by), "resid", [round(r, 4) for r in y_res])
print("marker radius %.2f px" % radius)
for rec in pts:
    eh = "NA" if rec["err_high"] is None else "%.4f" % rec["err_high"]
    el = "NA" if rec["err_low"] is None else "%.4f" % rec["err_low"]
    print("p=%.4f (nominal %.2f) value=%.4f [%s, %s]  Fig3I-normalized %.4f +- %.4f  diff %+.4f" % (
        rec["p"], rec["p_nominal"], rec["value"], el, eh, rec["fig3I_normalized_mean"],
        rec["fig3I_normalized_sem"], rec["diff_vs_fig3I"]))
print("RMSD (paper vs recomputed):")
for k, v in rms_table.items():
    print("  %-48s paper %-6s | F3I-norm 8pt %.3f 6pt %.3f | F4B 8pt %.3f 6pt %.3f" % (
        k, v["paper"], v["fig3I_norm_8pts"], v["fig3I_norm_6interior"], v["fig4B_digitized_8pts"], v["fig4B_digitized_6interior"]))
for c in curve_checks:
    print("curve check p=%.3f ink %s pred 2*=%.3f n2=%.3f 2adj(4/6)=%.3f" % (
        c["p"], [round(v, 3) for v in c["black_ink_values_in_column"]], c["pred_2star"], c["pred_n2"], c["pred_2adj_4of6"]))

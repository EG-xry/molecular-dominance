"""Digitize Hosy & Vivaudou 2014 Fig. 3I (bar chart, raster only).

Usage: python3 -I hosy_fig3I_digitize.py <fig3_png> <out_json> <overlay_png> <zoom_png>

Image: the lossless (Flate) 945 x 1719 px Figure 3 embedded in the publisher
PDF (pdfimages -all, page 4, image 2). Pixel coordinates are (col, row) with
the origin at the top-left pixel; a pixel's centre is at integer coordinates.

Method
  * blackness = 255 - max(R, G, B): black strokes ~255, white 0, the grey bar
    fill (192,192,192) 63, red label text ~0.
  * Every stroke position (tick, bar outline, error-bar cap) is the midpoint of
    its two half-maximum edges, each edge measured against the background on
    its own side (white above a bar top, grey fill below it). This is robust to
    the fill and to anti-aliasing; a plain centroid is reported as a check.
  * y calibration: 7 tick marks (0 ... 1.2) left of the axis line; least-
    squares line row -> value; residuals reported. Checks: the dashed y = 1
    reference line and the bottom outlines of all 8 bars (which must sit at 0).
  * bar mean: top outline, measured in interior columns that avoid the side
    outlines and the error-bar cap.
  * error bar: only the upper half (stem + cap) is drawn. Cap row measured in
    cap columns that avoid the stem. For the 0::1 bar the cap overlaps the bar
    outline; its centre is then the cap's upper edge + half the cap stroke width
    measured on the other seven caps.
"""
import json
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

SRC, OUT, OVL, ZOOM = sys.argv[1:5]
img = Image.open(SRC).convert("RGB")
a = np.asarray(img).astype(float)
K = 255.0 - a.max(axis=2)


def stroke_centre(prof, coords):
    """midpoint of the two half-maximum edges of the single stroke in prof."""
    prof = np.asarray(prof, float)
    i = int(np.argmax(prof))
    peak = prof[i]
    b0 = float(np.median(prof[:3]))
    b1 = float(np.median(prof[-3:]))
    h0 = (peak + b0) / 2.0
    h1 = (peak + b1) / 2.0
    j = i
    while prof[j - 1] > h0:
        j -= 1
    e0 = coords[j - 1] + (h0 - prof[j - 1]) / (prof[j] - prof[j - 1])
    j = i
    while prof[j + 1] > h1:
        j += 1
    e1 = coords[j] + (prof[j] - h1) / (prof[j] - prof[j + 1])
    return float((e0 + e1) / 2.0), float(e1 - e0), float(e0), float(e1)


def centroid(prof, coords, bg=0.0, thresh=20.0):
    w = np.asarray(prof, float) - bg
    w = np.where(w > thresh, w, 0.0)
    return float((w * coords).sum() / w.sum())


def row_profile(r0, r1, cols):
    return K[r0:r1 + 1, cols].mean(axis=1), np.arange(r0, r1 + 1, dtype=float)


# ------------------------------ calibration ------------------------------
TICK_COLS = list(range(197, 209))   # tick strokes, left of the axis line (cols 210-211)
TICKS = [(0.0, 1555), (0.2, 1478), (0.4, 1402), (0.6, 1325), (0.8, 1248), (1.0, 1172), (1.2, 1095)]
cal_ticks = []
for val, r in TICKS:
    p, rr = row_profile(r - 6, r + 6, TICK_COLS)
    c, wdt, e0, e1 = stroke_centre(p, rr)
    cal_ticks.append({"value": val, "row": c, "stroke_px": wdt, "row_centroid_check": centroid(p, rr)})
rows = np.array([t["row"] for t in cal_ticks])
vals = np.array([t["value"] for t in cal_ticks])
B, A = np.polyfit(rows, vals, 1)
resid = vals - (A + B * rows)
for t, r in zip(cal_ticks, resid):
    t["fit_residual_units"] = float(r)
to_val = lambda row: float(A + B * row)
PX = -1.0 / B

# checks: dashed y=1 line (between bars) and bar bottoms
checks = {"dashed_line_y1": []}
for c0, c1 in ((380, 405), (470, 498), (560, 590)):
    # dashed: take, for each row, the max over the columns (dash or gap)
    p = K[1162:1182, c0:c1].max(axis=1)
    rr = np.arange(1162, 1182, dtype=float)
    c, w, _, _ = stroke_centre(p, rr)
    checks["dashed_line_y1"].append({"cols": [c0, c1], "row": c, "value": to_val(c)})

BARS = [(231, 233, 285, 287), (322, 324, 375, 377), (409, 412, 464, 467), (501, 504, 556, 559),
        (593, 595, 649, 651), (689, 691, 741, 743), (776, 778, 830, 832), (865, 867, 920, 922)]
LABELS = ["1::0", "0.8::0.2", "0.6::0.4", "0.5::0.5", "0.4::0.6", "0.2::0.8", "0.16::0.84", "0::1"]
WT = [1.0, 0.8, 0.6, 0.5, 0.4, 0.2, 0.16, 0.0]
N = [87, 10, 6, 19, 10, 12, 13, 38]

checks["bar_bottoms"] = []
for i, (l0, l1, r0, r1) in enumerate(BARS):
    cols = list(range(l1 + 4, r0 - 3))
    p, rr = row_profile(1546, 1564, cols)
    c, w, _, _ = stroke_centre(p, rr)
    checks["bar_bottoms"].append({"bar": i + 1, "row": c, "value": to_val(c), "stroke_px": w})

# --------------------------------- bars ----------------------------------
results = []
for i, (l0, l1, r0, r1) in enumerate(BARS):
    inner = list(range(l1 + 4, l1 + 13)) + list(range(r0 - 12, r0 - 3))
    colprof = K[1080:1560, inner].mean(axis=1)
    r = 1545 - 1080
    while colprof[r] < 150:          # walk up through the grey fill to the top outline
        r -= 1
    rt = r + 1080
    p, rr = row_profile(rt - 8, rt + 9, inner)
    top, top_w, te0, te1 = stroke_centre(p, rr)
    rec = {"bar": i + 1, "label": LABELS[i], "wt_fraction": WT[i], "n": N[i],
           "left_outline_cols": [l0, l1], "right_outline_cols": [r0, r1], "inner_cols": [inner[0], inner[-1]],
           "bar_top_row": top, "bar_top_stroke_px": top_w, "bar_top_upper_edge": te0,
           "mean": to_val(top)}
    # stem: darkest column just above the outline, between the side outlines
    above = K[int(round(top)) - 9:int(round(top)) - 3, l1 + 3:r0 - 2].mean(axis=0)
    stem_c = l1 + 3 + int(np.argmax(above))
    rec["stem_col"] = stem_c
    rec["stem_strength"] = float(above.max())
    if above.max() > 120:
        stem_cols = [c for c in range(stem_c - 1, stem_c + 2) if K[int(top) - 6, c] > 120]
        rs = int(round(top)) - 4
        while K[rs - 1, stem_cols].max() > 120:
            rs -= 1
        capcols = [c for c in range(stem_c - 9, stem_c + 10) if abs(c - stem_c) >= 3]
        p, rr = row_profile(rs - 6, rs + 7, capcols)
        cap, cap_w, ce0, ce1 = stroke_centre(p, rr)
        rec.update({"cap_row": cap, "cap_stroke_px": cap_w, "cap_upper_edge": ce0, "stem_top_row": rs,
                    "cap_cols": [capcols[0], capcols[-1]], "cap_method": "half-max edge midpoint of cap stroke"})
    else:
        rec.update({"cap_row": None, "cap_method": "pending"})
    results.append(rec)

cap_w_mean = float(np.mean([r["cap_stroke_px"] for r in results if r.get("cap_stroke_px")]))
cap_w_sd = float(np.std([r["cap_stroke_px"] for r in results if r.get("cap_stroke_px")]))
for rec in results:
    if rec["cap_row"] is not None:
        continue
    i = rec["bar"] - 1
    l0, l1, r0, r1 = BARS[i]
    bt = int(round(rec["bar_top_row"]))
    ue = {}
    for c in range(l1 + 3, r0 - 2):
        col = K[bt - 10:bt + 1, c]       # white above, stroke at the bottom of the window
        pk = col.max()
        k = int(np.argmax(col > pk / 2))
        ue[c] = bt - 10 + (k - 1) + (pk / 2 - col[k - 1]) / (col[k] - col[k - 1])
    cols = np.array(sorted(ue))
    v = np.array([ue[c] for c in cols])
    base = float(np.median(v))
    bump = cols[v < base - 0.8]
    edge = float(np.median(v[v < base - 0.8]))
    rec.update({"bump_cols": [int(bump.min()), int(bump.max())], "outline_upper_edge_median": base,
                "cap_upper_edge": edge, "cap_row": edge + cap_w_mean / 2,
                "cap_method": ("cap overlaps the bar outline (visible only as a bump above it): cap centre = "
                               f"bump upper edge + mean cap stroke/2 ({cap_w_mean:.2f}/2 px)")})

for rec in results:
    rec["err_high"] = to_val(rec["cap_row"])
    rec["sem"] = rec["err_high"] - rec["mean"]
    rec["err_low_by_symmetry"] = rec["mean"] - rec["sem"]

out = {"source_image": SRC, "image_size": list(img.size),
       "calibration": {"ticks": cal_ticks, "value_per_row": B, "value_at_row0": A, "px_per_unit": PX,
                       "max_abs_residual_units": float(np.abs(resid).max()), "scale": "linear"},
       "checks": checks, "mean_cap_stroke_px": cap_w_mean, "sd_cap_stroke_px": cap_w_sd, "bars": results}
json.dump(out, open(OUT, "w"), indent=1)

# ------------------------------- overlays --------------------------------
try:
    FONT = ImageFont.load_default(size=22)
    FONT_S = ImageFont.load_default(size=16)
except TypeError:
    FONT = FONT_S = ImageFont.load_default()
x0c, y0c, x1c, y1c = 120, 1080, 945, 1600
S = 3
crop = img.crop((x0c, y0c, x1c, y1c)).resize(((x1c - x0c) * S, (y1c - y0c) * S), Image.NEAREST)
d = ImageDraw.Draw(crop)
P = lambda x, y: ((x - x0c + 0.5) * S, (y - y0c + 0.5) * S)
for t in cal_ticks:
    d.line([P(186, t["row"]), P(222, t["row"])], fill=(0, 120, 255), width=2)
    xx, yy = P(122, t["row"])
    d.text((xx, yy - 30), f"tick {t['value']:.1f} @ row {t['row']:.2f}", fill=(0, 90, 220), font=FONT_S)
for rec in results:
    l0, l1, r0, r1 = BARS[rec["bar"] - 1]
    d.line([P(l0 - 5, rec["bar_top_row"]), P(r1 + 5, rec["bar_top_row"])], fill=(230, 0, 0), width=2)
    d.line([P(rec["stem_col"] - 14, rec["cap_row"]), P(rec["stem_col"] + 14, rec["cap_row"])], fill=(0, 170, 0), width=2)
    xx, yy = P(l0 + 3, rec["bar_top_row"] + 6)
    d.text((xx, yy), f"{rec['mean']:.3f}\n+{rec['sem']:.3f}", fill=(200, 0, 0), font=FONT)
d.text((12, (y1c - y0c) * S - 34),
       "blue: tick centres used for calibration | red: bar-top outline = mean | green: cap = mean + SEM (upper half only drawn)",
       fill=(0, 0, 0), font=FONT_S)
crop.save(OVL)

# zoom panel for the two smallest error bars (bars 7, 8)
tiles = []
for b in (7, 8):
    rec = results[b - 1]
    l0, l1, r0, r1 = BARS[b - 1]
    y = int(rec["bar_top_row"])
    zx0, zy0, zx1, zy1 = l0 - 6, y - 30, r1 + 6, y + 12
    Z = 10
    t = img.crop((zx0, zy0, zx1, zy1)).resize(((zx1 - zx0) * Z, (zy1 - zy0) * Z), Image.NEAREST)
    dz = ImageDraw.Draw(t)
    Q = lambda x, yy: ((x - zx0 + 0.5) * Z, (yy - zy0 + 0.5) * Z)
    dz.line([Q(zx0, rec["bar_top_row"]), Q(zx1, rec["bar_top_row"])], fill=(255, 0, 0), width=2)
    dz.line([Q(zx0, rec["cap_row"]), Q(zx1, rec["cap_row"])], fill=(0, 190, 0), width=2)
    dz.text((8, 8), f"bar {b} ({rec['label']}): mean {rec['mean']:.4f}, cap {rec['err_high']:.4f}",
            fill=(0, 0, 0), font=FONT)
    tiles.append(t)
Wt = max(t.width for t in tiles)
canvas = Image.new("RGB", (Wt, sum(t.height for t in tiles) + 10), "white")
yy = 0
for t in tiles:
    canvas.paste(t, (0, yy))
    yy += t.height + 10
canvas.save(ZOOM)

print(f"calibration: value = {A:.6f} + {B:.8f}*row ; {PX:.3f} px/unit ; max |resid| {np.abs(resid).max():.5f}")
for t in cal_ticks:
    print(f"  tick {t['value']:.1f}: row {t['row']:.3f} (centroid {t['row_centroid_check']:.3f}) "
          f"stroke {t['stroke_px']:.2f}px resid {t['fit_residual_units']:+.5f}")
print("dashed y=1 line:", [f"{c['value']:.4f}" for c in checks["dashed_line_y1"]])
print("bar bottoms:", [f"{c['value']:+.4f}" for c in checks["bar_bottoms"]])
print(f"cap stroke {cap_w_mean:.2f} +- {cap_w_sd:.2f} px")
for rec in results:
    print(f"bar {rec['bar']} {rec['label']:10s} n={rec['n']:3d} top_row={rec['bar_top_row']:.2f} "
          f"(stroke {rec['bar_top_stroke_px']:.2f}) mean={rec['mean']:.4f} cap_row={rec['cap_row']:.2f} "
          f"err_high={rec['err_high']:.4f} sem={rec['sem']:.4f} stem_col={rec['stem_col']}")

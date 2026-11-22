"""Cross-check of the Hosy & Vivaudou 2014 Fig. 3I readings against two other
renderings of the same data (same n): Fig. 1A (300 uM ADP bars: WT n=87,
DNN n=38) and Fig. 2D (WT 1::0 n=87, 0.5::0.5 n=19, 0::1 n=38).

Usage: python3 -I hosy_crosscheck_fig1A_fig2D.py <fig1.jpg> <fig2.jpg> <out_json> <crop_dir>

Images: JPEGs embedded in the publisher PDF (pdfimages -all): page 3 image 0
= Fig. 1 (1580x1085); page 4 image 1 = Fig. 2 (967x1140).
blackness = 255 - max(R, G, B) so coloured fills (light blue, red, grey) are
not counted as ink. Stroke positions = midpoint of the two half-maximum edges
(same estimator as hosy_fig3I_digitize.py).
"""
import json
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

F1, F2, OUT, CROPS = sys.argv[1:5]


def blackness(img):
    a = np.asarray(img.convert("RGB")).astype(float)
    return 255.0 - a.max(axis=2)


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
    return float((e0 + e1) / 2.0), float(e1 - e0), float(e0), float(e1)


def upper_edge(col, r0):
    """half-max upper edge of the first stroke in a column segment (white above)."""
    pk = col.max()
    j = int(np.argmax(col > pk / 2))
    return r0 + (j - 1) + (pk / 2 - col[j - 1]) / (col[j] - col[j - 1])


try:
    FONT = ImageFont.load_default(size=16)
except TypeError:
    FONT = ImageFont.load_default()
out = {}

# =============================== Fig. 1A ==================================
im1 = Image.open(F1)
k1 = blackness(im1)
tick_vals = [2.5, 2.0, 1.5, 1.0, 0.5, 0.0]
tick_rows = []
for r in (23, 133, 245, 355, 465, 576):
    rr = np.arange(r - 6, r + 7, dtype=float)
    tick_rows.append(stroke_centre(k1[r - 6:r + 7, 286:299].mean(axis=1), rr)[0])
b1, a1 = np.polyfit(tick_rows, tick_vals, 1)
res1 = np.array(tick_vals) - (a1 + b1 * np.array(tick_rows))
v1 = lambda row: float(a1 + b1 * row)
fig1 = {"ticks": [[v, r] for v, r in zip(tick_vals, tick_rows)], "px_per_unit": float(-1 / b1),
        "max_abs_residual_units": float(np.abs(res1).max())}
# dashed reference line y = 1, between the 300 uM DNN bar and the 1 mM bar
rr = np.arange(345, 366, dtype=float)
dl = stroke_centre(k1[345:366, 830:870].max(axis=1), rr)[0]
fig1["dashed_line_y1_value"] = v1(dl)
bars1 = {}
# WT 300 uM (light blue): side outlines cols 747-749 / 786-788, stem near col 766
cols = list(range(752, 760)) + list(range(777, 783))
rr = np.arange(338, 354, dtype=float)
top, w, _, _ = stroke_centre(k1[338:354, cols].mean(axis=1), rr)
rr = np.arange(326, 343, dtype=float)
capcols = [c for c in range(755, 778) if abs(c - 766) >= 3]
cap, cw, _, _ = stroke_centre(k1[326:343, capcols].mean(axis=1), rr)
bars1["WT 300 uM (n=87)"] = {"top_row": top, "cap_row": cap, "mean": v1(top), "err_high": v1(cap),
                             "sem": v1(cap) - v1(top), "cap_stroke_px": cw}
# DNN 300 uM (red): cap overlaps the top outline. Outline-only columns are the
# few columns at both ends of the bar outside the cap (790-792, 821-822); cap
# columns 795-818. Compare half-max upper edges; the outline centre comes from
# the outline-only columns.
ue_out = [upper_edge(k1[554:567, c], 554) for c in (790, 791, 792, 821, 822)]
ue_cap = [upper_edge(k1[554:567, c], 554) for c in range(795, 819)]
rr = np.arange(556, 575, dtype=float)
top_d = []
for c in (790, 791, 792):          # left end only: the right end is 1-2 px from the side outline
    top_d.append(stroke_centre(k1[556:575, c], rr)[0])
top_d = float(np.mean(top_d))
cap_d = float(np.median(ue_cap)) + cw / 2       # cap centre from its upper edge + half stroke
bars1["DNN 300 uM (n=38)"] = {"top_row_outline_only_cols": top_d, "outline_upper_edge": float(np.median(ue_out)),
                              "cap_upper_edge": float(np.median(ue_cap)), "cap_row_est": cap_d,
                              "mean": v1(top_d), "err_high": v1(cap_d), "sem": v1(cap_d) - v1(top_d),
                              "sem_from_edge_offset": (float(np.median(ue_out)) - float(np.median(ue_cap))) / float(-1 / b1)}
fig1["bars"] = bars1
out["fig1A"] = fig1

# =============================== Fig. 2D ==================================
im2 = Image.open(F2)
k2 = blackness(im2)
tick_vals2 = [0.0, 0.5, 1.0]
tick_cols2 = []
for c in (224, 409, 579):
    cc = np.arange(c - 6, c + 7, dtype=float)
    tick_cols2.append(stroke_centre(k2[946:958, c - 6:c + 7].mean(axis=0), cc)[0])
sp1 = tick_cols2[1] - tick_cols2[0]; sp2 = tick_cols2[2] - tick_cols2[1]
px2 = tick_cols2[2] - tick_cols2[0]
v2 = lambda col: float((col - tick_cols2[0]) / px2)          # 0 and 1 ticks
B3, A3 = np.polyfit(tick_cols2, tick_vals2, 1)
v3 = lambda col: float(A3 + B3 * col)                        # all three ticks
cc = np.arange(218, 234, dtype=float)
base = stroke_centre(k2[986:1020, 218:234].mean(axis=0), cc)[0]
fig2 = {"ticks": [[v, c] for v, c in zip(tick_vals2, tick_cols2)], "spacing_0_to_0.5_px": sp1,
        "spacing_0.5_to_1_px": sp2, "px_per_unit_0_1_ticks": px2,
        "bar_base_col": base, "bar_base_value_0_1_ticks": v2(base),
        "where_0.5_tick_should_be_if_linear": (tick_cols2[0] + tick_cols2[2]) / 2}
spec = {"1::0 (n=87)": ((986, 1020), (592, 607), (608, 620), (1001, 1004)),
        "0.5::0.5 (n=19)": ((1040, 1074), (415, 429), (428, 440), (1055, 1058)),
        "0::1 (n=38)": ((1092, 1127), (234, 244), (243, 251), (1109, 1111))}
bars2 = {}
for name, ((r0, r1), (o0, o1), (q0, q1), (s0, s1)) in spec.items():
    rows_ol = [r for r in range(r0, r1 + 1) if not (s0 - 9 <= r <= s1 + 9)]
    cc = np.arange(o0, o1 + 1, dtype=float)
    ol, olw, _, _ = stroke_centre(k2[rows_ol, o0:o1 + 1].mean(axis=0), cc)
    rows_cap = [r for r in range(s0 - 5, s1 + 6) if not (s0 - 1 <= r <= s1 + 1)]
    cc = np.arange(q0, q1 + 1, dtype=float)
    cap, capw, _, _ = stroke_centre(k2[rows_cap, q0:q1 + 1].mean(axis=0), cc)
    bars2[name] = {"outline_col": ol, "cap_col": cap, "mean": v2(ol), "err_high": v2(cap),
                   "sem": v2(cap) - v2(ol), "mean_3tick_fit": v3(ol), "sem_3tick_fit": v3(cap) - v3(ol)}
fig2["bars"] = bars2
out["fig2D"] = fig2
json.dump(out, open(OUT, "w"), indent=1)

# ================================ crops ===================================
S = 2
c1 = im1.convert("RGB").crop((280, 0, 960, 600)).resize((680 * S, 600 * S), Image.NEAREST)
d = ImageDraw.Draw(c1)
P = lambda x, y: ((x - 280 + 0.5) * S, (y + 0.5) * S)
for v, r in zip(tick_vals, tick_rows):
    d.line([P(282, r), P(312, r)], fill=(0, 120, 255), width=2)
b = bars1["WT 300 uM (n=87)"]
d.line([P(742, b["top_row"]), P(793, b["top_row"])], fill=(255, 0, 0), width=2)
d.line([P(752, b["cap_row"]), P(781, b["cap_row"])], fill=(0, 170, 0), width=2)
b = bars1["DNN 300 uM (n=38)"]
d.line([P(786, b["top_row_outline_only_cols"]), P(830, b["top_row_outline_only_cols"])], fill=(255, 0, 0), width=2)
d.line([P(792, b["cap_row_est"]), P(822, b["cap_row_est"])], fill=(0, 170, 0), width=2)
d.text((10, 10), "Fig 1A cross-check: blue ticks; red bar tops; green caps", fill=(0, 0, 0), font=FONT)
c1.save(f"{CROPS}/hosy2014_fig1A_crosscheck_overlay.png")
S = 3
c2 = im2.convert("RGB").crop((200, 930, 700, 1140)).resize((500 * S, 210 * S), Image.NEAREST)
d = ImageDraw.Draw(c2)
P = lambda x, y: ((x - 200 + 0.5) * S, (y - 930 + 0.5) * S)
for v, c in zip(tick_vals2, tick_cols2):
    d.line([P(c, 938), P(c, 963)], fill=(0, 120, 255), width=2)
lin = (tick_cols2[0] + tick_cols2[2]) / 2
d.line([P(lin, 936), P(lin, 946)], fill=(255, 0, 255), width=2)
for name, ((r0, r1), _, _, (s0, s1)) in spec.items():
    b = bars2[name]
    d.line([P(b["outline_col"], r0 - 4), P(b["outline_col"], r1 + 4)], fill=(255, 0, 0), width=2)
    d.line([P(b["cap_col"], s0 - 9), P(b["cap_col"], s1 + 9)], fill=(0, 170, 0), width=2)
d.text((8, 600), "Fig 2D: blue ticks; magenta = linear position of 0.5 from 0,1 ticks; red bar ends; green caps",
       fill=(0, 0, 0), font=FONT)
c2.save(f"{CROPS}/hosy2014_fig2D_crosscheck_overlay.png")

print(json.dumps(out, indent=1))

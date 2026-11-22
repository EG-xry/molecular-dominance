"""Overlay the vector-extracted Clatot 2017 points on 600-dpi renderings of the
PDF pages, and convert the calibration ticks to pixel coordinates of those
renderings.

Usage: python3 -I clatot_overlay.py <points_json> <fig1c_png> <figS3C_png> <crop_dir> <out_json>

The renderings were made with
  pdftoppm -r 600 -f 2 -l 2 -x 2458 -y 3042 -W 1417 -H 1125   (main PDF, Fig. 1c)
  pdftoppm -r 600 -f 3 -l 3 -x 4125 -y 2500 -W 1750 -H 1375   (SI PDF, Supp. Fig. 3C)
so pixel = pt * 600/72 - crop origin.
"""
import json
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

PTS, F1C, FS3C, CROPS, OUT = sys.argv[1:6]
S = 600.0 / 72.0
ORIGIN = {"fig1c": (2458, 3042), "figS3C": (4125, 2500)}
SERIES = {"fig1c": {"circle": "WT + L325R", "triangle": "WT + L325R + difopein", "square": "S460A + L325R"},
          "figS3C": {"circle": "WT + R104W", "square": "S460A + R104W", "triangle": "WT + R878C"}}
COL = {"circle": (220, 0, 0), "triangle": (0, 120, 255), "square": (0, 160, 0)}
try:
    FONT = ImageFont.load_default(size=22)
    FONT_S = ImageFont.load_default(size=18)
except TypeError:
    FONT = FONT_S = ImageFont.load_default()

data = json.load(open(PTS))
pix_cal = {}
for key, png in (("fig1c", F1C), ("figS3C", FS3C)):
    im = Image.open(png).convert("RGB")
    d = ImageDraw.Draw(im)
    ox, oy = ORIGIN[key]
    P = lambda x, y: (x * S - ox, y * S - oy)
    cal = data[key]["calibration"]
    # calibration ticks
    xt_px = []
    for x_pt, v in cal["xticks_pt"]:
        x, y = P(x_pt, 0)
        ybase = P(0, cal["yticks_pt"][0][0])[1]
        d.line([(x, ybase + 14), (x, ybase + 40)], fill=(255, 0, 255), width=3)
        xt_px.append([v, round(x, 2), round(x_pt, 4)])
    yt_px = []
    for y_pt, v in cal["yticks_pt"]:
        x0 = P(cal["xticks_pt"][0][0], y_pt)[0]
        x, y = x0, P(0, y_pt)[1]
        d.line([(x - 45, y), (x - 15, y)], fill=(255, 0, 255), width=3)
        yt_px.append([v, round(y, 2), round(y_pt, 4)])
    pix_cal[key] = {"render": png, "dpi": 600, "crop_origin_px": [ox, oy],
                    "x_ticks_value_px_pt": xt_px, "y_ticks_value_px_pt": yt_px,
                    "px_per_x_unit": cal["pt_per_x_unit"] * S, "px_per_y_unit": cal["pt_per_y_unit"] * S}
    # recomputed curves for comparison: (1-x)^3 (the Methods' trimer formula)
    xs = np.linspace(0, 1, 200)
    def to_px(xv, yv):
        x_pt = (xv - cal["x_a"]) / cal["x_b"]
        y_pt = (yv - cal["y_a"]) / cal["y_b"]
        return P(x_pt, y_pt)
    pts3 = [to_px(x, (1 - x) ** 3) for x in xs]
    for i in range(0, len(pts3) - 1, 2):
        d.line([pts3[i], pts3[i + 1]], fill=(255, 140, 0), width=2)
    # points
    for r in data[key]["points"]:
        cx, cy = P(r["cx"], r["junction_pt"] if r.get("junction_pt") else r["cy"])
        c = COL[r["kind"]]
        d.ellipse([cx - 13, cy - 13, cx + 13, cy + 13], outline=c, width=3)
        d.line([(cx - 20, cy), (cx + 20, cy)], fill=c, width=2)
        if r.get("top_pt") is not None:
            for yy in (r["top_pt"], r["bottom_pt"]):
                _, py = P(0, yy)
                d.line([(cx - 16, py), (cx + 16, py)], fill=(255, 0, 255), width=2)
    # legend in a white band added below the rendering
    band = 30 * (len(SERIES[key]) + 3)
    canvas = Image.new("RGB", (im.width, im.height + band), "white")
    canvas.paste(im, (0, 0))
    im = canvas
    d = ImageDraw.Draw(im)
    lx, ly = 20, im.height - band + 8
    for i, (k, name) in enumerate(SERIES[key].items()):
        d.text((lx, ly + 30 * i), f"{k}: {name}", fill=COL[k], font=FONT)
    d.text((lx, ly + 30 * len(SERIES[key])), "magenta ticks: calibration ticks; magenta bars: located error-bar ends",
           fill=(160, 0, 160), font=FONT_S)
    d.text((lx, ly + 30 * len(SERIES[key]) + 26), "orange dashes: (1-x)^3 recomputed (Methods 'trimer'); compare with plotted dash-dot curve",
           fill=(200, 100, 0), font=FONT_S)
    im.save(f"{CROPS}/clatot2017_{key}_overlay.png")

json.dump(pix_cal, open(OUT, "w"), indent=1)
for k, v in pix_cal.items():
    print(k, "x ticks (value, px, pt):", v["x_ticks_value_px_pt"])
    print(k, "y ticks (value, px, pt):", v["y_ticks_value_px_pt"])
    print(k, "px per x unit %.2f, px per y unit %.2f" % (v["px_per_x_unit"], v["px_per_y_unit"]))

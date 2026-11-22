"""Read Clatot et al. 2017 Fig. 1c and Supplementary Fig. 3C from the vector
drawing in the PDFs (paths exported by svg_paths.py).

Usage: python3 -I clatot_vector_extract.py <work_dir> <out_json>

Coordinates are PDF points in a y-down page frame, as pdftocairo writes them.
For each panel the script:
  1. calibrates both axes from the drawn tick marks (the labels are checked
     separately against pdftotext -bbox word boxes, see report);
  2. finds the markers by shape: filled circle (Bezier outline), filled square
     (4-vertex rectangle), open triangle (3 vertices, white fill);
  3. finds the error bars: the vertical segments on the marker's x that
     meet at the mean, and the horizontal caps at their far ends;
  4. converts everything to data units.
"""
import json
import sys
from collections import defaultdict

WORK = sys.argv[1]
OUT = sys.argv[2]

BLACK = ("rgb(0%, 0%, 0%)",)
WHITE = ("rgb(100%, 100%, 100%)",)


def load(fn):
    return json.load(open(f"{WORK}/{fn}"))


def verts(p):
    return [(a, b) for sp in p["subpaths"] for a, b, c in sp["pts"] if c == 0]


def is_seg(p):
    v = verts(p)
    return len(p["subpaths"]) == 1 and len(v) == 2


def calibrate(paths, axis_y_line, axis_x_line, xticks_val, yticks_val, tol=0.05):
    """axis_y_line: y of the x-axis line; axis_x_line: x of the y-axis line."""
    xt, yt = [], []
    for p in paths:
        if p["style"].get("stroke") not in BLACK or not is_seg(p):
            continue
        (xa, ya), (xb, yb) = verts(p)
        # x-axis tick: vertical, starting on the x-axis line, pointing down
        if abs(xa - xb) < tol and min(ya, yb) == min(ya, yb) and abs(min(ya, yb) - axis_y_line) < tol \
                and 0.5 < abs(ya - yb) < 3:
            xt.append(round(xa, 4))
        # y-axis tick: horizontal, ending on the y-axis line, pointing left
        if abs(ya - yb) < tol and abs(max(xa, xb) - axis_x_line) < tol and 0.5 < abs(xa - xb) < 3:
            yt.append(round(ya, 4))
    xt = sorted(set(xt))
    yt = sorted(set(yt), reverse=True)  # bottom (largest y) first
    assert len(xt) == len(xticks_val), (xt, xticks_val)
    assert len(yt) == len(yticks_val), (yt, yticks_val)
    # least-squares linear fit value = a + b * coord
    def fit(c, v):
        n = len(c); mc = sum(c) / n; mv = sum(v) / n
        b = sum((ci - mc) * (vi - mv) for ci, vi in zip(c, v)) / sum((ci - mc) ** 2 for ci in c)
        a = mv - b * mc
        res = max(abs(a + b * ci - vi) for ci, vi in zip(c, v))
        return a, b, res
    ax, bx, rx = fit(xt, xticks_val)
    ay, by, ry = fit(yt, yticks_val)
    return {
        "xticks_pt": list(zip(xt, xticks_val)),
        "yticks_pt": list(zip(yt, yticks_val)),
        "x_a": ax, "x_b": bx, "x_max_tick_residual": rx,
        "y_a": ay, "y_b": by, "y_max_tick_residual": ry,
        "pt_per_x_unit": 1 / bx, "pt_per_y_unit": -1 / by,
    }


def classify_marker(p):
    st = p["style"]
    v = verts(p)
    nctrl = sum(1 for sp in p["subpaths"] for a, b, c in sp["pts"] if c == 1)
    w = p["bbox"][2] - p["bbox"][0]
    h = p["bbox"][3] - p["bbox"][1]
    if not (1.5 < w < 6.5 and 1.5 < h < 6.5):
        return None
    fill = st.get("fill")
    if fill in BLACK and nctrl >= 8:
        return "circle"
    uniq = []
    for q in v:
        if all(abs(q[0] - u[0]) > 1e-3 or abs(q[1] - u[1]) > 1e-3 for u in uniq):
            uniq.append(q)
    if fill in BLACK and nctrl == 0 and len(uniq) == 4:
        return "square"
    if fill in WHITE and nctrl == 0 and len(uniq) == 3:
        return "triangle"
    return None


def marker_centre(p, kind):
    v = verts(p)
    uniq = []
    for q in v:
        if all(abs(q[0] - u[0]) > 1e-3 or abs(q[1] - u[1]) > 1e-3 for u in uniq):
            uniq.append(q)
    if kind == "triangle":  # centroid of the 3 vertices
        return (sum(q[0] for q in uniq) / 3, sum(q[1] for q in uniq) / 3)
    bb = p["bbox"]
    return ((bb[0] + bb[2]) / 2, (bb[1] + bb[3]) / 2)


def extract(paths, cal, plot_box, legend_box, tol=0.06):
    x0, y0, x1, y1 = plot_box
    lx0, ly0, lx1, ly1 = legend_box
    markers = []
    seen = set()
    for p in paths:
        kind = classify_marker(p)
        if not kind:
            continue
        cx, cy = marker_centre(p, kind)
        if not (x0 - 3 <= cx <= x1 + 3 and y0 - 3 <= cy <= y1 + 3):
            continue
        if lx0 <= cx <= lx1 and ly0 <= cy <= ly1:
            continue  # legend key symbols (drawn inside the axes box)
        key = (kind, round(cx, 1), round(cy, 1))
        if key in seen:
            continue
        seen.add(key)
        markers.append({"kind": kind, "cx": cx, "cy": cy, "bbox": p["bbox"], "pid": p["id"],
                        "w": p["bbox"][2] - p["bbox"][0], "h": p["bbox"][3] - p["bbox"][1]})
    segs = [p for p in paths if p["style"].get("stroke") in BLACK and is_seg(p)]
    vsegs = []
    hsegs = []
    for p in segs:
        (xa, ya), (xb, yb) = verts(p)
        if abs(xa - xb) < 1e-3 and abs(ya - yb) > 0.2:
            vsegs.append({"x": xa, "ya": min(ya, yb), "yb": max(ya, yb), "pid": p["id"]})
        elif abs(ya - yb) < 1e-3 and 1.0 < abs(xa - xb) < 4.0:
            hsegs.append({"y": ya, "xa": min(xa, xb), "xb": max(xa, xb), "pid": p["id"]})
    out = []
    for m in markers:
        # error-bar halves are drawn just before the marker in the content
        # stream; take the closest preceding pair on the same x that meet.
        cand = [s for s in vsegs if abs(s["x"] - m["cx"]) < tol and s["pid"] < m["pid"]
                and m["pid"] - s["pid"] <= 6]
        # an error bar = an upper half whose lower end coincides with the
        # upper end of a lower half (the shared point is the plotted mean);
        # choose the pair whose junction is closest to the marker centre.
        pairs = [(u, d) for u in cand for d in cand
                 if u is not d and abs(u["yb"] - d["ya"]) < 0.01]
        rec = dict(m)
        rec["x_val"] = cal["x_a"] + cal["x_b"] * m["cx"]
        rec["y_marker_val"] = cal["y_a"] + cal["y_b"] * m["cy"]
        if pairs:
            u, d = min(pairs, key=lambda ud: abs((ud[0]["yb"] + ud[1]["ya"]) / 2 - m["cy"]))
            junction = (u["yb"] + d["ya"]) / 2
            rec["junction_pt"] = junction
            rec["junction_mismatch_pt"] = abs(u["yb"] - d["ya"])
            rec["top_pt"] = u["ya"]; rec["bottom_pt"] = d["yb"]
            capt = [h for h in hsegs if abs(h["y"] - u["ya"]) < 1e-3 and h["xa"] - tol <= m["cx"] <= h["xb"] + tol]
            capb = [h for h in hsegs if abs(h["y"] - d["yb"]) < 1e-3 and h["xa"] - tol <= m["cx"] <= h["xb"] + tol]
            rec["top_cap"] = bool(capt); rec["bottom_cap"] = bool(capb)
            rec["mean_val"] = cal["y_a"] + cal["y_b"] * junction
            rec["top_val"] = cal["y_a"] + cal["y_b"] * u["ya"]
            rec["bottom_val"] = cal["y_a"] + cal["y_b"] * d["yb"]
            rec["mean_source"] = "error-bar junction"
            rec["errbar_pids"] = [u["pid"], d["pid"]]
        else:
            rec["mean_val"] = rec["y_marker_val"]
            rec["top_val"] = rec["bottom_val"] = None
            rec["top_cap"] = rec["bottom_cap"] = False
            rec["mean_source"] = "marker centre (no error bar drawn)"
            rec["errbar_pids"] = []
        out.append(rec)
    out.sort(key=lambda r: (r["kind"], r["x_val"]))
    return out


def main():
    res = {}
    # ---------------- Fig. 1c (main article PDF page 2) -----------------
    P = load("clatot_fig1c_paths.json")
    cal = calibrate(P, axis_y_line=481.44, axis_x_line=319.44,
                    xticks_val=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
                    yticks_val=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # y-axis line top end (axis maximum)
    yaxis = [p for p in P if p["style"].get("stroke") in BLACK and is_seg(p)
             and abs(verts(p)[0][0] - 319.44) < 0.02 and abs(verts(p)[1][0] - 319.44) < 0.02
             and abs(verts(p)[0][1] - verts(p)[1][1]) > 50]
    top = min(min(v[1] for v in verts(p)) for p in yaxis)
    cal["y_axis_top_pt"] = top
    cal["y_axis_top_val"] = cal["y_a"] + cal["y_b"] * top
    # legend key symbols sit at x 404-410 pt, y 397-424 pt (next to the words WT, WT+difopein, S460A)
    pts = extract(P, cal, plot_box=(319.44, top, 449.39, 481.44), legend_box=(403.0, 396.0, 411.0, 425.0))
    res["fig1c"] = {"calibration": cal, "points": pts}
    # ---------------- Supplementary Fig. 3C (SI PDF page 3) -------------
    P = load("clatot_S3C_paths.json")
    cal = calibrate(P, axis_y_line=431.46, axis_x_line=524.45,
                    xticks_val=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
                    yticks_val=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
    yaxis = [p for p in P if p["style"].get("stroke") in BLACK and is_seg(p)
             and abs(verts(p)[0][0] - 524.45) < 0.02 and abs(verts(p)[1][0] - 524.45) < 0.02
             and abs(verts(p)[0][1] - verts(p)[1][1]) > 50]
    top = min(min(v[1] for v in verts(p)) for p in yaxis)
    cal["y_axis_top_pt"] = top
    cal["y_axis_top_val"] = cal["y_a"] + cal["y_b"] * top
    # legend key symbols sit at x 613-620 pt, y 333-366 pt (next to R104W, R104W + S460A, R878C)
    pts = extract(P, cal, plot_box=(524.45, top, 699.83, 431.46), legend_box=(611.0, 331.0, 622.0, 367.0))
    res["figS3C"] = {"calibration": cal, "points": pts}
    json.dump(res, open(OUT, "w"), indent=1)
    for k, v in res.items():
        c = v["calibration"]
        print(f"== {k}: pt per x unit {c['pt_per_x_unit']:.3f}, pt per y unit {c['pt_per_y_unit']:.3f}, "
              f"tick residuals x {c['x_max_tick_residual']:.5f} y {c['y_max_tick_residual']:.5f}, "
              f"y-axis top = {c['y_axis_top_val']:.4f}")
        for r in v["points"]:
            tv = "NA" if r["top_val"] is None else f"{r['top_val']:.4f}"
            bv = "NA" if r["bottom_val"] is None else f"{r['bottom_val']:.4f}"
            print(f"  {r['kind']:8s} x={r['x_val']:.4f} mean={r['mean_val']:.4f} "
                  f"(marker {r['y_marker_val']:.4f}) top={tv} bottom={bv} "
                  f"caps={r['top_cap']},{r['bottom_cap']} src={r['mean_source']}")


if __name__ == "__main__":
    main()

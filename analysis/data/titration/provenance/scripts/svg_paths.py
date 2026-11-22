"""Extract vector paths (page coordinates, y down, PDF points) from a
pdftocairo-generated SVG page.

Usage (always with python3 -I):
    python3 -I svg_paths.py <page.svg> <out.json> [x0 y0 x1 y1]

Only paths outside <defs> are exported (glyph outlines live in <defs> and are
referenced by <use>; text positions are taken from pdftotext -bbox instead).
Each exported path records: style (fill, stroke, stroke-width, dasharray),
the cumulative transform, the list of sub-paths as lists of (x, y) vertices
in page coordinates (Bezier control points included, flagged), and a bbox.
The optional clip box keeps only paths whose bbox intersects it.
"""
import json
import math
import re
import sys
import xml.etree.ElementTree as ET

SVGNS = "{http://www.w3.org/2000/svg}"


def parse_transform(t):
    m = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    if not t:
        return m
    for name, args in re.findall(r"(\w+)\s*\(([^)]*)\)", t):
        v = [float(a) for a in re.split(r"[ ,]+", args.strip()) if a]
        if name == "matrix":
            n = v
        elif name == "translate":
            n = [1, 0, 0, 1, v[0], v[1] if len(v) > 1 else 0.0]
        elif name == "scale":
            n = [v[0], 0, 0, v[1] if len(v) > 1 else v[0], 0, 0]
        elif name == "rotate":
            a = math.radians(v[0])
            n = [math.cos(a), math.sin(a), -math.sin(a), math.cos(a), 0, 0]
        else:
            raise ValueError("unsupported transform " + name)
        m = mul(m, n)
    return m


def mul(a, b):
    # compose: first b then a  (SVG: parent * child)
    return [
        a[0] * b[0] + a[2] * b[1],
        a[1] * b[0] + a[3] * b[1],
        a[0] * b[2] + a[2] * b[3],
        a[1] * b[2] + a[3] * b[3],
        a[0] * b[4] + a[2] * b[5] + a[4],
        a[1] * b[4] + a[3] * b[5] + a[5],
    ]


def apply(m, x, y):
    return (m[0] * x + m[2] * y + m[4], m[1] * x + m[3] * y + m[5])


TOK = re.compile(r"[MmLlHhVvCcZzQqSsTtAa]|[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?")


def parse_d(d):
    toks = TOK.findall(d)
    i = 0
    subpaths = []
    cur = None
    x = y = 0.0
    sx = sy = 0.0
    cmd = None
    while i < len(toks):
        t = toks[i]
        if re.match(r"[A-Za-z]", t):
            cmd = t
            i += 1
            if cmd in "Zz":
                if cur is not None:
                    cur["closed"] = True
                    subpaths.append(cur)
                    cur = None
                x, y = sx, sy
                continue
        def num():
            nonlocal i
            v = float(toks[i]); i += 1
            return v
        if cmd in "Mm":
            nx, ny = num(), num()
            if cmd == "m":
                nx += x; ny += y
            if cur is not None:
                subpaths.append(cur)
            cur = {"pts": [(nx, ny, 0)], "closed": False}
            x, y = nx, ny
            sx, sy = nx, ny
            cmd = "L" if cmd == "M" else "l"
        elif cmd in "Ll":
            nx, ny = num(), num()
            if cmd == "l":
                nx += x; ny += y
            cur["pts"].append((nx, ny, 0)); x, y = nx, ny
        elif cmd in "Hh":
            nx = num()
            if cmd == "h":
                nx += x
            cur["pts"].append((nx, y, 0)); x = nx
        elif cmd in "Vv":
            ny = num()
            if cmd == "v":
                ny += y
            cur["pts"].append((x, ny, 0)); y = ny
        elif cmd in "Cc":
            v = [num() for _ in range(6)]
            if cmd == "c":
                v = [v[0] + x, v[1] + y, v[2] + x, v[3] + y, v[4] + x, v[5] + y]
            cur["pts"].append((v[0], v[1], 1))
            cur["pts"].append((v[2], v[3], 1))
            cur["pts"].append((v[4], v[5], 0))
            x, y = v[4], v[5]
        else:
            raise ValueError("unsupported path command " + cmd)
    if cur is not None:
        subpaths.append(cur)
    return subpaths


def style_of(el, inherited):
    st = dict(inherited)
    for k in ("fill", "stroke", "stroke-width", "stroke-dasharray", "fill-opacity",
              "stroke-opacity", "fill-rule", "stroke-linecap"):
        if el.get(k) is not None:
            st[k] = el.get(k)
    if el.get("style"):
        for kv in el.get("style").split(";"):
            if ":" in kv:
                k, v = kv.split(":", 1)
                st[k.strip()] = v.strip()
    return st


def walk(el, m, st, out, in_defs, clips):
    tag = el.tag.replace(SVGNS, "")
    if tag == "defs":
        in_defs = True
    m2 = mul(m, parse_transform(el.get("transform")))
    st2 = style_of(el, st)
    if tag == "path" and not in_defs:
        sps = parse_d(el.get("d", ""))
        allp = []
        res = []
        for sp in sps:
            pts = [apply(m2, px, py) + (c,) for px, py, c in sp["pts"]]
            res.append({"pts": [[round(a, 4), round(b, 4), c] for a, b, c in pts],
                        "closed": sp["closed"]})
            allp += pts
        if allp:
            xs = [p[0] for p in allp]; ys = [p[1] for p in allp]
            out.append({
                "id": len(out),
                "style": st2,
                "ctm": [round(v, 6) for v in m2],
                "subpaths": res,
                "bbox": [min(xs), min(ys), max(xs), max(ys)],
                "clip": el.get("clip-path") or st2.get("clip-path-inherited"),
            })
    for ch in el:
        st_child = dict(st2)
        if el.get("clip-path"):
            st_child["clip-path-inherited"] = el.get("clip-path")
        walk(ch, m2, st_child, out, in_defs, clips)


def main():
    src, dst = sys.argv[1], sys.argv[2]
    box = [float(v) for v in sys.argv[3:7]] if len(sys.argv) >= 7 else None
    root = ET.parse(src).getroot()
    out = []
    walk(root, [1, 0, 0, 1, 0, 0], {}, out, False, {})
    if box:
        x0, y0, x1, y1 = box
        out = [p for p in out if not (p["bbox"][2] < x0 or p["bbox"][0] > x1 or
                                      p["bbox"][3] < y0 or p["bbox"][1] > y1)]
    json.dump(out, open(dst, "w"), indent=0)
    print(f"{len(out)} paths written to {dst}")


if __name__ == "__main__":
    main()

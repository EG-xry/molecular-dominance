#!/usr/bin/env python3
"""Design schematic: paired assay, assembly model, and requirements.

    python3 run_schematic.py                  # writes Schematic.png and Schematic.pdf

Draws the conceptual figure submitted as Fig. 1: (A) the two arms of a paired
assay, their references, and the panel slope that summarizes them, (B) the
random assembly model that links the arms for a dimer, and (C) the information
each conclusion requires. A dimer is drawn end on, as a disc of two half disc
subunits: wild type blue, mutant gold. The scatter in panel A and the bar
heights in panel B are illustrative values, not data.

Layout: each panel is a tinted card divided into white tiles on a common grid,
with equal gutters, so that every block of content fills its cell. Equations
are sized from their measured width to fill the tile that holds them. Drawn at
its printed width (palette.FIG_WIDTH); layout units are inches before scaling.
Schematic.tiff is written by export_figure_formats.py.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle, Wedge

from palette import (AXIS, BAND, BLUE, BLUE_LIGHT, BLUES, CONTEXT, FIG_WIDTH, GOLD, GRID, INK,
                     MUTED, RULE, TEAL, panel_label, rich_text, tint, use_figure_style)

WT, MUT = BLUE, GOLD             # wild type subunits blue, mutant subunits gold
W, H = 7.0, 4.5                  # layout units; the figure is printed FIG_WIDTH wide
SCALE = FIG_WIDTH / W
LEFT, MID, RIGHT = 0.06, 3.38, W - 0.05      # column edges
TOP, TITLE_Y = 4.24, 4.36        # top edge of the first cards; baseline of panel titles
BOTTOM = 0.08                    # bottom edge of the last cards
PAD, GUTTER = 0.08, 0.08         # card padding and the gap between tiles
R = 0.07                         # radius of a dimer
HEAD = dict(fontsize=8.5, weight="bold", color=INK)
BODY = dict(fontsize=8, color=INK)
SMALL = dict(fontsize=8, color=MUTED)
MATH = dict(fontsize=8.5, color=INK)
KEY_FILL = "#e4e8ed"             # the grey chip behind the colour key
# Panel C row headings: Arial Bold thickened by a thin stroke of the same colour,
# so they stand apart from the grey dataset names on the right. A larger size
# would crowd the last row's heading against its dataset names.
ROW_HEAD = dict(fontsize=8.5, weight="bold", color=INK)
ROW_STROKE = 0.3                 # points


# ------------------------------------------------------------------ drawing helpers
def half(ax, x, y, kind, side, r=R):
    """One subunit of a dimer: the left (side -1) or right (side 1) half of a disc."""
    start = 90 if side < 0 else -90
    ax.add_patch(Wedge((x, y), r, start, start + 180, facecolor=WT if kind == "W" else MUT,
                       edgecolor="white", lw=0.7, zorder=5))


def dimer(ax, x, y, kinds, r=R):
    """A dimer seen end on: two half disc subunits, left and right."""
    half(ax, x, y, kinds[0], -1, r)
    half(ax, x, y, kinds[1], 1, r)


def dimer_row(ax, x_centre, y, pairs, step=0.26, r=R):
    start = x_centre - step * (len(pairs) - 1) / 2
    for k, kinds in enumerate(pairs):
        dimer(ax, start + k * step, y, kinds, r)


def card(ax, x0, y0, x1, y1, fc=BAND, ec="none", lw=0.6):
    """A tinted panel background."""
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle="round,pad=0,rounding_size=0.07",
                                facecolor=fc, edgecolor=ec, lw=lw, zorder=1))


def tile(ax, x0, y0, x1, y1, radius=0.05):
    """A white cell of a card; its content is laid out to fill it."""
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle=f"round,pad=0,rounding_size={radius}",
                                facecolor="white", edgecolor=RULE, lw=0.6, zorder=2))


def rule(ax, x0, x1, y, color=RULE):
    ax.plot([x0, x1], [y, y], color=color, lw=0.6, zorder=3, solid_capstyle="butt")


def arrow(ax, x0, y0, x1, y1, color=MUTED, lw=0.9):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7,
                                 color=color, lw=lw, zorder=4, shrinkA=0, shrinkB=0))


def title(ax, x, y, letter, text):
    panel_label(ax, letter, text, x=x, y=y, pad=0, transform=ax.transData, weight="bold")


def inset(fig, x0, y0, x1, y1):
    """An axes placed in layout units on the figure."""
    new = fig.add_axes([x0 / W, y0 / H, (x1 - x0) / W, (y1 - y0) / H])
    new.patch.set_alpha(0)                  # let the tile show through
    return new


def text_width(ax, text, **props):
    """Width of a string in layout units, measured with the renderer."""
    fig = ax.figure
    screen_dpi = fig.dpi
    fig.set_dpi(600)
    probe = ax.text(0, 0, text, **props)
    width = probe.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi / SCALE
    probe.remove()
    fig.set_dpi(screen_dpi)
    return width


def fitted(ax, texts, width, largest, smallest=8.0, **props):
    """The largest font size, up to `largest`, at which every string fits in `width`.

    Sizes step by a quarter point. PLOS ONE asks for at least 8 pt, so a string
    that would need less raises an error instead of shrinking.
    """
    widest = max(text_width(ax, text, fontsize=largest, **props) for text in texts)
    size = np.floor(largest * min(1.0, width / widest) * 4) / 4
    if size < smallest:
        raise ValueError(f"{texts} need more than {width:.2f} layout units at {smallest} pt")
    return float(size)


# ------------------------------------------------------------------ panel A
ROW_STEP = 0.21                             # distance between two rows of dimers
A_SPLIT = 2.46                              # top of the lower card of panel A


def panel_a(fig, ax):
    title(ax, LEFT, TITLE_Y, "A", "")
    # Colour key on the title line: one mixed dimer with each half named, set on a
    # small grey chip whose right edge lines up with the cards below.
    key_y, gap, pad, half_h = TITLE_Y + 0.035, 0.1, 0.075, 0.075
    key_text = dict(fontsize=8, color=INK)
    right = MID - 0.1
    key_x = right - pad - text_width(ax, "mutant", **key_text) - gap
    left = key_x - gap - text_width(ax, "wild type", **key_text) - pad
    ax.add_patch(FancyBboxPatch((left, key_y - half_h), right - left, 2 * half_h,
                                boxstyle="round,pad=0,rounding_size=0.05",
                                facecolor=KEY_FILL, edgecolor="none", zorder=1))
    dimer(ax, key_x, key_y, "WM", r=0.06)
    ax.text(key_x - gap, key_y, "wild type", ha="right", va="center", **key_text)
    ax.text(key_x + gap, key_y, "mutant", ha="left", va="center", **key_text)

    # Design table: what each arm expresses, what it is compared with, and the
    # effect it yields. Each white dish holds one arm; twice as many complexes
    # mark the doubled total dose of the mixed arm.
    card(ax, LEFT, A_SPLIT + 0.2, MID - 0.1, TOP)
    c_expr, c_ref, c_eff = 1.13, 2.36, 3.07
    for x, label in ((c_expr, "expressed"), (c_ref, "reference"), (c_eff, "effect")):
        ax.text(x, TOP - 0.155, label, ha="center", va="center", **SMALL)
    arms = [
        ("Alone\narm", TOP - 0.425, [["MM"] * 3], [["WW"] * 3],
         "variant only", "wild type, one dose", "$x$"),
        ("Mixed\narm", TOP - 1.04, [["WW", "WM", "MM"], ["WM", "WW", "WM"]], [["WW"] * 3] * 2,
         "variant + wild type", "wild type, two doses", "$y$"),
    ]
    for name, yc, expressed, reference, note_e, note_r, effect in arms:
        rows = len(expressed)
        height = 0.27 + ROW_STEP * (rows - 1)
        ax.text(LEFT + 0.1, yc - 0.05, name, va="center", linespacing=1.1, **HEAD)
        for xc, grid, note in ((c_expr, expressed, note_e), (c_ref, reference, note_r)):
            tile(ax, xc - 0.43, yc - height / 2, xc + 0.43, yc + height / 2, radius=0.13)
            for j, pairs in enumerate(grid):
                dimer_row(ax, xc, yc + ROW_STEP * ((rows - 1) / 2 - j), pairs)
            ax.text(xc, yc - height / 2 - 0.095, note, ha="center", va="center", **BODY)
        ax.text((c_expr + c_ref) / 2, yc, "vs", ha="center", va="center", fontsize=8,
                style="italic", color=MUTED)
        arrow(ax, c_ref + 0.47, yc, c_eff - 0.125, yc)
        ax.add_patch(Circle((c_eff, yc), 0.105, facecolor="white", edgecolor=RULE, lw=0.7, zorder=2))
        ax.text(c_eff, yc + 0.003, effect, ha="center", va="center", fontsize=10, color=INK)
    rule(ax, LEFT + 0.1, MID - 0.2, TOP - 0.75)

    # One (x, y) pair per variant feeds the panel below.
    flow_x = (LEFT + MID - 0.1) / 2
    ax.add_patch(FancyArrowPatch((flow_x, A_SPLIT + 0.175), (flow_x, A_SPLIT + 0.025),
                                 color=CONTEXT, lw=0, zorder=3,
                                 arrowstyle="simple,head_width=7,head_length=5,tail_width=2.6"))

    # Lower card: a plot tile and a formula tile, side by side under one heading.
    card(ax, LEFT, BOTTOM, MID - 0.1, A_SPLIT)
    ax.text(LEFT + 0.1, A_SPLIT - 0.2, "Across a panel of variants", **HEAD)
    y0, y1 = BOTTOM + PAD, A_SPLIT - 0.31
    x0, x1 = LEFT + PAD, MID - 0.1 - PAD
    split = x1 - 1.2                        # the formula tile is 1.2 wide

    # Plot tile: illustrative variants, marker area proportional to the weight x^2.
    tile(ax, x0, y0, split - GUTTER / 2, y1)
    sc = inset(fig, x0 + 0.33, y0 + 0.33, split - GUTTER / 2 - 0.13, y1 - 0.14)
    rng = np.random.default_rng(7)
    xs = np.array([0.12, 0.2, 0.28, 0.35, 0.45, 0.55, 0.62, 0.74, 0.86, 0.95])
    ys = 0.48 * xs + rng.normal(0, 0.05, xs.size) * (0.4 + xs)
    sc.scatter(xs, ys, s=10 + 78 * xs ** 2, color=MUT, edgecolor="white", linewidth=0.5, zorder=3)
    b = float(np.sum(xs * ys) / np.sum(xs ** 2))
    sc.plot([0, 1.02], [0, 1.02 * b], color=TEAL, lw=1.3, zorder=2)
    sc.set_xlim(0, 1.08)
    sc.set_ylim(0, 0.75)
    sc.set_xticks([])
    sc.set_yticks([])
    for side in ("left", "bottom"):
        sc.spines[side].set_visible(False)
    for tip in ((1.0, 0.0), (0.0, 1.0)):    # axes drawn as arrows: directions, not scales
        sc.annotate("", xy=tip, xytext=(0, 0), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=AXIS, lw=0.8, mutation_scale=7,
                                    shrinkA=0, shrinkB=0))
    sc.set_xlabel("alone state effect $x$", fontsize=8, labelpad=3)
    sc.set_ylabel("mixed state effect $y$", fontsize=8, labelpad=3)
    sc.annotate("slope $b$", (1.02, 1.02 * b), xytext=(2, 5), textcoords="offset points",
                ha="right", va="bottom", fontsize=8, color=INK)

    # Formula tile: three cells that share its width and together fill its height.
    fx0 = split + GUTTER / 2
    tile(ax, fx0, y0, x1, y1)
    xc, inner = (fx0 + x1) / 2, x1 - fx0 - 0.2
    slope = [r"$b=\sum_i x_iy_i\,/\sum_i x_i^2$", r"$\;\;=\sum_i w_i\,d_i$"]
    # The weights written out in full: Arial has no proportional-to sign, and
    # PLOS ONE accepts only Arial, Times or Symbol in figures.
    weight = r"$w_i=x_i^2/\sum_j x_j^2$"
    size = fitted(ax, slope + [weight], inner, 10.5, color=INK)
    left = xc - text_width(ax, slope[0], fontsize=size) / 2
    span = y1 - y0
    ax.text(xc, y1 - 0.075 * span, "panel slope", ha="center", va="center", **HEAD)
    ax.text(left, y1 - 0.215 * span, slope[0], va="baseline", fontsize=size, color=INK)
    ax.text(left, y1 - 0.35 * span, slope[1], va="baseline", fontsize=size, color=INK)
    rule(ax, fx0 + 0.1, x1 - 0.1, y1 - 0.44 * span)
    ax.text(xc, y1 - 0.52 * span, "ratio of each variant", ha="center", va="center", **SMALL)
    ax.text(xc, y1 - 0.625 * span, r"$d_i=y_i/x_i$", ha="center", va="center", fontsize=size, color=INK)
    rule(ax, fx0 + 0.1, x1 - 0.1, y1 - 0.72 * span)
    ax.text(xc, y1 - 0.80 * span, "weight (marker area)", ha="center", va="center", **SMALL)
    ax.text(xc, y1 - 0.905 * span, weight, ha="center", va="center", fontsize=size, color=INK)


# ------------------------------------------------------------------ panel B
B_BOTTOM = 2.30
B_ROW = 3.14                                # centre of the gutter between the two tile rows


def complexes_table(ax, x0, x1, y_top):
    """Complex, probability, and activity, as a table that fills its tile."""
    width = x1 - x0
    cols = (x0 + 0.16 * width, x0 + 0.50 * width, x0 + 0.835 * width)
    for x, label in zip(cols, ("complex", "probability", "activity")):
        ax.text(x, y_top - 0.125, label, ha="center", va="center", **SMALL)
    rule(ax, x0 + 0.08, x1 - 0.08, y_top - 0.235)
    rows = [("WW", r"$(1-p)^2$", r"$a_0=1$"), ("WM", r"$2p(1-p)$", r"$a_1$"),
            ("MM", r"$p^2$", r"$a_2=m$")]
    row_h = 0.25
    for k, (kinds, prob, act) in enumerate(rows):
        yy = y_top - 0.375 - k * row_h
        dimer(ax, cols[0], yy, kinds)
        ax.text(cols[1], yy, prob, ha="center", va="center", fontsize=9, color=INK)
        ax.text(cols[2], yy, act, ha="center", va="center", fontsize=9, color=INK)
        if k:                               # hairline between rows
            rule(ax, x0 + 0.08, x1 - 0.08, yy + row_h / 2, color=GRID)


def formula_tile(ax, x0, y0, x1, y1, heading, equation, note, size):
    """A heading, one display equation sized to the tile, and the symbol it needs defined."""
    tile(ax, x0, y0, x1, y1)
    xc = (x0 + x1) / 2
    ax.text(xc, y1 - 0.13, heading, ha="center", va="center", **HEAD)
    ax.text(xc, (y0 + y1) / 2, equation, ha="center", va="center", fontsize=size, color=INK)
    ax.text(xc, y0 + 0.115, note, ha="center", va="center", **SMALL)


def panel_b(fig, ax):
    """The mixed readout as the sum of its three contributions, and what follows from it.

    Illustrative activities a1 = 0.4 and m = 0.2 at p = 0.5. Each bar segment
    takes the colours of its complex, so the mixed complex is half blue and
    half gold. The wild type complexes alone contribute (1-p)^2 of the two dose
    reference: that is the floor, whatever the other two activities are.
    """
    title(ax, MID + 0.06, TITLE_Y, "B", "")
    card(ax, MID + 0.06, B_BOTTOM, RIGHT, TOP, fc=tint(BLUE_LIGHT, 0.24))
    x0, x1 = MID + 0.06 + PAD, RIGHT - PAD
    y_top = TOP - PAD
    top_row = B_ROW + GUTTER / 2            # bottom edge of the upper tiles
    split = x0 + 1.95                       # the table tile is wider than the bar tile

    tile(ax, x0, top_row, split - GUTTER / 2, y_top)
    complexes_table(ax, x0, split - GUTTER / 2, y_top)

    bx0 = split + GUTTER / 2
    tile(ax, bx0, top_row, x1, y_top)
    ax.text((bx0 + x1) / 2, y_top - 0.125, "readout at $p=0.5$", ha="center", va="center", **SMALL)
    bars = inset(fig, bx0 + 0.10, top_row + 0.2, x1 - 0.34, y_top - 0.23)
    p, a1, m = 0.5, 0.4, 0.2
    width = 0.6
    bottom = 0.0
    for kinds, height in (("WW", (1 - p) ** 2), ("WM", 2 * p * (1 - p) * a1), ("MM", p * p * m)):
        if kinds[0] == kinds[1]:
            bars.bar(0, height, bottom=bottom, width=width, color=WT if kinds[0] == "W" else MUT, lw=0)
        else:                               # left and right halves, as in the dimer
            for side, kind in zip((-1, 1), kinds):
                bars.bar(side * width / 4, height, bottom=bottom, width=width / 2,
                         color=WT if kind == "W" else MUT, lw=0)
        bottom += height
        if kinds != "MM":
            bars.plot([-width / 2, width / 2], [bottom, bottom], color="white", lw=0.9)
    bars.bar(1, 1.0, width=width, color=tint(BLUE, 0.45), lw=0)
    floor = (1 - p) ** 2
    bars.plot([-0.4, 1.4], [floor, floor], color=INK, lw=0.8, ls=(0, (4, 2.5)))
    bars.text(1.47, floor, "floor", va="center", ha="left", fontsize=8, color=INK)
    bars.set_xlim(-0.45, 1.45)
    bars.set_ylim(0, 1.04)
    bars.set_xticks([0, 1], ["mixed", "WT + WT"], fontsize=8)
    bars.tick_params(axis="x", length=0, pad=2)
    bars.set_yticks([])
    bars.spines["left"].set_visible(False)

    # Two equal tiles; both display equations take the size at which the wider one fills its tile.
    mid = (x0 + x1) / 2
    transfer = r"$d(p)=2p(1-p)\,g_1+p^2$"
    floor_eq = r"$I_{\mathrm{mixed}}\geq I_{\mathrm{WT+WT}}\,(1-p)^2$"
    size = fitted(ax, [transfer, floor_eq], mid - GUTTER / 2 - x0 - 0.28, 11, color=INK)
    low_row = B_ROW - GUTTER / 2            # top edge of the lower tiles
    formula_tile(ax, x0, B_BOTTOM + PAD, mid - GUTTER / 2, low_row, "transfer ratio", transfer,
                 r"with $g_1=(1-a_1)/(1-m)$", size)
    formula_tile(ax, mid + GUTTER / 2, B_BOTTOM + PAD, x1, low_row, "floor, if activities \u2265 0",
                 floor_eq, r"$I_{\mathrm{WT+WT}}$: wild type at two doses", size)


# ------------------------------------------------------------------ panel C
def panel_c(fig, ax):
    title(ax, MID + 0.06, B_BOTTOM - 0.2, "C", "")
    steps = [
        ("Descriptive panel slope", "both arms on stated scales and origins",
         [("Kv7 (origin assumed), ", False), ("TP53", True)]),
        ("Corrected or cross screen slope", "+ same arm controls; anchors suited to the response",
         [("KCNQ2", True), (", ", False), ("TP53", True)]),
        ("Conditional dimer floor", r"+ achieved fraction $p$ and abundance per arm",
         [("SCN5A", True), (" (", False), ("p", True), (" assumed)", False)]),
        ("Unrestricted profile", r"+ $\geq n-1$ interior fractions; one profile across them",
         [("titrations (assembly fractions assumed)", False)]),
    ]
    x0, y_top, gap, indent = MID + 0.06, B_BOTTOM - 0.33, 0.05, 0.1
    h = (y_top - BOTTOM - gap * (len(steps) - 1)) / len(steps)
    for k, (heading, need, data) in enumerate(steps):
        y1 = y_top - k * (h + gap)
        xs = x0 + k * indent
        # Cumulative requirements: the blue ramp darkens as each row adds information.
        card(ax, xs, y1 - h, RIGHT, y1, fc=tint(BLUES[k], 0.16), ec=tint(BLUES[k], 0.6))
        ax.add_patch(Rectangle((xs + 0.04, y1 - h + 0.06), 0.035, h - 0.12, facecolor=BLUES[k],
                               edgecolor="none", zorder=2))
        ax.text(xs + 0.13, y1 - 0.165, heading,
                path_effects=[patheffects.withStroke(linewidth=ROW_STROKE, foreground=INK)],
                **ROW_HEAD)
        ax.text(xs + 0.13, y1 - 0.34, need, **BODY)
        rich_text(ax, RIGHT - 0.09, y1 - 0.165, data, transform=ax.transData, ha="right",
                  fontsize=8, color=MUTED)
        # The heading and the dataset names share a line; keep a clear gap between them.
        heading_end = xs + 0.13 + text_width(ax, heading, **ROW_HEAD)
        names = "".join(text for text, _ in data)
        if RIGHT - 0.09 - text_width(ax, names, fontsize=8) - heading_end < 0.15:
            raise ValueError(f"Row heading {heading!r} runs into its dataset names")


# ------------------------------------------------------------------ assembly
def draw():
    use_figure_style()
    fig = plt.figure(figsize=(W * SCALE, H * SCALE))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    panel_a(fig, ax)
    panel_b(fig, ax)
    panel_c(fig, ax)
    return fig


def main():
    fig = draw()
    fig.savefig("Schematic.png", dpi=300)
    fig.savefig("Schematic.pdf")
    print("wrote Schematic.png, Schematic.pdf")
    return fig


if __name__ == "__main__":
    fig = main()

#!/usr/bin/env python3
"""Figure 2: conditional origin sensitivity and descriptive TP53 score slopes.

Run run_deconvolution.py, run_tp53.py and run_tp53_hotspot.py before this
driver. The existing kv7_controls.csv contains the all controls residue
bootstrap intervals for panel A. No additive reference is defined on the
separately normalized TP53 score axes. Panel B contains the two mixed state
screens only: slopes on the deposited score axes and, as a conditional
calibration sensitivity, after anchoring every axis at the synonymous median
(0) and the selected hotspot median (1).

The figure is drawn at its printed width (palette.FIG_WIDTH). A soft band
groups the series of each readout or calibration in place of a category axis
line. Legends show marker shape and colour only; interval types and cohort
sizes are stated in the caption.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import blended_transform_factory

from palette import (BLUE_DARK, BLUE_MID, CONTEXT, FIG_WIDTH, GOLD, GOLD_LIGHT, INK, PILL, TEAL,
                     TEAL_LIGHT, group_bands, marker_handle, panel_label, soft_grid,
                     use_figure_style)

RING = dict(mec="white", mew=0.5)           # white ring that lifts a marker off its interval
FIG_HEIGHT = 3.95
BAND_HALF = 0.43                            # half width of each calibration band in panel B


def below(ax, points):
    """Axes fraction for a distance in points below the x axis."""
    height = ax.get_position().height * ax.figure.get_figheight() * 72
    return -points / height


def categorical_x(ax, positions, labels):
    ax.set_xticks(positions, labels)
    ax.tick_params(axis="x", length=0, pad=5)
    ax.spines["bottom"].set_visible(False)


def panel_a(ax, dec, controls):
    axes = ["dV1/2 (mV)", "current at +40 mV", "current at -40 mV"]
    short = ["$\\Delta V_{1/2}$", "current\n+40 mV", "current\n$-$40 mV"]
    xs = np.arange(len(axes))
    dec, controls = dec.loc[axes], controls.loc[axes]
    if not np.allclose(dec.b_control, controls.b, atol=1e-4):
        raise ValueError("control intervals and moment correction point estimates disagree")
    group_bands(ax, xs)
    soft_grid(ax, "y")

    # Raw slopes in the gold pair, candidate adjustments in blue and teal; each
    # series also has its own marker and horizontal offset.
    raw, adjusted = [], []
    ax.errorbar(xs - 0.26, dec.b_control,
                yerr=[dec.b_control - controls.ci_lo, controls.ci_hi - dec.b_control],
                fmt="^", ms=5.6, color=GOLD, mfc=GOLD_LIGHT, mec=GOLD, mew=0.7, lw=1.0)
    raw.append(marker_handle("^", GOLD_LIGHT, 5.6, mec=GOLD, mew=0.7, label="benign controls"))
    ax.plot(xs - 0.13, dec.b_panel_raw, "o", ms=5.4, color=GOLD, **RING)
    raw.append(marker_handle("o", GOLD, 5.4, label="disease panel", **RING))
    top = float(controls.ci_hi.max())
    for offset, key, marker, color, size, label in [
        (0.00, "as_published", "D", BLUE_DARK, 4.8, "published (offset biased)"),
        (0.13, "control_centred", "s", TEAL, 5.0, "control centred"),
        (0.26, "intercept", "v", BLUE_MID, 5.6, "centred with intercept"),
    ]:
        point = dec[f"d_{key}"]
        ax.errorbar(xs + offset, point,
                    yerr=[point - dec[f"ci_lo_{key}"], dec[f"ci_hi_{key}"] - point],
                    fmt=marker, ms=size, color=color, lw=1.0, **RING)
        adjusted.append(marker_handle(marker, color, size, label=label, **RING))
        top = max(top, float(dec[f"ci_hi_{key}"].max()))
    if top >= 1.0:
        raise ValueError("an interval leaves the plotted range; raise the y limit")
    categorical_x(ax, xs, short)
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("slope or moment specification")
    panel_label(ax, "A", [("KCNQ2", True), (": sensitivity to the origin", False)])
    # Two titled groups below the axis keep the plotting area free of text.
    style = dict(loc="upper left", alignment="left", borderaxespad=0, borderpad=0,
                 labelspacing=0.35, handlelength=1.2, handletextpad=0.5,
                 title_fontproperties={"weight": "bold", "size": 8})
    first = ax.legend(handles=raw, title="raw slope", bbox_to_anchor=(0.0, below(ax, 36)), **style)
    ax.add_artist(first)
    ax.legend(handles=adjusted, title="moment specification",
              bbox_to_anchor=(0.43, below(ax, 36)), **style)


def panel_b(ax, tp, hot, labels):
    """x groups are the two calibrations; the two screens sit side by side in each.

    A dotted line joins each screen across the two calibrations, and the
    Boettcher/Giacomelli ratio is printed under the pair it compares.
    """
    screens = list(tp.screen)
    first, second = screens
    marker = {first: ("s", 5.8), second: ("^", 6.6)}
    dep = dict(zip(tp.screen, tp.b))
    dep_lo = dict(zip(tp.screen, tp.ci_lo))
    dep_hi = dict(zip(tp.screen, tp.ci_hi))
    off = {first: -0.14, second: 0.14}
    group_bands(ax, [0, 1], half=BAND_HALF)
    soft_grid(ax, "y")
    for n in screens:
        shape, size = marker[n]
        a = hot["anchored_slopes"][n]
        lo, hi = hot["anchored_bootstrap_ci"][n]
        ax.plot([off[n], 1 + off[n]], [dep[n], a], color=CONTEXT, lw=0.8, ls=(0, (1, 2.2)), zorder=1)
        ax.errorbar(0 + off[n], dep[n], yerr=[[dep[n] - dep_lo[n]], [dep_hi[n] - dep[n]]],
                    fmt=shape, ms=size, color=GOLD, lw=1.1, **RING)
        ax.errorbar(1 + off[n], a, yerr=[[a - lo], [hi - a]], fmt=shape, ms=size,
                    color=TEAL, lw=1.1, **RING)
        l1, l2 = hot["leave_one_out_range"][n]
        ax.plot([1 + off[n] + 0.085] * 2, [l1, l2], color=TEAL_LIGHT, lw=3,
                solid_capstyle="butt")
    rr = hot["raw_ratio_second_over_first"]
    ar = hot["anchored_ratio_second_over_first"]
    al, ah = hot["anchored_ratio_bootstrap_ci"]
    categorical_x(ax, [0, 1], ["deposited\nscore axes", "hotspot\nanchored axes"])
    # The ratio row sits under the pair it compares and carries its own axis
    # style label, so no text is placed inside the plotting area.
    for x, text in ((0, f"{rr:.2f}"), (1, f"{ar:.2f} (95% CI {al:.2f}–{ah:.2f})")):
        ax.annotate(text, (x, 0), xycoords=("data", "axes fraction"), xytext=(0, -36),
                    textcoords="offset points", ha="center", va="top", color=INK, bbox=PILL)
    ax.annotate("slope ratio, Boettcher ÷ Giacomelli", (0.5, 0), xycoords="axes fraction",
                xytext=(0, -56), textcoords="offset points", ha="center", va="top", color=INK)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, 0.85)
    ax.set_ylabel("through origin slope on the common\nalone state axis")
    panel_label(ax, "B", [("TP53", True), (": two screens, two calibrations", False)])
    # The legend sits inside the hotspot anchored band, above its data. PLOS ONE
    # asks for at least 8 pt, so the labels wrap onto two lines instead of shrinking.
    handles = [marker_handle(marker[first][0], INK, 5.4, label=labels[first].replace(", ", ",\n"), **RING),
               marker_handle(marker[second][0], INK, 6.2, label=labels[second].replace(", ", ",\n"), **RING),
               plt.Line2D([], [], color=TEAL_LIGHT, lw=3, label="leave one hotspot\nout range")]
    legend = ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(1.0, 0.985),
                       bbox_transform=blended_transform_factory(ax.transData, ax.transAxes),
                       borderaxespad=0, borderpad=0, labelspacing=0.55, handlelength=1.1,
                       handletextpad=0.5)
    ax.figure.canvas.draw()
    box = legend.get_window_extent()
    (band_left, _), (band_right, _) = ax.transData.transform([(1 - BAND_HALF, 0), (1 + BAND_HALF, 0)])
    margin = 4 * ax.figure.dpi / 72
    if box.x0 < band_left + margin or box.x1 > band_right - margin:
        raise ValueError("the panel B legend no longer fits inside the hotspot anchored band")
    highest = max(max(hot["anchored_bootstrap_ci"][n][1], *hot["leave_one_out_range"][n]) for n in screens)
    if box.y0 < ax.transData.transform([(1, highest)])[0][1] + margin:
        raise ValueError("the panel B legend overlaps the anchored estimates")


def main():
    use_figure_style()
    controls = pd.read_csv("kv7_controls.csv")
    controls = controls[controls.panel == "all controls"].set_index("axis")
    dec = pd.read_csv("deconvolution.csv").set_index("axis")
    tp = pd.read_csv("tp53_portability.csv")
    tp = tp[(tp.scope == "common variant set") & (tp.state == "mixed")]
    if len(tp) != 2 or tp.n.nunique() != 1:
        raise ValueError("Figure 2B requires two mixed screens on the same variant cohort")
    tp = tp.sort_values("b", ascending=False)
    labels = {"Giacomelli, WT background + nutlin-3": "Giacomelli, WT + nutlin-3",
              "Boettcher, K562 WT-activity reporter": "Boettcher, K562 reporter"}
    hot = json.load(open("tp53_hotspot.json"))

    # Fixed axes positions leave room below each panel for its legend or ratio row.
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
    ax1 = fig.add_axes([0.075, 0.31, 0.40, 0.61])
    ax2 = fig.add_axes([0.585, 0.31, 0.40, 0.61])
    panel_a(ax1, dec, controls)
    panel_b(ax2, tp, hot, labels)
    fig.savefig("Figure2.png", dpi=300)
    fig.savefig("Figure2.pdf")
    fig.savefig("Figure2.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("wrote Figure2.png, Figure2.pdf, Figure2.tiff")


if __name__ == "__main__":
    main()

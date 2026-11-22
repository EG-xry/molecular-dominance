#!/usr/bin/env python3
"""Draw the main titration figure from the existing, frozen fit results.

No fitting or simulation occurs here. Source-stated input fractions, all
estimates, uncertainties and model tests are read directly from titration.json.
The complete seven-series display and literal-ratio sensitivity remain in S2.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image

from palette import (BLUE_DARK, CONTEXT, FIG_WIDTH, GOLD, MUTED,
                     RULE, TEAL, halo, panel_label, soft_grid, tint, use_figure_style)
from run_titration import design

HERE = Path(__file__).resolve().parent
CODING = "plotted_positions"
FIG_HEIGHT = 6.4
# Panel C: the two partner series share input fractions, so their points are
# offset horizontally by this much (in p) to keep the intervals apart. The
# fitted curves are drawn at the true fractions. The independent-subunit
# diagonal is a wide pale line under the data, because the S460A fit lies
# almost on it and a thin dotted line would disappear under the dashed fit.
DODGE = .02
REFERENCE = tint(CONTEXT, .42)


def profile_curve(fit, p):
    n = fit["n_subunits"]
    return p ** n + design(p, n) @ np.asarray(fit["profile_estimate"])


def fit_for(results, suffix):
    matches = [v["codings"][CODING] for k, v in results.items() if k.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"Expected one source-stated series for {suffix!r}")
    return matches[0]


def observations(ax, fit, color, marker, shift=0.0, hollow=False):
    # These intervals use the same propagated endpoint covariance as S2.
    face = dict(markerfacecolor="white", markeredgecolor=color, markeredgewidth=.9) if hollow \
        else dict(markeredgecolor="white", markeredgewidth=.4)
    return ax.errorbar(np.asarray(fit["mutant_fractions"]) + shift, fit["d_observed"],
                       yerr=1.96 * np.asarray(fit["d_se"]), fmt=marker,
                       color=color, ms=4.5, capsize=2, lw=.85, elinewidth=.85,
                       zorder=5, **face)


def axes_at(fig, left, bottom, width, height):
    """An axes placed in inches on the figure."""
    w, h = fig.get_size_inches()
    return fig.add_axes([left / w, bottom / h, width / w, height / h])


def main():
    data = json.loads((HERE / "titration.json").read_text())
    if data["primary_coding"] != CODING:
        raise ValueError("Main figure requires the documented input design")
    results = data["series"]
    katp = fit_for(results, "SUR2A wild type + DNN")
    wt = fit_for(results, ": WT + R104W")
    s460a = fit_for(results, ": S460A + R104W")
    p = np.linspace(0, 1, 401)
    use_figure_style()
    # Panel letters only: what each panel shows is stated in the caption.
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
    axa = axes_at(fig, .675, 3.75, 2.66, 2.3)
    axb = axes_at(fig, 4.54, 3.75, 2.66, 2.3)
    axc = axes_at(fig, .675, .55, 2.66, 2.3)
    axd = axes_at(fig, 5.14, .55, 2.06, 2.3)

    # A: a sufficient tetramer design estimates a curve and rejects all four
    # simple activation thresholds. Candidate P values use the full 26-test family.
    thresholds = [(name, candidate) for name, candidate in katp["candidates"].items()
                  if name.startswith("requires >=") and "wild-type subunits" in name]
    if len(thresholds) != 4:
        raise ValueError("Expected the four WT-count threshold candidates")
    for i, (_, candidate) in enumerate(thresholds, 1):
        if max(candidate["p_bonferroni_26"], candidate["p_holm_26"]) >= .05:
            raise ValueError("The main-figure threshold finding changed")
        curve = p ** 4 + design(p, 4) @ np.asarray(candidate["profile"])
        axa.plot(p, curve, color=GOLD, lw=.85, ls=(0, (3, 2)), alpha=.72)
        q = .66
        y = q ** 4 + (design([q], 4) @ np.asarray(candidate["profile"]))[0]
        axa.text(q + .018, y, str(i), color=GOLD, fontsize=8, va="center",
                 path_effects=halo(2), zorder=7)
    axa.plot(p, profile_curve(katp, p), color=BLUE_DARK, lw=1.8)
    observations(axa, katp, BLUE_DARK, "o")
    axa.set(xlim=(-.02, 1.02), ylim=(-.2, 1.1),
            xlabel="Input mutant fraction $p$", ylabel="Normalized activity loss $d(p)$")
    axa.set_xticks([0, .25, .5, .75, 1], ["0", ".25", ".50", ".75", "1"])
    axa.set_yticks([0, .5, 1])
    soft_grid(axa)
    panel_label(axa, "A")
    axa.legend(handles=[Line2D([], [], color=BLUE_DARK, marker="o", ms=4, label="Fitted profile"),
                        Line2D([], [], color=GOLD, ls="--", label="WT thresholds 1–4")],
               loc="upper left", bbox_to_anchor=(0, 1.005), borderaxespad=.1,
               fontsize=8, handlelength=1.6, labelspacing=.3)

    # B: the three recovered coefficients are unconstrained and their pointwise
    # intervals therefore extend outside [0,1]. Do not clip those intervals.
    g = np.asarray(katp["profile_estimate"])
    ci = np.asarray(katp["profile_ci"])
    axb.axhspan(0, 1, color=TEAL, alpha=.065, zorder=0)
    axb.plot([1, 2, 3], [0, .5, 1], color=GOLD, marker="D", mfc="white",
             ms=5, lw=1.0, ls="--", label="Equal-arrangement model")
    axb.errorbar([1, 2, 3], g, yerr=np.vstack([g-ci[:, 0], ci[:, 1]-g]),
                 fmt="o", color=BLUE_DARK, ms=5.5, lw=1.2, capsize=3,
                 label="Estimated profile", zorder=5)
    for k, value in enumerate(g, 1):
        axb.annotate(f"{value:.2f}".replace("-", "−"), (k, value), xytext=(7, 0),
                     textcoords="offset points", va="center", fontsize=8,
                     color=BLUE_DARK, path_effects=halo(2))
    axb.set(xlim=(.65, 3.55), ylim=(-.85, 1.55),
            xlabel="Mutant SUR2A subunits per channel $k$", ylabel="Activity-loss coefficient $g_k$")
    axb.set_xticks([1, 2, 3])
    axb.set_yticks([-.5, 0, .5, 1])
    soft_grid(axb)
    panel_label(axb, "B")
    handles, labels = axb.get_legend_handles_labels()
    axb.legend([handles[1], handles[0]], [labels[1], labels[0]], loc="upper left",
               bbox_to_anchor=(-.01, 1.0), borderaxespad=.2, labelspacing=.35,
               handlelength=1.6, fontsize=8)

    # C: two independently fitted series, not a new between-partner test.
    axc.plot(p, p, color=REFERENCE, lw=5, solid_capstyle="butt", zorder=1)
    axc.plot(p, profile_curve(wt, p), color=BLUE_DARK, lw=1.6, zorder=3)
    axc.plot(p, profile_curve(s460a, p), color=TEAL, lw=1.6, ls=(0, (4, 1.6)), zorder=3)
    observations(axc, wt, BLUE_DARK, "o", shift=-DODGE)
    observations(axc, s460a, TEAL, "s", shift=DODGE, hollow=True)
    axc.set(xlim=(-.02, 1.02), ylim=(-.06, 1.13),
            xlabel="Input mutant fraction $p$", ylabel="Normalized activity loss $d(p)$")
    axc.set_xticks([0, .25, .5, .75, 1], ["0", ".25", ".50", ".75", "1"])
    axc.set_yticks([0, .5, 1])
    soft_grid(axc)
    panel_label(axc, "C")
    axc.legend(handles=[Line2D([], [], color=BLUE_DARK, marker="o", ms=4, lw=1.6,
                              label="WT partner: $g_1=0.88$"),
                        Line2D([], [], color=TEAL, marker="s", ms=4, lw=1.6, ls=(0, (4, 1.6)),
                               mfc="white", mec=TEAL, mew=.9, label="S460A partner: $g_1=0.47$"),
                        Line2D([], [], color=REFERENCE, lw=5, solid_capstyle="butt",
                               label="Independent subunits")],
               loc="upper left", bbox_to_anchor=(0, 1.02), borderaxespad=.1,
               handlelength=1.5, labelspacing=.35, fontsize=8)

    # D: all six source-fraction dimer estimates, grouped by variant.
    rows = [(": WT + L325R", "L325R · WT", BLUE_DARK, "o"),
            (": S460A + L325R", "L325R · S460A", TEAL, "s"),
            (": WT + L325R + difopein", "L325R · difopein", TEAL, "^"),
            (": WT + R104W", "R104W · WT", BLUE_DARK, "o"),
            (": S460A + R104W", "R104W · S460A", TEAL, "s"),
            (": WT + R878C", "R878C · WT", GOLD, "D")]
    axd.axvline(.5, color=MUTED, lw=.95, ls=(0, (2, 2)))
    axd.axvline(1, color=RULE, lw=1.05)
    ys = [5.5, 4.5, 3.5, 2, 1, -.5]
    axd.axhspan(.55, 2.45, color=TEAL, alpha=.065, zorder=0)
    for y, (key, label, color, marker) in zip(ys, rows):
        fit = fit_for(results, key)
        value = fit["profile_estimate"][0]
        low, high = fit["profile_ci"][0]
        axd.errorbar(value, y, xerr=[[value-low], [high-value]], fmt=marker,
                     color=color, ms=5, lw=1.2, capsize=2.5, zorder=4)
    axd.set_yticks(ys, [row[1] for row in rows])
    axd.set(xlim=(-.1, 1.08), ylim=(-1.1, 6.15), xlabel="Heterodimer loss $g_1$")
    axd.set_xticks([0, .5, 1], ["0", "0.5", "1"])
    axd.spines["left"].set_visible(False)
    axd.tick_params(axis="y", length=0, pad=6, labelsize=8)
    panel_label(axd, "D", x=-.29)
    # Input fractions, interval type and the meaning of the reference lines are
    # stated in the caption rather than printed on the figure.
    fig.savefig(HERE / "TitrationMain.png", dpi=300)
    fig.savefig(HERE / "TitrationMain.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(HERE / "TitrationMain.tiff", dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(HERE / "TitrationMain.tiff") as rendered:
        rgb = rendered.convert("RGB")
    rgb.save(HERE / "TitrationMain.tiff", dpi=(600, 600), compression="tiff_lzw")
    rgb.close()
    plt.close(fig)
    print(f"Wrote TitrationMain.png/.pdf/.tiff (7.5 × {FIG_HEIGHT} in; source-stated fits, no refitting).")


if __name__ == "__main__":
    main()

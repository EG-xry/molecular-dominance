"""Figure palette and shared figure style for every manuscript and supplementary figure.

The author supplied the blue ramp and the teal pair, and asked for a golden
brown accent in place of the earlier brown.

Blue ramp, light to dark: #9ecae1, #4292c6, #2171b5, #08306b
Teal pair: #0b7a85 (dark), #3ba0a8 (light)
Gold pair: #b8860b (dark goldenrod), #d4a72c (light gold)

Roles, kept the same in every figure
- blue: wild-type subunits, model quantities such as the floor, and ordered
  families (the ramp runs from the most permissive to the strictest test)
- gold: mutant subunits, raw or deposited quantities, and flagged variants
- teal: adjusted, anchored, or fitted quantities and reference thresholds
- slate grey: context marks that carry no comparison

Checked with the palette validator of the dataviz guidance (OKLab distances,
also under simulated colour-vision deficiency). Gold separates from every blue,
from teal, and from the slate context grey for all viewers, and has more than
3:1 contrast on white. Teal and the darker blues are close, so they never mark
two categories in one panel without a second encoding (marker shape, line
style, or position). The light tones (#9ecae1 and #d4a72c) have low contrast on
white, so they carry fills, bands, bars, and wide lines rather than small
markers.
"""

BLUE_LIGHT = "#9ecae1"
BLUE_MID = "#4292c6"
BLUE = "#2171b5"
BLUE_DARK = "#08306b"
BLUES = (BLUE_LIGHT, BLUE_MID, BLUE, BLUE_DARK)

TEAL = "#0b7a85"
TEAL_LIGHT = "#3ba0a8"
GOLD = "#b8860b"
GOLD_LIGHT = "#d4a72c"
GOLD_DARK = "#7a5a06"       # an interval drawn over a gold marker

INK = "#222222"
MUTED = "#5f6368"
CONTEXT = "#8592a3"         # slate grey for marks that carry no comparison
RULE = "#c4ccd4"            # hairlines and chip borders
GRID = "#e9edf1"            # value gridlines, drawn behind the data
BAND = "#f3f5f8"            # neutral cards and alternating group bands
AXIS = "#3b4148"            # axis lines and ticks

# Figures are drawn at their printed size: 7.5 inches is the PLOS ONE maximum
# width (2250 px at 300 dpi), and PLOS asks for Arial text of 8 to 12 pt, so the
# type sizes below are the printed sizes.
FIG_WIDTH = 7.5
LETTER_SIZE = 11            # bold panel letters
TITLE_SIZE = 9              # panel titles

# A white rounded box for a value printed on or beside a mark.
PILL = dict(boxstyle="round,pad=0.3,rounding_size=0.9", facecolor="white", edgecolor=RULE,
            linewidth=0.6)


def use_figure_style():
    """Arial at printed sizes of 8 to 11 pt, open axes, and quiet axis lines."""
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "text.color": INK,
        "axes.titlesize": TITLE_SIZE,
        "axes.labelsize": 8,
        "axes.labelcolor": INK,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.axisbelow": True,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.color": AXIS,
        "ytick.color": AXIS,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "lines.linewidth": 1.2,
        "errorbar.capsize": 0,
        "hatch.linewidth": 0.6,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.sf": "Arial",
        "mathtext.fallback": "stixsans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def rich_text(ax, x, y, parts, transform=None, ha="left", **kwargs):
    """Text set run by run, so that one line can mix upright, italic, and bold.

    parts is a list of (string, style) pairs. style is False for upright, True
    for italic, or a dict of text properties for that run (for example
    {"weight": "bold", "fontsize": 11}). Mathtext italics leave digits upright,
    so gene symbols like KCNQ2 are set as separate italic runs; each run starts
    where the run before it ends, measured with the renderer in its own style.
    """
    from matplotlib.transforms import offset_copy
    fig = ax.figure
    trans = transform or ax.transAxes

    def properties(style):
        if isinstance(style, dict):
            return {**kwargs, **style}
        return {**kwargs, "fontstyle": "italic" if style else "normal"}

    # Measure at a high resolution so that pixel rounding of glyph widths at the
    # screen dpi does not open visible gaps between runs in the exported file.
    screen_dpi = fig.dpi
    fig.set_dpi(600)
    renderer = fig.canvas.get_renderer()

    def width(text, props):
        # Advance width, measured against a guard character so that leading and
        # trailing spaces count.
        probe = ax.text(0, 0, "|" + text + "|", transform=trans, **props)
        guard = ax.text(0, 0, "||", transform=trans, **props)
        w = probe.get_window_extent(renderer).width - guard.get_window_extent(renderer).width
        probe.remove()
        guard.remove()
        return w * 72.0 / fig.dpi

    runs = [(text, properties(style)) for text, style in parts]
    widths = [width(text, props) for text, props in runs]
    fig.set_dpi(screen_dpi)
    total = sum(widths)
    shift = {"left": 0.0, "center": total / 2, "right": total}[ha]
    # Offsets are kept in points, so later layout changes cannot stretch them.
    artists, offset = [], 0.0
    for (text, props), w in zip(runs, widths):
        run = offset_copy(trans, fig=fig, x=offset - shift, y=0, units="points")
        artists.append(ax.text(x, y, text, transform=run, ha="left", **props))
        offset += w
    return artists


def panel_label(ax, letter, title="", x=0.0, y=1.0, pad=6.0, transform=None, weight="normal",
                size=TITLE_SIZE, letter_size=LETTER_SIZE):
    """A bold panel letter followed by the panel title, on one baseline.

    title is a string or a list of rich_text parts (for italic gene symbols).
    The baseline sits pad points above (x, y), which default to the top left
    corner of the axes.
    """
    from matplotlib.transforms import offset_copy
    trans = offset_copy(transform or ax.transAxes, fig=ax.figure, x=0, y=pad, units="points")
    parts = [(title, False)] if isinstance(title, str) else list(title)
    runs = [(letter, {"weight": "bold", "fontsize": letter_size})]
    if parts and parts[0][0]:
        runs += [("  ", False)] + parts
    return rich_text(ax, x, y, runs, transform=trans, fontsize=size, weight=weight,
                     color=INK, va="baseline")


def soft_grid(ax, axis="y"):
    """Light value gridlines behind the data."""
    ax.grid(True, axis=axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def group_bands(ax, centres, half=0.43, horizontal=False, color=BAND):
    """A soft band behind each categorical group, in place of a category axis line."""
    span = ax.axhspan if horizontal else ax.axvspan
    for centre in centres:
        span(centre - half, centre + half, color=color, linewidth=0, zorder=0)


def rounded_bar(ax, base, value, centre, thickness, radius=2.5, horizontal=True, both=False,
                **kwargs):
    """A bar from base to value whose value end has rounded corners.

    With both=True the base end is rounded too, for an interval that floats
    free of an axis. The ends stay exactly at base and value. radius is in
    points; it is converted to data units with the current axes transform, so
    the axes position and limits must be final before the call.
    """
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path
    if value == base:
        return None
    (x0, y0), (x1, y1) = ax.transData.transform([(0, 0), (1, 1)])
    per_point = ax.figure.dpi / 72.0
    along, across = (abs(x1 - x0), abs(y1 - y0)) if horizontal else (abs(y1 - y0), abs(x1 - x0))
    ra = min(radius * per_point / along, abs(value - base) / (2 if both else 1))
    rc = min(radius * per_point / across, thickness / 2)
    sign = 1 if value > base else -1
    lo, hi = centre - thickness / 2, centre + thickness / 2
    k = 1 - 0.5523                       # Bezier offset for a quarter circle
    end = value
    points = [(end - sign * ra, lo),
              (end - sign * ra * k, lo), (end, lo + rc * k), (end, lo + rc),
              (end, hi - rc),
              (end, hi - rc * k), (end - sign * ra * k, hi), (end - sign * ra, hi)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4]
    if both:
        points += [(base + sign * ra, hi),
                   (base + sign * ra * k, hi), (base, hi - rc * k), (base, hi - rc),
                   (base, lo + rc),
                   (base, lo + rc * k), (base + sign * ra * k, lo), (base + sign * ra, lo)]
        codes += [Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
                  Path.CURVE4, Path.CURVE4, Path.CURVE4]
    else:
        points += [(base, hi), (base, lo)]
        codes += [Path.LINETO, Path.LINETO]
    points.append(points[0])
    codes.append(Path.CLOSEPOLY)
    if not horizontal:
        points = [(b, a) for a, b in points]
    kwargs.setdefault("linewidth", 0)
    patch = PathPatch(Path(points, codes), **kwargs)
    ax.add_patch(patch)
    return patch


def halo(width=2.4):
    """A white outline that lifts a line or text off the marks behind it."""
    from matplotlib import patheffects
    return [patheffects.withStroke(linewidth=width, foreground="white")]


def marker_handle(marker, color, size=6, **kwargs):
    """A legend entry that shows only the marker shape and colour, no error bar."""
    from matplotlib.lines import Line2D
    return Line2D([], [], marker=marker, color=color, linestyle="none", markersize=size, **kwargs)


def tint(hex_colour, amount):
    """Mix a colour with white; amount 0 gives white and 1 the colour itself."""
    hex_colour = hex_colour.lstrip("#")
    rgb = [int(hex_colour[i:i + 2], 16) for i in (0, 2, 4)]
    mixed = [round(255 - (255 - c) * amount) for c in rgb]
    return "#" + "".join(f"{c:02x}" for c in mixed)

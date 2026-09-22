# VENDORED COPY -- do not edit the original.
# source: /home/hlamm/Desktop/QC/dg_qudits/note/figures/fnalfig.py
# copied: 2026-09-22  (for fh_iqus; original is read-only)

"""Shared figure style for the FNAL 2026 QC applications report.

House style follows ~/Desktop/QC/lhc_combo_breaker/scripts/: serif/cm at 8pt,
Okabe-Ito CVD-safe hues, ticks inward on all four sides, no on-plot titles (the
LaTeX caption does that work), no boxed legends (labels sit tangent along
curves or beside marks).

One adaptation: that repo targets a 3.4in PRL column. This report's text block
is 6.5in (letter, 1in margins, set in fnal-report.cls), so FULL/HALF below are
sized for it while keeping the 8pt/0.6pt weights, so on-page text matches.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Okabe-Ito, CVD-safe. Same hex values as the reference scripts.
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
YELLOW = "#E69F00"

# Greys, darkest to lightest, as used for rules/bands/annotation in the reference.
INK = "#333333"
DARKGREY = "#555555"
MIDGREY = "#666666"
GREY = "#999999"
LIGHTGREY = "#bbbbbb"
PALEGREY = "#cccccc"

# Fermilab identity, from fnal-report.cls. Reserved for envelope/branding
# elements so it never collides with the Okabe-Ito data hues.
FNAL_BLUE = "#002147"

FULL = (6.5, 3.6)      # full text width
HALF = (3.2, 2.6)      # half width, two side by side
WIDE = (6.5, 2.6)      # full width, short
TALL = (6.5, 4.6)      # full width, tall

FIGDIR = Path(__file__).resolve().parent.parent / "figures"


def setup():
    """Apply the house rcParams. Call once at the top of each figure script."""
    plt.rcParams.update(
        {
            "font.size": 8,
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "axes.linewidth": 0.6,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.minor.width": 0.5,
            "ytick.minor.width": 0.5,
            "legend.frameon": False,
            "axes.unicode_minus": False,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )


def save(fig, name):
    """Write <name>.pdf into doc/Figures/ and report the path."""
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"{name}.pdf"
    fig.savefig(out)
    print(f"wrote {out}")
    return out


def label_along(ax, x, y_of, text, color, side=1, pad=7.0, fontsize=7):
    """Place text tangent to a curve, offset perpendicular to it.

    Lifted from breakeven.py in the reference repo. Must be called after
    fig.tight_layout() and fig.canvas.draw(), since the tangent angle is
    computed in display coordinates and needs the final axes box.
    """
    x0, x1 = x * 0.9, x * 1.1
    p0 = ax.transData.transform((x0, y_of(x0)))
    p1 = ax.transData.transform((x1, y_of(x1)))
    ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
    nx, ny = -np.sin(np.radians(ang)), np.cos(np.radians(ang))
    ax.annotate(
        text,
        xy=(x, y_of(x)),
        xytext=(side * pad * nx, side * pad * ny),
        textcoords="offset points",
        rotation=ang,
        rotation_mode="anchor",
        ha="center",
        va="center",
        color=color,
        fontsize=fontsize,
    )


def despine(ax, keep=("left", "bottom")):
    """Drop spines not in `keep`. For cartoons, pass keep=() for a bare canvas."""
    for side in ("top", "right", "bottom", "left"):
        if side not in keep:
            ax.spines[side].set_visible(False)


def blank(ax):
    """Strip an axes to a bare drawing canvas, for the cartoon figures."""
    despine(ax, keep=())
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")

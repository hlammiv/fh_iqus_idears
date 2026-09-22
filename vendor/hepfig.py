# VENDORED COPY -- do not edit the original.
# source: /home/hlamm/Desktop/QC/hepware_qec/figures_src/hepfig.py
# copied: 2026-09-22  (for fh_iqus; original is read-only)

"""Shared style for the HEP-aware QEC deck figures.

A copy of fnalqc/PAC_meeting/pacfig.py with the output directory moved to
hepware_qec/figures and a manifest entry written on every save, so that
figures/manifest.json records which script produced which PNG.
"""

import json
import os
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

# Deck palette, from IBM_meeting/deck_template.py
BLUE = "#004C97"
ORANGE = "#F88F2E"
INK = "#1A1A1A"
GREY = "#5A5A5A"
RULE = "#D8DDE4"
BAND = "#F2F5F9"

# Accents for the figures only (pass the data-viz contrast checks; the brand
# colors are text colors)
S1, S2, S3 = "#1A4C9A", "#E47818", "#1B8F6E"
RED = "#A33A0E"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "figures")
MANIFEST = os.path.join(OUT, "manifest.json")

sys.path.insert(0, ROOT)          # so figure scripts can `import resources`


def setup():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Carlito", "Calibri", "DejaVu Sans"],
        "axes.edgecolor": RULE, "axes.linewidth": 0.8,
        "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": GREY, "ytick.color": GREY,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.labelsize": 10, "legend.fontsize": 9,
        "legend.frameon": False,
        "figure.dpi": 200, "savefig.dpi": 200,
        "savefig.bbox": "tight", "savefig.facecolor": "white",
        "axes.unicode_minus": False,
        "mathtext.fontset": "dejavusans",
    })


def manifest_update(key, entry):
    os.makedirs(OUT, exist_ok=True)
    data = {}
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as fh:
            data = json.load(fh)
    data[key] = entry
    with open(MANIFEST, "w") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)


def save(fig, name, script=None):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"{name}.png")
    fig.savefig(path)
    plt.close(fig)
    manifest_update(name, {"kind": "mpl",
                           "script": os.path.relpath(script or sys.argv[0], ROOT),
                           "path": os.path.relpath(path, ROOT)})
    print(f"wrote {path}")
    return path


def clean_axes(ax, keep=("left", "bottom")):
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(sp in keep)

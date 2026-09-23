#!/usr/bin/env python3
"""Slide version of the resource estimate -- one panel, one message.

Same models as make_figure.py (fhcost/), but in the slide house style
(vendor/hepfig.py: sans-serif Carlito, dpi 200, PNG + manifest.json). The paper
figure carries three panels; a slide gets the main panel only, with the
conclusion written on it.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt

from vendor import hepfig as H
from fhcost.budget import DEFAULT
from fhcost import curves, classical

H.setup()
C = curves.COLORS
BAND_LO, BAND_HI = classical.band(DEFAULT)["band"]
n_grid = np.logspace(0, 8, 130)
Y = curves.evaluate(n_grid, DEFAULT)
msk = lambda y: np.where(y > 0, y, np.nan)

fig, ax = plt.subplots(figsize=(9.0, 5.0))
ax.axhspan(BAND_LO, BAND_HI, color=H.BAND, zorder=0)
ax.axhline(BAND_HI, color=H.GREY, lw=1.0, ls=(0, (5, 3)), zorder=1)

ax.plot(n_grid, msk(Y["ideal"]), color=C["ideal"], lw=1.6, ls=(0, (6, 3)), zorder=3)
ax.plot(n_grid, msk(Y["nisq_pec"]), color=C["nisq_pec"], lw=2.6, zorder=4)
ax.plot(n_grid, msk(Y["nisq_pec_mpf"]), color=C["nisq_pec"], lw=2.0, ls=(0, (4, 2)), zorder=4)
ax.plot(n_grid, msk(Y["star"]), color=C["star"], lw=2.6, zorder=4)
ax.fill_between(n_grid, msk(Y["surface_willow"]), msk(Y["surface"]),
                color=C["surface"], alpha=0.18, lw=0, zorder=2)
ax.plot(n_grid, msk(Y["surface"]), color=C["surface"], lw=2.8, zorder=5)
ax.plot(n_grid, msk(Y["surface_willow"]), color=C["surface"], lw=1.6, ls=(0, (2, 2)), zorder=5)

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 1e8); ax.set_ylim(2, 1.2e3)
ax.set_xlabel("n  (physical qubits)")
ax.set_ylabel("m  (resolvable lattice sites)")
ax.tick_params(which="both", direction="out", length=3)
H.clean_axes(ax)

def at(key, x):
    lg = np.log10(np.where(Y[key] > 0, Y[key], np.nan))
    return 10 ** np.interp(np.log10(x), np.log10(n_grid), lg)

lab = [("ideal", 3.0e2, "p = 0  (qubit-limited)", C["ideal"], 1.45),
       ("nisq_pec", 1.2e7, "NISQ + PEC", C["nisq_pec"], 0.62),
       ("nisq_pec_mpf", 1.2e7, "NISQ + PEC + MPF-4", C["nisq_pec"], 1.22),
       ("star", 1.2e7, "STAR", C["star"], 0.45),
       ("surface", 2.2e7, "surface-code FT", C["surface"], 1.42)]
for key, x, txt, col, dy in lab:
    ax.text(x, at(key, x) * dy, txt, color=col, fontsize=11, fontweight="bold",
            ha="center", va="center")

ax.text(1.4, BAND_HI * 1.25, f"estimated ED capacity   ({BAND_LO:.0f} – {BAND_HI:.0f} sites)",
        color=H.GREY, fontsize=10, va="bottom", ha="left")
ax.text(1.4, 2.6, "unmitigated p = 10$^{-3}$ never reaches m = 4",
        color=H.GREY, fontsize=9.5, va="bottom", ha="left")

xc = curves.summary()["n_ft_clears_classical_hi"]
ax.plot([xc], [BAND_HI], marker="o", ms=11, mfc="none", mec=C["surface"], mew=2.2, zorder=6)
ax.annotate("only fault tolerance clears\nthe estimated ED capacity",
            xy=(xc, BAND_HI), xytext=(3.0, 260), fontsize=11, color=C["surface"],
            fontweight="bold", ha="left",
            arrowprops=dict(arrowstyle="->", color=C["surface"], lw=1.8,
                            connectionstyle="arc3,rad=-0.25"))

from fhcost.record import model_id
fig.text(0.995, 0.004, f"model {model_id()}", ha="right", va="bottom",
         fontsize=5.5, color="0.6")

fig.tight_layout(pad=0.6)
H.save(fig, "fh_resource_estimate_slide", script=__file__)

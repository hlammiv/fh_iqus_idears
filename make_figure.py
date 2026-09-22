#!/usr/bin/env python3
"""Quantitative replacement for the slide-2 whiteboard sketch of Q1_Fermi_Hubbard-1.pdf.

VALUE PROVENANCE
  Every number is computed by fhcost/ (see fhcost/budget.py for the input ledger
  and METHODS.md for derivations and citations). Nothing here is hand-placed.

PANELS
  (a) resolvable lattice size m vs physical qubits n, one week, p = 1e-3
  (b) what each extrapolation buys at fixed n = 1e6  -- the user's question
  (c) why the time window decides everything: classical reach vs t_max

HOUSE STYLE: vendor/fnalfig.py -- serif/cm 8pt, Okabe-Ito, ticks in on all four
sides, no on-plot titles, no boxed legends (curves are labelled tangentially, so
identity never rests on colour alone). Palette CVD-validated by check_palette.py.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from vendor import fnalfig as F
from fhcost.budget import DEFAULT
from fhcost import curves, classical, nisq, ftqc

F.setup()
C, L = curves.COLORS, curves.LABELS
BAND_LO, BAND_HI = classical.band(DEFAULT)["band"]

n_grid = np.logspace(0, 8, 130)   # from a single physical qubit
Y = curves.evaluate(n_grid, DEFAULT)

fig = plt.figure(figsize=(7.0, 3.05))
gs = fig.add_gridspec(1, 3, width_ratios=[1.40, 0.98, 0.82], wspace=0.55)
axA, axB, axC = (fig.add_subplot(g) for g in gs)

# ---------------------------------------------------------------- (a) m vs n
ax = axA
ax.axhspan(BAND_LO, BAND_HI, color=F.PALEGREY, alpha=0.55, lw=0, zorder=0)
ax.axhline(BAND_HI, color=F.GREY, lw=0.6, ls=(0, (4, 2)), zorder=1)
ax.text(1.4, BAND_HI * 1.18, f"classically easy  ($m\\leq{BAND_HI:.0f}$)",
        fontsize=6.0, color=F.DARKGREY, va="bottom", ha="left")

def msk(y):
    return np.where(y > 0, y, np.nan)

ax.plot(n_grid, msk(Y["ideal"]), color=C["ideal"], lw=1.0, ls=(0, (5, 2)), zorder=3)
if np.nanmax(Y["nisq_none"]) > 0:
    ax.plot(n_grid, msk(Y["nisq_none"]), color=C["nisq_none"], lw=1.4, zorder=3)
else:
    # bias-limited below the smallest real lattice: there is no curve to draw
    ax.text(1.4, 1.35, r"unmitigated: $m<4$ at every $n$ (bias-limited)",
            fontsize=5.8, color=C["nisq_none"], va="bottom", ha="left")
ax.fill_between(n_grid, msk(Y["nisq_pec"]), msk(Y["nisq_pec_mpf"]),
                color=C["nisq_pec"], alpha=0.10, lw=0, zorder=2)
ax.plot(n_grid, msk(Y["nisq_pec"]), color=C["nisq_pec"], lw=1.6, zorder=4)
ax.plot(n_grid, msk(Y["nisq_pec_mpf"]), color=C["nisq_pec"], lw=1.3,
        ls=(0, (3, 1.6)), zorder=4)
ax.plot(n_grid, msk(Y["star"]), color=C["star"], lw=1.6, zorder=4)
ax.fill_between(n_grid, msk(Y["surface_willow"]), msk(Y["surface"]),
                color=C["surface"], alpha=0.16, lw=0, zorder=2)
ax.plot(n_grid, msk(Y["surface"]), color=C["surface"], lw=1.7, zorder=5)
ax.plot(n_grid, msk(Y["surface_willow"]), color=C["surface"], lw=1.0,
        ls=(0, (1.6, 1.4)), zorder=5)
ax.plot(n_grid, msk(Y["pinnacle"]), color=C["pinnacle"], lw=1.9, zorder=6)

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 1e8); ax.set_ylim(1.2, 1.2e3)
ax.set_xticks([10.0 ** e for e in range(0, 9, 2)])
ax.set_xlabel(r"$n$ (physical qubits)")
ax.set_ylabel(r"$m$ (resolvable lattice sites)")
ax.tick_params(which="both", direction="in", top=True, right=True, length=2.6)

def safe_label(ax, x, key, text, color, **kw):
    """label_along, but refuse to place a label outside the axes.

    An off-axis annotation silently drags tight_layout into expanding the canvas
    (this figure once came out 440 inches wide because the p=0 label sat at
    m ~ 6700 with ylim topping out at 1200).
    """
    yf = y_of(key)
    y = yf(x)
    y0, y1 = ax.get_ylim()
    x0, x1 = ax.get_xlim()
    if not np.isfinite(y) or not (y0 < y < y1) or not (x0 < x < x1):
        raise ValueError(f"label {text!r} for {key!r} at n={x:.2g} falls at m={y:.4g}, "
                         f"outside the axes y=({y0:g},{y1:g})")
    F.label_along(ax, x, yf, text, color, **kw)


def y_of(key):
    lg = np.log10(np.where(Y[key] > 0, Y[key], np.nan))
    return lambda x: 10 ** np.interp(np.log10(x), np.log10(n_grid), lg)

fig.tight_layout(pad=0.35)
fig.canvas.draw()
safe_label(ax, 1.3e3, "ideal", r"$p=0$", C["ideal"], side=-1, pad=6.5, fontsize=6.2)
safe_label(ax, 1.1e3, "nisq_pec", "NISQ + PEC", C["nisq_pec"], side=-1, pad=5.5, fontsize=6.2)
safe_label(ax, 2.2e2, "nisq_pec_mpf", "NISQ + PEC + MPF-4", C["nisq_pec"], side=1, pad=5.0, fontsize=6.2)

safe_label(ax, 4e7, "star", "STAR", C["star"], side=-1, pad=5.5, fontsize=6.2)
safe_label(ax, 2.2e7, "surface", "surface-code FT", C["surface"], side=-1, pad=6, fontsize=6.2)
safe_label(ax, 3.2e7, "pinnacle", "Pinnacle (QLDPC)", C["pinnacle"], side=1, pad=6, fontsize=6.2)
safe_label(ax, 1.3e7, "surface_willow", "measured $p_L$", C["surface"], side=-1, pad=5.5, fontsize=5.6)
ax.text(0.965, 0.955, "(a)", transform=ax.transAxes, fontsize=7.5, va="top", ha="right")

xc = curves.summary()["n_pinnacle_clears_classical_hi"]
if xc:
    ax.plot([xc], [BAND_HI], marker="o", ms=5.5, mfc="none",
            mec=C["pinnacle"], mew=1.2, zorder=7)
    ax.annotate("QLDPC clears it\n8$\\times$ earlier", xy=(xc, BAND_HI),
                xytext=(1.5, 520), fontsize=6.2, color=C["pinnacle"], ha="left",
                arrowprops=dict(arrowstyle="->", color=C["pinnacle"], lw=0.8,
                                connectionstyle="arc3,rad=-0.25"))

# ------------------------------------------------- (b) extrapolation ladder
ax = axB
rows = [r for r in curves.extrapolation_table(1e6) if not r[0].startswith("Surface")]
ft_rows = [r for r in curves.extrapolation_table(1e6) if r[0].startswith("Surface")]
rows = rows + ft_rows + [r for r in curves.extrapolation_table(1e6)
                         if r[0].startswith("  + MPF")]
seen, ordered = set(), []
for r in curves.extrapolation_table(1e6):
    if r[0] not in seen:
        seen.add(r[0]); ordered.append(r)
SHORT = {"NISQ, unmitigated": "unmitigated", "  + ZNE order 1": "+ ZNE-1",
         "  + ZNE order 2": "+ ZNE-2", "  + PEC (optimal mitigation)": "+ PEC",
         "  + PEC + MPF order 4": "+ PEC + MPF-4",
         "  + PEC + MPF order 6": "+ PEC + MPF-6",
         "  + PEC + MPF order 8": "+ PEC + MPF-8",
         "Surface-code FT": "surface FT", "  + MPF order 4": "+ MPF-4",
         "  + MPF order 6": "+ MPF-6", "  + MPF order 8": "+ MPF-8"}
names = [SHORT.get(r[0], r[0]) for r in ordered]
vals = [r[1] for r in ordered]
is_ft = [i for i, r in enumerate(ordered)
         if r[0].startswith("Surface") or (i > 0 and ordered[0][0] and False)]
ft_start = next(i for i, r in enumerate(ordered) if r[0].startswith("Surface"))
cols = [C["nisq_pec"]] * ft_start + [C["surface"]] * (len(ordered) - ft_start)
cols[0] = C["nisq_none"]

ypos = np.arange(len(ordered))[::-1]
ax.axvspan(BAND_LO, BAND_HI, color=F.PALEGREY, alpha=0.55, lw=0, zorder=0)
# a zero-length bar would put its value text at log(0) = -inf and drag
# tight_layout into a 440-inch canvas, so clamp to the axis floor
XFLOOR = 1.0
ax.barh(ypos, [max(v, XFLOOR) for v in vals], height=0.62, color=cols, zorder=3)
for y, v in zip(ypos, vals):
    ax.text(max(v, XFLOOR) * 1.12, y, f"{v:.0f}" if v > 0 else r"$<4$",
            va="center", ha="left", fontsize=5.8,
            color=F.INK if v > 0 else F.GREY)
ax.set_yticks(ypos); ax.set_yticklabels(names, fontsize=6.0)
ax.set_xscale("log"); ax.set_xlim(1, 600)
ax.set_xlabel(r"$m$ at $n=10^6$")
ax.tick_params(which="both", direction="in", top=True, length=2.6)
ax.tick_params(axis="y", length=0)
ax.text(0.04, 0.985, "(b)", transform=ax.transAxes, fontsize=7.5, va="top")
ax.text(np.sqrt(BAND_LO * BAND_HI), len(ordered) - 1.25, "classical",
        fontsize=5.8, color=F.DARKGREY, ha="center", va="center", rotation=90)

# ------------------------------------------- (c) the time window decides it
ax = axC
ts = np.logspace(-1.7, 1.0, 90)
mc = np.array([classical.max_m_fixed_t(t) for t in ts])
unb = mc >= classical.UNBOUNDED
ax.plot(ts[~unb], mc[~unb], color=F.DARKGREY, lw=1.6, zorder=3)
if unb.any():
    ax.fill_between(ts[unb], 1, 1e4, color=F.PALEGREY, alpha=0.75, lw=0, zorder=1)
    ax.text(ts[unb][len(ts[unb]) // 2], 260, "classical\nreach\nunbounded",
            fontsize=5.8, color=F.DARKGREY, ha="center", va="center")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(ts[0], ts[-1]); ax.set_ylim(1, 1e4)
ax.set_xlabel(r"$t_{\max}$  $[\hbar/J]$")
ax.set_ylabel(r"largest classically easy $m$")
ax.tick_params(which="both", direction="in", top=True, right=True, length=2.6)
ax.text(0.05, 0.955, "(c)", transform=ax.transAxes, fontsize=7.5, va="top")

for tv, lab, col in ((1.0 / 26.0, r"slide 1: $t \leq 1/m$", F.ORANGE),
                     (np.sqrt(26.0) / DEFAULT.v, r"adopted: $\sqrt{m}/v_B$", F.BLUE)):
    if tv >= ts[0]:
        ax.axvline(tv, color=col, lw=0.9, ls=(0, (3, 2)), zorder=4)
        ax.text(tv * 1.12, 1.6, lab, rotation=90, fontsize=5.8, color=col,
                va="bottom", ha="left")
    else:
        ax.annotate(lab, xy=(ts[0], 3.0), xytext=(ts[0] * 1.5, 3.0), fontsize=5.8,
                    color=col, va="center", ha="left",
                    arrowprops=dict(arrowstyle="->", color=col, lw=0.8))

out = pathlib.Path(__file__).parent / "figures"
out.mkdir(exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(out / f"fh_resource_estimate.{ext}", dpi=300)
print(f"wrote {out}/fh_resource_estimate.pdf and .png")

# ------------------------------------------------ console spot-check (house convention)
print(f"\nclassical band: m = {BAND_LO:.0f} .. {BAND_HI:.0f}")
print(f"{'n':>9} " + " ".join(f"{k:>14}" for k in
      ("ideal", "nisq_none", "nisq_pec", "nisq_pec_mpf", "star", "surface", "surface_willow")))
for n in (1e3, 1e4, 1e5, 1e6, 1e7, 1e8):
    row = curves.evaluate([n], DEFAULT)
    print(f"{n:>9.0e} " + " ".join(
        f"{row[k][0]:>14.4g}" for k in
        ("ideal", "nisq_none", "nisq_pec", "nisq_pec_mpf", "star", "surface", "surface_willow")))

#!/usr/bin/env python3
"""Platform comparison: connectivity bought against clock speed.

    python3 make_platforms.py   ->  figures/fh_platforms.pdf / .png

THE QUESTION
  Every other figure here assumes one machine: a nearest-neighbour superconducting
  grid, 10 ns gates, unlimited gate parallelism. Trapped-ion and neutral-atom
  machines offer the connectivity that qLDPC codes need and that would delete the
  fermionic swap network entirely. Is that trade worth taking?

WHY IT IS NOT OBVIOUS A PRIORI
  All-to-all removes real cost: at m = 256 the Jordan-Wigner swap network is 44
  layers per Trotter step against 12 without it, and the routing penalty goes to
  zero. Helios also has a BETTER two-qubit gate than the slide assumes,
  7.9e-4 against 1e-3. Against that, its layer takes 55 ms.

  Both sides are computed. Neither is asserted.

PANEL (b) IS THE POINT
  Which constraint binds -- qubits, or the one-week clock -- decides whether gate
  speed matters at all. On the superconducting grid the surface code is
  qubit-limited until n ~ 1e7 and gate speed is nearly free; past that it is
  clock-limited and a slower machine loses linearly.
"""
import sys, pathlib, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
import matplotlib.pyplot as plt
from vendor import fnalfig as F
from fhcost.budget import DEFAULT
from fhcost import curves, classical, nisq, ftqc, hubbard
from fhcost import platform as plat
from fhcost.record import model_id

F.setup()
BAND_LO, BAND_HI = classical.band(DEFAULT)["band"]
n_grid = np.logspace(3, 12, 110)   # out to 1e12: the qLDPC arms need it

# (label, config, colour, style). Each arm on the hardware it actually needs.
# (label, cfg, arm, colour, linestyle, label position along the curve)
ARMS = [
    ("SC grid, surface code", DEFAULT.but(platform="superconducting",
                                          use_platform_clock=True),
     "surface", curves.COLORS["surface"], "-", 0.72),
    ("SC x2 layers, QLDPC", DEFAULT.but(platform="sc_long_range",
                                        use_platform_clock=True),
     "pinnacle", curves.COLORS["pinnacle"], "-", 0.55),
    ("Helios, NISQ+PEC", DEFAULT.but(platform="helios", use_platform_clock=True,
                                     encoding="jw", p=7.9e-4),
     "nisq", curves.COLORS["nisq_none"], (0, (4, 2)), 0.45),
    ("atoms, NISQ+PEC ($p{=}5{\\times}10^{-3}$)", DEFAULT.but(platform="neutral_atom",
                                    use_platform_clock=True, encoding="jw", p=5e-3),
     "nisq", curves.COLORS["star"], (0, (1.5, 1.5)), 0.80),
    ("SC grid, NISQ+PEC", DEFAULT.but(platform="superconducting",
                                      use_platform_clock=True),
     "nisq", curves.COLORS["nisq_pec"], "-", 0.30),
    ("atoms, STAR", DEFAULT.but(platform="neutral_atom", use_platform_clock=True,
                                encoding="jw", p=5e-3),
     "star", curves.COLORS["star"], "-", 0.55),
    ("atoms, surface code", DEFAULT.but(platform="neutral_atom",
                                        use_platform_clock=True,
                                        encoding="jw", p=5e-3),
     "surface", curves.COLORS["star"], (0, (3, 1.5)), 0.5),
    # O14: the mobile-qubit platforms propose HIGH-RATE qLDPC codes, not the
    # surface code. Crediting them with it is what gives Helios an FT arm at all.
    ("Helios, QLDPC", DEFAULT.but(platform="helios", use_platform_clock=True,
                                  encoding="jw", p=7.9e-4),
     "pinnacle", curves.COLORS["nisq_none"], (0, (5, 2)), 0.72),
]


def reach(cfg, arm, n):
    if arm == "surface":
        return ftqc.max_m_surface(n, cfg.but(pl_model="fowler"))
    if arm == "star":
        return ftqc.max_m_star(n, cfg)
    if arm == "pinnacle":
        return ftqc.max_m_pinnacle(n, cfg)
    return nisq.max_m(n, cfg, "pec")


def binding(cfg, arm, n):
    """0 = nothing feasible, 1 = qubit-limited, 2 = clock-limited."""
    m = reach(cfg, arm, n)
    if m <= 0:
        return 0
    if arm == "star":
        pt = ftqc.star_point(m, cfg)
        cp = math.floor(n / pt["phys"])
        frac = ftqc.n_shots_total(cfg, m) * pt["t_shot"] / max(cp, 1) / cfg.budget_s
        return 2 if frac > 0.5 else 1
    if arm == "surface":
        pt = ftqc.surface_point(m, cfg.but(pl_model="fowler"))
        cp = math.floor(n / pt["phys"])
        frac = ftqc.n_shots_total(cfg, m) * pt["t_shot"] / max(cp, 1) / cfg.budget_s
    elif arm == "pinnacle":
        pt = ftqc.pinnacle_point(m, cfg)
        cp = math.floor(n / pt["phys"])
        frac = ftqc.n_shots_total(cfg, m) * pt["t_shot"] / max(cp, 1) / cfg.budget_s
    else:
        cp = math.floor(n / hubbard.counts(m, cfg)["q_per_copy"])
        frac = nisq.time_required(m, cfg, "pec") / max(cp, 1) / cfg.budget_s
    return 2 if frac > 0.5 else 1


fig = plt.figure(figsize=(7.0, 3.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.30, 0.92], wspace=0.06)
axA, axB = (fig.add_subplot(g) for g in gs)

# ---------------------------------------------------------------- (a) reach
ax = axA
ax.axhspan(BAND_LO, BAND_HI, color=F.PALEGREY, alpha=0.55, lw=0, zorder=0)
ax.text(1.2e3, BAND_HI * 1.25, f"estimated ED capacity ($m\\leq{BAND_HI:.0f}$)",
        fontsize=5.6, color=F.DARKGREY, va="bottom", ha="left")
for label, cfg, arm, col, ls, frac in ARMS:
    y = np.array([reach(cfg, arm, n) for n in n_grid], float)
    if not (y > 0).any():
        continue
    ax.plot(n_grid, np.where(y > 0, y, np.nan), color=col, lw=1.5, ls=ls, zorder=3)
    # label at ~62% along the curve's own visible range, inside the axes, so
    # nothing collides with panel (b)'s row labels on the right
    live = np.flatnonzero(y > 0)
    k = live[int(frac * (len(live) - 1))]
    ax.annotate(label, xy=(n_grid[k], y[k]), xytext=(0, 4),
                textcoords="offset points", fontsize=5.2, color=col,
                ha="center", va="bottom")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1e3, 1e12); ax.set_ylim(2.5, 3e3)
ax.set_xlabel(r"$n$ (physical qubits)")
ax.set_ylabel(r"$m$ (resolvable lattice sites)")
ax.text(0.97, 0.94, "(a)", transform=ax.transAxes, ha="right", va="top", fontsize=7)

# arms that never appear, and why
missing = [(lb, cf, a) for lb, cf, a, _c, _s, _f in ARMS
           if not any(reach(cf, a, n) > 0 for n in n_grid)]
notes = ["transversal FT gives slow clocks $\\Theta(1)$ syndrome rounds,",
         "not $\\Theta(d)$: ${\\sim}19{\\times}$ here — enough for STAR, not the surface code"]
ax.text(1.2e3, 2.75, "\n".join(notes), fontsize=4.9, color=F.DARKGREY,
        va="bottom", ha="left", linespacing=1.3,
        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=0.9))

# ------------------------------------------------- (b) which constraint binds
ax = axB
rows = [(lb, cf, a, col) for lb, cf, a, col, _s, _f in ARMS]
for i, (label, cfg, arm, col) in enumerate(rows):
    b = np.array([binding(cfg, arm, n) for n in n_grid])
    for val, hatch, alpha in ((1, "", 0.30), (2, "///", 0.75)):
        mask = b == val
        if mask.any():
            ax.fill_between(n_grid, i - 0.34, i + 0.34, where=mask, color=col,
                            alpha=alpha, lw=0, hatch=hatch, edgecolor=col)
ax.set_xscale("log"); ax.set_xlim(1e3, 1e12)
ax.set_yticks(range(len(rows)))
ax.yaxis.tick_right()
ax.set_yticklabels([r[0] for r in rows], fontsize=5.2)
ax.tick_params(axis="y", length=0, pad=2)
ax.set_ylim(-0.85, len(rows) - 0.4)
ax.set_xlabel(r"$n$ (physical qubits)")
ax.text(0.03, 0.94, "(b)", transform=ax.transAxes, ha="left", va="top", fontsize=7)
ax.text(1.2e3, -0.52, "solid: qubit-limited    hatched: clock-limited    "
                      "blank: nothing feasible",
        fontsize=4.9, color=F.DARKGREY, va="center", ha="left")

fig.text(0.995, 0.004, f"model {model_id()}", ha="right", va="bottom",
         fontsize=4.2, color="0.55")
out = pathlib.Path(__file__).parent / "figures"
out.mkdir(exist_ok=True)
fig.tight_layout(pad=0.5)
for ext in ("pdf", "png"):
    fig.savefig(out / f"fh_platforms.{ext}", dpi=300)
print(f"wrote {out}/fh_platforms.pdf and .png\n")

print(f"{'arm':<24}" + "".join(f"{f'1e{e}':>9}" for e in (6, 8, 10, 12)))
for label, cfg, arm, _c, _s, _f in ARMS:
    print(f"{label:<24}" + "".join(f"{reach(cfg, arm, 10.0**e):>9.1f}"
                                   for e in (6, 8, 10, 12)))
print(f"\n{'arm':<24}binding constraint by n (1=qubits 2=clock 0=infeasible)")
for label, cfg, arm, _c, _s, _f in ARMS:
    print(f"{label:<24}" + "".join(f"{binding(cfg, arm, 10.0**e):>9d}"
                                   for e in (6, 8, 10, 12)))

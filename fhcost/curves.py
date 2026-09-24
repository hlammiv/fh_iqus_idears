"""Evaluates every curve on a common n grid, and locates the crossings."""
from __future__ import annotations
import math
import numpy as np
from .budget import Config, DEFAULT
from . import nisq, ftqc, classical, hubbard

# colours validated for CVD by check_palette.py -- all 10 pairs pass
COLORS = {
    "ideal":        "#333333",   # INK   (dashed reference, not a data series)
    "nisq_none":    "#56B4E9",   # SKY
    "nisq_pec":     "#0072B2",   # BLUE
    "nisq_pec_mpf": "#0072B2",   # BLUE, dashed
    "star":         "#D55E00",   # ORANGE
    "surface":      "#009E73",   # GREEN
    "pinnacle":     "#E69F00",   # AMBER (Okabe-Ito) -- QLDPC arm; CVD-validated
    "surface_mpf":  "#009E73",   # GREEN, dashed
}

LABELS = {
    "ideal":        r"$p=0$ (qubit-limited)",
    "nisq_none":    r"NISQ, unmitigated",
    "nisq_pec":     r"NISQ + PEC",
    "nisq_pec_mpf": r"NISQ + PEC + MPF-4",
    "star":         r"STAR",
    "surface":      r"Surface-code FT",
    "pinnacle":     "Pinnacle (QLDPC)\n2 coupler layers",
    "surface_mpf":  r"Surface-code FT + MPF-4",
}

MPF_ORDER = 2      # k=2 -> order-4 multiproduct

# Pinnacle is costed on DIFFERENT HARDWARE from every other curve, and now says
# so instead of burying it in a docstring. Generalised bicycle codes are not
# embeddable in the slide's nearest-neighbour grid; Bravyi et al.
# (arXiv:2308.07915) put the requirement at degree-6 with two edge-disjoint
# planar subgraphs, i.e. two coupler layers. On the slide's grid this arm
# returns 0 -- not "small", unavailable. Keeping it on the figure under its own
# hardware, clearly labelled, is more informative than dropping it and more
# honest than pretending it runs on the grid.
PINNACLE_PLATFORM = "sc_long_range"

# The Trotter step count is the dominant uncertainty, so every arm is reported
# across TWO SCENARIOS rather than as a point value:
#   "bound"     Campbell's worst-case commutator bound
#   "measured"  the exact-diagonalisation calibration, extrapolated beyond its
#               domain of m <= 12 and tau <= 2
# These are SCENARIOS, not an uncertainty interval: neither is a guaranteed edge,
# the truth is not required to lie between them, and the shaded region carries no
# confidence level. See OPEN_ITEMS.md O9 and hubbard.calibration_status.
TROTTER_EDGES = ("extensive", "measured")


def evaluate_band(n_grid, cfg: Config = DEFAULT) -> dict:
    """{arm: (lo, hi)} over the Trotter-model uncertainty."""
    out = {}
    lo = evaluate(n_grid, cfg.but(trotter=TROTTER_EDGES[0]))
    hi = evaluate(n_grid, cfg.but(trotter=TROTTER_EDGES[1]))
    for k in lo:
        out[k] = (np.minimum(lo[k], hi[k]), np.maximum(lo[k], hi[k]))
    return out


def evaluate(n_grid, cfg: Config = DEFAULT) -> dict:
    mpf = cfg.but(trotter_order_k=MPF_ORDER)
    fow, wil = cfg.but(pl_model="fowler"), cfg.but(pl_model="willow")
    mpf_fow = mpf.but(pl_model="fowler")
    out = {k: [] for k in
           ("ideal", "nisq_none", "nisq_pec", "nisq_pec_mpf", "star",
            "surface", "surface_willow", "surface_mpf", "pinnacle")}
    for n in n_grid:
        out["ideal"].append(nisq.max_m_ideal(n, cfg))
        out["nisq_none"].append(nisq.max_m(n, cfg, "none"))
        out["nisq_pec"].append(nisq.max_m(n, cfg, "pec"))
        out["nisq_pec_mpf"].append(nisq.max_m(n, mpf, "pec"))
        out["star"].append(ftqc.max_m_star(n, cfg))
        out["surface"].append(ftqc.max_m_surface(n, fow))
        out["surface_willow"].append(ftqc.max_m_surface(n, wil))
        out["surface_mpf"].append(ftqc.max_m_surface(n, mpf_fow))
        out["pinnacle"].append(
            ftqc.max_m_pinnacle(n, cfg.but(platform=PINNACLE_PLATFORM)))
    return {k: np.array(v, float) for k, v in out.items()}


M_FLOOR = 4.0      # smallest real lattice: a 2x2. Below it there is nothing to run


def band_for_plot(lo, hi, y_floor: float, m_floor: float = M_FLOOR):
    """(lo_drawn, hi_drawn, off_scale) for a SCENARIO band on a log axis.

    A zero edge is not missing data. It says the pessimistic scenario does not
    reach even a 2x2 lattice -- which is the most informative thing the band has
    to say, and NaN-masking it made the band vanish precisely there. On the
    default configuration the NISQ+PEC band had a zero lower edge at EVERY
    plotted n, so no shading was drawn at all (review #9).

    The lower edge is therefore clipped to the axis floor and the clipped region
    returned as `off_scale`, so the figure can hatch it and say what it means.
    Where the UPPER edge is zero too, nothing is feasible and both edges come
    back NaN -- there is genuinely no band, and the figure annotates instead.
    """
    lo = np.asarray(lo, float)
    hi = np.asarray(hi, float)
    dead = hi <= 0.0                       # nothing feasible at all
    off = (lo <= 0.0) & ~dead              # pessimistic edge below m_floor
    lo_d = np.where(off, y_floor, lo)
    lo_d = np.where(dead, np.nan, lo_d)
    hi_d = np.where(dead, np.nan, hi)
    return lo_d, hi_d, off


CROSSING_LIMIT = 1e10


def crossing(f, target: float, lo: float = 1e2, hi: float = CROSSING_LIMIT) -> float | None:
    """Smallest n with f(n) >= target, bisected in log n. None if never."""
    if f(hi) < target:
        return None            # NOT reached within the search range -- see fmt_crossing
    if f(lo) >= target:
        return lo
    for _ in range(70):
        mid = math.sqrt(lo * hi)
        lo, hi = (lo, mid) if f(mid) >= target else (mid, hi)
    return math.sqrt(lo * hi)


def extrapolation_table(n: float = 1e6, cfg: Config = DEFAULT) -> list[tuple]:
    """The user's question: what does each extrapolation buy, at fixed n?

    (a) noise extrapolation  -- ZNE(k) and its optimal limit, PEC
    (b) Trotter-step (multiproduct/Richardson-in-dt) extrapolation, NET of ||c||_1^2
    (c) finite-size extrapolation -- handled separately; see METHODS.md, it changes
        WHICH m you need rather than which m you can run.
    """
    rows = []
    base = nisq.max_m(n, cfg, "none")
    rows.append(("NISQ, unmitigated", base, "baseline"))
    for k in (1, 2):
        rows.append((f"  + ZNE order {k}", nisq.max_m(n, cfg, f"zne{k}"), "(a) noise"))
    pec = nisq.max_m(n, cfg, "pec")
    rows.append(("  + PEC (as implemented)", pec, "(a) noise"))
    for k in (2, 3, 4):
        c = cfg.but(trotter_order_k=k)
        rows.append((f"  + PEC + MPF order {2*k}", nisq.max_m(n, c, "pec"),
                     f"(b) Trotter, ||c||_1={hubbard.multiproduct_l1(k):.2f}"))
    ft = ftqc.max_m_surface(n, cfg.but(pl_model="fowler"))
    rows.append(("Surface-code FT", ft, "baseline"))
    for k in (2, 3, 4):
        c = cfg.but(trotter_order_k=k, pl_model="fowler")
        rows.append((f"  + MPF order {2*k}", ftqc.max_m_surface(n, c),
                     f"(b) Trotter, ||c||_1={hubbard.multiproduct_l1(k):.2f}"))
    return rows


def fmt_crossing(v) -> str:
    """Report an unfound crossing honestly: the search was bounded."""
    return f"not reached below n = {CROSSING_LIMIT:.0e}" if v is None else f"{v:.2e}"


# ---------------------------------------------------------------------------
# WHICH INTEGER LATTICES CAN ACTUALLY HOLD THE STATE (second-pass review #10)
#
# floor(sqrt(m)) is not enough. The specified initial state is half-filled with
# S^z_tot = 0: one holon, one doublon, and S^z = 0 triplets covering every
# remaining site in pairs. Two constraints follow, and a 5x5 satisfies neither:
#
#   1. N_sites must be EVEN. Each triplet contributes one up and one down, and
#      so does the doublon, giving N_up = (N-2)/2 + 1 = N/2. On 25 sites that is
#      not an integer.
#   2. The N-2 sites left after the defects must admit a PERFECT DIMER COVERING.
#      On 25 sites, 23 remain -- odd, so no covering exists at all.
#
# Both are satisfiable exactly when N is even. A grid graph is bipartite, and
# N even forces at least one side even, which makes the two colour classes
# equal; putting the holon and the doublon on OPPOSITE sublattices then leaves
# equal classes, and Gomory's theorem gives a perfect matching of the rest. So
# the admissibility rule is "N even, defects on opposite sublattices", and the
# defect placement is ours to choose.
#
# Restricting to squares would be needlessly harsh -- the experiment itself uses
# 7x4 -- so rectangles are admitted and the most square one wins on a tie.
LATTICE_SPEC = ("half filling, S^z_tot = 0, one holon + one doublon, "
                "S^z = 0 triplets on a perfect dimer covering of the rest; "
                "nearest-neighbour grid, doubly periodic")


MAX_ASPECT = 2.0   # beyond this a grid is a quasi-1D ribbon, not a 2D lattice
                   # -- and a materially easier classical problem


def lattice_admissible(lx: int, ly: int,
                       max_aspect: float = MAX_ASPECT) -> tuple[bool, str]:
    """Can an lx x ly grid hold the specified initial state, and is it 2D?"""
    n = lx * ly
    if lx < 2 or ly < 2:
        return False, "a lattice needs at least two sites on a side to have a bond"
    if n % 2:
        return False, (f"{n} sites is odd: half filling with S^z_tot = 0 needs "
                       f"N_up = N/2, and {n} - 2 = {n - 2} sites cannot be "
                       f"perfectly dimer-covered either")
    a = max(lx, ly) / min(lx, ly)
    if a > max_aspect + 1e-12:
        return False, (f"aspect {a:.1f} > {max_aspect}: a ribbon this thin is "
                       f"quasi-1D, and a different -- easier -- problem")
    return True, ""


def admissible_lattices(m: float, square_only: bool = False,
                        max_aspect: float = MAX_ASPECT):
    """[(lx, ly, sites)] fitting capacity m, largest first, most square first."""
    lim = int(math.floor(m))
    out = []
    for lx in range(2, lim + 1):
        if lx * lx > lim:
            break
        for ly in range(lx, lim // lx + 1):
            if square_only and lx != ly:
                continue
            ok, _ = lattice_admissible(lx, ly, max_aspect)
            if ok and lx * ly <= lim:
                out.append((lx, ly, lx * ly))
    out.sort(key=lambda t: (-t[2], t[1] - t[0]))
    return out


def best_lattice(m: float, square_only: bool = False,
                 max_aspect: float = MAX_ASPECT):
    """(lx, ly, sites) -- the largest admissible lattice, or None."""
    a = admissible_lattices(m, square_only, max_aspect)
    return a[0] if a else None


def fmt_lattice(m: float, square_only: bool = False) -> str:
    b = best_lattice(m, square_only)
    return "none" if b is None else f"{b[0]}x{b[1]}"


def max_integer_L(m: float) -> int:
    """DEPRECATED: floor(sqrt(m)) ignores the initial state. Use best_lattice.

    Kept only so an old caller fails loudly rather than silently reporting an
    inadmissible size; selftest asserts it disagrees with best_lattice where the
    review said it does.
    """
    return int(math.floor(math.sqrt(max(m, 0.0))))


def summary(cfg: Config = DEFAULT) -> dict:
    b = classical.band(cfg)
    lo, hi = b["band"]
    mpf = cfg.but(trotter_order_k=MPF_ORDER)
    fow = cfg.but(pl_model="fowler")
    wil = cfg.but(pl_model="willow")
    pin = cfg.but(platform=PINNACLE_PLATFORM)
    return {
        "classical_band": (lo, hi),
        "n_ft_clears_classical_lo": crossing(lambda n: ftqc.max_m_surface(n, fow), lo),
        "n_ft_clears_classical_hi": crossing(lambda n: ftqc.max_m_surface(n, fow), hi),
        "n_ft_willow_clears_hi": crossing(lambda n: ftqc.max_m_surface(n, wil), hi),
        "n_pec_clears_classical_hi": crossing(lambda n: nisq.max_m(n, cfg, "pec"), hi),
        "n_pec_mpf_clears_classical_hi": crossing(lambda n: nisq.max_m(n, mpf, "pec"), hi),
        "n_star_clears_classical_hi": crossing(lambda n: ftqc.max_m_star(n, cfg), hi),
        "n_ft_overtakes_pec": crossing(
            lambda n: (ftqc.max_m_surface(n, fow) - nisq.max_m(n, cfg, "pec"))
            if ftqc.max_m_surface(n, fow) > 0 else -1.0, 1e-9),
        "pec_at_1e6": nisq.max_m(1e6, cfg, "pec"),
        "pec_mpf_at_1e6": nisq.max_m(1e6, mpf, "pec"),
        "star_at_1e6": ftqc.max_m_star(1e6, cfg),
        "ft_fowler_at_1e6": ftqc.max_m_surface(1e6, fow),
        "ft_willow_at_1e6": ftqc.max_m_surface(1e6, wil),
        "ft_mpf_at_1e6": ftqc.max_m_surface(1e6, mpf.but(pl_model="fowler")),
        "pinnacle_at_1e6": ftqc.max_m_pinnacle(1e6, pin),
        "n_pinnacle_clears_classical_hi": crossing(
            lambda n: ftqc.max_m_pinnacle(n, pin), hi),
        "pinnacle_on_slide_grid": ftqc.max_m_pinnacle(1e6, cfg),
    }


if __name__ == "__main__":
    s = summary()
    print("CLASSICAL BAND      m = %.0f .. %.0f\n" % s["classical_band"])
    for k, v in s.items():
        if k == "classical_band":
            continue
        print(f"  {k:<32} {'never' if v is None else f'{v:.4g}'}")
    print(f"\nEXTRAPOLATION LADDER at n = 1e6 (classical band "
          f"{s['classical_band'][0]:.0f}-{s['classical_band'][1]:.0f})\n")
    print(f"  {'strategy':<32}{'m':>8}   which extrapolation")
    for name, m, tag in extrapolation_table(1e6):
        flag = "  <-- clears classical" if m > s["classical_band"][1] else ""
        print(f"  {name:<32}{m:>8.1f}   {tag}{flag}")


# ---------------------------------------------------------------------------
# Which arm leads, and where they cross. check_docs.py verifies NUMBERS quoted
# in prose; it cannot see a claim like "STAR beats NISQ at every n". Those went
# stale silently when c_rot was measured -- STAR's crossover with NISQ moved from
# 5e4 to 1.3e5, its crossover with the surface code from 2e5 to 1.6e7, and
# Pinnacle stopped being the strongest arm above 1e10. So the orderings are
# computed here and asserted, like everything else.

ARMS = ("nisq_pec", "star", "surface", "pinnacle")


def arm_reach(name: str, n: float, cfg: Config = DEFAULT,
              own_hardware: bool = True) -> float:
    """One arm's m at budget n.

    `own_hardware` puts each arm on the hardware it actually requires, which is
    how the figure plots them: Fowler p_L for the surface code, and the
    two-coupler-layer chip the generalised bicycle codes need for Pinnacle.
    Pass False when the CALLER has already chosen a platform -- otherwise a
    Helios study silently gets its Pinnacle arm costed on sc_long_range, which
    is what the first version of this function did.
    """
    from . import nisq as _n, ftqc as _f
    if name == "nisq_pec":
        return _n.max_m(n, cfg, "pec")
    if name == "star":
        return _f.max_m_star(n, cfg)
    if name == "surface":
        return _f.max_m_surface(n, cfg.but(pl_model="fowler")
                                if own_hardware else cfg)
    if name == "pinnacle":
        return _f.max_m_pinnacle(n, cfg.but(platform=PINNACLE_PLATFORM)
                                 if own_hardware else cfg)
    raise ValueError(f"unknown arm {name!r}")


def crossover(a: str, b: str, cfg: Config = DEFAULT,
              lo: float = 1e2, hi: float = 1e13, rising: bool = True,
              own_hardware: bool = True):
    """Smallest n where `a` overtakes `b` (or falls behind, if rising=False).

    Bisection is not safe here -- the ordering is NOT monotone (surface leads
    Pinnacle below ~1e5, trails it to ~1e10, leads again after) -- so this scans.
    """
    import math
    want = (lambda x, y: x > y) if rising else (lambda x, y: x < y)
    prev = None
    e, e_hi, step = math.log10(lo), math.log10(hi), 0.02
    while e <= e_hi:
        n = 10.0 ** e
        cur = want(arm_reach(a, n, cfg, own_hardware),
                   arm_reach(b, n, cfg, own_hardware))
        if prev is False and cur:
            return n
        prev = cur
        e += step
    return None


def crossings(a: str, b: str, cfg: Config = DEFAULT, lo: float = 1e2,
              hi: float = 1e13, own_hardware: bool = True) -> list[tuple]:
    """EVERY n where the `a`/`b` ordering flips, as (n, "a>b" | "a<b").

    Reporting only the first is how "the surface code overtakes Pinnacle by
    n = 1e8" survived: surface leads below ~2.5e4, trails to ~1e10, then leads
    again. A single crossover number cannot describe that.
    """
    import math
    out, prev = [], None
    e, e_hi = math.log10(lo), math.log10(hi)
    while e <= e_hi:
        n = 10.0 ** e
        cur = arm_reach(a, n, cfg, own_hardware) > arm_reach(b, n, cfg, own_hardware)
        if prev is not None and cur != prev:
            out.append((n, f"{a}>{b}" if cur else f"{a}<{b}"))
        prev = cur
        e += 0.02
    return out


def leader(n: float, cfg: Config = DEFAULT) -> str:
    """Which arm reaches furthest at this budget."""
    return max(ARMS, key=lambda k: arm_reach(k, n, cfg))


def clears_at(name: str, cfg: Config = DEFAULT, lo: float = 1e2,
              hi: float = 1e13):
    """Smallest n at which this arm exceeds the top of the classical band."""
    import math
    from . import classical as _c
    top = _c.band(cfg)["band"][1]
    e, e_hi = math.log10(lo), math.log10(hi)
    while e <= e_hi:
        n = 10.0 ** e
        if arm_reach(name, n, cfg) > top:
            return n
        e += 0.02
    return None

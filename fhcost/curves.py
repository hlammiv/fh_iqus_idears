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
    "pinnacle":     r"Pinnacle (QLDPC)",
    "surface_mpf":  r"Surface-code FT + MPF-4",
}

MPF_ORDER = 2      # k=2 -> order-4 multiproduct


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
        out["pinnacle"].append(ftqc.max_m_pinnacle(n, cfg))
    return {k: np.array(v, float) for k, v in out.items()}


def crossing(f, target: float, lo: float = 1e2, hi: float = 1e10) -> float | None:
    """Smallest n with f(n) >= target, bisected in log n. None if never."""
    if f(hi) < target:
        return None
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


def summary(cfg: Config = DEFAULT) -> dict:
    b = classical.band(cfg)
    lo, hi = b["band"]
    mpf = cfg.but(trotter_order_k=MPF_ORDER)
    fow = cfg.but(pl_model="fowler")
    wil = cfg.but(pl_model="willow")
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
        "pinnacle_at_1e6": ftqc.max_m_pinnacle(1e6, cfg),
        "n_pinnacle_clears_classical_hi": crossing(
            lambda n: ftqc.max_m_pinnacle(n, cfg), hi),
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

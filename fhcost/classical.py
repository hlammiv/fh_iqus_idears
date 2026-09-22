"""How big an m the classical data center resolves -- the band the quantum curves must clear.

WHY THIS IS A HORIZONTAL BAND AND NOT A CURVE
  The x axis of the figure is physical QUBITS. Classical compute does not live
  on that axis, so the classical frontier is drawn as a horizontal band: the
  largest m a one-week run on a 1-100 PB data center resolves, independent of n.
  Everything below the band is classically easy and cannot be an advantage.

WHY EVERY CLASSICAL METHOD IS EXPONENTIAL IN m HERE
  At t_max = sqrt(m)/v the light cone has crossed the lattice, and that kills
  every structural shortcut at once:
    * Krylov / state vector:  dim = C(m, m/2)^2 ~ 4^m / m. Memory-bound.
    * snake MPS / PEPS:       entanglement across a length-L cut grows as
                              ~ ent_rate * v * t * L; at t = L/v that is
                              S ~ ent_rate * m bits, so chi ~ 2^{c m}. No
                              asymptotic win over ED -- a huge win at fixed
                              SHORT t, which is why the t_max panel exists.
    * light-cone cluster:     cone has (2 v t + 1)^2 = ~4m sites > m. No win.
  So the band is set by exact diagonalization, and tensor networks only matter
  in the short-time panel.

NOISE-INDUCED CLASSICAL SIMULABILITY (the second, dashed band)
  At fixed per-gate noise p, the output of an UNMITIGATED noisy circuit is
  classically simulable far more cheaply than the ideal one: Pauli-path /
  sparse-Pauli-dynamics methods truncate the Pauli weight because each gate
  damps high-weight paths by (1-p). Cost is quasi-polynomial in the circuit
  size at fixed p and fixed target accuracy [Aharonov, Gao, Landau, Liu &
  Vazirani, STOC 2023; Angrisani et al., arXiv:2403.13927; Schuster, Yin, Gao &
  Yao, arXiv:2407.12768]. This band therefore swallows the unmitigated-NISQ
  curve entirely. It is drawn DASHED because the published statements are
  asymptotic and the constants at these depths are not pinned down -- the
  honest position is "strongly suggestive, not a proof at m ~ 3".
"""
from __future__ import annotations
import math
from math import comb
from .budget import Config, DEFAULT
from .hubbard import t_max, cone_sites

BYTES_PER_AMP = 16.0      # complex128


def hilbert_dim(m: int) -> float:
    """Half-filled, Sz-resolved Hubbard sector: C(m, m/2)^2."""
    h = comb(m, m // 2)
    return float(h) * float(h)


def ed_max_m(cfg: Config = DEFAULT) -> int:
    """Largest m whose Krylov vectors fit in memory."""
    best = 0
    for m in range(2, 80, 2):
        if BYTES_PER_AMP * hilbert_dim(m) * cfg.n_krylov_vec > cfg.ram_bytes:
            break
        best = m
    return best


def ed_runtime_s(m: int, cfg: Config = DEFAULT) -> float:
    """One week of exaflops? Krylov matvecs x dim x nonzeros-per-row."""
    dim = hilbert_dim(m)
    t = t_max(m, cfg)
    matvecs = max(1.0, 2.0 * t * m * (4.0 + cfg.U_over_J)) * cfg.n_times
    return matvecs * dim * 4.0 * m * 2.0 / cfg.flops


def ed_frontier(cfg: Config = DEFAULT) -> int:
    """Largest m satisfying BOTH memory and the one-week clock."""
    m = ed_max_m(cfg)
    while m >= 2 and ed_runtime_s(m, cfg) > cfg.budget_s:
        m -= 2
    return max(m, 0)


def mps_bond_bits(m: float, cfg: Config = DEFAULT) -> float:
    """log2(chi) needed for a snake MPS of an L x L lattice at t_max."""
    L = math.sqrt(m)
    return cfg.ent_rate * cfg.v * t_max(m, cfg) * L


def mps_max_m(cfg: Config = DEFAULT) -> float:
    """Largest m a snake-MPS run finishes in a week, in memory AND in flops.

    Runtime is the binding constraint, not memory: a bond update costs
    O(chi^3), and there are ~ (Trotter steps) x m of them. Checking memory
    alone overstates the reachable m by nearly 2x.
    """
    from .hubbard import trotter_steps
    best = 0.0
    m = 2.0
    log2_ram = math.log2(cfg.ram_bytes)
    log2_flop_budget = math.log2(cfg.flops * cfg.budget_s)
    while m < 1e5:
        chi_bits = mps_bond_bits(m, cfg)
        log2_mem = 2 * chi_bits + math.log2(math.sqrt(m) * m * 4.0 * BYTES_PER_AMP)
        n_updates = trotter_steps(m, cfg) * m * cfg.n_times
        log2_flops = 3 * chi_bits + math.log2(max(n_updates, 1.0)) + math.log2(64.0)
        if log2_mem > log2_ram or log2_flops > log2_flop_budget:
            break
        best = m
        m += 2.0
    return best


def cluster_max_m(cfg: Config = DEFAULT) -> float:
    """Light-cone / cluster truncation: exact on (2vt+1)^2 sites, 4 states each."""
    best = 0.0
    m = 2.0
    while m < 1e5:
        nc = cone_sites(m, t_max(m, cfg), cfg)
        if nc * 2.0 > math.log2(cfg.ram_bytes / BYTES_PER_AMP):   # 4^nc amplitudes
            break
        best = m
        m += 2.0
    return best


# Entanglement density across a cut at the light-cone-crossing time, in bits per
# site. The state is near-volume-law there; Hubbard saturates at 2 bits/site, and
# a generic quench reaches an O(1) fraction of that. This is the single most
# sensitive classical input, so the band is reported ACROSS the range.
ENT_RANGE = (0.3, 1.5)


def band(cfg: Config = DEFAULT) -> dict:
    """The classical band: the union over methods, RAM and entanglement density.

    Low edge = the most pessimistic defensible classical machine (1 PB ED, or a
    strongly entangled state for the tensor network). High edge = the most
    optimistic (100 PB, weakly entangled). Quantum only counts as advantage
    above the HIGH edge.
    """
    ed_lo = ed_frontier(cfg.but(ram_bytes=1e15))
    ed_hi = ed_frontier(cfg.but(ram_bytes=100e15))
    mps_lo = mps_max_m(cfg.but(ent_rate=ENT_RANGE[1]))
    mps_hi = mps_max_m(cfg.but(ent_rate=ENT_RANGE[0]))
    # A classical attacker picks the BEST method available, so each edge is a
    # max over methods; the edges differ only in how favourable the assumptions
    # are (RAM, entanglement density).
    lo = max(ed_lo, mps_lo)
    hi = max(ed_hi, mps_hi)
    return {"ed_lo_1PB": ed_lo, "ed_mid": ed_frontier(cfg), "ed_hi_100PB": ed_hi,
            "mps_strong_ent": mps_lo, "mps_mid": mps_max_m(cfg),
            "mps_weak_ent": mps_hi, "cluster": cluster_max_m(cfg),
            "band": (lo, hi)}


if __name__ == "__main__":
    print("exact diagonalization: half-filled Hubbard sector\n")
    print(f"{'m':>4} {'dim':>12} {'RAM 4 vec':>13} {'1-week runtime':>16}")
    for m in range(16, 34, 2):
        d = hilbert_dim(m)
        print(f"{m:>4} {d:>12.3e} {BYTES_PER_AMP*d*4/2**40:>10.3g} TB "
              f"{ed_runtime_s(m):>13.3g} s")
    print()
    for label, ram in (("1 PB", 1e15), ("10 PB (Frontier-class)", 10e15), ("100 PB", 100e15)):
        c = DEFAULT.but(ram_bytes=ram)
        print(f"  {label:24} memory-limited m = {ed_max_m(c):>3}   "
              f"after the 1-week clock: {ed_frontier(c):>3}")
    print("\nsnake MPS at t_max, vs entanglement density (bits/site)")
    for e in (0.3, 0.6, 1.0, 1.5):
        c = DEFAULT.but(ent_rate=e)
        print(f"  ent_rate = {e:.1f}  ->  m = {mps_max_m(c):>4.0f}   "
              f"(log2 chi = {mps_bond_bits(max(mps_max_m(c),2), c):.1f})")
    b = band()
    print(f"\nlight-cone cluster at t_max:   m = {b['cluster']:.0f}   (cone is ~4m sites > m)")
    print(f"exact diagonalization, 1-100 PB: m = {b['ed_lo_1PB']} .. {b['ed_hi_100PB']}")
    print(f"CLASSICAL BAND (union):          m = {b['band'][0]:.0f} .. {b['band'][1]:.0f}")


# --------------------------------------------------------------- fixed-t panel
UNBOUNDED = 1e9   # sentinel: m is not the limiting quantity at this t


def max_m_fixed_t(t: float, cfg: Config = DEFAULT) -> float:
    """Best classical reach at a FIXED evolution time t (not t = sqrt(m)/v).

    This is the panel that shows why the time-window convention decides the
    whole question. At small t the causal cone of <Z_i(t) Z_j(0)> contains only
    a handful of sites, so a cluster expansion is exact at a cost that does not
    grow with m AT ALL -- classical reach is unbounded. Slide 1's literal
    "t in [0, 1/m]" sits deep in that regime, which is why it cannot be the
    intended question.
    """
    from .hubbard import trotter_steps
    ncone = (2.0 * cfg.v * t + 1.0) ** 2
    # cluster expansion: 4^ncone amplitudes, independent of m
    if 2.0 * ncone <= math.log2(cfg.ram_bytes / BYTES_PER_AMP):
        return UNBOUNDED
    best = float(ed_frontier(cfg))
    c = cfg.but(tmax_mode="const", tmax_const=t)
    m = 2.0
    log2_ram = math.log2(cfg.ram_bytes)
    log2_flops = math.log2(cfg.flops * cfg.budget_s)
    while m < 1e5:
        chi_bits = cfg.ent_rate * cfg.v * t * math.sqrt(m)
        mem = 2 * chi_bits + math.log2(math.sqrt(m) * m * 4.0 * BYTES_PER_AMP)
        ops = 3 * chi_bits + math.log2(max(trotter_steps(m, c) * m * cfg.n_times, 1.0)) + 6.0
        if mem > log2_ram or ops > log2_flops:
            break
        best = max(best, m)
        m += 2.0
    return best

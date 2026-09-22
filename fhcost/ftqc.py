"""Surface-code fault tolerance and the STAR architecture.

THE CONSTRAINT THAT IS EASY TO FORGET
  An FT resource estimate is not a qubit count. Relative error eps on a signal
  of size s needs N_tot = T/(s*eps)^2 shots, and EVERY shot is a full circuit of
  r ~ 10^3-10^4 Trotter steps at d rounds per logical layer. Counting only
  qubits, and checking only a single shot against the week, overstates the
  reachable m by more than an order of magnitude and wrongly makes the FT curve
  slope 1 forever. With the shot budget imposed, the curve is

      m = min[ qubit-limited (slope 1) ,  shot/time-limited (slope 4/9) ]

  and the bend lands inside the plotted range.

p_L IS THE BIGGEST UNCERTAINTY ON THIS CURVE -- so it is a BAND, not a line
  "fowler"  p_L = 0.1 (p/1e-2)^{(d+1)/2}        idealised; suppression x10 per
            Fowler, Mariantoni, Martinis & Cleland, PRA 86, 032324 (2012)
  "willow"  p_L = 1.43e-3 * 2.14^{-(d-7)/2}     MEASURED; suppression x2.14 per
            Delta d = 2, anchored at 0.143%/round at d=7
            Google Quantum AI, Nature 638, 920 (2025)
  The two differ by ~4x in the reachable m. Plot both edges.

OTHER CORRECTIONS FOLDED IN
  * Encoding: FT uses JORDAN-WIGNER (2m+1 logical), not Derby-Klassen. Under
    lattice surgery the circuit is a sequence of pi/8 Pauli-PRODUCT rotations
    whose cost is independent of Pauli weight, and fermionic swap networks are
    Clifford, hence free in T-count. Paying 1.5 qubits/mode for locality buys
    nothing here. NISQ and STAR still use DK.  [Litinski, Quantum 3, 128 (2019)]
  * Hamming-weight phasing: k identical-angle rotations cost 4(k-1) T plus
    log2(k) syntheses instead of k full syntheses -- a ~30x T-count collapse,
    which is why the magic-state factory is NOT the bottleneck.
    [Gidney, Quantum 2, 74 (2018); used in Campbell, QST 7, 015007 (2021)]
  * d is chosen against the PER-SHOT logical volume. A logical fault biases the
    shot it occurs in; it does not compound across independent shots. Summing
    the volume over all ~10^5 shots over-provisions d badly.
  * STAR rotation error is P_Z,1 = 2p/15 per injected state and is INDEPENDENT
    of d -- raising d cannot fix it. gamma^2 = exp(8 P_Z,1 N) gives
    Lambda_STAR = 0.533 p N_rot.  [Akahoshi et al., PRX Quantum 5, 010337 (2024)]
"""
from __future__ import annotations
import math
from .budget import Config, DEFAULT
from .hubbard import (counts, step_depth, eps_absolute, multiproduct_l1,
                      multiproduct_branches, t_max, signal_at, conf_z)
from .nisq import M_MIN

ROSS_SELINGER = 3.0        # T gates per Rz = 3 log2(1/eps) (Ross & Selinger 2016)
FACTORY_QUBITS = 4620.0    # (15-to-1)_{17,7,7} at p=1e-3  [Litinski, Quantum 3, 205 (2019)]
FACTORY_CYCLES = 42.6      # rounds per output T state, same source
CULT_QUBITS = 1024.0       # magic-state cultivation source, 2e-9 output at p=1e-3
CULT_CYCLES = 100.0        # rounds per state  [Gidney, Shutty & Jones, arXiv:2409.17595]


def magic_cost(cfg: Config = DEFAULT) -> tuple[float, float]:
    """(qubits, rounds) per magic-state source. Cultivation is ~4.5x smaller."""
    if cfg.magic_source == "cultivation":
        return CULT_QUBITS, CULT_CYCLES
    if cfg.magic_source == "litinski":
        return FACTORY_QUBITS, FACTORY_CYCLES
    raise ValueError(f"unknown magic_source {cfg.magic_source!r}")


def star_theta(m: float, cfg: Config = DEFAULT) -> tuple[float, float]:
    """(sum of |rotation angles|, rotation count), angle-weighted convention.

    Termwise second-order evolution at U/J = 4: Theta = (5L^2 - 4L) tau and
    R = [16 L(L-1) + 2L^2 - 2] r + 1.  Their angle-dependent STAR law charges
    p per unit ANGLE rather than per rotation, which is the mechanism that makes
    STAR worth building -- most Trotter rotations are small-angle.
    """
    L = math.sqrt(m)
    r = counts(m, cfg)["steps"]
    theta = (5.0 * m - 4.0 * L) * t_max(m, cfg)
    rot = (16.0 * L * (L - 1.0) + 2.0 * m - 2.0) * r + 1.0
    return theta, max(rot, 1.0)


def p_logical(d: int, cfg: Config = DEFAULT) -> float:
    if cfg.pl_model == "fowler":
        return 0.1 * (cfg.p / cfg.p_th) ** ((d + 1) / 2.0)
    if cfg.pl_model == "willow":
        return 1.43e-3 * 2.14 ** (-(d - 7) / 2.0)
    if cfg.pl_model == "star_fit":      # Akahoshi et al. Eq. 18
        return 0.0679 * (cfg.p / 0.00385) ** ((d + 1) / 2.0)
    raise ValueError(f"unknown pl_model {cfg.pl_model!r}")


def storage_per_logical(d: int, cfg: Config = DEFAULT) -> float:
    return cfg.routing * 2.0 * d * d


def n_shots_total(cfg: Config = DEFAULT, m: float | None = None) -> float:
    """Shots for relative error eps on a signal of magnitude s, over T times.

    Multiproduct extrapolation amplifies shot noise by ||c||_1^2; that cost is
    charged here so the depth saving is reported NET, not gross.
    """
    # per-branch: no PEC here, so v_i = 1, but branch i still costs k_i x the
    # runtime, so the weight is (sum_i |c_i| sqrt(k_i))^2 rather than ||c||_1^2
    w = sum(abs(ci) * math.sqrt(ki) for ki, ci in
            multiproduct_branches(cfg.trotter_order_k))
    # same ledger as NISQ: a frac_stat share, as a confidence half-width.
    # This previously used the FULL tolerance, which is why FT ran 1.9x over.
    sig = cfg.s_sig if m is None else signal_at(m, cfg)
    delta = cfg.frac_stat * cfg.eps * max(sig, cfg.s_res_min) / conf_z(cfg)
    return cfg.n_times * w * w / delta ** 2


def t_counts(m: float, cfg: Config = DEFAULT) -> tuple[float, float]:
    """(total T gates, sequential T-layers) with Hamming-weight phasing."""
    c = counts(m, cfg)
    r = c["steps"]
    n_distinct = max(c["n_rot"] / max(m, 1.0), 1.0)      # distinct angles
    eps_syn = cfg.frac_syn * eps_absolute(cfg, m) / max(n_distinct, 1.0)
    n_syn = ROSS_SELINGER * math.log2(1.0 / eps_syn)
    lg = math.log2(max(m, 2.0))
    n_t = r * cfg.c_rot * (4.0 * (m - 1.0) + lg * n_syn)
    d_t = r * cfg.c_rot * (2.0 * lg + n_syn)
    return max(n_t, 1.0), max(d_t, 1.0)


def surface_point(m: float, cfg: Config = DEFAULT) -> dict | None:
    """Footprint and per-shot runtime for lattice size m under full FT."""
    q_L = 2.0 * m + cfg.n_ancilla              # Jordan-Wigner; see module docstring
    n_t, d_t = t_counts(m, cfg)
    eps_L = cfg.frac_logical * eps_absolute(cfg, m)
    for d in range(3, cfg.d_max, 2):
        rounds = d_t * d                        # PER SHOT
        if q_L * rounds * p_logical(d, cfg) > eps_L:
            continue
        # factories sized for throughput: n_t states in `rounds` rounds
        fq, fc = magic_cost(cfg)
        need_fac = max(1.0, math.ceil(n_t * fc / rounds))
        n_fac = max(float(cfg.n_factories), need_fac)
        phys = q_L * storage_per_logical(d, cfg) + n_fac * fq
        return {"m": m, "d": d, "q_L": q_L, "n_t": n_t, "d_t": d_t,
                "rounds": rounds, "n_fac": n_fac, "phys": phys,
                "t_shot": rounds * cfg.t_round}
    return None


def max_m_surface(n: float, cfg: Config = DEFAULT, m_hi: float = 1e6) -> float:
    def ok(m):
        n_tot = n_shots_total(cfg, m)
        pt = surface_point(m, cfg)
        if pt is None or pt["phys"] > n:
            return False
        copies = n / pt["phys"]
        if copies < 1:
            return False
        return n_tot * pt["t_shot"] / copies <= cfg.budget_s

    if not ok(M_MIN):
        return 0.0
    lo, hi = M_MIN, m_hi
    if ok(hi):
        return hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


# ----------------------------------------------------------------- STAR
def _star_lam(m: float, cfg: Config = DEFAULT) -> float:
    """Half the log sampling overhead: Gamma^2 = exp(2 * lam)."""
    if cfg.star_law == "angle":
        theta, rot = star_theta(m, cfg)
        return 0.5 * (4.0 * cfg.star_alpha * cfg.p * theta + 4.0 * cfg.star_floor * rot)
    if cfg.star_law == "count":
        return cfg.star_kappa * cfg.p * max(counts(m, cfg)["n_rot_cone"], 1.0)
    raise ValueError(f"unknown star_law {cfg.star_law!r}")


def star_point(m: float, cfg: Config = DEFAULT) -> dict | None:
    c = counts(m, cfg)
    q_L = 2.0 * m + cfg.n_ancilla
    n_rot = max(c["n_rot"], 1.0)
    eps_L = cfg.frac_logical * eps_absolute(cfg, m)
    for d in range(3, cfg.d_max, 2):
        rounds = c["steps"] * cfg.star_rounds_per_step * d
        if q_L * rounds * p_logical(d, cfg) > eps_L:
            continue
        return {"m": m, "d": d, "q_L": q_L, "n_rot": n_rot,
                # d-INDEPENDENT either way: raising the code distance cannot fix
                # an injected rotation error.
                "lam": _star_lam(m, cfg),
                "phys": q_L * storage_per_logical(d, cfg),      # no factories
                "t_shot": rounds * cfg.t_round}
    return None


def max_m_star(n: float, cfg: Config = DEFAULT, m_hi: float = 1e6) -> float:
    def ok(m):
        pt = star_point(m, cfg)
        if pt is None or pt["phys"] > n:
            return False
        copies = n / pt["phys"]
        if copies < 1:
            return False
        need_log = math.log(n_shots_total(cfg, m)) + 2 * pt["lam"]
        have_log = math.log(cfg.budget_s * copies / pt["t_shot"])
        return need_log <= have_log
    if not ok(M_MIN):
        return 0.0
    lo, hi = M_MIN, m_hi
    if ok(hi):
        return hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


if __name__ == "__main__":
    print(f"shots needed for eps_rel={DEFAULT.eps} on s={DEFAULT.s_sig}, "
          f"T={DEFAULT.n_times}: N = {n_shots_total():,.0f}\n")
    for model in ("fowler", "willow"):
        cfg = DEFAULT.but(pl_model=model)
        print(f"--- p_L model: {model} ---")
        print(f"{'n':>10} {'surface FT':>11} {'d':>4} {'STAR':>9} {'binding':>14}")
        for n in (1e4, 1e5, 1e6, 1e7, 1e8):
            ms = max_m_surface(n, cfg)
            pt = surface_point(ms, cfg) if ms else None
            if pt:
                copies = math.floor(n / pt["phys"])
                frac = n_shots_total(cfg) * pt["t_shot"] / max(copies, 1) / cfg.budget_s
                bind = "shots/clock" if frac > 0.5 else "qubits"
            else:
                bind = "-"
            print(f"{n:>10.0e} {ms:>11.4g} {(pt['d'] if pt else 0):>4} "
                  f"{max_m_star(n, cfg):>9.4g} {bind:>14}")
        print()


# ------------------------------------------------- Pinnacle (QLDPC) architecture
# Webster, Berent, Chandra, Hockings, Baspin, Thomsen, Smith & Cohen (Iceberg
# Quantum), "The Pinnacle Architecture: Reducing the cost of breaking RSA-2048 to
# 100 000 physical qubits using quantum LDPC codes", arXiv:2602.11457v2 (2026).
# refs/2602.11457_pinnacle.pdf
#
# Generalised bicycle codes [[n, k, d]], dt = d + 2 code cycles per logical cycle,
# and n_pb physical qubits per processing block (code block + 4 gadgets + 4 bridges).
# Their Table I:
GB_CODES = [   # (n,    k,  d,  dt, n_pb)
    (30,   8,  4,  6,  140),
    (62,  10,  6,  8,  244),
    (126, 12, 10, 12,  452),
    (254, 14, 16, 18,  860),
    (510, 16, 24, 26, 1620),
]
# Logical error per logical qubit per logical cycle. Their Eq. (13) is
# p_L = (A/k)(p/B)^(d/2+C); fitting their Table III gives C = 1/2 (the same form
# as the surface code) with A = 5.84 and a threshold B = 1.58%. The fit reproduces
# every entry of their Table III to within 26% across twelve orders of magnitude.
PIN_A, PIN_B = 5.84, 0.0158
# One magic engine per processing unit, itself a GB code block, delivering one
# |T> per logical cycle. Calibrated so the footprint reproduces their Table IV
# L = 8 point (19 kq): see pinnacle_validate().
PIN_ENGINE = 1.2


def p_logical_gb(k: int, d: int, cfg: Config = DEFAULT) -> float:
    return (PIN_A / k) * (cfg.p / PIN_B) ** ((d + 1) / 2.0)


def pinnacle_point(m: float, cfg: Config = DEFAULT) -> dict | None:
    """Our dynamics workload run on the Pinnacle architecture.

    Their own Fermi-Hubbard application (their Table IV) is GROUND-STATE ENERGY
    at 0.5% relative energy error -- a different workload from this dynamics
    problem. What is borrowed here is the architecture: the code family, the
    logical error model, the (d+2) code-cycle logical clock, and one magic engine
    per processing unit. Their logical-qubit count N = 2L^2 + 2 is the same
    Jordan-Wigner count we already use.
    """
    q_L = 2.0 * m + cfg.n_ancilla
    n_t, d_t = t_counts(m, cfg)
    eps_L = cfg.frac_logical * eps_absolute(cfg, m)
    for (n_code, k, d, dt, n_pb) in GB_CODES:
        units = math.ceil(q_L / k)
        # each unit consumes at most one T per logical cycle
        cycles = max(d_t, n_t / units) / max(cfg.lanes, 1.0)
        if q_L * cycles * p_logical_gb(k, d, cfg) > eps_L:
            continue
        phys = units * n_pb * (1.0 + PIN_ENGINE)
        return {"m": m, "d": d, "k": k, "units": units, "n_t": n_t,
                "cycles": cycles, "phys": phys,
                "t_shot": cycles * dt * cfg.t_round}
    return None


PIN_N_VALID = 1e10   # above this the 5-code family is exhausted; see pinnacle_point


def max_m_pinnacle(n: float, cfg: Config = DEFAULT, m_hi: float = 1e6) -> float:
    """Largest m in a week.

    CAVEAT: their published family has only five codes, topping out at d = 24.
    Past n ~ 1e10 the model is pinned at that code and the curve saturates
    artificially (m = 6178 at both n = 1e11 and 1e13). That is this model
    running out of table rows, not physics -- real GB families extend further.
    The plotted range (n <= 1e8) is well inside the valid region.
    """
    def ok(m):
        n_tot = n_shots_total(cfg, m)
        pt = pinnacle_point(m, cfg)
        if pt is None or pt["phys"] > n:
            return False
        copies = n / pt["phys"]
        return copies >= 1 and n_tot * pt["t_shot"] / copies <= cfg.budget_s

    if not ok(M_MIN):
        return 0.0
    lo, hi = M_MIN, m_hi
    if ok(hi):
        return hi
    for _ in range(34):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def pinnacle_validate() -> None:
    """Footprint cross-check against their Table IV (ground-state workload)."""
    print("Pinnacle footprint vs their Table IV (p = 1e-3, N = 2L^2 + 2 logical)")
    print(f"{'L':>4}{'N_log':>7}{'code':>16}{'units':>7}{'ours':>9}{'paper':>8}{'ratio':>7}")
    paper = {8: 19e3, 10: 25e3, 12: 35e3, 14: 45e3, 16: 58e3, 18: 71e3, 20: 87e3}
    for L, want in paper.items():
        N = 2 * L * L + 2
        n_code, k, d, dt, n_pb = GB_CODES[3]        # d = 16, their p=1e-3 choice
        units = math.ceil(N / k)
        got = units * n_pb * (1.0 + PIN_ENGINE)
        print(f"{L:>4}{N:>7}{f'[[{n_code},{k},{d}]]':>16}{units:>7}"
              f"{got/1e3:>8.0f}k{want/1e3:>7.0f}k{got/want:>7.2f}")

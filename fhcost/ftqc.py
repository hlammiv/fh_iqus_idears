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

# Cultivation ends by escaping into a d = 15 grafted matchable code, so one unit
# occupies a d = 15 patch. Its expected volume per ACCEPTED state, retries
# included, is read off Fig. 1 of arXiv:2409.17595 at ~3e4 qubit-rounds; that is
# a value read from a log-log scatter plot, good to about a factor of two, and
# it is the weakest number in this file. Cycles = volume / footprint.
CULT_FOOTPRINT = 2.0 * 15 * 15        # 450 qubits, one d = 15 patch
CULT_VOLUME_1E3 = 3.0e4               # qubit-rounds per accepted state, p = 1e-3
CULT_CYCLES_1E3 = CULT_VOLUME_1E3 / CULT_FOOTPRINT
# "a 2x noise strength improvement ... becomes a 50x logical error rate
# improvement and a 10x cost reduction" -- same construction, so the footprint
# is held and the volume divided.
CULT_CYCLES_5E4 = CULT_CYCLES_1E3 / 10.0

# ---------------------------------------------------------------------------
# MAGIC-STATE SOURCES -- published operating points, not a formula
#
# The previous ladder gave a two-level factory an output error of 35 (4.5e-8)^3
# = 3.2e-21 and charged it 42.6 cycles. Both are wrong, and wrong in the
# optimistic direction. The cubic law describes the suppression of INPUT-state
# error only; it says nothing about faults in the distillation circuitry, whose
# footprint and cycle count do not shrink when the inputs get cleaner. Litinski
# determines p_out NUMERICALLY (5-qubit density-matrix simulation including
# storage errors and faulty T measurements) and reports 4.5e-20 at 128 cycles
# for the comparable two-level protocol -- 14x worse in error and 3x slower.
#
# So the table below IS the model: every row is a published operating point at a
# stated physical error rate. Nothing is extrapolated between rows, and nothing
# is extrapolated in p.
#
#   LITINSKI  Quantum 3, 205 (2019), arXiv:1905.06903, Table 1.
#             "Cycles" is already per OUTPUT state, so the 20-to-4 protocols'
#             four outputs are divided in. The rejection rate is already inside
#             it too: the time cost is 6 d_m / (1 - p_fail), and the published
#             counts exceed 6 d_m by 0-2%, which IS p_fail.
#             Footprint cross-check: a (15-to-1)_{dX,dZ,dm} block costs
#             2 (dX + 4 dZ) 3 dX + 4 dm qubits, which reproduces 810, 1150,
#             2070 and 4620 to the table's rounding. selftest asserts this.
#
#   CULTIVATION  Gidney, Shutty & Jones, arXiv:2409.17595, Fig. 1-2. End-to-end
#             (grown) error 2e-9 at p = 1e-3 with a 99% DISCARD rate, and 4e-11
#             at p = 5e-4 with 90%. Their cost axis is expected volume in
#             qubit-rounds INCLUDING retries, so the discard rate is paid there.
#
# (name, p_phys, p_out, qubits, cycles per ACCEPTED state, family)
MAGIC_SOURCES = [
    # --- Litinski Table 1, p_phys = 1e-4 ---
    ("(15-to-1)_7,3,3",                      1e-4, 4.4e-8,    810.0,  18.1, "litinski"),
    ("(15-to-1)_9,3,3 small-footprint",      1e-4, 1.5e-9,    762.0,  36.2, "litinski"),
    ("(15-to-1)_9,3,3",                      1e-4, 9.3e-10,  1150.0,  18.1, "litinski"),
    ("(15-to-1)_11,5,5",                     1e-4, 1.9e-11,  2070.0,  30.0, "litinski"),
    ("(15-to-1)^4_9,3,3 x (20-to-4)_15,7,9", 1e-4, 2.4e-15, 16400.0,  90.3, "litinski"),
    ("(15-to-1)^4_9,3,3 x (15-to-1)_25,9,9", 1e-4, 6.3e-25, 18600.0,  67.8, "litinski"),
    # --- Litinski Table 1, p_phys = 1e-3 ---
    ("(15-to-1)_17,7,7",                     1e-3, 4.5e-8,   4620.0,  42.6, "litinski"),
    ("(15-to-1)_9,5,5 x (15-to-1)_21,9,11",  1e-3, 6.1e-10,  7780.0, 469.0, "litinski"),
    ("(15-to-1)^6_13,5,5 x (20-to-4)_23,11,13", 1e-3, 1.4e-10, 43300.0, 130.0, "litinski"),
    ("(15-to-1)^4_13,5,5 x (20-to-4)_27,13,15", 1e-3, 2.6e-11, 46800.0, 157.0, "litinski"),
    ("(15-to-1)^6_11,5,5 x (15-to-1)_25,11,11", 1e-3, 2.7e-12, 30700.0,  82.5, "litinski"),
    ("(15-to-1)^6_13,5,5 x (15-to-1)_29,11,13", 1e-3, 3.3e-14, 39100.0,  97.5, "litinski"),
    ("(15-to-1)^6_17,7,7 x (15-to-1)_41,17,17", 1e-3, 4.5e-20, 73400.0, 128.0, "litinski"),
    # --- cultivation ---
    ("cultivation d1=5",                     1e-3, 2.0e-9,  CULT_FOOTPRINT, CULT_CYCLES_1E3, "cultivation"),
    ("cultivation d1=5",                     5e-4, 4.0e-11, CULT_FOOTPRINT, CULT_CYCLES_5E4, "cultivation"),
]

LADDER_P = sorted({row[1] for row in MAGIC_SOURCES})     # 5e-4, 1e-3 ... and 1e-4


def factory_ladder(cfg: Config = DEFAULT) -> tuple[list, float | None, str]:
    """(rows, ladder p_phys, status) for this machine's physical error rate.

    There is no honest interpolation of a numerically simulated p_out in p, so
    nothing is interpolated: every row tabulated at a p_phys AT OR ABOVE cfg.p is
    admissible, because a factory characterised on noisier hardware also works on
    quieter hardware. That makes the ladder pessimistic off the tabulated points
    and exact on them. Above the largest tabulated p_phys there is no published
    operating point at all and the model REFUSES rather than extrapolating --
    which is why the FT arms vanish above p = 1e-3, and that is the honest answer
    rather than a modelling gap hidden behind a fitted curve.
    """
    rows = [r for r in MAGIC_SOURCES if r[1] >= cfg.p]
    if cfg.magic_source in ("litinski", "cultivation"):
        rows = [r for r in rows if r[5] == cfg.magic_source]
    if not rows:
        if cfg.p > max(LADDER_P):
            return [], None, (f"no published factory at p = {cfg.p:.1e}; Litinski "
                              f"and cultivation both stop at {max(LADDER_P):.0e}")
        return [], None, f"no {cfg.magic_source} source tabulated at p >= {cfg.p:.1e}"
    q = min(r[1] for r in rows)
    status = ("exact" if any(abs(r[1] - cfg.p) < 1e-18 for r in rows) else
              f"conservative: nearest tabulated p_phys is {q:.0e}, machine is "
              f"{cfg.p:.1e}")
    return rows, q, status


def select_factory(p_target: float, cfg: Config = DEFAULT,
                   n_t: float | None = None, rounds: float | None = None):
    """The source minimising TOTAL magic qubits, not the smallest single unit.

    Two corrections over the previous version, both from review #4:

      * it consults cfg.p. Before, the same specification came out at p = 1e-5,
        1e-3 and 3e-3, which cannot be read as a hardware sensitivity.
      * it optimises the whole magic plant. Picking the smallest unit and only
        then asking how many are needed is the wrong order: cultivation's 450
        qubits beat a 4620-qubit 15-to-1 block per unit, but at 67 cycles per
        accepted state against 42.6 it needs more units, and which wins depends
        on the T rate the circuit actually demands. With n_t and rounds given,
        the cost compared is units x footprint; without them it falls back to
        the unit footprint and says so.

    Returns (name, qubits, cycles, p_out, units) or None.
    """
    rows, _, _ = factory_ladder(cfg)
    ok = [r for r in rows if r[2] <= p_target]
    if not ok:
        return None

    def plant(r):
        if n_t is None or rounds is None or rounds <= 0:
            return r[3], 1.0                      # no throughput information
        u = max(1.0, math.ceil(n_t * r[4] / rounds))
        return u * r[3], u

    best = min(ok, key=lambda r: plant(r)[0])
    return (best[0], best[3], best[4], best[2], plant(best)[1])


def magic_cost(cfg: Config = DEFAULT) -> tuple[float, float]:
    """(qubits, cycles) of the cheapest source in this ladder, ignoring fidelity.
    Diagnostic only; select_factory() is what the model uses."""
    rows, _, _ = factory_ladder(cfg)
    if not rows:
        raise ValueError(f"no magic-state source tabulated at p = {cfg.p:.1e}")
    r = min(rows, key=lambda x: x[3])
    return r[3], r[4]


def hwp_workspace(m: float, cfg: Config = DEFAULT) -> float:
    """Logical qubits Hamming-weight phasing needs, which were never charged.

    Phasing k identical-angle rotations computes their Hamming weight into a
    ceil(log2 k)-qubit register and applies one rotation per register bit, so the
    T-count saving comes with a workspace cost. A phase-gradient register of
    n_syn qubits is also needed; it is catalytic, so it is paid once.
    [Gidney, Quantum 2, 74 (2018)]
    """
    c = counts(m, cfg)
    n_distinct = max(c["n_rot"] / max(m, 1.0), 1.0)
    eps_syn = cfg.frac_syn * eps_absolute(cfg, m) / max(n_distinct, 1.0)
    n_syn = ROSS_SELINGER * math.log2(1.0 / eps_syn)
    return math.ceil(math.log2(max(m, 2.0))) + math.ceil(n_syn)


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
        # A FIXED EXPERIMENTAL ANCHOR at the measured device error, not a
        # p-dependent family: this expression has no cfg.p in it. Sweeping p with
        # pl_model="willow" therefore moves nothing, and a p-sensitivity column
        # computed that way is fixed by construction rather than by physics.
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
    """Footprint and per-shot runtime for lattice size m under full FT.

    Review #7 added three things that were missing:
      * Hamming-weight phasing was taking its T-count saving without paying for
        its workspace or its phase-gradient register,
      * magic-state infidelity was never budgeted, so the model consumed T states
        up to 320x too noisy at large n with a fixed factory spec, and
      * a logical failure flips a +-1 outcome, biasing the estimator by up to
        TWICE the failure probability, not once.

    Second-pass review #4 added two more:
      * the SAME factor of two now applies to the magic budget. A faulty T state
        corrupts a +-1 measurement exactly as a logical failure does, so charging
        the logical channel 2 p_L and the magic channel 1 p_T was an inconsistency
        in the ledger, not a modelling choice. It costs a factor of two in the
        per-state target, which moves the factory one rung up the ladder.
      * the factory is chosen INSIDE the distance loop and by total plant size.
        The number of units needed depends on `rounds`, which depends on d, so
        choosing the source first and counting units afterwards optimises the
        wrong quantity.
    """
    q_L = 2.0 * m + cfg.n_ancilla + hwp_workspace(m, cfg)
    n_t, d_t = t_counts(m, cfg)
    eps_L = cfg.frac_logical * eps_absolute(cfg, m)
    # per-state magic target: union bound over n_t states, x2 for the sign flip
    p_target = cfg.frac_magic * eps_absolute(cfg, m) / (2.0 * max(n_t, 1.0))
    for d in range(3, cfg.d_max, 2):
        # PER SHOT. No `lanes` divisor here: the factory bank below is already
        # SIZED for throughput (the unit count delivers n_t states within
        # `rounds`), so dividing again would count the same parallelism twice.
        # `lanes` applies only where supply is genuinely serialised -- Pinnacle's
        # single engine.
        rounds = d_t * d
        if 2.0 * q_L * rounds * p_logical(d, cfg) > eps_L:  # failure -> bias is x2
            continue
        fac = select_factory(p_target, cfg, n_t=n_t, rounds=rounds)
        if fac is None:
            return None                 # no published source is clean enough
        fname, fq, fc, f_pT, units = fac
        n_fac = max(float(cfg.n_factories), units)
        phys = q_L * storage_per_logical(d, cfg) + n_fac * fq
        return {"m": m, "d": d, "q_L": q_L, "n_t": n_t, "d_t": d_t,
                "rounds": rounds, "n_fac": n_fac, "phys": phys,
                "factory": fname, "p_T": f_pT, "p_T_target": p_target,
                "magic_qubits": n_fac * fq,
                "workspace": hwp_workspace(m, cfg),
                "t_shot": rounds * cfg.t_round}
    return None


def max_m_surface(n: float, cfg: Config = DEFAULT, m_hi: float = 1e6) -> float:
    def ok(m):
        n_tot = n_shots_total(cfg, m)
        pt = surface_point(m, cfg)
        if pt is None or pt["phys"] > n:
            return False
        copies = math.floor(n / pt["phys"])
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
        copies = math.floor(n / pt["phys"])
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


def p_logical_gb(k: int, d: int, cfg: Config = DEFAULT) -> float:
    return (PIN_A / k) * (cfg.p / PIN_B) ** ((d + 1) / 2.0)


ENGINE_QUBITS = 4410.0     # one magic engine, from their Hubbard footprint formula


def pinnacle_footprint(m: float, cfg: Config = DEFAULT, d: int = 24) -> float:
    """Their published formula, exactly: n = 1620 ceil((L^2+1)/8) + 4410.

    Decoded: 1620 is n_pb for the d = 24 code [[510,16,24]]; ceil((L^2+1)/8) is
    ceil((2m+2)/16), the number of processing BLOCKS at k = 16 logical qubits
    each; and 4410 is ONE magic engine for the whole machine. One processing
    unit, one engine, no memory. Reproduces all seven rows of their Table IV to
    within 2% with no free parameter.
    """
    k, n_pb = {4: (8, 140), 6: (10, 244), 10: (12, 452), 16: (14, 860),
               24: (16, 1620)}[d]
    blocks = math.ceil((2.0 * m + 2.0) / k)
    return blocks * n_pb + cfg.pin_engines * ENGINE_QUBITS


def pinnacle_point(m: float, cfg: Config = DEFAULT) -> dict | None:
    """Our dynamics workload on the Pinnacle architecture.

    CORRECTED against the paper (review #6). Two structural errors before:
      * the magic engine was scaled with the block count; there is exactly ONE,
      * so T states arrive at one per logical cycle for the WHOLE machine, not
        one per block. The previous model divided the T-count by the block count
        and was therefore ~m times too fast on a T-heavy workload.

    CONNECTIVITY CAVEAT: generalised bicycle codes need non-local qLDPC
    connectivity. The slide specifies a nearest-neighbour 2D grid, which does not
    provide it. This arm is therefore costed under a different hardware
    assumption from every other curve on the figure, and the comparison is not
    like-for-like. See OPEN_ITEMS.md O11.
    """
    q_L = 2.0 * m + cfg.n_ancilla
    n_t, d_t = t_counts(m, cfg)
    eps_L = cfg.frac_logical * eps_absolute(cfg, m)
    for (n_code, k, d, dt, n_pb) in GB_CODES:
        blocks = math.ceil(q_L / k)
        # ONE engine (by default): T consumption is serialised across the machine
        cycles = max(d_t, n_t / max(cfg.pin_engines, 1)) / max(cfg.lanes, 1.0)
        if q_L * cycles * p_logical_gb(k, d, cfg) > eps_L:
            continue
        phys = blocks * n_pb + cfg.pin_engines * ENGINE_QUBITS
        return {"m": m, "d": d, "k": k, "blocks": blocks, "n_t": n_t,
                "cycles": cycles, "phys": phys,
                "t_shot": cycles * dt * cfg.t_round}
    return None


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
        copies = math.floor(n / pt["phys"])
        return copies >= 1 and n_tot * pt["t_shot"] / copies <= cfg.budget_s

    if not ok(M_MIN):
        return 0.0
    lo, hi = M_MIN, m_hi
    if ok(hi):
        return hi
    for _ in range(80):          # same iteration count as every other arm, so a
        mid = 0.5 * (lo + hi)    # parity check against the explorer is exact
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def pinnacle_validate(cfg: Config = DEFAULT) -> None:
    """Footprint check against their Table IV -- NO free parameters."""
    print("Pinnacle footprint vs Table IV, d = 24, one engine (p = 1e-3)")
    print(f"{'L':>4}{'m':>6}{'blocks':>8}{'ours':>9}{'paper':>8}{'ratio':>7}")
    paper = {8: 19e3, 10: 25e3, 12: 35e3, 14: 45e3, 16: 58e3, 18: 71e3, 20: 87e3}
    worst = 0.0
    for L, want in paper.items():
        got = pinnacle_footprint(L * L, cfg, d=24)
        worst = max(worst, abs(got / want - 1))
        print(f"{L:>4}{L*L:>6}{math.ceil((2*L*L+2)/16):>8}"
              f"{got/1e3:>8.0f}k{want/1e3:>7.0f}k{got/want:>7.3f}")
    print(f"worst deviation {100*worst:.1f}%  (rounding in their published table)")

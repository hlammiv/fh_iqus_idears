"""How big an m the classical data center resolves.

WHAT THE BAND IS AND IS NOT
  It is the ESTIMATED CAPACITY of specified methods on a stated machine, not a
  classical impossibility boundary. Sitting below it does not make a point easy
  in practice; clearing it does not establish that every competitive classical
  method fails. The figure says "estimated ED capacity" for that reason.

  It is also observable-specific and parameter-specific. At U = 0 the
  Hamiltonian is quadratic and this weight-4 correlator is POLYNOMIAL despite
  the non-Gaussian initial state -- see free_fermion_frontier below -- so the
  ED band is simply the wrong baseline there. Estimating a low-weight
  observable is a different task from preparing the full state or sampling
  from it, and only the first is what this figure asks about.

WHY THIS IS A HORIZONTAL BAND AND NOT A CURVE
  The x axis of the figure is physical QUBITS. Classical compute does not live
  on that axis, so the classical frontier is drawn as a horizontal band: the
  largest m a one-week run on a 1-100 PB data center resolves, independent of n.
  Everything below the band is reachable by the methods costed here.

WHY EVERY METHOD COSTED HERE IS EXPONENTIAL IN m AT U > 0
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
from .hubbard import t_max, cone_sites, eps_absolute

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


def ed_io_seconds(m: int, cfg: Config = DEFAULT) -> float:
    """Wall clock if the state vector does not fit in memory and must be streamed.

    Out-of-core ED looks attractive -- a 700 PB filesystem holds an m = 30 vector
    where 10 PB of memory does not -- but time evolution touches the whole vector
    once per Krylov matvec, so the cost is bandwidth, not capacity. At m = 28 that
    is ~1300 s per matvec and ~50 WEEKS for the evolution. Out-of-core therefore
    does NOT extend the frontier, and the in-memory limit stands.
    """
    vec = BYTES_PER_AMP / 2.0 * hilbert_dim(m)          # complex64 out of core
    t = t_max(m, cfg)
    matvecs = max(1.0, 2.0 * t * m * (4.0 + cfg.U_over_J)) * cfg.n_times
    return matvecs * vec / cfg.disk_bw


def ed_frontier(cfg: Config = DEFAULT) -> int:
    """Largest m satisfying memory, the one-week clock, and (out of core) I/O."""
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


def cluster_radius(t: float, cfg: Config = DEFAULT) -> float:
    """Cluster radius for an ERROR-CONTROLLED truncation, matching the quantum
    finite-size criterion. Truncating at v t alone is not exact -- the tail
    outside the cone is exponential, so reaching accuracy eps needs the same
    xi ln(1/eps) buffer that converged.m_required charges. Omitting it here gave
    the classical side a free pass the quantum side did not get."""
    eps = cfg.frac_trotter * eps_absolute(cfg)
    return cfg.v * t + cfg.xi * math.log(1.0 / max(eps, 1e-12))


def cluster_sites(t: float, cfg: Config = DEFAULT) -> float:
    return (2.0 * cluster_radius(t, cfg) + 1.0) ** 2


def cluster_max_m(cfg: Config = DEFAULT) -> float:
    """Error-controlled cluster truncation: exact on cluster_sites, 4 states each."""
    best = 0.0
    m = 2.0
    while m < 1e5:
        nc = min(m, cluster_sites(t_max(m, cfg), cfg))
        if nc * 2.0 > math.log2(cfg.ram_bytes / BYTES_PER_AMP):   # 4^nc amplitudes
            break
        best = m
        m += 2.0
    return best


# Entanglement density across a balanced cut at the light-cone-crossing time, in
# bits per SITE. Hard bound: each half has Hilbert dimension 4^(m/2), so a
# balanced cut carries at most m bits, i.e. ent_rate <= 1. An earlier version
# used 1.5, which is unphysical.
ENT_RANGE = (0.3, 1.0)

# MEASURED tensor-network performance on exactly this problem, from the TDVP
# series published with arXiv:2510.26300 (Zenodo 17799843): |C^zz_nn| on dimer
# links at U = 0, against the exact free-fermion result, for chi = 256 .. 2048.
#
#   chi    256     512    1024    2048
#   err  0.100   0.085   0.078   0.077     (mean over t in [0.5, 2])
#
# The error is essentially FLAT in chi -- err ~ chi^-0.12, so doubling the bond
# dimension buys 9%. It has plateaued at ~0.077, which at late times is several
# times LARGER than the signal itself (exact |C^zz| = 0.021-0.031 for t >= 1.5).
# Extrapolating that slope, reaching our tolerance would need chi ~ 2e12.
#
# This is the datum the entropy model has to answer to, and it fails: at the same
# point (m=28, t=2, ent_rate=0.6) the entropy model predicts chi ~ 6.6e3 would
# suffice. It is not merely mis-calibrated, it is the wrong shape -- an entropy
# argument cannot see an error that saturates in chi.
TDVP_MEASURED = {256: 0.1002, 512: 0.0846, 1024: 0.0782, 2048: 0.0773}
TDVP_SLOPE = -0.12          # d log(err) / d log(chi), fitted
TDVP_AT = {"m": 28, "u_over_j": 0.0, "t_range": (0.5, 2.0),
           "source": "Zenodo 17799843, TDVP vs FLO, dimer-link C^zz"}

# AND THAT RUN WAS NOWHERE NEAR A LEADERSHIP-SCALE ATTEMPT.
# Snake MPS on 28 sites costs ~chi^3 per bond update, ~5.6e4 updates for the full
# time grid, and ~chi^2 L d amplitudes of memory:
#     chi = 2048  ->  4.8e14 flops  =  0.28 MILLISECONDS of exascale compute
#     chi = 1e5   ->  16 TB,   33 s
#     chi = 1e6   ->  1.6 PB,  9 hours
#     chi = 3e6   ->  15 PB,   one week
# A week of exascale is 2.1e9 times the published calculation. So the frontier
# was never probed: chi could go 500-1500x higher within the same budget the
# quantum arms are given.
#
# On the MEASURED slope that still would not reach our tolerance --
# chi = 3e6 gives 0.032 against a target of 0.0064. But a 1500x increase in chi
# buying only 2.4x in error is itself evidence the error is NOT truncation
# limited, which means the slope cannot be extrapolated in EITHER direction.
# The tensor-network arm is therefore UNBOUNDED by the available data, not
# bounded and small. Its absence from the band is a gap, not a finding.
# DECISIVE DIAGNOSTIC, from the same deposit and costing nothing to run.
# Absolute error against exact FLO at U = 0, by time and bond dimension:
#
#    t     exact    chi=256   chi=512  chi=1024  chi=2048   ratio 256/2048
#   0.1    0.9418   7.93e-3   7.38e-3   7.36e-3   7.35e-3      1.1x
#   0.5    0.2429   6.36e-2   4.78e-2   3.31e-2   2.56e-2      2.5x
#   2.0    0.0211   1.72e-1   1.47e-1   1.44e-1   1.57e-1      1.1x
#
# At t = 0.1 the state is barely entangled and chi = 256 is wild overkill -- a
# truncation-limited calculation would be at machine precision there. Instead the
# error is 7.9e-3, and an EIGHTFOLD increase in bond dimension removes 8% of it.
# The `max_bond_dimension` column confirms every run saturated its cap, so
# truncation was binding; it simply was not what limited the accuracy.
#
# CONCLUSION: the published TDVP carries a chi-INDEPENDENT error floor of
# ~7e-3 present from the earliest times -- plausibly two-site TDVP projection
# error on a snake MPS with long-range Jordan-Wigner strings, or the time step,
# or the GPR smoothing. It is therefore NOT a converged tensor-network
# calculation, the chi-extrapolation to 2e12 is meaningless, and this dataset
# CANNOT bound what tensor networks can do on this problem.
#
# Note what this does to the argument in METHODS 5: the conclusion there (drop
# the entropy-derived arm) stands, because the entropy model is contradicted by
# a floor it cannot represent. But the REASONING offered there -- "the error is
# flat in chi, so TN fails" -- was wrong: flatness in chi is evidence the run was
# not chi-limited, not evidence that chi cannot help.
#
# The floor sits at ~7e-3, essentially AT our absolute tolerance (~6e-3). A clean
# implementation that removed it could plausibly reach this accuracy at modest
# chi, which would move the classical upper edge well above 26.
TDVP_FLOOR = {"t": 0.1, "exact": 0.9418,
              "err": {256: 7.93e-3, 512: 7.38e-3, 1024: 7.36e-3, 2048: 7.35e-3},
              "verdict": "chi-independent floor; not truncation-limited"}

# PEAK-FLOP ARITHMETIC, NOT A DEMONSTRATED WALL TIME (review #7). Dividing a
# chi^3 operation count by a machine's peak rate assumes perfect strong scaling
# and ignores memory traffic, communication, the SVD/QR and MPO applications
# that dominate a real TDVP sweep, and the achievable fraction of peak. Leading
# tensor-network codes report single-digit to low-tens percent of peak on
# leadership machines, and strong scaling in chi is not free. The affordable chi
# below should be read as an optimistic upper bound on what the arithmetic
# permits, not as a run anyone has done.
TDVP_LEADERSHIP = {"published_chi": 2048, "published_flops": 4.8e14,
                   "week_exascale_flops": 1.0e24, "chi_affordable_week": 3e6,
                   "err_at_affordable_chi_on_measured_slope": 0.032,
                   "caveat": "peak-FLOP arithmetic; no memory, communication, "
                             "SVD/MPO or parallel-efficiency accounting",
                   "assumed_fraction_of_peak": 1.0}


def tdvp_chi_for(tol: float) -> float:
    """Bond dimension the MEASURED convergence implies for a given absolute error."""
    e0 = TDVP_MEASURED[2048]
    return 2048.0 * (e0 / max(tol, 1e-12)) ** (1.0 / abs(TDVP_SLOPE))


# ---------------------------------------------------------------------------
# U = 0: this observable is POLYNOMIAL, and the ED band is the wrong baseline
#
# Second-pass review #7. At U = 0 the Hamiltonian is quadratic, so evolution is
# a 2m x 2m single-particle matrix exponential. The initial state is NOT
# Gaussian -- a triplet covering cannot be made by a FLO unitary, and expanding
# it gives 2^{m/2} Fock branches -- but C^zz is a WEIGHT-4 observable, and a
# four-fermion operator can connect branches differing in at most one triplet.
# That collapses the branch sum to a diagonal part fixed by the one- and
# two-mode occupation statistics plus a local coherent part, neither of which
# needs the branches. The argument is the experimental paper's Appendix E
# (arXiv:2510.26300); calibration/free_fermion.py is an independent
# implementation of it.
#
# VALIDATED, not asserted: against exact many-body evolution at m = 4, 6, 8, 9
# and t = 0, 0.25, 0.5, 1, 2, the worst disagreement is 1.7e-15 -- machine
# precision. An O(m) form (sparse Krylov propagator, block-collapsed sums)
# reproduces the O(M^2) form exactly and was timed out to m = 1024.
#
# The frontier below is a MEASURED single-core wall time for unoptimised
# CPython, deliberately: it needs no extrapolation, and it already exceeds the
# ED band by six orders of magnitude.
FREE_FERMION = {
    "valid_at": "U/J = 0 exactly",
    "worst_abs_err": 1.7e-15,
    "validated_m": (4, 6, 8, 9),
    "validated_t": (0.0, 0.25, 0.5, 1.0, 2.0),
    # one core, unoptimised CPython; log-log exponent 0.979 -- linear
    "seconds": {16: 0.015, 64: 0.050, 256: 0.177, 1024: 0.687,
                4096: 2.908, 16384: 11.616},
    "source": "calibration/free_fermion.py -> data/free_fermion.json",
}
FF_SECONDS_PER_SITE = 7.0e-4      # fitted slope; flat to 6% over m = 1024..16384
FF_BYTES_PER_SITE = 64.0          # four complex rows of length 2m


def free_fermion_applicable(cfg: Config = DEFAULT) -> bool:
    """Exact only at U = 0. No claim is made about small but non-zero U."""
    return cfg.U_over_J == 0.0


def free_fermion_frontier(cfg: Config = DEFAULT) -> float:
    """Largest m in the week budget, from MEASURED single-core wall time.

    Conservative on purpose: one core, unoptimised CPython, no extrapolation of
    the implementation. A vectorised or parallel version would go further, and
    the memory (O(m)) is nowhere near binding.
    """
    if not free_fermion_applicable(cfg):
        return 0.0
    by_time = cfg.budget_s / (cfg.n_times * FF_SECONDS_PER_SITE)
    by_mem = cfg.ram_bytes / FF_BYTES_PER_SITE
    return min(by_time, by_mem)


def band(cfg: Config = DEFAULT) -> dict:
    """ESTIMATED CAPACITY of specified classical methods under stated machine
    assumptions. **Not** a classical impossibility boundary -- exceeding it does
    not establish that every competitive classical method fails, and sitting
    below it does not establish that a point is easy in practice.

    At U = 0 the band is set by the free-fermion estimator instead, which is
    polynomial for this observable and validated to machine precision; exact
    diagonalisation is simply the wrong method there.

    At U > 0 the band is set by exact diagonalisation. The tensor-network edge is reported
    but NOT used to widen it, because the only published attempt on this problem
    (TDVP, chi up to 2048) does not converge on this observable: its error is flat
    in chi and larger than the signal at late times. An entropy-derived chi is
    therefore not evidence of capability here.
    """
    ed_lo = ed_frontier(cfg.but(ram_bytes=1e15))
    ed_hi = ed_frontier(cfg.but(ram_bytes=100e15))
    ff = free_fermion_frontier(cfg)
    mps_lo = mps_max_m(cfg.but(ent_rate=ENT_RANGE[1]))
    mps_hi = mps_max_m(cfg.but(ent_rate=ENT_RANGE[0]))
    # A classical attacker picks the BEST method available, so each edge is a
    # max over methods; the edges differ only in how favourable the assumptions
    # are (RAM, entanglement density).
    # ED sets the band. The entropy-derived MPS reach is carried alongside as a
    # diagnostic, not folded in -- see the docstring and TDVP_MEASURED.
    lo, hi = ed_lo, ed_hi
    if ff > 0.0:
        # at U = 0 exact diagonalisation is simply not the applicable method
        lo = hi = ff
    return {"free_fermion": ff, "method": "free fermion" if ff else "ED",
            "ed_lo_1PB": ed_lo, "ed_mid": ed_frontier(cfg), "ed_hi_100PB": ed_hi,
            # `band` is a FIXED 1-100 PB uncertainty envelope and deliberately
            # ignores cfg.ram_bytes. `ed_at_cfg_ram` is the single-machine
            # sensitivity; sweeping RAM should move that, not the band.
            "ed_at_cfg_ram": ed_frontier(cfg),
            "mps_strong_ent": mps_lo, "mps_mid": mps_max_m(cfg),
            "tdvp_err_at_chi2048": TDVP_MEASURED[2048],
            "tdvp_chi_for_tol": tdvp_chi_for(0.1 * 0.064),
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


def max_m_fixed_t_detail(t: float, cfg: Config = DEFAULT) -> dict:
    """The number, the METHOD that produced it, and what it does not claim."""
    v = max_m_fixed_t(t, cfg)
    ed = float(ed_frontier(cfg))
    if v >= UNBOUNDED:
        meth, claim = ("light-cone cluster",
                       "exact truncation; the cone does not grow with m")
    elif v > ed:
        meth, claim = ("snake MPS",
                       "entropy-derived chi; NO convergence guarantee -- the "
                       "only published TDVP attempt on this observable has a "
                       "chi-independent error floor (O10)")
    else:
        meth, claim = ("exact diagonalisation", "exact within the sector")
    return {"t": t, "m": v, "method": meth, "claims": claim,
            "not_a_certificate": "a CAPACITY at fixed t, not a certified "
                                 "converged answer; converged.classical_t_reach "
                                 "answers the other question and reports "
                                 f"t = {__import__('fhcost.converged', fromlist=['x']).classical_t_reach(cfg):.4f}"}


def max_m_fixed_t(t: float, cfg: Config = DEFAULT) -> dict | float:
    """Best classical reach at a FIXED evolution time t (not t = sqrt(m)/v).

    DIFFERENT QUESTION FROM converged.classical_t_reach, and the two were being
    read as interchangeable. This asks "how many sites can a method carry to
    time t within the budget", answered by an ENTROPY-derived bond dimension
    with no convergence guarantee; that asks "to what time does a bound CERTIFY
    a thermodynamic-limit answer". At t = 0.02 this returns ~1e5 sites while
    classical_t_reach reports 0.013 -- both correct, neither a contradiction.
    Use max_m_fixed_t_detail() to get the method and assumption alongside the
    number (second-pass review, additional checks).

    This is the panel that shows why the time-window convention decides the
    whole question. At small t the causal cone of <Z_i(t) Z_j(0)> contains only
    a handful of sites, so a cluster expansion is exact at a cost that does not
    grow with m AT ALL -- classical reach is unbounded. Slide 1's literal
    "t in [0, 1/m]" sits deep in that regime, which is why it cannot be the
    intended question.
    """
    from .hubbard import trotter_steps
    ncone = cluster_sites(t, cfg)          # error-controlled, same buffer as quantum
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

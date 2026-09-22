"""Shared constants and the Config object every cost model takes.

VALUE PROVENANCE (all inputs, with source; ESTIMATE unless marked SPEC)
  SPEC (from Q1_Fermi_Hubbard-1.pdf, "HARDWARE 1"):
    p = 1e-3 two-qubit depolarizing; dt_gate = 10 ns; dt_meas = 100 ns;
    nearest-neighbour CNOT + Z(theta) + H on a sqrt(n) x sqrt(n) grid;
    eps = 1e-1 RELATIVE error on <Z_i(t) Z_j(0)>; one week of wall clock.
  budget_s = 6.048e5 s            SPEC ("if we have the hardware for a week").
  v = 4.0  (units J a / hbar)     ESTIMATE. The 2D tight-binding band
    eps(k) = -2J(cos kx + cos ky) has max group velocity |d eps/dk| = 2J along a
    lattice axis. A TWO-point correlator spreads on the two-particle light cone,
    i.e. at ~2 x that, so v = 4J. Knob: Config.v.
  t_max = sqrt(m)/v               USER DECISION. Evolve until the light cone
    crosses the lattice. This supersedes the literal "t in [0, 1/m]" printed on
    slide 1, which puts the light cone inside one lattice spacing and makes the
    observable classically trivial at any m.  See METHODS.md.
  c_g = 15 two-qubit gates / site / 2nd-order Trotter step   ESTIMATE.
    2m bonds x 2 spins hopping (~3-6 2q gates each as a Givens rotation) + m
    on-site ZZ (2 CNOT + 1 Rz), divided by m sites.
  c_rot = 5 arbitrary-angle rotations / site / step          ESTIMATE, same count.
  c_w = 1.0                       ESTIMATE. Absorbed O(1) prefactor in the
    2nd-order Trotter commutator sum W2 = c_w * m  (units J=1, U/J=4).
  f_emp = 5.0                     ESTIMATE. Worst-case Trotter bounds overestimate
    the true error by 1-2 orders for Hubbard [Childs, Su, Tran, Wiebe, Zhu,
    PRX 11, 011020 (2021)]; since r ~ sqrt(W), a 25x error overestimate is a 5x
    step-count overestimate. Only used when trotter="empirical".
  s_sig = 0.1                     ESTIMATE. Magnitude of the correlator being
    resolved. Relative error eps on a signal of size s costs (s*eps)^-2 shots.
  n_times = 20                    ESTIMATE for the "list of T times" on slide 1.
  p_th = 1e-2, p_L = 0.1 (p/p_th)^(d/2)   Fowler et al., PRA 86, 032324 (2012).
  t_round = 1e-6 s                Superconducting surface-code round.
  ross_selinger = 1.15 * log2(1/eps_rot) T gates per arbitrary Rz
                                  Ross & Selinger, arXiv:1403.2975.
  eps is RELATIVE: eps_abs = eps_rel * s_sig. The Trotter budget is absolute, so
    getting this wrong inflates r by sqrt(10) and costs 1.67x in m.
  v_B = 2J   ESTIMATE. Free-fermion max axial group velocity of
    eps(k) = -2J(cos kx + cos ky); this is the operator-spreading / butterfly
    speed that sets the causal cone. The correlation FRONT moves at ~2 v_B = 4J
    (quasiparticle-pair picture) and is used only for finite-size reasoning.
    The rigorous Lieb-Robinson bound (~2e J z) is ~10x looser and is NOT used.
  P_Z,1 = 2p/15 per injected STAR rotation, RUS total 4p/15, d-INDEPENDENT
                                  Akahoshi et al., PRX Quantum 5, 010337 (2024).
  Willow-calibrated p_L: 0.143%/round at d=7, suppression Lambda = 2.14 per
    Delta d = 2                   Google Quantum AI, Nature 638, 920 (2025).
  factory_tiles = 72 d^2 per 15-to-1 T factory (a 12d x 6d patch, one state per
    5.5d rounds)                  Litinski, Quantum 3, 128 (2019).
"""
from __future__ import annotations
from dataclasses import dataclass, replace

WEEK_S = 7 * 24 * 3600.0          # 6.048e5


@dataclass(frozen=True)
class Config:
    # --- hardware (SPEC) ---
    p: float = 1e-3               # two-qubit depolarizing rate
    dt_gate: float = 10e-9
    dt_meas: float = 100e-9
    budget_s: float = WEEK_S

    # --- problem (SPEC + user decision) ---
    eps: float = 0.1              # relative error target
    n_times: int = 20             # T time points
    U_over_J: float = 4.0
    v: float = 2.0                # BUTTERFLY velocity v_B, sets t_max and the causal cone
    v_corr: float = 4.0           # correlation-front speed (Calabrese-Cardy pair picture),
                                  # sets finite-size contamination only
    tmax_mode: str = "sqrt_m"     # "sqrt_m" | "const" | "inv_m"
    tmax_const: float = 1.0       # used when tmax_mode == "const"

    # --- circuit model (estimates) ---
    encoding: str = "compact"     # "compact" (Derby-Klassen) | "jw"
    c_g: float = 15.0             # 2q gates / site / step
    c_rot: float = 5.0            # rotations / site / step
    c_w: float = 1.0              # Trotter commutator prefactor
    trotter: str = "measured"     # "measured" | "extensive" | "lightcone" | "empirical"
                                  # "measured" uses W_MEASURED below: exact-diagonalisation
                                  # calibration of THIS observable and initial state on
                                  # patches of 4-12 sites, matching the collaborator's
                                  # geometry. See hubbard.w_commutator.
    f_emp: float = 5.0
    trotter_order_k: int = 1      # 1 = plain 2nd order; k>1 = 2k-order multiproduct

    # --- statistics ---
    s_sig: float = 0.1            # |<Z_i(t)Z_j(0)>_c| at the light-cone front.
                                  # THE largest un-pinned number in the model: it enters
                                  # twice (shots ~ 1/s^2, Trotter r ~ s^{-1/2}), net m ~ s^{2/9}.
    # WHICH GATES ACTUALLY DAMP THE OBSERVABLE
    #   "cone"    -- every gate in the causal cone, weighted by lightcone_frac.
    #                Geometric; what this model originally assumed.
    #   "support" -- only gates whose qubits overlap the Heisenberg-evolved
    #                observable's SUPPORT. Measured: on the Phasecraft/Quantinuum
    #                circuit (2415 two-qubit gates, their quoted p = 1e-3) the
    #                observed attenuation is Lambda = 0.20, against 2.58 from the
    #                cone model -- a 12.9x overcharge. The mechanism was checked
    #                by a weight test on their own data: Lambda(ZZ)/Lambda(Z) =
    #                1.91, where damping-proportional-to-operator-weight predicts
    #                2 and uniform-per-gate predicts 1.
    damping_model: str = "cone"   # "cone" | "support"  -- see OPEN_ITEMS.md O5
    w_obs0: float = 4.0           # qubits in the observable's support
                                  # (C^zz on nearest neighbours: 2 sites x 2 spins)
    support_growth: float = 0.0   # extra support per gate layer. FITTED to a SINGLE
                                  # circuit and consistent with zero there; the true
                                  # value cannot be zero at long times, when the
                                  # operator must eventually fill the lattice. This is
                                  # the least-constrained input in the model.
    lightcone_frac: float = 1.0 / 3.0   # "cone" model only
    depol_factor: float = 16.0 / 15.0   # ATTENUATION only: a non-identity Pauli is damped
                                  # by (1 - 16p/15) per 2-qubit depolarizing gate (7 of the
                                  # 15 non-identity Paulis commute, 8 anticommute). This is
                                  # a physical fact and must NOT move when the PEC costing
                                  # convention changes.
    n_ancilla: int = 1            # Hadamard-test ancilla for the two-time correlator

    # --- fault tolerance ---
    p_th: float = 1e-2
    t_round: float = 1e-6
    routing: float = 2.0          # 2 tiles per algorithmic qubit -> ~4d^2 (Litinski data block)
    pl_model: str = "fowler"      # "fowler" (idealised) | "willow" (measured); see ftqc.py
    n_factories: int = 8
    d_max: int = 101

    # --- STAR (least-pinned; see METHODS.md) ---
    pin_engines: int = 1          # magic engines. The paper's Hubbard calculation uses
                                  # exactly ONE for the whole machine, so T states arrive
                                  # at one per logical cycle regardless of lattice size.
                                  # Exposed because adding engines is the obvious scaling
                                  # lever, but >1 is OUR extrapolation, not theirs.
    pin_nonlocal: bool = True     # Pinnacle needs qLDPC connectivity, which the slide's
                                  # nearest-neighbour 2D grid does NOT provide. Keeping
                                  # this True means the arm is costed under a DIFFERENT
                                  # hardware assumption from every other curve.
    star_rounds_per_step: float = 20.0   # sequential logical layers per Trotter step, in
                                  # units of d. STAR's largest tax: a us*d code cycle
                                  # instead of a 10 ns gate clock. Least-pinned STAR input.
    star_kappa: float = 0.533     # Lambda_STAR = 0.533 * p * n_rot, from P_Z,1 = 2p/15 and
                                  # gamma^2 = exp(8 P_Z,1 N)  [Akahoshi et al. Eqs. 6-7, 15]

    # --- classical ---
    ram_bytes: float = 10e15      # 10 PB, Frontier-class in-memory
    disk_bytes: float = 700e15    # Orion-class parallel filesystem
    disk_bw: float = 10e12        # ~10 TB/s. Out-of-core ED is memory-feasible but
                                  # BANDWIDTH-infeasible: see classical.ed_frontier.
    n_krylov_vec: int = 4
    flops: float = 1.7e18         # exascale
    ent_rate: float = 0.6         # bits of entanglement per site at t_max (see classical.py)

    # ---- reconciliation knobs: every assumption that differs from
    # ---- an independent parallel resource model. See fhcost/presets.py.
    trotter_mode: str = "bound"   # "bound" = Campbell commutator bound (ours)
                                  # "fixed_density" = r = steps_per_tau * t (theirs)
    steps_per_tau: float = 4.0    # only used by "fixed_density"
    pec_model: str = "exact"      # "exact"  = gamma = (15+14p)/(15-16p) per gate, from the
                                  #            signed Pauli inverse of the depolarizing
                                  #            channel (verified against its Pauli transfer
                                  #            matrix); log(Gamma^2) = 2 G ln(gamma).
                                  # "linear" = pec_coeff * p * G, for comparison only.
    pec_coeff: float = 4.0        # only used by pec_model="linear". 2 ln(gamma)/p = 4.000268
                                  # at p=1e-3, so this is the correct linearisation.
                                  # ATTENUATION and CANCELLATION ONE-NORM ARE DIFFERENT
                                  # QUANTITIES: an earlier version used one number (32/15)
                                  # for both, which made the PEC cost 1.9x too cheap.
    noise_channels: float = 1.0   # multiplier on Lambda for 1q + idle + SPAM.
                                  # ours 1.0 (2q only); theirs 20.02/10.59 = 1.89.
    star_law: str = "count"       # "count" = 0.533*p*N_rot (ours)
                                  # "angle" = 4*alpha*p*Theta + 4*floor*R (theirs)
    star_alpha: float = 1.5       # angle law slope   [Toshio et al., PRX 15, 021057]
    star_floor: float = 1e-5      # residual per-rotation error that angle cannot remove
    magic_source: str = "litinski"  # "litinski" 15-to-1 | "cultivation" [Gidney 2409.17595]
    lanes: float = 1.0            # concurrent logical delivery lanes (theirs: 4)
    routing_power: float = 0.0    # EXTRA powers of L = sqrt(m) in the 2-qubit gate
                                  # count from routing/SWAP on a NN grid.
                                  # ours 0 (G ~ m*r, per-site estimate);
                                  # theirs 1 (G ~ m^1.5*r, from the compiled
                                  # Willow anchor G2 = 4372 (r/3) (L/6)^3).
    fs_speed: str = "butterfly"   # which speed carries BOUNDARY contamination inward:
                                  # "butterfly"  v_B   -- a boundary perturbation spreads
                                  #                       at the operator-spreading speed
                                  # "correlation" 2 v_B -- if you (over)charge the
                                  #                       quasiparticle-PAIR front
                                  # "lieb_robinson" -- the rigorous bound, ~10x v_B
    v_lr: float = 20.0            # rigorous Lieb-Robinson speed, ~2 e J z. Far looser
                                  # than anything observed; included only as an upper edge.
    xi: float = 1.0               # correlation length, lattice units. Sets how far past
                                  # the light cone L must reach for finite-size
                                  # convergence: error ~ exp(-(L - 2 v_corr t)/xi).
    observable: str = "czz_nn"    # "czz_nn"   = equal-time CONNECTED nearest-neighbour spin
                                  #              correlation C^zz_ij(t), the observable the
                                  #              Phasecraft/Quantinuum experiment reports.
                                  #              Measured in the occupation basis: NO ancilla.
                                  # "two_time" = <Z_i(t)Z_j(0)>, needs a Hadamard ancilla.
                                  #              NOTE: its CONNECTED part vanishes identically
                                  #              from a Z-eigenstate, so it is only meaningful
                                  #              with a non-eigenstate start.
                                  # "local"    = <n_up n_dn> doublon density, no ancilla
    # --- signal model for C^zz_nn(t), from arXiv:2510.26300 ---
    # The initial state is a dimer covering of S^z_tot = 0 TRIPLETS (not a Neel/Z
    # eigenstate -- which is what makes the connected correlator nonzero at all).
    # On a dimer link that triplet gives <S^z_i S^z_j> = -1/4 and <S^z_i> = 0, so
    # |C^zz| = 1 exactly at t = 0. The order then melts and leaves a small residual
    # antiferromagnetic correlation. The paper reports melting at t ~ 0.4-0.7, that
    # melting is SLOWER for larger U, and that the residual is LARGER for larger U.
    # The residual magnitude is not quoted numerically there, so it is parameterised.
    # Experiment-calibrated mitigation (TFLO + GPR), anchored to arXiv:2510.26300.
    # They reach ~0.005 absolute on C^zz with 160 shots x 20 time points = 3200
    # shots, against 1/delta^2 = 4e4 for naive shot noise on a bounded observable.
    # So the effective sampling overhead is ~0.08, i.e. BELOW one: GPR borrows
    # statistics across correlated time points and TFLO removes bias by training
    # on a classically simulable free-fermion point. There is no exponential cost
    # here at all -- but it is bias-limited, and the bias is characterised only out
    # to the Lambda they actually ran at.
    exp_cal_overhead: float = 0.08   # measured effective sampling overhead
    exp_cal_lambda_max: float = 0.20  # largest Lambda demonstrated. Beyond this the
                                  # method is UNCHARACTERISED, not known to fail.
    # ---- ERROR LEDGER -------------------------------------------------------
    # Every contribution is a share of the SAME absolute tolerance and the shares
    # must sum to at most 1. They previously did not: Trotter, synthesis and
    # logical each took 0.30 while sampling took 0.50 (NISQ) or the whole budget
    # (FT), i.e. 1.4x and 1.9x over. Synthesis and logical get the small shares
    # because both are only logarithmic in the resource, so buying them down is
    # cheap; Trotter and statistics get the rest.
    frac_trotter: float = 0.25
    frac_syn: float = 0.10
    frac_logical: float = 0.05
    frac_magic: float = 0.05      # magic-state infidelity, union-bounded over the T count.
                                  # Previously UNBUDGETED: the model consumed T states up
                                  # to 320x too noisy at large n and never checked.
    frac_mitig: float = 0.05      # residual mitigation bias (ZNE); unused by PEC
    frac_stat: float = 0.50       # statistical HALF-WIDTH at confidence conf_z
    conf_z: float = 1.9600        # per-time two-sided 95%
    simultaneous: bool = False    # True -> Bonferroni over n_times (z = 3.02 at T=20)
    # FITTED to the published data (Zenodo 17799843), dimer-link C^zz, TFLO+GPR
    # mitigated, t in [0.1,2], 20 points, both U. The decay is GAUSSIAN, not
    # exponential: a free stretch exponent fits beta = 2.15 +- 0.09 (U=0) and
    # 1.98 +- 0.05 (U=4), i.e. consistent with exactly 2, as a quench requires
    # (dC/dt = 0 at t = 0). Forcing beta = 1 drives s_short to 1.22-1.25, which is
    # 7-8 sigma above the exact triplet value of 1 -- a clean sign the form, not
    # the data, is wrong.
    signal_regime: str = "curve"  # "short" | "long" | "curve" | "fixed" (use s_sig)
    signal_beta: float = 2.0      # stretch exponent; 2 = Gaussian    FITTED
    s_short: float = 1.0          # EXACT from the triplet algebra; fit gives 1.000 +- 0.004
    s_res_min: float = 0.043      # residual at U = 0                 FITTED +- 0.006
    s_res_slope: float = 0.021    # residual growth per (U/J)/4       FITTED +- 0.003
    t_melt_base: float = 0.426    # melting time at U = 0, in 1/J     FITTED +- 0.007
    t_melt_slope: float = 0.021   # melting-time growth per (U/J)/4   FITTED +- 0.008

    def but(self, **kw) -> "Config":
        return replace(self, **kw)


DEFAULT = Config()

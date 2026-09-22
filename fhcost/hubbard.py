"""m -> circuit: qubits, Trotter steps, gate counts, depth, wall-clock per shot.

Shared by every quantum curve. See budget.py for the value-provenance ledger.

THE ONE EXPONENT THAT MATTERS
  t_max = sqrt(m)/v, and 2nd-order Trotter needs r ~ t^{3/2} sqrt(W2/eps_T) with
  W2 = c_w * m, so

      r(m)       ~ (m^{1/2}/v)^{3/2} * sqrt(c_w m / eps_T)  ~  m^{5/4}
      G_total(m) = c_g * m * r(m)                           ~  m^{9/4}

  i.e. alpha = 9/4. Everything downstream -- the NISQ saturation exponent
  1/alpha = 4/9, the FT qubit budget, the crossings -- rides on this.

NOTE ON THE LOCAL-OBSERVABLE REFINEMENT
  <Z_i(t) Z_j(0)> is local, so only the causal cone contributes to the Trotter
  error: W2_eff = c_w * min(m, N_cone) with N_cone = (2 v t)^2. At the
  light-cone-crossing time t = sqrt(m)/v we have N_cone = 4m > m, so the
  refinement buys NOTHING -- the cone IS the lattice. It only bites at fixed
  short t, which is why the t_max panel of the figure exists.
"""
from __future__ import annotations
import math
from .budget import Config, DEFAULT


# Measured second-order Trotter error constant, W_eff = err * r^2 / tau^3, from
# exact evolution of the dimerised-triplet state and the C^zz_nn observable on
# patches of 4-12 sites (the collaborator's cross5/square2/square3/rectangle2x3
# geometry plus 3x4). Medians over n >= 9 and r >= 8; n = 6 is excluded because
# the lattice clips the causal cone there.
#
#   * W_eff is constant in r to four digits, so the tau^3 / r^2 form is exact.
#   * W_eff is INDEPENDENT of m. Campbell's extensive bound W = 9.5 m is therefore
#     wrong in kind for a two-site observable, and its overestimate grows without
#     bound: 150x at m = 9, ~450x at the operating point, and larger still beyond.
#   * W_eff FALLS with tau, because the absolute error on an observable that melts
#     to a small residual cannot keep growing as tau^3.
W_MEASURED = {          # W_eff[U/J][tau], medians over n >= 9, r >= 8
    0.0: {0.25: 0.0682, 0.5: 0.0466, 1.0: 0.0083, 2.0: 0.0067},
    4.0: {0.25: 1.4500, 0.5: 1.2896, 1.0: 0.1899, 2.0: 0.1190},
    8.0: {0.25: 4.6275, 0.5: 1.6841, 1.0: 0.2368, 2.0: 0.0701},
}
W_TAU_MIN, W_TAU_MAX = 0.25, 2.0
W_U_MIN, W_U_MAX = 0.0, 8.0


def _loginterp(x, xs, ys):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for a, b, ya, yb in zip(xs, xs[1:], ys, ys[1:]):
        if a <= x <= b:
            f = (math.log(x) - math.log(a)) / (math.log(b) - math.log(a))
            return math.exp((1 - f) * math.log(ya) + f * math.log(yb))
    raise AssertionError


def w_measured(t: float, u: float = 4.0) -> float:
    """Measured W_eff at evolution time t and coupling U/J. Clamped outside range.

    U-dependence is strong and was invisible while the calibration ran at U = 4
    only: W_eff spans 68x from U = 0 to U = 8 at tau = 0.25. Most of that is the
    step from U = 0 (where the model is free and the two hopping colours nearly
    commute, so the on-site term carries all the Trotter error) to U = 4; beyond
    that it flattens, consistent with the dynamics slowing as 4J^2/U.

    The U = 8, tau = 2 entry sits BELOW U = 4 -- non-monotonic, and the same
    accidental-zero-crossing artefact that shows up in the raw error there. It is
    kept rather than smoothed, but it should not be read as physics.
    """
    taus = sorted(W_MEASURED[0.0])
    us = sorted(W_MEASURED)
    per_u = [_loginterp(t, taus, [W_MEASURED[uu][x] for x in taus]) for uu in us]
    if u <= us[0]:
        return per_u[0]
    if u >= us[-1]:
        return per_u[-1]
    for a, b, ya, yb in zip(us, us[1:], per_u, per_u[1:]):
        if a <= u <= b:
            f = (u - a) / (b - a)            # linear in U, log in W
            return math.exp((1 - f) * math.log(ya) + f * math.log(yb))
    raise AssertionError


def w_commutator(cfg: Config = DEFAULT) -> float:
    """Second-order Trotter commutator norm per site, W2 = w * m.

    Campbell's PLAQ bound for 2D Fermi-Hubbard is the tightest published closed
    form: w ~ 9.5/site at U/J = 4 and ~24.3/site at U/J = 8 [Campbell, QST 7,
    015007 (2021)]. The nested commutators are dominated by terms linear and
    quadratic in U, so interpolate in U/J through those two anchors.
    """
    u = cfg.U_over_J
    w = 9.5 + (24.3 - 9.5) * (u - 4.0) / 4.0
    return cfg.c_w * max(w, 1.0)


def t_max(m: float, cfg: Config = DEFAULT) -> float:
    """Longest evolution time required, in units hbar/J."""
    if cfg.tmax_mode == "sqrt_m":
        return math.sqrt(m) / cfg.v
    if cfg.tmax_mode == "const":
        return cfg.tmax_const
    if cfg.tmax_mode == "inv_m":
        return 1.0 / m
    raise ValueError(f"unknown tmax_mode {cfg.tmax_mode!r}")


def cone_sites(m: float, t: float, cfg: Config = DEFAULT) -> float:
    """HAMILTONIAN light cone: sites within v t of the observable, capped.

    This is an approximate physical cone, NOT the circuit's exact backward
    dependency cone (see circuit_cone_sites). Hamiltonian evolution also has
    tails outside it, so truncation at v t is not exact -- it is accurate only to
    the exponential tail, which is why the finite-size criterion in converged.py
    carries an explicit xi ln(1/eps) buffer and this does not.
    """
    return min(float(m), (2.0 * cfg.v * t + 1.0) ** 2)


def circuit_cone_sites(m: float, cfg: Config = DEFAULT) -> float:
    """Backward dependency cone of the COMPILED circuit, in sites.

    Different object from the Hamiltonian cone: it is set by gate connectivity
    and depth, not by a velocity. A second-order Trotter step applies four
    nearest-neighbour hopping colour classes, so the observable's support grows
    by roughly two sites per step. With r in the tens to hundreds this saturates
    at the whole lattice almost immediately -- which is precisely why a cone
    argument cannot explain the measured damping, and why the support model
    (O5) does. Kept as a diagnostic.
    """
    r = trotter_steps(m, cfg)
    return min(float(m), (2.0 * 2.0 * r + 1.0) ** 2)


def mean_cone_fraction(m: float, cfg: Config = DEFAULT) -> float:
    """Time-averaged cone fraction over the evolution, with the lattice cap.

    The cone grows as (2 v u + 1)^2 until it hits m, then stays there. Averaging
    that over u in [0, t_max] gives 2/3 + 1/(2 sqrt m) - 1/(6 m^(3/2)), which
    tends to 2/3 -- NOT the 1/3 of an uncapped cone, which this model used. The
    closed form is checked against numerical integration in selftest.
    """
    if cfg.tmax_mode != "sqrt_m":
        t = t_max(m, cfg)
        if t <= 0:
            return 1.0
        n = 256
        acc = sum(cone_sites(m, t * (i + 0.5) / n, cfg) for i in range(n)) / n
        return acc / max(m, 1e-12)
    rm = math.sqrt(max(m, 1.0))
    return 2.0 / 3.0 + 1.0 / (2.0 * rm) - 1.0 / (6.0 * rm ** 3)


def qubits_per_site_total(m: float, cfg: Config = DEFAULT) -> float:
    """Qubits in one copy of the register (used as the denominator for support)."""
    return qubits_per_copy(m, cfg)


def qubits_per_copy(m: float, cfg: Config = DEFAULT) -> float:
    """Physical/logical qubits to hold one copy of the lattice."""
    # + n_ancilla for the Hadamard test: <Z_i(t) Z_j(0)> is a TWO-TIME correlator
    # and cannot be read off a prepare-evolve-measure circuit.
    # only the two-time correlator needs a Hadamard-test ancilla; the experiment's
    # equal-time C^zz is read straight out of the occupation basis
    anc = cfg.n_ancilla if cfg.observable == "two_time" else 0
    if cfg.encoding == "compact":       # Derby-Klassen, ~1.5 qubits per mode
        return 3.0 * m + anc
    if cfg.encoding == "jw":            # 2 modes per site, no ancillas
        return 2.0 * m + anc
    raise ValueError(f"unknown encoding {cfg.encoding!r}")


def t_melt(cfg: Config = DEFAULT) -> float:
    """Melting time of the initial triplet order. Slower at larger U, but only
    just: the fit gives 0.426 -> 0.447 between U/J = 0 and 4, a 5% effect. The
    paper's much larger visual difference comes from the residual, not the
    timescale."""
    return cfg.t_melt_base + cfg.t_melt_slope * (cfg.U_over_J / 4.0)


def s_residual(cfg: Config = DEFAULT) -> float:
    """Late-time antiferromagnetic residual of C^zz_nn. Larger at larger U."""
    return cfg.s_res_min + cfg.s_res_slope * (cfg.U_over_J / 4.0)


def signal(t: float, cfg: Config = DEFAULT) -> float:
    """|C^zz_nn(t)|: 1 at t=0 (exact triplet value), decaying to the AF residual.

    This is the knob that makes short- and long-time runs genuinely different
    problems, and their relative difficulty moves with U/J: raising U both slows
    the melt (spin dynamics set by the superexchange 4J^2/U) and raises the
    residual, so long-time targets get EASIER in signal while the circuit gets
    harder. The two effects pull opposite ways.
    """
    if cfg.signal_regime == "fixed":
        return cfg.s_sig
    if cfg.signal_regime == "short":
        return cfg.s_short
    res = s_residual(cfg)
    if cfg.signal_regime == "long":
        return res
    # Gaussian kernel (beta = 2), not exponential: the quench has dC/dt = 0 at t = 0.
    return res + (cfg.s_short - res) * math.exp(
        -((t / max(t_melt(cfg), 1e-9)) ** cfg.signal_beta))


def signal_at(m: float, cfg: Config = DEFAULT) -> float:
    """The signal at the longest time this lattice is evolved to."""
    return max(signal(t_max(m, cfg), cfg), 1e-6)


def conf_z(cfg: Config = DEFAULT) -> float:
    """Confidence factor. Bonferroni across the time grid if `simultaneous`."""
    if not cfg.simultaneous:
        return cfg.conf_z
    from scipy.stats import norm            # only needed for the stronger claim
    return float(norm.ppf(1.0 - 0.05 / (2.0 * max(cfg.n_times, 1))))


def error_ledger(cfg: Config = DEFAULT) -> dict:
    """The shares, and a hard check that they fit inside the tolerance."""
    d = {"trotter": cfg.frac_trotter, "finite size": cfg.frac_finite,
         "synthesis": cfg.frac_syn, "logical": cfg.frac_logical,
         "magic states": cfg.frac_magic, "mitigation bias": cfg.frac_mitig,
         "statistical": cfg.frac_stat}
    d["TOTAL"] = sum(d.values())
    if d["TOTAL"] > 1.0 + 1e-9:
        raise ValueError(f"error budget over-allocated: {d['TOTAL']:.3f} > 1")
    return d


def eps_absolute(cfg: Config = DEFAULT, m: float | None = None) -> float:
    """eps is RELATIVE; the Trotter/synthesis budgets are absolute.

    An absolute floor is applied so the target stays finite where the signal
    passes through zero (review finding #10).
    """
    error_ledger(cfg)          # validate at the entry point, not only on request
    s = cfg.s_sig if m is None else signal_at(m, cfg)
    return cfg.eps * max(s, cfg.s_res_min)


def trotter_steps(m: float, cfg: Config = DEFAULT, eps_trot: float | None = None) -> float:
    """Second-order (or 2k-order multiproduct) Trotter step count."""
    if eps_trot is None:
        eps_trot = cfg.frac_trotter * eps_absolute(cfg, m)
    t = t_max(m, cfg)
    if t <= 0:
        return 1.0
    if cfg.trotter_mode == "fixed_density":
        # A step-count CONVENTION, not an error model, so it outranks cfg.trotter:
        # r = ceil(4*tau) regardless of how the error would otherwise be estimated.
        # Not calibrated to the accuracy target; treat it as a floor.
        return max(1.0, math.ceil(cfg.steps_per_tau * t))
    if cfg.trotter == "measured":
        # W_eff is a property of the OBSERVABLE, not the lattice, so no factor of m
        W2 = w_measured(t, cfg.U_over_J)
        k = max(1, int(cfg.trotter_order_k))
        r = (t ** 1.5 * math.sqrt(W2 / eps_trot) if k == 1
             else t ** (1.0 + 1.0 / (2 * k)) * (W2 / eps_trot) ** (1.0 / (2 * k)))
        return max(1.0, r)
    W2 = w_commutator(cfg) * (cone_sites(m, t, cfg) if cfg.trotter == "lightcone" else m)
    k = max(1, int(cfg.trotter_order_k))
    if k == 1:
        r = t ** 1.5 * math.sqrt(W2 / eps_trot)
    else:
        # order-2k multiproduct / Richardson-in-dt: r ~ t^{1+1/2k} (W/eps)^{1/2k}
        r = t ** (1.0 + 1.0 / (2 * k)) * (W2 / eps_trot) ** (1.0 / (2 * k))
    if cfg.trotter == "empirical":
        r /= cfg.f_emp
    return max(1.0, r)


def step_depth(m: float, cfg: Config = DEFAULT) -> float:
    """Two-qubit-gate layers per Trotter step."""
    base = 12.0                                  # 4 hopping colour classes x 2 spins + onsite
    if cfg.encoding == "jw":
        base += 2.0 * math.sqrt(m)               # Kivlichan fermionic swap network
    return base


def counts(m: float, cfg: Config = DEFAULT) -> dict:
    """Everything the downstream models need, in one dict."""
    r = trotter_steps(m, cfg)
    t = t_max(m, cfg)
    depth = r * step_depth(m, cfg)
    # routing_power adds powers of L = sqrt(m): a compiled NN-grid circuit needs
    # SWAP networks that a per-site gate estimate does not see.
    route = m ** (0.5 * cfg.routing_power)
    g_total = cfg.c_g * m * r * route
    # Fraction of the circuit's gates that actually damp the observable.
    if cfg.damping_model == "support":
        # only gates overlapping the Heisenberg-evolved operator's support
        frac = min(1.0, (cfg.w_obs0 + cfg.support_growth * depth)
                   / qubits_per_site_total(m, cfg))
    elif cfg.damping_model == "cone":
        # the TIME-AVERAGED cone fraction, not a constant: with the lattice cap
        # the average is ~2/3, not the 1/3 of an uncapped cone
        frac = mean_cone_fraction(m, cfg)
    else:
        raise ValueError(f"unknown damping_model {cfg.damping_model!r}")
    g_cone = cfg.c_g * m * r * route * frac
    return {
        "m": m,
        "t_max": t,
        "steps": r,
        "q_per_copy": qubits_per_copy(m, cfg),
        "g_total": g_total,
        "g_cone": g_cone,
        "n_rot": cfg.c_rot * m * r * route,
        # rotations inside the causal cone. STAR's noise, like NISQ's, only
        # matters where it can reach the observable, so this must carry the SAME
        # light-cone factor as g_cone -- applying it to one and not the other
        # penalises STAR by exactly 3x.
        "n_rot_cone": cfg.c_rot * m * r * route * frac,
        "damp_frac": frac,
        "depth": depth,
        "t_circuit": depth * cfg.dt_gate + cfg.dt_meas,
    }


def multiproduct_branches(k: int) -> list[tuple[int, float]]:
    """[(step-count multiplier, Lagrange weight)] for an order-2k multiproduct.

    CLASSICAL EXTRAPOLATION OF EXPECTATION VALUES: branch i is its own circuit,
    run at k_i times the base step count, so it has k_i times the gates, k_i
    times the depth, and -- under PEC -- an overhead exponential in k_i. Charging
    only ||c||_1^2 ignores all of that; the deepest branch dominates.

    The convergence order is VALIDATED against exact diagonalisation: measured
    3.5-4.1 for the order-4 formula and 5.7-6.5 for order-6, with error gains up
    to 5e6 over plain Trotter at the same base step count.
    """
    if k <= 1:
        return [(1, 1.0)]
    x = [(1.0 / i) ** 2 for i in range(1, k + 1)]
    out = []
    for i, xi in enumerate(x):
        c = 1.0
        for j, xj in enumerate(x):
            if i != j:
                c *= xj / (xj - xi)
        out.append((i + 1, c))
    return out


def multiproduct_l1(k: int) -> float:
    """1-norm of the multiproduct coefficients for an order-2k formula.

    Richardson extrapolation of 2nd-order Trotter on step counts (k_1..k_j) =
    (1,2,...,j) with j = k: the coefficients that cancel orders 2..2k-2 are the
    Lagrange weights at nodes h_i = 1/k_i, evaluated at h = 0. Their 1-norm is
    the shot-noise amplification, and it grows fast with k -- which is exactly
    the cost that must be charged against the depth saving.
    """
    if k <= 1:
        return 1.0
    nodes = [(1.0 / i) ** 2 for i in range(1, k + 1)]     # 2nd order -> error in h^2
    tot = 0.0
    for i, xi in enumerate(nodes):
        c = 1.0
        for j, xj in enumerate(nodes):
            if i != j:
                c *= xj / (xj - xi)
        tot += abs(c)
    return tot


if __name__ == "__main__":
    print(f"{'m':>6} {'t_max':>7} {'steps':>10} {'G_total':>11} {'depth':>10} "
          f"{'t_circ(s)':>11} {'q/copy':>8}")
    for m in (4, 16, 36, 64, 100, 256, 1024):
        c = counts(m)
        print(f"{m:>6} {c['t_max']:>7.2f} {c['steps']:>10.1f} {c['g_total']:>11.3g} "
              f"{c['depth']:>10.3g} {c['t_circuit']:>11.3g} {c['q_per_copy']:>8.0f}")
    import numpy as np
    ms = np.array([64.0, 256.0, 1024.0, 4096.0])
    gs = np.array([counts(x)["g_total"] for x in ms])
    print(f"\nalpha in G ~ m^alpha: {np.polyfit(np.log(ms), np.log(gs), 1)[0]:.4f}  (expect 2.25)")
    print("multiproduct 1-norm:", {k: round(multiproduct_l1(k), 2) for k in (1, 2, 3, 4)})

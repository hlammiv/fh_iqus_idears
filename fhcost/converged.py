"""Finite-size extrapolation: what happens when you only need a CONVERGED answer.

THE POINT
  Asking "what is the largest m?" assumes m is the figure of merit. For a LOCAL
  observable it is not. Finite-size error falls as exp(-(L - 2 v_corr t)/xi), so
  once the lattice is bigger than the correlation front has travelled, the finite
  answer IS the infinite answer, and running a larger lattice buys nothing:

      m_required(t) = (2 v_corr t + xi ln(1/eps_fs))^2

  Below that L the answer is a finite-lattice quantity; above it, it is the
  thermodynamic limit. The crossover is sharp -- exponentially so in L.

WHY THE CURRENT CONVENTION CANNOT SEE THIS
  With t_max = sqrt(m)/v_B we have 2 v_corr t_max = 4 sqrt(m) = 4L > L always.
  The front has hit the boundary by construction, so no lattice on the plot is
  ever converged and finite-size extrapolation buys exactly nothing. That is a
  property of the chosen time window, not of the method.

WHERE THE EXPONENTIAL TAKEOFF ACTUALLY IS
  In the converged regime the contest is no longer "how big a lattice" but "how
  far in time", and there the scaling separates violently:

      classical (light-cone / cluster):  cost ~ 4^{(2 v_B t + 1)^2}   -- exp(t^2)
      quantum:                           cost ~ poly(t) on m_required(t) sites

  The classical cost is exponential in the cone AREA. Every unit of t you add
  multiplies it by a large constant; the quantum cost grows polynomially. So the
  advantage takes off like exp(t^2) in TIME, at fixed (modest) lattice size --
  not like anything in m.
"""
from __future__ import annotations
import math
from .budget import Config, DEFAULT
from . import nisq, ftqc
from .hubbard import eps_absolute

BYTES_PER_AMP = 16.0


def fs_speed(cfg: Config = DEFAULT) -> float:
    """Speed at which the boundary contaminates the observable.

    WHAT THIS CRITERION ASSUMES, AND WHERE IT COMES FROM
      The form L >~ v t + xi ln(1/eps) is a ballistic front plus an exponential
      tail. Both pieces are state- and observable-dependent, and the usual
      quoted version is NOT generic:

      * The exponential tail outside the cone is a QUENCH statement -- it is
        clean for evolution from a low-entanglement product/Neel state, where
        Calabrese-Cardy quasiparticle pairs are emitted at t = 0. From a ground
        or thermal state the lattice already carries equilibrium correlations of
        range xi_eq BEFORE any evolution, and near half filling at low T the 2D
        Hubbard antiferromagnetic xi_eq grows -- that term then dominates and
        this criterion understates m_required badly.
      * The SPEED is the operator-spreading (butterfly) velocity v_B, not the
        correlation-front speed 2 v_B. A boundary perturbation is a single
        disturbance propagating inward; it does not need a quasiparticle pair.
        Charging 2 v_B (as an earlier version of this model did) makes
        m_required 4x too large.
      * The rigorous Lieb-Robinson velocity is ~10x v_B and would blow this up
        by ~100x in m. It is an upper edge, not an estimate.

      So this is a band, not a number, and it is conditional on a product-state
      quench. The initial state is not pinned by slide 1 -- see METHODS.md 0(b).
    """
    if cfg.fs_speed == "butterfly":
        return cfg.v
    if cfg.fs_speed == "correlation":
        return cfg.v_corr
    if cfg.fs_speed == "lieb_robinson":
        return cfg.v_lr
    raise ValueError(f"unknown fs_speed {cfg.fs_speed!r}")


def m_required(t: float, cfg: Config = DEFAULT) -> float:
    """Sites needed for a thermodynamic-limit answer at evolution time t.

    The boundary sits L/2 from the observable, so the condition is
    (L/2) > v t + xi ln(1/eps_fs), i.e. L > 2 v t + 2 xi ln(1/eps_fs).
    """
    eps_fs = cfg.frac_trotter * eps_absolute(cfg)
    L = 2.0 * (fs_speed(cfg) * t + cfg.xi * math.log(1.0 / max(eps_fs, 1e-12)))
    return L * L


def m_required_extrapolated(t: float, cfg: Config = DEFAULT, n_sizes: int = 4,
                            gain: float = 10.0) -> float:
    """m needed when you FIT the finite-size trend instead of waiting it out.

    Running n_sizes lattices and fitting exp(-(L-2 v_corr t)/xi) removes a factor
    `gain` of the residual, so the required L drops by xi*ln(gain) -- a constant
    shift, not a change of scaling. It is paid for with n_sizes x the shots.
    """
    eps_fs = cfg.frac_trotter * eps_absolute(cfg)
    L = 2.0 * cfg.v_corr * t + cfg.xi * math.log(1.0 / max(eps_fs * gain, 1e-12))
    return max(L, 1.0) ** 2


def classical_t_reach(cfg: Config = DEFAULT) -> float:
    """Largest t a light-cone/cluster expansion handles, at ANY lattice size.

    The cone holds N_c = (2 v_B t + 1)^2 sites at 4 states each, and its cost is
    independent of m. This is the number the quantum machine has to beat.
    """
    cap = math.log2(cfg.ram_bytes / BYTES_PER_AMP)
    lo, hi = 0.0, 50.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if 2.0 * (2.0 * cfg.v * mid + 1.0) ** 2 <= cap else (lo, mid)
    return lo


def quantum_t_reach(n: float, cfg: Config = DEFAULT, arm: str = "pec") -> float:
    """Largest t reachable in a week while ALSO holding m_required(t) sites."""
    def ok(t):
        if t <= 0:
            return True
        c = cfg.but(tmax_mode="const", tmax_const=t)
        need = m_required(t, cfg)
        if arm == "pinnacle":
            got = ftqc.max_m_pinnacle(n, c)
        elif arm == "surface":
            got = ftqc.max_m_surface(n, c.but(pl_model=cfg.pl_model))
        elif arm == "star":
            got = ftqc.max_m_star(n, c)
        else:
            got = nisq.max_m(n, c, arm)
        return got >= need
    if not ok(0.05):
        return 0.0
    lo, hi = 0.05, 60.0
    for _ in range(44):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def classical_cone_cost_log2(t: float, cfg: Config = DEFAULT) -> float:
    """log2 of the amplitudes a converged classical cluster calculation must hold."""
    return 2.0 * (2.0 * cfg.v * t + 1.0) ** 2


if __name__ == "__main__":
    cfg = DEFAULT
    tc = classical_t_reach(cfg)
    print("CONVERGED REGIME -- the contest is time, not lattice size\n")
    print(f"classical light-cone reach (any m): t = {tc:.2f}  "
          f"[cone {(2*cfg.v*tc+1)**2:.0f} sites, 10 PB]")
    print(f"m required at that t:               {m_required(tc, cfg):.0f} sites\n")
    print(f"{'t':>6} {'m_required':>11} {'classical cone cost':>21} {'vs 10 PB':>10}")
    for t in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
        c2 = classical_cone_cost_log2(t, cfg)
        print(f"{t:>6.2f} {m_required(t, cfg):>11.0f} {'2^%.0f amplitudes' % c2:>21} "
              f"{'fits' if c2 <= math.log2(cfg.ram_bytes/16) else '10^%.0f x over' % ((c2-math.log2(cfg.ram_bytes/16))*0.301):>10}")
    print(f"\n{'n':>9} {'NISQ+PEC':>10} {'STAR':>8} {'surface FT':>11}   (largest converged t)")
    for n in (1e4, 1e6, 1e8, 1e10):
        print(f"{n:>9.0e} {quantum_t_reach(n, cfg, 'pec'):>10.2f} "
              f"{quantum_t_reach(n, cfg, 'star'):>8.2f} "
              f"{quantum_t_reach(n, cfg, 'surface'):>11.2f}")

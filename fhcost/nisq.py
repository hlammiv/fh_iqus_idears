"""The mitigation ceiling: how big an m a NOISY, UNCORRECTED device reaches in a week.

THE LAW (re-implemented here; original: ~/Desktop/QC/logdet/nighthawk_review/
mitigation_ceiling.py, read-only source, not imported)
  A circuit with G two-qubit gates at rate p suffers Lambda = G*p expected
  faults and its signal is attenuated by exp(-Lambda). Beating that with
  statistics costs

      N_shots ~ Gamma(Lambda)^2 / (s*eps)^2 ,   Gamma >= exp(Lambda)

  where Gamma = exp(Lambda) is the probabilistic-error-cancellation sampling
  overhead, and Richardson ZNE is strictly worse because it pays the LARGEST
  noise-scaling node. Takagi, Endo, Benjamin & Mitarai, npj QI 8, 114 (2022)
  prove the exponential is unavoidable for ANY mitigation strategy, so
  Gamma = exp(Lambda) is an optimistic envelope, not an artifact of PEC.

WHY THIS SATURATES IN n
  More qubits do not make gates better. They buy PARALLEL COPIES of the circuit,
  P = floor(n / q_per_copy), hence P times more shots per week. But shots enter
  through a logarithm:

      Lambda_max = (1/2) ln( budget * P / t_circuit * (s*eps)^2 / T )

  so Lambda_max ~ (1/2) ln n, and since G ~ m^alpha with alpha = 9/4,

      m_max ~ (ln n / p)^{4/9}

  -- SLOWER than logarithmic. Ten decades of n buy ~1.5x in m. That is the
  quantitative content of the slide's "NISQ saturates" box, and it is stronger
  than the slide claims: it holds for every mitigation strategy at once.
"""
from __future__ import annotations
import math
from .budget import Config, DEFAULT
from .hubbard import counts, qubits_per_copy, multiproduct_l1

M_MIN = 4.0      # smallest real lattice is 2x2; m=1 has no hopping term


def richardson_nodes(k: int) -> list[float]:
    """Noise-amplification factors for order-k Richardson ZNE: 1, 3, 5, ..."""
    return [1.0 + 2.0 * i for i in range(k + 1)]


def richardson_coeffs(nodes: list[float]) -> list[float]:
    """Lagrange weights extrapolating to zero noise."""
    out = []
    for i, li in enumerate(nodes):
        c = 1.0
        for j, lj in enumerate(nodes):
            if i != j:
                c *= lj / (lj - li)
        out.append(c)
    return out


def overhead(lam: float, strategy: str = "pec") -> float:
    """Gamma(Lambda): the factor by which the shot count's sqrt is inflated."""
    if strategy == "none":
        return 1.0
    if strategy == "pec":
        return math.exp(min(lam, 700.0))
    if strategy.startswith("zne"):
        k = int(strategy[3:]) if len(strategy) > 3 else 2
        nodes = richardson_nodes(k)
        cs = richardson_coeffs(nodes)
        # optimal shot allocation across nodes -> (sum_i |c_i| e^{l_i Lambda})
        return sum(abs(c) * math.exp(min(l * lam, 700.0)) for c, l in zip(cs, nodes))
    raise ValueError(f"unknown strategy {strategy!r}")


def shots_required(m: float, cfg: Config = DEFAULT, strategy: str = "pec") -> float:
    c = counts(m, cfg)
    # Gamma^2 = exp(pec_coeff * noise_channels * p * G_cone); lam is half that
    # exponent so gamma = exp(lam) keeps the ZNE node algebra unchanged.
    lam = 0.5 * cfg.pec_coeff * cfg.noise_channels * cfg.p * c["g_cone"]
    if strategy == "none":
        # unmitigated: the attenuation bias itself must sit under eps
        if 1.0 - math.exp(-lam) > cfg.eps:
            return math.inf
    g = overhead(lam, strategy)
    if not math.isfinite(g):
        return math.inf
    # multiproduct (Richardson-in-dt) extrapolation cuts the DEPTH but amplifies
    # shot noise by the coefficient 1-norm. ||c||_1 = O(log k) -- polylogarithmic,
    # not exponential [Low, Kliuchnikov & Wiebe arXiv:1907.11679; Vazquez et al.,
    # Quantum 7, 1067 (2023)] -- so the depth saving wins easily.
    l1 = multiproduct_l1(cfg.trotter_order_k)
    return cfg.n_times * l1 * l1 * g * g / (cfg.s_sig * cfg.eps) ** 2


def shots_available(m: float, n: float, cfg: Config = DEFAULT) -> float:
    c = counts(m, cfg)
    # fractional copies: with >=1e5 shots the chip is time-shared, so rounding the
    # packing down is an artifact that puts false staircases on the curve
    copies = n / c["q_per_copy"]
    if copies < 1:
        return 0.0
    return cfg.budget_s * copies / c["t_circuit"]


def feasible(m: float, n: float, cfg: Config = DEFAULT, strategy: str = "pec") -> bool:
    return shots_required(m, cfg, strategy) <= shots_available(m, n, cfg)


def max_m(n: float, cfg: Config = DEFAULT, strategy: str = "pec",
          m_hi: float = 1e6) -> float:
    """Largest m satisfying the week budget. 0 if even m=1 fails."""
    if not feasible(M_MIN, n, cfg, strategy):
        return 0.0
    lo, hi = M_MIN, m_hi
    if feasible(hi, n, cfg, strategy):
        return hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if feasible(mid, n, cfg, strategy):
            lo = mid
        else:
            hi = mid
    return lo


def max_m_ideal(n: float, cfg: Config = DEFAULT) -> float:
    """p = 0 reference: qubit-limited only, plus the week clock on ~T/eps^2 shots."""
    m = n / (3.0 if cfg.encoding == "compact" else 2.0)
    if m < M_MIN:
        return 0.0
    need = cfg.n_times / (cfg.s_sig * cfg.eps) ** 2
    lo, hi = M_MIN, m
    if need * counts(hi, cfg)["t_circuit"] <= cfg.budget_s * max(1.0, math.floor(n / qubits_per_copy(hi, cfg))):
        return hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        ok = need * counts(mid, cfg)["t_circuit"] <= cfg.budget_s * max(
            1.0, math.floor(n / qubits_per_copy(mid, cfg)))
        lo, hi = (mid, hi) if ok else (lo, mid)
    return lo


if __name__ == "__main__":
    print("Gamma(Lambda) by strategy -- PEC is the optimistic envelope\n")
    print(f"{'Lambda':>7} {'none':>10} {'PEC':>12} {'ZNE k=1':>12} {'ZNE k=2':>12}")
    for lam in (0.1, 1, 3, 5, 10):
        row = [overhead(lam, s) for s in ("none", "pec", "zne1", "zne2")]
        print(f"{lam:>7.1f} " + " ".join(f"{x:>11.3g}" for x in row))

    print("\nresolvable m vs n (one week, p=1e-3)\n")
    print(f"{'n':>10} {'ideal p=0':>10} {'unmitig.':>10} {'ZNE k=2':>10} {'PEC':>10} {'Lambda_max':>11}")
    for n in (1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e9, 1e12):
        mi = max_m_ideal(n)
        m0 = max_m(n, strategy="none")
        mz = max_m(n, strategy="zne2")
        mp = max_m(n, strategy="pec")
        lam = DEFAULT.p * counts(max(mp, 1.0))["g_cone"]
        print(f"{n:>10.0e} {mi:>10.4g} {m0:>10.3g} {mz:>10.3g} {mp:>10.3g} {lam:>11.2f}")

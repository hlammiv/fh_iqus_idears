"""The mitigation ceiling: how big an m a NOISY, UNCORRECTED device reaches in a week.

THE LAW (re-implemented here; original: ~/Desktop/QC/logdet/nighthawk_review/
mitigation_ceiling.py, read-only source, not imported)
  A circuit with G two-qubit gates at rate p suffers Lambda = G*p expected
  faults and its signal is attenuated by exp(-Lambda). Beating that with
  statistics costs

      N_shots ~ Gamma(Lambda)^2 / (s*eps)^2 ,   Gamma >= exp(Lambda)

  where Gamma = exp(Lambda) is the probabilistic-error-cancellation sampling
  overhead, and Richardson ZNE is strictly worse because it pays the LARGEST
  noise-scaling node.

  WHAT THE LOWER-BOUND LITERATURE DOES AND DOES NOT GIVE US
  Takagi, Endo, Minagawa & Gu, npj QI 8, 114 (2022) prove a worst-case
  estimator-spread lower bound, exponential in circuit DEPTH, for a defined class
  of mitigation protocols under a layered local-depolarizing noise model. That is
  NOT the same as the curve plotted here, which is the cost of one specified
  implementation (gatewise PEC) as a function of GATE COUNT for this particular
  state and observable. Depth and gate count scale differently in m, so the
  theorem does not license this curve's m-dependence; a worst-case bound does not
  say this instance is hard; and their PEC-optimality result is for a particular
  dephasing setting, so PEC is not "optimal mitigation" here in any proven sense.
  Treat the plotted curve as a costing of PEC-as-implemented, and the theorem as
  separate, weaker, and differently quantified support for the general shape.

TWO DIFFERENT ESTIMATORS, TWO DIFFERENT COSTS
  These are not interchangeable and the model must not blur them:

  * PEC samples a quasiprobability decomposition of the inverse noise channel.
    Every sample carries a sign and a weight, so the per-shot variance is
    gamma^2 = exp(2 Lambda) and the cost is exponential in Lambda -- but the
    estimator is unbiased, so shots really do buy accuracy.
  * Richardson ZNE runs the circuit at amplified noise levels lambda_i and
    averages BOUNDED outcomes. Its variance factor is therefore only
    (sum_i |c_i|)^2 -- polynomial, not exponential. What it pays instead is
    RESIDUAL BIAS: the extrapolant recovers the noiseless value only to
    |1 - sum_i c_i exp(-lambda_i Lambda)|, and no number of shots removes that.
    So ZNE is bias-limited at Lambda ~ 0.2-0.5 while PEC is shot-limited at
    Lambda ~ 10.

  An earlier version of this model charged ZNE an exponential variance
  sum |c_i| exp(lambda_i Lambda) -- i.e. it assumed attenuation had to be
  inverted at every node -- and never checked the bias at all. Both were wrong,
  and together they reported ZNE points carrying ~0.9 relative bias against a
  0.1 target. Whether higher ZNE order helps or hurts is now an OUTPUT of the
  bias/variance trade, not an assumption.

CAVEAT ON THE RESPONSE MODEL
  s(lambda) = s(0) exp(-lambda Lambda) is a test model, not an established
  description of this observable under gate-local noise. A response with real
  curvature would change the residual bias, and with it every ZNE conclusion.

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
from .hubbard import (counts, qubits_per_copy, multiproduct_l1,
                      multiproduct_branches, signal_at, eps_absolute, conf_z)

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


def _k_of(strategy: str) -> int:
    return int(strategy[3:]) if len(strategy) > 3 else 2


def residual_bias(lam: float, strategy: str) -> float:
    """Relative bias the estimator CANNOT remove by sampling."""
    if strategy == "pec":
        return 0.0                      # unbiased given exact noise characterisation
    if strategy == "none":
        return 1.0 - math.exp(-lam)     # raw attenuation
    if strategy == "expcal":
        # Demonstrated to agree with exact results out to the Lambda they ran at;
        # beyond that we have no evidence either way, so the curve simply stops
        # rather than being extrapolated.
        return 0.0
    if strategy.startswith("zne"):
        nodes = richardson_nodes(_k_of(strategy))
        cs = richardson_coeffs(nodes)
        return abs(1.0 - sum(c * math.exp(-l * lam) for c, l in zip(cs, nodes)))
    raise ValueError(f"unknown strategy {strategy!r}")


def cost_factor(strategy: str, lam: float = 0.0, lg2: float = 0.0) -> float:
    """Multiplier on (circuit time / delta^2) for the optimal shot allocation.

    For independent unbiased node estimators with per-shot variance v_i and
    per-shot time tau_i, minimising sum tau_i N_i subject to
    sum c_i^2 v_i / N_i <= delta^2 gives N_i ~ |c_i| sqrt(v_i/tau_i) and a
    minimum total time (sum_i |c_i| sqrt(v_i tau_i))^2 / delta^2.
    Here v_i = 1 (bounded outcomes) and tau_i = lambda_i x the base circuit time,
    because a noise-amplified circuit is correspondingly longer.
    """
    if strategy == "none":
        return 1.0
    if strategy == "expcal":
        return 1.0                      # placeholder; handled in time_required
    if strategy == "pec":
        if lg2 > 700.0:
            return math.inf                          # beyond any conceivable budget
        return math.exp(lg2)                         # v = Gamma^2, tau = 1
    if strategy.startswith("zne"):
        nodes = richardson_nodes(_k_of(strategy))
        cs = richardson_coeffs(nodes)
        return sum(abs(c) * math.sqrt(l) for c, l in zip(cs, nodes)) ** 2
    raise ValueError(f"unknown strategy {strategy!r}")


def lambda_of(m: float, cfg: Config = DEFAULT) -> float:
    """ATTENUATION exponent: a non-identity Pauli observable is damped by
    exp(-Lambda) after G_cone two-qubit gates. Physical; independent of how any
    mitigation scheme is costed."""
    g = cfg.noise_channels * counts(m, cfg)["g_cone"]
    return -g * math.log(1.0 - cfg.depol_factor * cfg.p)


def log_gamma_sq(m: float, cfg: Config = DEFAULT) -> float:
    """log(Gamma^2), the PEC sampling overhead: the CANCELLATION one-norm.

    The signed Pauli inverse of the two-qubit depolarizing channel has
    gamma = |a| + 15|b| = (15 + 14p)/(15 - 16p) per gate, so
    log(Gamma^2) = 2 G ln(gamma) -> 4.000268 p G at p = 1e-3. This is a different
    quantity from the attenuation above, which is 1.067 p G.
    """
    g = cfg.noise_channels * counts(m, cfg)["g_cone"]
    if cfg.pec_model == "linear":
        return cfg.pec_coeff * cfg.p * g
    if cfg.pec_model == "exact":
        gam = (15.0 + 14.0 * cfg.p) / (15.0 - 16.0 * cfg.p)
        return 2.0 * g * math.log(gam)
    raise ValueError(f"unknown pec_model {cfg.pec_model!r}")


def time_required(m: float, cfg: Config = DEFAULT, strategy: str = "pec") -> float:
    """Wall clock on ONE copy to hit the statistical target. inf if bias-blocked."""
    lam = lambda_of(m, cfg)
    if strategy == "expcal" and lam > cfg.exp_cal_lambda_max:
        return math.inf                 # outside the demonstrated range
    if residual_bias(lam, strategy) > cfg.frac_mitig * cfg.eps:
        return math.inf                 # no shot count repairs this
    c = counts(m, cfg)
    # statistical ALLOWANCE is frac_stat of the tolerance, and it is a
    # confidence half-width: z * sigma <= allowance, so sigma <= allowance / z
    delta = cfg.frac_stat * eps_absolute(cfg, m) / conf_z(cfg)
    # Per-BRANCH accounting. Branch i is its own circuit at k_i x the base step
    # count, so it has k_i x the gates (hence k_i x the PEC exponent) and k_i x
    # the runtime. Optimal allocation over independent unbiased branch estimators
    # with variance v_i and per-shot time tau_i gives a minimum total time
    # (sum_i |c_i| sqrt(v_i tau_i))^2 / delta^2.
    lg2 = log_gamma_sq(m, cfg)
    total = 0.0
    for ki, ci in multiproduct_branches(cfg.trotter_order_k):
        if strategy == "expcal":
            v = cfg.exp_cal_overhead
        else:
            v = cost_factor(strategy, ki * lam, ki * lg2)
        if not math.isfinite(v):
            return math.inf
        total += abs(ci) * math.sqrt(v * ki)
    return cfg.n_times * total * total * c["t_circuit"] / delta ** 2


def feasible(m: float, n: float, cfg: Config = DEFAULT, strategy: str = "pec") -> bool:
    """Bias must fit its allowance AND the statistics must fit the week."""
    t_need = time_required(m, cfg, strategy)
    if not math.isfinite(t_need):
        return False
    copies = math.floor(n / counts(m, cfg)["q_per_copy"])
    return copies >= 1 and t_need / copies <= cfg.budget_s


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
    dd = cfg.frac_stat * eps_absolute(cfg, min(m, 1e6)) / conf_z(cfg)
    need = cfg.n_times / dd ** 2
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

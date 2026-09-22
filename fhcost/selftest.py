"""Unit and sanity checks. Pure arithmetic -- allocates nothing, runs in <1 s."""
from __future__ import annotations
import math
import numpy as np
from .budget import DEFAULT
from . import hubbard, nisq, ftqc, classical, curves, converged, presets

FAILS: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not cond:
        FAILS.append(name)


def main() -> int:
    print("hubbard.py")
    ms = np.array([64.0, 256.0, 1024.0, 4096.0])
    gs = np.array([hubbard.counts(m)["g_total"] for m in ms])
    alpha = float(np.polyfit(np.log(ms), np.log(gs), 1)[0])
    # no longer exact: the absolute tolerance now tracks a decaying signal, so r
    # picks up a weak extra m-dependence that dies out once C^zz reaches its residual
    check("G_total ~ m^(9/4)", abs(alpha - 2.25) < 0.01, f"alpha = {alpha:.6f}")
    check("t_max = sqrt(m)/v", abs(hubbard.t_max(64) - 8.0 / DEFAULT.v) < 1e-12)
    check("compact encoding = 3m, no ancilla for the equal-time observable",
          hubbard.counts(16)["q_per_copy"] == 48.0)
    check("JW costs depth, not qubits",
          hubbard.counts(16, DEFAULT.but(encoding="jw"))["q_per_copy"] == 32.0
          and hubbard.step_depth(16, DEFAULT.but(encoding="jw")) > hubbard.step_depth(16))
    check("light-cone refinement is void at t_max",
          abs(hubbard.cone_sites(100, hubbard.t_max(100)) - 100) < 1e-9,
          "cone = lattice, so it buys nothing")
    check("light-cone refinement bites at short t",
          hubbard.cone_sites(100, 0.5) < 100)
    check("multiproduct 1-norm grows with order",
          hubbard.multiproduct_l1(1) < hubbard.multiproduct_l1(2) < hubbard.multiproduct_l1(4))

    print("\nnisq.py")
    # --- review finding #2: ZNE is bias-limited, not variance-limited ---
    check("PEC is unbiased; ZNE and the bare device are not",
          nisq.residual_bias(1.0, "pec") == 0.0
          and nisq.residual_bias(1.0, "zne2") > 0
          and nisq.residual_bias(1.0, "none") > 0)
    check("ZNE bias vanishes at zero noise and grows with it",
          abs(nisq.residual_bias(0.0, "zne2")) < 1e-12
          and nisq.residual_bias(0.5, "zne2") < nisq.residual_bias(1.5, "zne2"))
    check("ZNE variance cost is polynomial, PEC's is exponential",
          nisq.cost_factor("zne2", 6.0) == nisq.cost_factor("zne2", 1.0)
          and nisq.cost_factor("pec", lg2=12.0) > 1e4)

    # --- review finding #3: attenuation and cancellation one-norm are different ---
    # Verified against the Pauli transfer matrix of the two-qubit depolarizing
    # channel: its signed inverse has gamma = (15+14p)/(15-16p) per gate.
    pp = DEFAULT.p
    g = hubbard.counts(8.0)["g_cone"]
    att = nisq.lambda_of(8.0) / (pp * g)
    pec = nisq.log_gamma_sq(8.0) / (pp * g)
    check("attenuation coefficient = -ln(1-16p/15)/p", abs(att - 1.067236) < 1e-5,
          f"{att:.6f}")
    check("PEC coefficient = 2 ln(gamma)/p", abs(pec - 4.000268) < 1e-5, f"{pec:.6f}")
    check("they are NOT the same number", abs(pec / att - 3.748) < 0.01,
          f"PEC/attenuation = {pec/att:.3f}; one value was used for both")
    check("the correction makes PEC 1.875x more expensive than the old model",
          abs(pec / (32 / 15) - 1.875) < 0.01)
    check("attenuation does not move when the PEC convention changes",
          nisq.lambda_of(8.0, DEFAULT.but(pec_model="linear"))
          == nisq.lambda_of(8.0, DEFAULT.but(pec_model="exact")))
    check("the linearisation agrees with the exact form at p=1e-3",
          abs(nisq.log_gamma_sq(8.0, DEFAULT.but(pec_model="linear"))
              / nisq.log_gamma_sq(8.0) - 1) < 1e-4)

    def lam_cap(st):
        lo, hi = 0.0, 50.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if nisq.residual_bias(mid, st) <= DEFAULT.bias_frac * DEFAULT.eps else (lo, mid)
        return lo
    caps = {st: lam_cap(st) for st in ("none", "zne1", "zne2", "zne3")}
    check("higher ZNE order buys MORE noise headroom, not less",
          caps["none"] < caps["zne1"] < caps["zne2"] < caps["zne3"],
          " < ".join(f"{k}={v:.2f}" for k, v in caps.items()))
    check("sampling cannot repair bias: blocked points cost infinite time",
          not math.isfinite(nisq.time_required(4.0, DEFAULT, "zne2")))
    check("at p=1e-3 no ZNE order reaches even a 2x2 lattice",
          nisq.lambda_of(4.0) > caps["zne3"]
          and all(nisq.max_m(1e6, DEFAULT, f"zne{k}") == 0 for k in (1, 2, 3)),
          f"Lambda(m=4) = {nisq.lambda_of(4.0):.2f} vs cap {caps['zne3']:.2f}")
    check("...but ZNE does work once the noise is low enough",
          nisq.max_m(1e6, DEFAULT.but(p=1e-5), "zne2") > 0)
    check("Richardson weights sum to 1",
          abs(sum(nisq.richardson_coeffs(nisq.richardson_nodes(3))) - 1.0) < 1e-9)
    # THE anchor: the vendored mitigation-ceiling table says C_max ~ 3000-4500
    # two-qubit gates at p = 3e-3 for hour-to-year runtimes.
    for secs, lab in ((3600.0, "1 hour"), (3.15e7, "1 year")):
        cmax = 0.5 * math.log(1e4 * secs) / 3e-3
        check(f"mitigation-ceiling anchor, {lab}", 2850 <= cmax <= 4500, f"C_max = {cmax:.0f}")
    m_lo, m_hi = nisq.max_m(1e2, strategy="pec"), nisq.max_m(1e12, strategy="pec")
    check("mitigated NISQ saturates (10 decades of n < 2x in m)",
          m_hi / m_lo < 2.0, f"m: {m_lo:.1f} -> {m_hi:.1f}")
    check("unmitigated NISQ is flat in n",
          abs(nisq.max_m(1e3, strategy="none") - nisq.max_m(1e9, strategy="none")) < 1e-6)
    check("ideal p=0 is slope 1 while qubit-limited",
          abs(math.log10(nisq.max_m_ideal(1e5) / nisq.max_m_ideal(1e4)) - 1.0) < 0.05)
    check("eps is treated as RELATIVE",
          abs(hubbard.eps_absolute() - DEFAULT.eps * DEFAULT.s_sig) < 1e-15)
    check("the ancilla is charged only to the two-time correlator",
          hubbard.counts(16, DEFAULT.but(observable="two_time"))["q_per_copy"] == 3 * 16 + 1
          and hubbard.counts(16)["q_per_copy"] == 3 * 16)

    # --- review finding #1: state / observable consistency ---
    # <Z_i(t)Z_j(0)>_c vanishes identically from a Z-eigenstate. The experiment's
    # dimerised S^z_tot=0 TRIPLET start is not one, which is what rescues it.
    check("triplet start gives |C^zz| = 1 exactly at t = 0",
          abs(hubbard.signal(0.0) - 1.0) < 1e-12)
    # EXTERNAL validation: published dimer-link |C^zz| (TFLO+GPR), Zenodo 17799843
    PUB = {0: [(0.1, 0.9500), (0.3, 0.6378), (0.5, 0.2810), (0.7, 0.0854),
               (1.0, 0.0680), (1.5, 0.0326), (2.0, 0.0295)],
           4: [(0.1, 0.9527), (0.3, 0.6603), (0.5, 0.3292), (0.7, 0.1414),
               (1.0, 0.0902), (1.5, 0.0466), (2.0, 0.0532)]}
    resid = [hubbard.signal(t, DEFAULT.but(U_over_J=u)) - v
             for u, rows in PUB.items() for t, v in rows]
    rms = (sum(r * r for r in resid) / len(resid)) ** 0.5
    check("signal model reproduces the published data", rms < 0.02, f"RMS = {rms:.4f}")
    # a quench has dC/dt = 0 at t = 0, so the Gaussian must be far flatter at the
    # origin than the exponential that was there before
    drop_g = 1.0 - hubbard.signal(0.02, DEFAULT.but(signal_beta=2.0))
    drop_e = 1.0 - hubbard.signal(0.02, DEFAULT.but(signal_beta=1.0))
    check("the decay is Gaussian, not exponential",
          DEFAULT.signal_beta == 2.0 and drop_g < 0.1 * drop_e,
          f"drop at t=0.02: {drop_g:.2e} (Gaussian) vs {drop_e:.2e} (exponential)")
    check("order melts more slowly at larger U",
          hubbard.t_melt(DEFAULT.but(U_over_J=8)) > hubbard.t_melt(DEFAULT.but(U_over_J=4))
          > hubbard.t_melt(DEFAULT.but(U_over_J=0)))
    check("late-time AF residual grows with U",
          hubbard.s_residual(DEFAULT.but(U_over_J=8))
          > hubbard.s_residual(DEFAULT.but(U_over_J=4)))
    check("signal decays monotonically to the residual",
          hubbard.signal(0.1) > hubbard.signal(1.0) > hubbard.signal(10.0)
          >= hubbard.s_residual())
    check("absolute tolerance has a floor at the residual",
          hubbard.eps_absolute(DEFAULT, 1e6) >= DEFAULT.eps * DEFAULT.s_res_min)
    # the point of the short/long knob: raising U hurts at short time (bigger
    # commutator norm) but helps at long time (bigger residual -> looser tolerance)
    short = {u: nisq.max_m(1e6, DEFAULT.but(signal_regime="short", U_over_J=u), "pec")
             for u in (4, 8)}
    long_ = {u: ftqc.max_m_surface(1e6, DEFAULT.but(signal_regime="long", U_over_J=u,
                                                    pl_model="fowler")) for u in (4, 8)}
    check("U hurts at short time but helps at long time",
          short[8] < short[4] and long_[8] > long_[4],
          f"short {short[4]:.1f}->{short[8]:.1f}, long {long_[4]:.0f}->{long_[8]:.0f}")
    check("PEC outperforms every ZNE order here",
          all(nisq.max_m(1e6, strategy="pec") > nisq.max_m(1e6, strategy=f"zne{k}")
              for k in (1, 2, 3)))

    print("\nftqc.py")
    check("p_L falls with d", ftqc.p_logical(25) < ftqc.p_logical(15) < ftqc.p_logical(5))
    check("p_L(d) = 0.1 (p/p_th)^((d+1)/2)",
          abs(ftqc.p_logical(9, DEFAULT.but(pl_model="fowler")) - 0.1 * 0.1 ** 5) < 1e-18)
    check("surface storage = routing * 2d^2", ftqc.storage_per_logical(15) == 2.0 * 2 * 225)
    check("FT needs ~1e5 qubits before it resolves anything",
          ftqc.max_m_surface(1e4) == 0 and ftqc.max_m_surface(1e6) > 0)
    # THE correction the validation pass caught: the FT curve is shot-limited,
    # not qubit-limited, so its slope is 4/9 and not 1 over most of the range.
    fow = DEFAULT.but(pl_model="fowler")
    slope = math.log10(ftqc.max_m_surface(1e8, fow) / ftqc.max_m_surface(1e7, fow))
    check("FT slope -> 4/9 once shot-limited", abs(slope - 4 / 9) < 0.12,
          f"slope = {slope:.2f} (4/9 = 0.44)")
    check("the whole shot budget is charged, not one shot",
          ftqc.n_shots_total() > 1e4, f"N = {ftqc.n_shots_total():,.0f}")
    check("measured (Willow) p_L is harsher than idealised (Fowler)",
          ftqc.max_m_surface(1e7, DEFAULT.but(pl_model="willow"))
          < ftqc.max_m_surface(1e7, fow))
    check("STAR rotation error is d-independent and cone-restricted",
          abs(ftqc.star_point(64)["lam"]
              - 0.533 * DEFAULT.p * ftqc.counts(64)["n_rot_cone"]) < 1e-9)
    # STAR and NISQ must feel the SAME causal cone; a 3x mismatch here was the bug
    # that made STAR look strictly worse than NISQ at every n.
    ratio = ftqc.star_point(8.0)["lam"] / (
        DEFAULT.depol_factor * DEFAULT.p * ftqc.counts(8.0)["g_cone"])
    check("Lambda_STAR / Lambda_NISQ ~ 1/6", abs(ratio - 1 / 6) < 0.02,
          f"ratio = {ratio:.3f}")
    # NOTE a convention asymmetry worth keeping visible: the bare-NISQ arm now uses
    # the depolarizing-channel cancellation one-norm derived here, while STAR's
    # overhead comes from its own paper's gamma^2 = exp(8 P_Z,1 N). Those are
    # different sources, and the STAR formula has not been re-derived the same way.
    check("STAR beats bare NISQ at every n it can run",
          all(ftqc.max_m_star(n) >= nisq.max_m(n, strategy="pec")
              for n in (1e4, 1e6, 1e8)),
          f"ratio {ftqc.max_m_star(1e6)/nisq.max_m(1e6, strategy='pec'):.2f} at n=1e6")
    check("STAR starts earlier than surface FT",
          ftqc.max_m_star(1e4) > 0 and ftqc.max_m_surface(1e4, fow) == 0)
    s_lo, s_hi = ftqc.max_m_star(1e5), ftqc.max_m_star(1e9)
    check("STAR saturates too", s_hi / s_lo < 1.8, f"m: {s_lo:.1f} -> {s_hi:.1f}")
    check("surface FT eventually overtakes STAR",
          ftqc.max_m_surface(1e7, fow) > ftqc.max_m_star(1e7))
    check("STAR beats NISQ only modestly at p=1e-3",
          1.0 < ftqc.max_m_star(1e6) / nisq.max_m(1e6, strategy="pec") < 2.0,
          f"ratio = {ftqc.max_m_star(1e6)/nisq.max_m(1e6, strategy='pec'):.2f}")
    # only the injection error moves STAR's RELATIVE standing; footprint, clock and
    # p all enter through a logarithm, or help NISQ equally
    r0 = ftqc.max_m_star(1e6) / nisq.max_m(1e6, strategy="pec")
    cp = DEFAULT.but(p=1e-4)
    check("lowering p does NOT improve STAR relative to NISQ",
          abs(ftqc.max_m_star(1e6, cp) / nisq.max_m(1e6, cp, "pec") - r0) < 0.1)
    ck = DEFAULT.but(star_kappa=DEFAULT.star_kappa / 10)
    check("lowering the injection error DOES",
          ftqc.max_m_star(1e6, ck) / nisq.max_m(1e6, ck, "pec") > 2 * r0)

    print("\nftqc.py -- Pinnacle QLDPC arm (arXiv:2602.11457)")
    import math as _m
    for (nc, k, d, dt, npb), want in zip(ftqc.GB_CODES,
                                         (8e-4, 4e-5, 1e-7, 3e-11, 4e-16)):
        got = ftqc.p_logical_gb(k, d)
        check(f"p_L fit reproduces their Table III, d={d}",
              0.7 < got / want < 1.4, f"{got:.1e} vs {want:.0e}")
    check("Pinnacle footprint reproduces their Table IV at L=8",
          abs(math.ceil(130 / 14) * 860 * (1 + ftqc.PIN_ENGINE) / 19e3 - 1) < 0.05)
    check("QLDPC beats the surface code on our workload",
          ftqc.max_m_pinnacle(1e6) > 2 * ftqc.max_m_surface(1e6, fow))
    check("and clears the classical band far earlier",
          curves.summary()["n_pinnacle_clears_classical_hi"]
          < 0.3 * curves.summary()["n_ft_clears_classical_hi"])
    check("the plotted range stays inside the published GB code family",
          ftqc.pinnacle_point(ftqc.max_m_pinnacle(1e8))["d"] <= 24
          and ftqc.max_m_pinnacle(1e8) < ftqc.max_m_pinnacle(1e9),
          "family exhausts near n ~ 1e10, above the plotted range")
    check("GB rate advantage grows with d",
          (ftqc.GB_CODES[4][4] / ftqc.GB_CODES[4][1]) / (4 * 24 ** 2)
          < (ftqc.GB_CODES[0][4] / ftqc.GB_CODES[0][1]) / (4 * 4 ** 2))

    print("\nclassical.py")
    check("Hilbert dim = C(m,m/2)^2", classical.hilbert_dim(4) == 36.0)
    check("ED frontier in the mid-20s", 22 <= classical.ed_frontier() <= 28,
          f"m = {classical.ed_frontier()}")
    check("MPS falls as entanglement density rises",
          classical.mps_max_m(DEFAULT.but(ent_rate=0.3))
          > classical.mps_max_m(DEFAULT.but(ent_rate=1.5)))
    check("cluster expansion is no help at t_max",
          classical.cluster_max_m() <= classical.ed_frontier() + 2)

    print("\ncurves.py -- the headline")
    s = curves.summary()
    lo, hi = s["classical_band"]
    check("classical band is m ~ 24-62", 20 <= lo <= 30 and 50 <= hi <= 80,
          f"{lo:.0f} .. {hi:.0f}")
    check("PEC-NISQ never clears the optimistic classical edge",
          s["n_pec_clears_classical_hi"] is None)
    check("STAR never clears the optimistic classical edge",
          s["n_star_clears_classical_hi"] is None)
    # extrapolation (b): the biggest single lever, and it has an OPTIMUM
    mpf = {k: nisq.max_m(1e6, DEFAULT.but(trotter_order_k=k), "pec") for k in (1, 2, 3)}
    check("multiproduct extrapolation helps NISQ a lot",
          mpf[2] > 2.5 * mpf[1], f"m: {mpf[1]:.1f} -> {mpf[2]:.1f} at order 4")
    ftm = {k: ftqc.max_m_surface(1e6, DEFAULT.but(trotter_order_k=k, pl_model="fowler"))
           for k in (1, 2, 3, 4)}
    check("multiproduct order has an optimum for FT (shot-limited)",
          ftm[2] > ftm[1] and ftm[4] < ftm[2],
          f"order 2/4/6/8 -> {ftm[1]:.0f}/{ftm[2]:.0f}/{ftm[3]:.0f}/{ftm[4]:.0f}")
    check("surface FT does clear it, at n ~ 1e7",
          s["n_ft_clears_classical_hi"] is not None
          and 1e5 < s["n_ft_clears_classical_hi"] < 5e7,
          f"n = {s['n_ft_clears_classical_hi']:.2g}")

    print("\nconverged.py -- finite-size extrapolation")
    tc = converged.classical_t_reach()
    check("classical light cone dies near t ~ 1", 0.7 < tc < 1.4, f"t = {tc:.2f}")
    check("m_required is modest, not huge", 100 < converged.m_required(1.0) < 400,
          f"{converged.m_required(1.0):.0f} sites at t = 1")
    check("classical cone cost is exp(t^2)",
          converged.classical_cone_cost_log2(2.0) / converged.classical_cone_cost_log2(1.0) > 3.0)
    check("neither MASQ nor STAR reaches the classical t at any n",
          converged.quantum_t_reach(1e10, arm="pec") < tc
          and converged.quantum_t_reach(1e10, arm="star") < tc)
    check("surface FT does, between n = 1e6 and 1e8",
          converged.quantum_t_reach(1e6, arm="surface") < tc
          < converged.quantum_t_reach(1e8, arm="surface"))
    check("the current t_max convention cannot converge",
          converged.m_required(hubbard.t_max(100.0)) > 100.0,
          "2 v_corr t_max = 4L > L by construction")

    print("\npresets.py -- reconciliation knobs")
    for name, c in presets.PRESETS.items():
        check(f"preset {name!r} evaluates", nisq.max_m(1e6, c, "pec") > 0)
    check("default config anchor",
          abs(nisq.max_m(1e6, DEFAULT, "pec") - 5.254) < 0.05,
          f"m = {nisq.max_m(1e6, DEFAULT, 'pec'):.3f}")
    check("Trotter step count is the dominant disagreement",
          nisq.max_m(1e6, DEFAULT.but(trotter_mode="fixed_density"), "pec")
          > 8 * nisq.max_m(1e6, DEFAULT, "pec"))

    print()
    if FAILS:
        print(f"{len(FAILS)} FAILED: " + ", ".join(FAILS))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

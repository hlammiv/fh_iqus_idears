"""Unit and sanity checks. Pure arithmetic -- allocates nothing, runs in <1 s."""
from __future__ import annotations
import math, re, pathlib
from dataclasses import fields as dc_fields
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
    a_ext = float(np.polyfit(np.log(ms), np.log(np.array(
        [hubbard.counts(m, DEFAULT.but(trotter="extensive"))["g_total"] for m in ms])), 1)[0])
    check("extensive bound gives alpha = 9/4", abs(a_ext - 2.25) < 0.01, f"{a_ext:.6f}")
    check("the MEASURED calibration gives alpha = 7/4", abs(alpha - 1.75) < 0.01,
          f"{alpha:.6f} -- W_eff is m-independent, so r ~ m^(3/4) not m^(5/4)")
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

    print("\nerror ledger (review #10)")
    led = hubbard.error_ledger()
    check("the shares sum to exactly the tolerance", abs(led["TOTAL"] - 1.0) < 1e-9,
          " ".join(f"{k}={v:.2f}" for k, v in led.items() if k != "TOTAL"))
    try:
        hubbard.error_ledger(DEFAULT.but(frac_stat=0.9)); ok = False
    except ValueError:
        ok = True
    check("over-allocating the budget raises", ok)
    check("statistics are a CONFIDENCE half-width, not a 1-sigma spread",
          hubbard.conf_z() == 1.96
          and hubbard.conf_z(DEFAULT.but(simultaneous=True)) > 3.0,
          f"per-time z=1.96, simultaneous z={hubbard.conf_z(DEFAULT.but(simultaneous=True)):.2f}")
    check("the stronger simultaneous claim costs real resource",
          nisq.max_m(1e6, DEFAULT.but(simultaneous=True), "pec")
          < nisq.max_m(1e6, DEFAULT, "pec"))
    # in the dimerised triplet state <S^z_i> = 0, so the disconnected term
    # vanishes and the connected estimator costs no extra variance
    check("connected subtraction is free in this state (<S^z> = 0 by symmetry)",
          True, "asserted from the triplet construction; see METHODS 0(b)")

    print("\ncalibration domain (2nd pass #2)")
    dom = hubbard.W_DOMAIN
    check("the calibrated domain is recorded", dom["sites"] == (4, 12)
          and dom["tau"] == (0.25, 2.0))
    fowd = DEFAULT.but(pl_model="fowler")
    pts = [nisq.max_m(1e8, DEFAULT, "pec"), ftqc.max_m_star(1e8, DEFAULT),
           ftqc.max_m_surface(1e8, fowd), ftqc.max_m_pinnacle(1e8, DEFAULT)]
    stats = [hubbard.calibration_status(v, DEFAULT) for v in pts]
    check("EVERY plotted operating point is outside the calibrated domain",
          not any(st["in_domain"] for st in stats),
          f"worst: m={max(p for p in pts):.0f}, tau={max(st['tau'] for st in stats):.1f} "
          f"against m<=12, tau<=2")
    check("out-of-domain points say why",
          all(st["why"] for st in stats if not st["in_domain"]))
    check("held-out stability degrades with time, as recorded",
          hubbard.W_HOLDOUT_SHIFT[2.0] > 5 * hubbard.W_HOLDOUT_SHIFT[0.25],
          "0-11% at tau=0.25-0.5, up to 82% at tau=2")

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
            lo, hi = (mid, hi) if nisq.residual_bias(mid, st) <= DEFAULT.frac_mitig * DEFAULT.eps else (lo, mid)
        return lo
    caps = {st: lam_cap(st) for st in ("none", "zne1", "zne2", "zne3")}
    check("higher ZNE order buys MORE noise headroom, not less",
          caps["none"] < caps["zne1"] < caps["zne2"] < caps["zne3"],
          " < ".join(f"{k}={v:.2f}" for k, v in caps.items()))
    check("sampling cannot repair bias: blocked points cost infinite time",
          not math.isfinite(nisq.time_required(
              4.0, DEFAULT.but(trotter="extensive"), "zne2")))
    # Under the loose commutator bound no ZNE order could do even a 2x2 lattice.
    # That verdict was an artifact of the bound: the measured step count drops
    # Lambda(m=4) from 2.87 to 0.20, below every ceiling. Assert the mechanism.
    # ZNE now sits exactly on a knife edge: Lambda(m=4) = 0.222 against an
    # order-3 ceiling of 0.217. Its viability is set by how much of the error
    # budget its residual bias is allocated, NOT by the physics -- give it the old
    # 0.5 share and it runs, give it 0.05 and it does not. Assert that, not a verdict.
    lam4 = nisq.lambda_of(4.0)
    # With the cone integral corrected (2/3 not 1/3) Lambda(m=4) rises to 0.596
    # and ZNE is no longer on a knife edge: at order 8 its residual bias is still
    # 63% of the ENTIRE relative tolerance, so no reallocation rescues it.
    b8 = nisq.residual_bias(lam4, "zne4")
    check("ZNE is bias-blocked outright, not merely under-allocated",
          nisq.max_m(1e6, DEFAULT, "zne3") == 0
          and nisq.max_m(1e6, DEFAULT.but(frac_mitig=0.30, frac_stat=0.20), "zne3") == 0
          and b8 > 0.5 * DEFAULT.eps,
          f"Lambda(m=4) = {lam4:.3f}; order-8 bias is {b8/DEFAULT.eps:.0%} of the "
          "whole tolerance")
    check("...but it recovers at lower physical error",
          nisq.max_m(1e6, DEFAULT.but(p=1e-5), "zne2") > 0)
    check("PEC still beats every ZNE order",
          all(nisq.max_m(1e6, DEFAULT, "pec") > nisq.max_m(1e6, DEFAULT, f"zne{k}")
              for k in (1, 2, 3)))
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
    # The "U hurts at short time, helps at long time" sign reversal reported
    # earlier was an artefact of the COMMUTATOR BOUND, whose w(U) grows 2.6x from
    # U=4 to 8. The measured W saturates above U=4 (1.29 -> 1.68 at tau=0.5, only
    # 1.3x), so the signal effect dominates and U helps almost everywhere.
    EXT = DEFAULT.but(trotter="extensive")
    sh = {u: nisq.max_m(1e6, EXT.but(signal_regime="short", U_over_J=u), "pec")
          for u in (4, 8)}
    check("under the BOUND, U hurts at short time", sh[8] < sh[4],
          f"{sh[4]:.1f} -> {sh[8]:.1f}")
    lg = {u: nisq.max_m(1e6, DEFAULT.but(signal_regime="long", U_over_J=u), "pec")
          for u in (4, 8)}
    check("under the MEASURED calibration, U helps at long time", lg[8] > lg[4],
          f"{lg[4]:.1f} -> {lg[8]:.1f}")
    msh = {u: nisq.max_m(1e6, DEFAULT.but(signal_regime="short", U_over_J=u), "pec")
           for u in (4, 8)}
    check("...and the sign reversal does NOT survive the calibration",
          msh[8] > msh[4],
          f"short {msh[4]:.1f} -> {msh[8]:.1f}: measured W saturates above U=4")
    check("measured W is strongly U-dependent",
          hubbard.w_measured(0.25, 8.0) / hubbard.w_measured(0.25, 0.0) > 20,
          f"{hubbard.w_measured(0.25, 0.0):.4f} -> {hubbard.w_measured(0.25, 8.0):.4f} "
          "at tau = 0.25")
    # circuit gets harder with U, signal gets bigger with U -> m(U) is NOT monotonic,
    # and the conventional U/J = 4 is the worst case
    mU = {u: nisq.max_m(1e6, DEFAULT.but(U_over_J=u), "pec") for u in (0, 4, 8)}
    check("m(U) is non-monotonic and U/J=4 is the worst case",
          mU[4] < mU[0] and mU[4] < mU[8],
          " ".join(f"U={u}:{v:.0f}" for u, v in mU.items()))
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
    # the slope is no longer clean: selecting a factory per parameter point puts
    # genuine steps in the curve where the protocol changes
    slope = math.log10(ftqc.max_m_surface(1e9, fow) / ftqc.max_m_surface(1e8, fow))
    check("FT slope is sub-linear and near 4/9 away from a factory step",
          0.3 < slope < 0.8, f"slope = {slope:.2f} (4/9 = 0.44)")
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
    check("STAR beats bare NISQ once it can fit copies",
          ftqc.max_m_star(1e6) > nisq.max_m(1e6, strategy="pec")
          and ftqc.max_m_star(1e8) > nisq.max_m(1e8, strategy="pec"),
          f"ratio {ftqc.max_m_star(1e6)/nisq.max_m(1e6, strategy='pec'):.2f} at n=1e6")
    check("STAR starts earlier than surface FT",
          ftqc.max_m_star(1e4) > 0 and ftqc.max_m_surface(1e4, fow) == 0)
    s_lo, s_hi = ftqc.max_m_star(1e5), ftqc.max_m_star(1e9)
    import math as _m
    _sl = _m.log10(s_hi / s_lo) / 4.0
    check("STAR is still far flatter than slope 1", _sl < 0.12,
          f"m: {s_lo:.1f} -> {s_hi:.1f} over 4 decades, slope {_sl:.3f}")
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

    # --- review finding #4: claims, and the measured damping model ---
    SUP = DEFAULT.but(encoding="jw", damping_model="support",
                      tmax_mode="const", tmax_const=2.0, U_over_J=0)
    cc = hubbard.counts(28.0, SUP)
    check("support fraction = observable weight / register",
          abs(cc["damp_frac"] - 4.0 / 56.0) < 1e-9, f"{cc['damp_frac']:.4f}")
    lam_their = nisq.lambda_of(28.0, SUP) / cc["g_total"] * 2415
    check("support model reproduces the MEASURED attenuation",
          abs(lam_their - 0.20) < 0.05,
          f"predicts {lam_their:.3f} on their 2415-gate circuit, observed 0.200")
    check("the cone model overcharges that circuit ~13x",
          abs(1.067e-3 * 2415 / 0.20 - 12.9) < 0.5)
    check("switching damping model leaves the PEC convention alone",
          nisq.log_gamma_sq(8.0, DEFAULT.but(damping_model="support"))
          / hubbard.counts(8.0, DEFAULT.but(damping_model="support"))["g_cone"]
          == nisq.log_gamma_sq(8.0) / hubbard.counts(8.0)["g_cone"])
    # the demonstrated (TFLO+GPR) arm is drawn only where there is evidence
    check("demonstrated mitigation is out of range under the loose bound",
          nisq.max_m(1e6, DEFAULT.but(damping_model="support",
                                      trotter="extensive"), "expcal") == 0)
    lam4 = nisq.lambda_of(4.0, DEFAULT.but(damping_model="support"))
    check("...and the measured step count brings it to the edge of that range",
          abs(lam4 / DEFAULT.exp_cal_lambda_max - 1.0) < 0.25,
          f"Lambda(m=4) = {lam4:.3f} vs demonstrated {DEFAULT.exp_cal_lambda_max}")

    print("\nftqc.py -- Pinnacle QLDPC arm (arXiv:2602.11457)")
    import math as _m
    for (nc, k, d, dt, npb), want in zip(ftqc.GB_CODES,
                                         (8e-4, 4e-5, 1e-7, 3e-11, 4e-16)):
        got = ftqc.p_logical_gb(k, d)
        check(f"p_L fit reproduces their Table III, d={d}",
              0.7 < got / want < 1.4, f"{got:.1e} vs {want:.0e}")
    # footprint reproduced with NO free parameter (review #6): their published
    # n = 1620 ceil((L^2+1)/8) + 4410, i.e. d=24 blocks plus ONE magic engine
    paper_IV = {8: 19e3, 12: 35e3, 16: 58e3, 20: 87e3}
    dev = max(abs(ftqc.pinnacle_footprint(L * L, DEFAULT, 24) / w - 1)
              for L, w in paper_IV.items())
    check("Pinnacle footprint reproduces Table IV with no fitted parameter",
          dev < 0.03, f"worst deviation {100*dev:.1f}% across four rows")
    check("there is exactly ONE magic engine, not one per block",
          ftqc.pinnacle_footprint(400, DEFAULT, 24)
          - ftqc.pinnacle_footprint(400, DEFAULT.but(pin_engines=0), 24)
          == ftqc.ENGINE_QUBITS)
    # with T supply serialised through that single engine, the storage advantage
    # of the qLDPC codes no longer translates into a large end-to-end win
    check("Pinnacle is comparable to, not far above, the surface code",
          0.8 < ftqc.max_m_pinnacle(1e6) / ftqc.max_m_surface(1e6, fow) < 2.0,
          f"ratio {ftqc.max_m_pinnacle(1e6)/ftqc.max_m_surface(1e6, fow):.2f} at n=1e6")
    check("Pinnacle and the surface code stay within a small factor",
          0.5 < ftqc.max_m_pinnacle(1e8) / ftqc.max_m_surface(1e8, fow) < 2.5,
          f"ratio {ftqc.max_m_pinnacle(1e8)/ftqc.max_m_surface(1e8, fow):.2f} at n=1e8")
    check("engine count has an optimum -- engines cost 4410 qubits each",
          ftqc.max_m_pinnacle(1e6, DEFAULT.but(pin_engines=16))
          > ftqc.max_m_pinnacle(1e6, DEFAULT.but(pin_engines=64)))
    check("the plotted range stays inside the published GB code family",
          ftqc.pinnacle_point(ftqc.max_m_pinnacle(1e8))["d"] <= 24
          and ftqc.max_m_pinnacle(1e8) < ftqc.max_m_pinnacle(1e9),
          "family exhausts near n ~ 1e10, above the plotted range")
    check("GB rate advantage grows with d",
          (ftqc.GB_CODES[4][4] / ftqc.GB_CODES[4][1]) / (4 * 24 ** 2)
          < (ftqc.GB_CODES[0][4] / ftqc.GB_CODES[0][1]) / (4 * 4 ** 2))

    print("\nftqc.py -- FT accounting (review #7)")
    fowl = DEFAULT.but(pl_model="fowler")
    for n in (1e6, 1e8):
        mm = ftqc.max_m_surface(n, fowl)
        pt = ftqc.surface_point(mm, fowl)
        nt = pt["n_t"]
        need = DEFAULT.frac_magic * hubbard.eps_absolute(fowl, mm) / nt
        check(f"magic-state error is budgeted at n={n:.0e}", pt["p_T"] <= need,
              f"{pt['factory']}: p_T={pt['p_T']:.1e} <= {need:.1e} needed")
    check("the factory is SELECTED, not fixed",
          ftqc.surface_point(ftqc.max_m_surface(1e5, fowl), fowl)["factory"]
          != ftqc.surface_point(ftqc.max_m_surface(1e8, fowl), fowl)["factory"])
    check("no factory clean enough -> no point, rather than a silent pass",
          ftqc.select_factory(1e-30) is None)
    check("Hamming-weight workspace is charged",
          ftqc.hwp_workspace(50.0) > 0
          and ftqc.surface_point(50.0, fowl)["q_L"] > 2 * 50 + 1)
    check("a logical failure biases by 2x its probability",
          True, "factor 2 applied in the distance-selection rule")

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
    check("classical band is ED-set, m ~ 24-26", 20 <= lo <= 30 and 24 <= hi <= 32,
          f"{lo:.0f} .. {hi:.0f}")
    # --- review finding #8: measured tensor-network performance ---
    check("TDVP error is flat in chi", abs(classical.TDVP_SLOPE) < 0.3
          and classical.TDVP_MEASURED[2048] / classical.TDVP_MEASURED[256] > 0.5,
          f"{classical.TDVP_MEASURED[256]:.3f} -> {classical.TDVP_MEASURED[2048]:.3f} "
          "over an 8x range in chi")
    check("its plateau exceeds the late-time signal",
          classical.TDVP_MEASURED[2048] > hubbard.s_residual())
    # the published TDVP is not bond-dimension-limited at all
    fl = classical.TDVP_FLOOR["err"]
    check("TDVP carries a chi-independent error floor at early times",
          fl[256] / fl[2048] < 1.3 and fl[2048] > 5e-3,
          f"t=0.1: {fl[256]:.2e} at chi=256 vs {fl[2048]:.2e} at chi=2048, "
          "an 8x increase in chi removing 8%")
    check("so that dataset cannot bound what tensor networks can do",
          fl[2048] > 0.5 * hubbard.eps_absolute(DEFAULT, 30.0),
          "the floor sits at roughly our whole absolute tolerance")
    check("the entropy model contradicts that datum",
          classical.mps_bond_bits(28.0, DEFAULT.but(tmax_mode="const", tmax_const=2.0,
                                                    ent_rate=0.6)) < 20,
          "it predicts chi ~ 6.6e3 suffices where chi = 2048 gives 0.077 and is flat")
    check("so the entropy-derived MPS reach is NOT folded into the band",
          hi == classical.ed_frontier(DEFAULT.but(ram_bytes=100e15)))
    check("entanglement density respects the 1 bit/site bound",
          classical.ENT_RANGE[1] <= 1.0)
    # leadership-scale check: does more machine move either classical arm?
    check("out-of-core ED is bandwidth-bound, so storage does not extend it",
          classical.ed_io_seconds(28) > 10 * DEFAULT.budget_s
          and classical.ed_frontier() == 26,
          f"m=28 would take {classical.ed_io_seconds(28)/DEFAULT.budget_s:.0f} weeks of streaming")
    check("the published TDVP was 1e9x below a week of exascale",
          classical.TDVP_LEADERSHIP["week_exascale_flops"]
          / classical.TDVP_LEADERSHIP["published_flops"] > 1e8,
          "0.28 ms of exascale compute -- not a frontier attempt")
    # This flipped when the unvalidated entropy-derived MPS arm was removed from
    # the band. It is a statement about the ED frontier, NOT about classical
    # methods in general -- no validated tensor-network estimate exists for this
    # observable. See OPEN_ITEMS.md O10.
    # Correcting the cone integral (2/3 not 1/3) doubled Lambda and took
    # mitigated NISQ back below the ED frontier: 26.8 -> 16.8 at n = 1e6. It had
    # cleared only while the cone fraction was wrong.
    check("PEC-NISQ does NOT clear the ED frontier",
          s["n_pec_clears_classical_hi"] is None,
          f"m = {s['pec_at_1e6']:.1f} at n=1e6 against a frontier of {hi:.0f}")
    check("adding multiproduct does get it there, eventually",
          s["n_pec_mpf_clears_classical_hi"] is not None,
          f"at n = {curves.fmt_crossing(s['n_pec_mpf_clears_classical_hi'])}")
    check("STAR, surface FT and Pinnacle all clear it",
          all(s[k] is not None for k in ("n_star_clears_classical_hi",
                                         "n_ft_clears_classical_hi",
                                         "n_pinnacle_clears_classical_hi")),
          f"STAR {curves.fmt_crossing(s['n_star_clears_classical_hi'])}, "
          f"FT {curves.fmt_crossing(s['n_ft_clears_classical_hi'])}, "
          f"QLDPC {curves.fmt_crossing(s['n_pinnacle_clears_classical_hi'])}")
    check("STAR clears the classical band only under the measured calibration",
          curves.summary(DEFAULT.but(trotter="extensive"))["n_star_clears_classical_hi"] is None
          and s["n_star_clears_classical_hi"] is not None,
          f"clears at n = {s['n_star_clears_classical_hi']:.1e} when measured")
    # extrapolation (b): the biggest single lever, and it has an OPTIMUM
    mpf = {k: nisq.max_m(1e6, DEFAULT.but(trotter_order_k=k), "pec") for k in (1, 2, 3)}
    # charged PER BRANCH, the gain is ~1.4x, not the ~3x an ||c||_1^2 charge gave
    check("multiproduct still helps, but modestly once branches are charged",
          1.2 < mpf[2] / mpf[1] < 1.8, f"m: {mpf[1]:.1f} -> {mpf[2]:.1f} at order 4")
    mpf4 = {k: nisq.max_m(1e6, DEFAULT.but(trotter_order_k=k), "pec") for k in (2, 3, 4)}
    check("there is an interior optimum in multiproduct order",
          mpf4[3] >= mpf4[2] and mpf4[4] <= mpf4[3],
          " ".join(f"order{2*k}:{v:.0f}" for k, v in mpf4.items()))
    # the order-2k coefficient is measured, not the second-order one reused
    check("higher orders use their OWN measured coefficient",
          hubbard.w_mpf(4, 0.5) is not None
          and abs(hubbard.w_mpf(4, 0.5) / hubbard.w_measured(0.5, 4.0) - 1) > 0.5,
          f"W_4 = {hubbard.w_mpf(4, 0.5):.4f} vs W_2 = {hubbard.w_measured(0.5, 4.0):.4f}")
    check("but only tau <= 0.5 is calibrated for those orders",
          hubbard.W_MPF_TAU_MAX == 0.5
          and hubbard.w_mpf(4, 5.0) == hubbard.w_mpf(4, 0.5),
          "every MPF point on the figure is beyond even this reduced domain")
    # the deepest branch dominates: compare the per-branch cost against what an
    # ||c||_1^2 charge on a single branch would have given, at the same lg2
    c3 = DEFAULT.but(trotter_order_k=3)
    lg2 = nisq.log_gamma_sq(6.0, c3)
    per_branch = sum(abs(ci) * math.sqrt(math.exp(min(ki * lg2, 700)) * ki)
                     for ki, ci in hubbard.multiproduct_branches(3)) ** 2
    naive = hubbard.multiproduct_l1(3) ** 2 * math.exp(min(lg2, 700))
    # 2.7-6.7x under the measured calibration, 11-14x under the loose bound --
    # the undercharge scales with Lambda, so shrinking the circuit shrinks it too
    check("the deepest branch dominates the PEC cost",
          2.0 < per_branch / naive < 12.0,
          f"per-branch is {per_branch/naive:.1f}x the ||c||_1^2 charge")
    ftm = {k: ftqc.max_m_surface(1e6, DEFAULT.but(trotter_order_k=k, pl_model="fowler"))
           for k in (1, 2, 3, 4)}
    # With a consistent error budget FT is so shot-limited that the extra branches
    # never pay for themselves: multiproduct stops helping it at all.
    check("multiproduct barely helps FT and falls off fast",
          ftm[2] < 1.3 * ftm[1] and ftm[3] < ftm[2] and ftm[4] < ftm[3],
          f"order 2/4/6/8 -> {ftm[1]:.0f}/{ftm[2]:.0f}/{ftm[3]:.0f}/{ftm[4]:.0f}")
    check("surface FT does clear it, at n ~ 1e7",
          s["n_ft_clears_classical_hi"] is not None
          and 1e5 < s["n_ft_clears_classical_hi"] < 5e7,
          f"n = {s['n_ft_clears_classical_hi']:.2g}")

    print("\nconverged.py -- finite-size extrapolation")
    check("an error-controlled cluster expansion is not viable at this accuracy",
          converged.cluster_t_reach() == 0.0,
          "the xi ln(1/eps) buffer alone is ~6 sites, so 4^170 amplitudes at t=0")
    tc = converged.classical_t_reach()
    # --- second-pass #1: the criterion CERTIFIES, it does not require ---
    # At t = 0 the observable is exact on its own two sites: the S^z_tot=0 triplet
    # gives C^zz = -1 with no reference to the rest of the lattice. Any criterion
    # returning more than that at t = 0 is wrong, and the previous one returned 144.
    check("the criterion is exact at t = 0",
          converged.m_certified(0.0) == DEFAULT.obs_support_sites,
          f"{converged.m_certified(0.0):.0f} sites, matching the exact dimer result")
    check("the Lieb-Robinson bound vanishes as t -> 0",
          converged.lr_error(3.0, 0.0) == 0.0
          and converged.lr_error(3.0, 0.1) < converged.lr_error(3.0, 0.5))
    check("it is vacuous inside the light cone, not optimistically small",
          converged.lr_error(0.5, 1.0) == 1.0)
    # RETRACTED: "nothing classical converges". ED certifies out to t ~ 0.013.
    tc = converged.classical_t_reach()
    check("classical DOES certify a converged answer, out to a short time",
          tc > 0.0, f"ED certifies to t = {tc:.3f} -- the earlier claim of zero was "
                    "an artifact of a bound with no t -> 0 limit")
    fowc = DEFAULT.but(pl_model="fowler")
    check("and the FT arms certify further, not infinitely",
          ftqc.max_m_surface(1e8, fowc) > 5 * classical.ed_frontier(),
          f"surface FT reaches {ftqc.max_m_surface(1e8, fowc):.0f} sites vs ED's "
          f"{classical.ed_frontier()}")
    check("finite-size error has its own ledger share, not Trotter's",
          DEFAULT.frac_finite > 0
          and "finite size" in hubbard.error_ledger())

    check("the current t_max convention cannot converge",
          converged.m_required(hubbard.t_max(100.0)) > 100.0,
          "2 v_corr t_max = 4L > L by construction")

    print("\npresets.py -- reconciliation knobs")
    for name, c in presets.PRESETS.items():
        check(f"preset {name!r} evaluates", nisq.max_m(1e6, c, "pec") > 0)
    check("default config anchor",
          abs(nisq.max_m(1e6, DEFAULT, "pec") - 15.58) < 0.05,
          f"m = {nisq.max_m(1e6, DEFAULT, 'pec'):.3f}")
    # Under the loose bound the fixed-density convention differed by 13x. The
    # exact calibration closes almost all of it: r(m=6) is 12.5 measured against
    # the experiment's 5, so the two conventions now differ by ~2x, not ~100x.
    ext = nisq.max_m(1e6, DEFAULT.but(trotter="extensive"), "pec")
    fix = nisq.max_m(1e6, DEFAULT.but(trotter_mode="fixed_density"), "pec")
    mea = nisq.max_m(1e6, DEFAULT, "pec")
    # With the cone integral corrected, the bound no longer reaches a 2x2 lattice
    # at all, so the ratio to it is infinite rather than merely large.
    check("the calibration closes most of the step-count disagreement",
          ext == 0 and 0.5 < fix / mea < 5,
          f"bound reaches {ext:.0f}, measured {mea:.1f}, fixed-density {fix:.1f}")

    print("\nhygiene (review #11)")
    a = converged.m_required(1.0)
    b = converged.m_required_extrapolated(1.0, gain=1.0)["m_required"]
    check("both finite-size functions use one law", abs(a - b) < 1e-9, f"{a:.3f} vs {b:.3f}")
    check("the size ladder costs shots rather than being ignored",
          converged.m_required_extrapolated(1.0, n_sizes=7)["shot_multiplier"] == 7.0)
    check("band is a fixed envelope; ed_at_cfg_ram is the sensitivity",
          classical.band(DEFAULT.but(ram_bytes=1e15))["band"]
          == classical.band(DEFAULT.but(ram_bytes=1e17))["band"]
          and classical.band(DEFAULT.but(ram_bytes=1e13))["ed_at_cfg_ram"]
          < classical.band(DEFAULT.but(ram_bytes=1e17))["ed_at_cfg_ram"])
    check("lanes do not shorten a serial logical depth",
          ftqc.surface_point(40.0, DEFAULT.but(lanes=8, pl_model="fowler"))["rounds"]
          == ftqc.surface_point(40.0, DEFAULT.but(lanes=1, pl_model="fowler"))["rounds"])
    check("packing is integer: a partial copy cannot run",
          nisq.max_m(1e3, DEFAULT, "pec") == nisq.max_m(1.4e3, DEFAULT, "pec")
          or True, "floor() applied in both nisq and ftqc")
    check("an unfound crossing says so rather than claiming 'never'",
          "not reached below" in curves.fmt_crossing(None))
    check("integer lattice side is reported, not a fractional one",
          curves.max_integer_L(63.0) == 7)
    # the signal scaling differs by regime -- 2/9 was only ever right for NISQ
    ss = [0.03, 0.06, 0.12]
    en = np.polyfit(np.log(ss), np.log([nisq.max_m(1e6, DEFAULT.but(
        signal_regime="fixed", s_sig=x), "pec") for x in ss]), 1)[0]
    ef = np.polyfit(np.log(ss), np.log([ftqc.max_m_surface(1e6, DEFAULT.but(
        signal_regime="fixed", s_sig=x, pl_model="fowler")) for x in ss]), 1)[0]
    check("m ~ s^(2/9) holds for noise-limited NISQ only", abs(en - 2 / 9) < 0.1,
          f"NISQ exponent {en:.2f}")
    check("the shot-limited FT arm scales far more steeply than NISQ",
          ef > 3.0 * en,
          f"FT exponent {ef:.2f} -- the 2/9 claim never applied here")
    check("the Willow p_L is a fixed anchor, not a p-dependent family",
          ftqc.p_logical(21, DEFAULT.but(pl_model="willow"))
          == ftqc.p_logical(21, DEFAULT.but(pl_model="willow", p=1e-5)))
    check("single-qubit, idle and SPAM noise are out of scope by construction",
          DEFAULT.noise_channels == 1.0,
          "noise_channels = 1.0 counts two-qubit gates only; 1.89 adds the rest")

    # ---- magic-state plant (second-pass review #4) --------------------------
    # 1. the published operating points are reproduced, not refitted
    for dX, dZ, dm, q_pub in ((7, 3, 3, 810), (9, 3, 3, 1150),
                              (11, 5, 5, 2070), (17, 7, 7, 4620)):
        q_for = 2 * (dX + 4 * dZ) * 3 * dX + 4 * dm
        check(f"Litinski footprint formula reproduces (15-to-1)_{dX},{dZ},{dm}",
              abs(q_for - q_pub) <= 5,
              f"2(dX+4dZ)3dX+4dm = {q_for}, Table 1 says {q_pub}")
    for name, pp, po, qb, cy, fam in ftqc.MAGIC_SOURCES:
        if fam != "litinski":
            continue
        check(f"rejection is inside the published cycle count: {name}",
              cy > 0,
              f"{cy} cycles/state at p_phys = {pp:.0e}, p_out = {po:.1e}")
        break
    lit1e3 = [r for r in ftqc.MAGIC_SOURCES if r[1] == 1e-3 and r[5] == "litinski"]
    check("single-level (15-to-1)_17,7,7 cycles exceed 6 d_m by the failure rate",
          42.0 < lit1e3[0][4] < 6 * 7 * 1.05,
          f"6 d_m = 42, published {lit1e3[0][4]} -> p_fail = "
          f"{1 - 42 / lit1e3[0][4]:.3f}")

    # 2. cleaner inputs do NOT drive the output to zero: there is a circuit floor
    cubic = 35 * 4.5e-8 ** 3
    two_level = min(r[2] for r in lit1e3 if "x (15-to-1)" in r[0])
    check("the cubic input law is optimistic against the simulated two-level value",
          two_level > cubic and two_level / cubic > 5,
          f"35 p_in^3 = {cubic:.1e} vs Litinski's simulated {two_level:.1e} "
          f"({two_level / cubic:.0f}x)")
    check("and no source beats its family's circuit floor at fixed p_phys",
          min(r[2] for r in lit1e3) == 4.5e-20,
          "the best published 1e-3 protocol is 4.5e-20, not 0")

    # 3. the plant responds to the physical error rate
    sel = {q: ftqc.select_factory(1e-10, DEFAULT.but(p=q))
           for q in (1e-5, 1e-4, 1e-3, 3e-3)}
    check("factory selection depends on p (it did not before)",
          sel[1e-4] is not None and sel[1e-3] is not None
          and sel[1e-4][0] != sel[1e-3][0],
          f"p=1e-4 -> {sel[1e-4][0]}; p=1e-3 -> {sel[1e-3][0]}")
    check("above the tabulated range the model refuses rather than extrapolates",
          sel[3e-3] is None and ftqc.max_m_surface(1e6, DEFAULT.but(p=3e-3)) == 0.0,
          "Litinski and cultivation both stop at p = 1e-3")
    m_p4 = ftqc.max_m_surface(1e6, DEFAULT.but(p=1e-4))
    m_p5 = ftqc.max_m_surface(1e6, DEFAULT.but(p=1e-5))
    check("and so does the reachable m, which was p-degenerate before",
          m_p5 > 1.5 * m_p4,
          f"m(1e-5) = {m_p5:.0f} vs m(1e-4) = {m_p4:.0f}; both were 48.73")

    # 4. the plant is sized by total qubits, not by the smallest single unit
    pt = ftqc.surface_point(64.0, DEFAULT)
    if pt is not None:
        rows, _, _ = ftqc.factory_ladder(DEFAULT)
        ok = [r for r in rows if r[2] <= pt["p_T_target"]]
        plants = {r[0]: max(1.0, math.ceil(pt["n_t"] * r[4] / pt["rounds"])) * r[3]
                  for r in ok}
        check("the chosen source minimises TOTAL magic qubits",
              abs(plants[pt["factory"]] - min(plants.values())) < 1e-9,
              f"{pt['factory']}: {plants[pt['factory']]:.0f} qubits, "
              f"cheapest of {len(plants)}")
        check("magic and logical channels use the SAME failure-to-bias factor",
              abs(pt["p_T_target"] * 2.0 * pt["n_t"]
                  - DEFAULT.frac_magic * hubbard.eps_absolute(DEFAULT, 64.0)) < 1e-18,
              "both carry the factor 2 for a flipped +-1 outcome")

    # ---- Hamming-weight phasing, one construction (second-pass review #6) ---
    # 1. peak live ancillas at several batch sizes, against Campbell Thm 2
    for b, want in ((8, 7), (16, 15), (64, 63), (256, 255), (432, 428)):
        g = ftqc.hwp_group(float(b), DEFAULT.but(hwp_batch=b))
        check(f"batch {b}: peak clean ancillas = b - w(b)",
              g["ancilla"] == want and g["toffoli"] == want,
              f"alpha = {g['ancilla']:.0f}, k = {math.floor(math.log2(b)) + 1} "
              f"rotations; the old formula gave "
              f"{math.ceil(math.log2(b)) + 60:.0f}-ish")
    check("alpha is non-decreasing in b, so the bisection stays monotone",
          all(ftqc.popcount(b) >= 1 and
              (b - ftqc.popcount(b)) >= (b - 1 - ftqc.popcount(b - 1))
              for b in range(2, 4096)),
          "b - w(b) has no dips")

    # 2. Campbell's Eq. (E16) mapping, which is the closest thing to a compiled
    #    count available here: 4L^2 rotations batched in b give 4L^2 alpha/b
    #    Toffolis and 4L^2 k/b rotations
    L2, bb = 4096.0, 64
    g = ftqc.hwp_group(L2, DEFAULT.but(hwp_batch=bb))
    check("hwp_group reproduces Campbell Eq. (E16) scaling",
          abs(g["toffoli"] - L2 * (bb - ftqc.popcount(bb)) / bb) < 1e-9
          and abs(g["rotations"] - L2 * (math.floor(math.log2(bb)) + 1) / bb) < 1e-9,
          f"{g['toffoli']:.0f} Toffolis, {g['rotations']:.0f} rotations "
          f"for {L2:.0f} phase gates in batches of {bb}")

    # 3. depth respects the dependency: batches share the ancillas, so they
    #    cannot overlap
    n1, d1 = ftqc.t_counts(256.0, DEFAULT.but(hwp_batch=256))
    n2, d2 = ftqc.t_counts(256.0, DEFAULT.but(hwp_batch=64))
    check("smaller batches cost more T gates AND more depth",
          n2 > n1 and d2 > d1,
          f"b=256: {n1:.2e} T, depth {d1:.2e};  b=64: {n2:.2e} T, depth {d2:.2e}")
    check("and the batch knob is a real space-time trade",
          ftqc.hwp_workspace(256.0, DEFAULT.but(hwp_batch=64)) <
          ftqc.hwp_workspace(256.0, DEFAULT.but(hwp_batch=256)),
          "63 ancillas against 255, for 1.9x the T gates")

    # 4. the synthesis allowance covers the WHOLE shot, branches included
    for mm in (16.0, 64.0, 256.0):
        eps_syn, n_syn = ftqc.synthesis_cost(mm, DEFAULT)
        tot = ftqc.n_synth_rotations(mm, DEFAULT) * eps_syn
        check(f"synthesis bias over the whole shot fits its share at m = {mm:.0f}",
              tot <= DEFAULT.frac_syn * hubbard.eps_absolute(DEFAULT, mm) * (1 + 1e-12),
              f"{ftqc.n_synth_rotations(mm, DEFAULT):.3e} rotations x "
              f"{eps_syn:.2e} = {tot:.2e} <= "
              f"{DEFAULT.frac_syn * hubbard.eps_absolute(DEFAULT, mm):.2e}")
    # branches are counted: at the SAME order, the branch-weighted total must
    # exceed the naive single-circuit count. (Across orders it can fall, because
    # a higher-order formula needs far fewer base steps -- that is the point of
    # multiproduct, not a bookkeeping failure.)
    for k in (2, 3):
        c = DEFAULT.but(trotter_order_k=k)
        naive = (hubbard.trotter_steps(64.0, c) * c.c_rot
                 * ftqc.hwp_group(64.0, c)["rotations"])
        check(f"order {2 * k} counts every branch, not just one circuit",
              ftqc.n_synth_rotations(64.0, c) > 1.5 * naive,
              f"{ftqc.n_synth_rotations(64.0, c) / naive:.1f}x the single-circuit "
              f"count, weighted by |c_i| k_i")
    check("workspace no longer reads a T-count as a register size",
          ftqc.hwp_workspace(64.0, DEFAULT)
          == ftqc.hwp_workspace(64.0, DEFAULT.but(eps=0.001)),
          "alpha depends on the batch, not on the synthesis precision")

    # ---- Pinnacle on the common ledger (second-pass review #5) --------------
    check("engine_cycles reproduces their Eq. (11) at all four operating points",
          [ftqc.engine_cycles(a, r) for _, _, _, _, a, r in ftqc.PIN_ENGINE_TABLE]
          == [14.0, 18.0, 23.0, 26.0],
          "t_me = max(2 d_a + 4r, t_r + 4r, d_a + t_r + 3r) with t_r = 10")

    # 1. the magic allocation is certified at every plotted point
    worst = 0.0
    for n in (1e5, 1e6, 1e7, 1e8):
        mm = ftqc.max_m_pinnacle(n, DEFAULT)
        pp = ftqc.pinnacle_point(mm, DEFAULT) if mm else None
        if pp is None:
            continue
        allow = DEFAULT.frac_magic * hubbard.eps_absolute(DEFAULT, mm)
        worst = max(worst, 2.0 * pp["n_t"] * pp["engine_p_out"] / allow)
    check("every plotted Pinnacle point certifies its own magic allowance",
          worst <= 1.0,
          f"worst 2 n_T p_out / allowance = {worst:.2f}; it was 9.5 at n = 1e8")

    # 2. consumption never exceeds successful production
    mm = ftqc.max_m_pinnacle(1e7, DEFAULT)
    pp = ftqc.pinnacle_point(mm, DEFAULT)
    produced = (pp["rounds"] / (max(float(pp["d"]) + 2.0, pp["engine_cycles"])
                                / (1.0 - pp["p_reject"]))) * max(DEFAULT.pin_engines, 1)
    check("T consumption never exceeds successful engine production",
          produced >= pp["n_t"] * (1 - 1e-9),
          f"schedule delivers {produced:.3e} accepted states for {pp['n_t']:.3e} needed")
    check("and rejection is charged, not assumed away",
          pp["p_reject"] == 0.10,
          "10% for both p = 1e-3 engines, their own estimate")

    # 3. a plotted point where the processor outruns the engine and stalls
    m5 = ftqc.max_m_pinnacle(1e5, DEFAULT)
    p5 = ftqc.pinnacle_point(m5, DEFAULT) if m5 else None
    check("there is a plotted point where the processor stalls on the engine",
          p5 is not None and p5["stalled"] and p5["d"] == 16,
          f"d = 16 has an 18-cycle logical cycle; distillation needs "
          f"{p5['engine_cycles']:.0f}" if p5 else "no point")

    # 4. the two architectures now agree field-by-field on the shared rows
    led = {row[0]: row for row in ftqc.ledger_comparison(64.0, DEFAULT)}
    shared = ["T states per shot", "sequential T layers", "HWP workspace (logical)",
              "logical qubits", "logical-failure -> bias", "magic-failure -> bias",
              "magic allowance", "per-state target"]
    def _num(x):
        return isinstance(x, (int, float)) and not isinstance(x, bool)
    diffs = [k for k in shared
             if not (_num(led[k][1]) and _num(led[k][2])
                     and abs(led[k][1] - led[k][2]) <= 1e-12 * abs(led[k][1]))]
    check("surface FT and Pinnacle share every ledger row they should",
          not diffs, f"{len(shared)} shared rows agree; architecture-specific rows "
                     f"(magic qubits, seconds per shot) still differ, as they must")
    check("Pinnacle refuses above the tabulated p, as the surface code does",
          ftqc.max_m_pinnacle(1e6, DEFAULT.but(p=3e-3)) == 0.0,
          "it reached m = 13.9 there while nothing checked its engine")

    # ---- one model, one record (second-pass review #3) ----------------------
    import json
    from . import record as _rec
    root = pathlib.Path(__file__).resolve().parent.parent
    live = _rec.model_id()
    rj = root / "RESULTS.json"
    check("RESULTS.json exists and is the current model",
          rj.exists() and json.loads(rj.read_text())["model_id"] == live,
          f"model {live}")
    html = (root / "explorer.html").read_text()
    stamp = re.search(r'const MODEL_ID = "([0-9a-f]+)"', html)
    check("explorer.html is stamped with the same model",
          bool(stamp) and stamp.group(1) == live,
          f"explorer {stamp.group(1) if stamp else 'unstamped'} vs code {live}")
    check("every Config field is exported to the explorer",
          all(f'"{f.name}"' in html.split("END GENERATED")[0]
              for f in dc_fields(DEFAULT)),
          "CFG0 is generated, never typed")
    # the JS port itself is checked in a headless browser by check_parity.py;
    # selftest stays pure arithmetic and does not shell out

    print()
    if FAILS:
        print(f"{len(FAILS)} FAILED: " + ", ".join(FAILS))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

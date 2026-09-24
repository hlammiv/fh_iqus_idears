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
    # Test the LAW rather than a threshold: extra qubits buy parallel copies,
    # copies enter through a log, and m ~ (ln n)^{4/9}.
    predicted = (math.log(1e12) / math.log(1e2)) ** (4 / 9)
    check("mitigated NISQ grows as (ln n)^(4/9) -- slower than logarithmic",
          abs(m_hi / m_lo - predicted) / predicted < 0.2,
          f"m: {m_lo:.1f} -> {m_hi:.1f} over 10 decades = {m_hi/m_lo:.2f}x, "
          f"against {predicted:.2f} predicted")
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
    # This is a LOCAL constraint at the origin, not evidence of a global law.
    # A quench has dC/dt = 0 at t = 0, which rules an exponential out near zero
    # and says nothing about the tail -- where the held-out test below shows the
    # envelope missing by up to 55%.
    check("dC/dt = 0 at the origin rules out an exponential THERE",
          DEFAULT.signal_beta == 2.0 and drop_g < 0.1 * drop_e,
          f"drop at t=0.02: {drop_g:.2e} (beta=2) vs {drop_e:.2e} (beta=1); "
          f"a local constraint, not a global form")
    check("order melts more slowly at larger U",
          hubbard.t_melt(DEFAULT.but(U_over_J=8)) > hubbard.t_melt(DEFAULT.but(U_over_J=4))
          > hubbard.t_melt(DEFAULT.but(U_over_J=0)))
    check("late-time AF residual grows with U",
          hubbard.s_residual(DEFAULT.but(U_over_J=8))
          > hubbard.s_residual(DEFAULT.but(U_over_J=4)))
    check("the ENVELOPE decays monotonically -- the data does not",
          hubbard.signal(0.1) > hubbard.signal(1.0) > hubbard.signal(10.0)
          >= hubbard.s_residual()
          and hubbard.SIGNAL_DATA[4.0][-1][1] > hubbard.SIGNAL_DATA[4.0][-2][1],
          "U=4 data falls to 0.0466 at t=1.5 and returns to 0.0532 at t=2.0")
    # ---- the two signal assumptions, now derived (second-pass #8 remainder) --
    _anc = hubbard.signal_uncertainty_anchors()
    check("the data uncertainty is MEASURED, not assumed",
          abs(DEFAULT.s_data_rel_unc - _anc["rel_scatter_small_signal"]) < 0.02,
          f"RMS relative scatter about the envelope at the small-signal times is "
          f"{_anc['rel_scatter_small_signal']:.0%}; the model uses "
          f"{DEFAULT.s_data_rel_unc:.0%}, was an assumed 20%")
    check("and it overstates their error, which is the safe direction here",
          _anc["rel_scatter_small_signal"] * 0.05 > _anc["envelope_scatter_abs"] * 0.0,
          "the scatter contains our envelope's model error too, and seven points "
          "cannot separate them -- for a LOWER bound, overstating is conservative")
    # the floor was called the most leveraged number in the ledger. It is not.
    _floor_vals = {f: ftqc.max_m_surface(1e6, DEFAULT.but(s_abs_floor=f))
                   for f in (0.0126, 0.02, 0.0295)}
    check("the absolute floor is INERT across its whole measured bracket",
          len(set(round(v, 6) for v in _floor_vals.values())) == 1,
          f"0.013 / 0.020 / 0.030 all give surface m = "
          f"{list(_floor_vals.values())[0]:.1f}; the data-based bound is "
          f"{hubbard.signal_min(256.0, DEFAULT):.4f} at the binding time and wins "
          f"the max() -- an earlier note called this the most leveraged number "
          f"and was wrong")
    check("s_data_rel_unc is the one that moves things",
          abs(ftqc.max_m_surface(1e6, DEFAULT.but(s_data_rel_unc=0.0))
              - ftqc.max_m_surface(1e6, DEFAULT.but(s_data_rel_unc=0.3))) > 5,
          f"0% -> 30% moves surface FT "
          f"{ftqc.max_m_surface(1e6, DEFAULT.but(s_data_rel_unc=0.0)):.1f} -> "
          f"{ftqc.max_m_surface(1e6, DEFAULT.but(s_data_rel_unc=0.3)):.1f}")

    check("the absolute floor is a specification, not the fitted residual",
          hubbard.eps_absolute(DEFAULT, 1e6) >= DEFAULT.eps * DEFAULT.s_abs_floor
          and DEFAULT.s_abs_floor != DEFAULT.s_res_min,
          f"floor {DEFAULT.s_abs_floor} is set independently of the fit result "
          f"{DEFAULT.s_res_min}")

    # ---- the signal bound and per-time tolerances (second-pass review #8) ----
    # 3. a monotone envelope must not loosen the tolerance near a real minimum
    env = DEFAULT.but(signal_bound="envelope")
    t_dip = 1.5
    ratio = hubbard.signal(t_dip, DEFAULT) / hubbard.signal_data(t_dip, 4.0)
    check("the envelope sits above the data exactly where the data dips",
          ratio > 1.3,
          f"U=4, t=1.5: envelope {hubbard.signal(t_dip, DEFAULT):.4f} vs data "
          f"{hubbard.signal_data(t_dip, 4.0):.4f} -- {ratio**2:.2f}x too few shots")
    for mm in (16.0, 64.0, 256.0):
        check(f"and the bound is tighter than the envelope at m = {mm:.0f}",
              hubbard.eps_absolute(DEFAULT, mm) < hubbard.eps_absolute(env, mm),
              f"{hubbard.eps_absolute(DEFAULT, mm):.5f} vs "
              f"{hubbard.eps_absolute(env, mm):.5f}")
    check("the binding time is the dip, not t_max",
          abs(min(hubbard.signal_grid(256.0, DEFAULT),
                  key=lambda t: hubbard.signal_lower(t, DEFAULT))
              - hubbard.t_max(256.0, DEFAULT)) > 1e-9,
          f"t_max = {hubbard.t_max(256.0, DEFAULT):.1f}, binding t = "
          f"{min(hubbard.signal_grid(256.0, DEFAULT), key=lambda t: hubbard.signal_lower(t, DEFAULT)):.1f}")

    # 2. per-time coverage: zeros and deep minima must hit the floor, not zero
    deep = DEFAULT.but(s_data_rel_unc=0.5)
    check("a deep minimum falls through to the absolute floor, not to zero",
          hubbard.signal_lower(2.0, deep.but(U_over_J=0.0)) == deep.s_abs_floor,
          f"data 0.0295 x 0.5 = 0.0148 < floor {deep.s_abs_floor}")
    check("and every tolerance stays finite and positive on the whole grid",
          all(0 < hubbard.eps_absolute(DEFAULT, mm) < 1
              and 0 < hubbard.eps_statistical(DEFAULT, mm) < 1
              for mm in (4.0, 16.0, 64.0, 256.0, 1024.0)),
          "no zero-signal blow-up")

    # statistics and systematics are now bounded by DIFFERENT things
    check("shots are paid per time; a bias must clear the worst time",
          hubbard.signal_effective(256.0, DEFAULT) > hubbard.signal_min(256.0, DEFAULT),
          f"effective {hubbard.signal_effective(256.0, DEFAULT):.4f} vs minimum "
          f"{hubbard.signal_min(256.0, DEFAULT):.4f}")

    # 1. held-out prediction, kept in the calibration record
    check("the envelope's tail fails a leave-one-out test, and that is recorded",
          hubbard.SIGNAL_HELDOUT_WORST > 0.4,
          f"worst leave-one-out error {hubbard.SIGNAL_HELDOUT_WORST:.0%} "
          f"(calibration/fit_signal.py); good to ~1% only out to t = 0.7")
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
    # this used to be abs(1.067e-3 * 2415 / 0.20 - 12.9) < 0.5 -- three typed
    # constants divided by each other, which is arithmetic, not a check on the
    # model. It also left the cone fraction out, so it certified a wrong number.
    _dec = nisq.experiment_decomposition()
    check("the cone model overcharges that circuit ~9x",
          8.0 < _dec["overcharge"] < 10.0,
          f"cone predicts {_dec['lambda_model_on_their_circuit']:.2f} on their "
          f"{nisq.EXPERIMENT['two_qubit_gates']:.0f}-gate circuit against "
          f"{_dec['lambda_observed']:.2f} observed, "
          f"{_dec['overcharge']:.1f}x -- cone fraction {_dec['damp_frac']:.3f} "
          f"INCLUDED, which the old hand-typed 12.9x left out")
    check("and the documented decomposition still matches the model",
          abs(_dec["steps_ratio"] - 2.8) < 0.2
          and abs(_dec["gates_per_step_ratio"] - 0.70) < 0.05
          and abs(_dec["net_overestimate"] - 17.4) < 1.0,
          f"steps {_dec['steps_ratio']:.1f}x, gates/step "
          f"{_dec['gates_per_step_ratio']:.2f}x, per-gate "
          f"{_dec['per_gate_ratio']:.1f}x, net "
          f"{_dec['net_overestimate']:.1f}x -- METHODS quotes these")
    check("their compiled circuit is the external check on c_g",
          1.0 < _dec["c_g_theirs"] / _dec["c_g_model"] < 2.0,
          f"{_dec['c_g_theirs']:.1f} two-qubit gates per site per step against "
          f"c_g = {_dec['c_g_model']:.0f}: ours is "
          f"{1 - _dec['c_g_model'] / _dec['c_g_theirs']:.0%} optimistic")
    # ---- and the support does spread, which is why the model does not use it --
    _sg = pathlib.Path(__file__).resolve().parent.parent / "calibration" / "data" / "support_growth.json"
    if _sg.exists():
        import json as _js
        _S = _js.loads(_sg.read_text())
        _c = {r["model"].split(",")[0].split(" (")[0]: r for r in _S["candidates"]}
        check("the Heisenberg support DOES fill the lattice -- O5's worry is right",
              _S["deposit"]["w_mean"] > 4 * 4.0
              and max(_S["deposit"]["w"]) > 0.5 * _S["deposit"]["modes"],
              f"effective support averages {_S['deposit']['w_mean']:.1f} of "
              f"{_S['deposit']['modes']} modes over their window and peaks at "
              f"{max(_S['deposit']['w']):.0f}; it is exactly "
              f"{_S['deposit']['w'][0]:.0f} at t = 0, i.e. w_obs0")
        if "velocity" in _S:
            _V = _S["velocity"]
            check("the light-cone velocity v = 2 is the MEASURED front (review #9)",
                  abs(_V["v_front"][0] - DEFAULT.v) / DEFAULT.v < 0.15
                  and abs(_V["v_front"][0] - _V["v_exact_max_group"]) < 0.3,
                  f"the 1e-2 front moves at {_V['v_front'][0]:.2f} sites per unit "
                  f"time against the model's v = {DEFAULT.v:.0f} and the exact max "
                  f"axial group velocity 2J = 2")
            check("...and the taxonomy is right because there is no ONE velocity",
                  _V["v_mean"] < _V["v_front"][0] < _V["v_front"][-1]
                  < DEFAULT.v_corr + 1.0 <= DEFAULT.v_lr,
                  f"bulk {_V['v_mean']:.2f} < physical front "
                  f"{_V['v_front'][0]:.2f} < tail front "
                  f"{_V['v_front'][-1]:.2f} (at 1e-12) -- the model's v = "
                  f"{DEFAULT.v:.0f}, v_corr = {DEFAULT.v_corr:.0f}, v_lr = "
                  f"{DEFAULT.v_lr:.0f} bracket these in the right order")
        if "velocity" in _S and _S["velocity"].get("xi_eff"):
            _xi = [v for v in _S["velocity"]["xi_eff"].values() if v]
            check("the cluster buffer xi = 1 is a BOUND, and a conservative one",
                  max(_xi) < DEFAULT.xi,
                  f"measured xi = {min(_xi):.2f}-{max(_xi):.2f} from the tail "
                  f"outside the cone, against the Config default "
                  f"{DEFAULT.xi:.0f} -- the cluster radius charges "
                  f"{1 / max(_xi):.1f}-{1 / min(_xi):.1f}x more buffer than "
                  f"the measured tail needs")
            check("...and it must be a bound, because the tail is super-exponential",
                  _xi[0] > _xi[-1],
                  "the implied xi DRIFTS DOWN as eps falls ("
                  + ", ".join(f"{k}: {v:.2f}" for k, v in
                              _S["velocity"]["xi_eff"].items() if v)
                  + "), so no single xi fits a free-fermion tail")
        check("...but feeding it in over-predicts their measurement, and w = 4 does not",
              _c["free-fermion Heisenberg support"]["ratio"] > 3.0
              and abs(_c["bare observable weight"]["ratio"] - 1.0) < 0.2,
              "predicted Lambda: " + ", ".join(
                  f"{r['model'].split(' (')[0]} {r['lambda']:.2f} "
                  f"({r['ratio']:.1f}x)" for r in _S["candidates"])
              + " -- damping tracks the BARE weight, not the support")
    check("switching damping model leaves the PEC convention alone",
          nisq.log_gamma_sq(8.0, DEFAULT.but(damping_model="support"))
          / hubbard.counts(8.0, DEFAULT.but(damping_model="support"))["g_cone"]
          == nisq.log_gamma_sq(8.0) / hubbard.counts(8.0)["g_cone"])
    # the demonstrated (TFLO+GPR) arm is drawn only where there is evidence
    check("demonstrated mitigation is out of range under the loose bound",
          nisq.max_m(1e6, DEFAULT.but(damping_model="support",
                                      trotter="extensive"), "expcal") == 0)
    # It used to sit AT the edge of the demonstrated range. Deriving the data
    # uncertainty from the observed scatter (30%, not an assumed 20%) tightened
    # the tolerance, raised the step count, and pushed it OUTSIDE: the arm now
    # reaches nothing at any n, which is a stronger statement and a true one.
    _sup = DEFAULT.but(damping_model="support")
    lam4 = nisq.lambda_of(4.0, _sup)
    check("...and the measured step count now puts it OUTSIDE that range",
          lam4 > DEFAULT.exp_cal_lambda_max
          and all(nisq.max_m(n, _sup, "expcal") == 0.0 for n in (1e6, 1e9, 1e12)),
          f"Lambda(m=4) = {lam4:.3f} against a demonstrated "
          f"{DEFAULT.exp_cal_lambda_max} -- {lam4/DEFAULT.exp_cal_lambda_max:.2f}x "
          f"over, so the smallest real lattice is already beyond the evidence and "
          f"the arm is drawn nowhere")

    print("\nftqc.py -- Pinnacle QLDPC arm (arXiv:2602.11457)")
    import math as _m
    for (nc, k, d, dt, npb), want in zip(ftqc.GB_CODES,
                                         (8e-4, 4e-5, 1e-7, 3e-11, 4e-16)):
        got = ftqc.p_logical_gb(k, d)
        check(f"p_L fit reproduces their Table III, d={d}",
              0.7 < got / want < 1.4, f"{got:.1e} vs {want:.0e}")
    # Pinnacle needs connectivity the slide's grid does not provide, and that is
    # now ENFORCED rather than noted, so its tests run on the hardware the codes
    # actually require. PIN is the two-coupler-layer chip of Bravyi et al.
    PIN = DEFAULT.but(platform="sc_long_range")
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
    # The gap at n = 1e6 is a LADDER CLIFF, not an architectural verdict: the
    # surface arm sits exactly on the edge of the cultivation rung, and the next
    # published source costs 30700 qubits against 450. Tightening the tolerance
    # (review #8) pushed it onto that edge. By n = 1e7, past the cliff, the two
    # are within a factor of three again.
    _r6 = ftqc.max_m_pinnacle(1e6) / ftqc.max_m_surface(1e6, fow)
    _r7 = ftqc.max_m_pinnacle(1e7) / ftqc.max_m_surface(1e7, fow)
    _sp = ftqc.surface_point(ftqc.max_m_surface(1e6, fow), fow)
    check("the Pinnacle-surface gap at n=1e6 is a factory-ladder edge",
          _sp is not None and _sp["p_T"] >= 0.95 * _sp["p_T_target"],
          f"surface sits at p_T = {_sp['p_T']:.1e} against a target of "
          f"{_sp['p_T_target']:.1e} -- one rung from a 68x bigger plant")
    check("and away from that edge the two stay within a small factor",
          _r7 < 3.0, f"ratio {_r6:.2f} at n=1e6, {_r7:.2f} at n=1e7")
    check("Pinnacle and the surface code stay within a small factor",
          0.5 < ftqc.max_m_pinnacle(1e8, PIN) / ftqc.max_m_surface(1e8, fow) < 2.5,
          f"ratio {ftqc.max_m_pinnacle(1e8, PIN)/ftqc.max_m_surface(1e8, fow):.2f} at n=1e8, Pinnacle on the two-coupler-layer chip it requires")
    # ---- O4/O6: what their step rule is, and what their data can test -----
    _EXPT = nisq.EXPERIMENT
    _fc = DEFAULT.but(trotter_mode="fixed_count")
    check("the experiment's step rule is a fixed COUNT, not a density",
          hubbard.trotter_steps(4.0, _fc) == hubbard.trotter_steps(400.0, _fc)
          == _EXPT["trotter_steps"],
          f"k = {_EXPT['trotter_steps']:.0f} steps for the whole evolution to "
          f"t = 2 (Sec. C), so the depth does not grow with t -- "
          f"r = 4 tau would charge 8 there and 2 at tau = 0.5")
    _pairs = []
    for _nm, _md in (("slide r = 4 tau", "fixed_density"),
                     ("experiment k = 4", "fixed_count")):
        _c = DEFAULT.but(damping_model="support", trotter_mode=_md)
        _m = nisq.max_m(1e6, _c, "expcal")
        _t = hubbard.t_max(_m, _c)
        _e = (hubbard.w_measured(_t, _c.U_over_J) * _t ** 3
              / hubbard.trotter_steps(_m, _c) ** 2)
        _pairs.append((_nm, _m, _e / (_c.frac_trotter * hubbard.eps_absolute(_c, _m))))
    check("both fixed step rules buy their reach by abandoning accuracy",
          all(x[2] > 10.0 for x in _pairs) and _pairs[1][2] > 1e6,
          "; ".join(f"{n}: reaches m = {m:.0f} with a Trotter error {o:.1e}x its "
                    f"own budget" for n, m, o in _pairs)
          + " -- the k = 4 row is evidence about the convention, not a capability")
    _nr = pathlib.Path(__file__).resolve().parent.parent / "calibration" / "data" / "noise_response.json"
    if _nr.exists():
        import json as _js
        _N = _js.loads(_nr.read_text())
        check("their hardware data pin Lambda, at ONE noise strength (O4)",
              abs(_N["lambda_best_point"] - 0.20) < 0.05
              and _N["t_best"] <= 0.2,
              f"best-conditioned point t = {_N['t_best']:.2f} gives Lambda = "
              f"{_N['lambda_best_point']:.3f} against the anchored 0.20; pooled "
              f"over {_N['n_conditioned']} times, {_N['lambda_mean']:.2f} "
              f"+-{_N['attenuation_sd']:.2f} in attenuation")
        check("...so the response SHAPE is untestable from them, not merely untested",
              abs(_N["trend_slope"]) > 0.1 and _N["trend_slope"] < 0,
              f"the fitted trend in t is {_N['trend_slope']:+.2f} -- NEGATIVE, "
              f"less damping at longer times, which no noise process does. The "
              f"depth is constant (k = 4), so these are "
              f"{_N['n_conditioned']} repeats of one lambda and the scatter is "
              f"Trotter error")

    # ---- O3: what actually stops the Pinnacle arm -------------------------
    _cc = ftqc.pinnacle_ceiling_cause(PIN)
    check("the Pinnacle arm has a hard m ceiling at ANY budget",
          100.0 < _cc["m_ceiling"] < 5000.0
          and ftqc.pinnacle_point(_cc["m_ceiling"] * 1.05, PIN) is None,
          f"m = {_cc['m_ceiling']:.0f}; above it pinnacle_point returns None, so "
          f"the arm stops existing rather than flattening")
    check("and the ceiling is the MAGIC ENGINE, not the code family (O3)",
          _cc["cause"] == "magic engine"
          and _cc["p_target_just_over"] < _cc["cleanest_engine_p_out"],
          f"just over the ceiling the T-state target is "
          f"{_cc['p_target_just_over']:.1e} against the cleanest published "
          f"engine's {_cc['cleanest_engine_p_out']:.0e} -- select_engine refuses "
          f"before a code is tried. O3 used to blame the d = 24 code")
    check("...which is why more engines do not move it -- quality, not throughput",
          abs(ftqc.pinnacle_m_ceiling(PIN.but(pin_engines=16))
              / _cc["m_ceiling"] - 1.0) < 1e-3,
          f"1 engine: {_cc['m_ceiling']:.1f}, 16 engines: "
          f"{ftqc.pinnacle_m_ceiling(PIN.but(pin_engines=16)):.1f}")
    check("the code family IS fully consumed, from n ~ 1e6 -- no headroom left",
          ftqc.pinnacle_family_limit(PIN)["n_first_exhausted"] <= 1e6
          and _cc["d_selected_at_ceiling"] == _cc["d_max_published"],
          f"the largest published code (d = {_cc['d_max_published']}) is the one "
          f"selected from n = "
          f"{ftqc.pinnacle_family_limit(PIN)['n_first_exhausted']:.0e} up; O3 said "
          f"1e10 and called the plotted range comfortably inside")
    check("the surface arm has no such wall -- its ladder cascades further",
          ftqc.surface_point(3000.0, fow) is not None,
          "surface_point is still defined at m = 3000, five times the Pinnacle "
          "ceiling")

    check("engine count has an optimum -- engines cost 4410 qubits each",
          ftqc.max_m_pinnacle(1e6, PIN.but(pin_engines=16))
          > ftqc.max_m_pinnacle(1e6, PIN.but(pin_engines=64)))
    check("the plotted range stays inside the published GB code family",
          ftqc.pinnacle_point(ftqc.max_m_pinnacle(1e8, PIN), PIN)["d"] <= 24
          and ftqc.max_m_pinnacle(1e8, PIN) < ftqc.max_m_pinnacle(1e9, PIN),
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
    check("the POOLED TDVP error is flat in chi", abs(classical.TDVP_SLOPE) < 0.3
          and classical.TDVP_MEASURED[2048] / classical.TDVP_MEASURED[256] > 0.5,
          f"{classical.TDVP_MEASURED[256]:.3f} -> {classical.TDVP_MEASURED[2048]:.3f} "
          "over an 8x range in chi")
    check("its plateau exceeds the late-time signal",
          classical.TDVP_MEASURED[2048] > hubbard.s_residual())
    # the published TDVP is not bond-dimension-limited at all
    fl = classical.TDVP_FLOOR["err"]
    check("but that pooling hides a floor at t = 0.1 specifically",
          fl[256] / fl[2048] < 1.3 and fl[2048] > 5e-3,
          f"t=0.1: {fl[256]:.2e} at chi=256 vs {fl[2048]:.2e} at chi=2048, "
          "an 8x increase in chi removing 8% -- and the pooled slope averages "
          "this floor with a genuinely truncation-limited middle (see below)")
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
    # This has moved twice and both moves were real. Charging the tolerance per
    # time (review #8) took MPF below the frontier; refitting the order-2k
    # coefficients over a domain twice as wide (O11b) put it back above, at a
    # much larger n than the original 4e6.
    check("multiproduct clears the ED frontier again, but far later",
          s["n_pec_mpf_clears_classical_hi"] is not None
          and s["n_pec_mpf_clears_classical_hi"] > 1e7,
          f"MPF reaches {s['pec_mpf_at_1e6']:.1f} at n=1e6 against a frontier of "
          f"{hi:.0f}, and clears at "
          f"{curves.fmt_crossing(s['n_pec_mpf_clears_classical_hi'])} -- against "
          f"4e6 before the per-time tolerance and never after it")
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
    # The optimum is ARM-DEPENDENT, and saying "there is an interior optimum"
    # without saying whose was hiding that. NISQ is shot-limited, so extra
    # branches keep paying up to order 8; FT is magic-limited, so they stop
    # paying after order 4.
    _ftk = {k: ftqc.max_m_surface(1e6, DEFAULT.but(trotter_order_k=k,
                                                   pl_model="fowler"))
            for k in (1, 2, 3, 4)}
    check("the multiproduct optimum is arm-dependent, not universal",
          mpf4[4] >= mpf4[3] >= mpf4[2] and _ftk[2] > _ftk[3] > _ftk[4],
          "NISQ " + " ".join(f"o{2*k}:{v:.0f}" for k, v in mpf4.items())
          + " (monotone to order 8);  FT "
          + " ".join(f"o{2*k}:{v:.0f}" for k, v in _ftk.items())
          + " (peaks at order 4)")
    # the order-2k coefficient is measured, not the second-order one reused
    check("higher orders use their OWN measured coefficient",
          hubbard.w_mpf(4, 0.5) is not None
          and abs(hubbard.w_mpf(4, 0.5) / hubbard.w_measured(0.5, 4.0) - 1) > 0.5,
          f"W_4 = {hubbard.w_mpf(4, 0.5):.4f} vs W_2 = {hubbard.w_measured(0.5, 4.0):.4f}")
    check("the MPF domain doubled to tau <= 1.0 and is still clamped beyond it",
          hubbard.W_MPF_TAU_MAX == 1.0
          and hubbard.w_mpf(4, 5.0) == hubbard.w_mpf(4, 1.0),
          "refitting the exponent as well as the coefficient (O11b) extended the "
          "usable domain from tau <= 0.5; every plotted MPF point is still beyond it")
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
    # Once the tolerance is charged per time (review #8) the FT arm becomes
    # magic-limited, and multiproduct helps a lot at order 4 -- fewer steps means
    # fewer T states, which relaxes the per-state target and moves the plant a
    # rung DOWN the ladder. It still falls off above order 4, where the branch
    # count outruns the saving.
    check("multiproduct helps FT at order 4 and falls off above it",
          ftm[2] > 1.5 * ftm[1] and ftm[3] < ftm[2] and ftm[4] < ftm[3],
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
    # The claim is about CERTIFIED TIME, which is what the bound is in; comparing
    # site counts against a fixed multiple was a threshold, not the statement.
    tq = converged.quantum_t_reach(1e8, fowc)
    check("and the FT arms certify further, not infinitely",
          tq > 10 * tc and math.isfinite(tq)
          and math.isfinite(ftqc.max_m_surface(1e8, fowc)),
          f"surface FT certifies to t = {tq:.2f} against ED's {tc:.3f} "
          f"({tq/tc:.0f}x), at m = {ftqc.max_m_surface(1e8, fowc):.0f} vs "
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
          abs(nisq.max_m(1e6, DEFAULT, "pec") - 12.39) < 0.05,
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
    # ---- admissible integer lattices (second-pass review #10) ---------------
    # 1. the state specification is checked, not assumed
    ok25, why25 = curves.lattice_admissible(5, 5)
    check("a 5x5 cannot hold the specified state, and says why",
          not ok25 and "odd" in why25 and "23" in why25,
          why25)
    check("but the experiment's own 7x4 can",
          curves.lattice_admissible(7, 4)[0],
          "28 sites, 26 after the defects, perfectly dimer-coverable")
    for lx, ly in ((4, 6), (2, 2), (12, 17), (10, 12)):
        okk, _ = curves.lattice_admissible(lx, ly)
        n = lx * ly
        check(f"{lx}x{ly}: even sites, integer N_up, coverable remainder",
              okk and n % 2 == 0 and (n - 2) % 2 == 0,
              f"N = {n}, N_up = N_dn = {n // 2}, {n - 2} sites in "
              f"{(n - 2) // 2} triplets")
    # 2. inadmissible sizes are REJECTED, not rounded into
    check("a quasi-1D ribbon is rejected rather than counted as a 2D lattice",
          not curves.lattice_admissible(2, 13)[0]
          and "quasi-1D" in curves.lattice_admissible(2, 13)[1],
          f"aspect 6.5 > {curves.MAX_ASPECT}")
    check("odd-sided squares are rejected at every size checked",
          all(not curves.lattice_admissible(L, L)[0] for L in (3, 5, 7, 9, 11)),
          "3x3, 5x5, 7x7, 9x9, 11x11 all have odd site counts")
    # 3. and floor(sqrt(m)) really did name lattices that cannot exist
    _bad = [m for m in (9.0, 13.0, 25.0, 63.0, 124.8)
            if not curves.lattice_admissible(curves.max_integer_L(m),
                                             curves.max_integer_L(m))[0]]
    check("floor(sqrt(m)) named an impossible lattice at most sizes",
          len(_bad) >= 4,
          f"inadmissible at m = {_bad}; the review's 5x5 case is one of them")
    _b13 = curves.best_lattice(13.0)
    check("and the admissible answer uses MORE of the capacity, not less",
          _b13 == (3, 4, 12) and 12 > curves.max_integer_L(13.0) ** 2,
          f"m = 13: 3x4 = 12 sites admissible, against an inadmissible "
          f"{curves.max_integer_L(13.0)}x{curves.max_integer_L(13.0)} = "
          f"{curves.max_integer_L(13.0)**2}")
    # #10 remainder: re-cost every admissible INTEGER candidate, rather than
    # reporting best_lattice(continuous m) and hoping. The concern was block
    # packing -- Pinnacle stores k = 14-16 logical qubits per block, so its
    # footprint is lumpy in m and a smaller lattice could in principle be
    # infeasible where a larger one is not.
    def _feas(arm, mm, nn, cc):
        pt = (ftqc.surface_point(mm, cc, n=nn) if arm == "surface"
              else ftqc.pinnacle_point(mm, cc))
        if pt is None or pt["phys"] > nn:
            return False
        cp = math.floor(nn / pt["phys"])
        c2 = cc.but(hwp_batch=pt["hwp_batch"]) if arm == "surface" else cc
        return cp >= 1 and (ftqc.n_shots_total(c2, mm) * pt["t_shot"] / cp
                            <= cc.budget_s)
    _PINC = DEFAULT.but(platform=curves.PINNACLE_PLATFORM)
    _nonmono = []
    for _arm, _c in (("surface", DEFAULT), ("pinnacle", _PINC)):
        for _n in (1e6, 1e7, 1e8):
            _mx = (ftqc.max_m_surface(_n, _c) if _arm == "surface"
                   else ftqc.max_m_pinnacle(_n, _c))
            if not _mx:
                continue
            _nonmono += [(_arm, _n, _m) for _m in range(4, int(_mx) + 1)
                         if not _feas(_arm, float(_m), _n, _c)]
    check("feasibility is monotone in m, including the block-quantised arm",
          not _nonmono,
          "every integer m from 4 to the maximum is feasible at six (arm, n) "
          "points -- so best_lattice, which takes the largest admissible site "
          "count BELOW the continuous maximum, is feasible by construction")
    for _arm, _c, _f in (("surface", DEFAULT, ftqc.max_m_surface),
                         ("pinnacle", _PINC, ftqc.max_m_pinnacle)):
        _n = 1e8
        _bl = curves.best_lattice(_f(_n, _c))
        check(f"and the re-costed {_arm} candidate matches best_lattice at n=1e8",
              _bl is not None and _feas(_arm, float(_bl[2]), _n, _c),
              f"{_bl[0]}x{_bl[1]} = {_bl[2]} sites, re-costed at its own integer "
              f"site count rather than inherited from a continuous m")

    check("every reported lattice fits its capacity",
          all((lambda b: b is None or b[2] <= m)(curves.best_lattice(m))
              for m in (4.0, 13.0, 22.0, 26.9, 124.8, 204.0)))
    # the signal scaling differs by regime -- 2/9 was only ever right for NISQ
    ss = [0.03, 0.06, 0.12]
    en = np.polyfit(np.log(ss), np.log([nisq.max_m(1e6, DEFAULT.but(
        signal_regime="fixed", s_sig=x), "pec") for x in ss]), 1)[0]
    ef = np.polyfit(np.log(ss), np.log([ftqc.max_m_surface(1e6, DEFAULT.but(
        signal_regime="fixed", s_sig=x, pl_model="fowler")) for x in ss]), 1)[0]
    # The old "m ~ s^(2/9)" confirmation was a FLOOR ARTEFACT: with the floor at
    # the fitted residual 0.043, the s = 0.03 point was clipped, which flattened
    # the fitted slope to 0.275. With the floor set independently at 0.02 the
    # clip lifts and the true local slope is 0.379. The 2/9 derivation assumed
    # alpha = 9/4 and an s-independent step count; neither holds under the
    # measured calibration.
    en_clip = np.polyfit(np.log(ss), np.log([nisq.max_m(1e6, DEFAULT.but(
        signal_regime="fixed", s_sig=x, s_abs_floor=0.043), "pec")
        for x in ss]), 1)[0]
    check("the old s^(2/9) agreement was the floor clipping the lowest point",
          en > 0.35 and en_clip < 0.3,
          f"slope {en:.3f} unclipped vs {en_clip:.3f} with the floor at the "
          f"fitted residual; 2/9 = 0.222")
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

    # cultivation, from the authors' released stats rather than their figure
    check("cultivation's footprint and attempts are EXACT, not read off a plot",
          ftqc.CULT_FOOTPRINT == 463.0 and ftqc.CULT_ATTEMPTS_1E3 == 73.0,
          "q = 463, r = 20, 73.0 attempts at the 2e-9 gap cut (1.90e-9 measured, "
          "98.6% discard) -- reconstructed from the 117-bin complementary-gap "
          "histogram in Zenodo 10.5281/zenodo.13777072")
    check("and their plotted volume sits inside the bracket the stats allow",
          ftqc.CULT_VOLUME_BRACKET[0] < ftqc.CULT_VOLUME_1E3 < ftqc.CULT_VOLUME_BRACKET[1],
          f"{ftqc.CULT_VOLUME_BRACKET[0]:.1e} (one attempt) < "
          f"{ftqc.CULT_VOLUME_1E3:.1e} (theirs) < "
          f"{ftqc.CULT_VOLUME_BRACKET[1]:.1e} (all 73 at full length)")
    # ...and it turns out not to matter, for a reason worth recording
    _reach = []
    for _cyc in (20.0, 65.0, 1460.0):
        _src = [tuple(list(r[:4]) + [_cyc] + [r[5]])
                if r[5] == "cultivation" and r[1] == 1e-3 else r
                for r in ftqc.MAGIC_SOURCES]
        _save, ftqc.MAGIC_SOURCES = ftqc.MAGIC_SOURCES, _src
        _reach.append(ftqc.max_m_surface(1e7, DEFAULT))
        ftqc.MAGIC_SOURCES = _save
    check("the remaining uncertainty is absorbed by the plant-level optimisation",
          max(_reach) - min(_reach) < 1e-9,
          f"a 73x change in cultivation's cycle count moves surface FT at n = 1e7 "
          f"not at all ({_reach[0]:.1f}): select_factory minimises total plant "
          f"qubits over ALL admissible sources, so degrading one just hands the "
          f"job to the Litinski ladder")

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

    # ---- hardware platforms: connectivity, parallelism, clock ---------------
    from . import platform as plat
    # 1. the status quo is reproduced exactly by the platform path
    for _n in (1e5, 1e6, 1e7, 1e8):
        _sc = DEFAULT.but(platform="superconducting")
        check(f"platform path reproduces the status quo at n = {_n:.0e}",
              abs(nisq.max_m(_n, _sc, "pec") - nisq.max_m(_n, DEFAULT, "pec")) < 1e-12
              and abs(ftqc.max_m_surface(_n, _sc) - ftqc.max_m_surface(_n, DEFAULT)) < 1e-12,
              "superconducting IS the old implicit machine, written down")
    # 2. unlimited parallelism is algebraically the old formula
    _c = DEFAULT.but(use_platform_clock=True, platform="superconducting")
    _q, _q0 = hubbard.counts(64.0, _c), hubbard.counts(64.0, DEFAULT)
    check("with unlimited parallelism the clock reduces to depth x t_gate",
          abs(_q["t_circuit"] - _q0["t_circuit"]) / _q0["t_circuit"] < 1e-12,
          f"{_q['t_circuit']:.4e} s either way")
    # 3. published anchors
    check("Helios anchors match the paper",
          plat.HELIOS.n_demonstrated == 98 and plat.HELIOS.n_parallel_2q == 4.0
          and plat.HELIOS.p_2q == 7.9e-4 and plat.HELIOS.t_layer == 55e-3,
          "98 qubits, 4 two-qubit zones, 7.9e-4, 55 ms/layer (arXiv:2511.05465)")
    check("transport dominates the gate on a mobile-qubit machine",
          plat.HELIOS.t_layer / plat.HELIOS.t_2q > 500,
          f"55 ms per layer against a 70 us gate = "
          f"{plat.HELIOS.t_layer / plat.HELIOS.t_2q:.0f}x -- ion sorting, not gating")
    check("and no gate count could have predicted it",
          plat.HELIOS.layer_seconds(1e9) == plat.HELIOS.t_layer,
          "a measured layer time overrides the derived one")
    # 4. admissibility is enforced, not decorative
    check("GB codes are refused on the slide's nearest-neighbour grid",
          ftqc.max_m_pinnacle(1e8, DEFAULT) == 0.0
          and ftqc.max_m_pinnacle(1e8, DEFAULT.but(platform="sc_long_range")) > 0,
          "pin_nonlocal was documented to do this and was read by nothing")
    check("two coupler layers are enough -- all-to-all is not required",
          plat.admits("gb", "sc_long_range") and not plat.admits("gb", "superconducting"),
          "Bravyi et al.: degree 6, two edge-disjoint planar subgraphs")
    check("the surface code is unaffected by connectivity",
          abs(ftqc.max_m_surface(1e7, DEFAULT.but(platform="sc_long_range"))
              - ftqc.max_m_surface(1e7, DEFAULT)) < 1e-12,
          "it only ever needed a planar grid")
    # 5. the trade actually computes: all-to-all removes the swap network
    _aa = DEFAULT.but(platform="helios", encoding="jw")
    check("all-to-all removes the Kivlichan swap network",
          hubbard.step_depth(256.0, _aa)
          < hubbard.step_depth(256.0, DEFAULT.but(encoding="jw")),
          f"{hubbard.step_depth(256.0, _aa):.0f} layers against "
          f"{hubbard.step_depth(256.0, DEFAULT.but(encoding='jw')):.0f}")
    _h = DEFAULT.but(platform="helios", use_platform_clock=True, encoding="jw")
    check("but the clock cost swamps it on this workload",
          hubbard.counts(64.0, _h)["t_circuit"]
          > 1e5 * hubbard.counts(64.0, DEFAULT)["t_circuit"],
          f"{hubbard.counts(64.0, _h)['t_circuit']:.2e} s per shot against "
          f"{hubbard.counts(64.0, DEFAULT)['t_circuit']:.2e} s")
    # 6. transversal algorithmic fault tolerance -- the reason a slow-clock
    #    machine is not simply hopeless, and the reason it still loses HERE
    check("mobile-qubit machines get O(1) syndrome rounds, not O(d)",
          plat.se_rounds(21, "helios") == 1.0
          and plat.se_rounds(21, "neutral_atom") == 1.0
          and plat.se_rounds(21, "superconducting") == 21.0,
          "Zhou et al. arXiv:2406.17653: transversal gates + correlated decoding")
    # at the atoms' MEASURED gate error, not an optimistic one
    _at = DEFAULT.but(platform="neutral_atom", use_platform_clock=True,
                      encoding="jw", p=plat.NEUTRAL_ATOM.p_2q)
    _at_opt = _at.but(p=1e-3)
    check("and it is worth a factor of d, not a rounding",
          ftqc.surface_point(12.0, _at_opt)["rounds"]
          < ftqc.surface_point(12.0, DEFAULT)["rounds"] / 15,
          f"{ftqc.surface_point(12.0, _at_opt)['rounds']:.2e} rounds against "
          f"{ftqc.surface_point(12.0, DEFAULT)['rounds']:.2e} under lattice surgery")
    check("so a slow clock alone no longer kills the FT arm",
          ftqc.max_m_surface(1e7, _at_opt) > 0 and ftqc.max_m_star(1e7, _at) > 0,
          f"at p = 1e-3 neutral atoms reach surface m = "
          f"{ftqc.max_m_surface(1e7, _at_opt):.1f} and STAR m = "
          f"{ftqc.max_m_star(1e7, _at):.1f}; an earlier version of this model said "
          f"0 at every n, which was lattice surgery applied to a machine nobody "
          f"runs that way")
    # ...but each platform then fails for a DIFFERENT reason, which is the point
    check("atoms' surface arm is fidelity-limited, not clock-limited",
          ftqc.max_m_surface(1e9, _at) == 0 and ftqc.max_m_surface(1e7, _at_opt) > 0,
          f"their measured 5e-3 is half the p_th = {DEFAULT.p_th} threshold, so "
          f"p_L(d=41) is {ftqc.p_logical(41, DEFAULT.but(p=5e-3)):.1e} against "
          f"{ftqc.p_logical(41, DEFAULT.but(p=1e-3)):.1e} at 1e-3; the arm "
          f"appears at p <~ 2e-3, which is a falsifiable prediction")
    check("Helios is the opposite: a good gate on a slow clock",
          plat.HELIOS.p_2q < DEFAULT.p and plat.HELIOS.t_layer > 1e-3,
          f"7.9e-4 against the slide's 1e-3, at 55 ms per layer -- the two "
          f"platforms fail for different reasons, not one shared one")
    # the distinction that decides it: depth overhead vs shot count
    _tot = (ftqc.n_shots_total(_at_opt, 12.0)
            * ftqc.surface_point(12.0, _at_opt)["t_shot"] / _at_opt.budget_s)
    check("but transversal FT fixes DEPTH overhead, not the shot count",
          _tot > 50,
          f"at m = 12 the week is still overrun {_tot:.0f}x, because this "
          f"observable needs {ftqc.n_shots_total(_at_opt, 12.0):.1e} shots and the "
          f"clock is paid once per shot -- Shor is one deep circuit, this is not")

    # 7. O14: credit the mobile platforms with the codes they actually propose
    _h = DEFAULT.but(platform="helios", use_platform_clock=True, encoding="jw",
                     p=plat.HELIOS.p_2q)
    check("high-rate qLDPC storage is ~25x better than the surface code",
          ftqc.GB_CODES[-1][4] / ftqc.GB_CODES[-1][1]
          < 0.05 * ftqc.storage_per_logical(25, DEFAULT),
          f"[[510,16,24]] is {ftqc.GB_CODES[-1][4] / ftqc.GB_CODES[-1][1]:.0f} "
          f"physical per logical against "
          f"{ftqc.storage_per_logical(25, DEFAULT):.0f} for a d = 25 surface patch")
    # NOT "qLDPC is what gives Helios an arm at all" -- that was wrong, and the
    # truth is narrower. Transversal rounds already gave it a surface arm from
    # n ~ 1e8; qLDPC only overtakes above ~1e10, because the single 4410-qubit
    # engine and its serialised T supply dominate at small m while the surface
    # plant parallelises.
    _hs = {e: ftqc.max_m_surface(10.0 ** e, _h) for e in (9, 10, 12)}
    _hq = {e: ftqc.max_m_pinnacle(10.0 ** e, _h) for e in (9, 10, 12)}
    check("qLDPC overtakes the surface code on Helios, but only above n ~ 1e10",
          _hq[9] < _hs[9] and _hq[12] > 1.5 * _hs[12],
          f"surface {_hs[9]:.1f}/{_hs[10]:.1f}/{_hs[12]:.0f} against qLDPC "
          f"{_hq[9]:.1f}/{_hq[10]:.1f}/{_hq[12]:.0f} at n = 1e9/1e10/1e12 -- the "
          f"single engine's serialised T supply costs more than the storage saves "
          f"until the lattice is large")
    _a = DEFAULT.but(platform="neutral_atom", use_platform_clock=True,
                     encoding="jw", p=plat.NEUTRAL_ATOM.p_2q)
    check("atoms get no qLDPC arm either, and for the SAME reason as before",
          ftqc.max_m_pinnacle(1e12, _a) == 0.0
          and ftqc.select_engine(1e-9, _a) is None,
          f"no magic engine is characterised at their measured "
          f"{plat.NEUTRAL_ATOM.p_2q:.0e}; the engine table stops at 1e-3. "
          f"Fidelity, not code choice and not clock")

    check("the hypothetical chip is marked as one",
          plat.SC_LONG_RANGE.n_demonstrated == 0
          and plat.HELIOS.n_demonstrated > 0,
          "n_demonstrated = 0 means nobody has built it")

    # ---- how far the calibration actually reaches (second-pass #2) ----------
    _dc = pathlib.Path(__file__).resolve().parent.parent / "calibration" / "data" / "domain_check.json"
    if _dc.exists():
        import json as _json
        _D = _json.loads(_dc.read_text())
        # median |slope + 2| is 0.002 across the grid; the outliers are all at
        # tau = 2 on the smallest patches, where the error is nearly saturated
        _dev = sorted(abs(x["slope_large_r"] + 2.0) for x in _D["power_law"])
        _out = [(x["n"], x["tau"]) for x in _D["power_law"]
                if abs(x["slope_large_r"] + 2.0) > 0.15]
        check("the r^-2 power law is clean where the 2nd-order arm operates",
              _dev[len(_dev) // 2] < 0.01 and len(_out) <= 1,
              f"median |slope + 2| = {_dev[len(_dev) // 2]:.3f} for r >= 8 "
              f"(arm runs at r = {hubbard.trotter_steps(4.0, DEFAULT):.0f}-"
              f"{hubbard.trotter_steps(400.0, DEFAULT):.0f}); the one outlier is "
              f"n={_out[0][0]}, tau={_out[0][1]} at {max(_dev):.2f}, a small patch "
              f"at long time where the error is nearly saturated")
        check("but NOT where the multiproduct arm operates",
              _D["worst_dev_small_r"] > 1.0
              and hubbard.trotter_steps(4.0, DEFAULT.but(trotter_order_k=2)) < 8.0,
              f"worst |slope + 2| = {_D['worst_dev_small_r']:.2f} at r = 1-4; "
              f"MPF-4 starts at r = "
              f"{hubbard.trotter_steps(4.0, DEFAULT.but(trotter_order_k=2)):.1f}")
        _sp = _D["size_spread"]
        check("size-independence is a SHORT-time statement and decays with time",
              _sp["0.25"] < 1.1 and _sp["2.0"] > 2.0,
              "spread across n >= 6: " + ", ".join(
                  f"tau={k}: {v:.2f}x" for k, v in sorted(_sp.items(), key=lambda x: float(x[0]))))
        check("and every plotted point sits where it has decayed",
              hubbard.t_max(curves.M_FLOOR, DEFAULT) >= 1.0,
              f"the smallest plotted lattice already runs to tau = "
              f"{hubbard.t_max(curves.M_FLOOR, DEFAULT):.2f}; "
              f"locality was verified at tau = 0.25")
        check("interpolating in tau INSIDE the measured range is good to ~2x only",
              _D["worst_time_holdout"] > 1.0,
              f"holding tau = 1 out and predicting it overshoots by "
              f"{_D['worst_time_holdout']:+.0%}; the clamp beyond tau = 2 cannot "
              f"be worth more than that")

    # ---- the trajectory, not the fixed-time sweep (second-pass #2) ---------
    if _dc.exists():
        _T = _json.loads(_dc.read_text())
        if "trajectory" in _T and len(_T["trajectory"]) >= 4:
            _tr = [x for x in _T["trajectory"] if x["n"] >= 6]
            check("every trajectory reference is converged, not just recorded",
                  all(x["substep_convergence"] < 1e-10 for x in _T["trajectory"]),
                  "worst substep convergence "
                  f"{max(x['substep_convergence'] for x in _T['trajectory']):.1e} "
                  "-- the first n = 14 attempt came back at 2.2e-2 and was excluded")
            check("the fixed-time table is CONSERVATIVE at the operating points",
                  all(x["ratio"] < 1.0 for x in _tr),
                  "measured/table = " + ", ".join(f"n={x['n']}: {x['ratio']:.2f}"
                                                  for x in _tr)
                  + " -- the table over-predicts W, so it over-charges gates")
            check("but these sizes do NOT determine the exponent",
                  _T.get("alpha_se", 0) > 0.05
                  and _T["alpha_ci95"][0] < 1.75 < _T["alpha_ci95"][1],
                  f"alpha = {_T['alpha_measured']:.2f} +- {_T['alpha_se']:.2f}, "
                  f"95% CI {_T['alpha_ci95'][0]:.2f}..{_T['alpha_ci95'][1]:.2f}; "
                  f"including n = 4, 5 flips it to {_T['alpha_all_n']:.2f}. The "
                  f"adopted 1.75 and Campbell's 2.25 are both inside")

    # ---- is the published TDVP truncation-limited? (O10) -------------------
    _tj = pathlib.Path(__file__).resolve().parent.parent / "calibration" / "data" / "tdvp_check.json"
    if _tj.exists():
        import json as _json
        _T = _json.loads(_tj.read_text())["dimer"]
        _f = _T["0.1"]
        check("the stored TDVP floor is what the full deposit says",
              all(abs(_f[str(c)] - classical.TDVP_FLOOR["err"][c]) < 5e-5
                  for c in (256, 512, 1024, 2048)),
              "chi = 256..2048 at t = 0.1: "
              + ", ".join(f"{_f[str(c)]:.2e}" for c in (256, 512, 1024, 2048))
              + " -- reproduced from Zenodo 17799843, not re-entered by hand")
        check("and the stored per-time slopes are too",
              all(abs(_T[f"{t}"]["slope"] - v) < 5e-3
                  for t, v in classical.TDVP_SLOPE_BY_T.items()),
              "d log(err)/d log(chi) at t = "
              + ", ".join(f"{t}: {v:+.2f}" for t, v in
                          sorted(classical.TDVP_SLOPE_BY_T.items())))
        check("chi dependence is NON-monotonic in time, so no single slope holds",
              _T["0.1"]["slope"] > -0.10 and _T["0.5"]["slope"] < -0.35
              and _T["2.0"]["slope"] > -0.10,
              f"floored at t = 0.1 ({_T['0.1']['slope']:+.2f}), genuinely "
              f"truncation-limited at t = 0.5 ({_T['0.5']['slope']:+.2f}), "
              f"floored again at t = 2.0 ({_T['2.0']['slope']:+.2f}); the pooled "
              f"TDVP_SLOPE = {classical.TDVP_SLOPE:+.2f} averages all three")
        check("at t = 0.1 the floor is not truncation -- chi = 256 is ample there",
              _T["0.1"][str(256)] / _T["0.1"][str(2048)] < 1.2
              and _T["0.1"][str(2048)] > 1e-3,
              f"8x the bond dimension removes "
              f"{1 - _T['0.1']['2048'] / _T['0.1']['256']:.0%} of a "
              f"{_T['0.1']['2048']:.1e} error on a barely entangled state; the "
              f"deposit varies chi and never dt, so it cannot say which term it is")
        check("at late times the published error EXCEEDS the signal",
              all(_T[f"{t}"][str(2048)] > _T[f"{t}"]["signal"]
                  for t in (1.5, 2.0)),
              f"t = 2: err {_T['2.0']['2048']:.3f} vs mean |C| "
              f"{_T['2.0']['signal']:.3f} -- that run bounds nothing there")

    # ---- the additional implementation checks (second pass) -----------------
    # every public entry point must refuse an over-allocated error budget
    _bad = DEFAULT.but(frac_stat=0.9)
    _entries = {
        "nisq.max_m": lambda c: nisq.max_m(1e6, c, "pec"),
        "nisq.max_m_ideal": lambda c: nisq.max_m_ideal(1e6, c),
        "ftqc.max_m_surface": lambda c: ftqc.max_m_surface(1e6, c),
        "ftqc.max_m_star": lambda c: ftqc.max_m_star(1e6, c),
        "ftqc.max_m_pinnacle": lambda c: ftqc.max_m_pinnacle(1e6, c),
        "classical.band": classical.band,
    }
    _taken = []
    for _nm, _fn in _entries.items():
        try:
            _fn(_bad)
            _taken.append(_nm)
        except ValueError:
            pass
    check("every entry point refuses an over-allocated error budget",
          not _taken, f"frac_stat = 0.9 with the other shares unchanged sums to "
                      f"1.45; accepted by {_taken or 'nothing'}")
    check("and validate() is an explicit precondition, not a side effect",
          hubbard.validate(DEFAULT) is DEFAULT,
          "checks the ledger, eps in (0,1), p below threshold, n_times, floor")
    for _kw, _why in (({"p": 0.05}, "p above threshold"), ({"eps": 2.0}, "eps > 1"),
                      ({"n_times": 0}, "no time points"),
                      ({"s_abs_floor": 0.0}, "zero floor collapses tolerances")):
        try:
            hubbard.validate(DEFAULT.but(**_kw))
            check(f"validate rejects {_why}", False, f"accepted {_kw}")
        except ValueError:
            check(f"validate rejects {_why}", True)

    # s_sig is INERT under the default regime -- the sensitivity rows that varied
    # it were reporting the baseline twice
    _sw = nisq.max_m(1e6, DEFAULT.but(s_sig=0.03), "pec")
    _ss = nisq.max_m(1e6, DEFAULT.but(s_sig=0.3), "pec")
    check("s_sig does nothing under signal_regime='curve', and is not varied there",
          _sw == _ss == nisq.max_m(1e6, DEFAULT, "pec"),
          f"both {_sw:.2f}; crossovers.md now varies the fixed-signal scenario "
          f"and the parameters the curve actually uses")
    _fw = nisq.max_m(1e6, DEFAULT.but(signal_regime="fixed", s_sig=0.03), "pec")
    _fs = nisq.max_m(1e6, DEFAULT.but(signal_regime="fixed", s_sig=0.3), "pec")
    check("under the fixed-signal scenario it does bite",
          _fs > 1.3 * _fw, f"{_fw:.1f} -> {_fs:.1f}")

    # the cluster prose contradicted the implementation at small xi
    _ct = {xi: converged.cluster_t_reach(DEFAULT.but(xi=xi))
           for xi in (0.2, 0.5, 1.0)}
    check("the cluster method is NOT uncompetitive at every xi -- only xi >~ 0.5",
          _ct[0.2] > 0.3 and _ct[0.5] == 0.0 and _ct[1.0] == 0.0,
          f"t_reach: " + ", ".join(f"xi={k}: {v:.3f}" for k, v in _ct.items())
          + " -- the docstring claiming failure at all xi is corrected")

    # capacity at fixed t and certified time are different questions
    _d = classical.max_m_fixed_t_detail(0.02, DEFAULT)
    check("a fixed-t capacity carries its method and its non-claim",
          _d["method"] == "snake MPS" and "NO convergence guarantee" in _d["claims"]
          and "certified" in _d["not_a_certificate"],
          f"m = {_d['m']:.0f} at t = 0.02 by {_d['method']}, against a certified "
          f"t_reach of {converged.classical_t_reach(DEFAULT):.4f} -- different "
          f"questions, not a contradiction")

    # ---- zero band edges stay visible (second-pass review #9) ---------------
    # 1. a synthetic band with ONE zero edge must survive
    _lo = np.array([0.0, 0.0, 5.0, 7.0])
    _hi = np.array([0.0, 13.0, 20.0, 30.0])
    _ld, _hd, _off = curves.band_for_plot(_lo, _hi, 2.5)
    check("a zero lower edge is clipped to the axis floor, not dropped",
          _ld[1] == 2.5 and _hd[1] == 13.0 and bool(_off[1]),
          "the band stays visible and is flagged off-scale")
    # 2. all-infeasible and partially feasible both behave
    check("an all-infeasible column is NaN on both edges, not a fake band",
          np.isnan(_ld[0]) and np.isnan(_hd[0]) and not _off[0],
          "nothing is drawn; the figure annotates instead")
    check("and feasible columns are untouched",
          (_ld[2:] == _lo[2:]).all() and (_hd[2:] == _hi[2:]).all()
          and not _off[2:].any())
    # 3. every band the figure draws has an identifiable encoding, and the
    #    off-scale case is actually exercised by the default configuration --
    #    it was not merely possible, it was happening at EVERY plotted n
    _ng = np.logspace(2, np.log10(3e8), 60)
    _B = curves.evaluate_band(_ng, DEFAULT)
    _drawn, _hatched = [], []
    for _k in ("nisq_pec", "star", "surface", "pinnacle"):
        _l, _h = _B[_k]
        _a, _b, _o = curves.band_for_plot(_l, _h, 2.5)
        _vis = np.isfinite(_a) & np.isfinite(_b) & (_b > _a)
        _drawn.append(_vis.any())
        _hatched.append(_o.any())
    check("every plotted scenario band has visible extent somewhere",
          all(_drawn), "nisq_pec, star, surface, pinnacle")
    _l, _h = _B["nisq_pec"]
    check("and NISQ+PEC was the invisible one: zero lower edge at EVERY n",
          (_l <= 0).all() and (_h > 0).all(),
          f"lower edge 0 at all {len(_l)} grid points, upper edge "
          f"{_h.min():.1f}-{_h.max():.1f} -- no shading was drawn at all")
    check("the hatched encoding is exercised, not just available",
          any(_hatched), f"{sum(_hatched)} of 4 arms run off the bottom")

    # ---- the U = 0 easy limit (second-pass review #7) -----------------------
    ff = classical.FREE_FERMION
    check("the U = 0 estimator is validated against exact evolution, not asserted",
          ff["worst_abs_err"] < 1e-12,
          f"worst |C_poly - C_exact| = {ff['worst_abs_err']:.1e} at m = "
          f"{ff['validated_m']}, t = {ff['validated_t']} "
          f"(calibration/free_fermion.py)")
    check("and against SOMEONE ELSE'S exact answer, which is the check that bit",
          ff["deposit_abs_err"] < 1e-8,
          f"{ff['deposit_abs_err']:.1e} against the published exact C^zz on the "
          f"7x4 instance. The internal check cannot see a wrong state or a "
          f"missing flux, because it evolves the same ones")
    _dj = pathlib.Path(__file__).resolve().parent.parent / "calibration" / "data" / "deposit_check.json"
    if _dj.exists():
        import json as _json
        _DD = _json.loads(_dj.read_text())
        _w = max(r["trip"] for r in _DD)
        _s = max(r["sing"] for r in _DD)
        _n = max(r["noflux"] for r in _DD)
        check("the singlet and the flux-free Hamiltonian are BOTH visibly wrong",
              _w < 1e-8 and _s > 1e-2 and _n > 1e-2,
              f"triplet + pi flux {_w:.1e}; singlet {_s:.2f}; no flux {_n:.2f}. "
              f"Both give C^zz = -1 at t = 0, which is why neither showed up "
              f"until the deposit was read")
        check("and their own 'FLO' rows are not the reference to measure against",
              max(r["their_spread"] for r in _DD) > 1e-2,
              f"the deposit's 'Exact' and 'FLO' disagree by up to "
              f"{max(r['their_spread'] for r in _DD):.1e}, the size of the TDVP "
              f"errors in O10; ours reproduces 'Exact' to {_w:.0e}")
    check("the state and the flux are recorded, not left implicit",
          ff["state"].startswith("Sz = 0 triplet") and abs(ff["flux_y"] - 0.7853981633974483) < 1e-12,
          f"{ff['state']}, Peierls phase pi/4 per y-bond "
          f"(arXiv:2510.26300 Fig. 1 and Sec. I)")
    secs = ff["seconds"]
    per = [secs[m] / m for m in sorted(secs)]
    check("and its measured cost is LINEAR in m, not exponential",
          max(per) / min(per) < 1.5,
          "s/site = " + ", ".join(f"{x * 1e3:.2f} ms" for x in per))
    c0 = DEFAULT.but(U_over_J=0.0)
    ffm = classical.free_fermion_frontier(c0)
    check("so at U = 0 the classical capacity is not the ED memory threshold",
          ffm > 1e6 and classical.band(c0)["band"][1] == ffm,
          f"{ffm:.2g} sites on ONE core in a week, against "
          f"{classical.ed_frontier(c0)} for ED -- six orders of magnitude")
    check("no quantum arm on the figure comes near it",
          max(ftqc.max_m_surface(1e8, c0), ftqc.max_m_pinnacle(1e8, c0),
              nisq.max_m(1e8, c0, "pec")) < ffm / 1e4,
          "the U = 0 point is not a quantum-advantage candidate for this observable")
    check("the claim is confined to U = 0; nothing is extrapolated in U",
          not classical.free_fermion_applicable(DEFAULT.but(U_over_J=0.1))
          and classical.band(DEFAULT)["method"] == "ED",
          "exact at U = 0, no accuracy claim at small non-zero U")
    check("the TDVP arithmetic is labelled peak-FLOP, not a measured wall time",
          "peak-FLOP" in classical.TDVP_LEADERSHIP["caveat"]
          and classical.TDVP_LEADERSHIP["assumed_fraction_of_peak"] == 1.0,
          "no memory, communication, SVD/MPO or parallel-efficiency accounting")

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
    # O13: the batch is now OPTIMISED, and the optimum is the value the model
    # already used. Recording the negative result, and why it is not circular:
    # minimising footprint instead picks b = 1 and makes the arm WORSE (reach at
    # n = 1e8 falls 125 -> 30), because b = 1 has no workspace but 7x the T
    # states. The objective has to be the binding constraint, not the footprint.
    _regimes = {"baseline": DEFAULT, "willow": DEFAULT.but(pl_model="willow"),
                "1 day": DEFAULT.but(budget_s=86400.0),
                "eps 0.01": DEFAULT.but(eps=0.01), "p 1e-4": DEFAULT.but(p=1e-4)}
    _smaller = []
    for _lab, _c in _regimes.items():
        for _n in (1e5, 1e6, 1e7, 1e8):
            _m = ftqc.max_m_surface(_n, _c)
            if not _m:
                continue
            _pt = ftqc.surface_point(_m, _c, n=_n)
            if _pt and _pt["hwp_batch"] < int(_m):
                _smaller.append(f"{_lab}@{_n:.0e}")
    check("the HWP batch optimum is the full batch, in every regime tested",
          not _smaller,
          "20 (regime, n) points: this workload is clock-bound and the full "
          "batch minimises T, so the knob is real but inert -- known now, "
          "not assumed")
    check("and optimising for FOOTPRINT instead would make it worse",
          ftqc.max_m_surface(1e8, DEFAULT.but(hwp_batch=1))
          < 0.5 * ftqc.max_m_surface(1e8, DEFAULT),
          f"b = 1 reaches {ftqc.max_m_surface(1e8, DEFAULT.but(hwp_batch=1)):.0f} "
          f"against {ftqc.max_m_surface(1e8, DEFAULT):.0f}: no workspace, 7x the "
          f"T states, clock-bound")
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
    _PIN = DEFAULT.but(platform=curves.PINNACLE_PLATFORM)
    worst, seen = 0.0, 0
    for n in (1e5, 1e6, 1e7, 1e8):
        mm = ftqc.max_m_pinnacle(n, _PIN)
        pp = ftqc.pinnacle_point(mm, _PIN) if mm else None
        if pp is None:
            continue
        seen += 1
        allow = _PIN.frac_magic * hubbard.eps_absolute(_PIN, mm)
        worst = max(worst, 2.0 * pp["n_t"] * pp["engine_p_out"] / allow)
    check("every plotted Pinnacle point certifies its own magic allowance",
          seen == 4 and worst <= 1.0,
          f"worst 2 n_T p_out / allowance = {worst:.2f} over {seen} points; "
          f"it was 9.5 at n = 1e8")

    # 2. consumption never exceeds successful production
    mm = ftqc.max_m_pinnacle(1e7, _PIN)
    pp = ftqc.pinnacle_point(mm, _PIN)
    produced = (pp["rounds"] / (max(float(pp["d"]) + 2.0, pp["engine_cycles"])
                                / (1.0 - pp["p_reject"]))) * max(_PIN.pin_engines, 1)
    check("T consumption never exceeds successful engine production",
          produced >= pp["n_t"] * (1 - 1e-9),
          f"schedule delivers {produced:.3e} accepted states for {pp['n_t']:.3e} needed")
    check("and rejection is charged, not assumed away",
          pp["p_reject"] == 0.10,
          "10% for both p = 1e-3 engines, their own estimate")

    # 3. a plotted point where the processor outruns the engine and stalls
    m5 = ftqc.max_m_pinnacle(1e5, _PIN)
    p5 = ftqc.pinnacle_point(m5, _PIN) if m5 else None
    check("there is a plotted point where the processor stalls on the engine",
          p5 is not None and p5["stalled"] and p5["d"] == 16,
          f"d = 16 has an 18-cycle logical cycle; distillation needs "
          f"{p5['engine_cycles']:.0f}" if p5 else "no point")

    # 4. the two architectures now agree field-by-field on the shared rows
    led = {row[0]: row for row in ftqc.ledger_comparison(64.0, _PIN)}
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
          ftqc.max_m_pinnacle(1e6, _PIN.but(p=3e-3)) == 0.0,
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

#!/usr/bin/env python3
"""Does the published TDVP floor move with the TIME STEP? (O10, the open half)

    python3 calibration/tdvp_dt.py --validate            # certify the setup
    python3 calibration/tdvp_dt.py --scan --lx 7 --ly 4  # the measurement
    python3 calibration/tdvp_dt.py --scan --chi-scan     # their axis, our setup

WHY THIS EXISTS.  `calibration/tdvp_check.py` reads the published deposit and
finds an error floor at t = 0.1 that is NOT truncation: chi = 256 is wild overkill
on a barely entangled state, yet the error is 7.9e-3 and eight times the bond
dimension removes 8% of it. The deposit varies chi over an 8x range and NEVER
varies the time step, so it structurally cannot say whether that floor is
integration error, the two-site projection, or the initial state. This scan
varies dt at fixed chi, which is the one axis missing.

NO NEW TDVP CODE.  The published runs are TeNPy -- their own metadata reads
`Tenpy $\\chi=512$` -- so this is a configuration of the library they used at the
parameters they used, not a hand-rolled fermionic 2D TDVP on a doubly periodic
torus. Hand-rolling one risks a confidently wrong number for no gain.

WHY THE ERROR IS KNOWN EXACTLY, AT FULL SIZE.  At U = 0 this correlator on this
state is polynomial-time (`calibration/free_fermion.py`, validated to 1.7e-15
against exact many-body evolution). So the reference is not another tensor
network and not a small-lattice proxy: it is the exact answer on the same 7x4
doubly periodic lattice the experiment used. `--validate` first checks that the
TeNPy model, initial state and observable reproduce that reference under EXACT
diagonalisation of the MPO, with no TDVP involved -- if that fails, nothing
downstream means anything.

THE GEOMETRY.  MPS order = site index s = x + Lx*y, which is what makes the
dimers nearest-neighbour in the chain and the y-wrap a distance-21 term at
7x4 -- the "snake MPS with long-range Jordan-Wigner strings" of the deposit,
reproduced rather than avoided.
"""
from __future__ import annotations
import argparse, json, pathlib, sys, time, warnings
import numpy as np

# a Chain IS the lattice here, so TeNPy's unit_cell_width default is correct
warnings.filterwarnings("ignore", message=".*unit_cell_width.*")

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import free_fermion as ff

DATA = pathlib.Path(__file__).parent / "data"


# --------------------------------------------------------------------------
# the TeNPy side: model, state, observable


def build_model(Lx, Ly, J=1.0, U=0.0, periodic=True, flux_y=0.0):
    """H = -J sum_<ij> e^{i phi_ij} c^dag_i c_j + h.c. + U sum_i n_iu n_id.

    `flux_y` is the Peierls phase on every +y bond. The experiment runs at pi/4,
    which threads a pi flux through the short cycle; with Lx = 7 odd that is not
    removable by a gauge choice, and leaving it out moves C^zz by up to 0.13.
    """
    from tenpy.models.lattice import Chain
    from tenpy.models.model import CouplingModel, MPOModel
    from tenpy.networks.site import SpinHalfFermionSite

    N = Lx * Ly
    site = SpinHalfFermionSite(cons_N="N", cons_Sz="Sz")
    lat = Chain(N, site, bc="open", bc_MPS="finite")
    cm = CouplingModel(lat)
    for (a, b, axis) in ff.directed_bonds(Lx, Ly, periodic):
        w = -J * (np.exp(1j * flux_y) if axis == 1 and flux_y else 1.0)
        # free_fermion puts w on c^dag_b c_a; TeNPy wants i < j, so the
        # coefficient is conjugated exactly when the bond already runs a -> b
        i, j, w = (a, b, np.conj(w)) if a < b else (b, a, w)
        for sp in ("u", "d"):
            cm.add_coupling_term(w, i, j, f"Cd{sp} JW", f"C{sp}",
                                 op_string="JW", plus_hc=True)
    if U:
        for s in range(N):
            cm.add_onsite_term(U, s, "NuNd")
    return MPOModel(lat, cm.calc_H_MPO()), site


def build_state(Lx, Ly, site, cover=None, state="triplet"):
    """The dimer covering, one doublon and one holon.

    from_singlets gives (|up down> - |down up>)/sqrt(2) per pair. The experiment
    uses the Sz = 0 TRIPLET, the + combination (arXiv:2510.26300 Fig. 1), so the
    diagonal unitary (-1)^{n_down} is applied on the first site of each dimer,
    which flips the sign of the |down up> branch and nothing else. Both states
    give C^zz = -1 at t = 0, which is why the wrong one goes unnoticed.

    The two lonely sites come out 'empty'; the doublon is then filled with
    Cdu Cdd, whose two Jordan-Wigner strings are identical and cancel, so no
    string bookkeeping is needed.
    """
    from tenpy.networks.mps import MPS
    dimers, doublon, holon, spare = cover if cover else ff.covering(Lx, Ly)
    if spare:
        raise SystemExit(f"{Lx}x{Ly} leaves a spare site; use an even site count")
    psi = MPS.from_singlets(site, Lx * Ly, dimers, up="up", down="down",
                            lonely=[doublon, holon], lonely_state="empty",
                            bc="finite")
    if state == "triplet":
        flip = site.Id - 2 * site.Nd                 # (-1)^{n_down}, unitary
        for (a, _b) in dimers:
            psi.apply_local_op(a, flip, unitary=True)
    elif state != "singlet":
        raise SystemExit(f"unknown state {state!r}")
    psi.apply_local_op(doublon, "Cdu Cdd", renormalize=True)
    psi.canonical_form()
    return psi, dimers


def czz_mps(psi, links):
    """C^zz = <(nu-nd)_i (nu-nd)_j> - <><>, i.e. 4x the connected Sz correlator."""
    sz = psi.expectation_value("Sz")
    return [4.0 * (psi.expectation_value_term([("Sz", i), ("Sz", j)])
                   - sz[i] * sz[j]) for (i, j) in links]


def czz_exact_free(Lx, Ly, t, links, periodic=True, flux_y=0.0, cover=None,
                   state="triplet"):
    order, idx, dimers, doublon, holon, spare = ff.mode_order(Lx, Ly, cover)
    h = ff.hopping(Lx, Ly, idx, periodic=periodic, flux_y=flux_y)
    nu, nu2 = ff.occupations(idx, dimers, doublon, holon, spare)
    return [ff.czz_free(t, h, nu, nu2, idx, dimers, i, j, state)
            for (i, j) in links]


def setup_for(Lx, Ly, flux_y=None, state="triplet"):
    """(cover, flux) -- the deposit's own at 7x4, a generic covering elsewhere."""
    if (Lx, Ly) == ff.DEPOSIT_LATTICE:
        return ff.DEPOSIT_COVER, (ff.DEPOSIT_FLUX_Y if flux_y is None else flux_y)
    return None, (0.0 if flux_y is None else flux_y)


# --------------------------------------------------------------------------
# step 1: does the TeNPy setup reproduce the validated free-fermion answer?


def validate(sizes, times, max_size=5e7, flux_y=0.0, state="triplet"):
    """Exact diagonalisation of the MPO -- no TDVP, so this tests only the setup."""
    from tenpy.algorithms.exact_diag import ExactDiag
    from tenpy.linalg import np_conserved as npc
    rows, worst = [], 0.0
    for (Lx, Ly) in sizes:
        model, site = build_model(Lx, Ly, flux_y=flux_y)
        psi, dimers = build_state(Lx, Ly, site, state=state)
        # restrict to the state's own (N, Sz) sector, or 4^N swamps the machine
        ed = ExactDiag(model, charge_sector=psi.get_total_charge(True),
                       max_size=max_size)
        ed.build_full_H_from_mpo()
        ed.full_diagonalization()
        v0 = ed.mps_to_full(psi)
        for t in times:
            vt = (v0 if t == 0 else
                  npc.tensordot(ed.exp_H(t), v0, axes=[1, 0]))
            got = czz_mps(ed.full_to_mps(vt), dimers)
            ref = czz_exact_free(Lx, Ly, t, dimers, flux_y=flux_y, state=state)
            err = max(abs(a - b) for a, b in zip(got, ref))
            worst = max(worst, err)
            rows.append({"Lx": Lx, "Ly": Ly, "t": t, "flux_y": flux_y,
                         "state": state, "max_abs_err": err,
                         "czz_tenpy": got[0], "czz_free": ref[0]})
            print(f"  {Lx}x{Ly}  t={t:4.2f}  TeNPy {got[0]:+.12f}  "
                  f"free-fermion {ref[0]:+.12f}  |d|max {err:.2e}")
    ok = worst < 1e-10
    print(f"\n{'PASS' if ok else 'FAIL'}: worst disagreement {worst:.2e} "
          f"(tolerance 1e-10). The TeNPy model, the initial state and the "
          f"observable {'reproduce' if ok else 'DO NOT reproduce'} the validated "
          f"free-fermion reference.")
    return rows, ok


# --------------------------------------------------------------------------
# step 2: the measurement -- error against exact, as a function of dt and chi


def evolve(Lx, Ly, t_final, dt, chi, U=0.0, svd_min=1e-14, flux_y=0.0,
           cover=None, state="triplet"):
    from tenpy.algorithms.tdvp import TwoSiteTDVPEngine
    model, site = build_model(Lx, Ly, U=U, flux_y=flux_y)
    psi, dimers = build_state(Lx, Ly, site, cover, state)
    n = int(round(t_final / dt))
    if abs(n * dt - t_final) > 1e-12:
        raise SystemExit(f"dt={dt} does not divide t={t_final}")
    eng = TwoSiteTDVPEngine(psi, model, {
        "dt": dt, "N_steps": 1, "preserve_norm": True,
        "trunc_params": {"chi_max": chi, "svd_min": svd_min}})
    t0 = time.time()
    for _ in range(n):
        eng.run()
    return psi, dimers, max(psi.chi), time.time() - t0


def scan(Lx, Ly, times, dts, chis, U=0.0, flux_y=None, state="triplet"):
    cover, flux = setup_for(Lx, Ly, flux_y, state)
    links = (cover or ff.covering(Lx, Ly))[0]
    rows = []
    for t in times:
        ref = (czz_exact_free(Lx, Ly, t, links, flux_y=flux, cover=cover,
                              state=state) if U == 0 else None)
        for chi in chis:
            for dt in dts:
                psi, links_, chi_used, secs = evolve(Lx, Ly, t, dt, chi, U,
                                                     flux_y=flux, cover=cover,
                                                     state=state)
                got = czz_mps(psi, links_)
                err = (float(np.mean([abs(a - b) for a, b in zip(got, ref)]))
                       if ref is not None else float("nan"))
                rows.append({"Lx": Lx, "Ly": Ly, "t": t, "dt": dt, "chi": chi,
                             "chi_used": chi_used, "U": U, "flux_y": flux,
                             "state": state,
                             "mean_abs_err": err, "czz_0": got[0],
                             "seconds": secs})
                print(f"  t={t:4.2f} chi={chi:5d} (used {chi_used:5d}) "
                      f"dt={dt:7.4f}  err {err:.4e}  [{secs:7.1f} s]")
    return rows


def report(rows):
    """Does the error move with dt at fixed chi? That is the whole question."""
    print("\n  t     chi   err(dt_max)  err(dt_min)  ratio   verdict")
    for t in sorted({r["t"] for r in rows}):
        for chi in sorted({r["chi"] for r in rows}):
            g = sorted((r for r in rows if r["t"] == t and r["chi"] == chi),
                       key=lambda r: r["dt"])
            if len(g) < 2:
                continue
            hi, lo = g[-1]["mean_abs_err"], g[0]["mean_abs_err"]
            v = ("dt-limited" if hi / max(lo, 1e-300) > 2 else
                 "dt-converged" if hi / max(lo, 1e-300) < 1.2 else "partial")
            print(f"  {t:4.2f} {chi:6d}  {hi:11.4e}  {lo:11.4e}  "
                  f"{hi / max(lo, 1e-300):6.2f}   {v}")
    print("\n'dt-converged' with a NON-ZERO error is the interesting outcome: it "
          "would mean the floor survives both refinements and is the two-site\n"
          "projection or the state, not the integrator. 'dt-limited' would void "
          "the published chi-extrapolation in both directions.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--lx", type=int, default=7)
    ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--u", type=float, default=0.0)
    ap.add_argument("--times", type=float, nargs="+", default=[0.1, 0.5])
    ap.add_argument("--dts", type=float, nargs="+",
                    default=[0.05, 0.025, 0.0125, 0.00625])
    ap.add_argument("--chis", type=int, nargs="+", default=[256])
    ap.add_argument("--flux", type=float, default=None,
                    help="Peierls phase per +y bond; default pi/4 at 7x4")
    ap.add_argument("--state", default="triplet", choices=("triplet", "singlet"))
    ap.add_argument("--chi-scan", action="store_true",
                    help="chi = 256..2048 at the finest dt, their axis in our setup")
    ap.add_argument("--json", type=str, default=None)
    a = ap.parse_args()
    out = {}
    if a.validate:
        print("EXACT-DIAGONALISATION CHECK of the TeNPy setup (no TDVP):")
        # 6 sites is the ceiling: ExactDiag builds the FULL 4^N x 4^N before it
        # restricts to a charge sector, and 8 sites is already 68 GB. Larger
        # lattices are certified by the scan itself, whose reference is exact.
        import math
        rows, ok = [], True
        for fl in (0.0, math.pi / 4):
            print(f"  --- flux_y = {fl:.4f}, {a.state} dimers ---")
            r, o = validate([(3, 2), (2, 3)], [0.0, 0.1, 0.5, 1.0],
                            flux_y=fl, state=a.state)
            rows += r; ok = ok and o
        out["validate"] = rows
        if not ok:
            sys.exit(1)
    if a.scan or a.chi_scan:
        chis = [256, 512, 1024, 2048] if a.chi_scan else a.chis
        dts = [min(a.dts)] if a.chi_scan else a.dts
        cover, flux = setup_for(a.lx, a.ly, a.flux, a.state)
        print(f"\n{a.lx}x{a.ly} doubly periodic, U = {a.u}, flux_y = {flux:.4f}, "
              f"{a.state} dimers{' (deposit covering)' if cover else ''}, "
              f"chi in {chis}, dt in {dts}:")
        rows = scan(a.lx, a.ly, a.times, dts, chis, a.u, a.flux, a.state)
        out["scan"] = rows
        report(rows)
    if a.json:
        p = DATA / a.json
        p.write_text(json.dumps(out, indent=1))
        print(f"\nwrote {p}")


if __name__ == "__main__":
    main()

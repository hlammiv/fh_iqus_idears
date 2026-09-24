#!/usr/bin/env python3
"""How fast does the observable's support actually grow? (O5)

    python3 calibration/support_growth.py [--json]

O5 says the support damping model reproduces the ONE measured attenuation but
that its SCALING is a one-point fit: `support_growth = 0` cannot be right at long
times, "when the operator must eventually fill the lattice". That sentence is an
assumption. At U = 0 it is a calculation, and this file does it on the very
instance the attenuation was measured on -- 7x4 doubly periodic, pi flux, the
dimer covering read off the deposit (`free_fermion.py --deposit`, 2.9e-10).

WHAT IS COMPUTED.  In the Heisenberg picture at U = 0, n_p(t) = sum_ab R*_pa R_pb
c^dag_a c_b with R = exp(-i h t), so the single-particle weight of the evolved
operator on mode a is p_a(t) = mean_p |R_pa|^2 over the four modes of C^zz. That
is a probability distribution, and its inverse participation ratio

    w(t) = 1 / sum_a p_a(t)^2

is the EFFECTIVE number of modes -- equivalently qubits, under Jordan-Wigner --
the operator occupies. At t = 0 it is exactly 4, which is the model's `w_obs0`.

A gate sitting at circuit time s sees the observable evolved backwards from the
end by T - s, so the quantity that multiplies a whole circuit is the TIME AVERAGE
of w over [0, T], not its endpoint.

WHAT IT SETTLES.  The support does fill the lattice, and fast -- by t ~ 1 on this
instance. So the literal worry in O5 is correct. But feeding that measured
support into the damping model OVER-predicts the observed attenuation by about
four times, where the flat `w = 4` reproduces it. The conclusion is not that the
support stays small; it is that the damping does NOT track the support. Their own
weight test says the same thing from the other side: Lambda(ZZ)/Lambda(Z) = 1.91
tracks the BARE weights 2 and 1, which a spread-out operator could not do.

So `support_growth = 0` stops being an unconstrained fit and becomes the only one
of the three candidates that survives the measurement. The default stays "cone"
regardless -- this changes no headline, it removes a hole.
"""
from __future__ import annotations
import argparse, json, pathlib, sys
import numpy as np
from scipy.linalg import expm

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import free_fermion as ff

LAMBDA_OBSERVED = 0.200        # their 2415-gate circuit, arXiv:2510.26300
THEIR_GATES = 2415
T_WINDOW = 2.0                 # the evolution their circuit implements


def support_curve(h, idx, modes, times):
    """Effective number of modes (= qubits) the evolved observable occupies."""
    out = []
    for t in times:
        R = expm(-1j * h * t)
        p = np.mean([np.abs(R[q]) ** 2 for q in modes], axis=0)
        out.append(float(1.0 / np.sum(p ** 2)))
    return out


def deposit_support(n_t=81, t_max=T_WINDOW):
    o, idx, dm, D, H, sp, h = ff.deposit_setup()
    i, j = dm[0]
    modes = [idx[(i, 0)], idx[(i, 1)], idx[(j, 0)], idx[(j, 1)]]
    ts = np.linspace(0.0, t_max, n_t)
    w = support_curve(h, idx, modes, ts)
    return ts, np.array(w), h.shape[0]


def sizes_scan(sizes=((3, 2), (4, 2), (4, 4), (6, 4), (7, 4), (8, 6)),
               t_max=T_WINDOW, n_t=41):
    """When does the support saturate? The cone says t ~ L / (2 v)."""
    rows = []
    for (Lx, Ly) in sizes:
        o, idx, dm, D, H, spare = ff.mode_order(Lx, Ly)
        h = ff.hopping_sparse(Lx, Ly, idx, periodic=True).toarray()
        i, j = dm[0]
        modes = [idx[(i, 0)], idx[(i, 1)], idx[(j, 0)], idx[(j, 1)]]
        ts = np.linspace(0.0, t_max, n_t)
        w = np.array(support_curve(h, idx, modes, ts))
        M = 2 * Lx * Ly
        sat = next((t for t, x in zip(ts, w) if x > 0.5 * M), None)
        rows.append({"Lx": Lx, "Ly": Ly, "modes": M,
                     "t_half_fill": (float(sat) if sat is not None else None),
                     "w_mean": float(np.trapezoid(w, ts) / (ts[-1] - ts[0])),
                     "w_end": float(w[-1]),
                     "cone_at_tmax": 2.0 * 2.0 * float(ts[-1])})
        print(f"  {Lx}x{Ly}  {M:>3} modes   fills half by t = "
              f"{'%.2f' % sat if sat is not None else '  -- '}   "
              f"time-averaged w = {rows[-1]['w_mean']:6.2f}   w(T) = {w[-1]:6.2f}")
    return rows


def main(as_json=False):
    ts, w, M = deposit_support()
    w_mean = float(np.trapezoid(w, ts) / (ts[-1] - ts[0]))
    print(f"THE PUBLISHED INSTANCE: 7x4 doubly periodic, pi flux, {M} modes\n")
    print("     t   effective support w(t)   fraction of register   cone 2vt")
    for t in (0.0, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0):
        k = int(np.argmin(np.abs(ts - t)))
        print(f"  {t:4.2f}   {w[k]:20.2f}   {w[k] / M:20.3f}   {4.0 * t:8.1f}")
    print(f"\n  time-averaged over [0, {T_WINDOW}]: w = {w_mean:.2f} of {M} modes "
          f"({w_mean / M:.3f} of the register)")
    print(f"  the support DOES fill the lattice, by t ~ 1 -- O5's worry is right")

    # what each candidate damping fraction predicts for their measured circuit
    from fhcost import hubbard, nisq
    from fhcost.budget import DEFAULT
    SUP = DEFAULT.but(encoding="jw", damping_model="support",
                      tmax_mode="const", tmax_const=T_WINDOW, U_over_J=0)
    cc = hubbard.counts(28.0, SUP)
    per_gate = nisq.lambda_of(28.0, SUP) / cc["g_total"] / cc["damp_frac"]
    CONE = SUP.but(damping_model="cone")
    cands = [("cone, time-averaged", hubbard.counts(28.0, CONE)["damp_frac"]),
             ("free-fermion Heisenberg support (measured here)", w_mean / M),
             ("bare observable weight, w = 4", 4.0 / M)]
    print(f"\n  predicted attenuation on their {THEIR_GATES}-gate circuit:\n")
    print(f"  {'damping fraction from':>46}  {'frac':>6}  {'Lambda':>7}  vs 0.200")
    out = []
    for name, frac in cands:
        lam = per_gate * frac * THEIR_GATES
        out.append({"model": name, "frac": frac, "lambda": lam,
                    "ratio": lam / LAMBDA_OBSERVED})
        print(f"  {name:>46}  {frac:6.3f}  {lam:7.3f}  {lam / LAMBDA_OBSERVED:6.1f}x")
    print(f"\n  So the measured support over-predicts by "
          f"{out[1]['ratio']:.1f}x and the cone by {out[0]['ratio']:.0f}x, while "
          f"the BARE weight reproduces it.\n  Damping does not track the support. "
          f"Their own weight test agrees from the other side:\n  "
          f"Lambda(ZZ)/Lambda(Z) = 1.91 tracks the bare weights 2 and 1, which an "
          f"operator spread over\n  {w_mean:.0f} modes could not do.")

    print("\nWHEN DOES IT SATURATE?  (generic lattices, no flux, U = 0)\n")
    rows = sizes_scan()
    if as_json:
        p = HERE / "data" / "support_growth.json"
        p.write_text(json.dumps(
            {"deposit": {"t": ts.tolist(), "w": w.tolist(), "modes": M,
                         "w_mean": w_mean, "t_window": T_WINDOW},
             "candidates": out, "sizes": rows,
             "lambda_observed": LAMBDA_OBSERVED, "gates": THEIR_GATES}, indent=1))
        print(f"\nwrote {p}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    main(ap.parse_args().json)

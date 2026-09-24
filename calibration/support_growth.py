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


VELOCITY_THRESHOLDS = (1e-2, 1e-4, 1e-8, 1e-12)


def velocity_scan(L=29, times=(0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0),
                  thresholds=VELOCITY_THRESHOLDS):
    """How fast does the operator front move -- and is there ONE answer? (review #9)

    The model carries three velocities for three jobs: v = 2 sets the light cone
    and t_max, v_corr = 4 the correlation spread, v_lr = 20 the Lieb-Robinson
    bound. Review #9 called the taxonomy inconsistent. At U = 0 it can simply be
    measured, on an open lattice big enough that the front does not reach the
    edge, with the observable at the centre and the Chebyshev radius as the
    distance (the cone on a square lattice is a square).

    The answer is that the spreading velocity is THRESHOLD-DEPENDENT, which is
    why one number cannot serve: the exponentially small tail outruns the bulk.
    """
    o, idx, dm, D, H, sp = ff.mode_order(L, L)
    h = ff.hopping_sparse(L, L, idx, periodic=False).toarray()
    c = L // 2
    s0, s1 = c + L * c, c + 1 + L * c
    modes = [idx[(s0, 0)], idx[(s0, 1)], idx[(s1, 0)], idx[(s1, 1)]]
    site_of = {v: k[0] for k, v in idx.items()}
    dist = np.array([max(abs(site_of[a] % L - c), abs(site_of[a] // L - c))
                     for a in range(2 * L * L)], float)
    ds = sorted(set(dist))
    rows = []
    print(f"  {L}x{L} open lattice (half-width {c}), observable at the centre\n")
    print("     t   mean r" + "".join(f"  front {x:>7.0e}" for x in thresholds))
    for t in times:
        R = expm(-1j * h * t)
        p = np.mean([np.abs(R[q]) ** 2 for q in modes], axis=0)
        fr = [max([d for d in ds if p[dist == d].sum() > x] or [0.0])
              for x in thresholds]
        rows.append({"t": float(t), "mean_r": float(p @ dist),
                     "front": [float(x) for x in fr]})
        print(f"  {t:4.1f}  {rows[-1]['mean_r']:7.3f}"
              + "".join(f"  {x:12.1f}" for x in fr))
    T = np.array([r["t"] for r in rows])
    out = {"L": L, "half_width": c, "thresholds": list(thresholds), "rows": rows,
           "v_exact_max_group": 2.0}
    Y = np.array([r["mean_r"] for r in rows])
    out["v_mean"] = float(np.polyfit(T, Y, 1)[0])
    out["v_front"] = []
    print()
    print(f"  {'mean radius':16s} v = {out['v_mean']:.3f}")
    for i, x in enumerate(thresholds):
        Y = np.array([r["front"][i] for r in rows])
        ok = Y < c - 0.5
        v = float(np.polyfit(T[ok], Y[ok], 1)[0]) if ok.sum() >= 3 else None
        out["v_front"].append(v)
        print(f"  front at {x:<8.0e} v = {v:.3f}" if v else
              f"  front at {x:<8.0e} saturates too early")
    # ---- the CLUSTER BUFFER: xi ln(1/eps), the other half of review #9 ----
    # classical.cluster_radius charges v t + xi ln(1/eps) extra sites so that the
    # tail outside the cone is below eps. xi = 1 is a Config default with nothing
    # behind it. The same data measure it: buffer = front(eps) - v t.
    v0 = out["v_exact_max_group"]
    print(f"\n  cluster buffer, front(eps) - v t, with v = {v0:.0f}:\n")
    print("     t" + "".join(f"   eps {x:>7.0e}" for x in thresholds))
    buf = {x: [] for x in thresholds}
    for r in rows:
        line = f"  {r['t']:4.1f}"
        for i, x in enumerate(thresholds):
            b = r["front"][i] - v0 * r["t"]
            if r["front"][i] < c - 0.5:
                buf[x].append(b)
            line += f"  {b:11.1f}" + ("" if r["front"][i] < c - 0.5 else "*")
        print(line)
    print("   (* = front has reached the boundary, excluded)")
    out["buffer"] = {f"{x:.0e}": (float(np.mean(v)) if v else None)
                     for x, v in buf.items()}
    out["xi_eff"] = {f"{x:.0e}": (float(np.mean(v) / np.log(1.0 / x)) if v else None)
                     for x, v in buf.items()}
    print(f"\n  {'eps':>10} {'mean buffer':>12} {'ln(1/eps)':>10} {'implied xi':>11}")
    for x in thresholds:
        b, xi = out["buffer"][f"{x:.0e}"], out["xi_eff"][f"{x:.0e}"]
        if b is None:
            continue
        print(f"  {x:10.0e} {b:12.2f} {np.log(1/x):10.2f} {xi:11.3f}")
    xis = [v for v in out["xi_eff"].values() if v]
    print(f"\n  The buffer is flat in t at fixed eps, so the v t + xi ln(1/eps) "
          f"SHAPE is right.\n  The implied xi is {min(xis):.2f}-{max(xis):.2f}, "
          f"against the Config default xi = 1: the cluster\n  radius is charged "
          f"{1.0 / max(xis):.1f}-{1.0 / min(xis):.1f}x more buffer than the "
          f"measured tail needs. Conservative,\n  and it drifts DOWN as eps "
          f"falls because a free-fermion tail is super-exponential,\n  not "
          f"exponential -- so a single xi is a bound, not a fit.")

    print(f"\n  The exact max axial group velocity of -2J(cos kx + cos ky) is "
          f"2J = 2.000,\n  and the 1e-2 front measures "
          f"{out['v_front'][0]:.2f} -- that is the physical cone, and the model's "
          f"v = 2.\n  The tail is faster and keeps getting faster as the "
          f"threshold drops, which is\n  exactly why a Lieb-Robinson constant "
          f"has to exceed it. The three velocities the\n  model carries are "
          f"three different questions, and they are ordered correctly.")
    return out


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
    print("\nHOW FAST DOES THE FRONT MOVE?  (review #9, the velocity taxonomy)\n")
    vel = velocity_scan()
    if as_json:
        p = HERE / "data" / "support_growth.json"
        p.write_text(json.dumps(
            {"deposit": {"t": ts.tolist(), "w": w.tolist(), "modes": M,
                         "w_mean": w_mean, "t_window": T_WINDOW},
             "candidates": out, "sizes": rows, "velocity": vel,
             "lambda_observed": LAMBDA_OBSERVED, "gates": THEIR_GATES}, indent=1))
        print(f"\nwrote {p}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    main(ap.parse_args().json)

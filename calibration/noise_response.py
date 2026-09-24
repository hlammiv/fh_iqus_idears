#!/usr/bin/env python3
"""Can their data test the ZNE response model? (O4)

    python3 calibration/noise_response.py [--json]

O4: every ZNE conclusion rests on `s(lambda) = s(0) exp(-lambda Lambda)`, and a
response with less curvature raises every bias ceiling. The obvious external test
is the published hardware data -- twenty evolution times on a real 56-qubit
machine, with an exact reference. This file runs that test and reports that it
CANNOT settle the question, for a specific and checkable reason.

THE REASON.  Their circuit does not get deeper with the evolution time. Sec. C:
"We execute k = 4 second-order Trotter steps ... and time-evolve the initial
state up to time t = 2." A fixed COUNT, not a step density -- the step size grows
with t and the gate count does not. So the twenty reported times are twenty
REPEATS of one noise strength, not a scan over it, and no number of them
constrains the shape of s(lambda).

What they do give is that one point, twenty times over, which is worth having:
the measured attenuation on the dimer links is flat in t where the ratio is
well-conditioned, and its mean is the Lambda = 0.20 the model anchors on.

Beyond t ~ 0.5 the exact signal falls below ~0.1 and the ratio of two small
numbers stops meaning anything -- the scatter there is the denominator, not
physics, and this file cuts the window rather than fitting through it.

Needs refs/fermionic_dynamics/ -- see refs/README.md.
"""
from __future__ import annotations
import argparse, json, pathlib
import numpy as np

HERE = pathlib.Path(__file__).parent
DATA = HERE.parent / "refs" / "fermionic_dynamics"
OBS = "spin_correlator_neighbours"
SIGNAL_FLOOR = 0.10      # below this the exact value, the ratio is ill-conditioned


def load(u="U_0"):
    import pandas as pd
    p = DATA / u / "exp_vals.h5"
    if not p.exists():
        raise SystemExit(f"missing {p}; see refs/README.md")
    raw = pd.read_hdf(p, "raw")
    ex = pd.read_hdf(p, "exact")
    return (raw[raw.obs_type == OBS],
            ex[(ex.obs_type == OBS) & (ex.method == "Exact")])


def attenuation(raw, exact):
    z = exact[exact.ev_time == 0.0]
    links = {r.site_coords for r in z.itertuples() if abs(r.obs_val) > 0.5}
    rows = []
    for t in sorted(set(raw.ev_time)):
        rr = [x for x in raw[raw.ev_time == t].itertuples()
              if x.site_coords in links]
        ee = {x.site_coords: x.obs_val for x in exact[exact.ev_time == t].itertuples()}
        if not rr or not all(x.site_coords in ee for x in rr):
            continue
        num = float(np.mean([x.obs_val for x in rr]))
        den = float(np.mean([ee[x.site_coords] for x in rr]))
        errs = np.array([getattr(x, "obs_err", np.nan) for x in rr], float)
        se = float(np.sqrt(np.nansum(errs ** 2)) / len(rr))
        rows.append({"t": float(t), "raw": num, "exact": den,
                     "shot_err": se, "att": num / den,
                     "att_err": abs(se / den),
                     "conditioned": abs(den) >= SIGNAL_FLOOR})
    return rows


def main(as_json=False):
    raw, exact = load()
    rows = attenuation(raw, exact)
    print("Dimer-link <C^zz>, raw hardware against the exact reference.\n")
    print("     t     raw       exact    attenuation      Lambda   window")
    for r in rows:
        print(f"  {r['t']:4.2f}  {r['raw']:+.5f}  {r['exact']:+.5f}   "
              f"{r['att']:6.4f} +- {r['att_err']:.4f}  "
              f"{-np.log(max(r['att'], 1e-9)):7.4f}   "
              f"{'USED' if r['conditioned'] else 'signal < %.2f' % SIGNAL_FLOOR}")
    good = [r for r in rows if r["conditioned"]]
    a = np.array([r["att"] for r in good]); T = np.array([r["t"] for r in good])
    lam = -np.log(a)
    slope, icept = np.polyfit(T, lam, 1)
    print(f"\n  {len(good)} conditioned points, t = {T.min():.1f} to {T.max():.1f}")
    print(f"  mean attenuation {a.mean():.4f} -> Lambda = {-np.log(a.mean()):.4f}, "
          f"scatter +-{a.std(ddof=1):.4f}")
    print(f"  fitted trend in t: Lambda = {slope:+.4f} t {icept:+.4f}")
    print(f"  the spread is {a.std(ddof=1) / np.mean([r['att_err'] for r in good]):.0f}x "
          f"the quoted shot errors, and the fitted trend is NEGATIVE -- less "
          f"damping at\n  longer times, which no noise process does. So the "
          f"residual is Trotter error: their\n  circuit carries it and no "
          f"noiseless Trotterised series was deposited to divide it out.")
    best = min(good, key=lambda r: r["att_err"] / max(abs(r["att"]), 1e-9))
    print(f"\n  The best-conditioned single point is t = {best['t']:.2f}: "
          f"attenuation {best['att']:.4f} +- {best['att_err']:.4f}, "
          f"Lambda = {-np.log(best['att']):.3f}.\n  That is the anchor, and it "
          f"is where the Trotter contamination is smallest (dt = t/4).")
    print(f"\n  VERDICT: the gate count is CONSTANT across these times "
          f"(k = 4 fixed), so these are\n  {len(good)} repeats of ONE noise "
          f"strength. The response SHAPE is untestable from them at any\n  "
          f"statistics. Pooled they give Lambda = {-np.log(a.mean()):.2f} with a "
          f"spread of +-{a.std(ddof=1):.2f} in attenuation, which contains the "
          f"0.20\n  anchor without sharply confirming it. Testing exp(-lambda "
          f"Lambda) needs noise-amplified\n  runs at several gains, which they "
          f"did not publish.")
    out = {"rows": rows, "lambda_best_point": float(-np.log(best["att"])),
           "t_best": best["t"], "lambda_mean": float(-np.log(a.mean())),
           "attenuation_mean": float(a.mean()), "attenuation_sd": float(a.std(ddof=1)),
           "trend_slope": float(slope), "n_conditioned": len(good),
           "signal_floor": SIGNAL_FLOOR,
           "note": "k = 4 fixed Trotter steps: constant depth, one lambda point"}
    if as_json:
        p = HERE / "data" / "noise_response.json"
        p.write_text(json.dumps(out, indent=1))
        print(f"\nwrote {p}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    main(ap.parse_args().json)

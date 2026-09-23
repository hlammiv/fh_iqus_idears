#!/usr/bin/env python3
"""The C^zz signal: the envelope fit, its residuals, and a held-out test.

Second-pass review #8. The costing model set its absolute tolerance from a
smooth monotone Gaussian envelope. The published data is neither smooth nor
monotone: at U/J = 4 it falls to 0.0466 at t = 1.5 and comes back up to 0.0532
at t = 2.0. The envelope reads 0.0640 there -- 37% high, which is 1.89x too few
shots at the point that actually binds.

This script keeps what the review asked to be preserved: the data, the fit, the
residuals, and a held-out prediction test. It is the single source of the data
table; `fhcost.hubbard` imports SIGNAL_DATA from here in spirit (the values are
mirrored there with this file named as their provenance).

    python3 calibration/fit_signal.py            # fit + residuals + held-out
    python3 calibration/fit_signal.py --json     # also write data/signal_fit.json

WHAT THE FIT IS FOR, AND WHAT IT IS NOT FOR
  The envelope is a physical description: |C^zz| starts at 1 (exact triplet
  algebra), melts on a timescale t_melt, and settles to an antiferromagnetic
  residual. It is good to ~1% out to t = 0.7.

  It is NOT a lower bound on the signal, and the costing model needs one. In the
  tail the envelope wanders from 0.69x to 1.46x the data -- it both over- and
  under-predicts -- so using it to set a tolerance silently loosens the target
  near a minimum. The costing path therefore uses the DATA, reduced by a stated
  uncertainty, and the envelope only for extrapolation beyond the measured
  range.
"""
from __future__ import annotations
import argparse, json, math, pathlib
import numpy as np

# Published dimer-link |C^zz| (TFLO + GPR), Zenodo 17799843, accompanying
# arXiv:2510.26300. 7x4 double-periodic lattice, flux pi, 160 shots per point.
# These are MITIGATED experimental values, not exact ground truth.
SIGNAL_DATA = {
    0.0: [(0.1, 0.9500), (0.3, 0.6378), (0.5, 0.2810), (0.7, 0.0854),
          (1.0, 0.0680), (1.5, 0.0326), (2.0, 0.0295)],
    4.0: [(0.1, 0.9527), (0.3, 0.6603), (0.5, 0.3292), (0.7, 0.1414),
          (1.0, 0.0902), (1.5, 0.0466), (2.0, 0.0532)],
}
DATA_RANGE = (0.1, 2.0)


def envelope(t, s_res, t_melt, beta, s_short=1.0):
    return s_res + (s_short - s_res) * np.exp(-((np.asarray(t) / t_melt) ** beta))


def fit_one(rows, free_beta=True):
    """Least squares on (s_res, t_melt, beta). s_short is FIXED at the exact 1."""
    from scipy.optimize import curve_fit
    t = np.array([r[0] for r in rows])
    y = np.array([r[1] for r in rows])
    if free_beta:
        f = lambda tt, sr, tm, b: envelope(tt, sr, tm, b)
        p0, bounds = [0.04, 0.43, 2.0], ([0, 0.1, 0.5], [0.5, 2.0, 4.0])
    else:
        f = lambda tt, sr, tm: envelope(tt, sr, tm, 2.0)
        p0, bounds = [0.04, 0.43], ([0, 0.1], [0.5, 2.0])
    p, cov = curve_fit(f, t, y, p0=p0, bounds=bounds, maxfev=20000)
    pred = f(t, *p)
    return p, np.sqrt(np.diag(cov)), y - pred


def held_out(rows):
    """Leave-one-out: refit without each point and predict it.

    Review #8, test 1. Seven points and three parameters, so this is a weak test
    by construction -- it is reported with that caveat rather than as validation.
    """
    out = []
    for k in range(len(rows)):
        sub = rows[:k] + rows[k + 1:]
        try:
            p, _, _ = fit_one(sub)
        except Exception:
            continue
        t, y = rows[k]
        pred = float(envelope(t, *p))
        out.append({"t": t, "data": y, "predicted": pred,
                    "rel_err": (pred - y) / max(y, 1e-12)})
    return out


def report():
    res = {"data": {str(u): v for u, v in SIGNAL_DATA.items()}, "fits": {}}
    for u, rows in SIGNAL_DATA.items():
        p, sd, r = fit_one(rows)
        p2, sd2, r2 = fit_one(rows, free_beta=False)
        ho = held_out(rows)
        res["fits"][str(u)] = {
            "free_beta": {"s_res": p[0], "t_melt": p[1], "beta": p[2],
                          "sigma": list(map(float, sd)),
                          "residuals": list(map(float, r))},
            "beta_fixed_2": {"s_res": p2[0], "t_melt": p2[1],
                             "sigma": list(map(float, sd2)),
                             "residuals": list(map(float, r2))},
            "held_out": ho,
            "data_min": min(v for _, v in rows),
            "data_min_t": min(rows, key=lambda x: x[1])[0],
            "monotone": all(rows[i][1] >= rows[i + 1][1]
                            for i in range(len(rows) - 1)),
        }
        print(f"U/J = {u}")
        print(f"  free beta : s_res {p[0]:.4f}+-{sd[0]:.4f}  t_melt {p[1]:.4f}"
              f"+-{sd[1]:.4f}  beta {p[2]:.3f}+-{sd[2]:.3f}")
        print(f"  beta = 2  : s_res {p2[0]:.4f}+-{sd2[0]:.4f}  t_melt {p2[1]:.4f}"
              f"+-{sd2[1]:.4f}   RMS {np.sqrt((r2**2).mean()):.4f}")
        print(f"  data is monotone: {res['fits'][str(u)]['monotone']}   "
              f"minimum {res['fits'][str(u)]['data_min']:.4f} at "
              f"t = {res['fits'][str(u)]['data_min_t']}")
        print("  residuals (beta=2), model - data:")
        for (t, y), rr in zip(rows, r2):
            print(f"     t={t:4.1f}  data {y:.4f}  resid {rr:+.4f}  "
                  f"ratio {(y + rr) / y:5.2f}")
        print("  leave-one-out:")
        for h in ho:
            print(f"     t={h['t']:4.1f}  data {h['data']:.4f}  "
                  f"predicted {h['predicted']:.4f}  rel {h['rel_err']:+.1%}")
        print()

    # the number that matters for costing: how far the envelope can sit ABOVE
    # the data, since that is what loosens a tolerance
    worst = 0.0
    for u, rows in SIGNAL_DATA.items():
        p2, _, _ = fit_one(rows, free_beta=False)
        for t, y in rows:
            worst = max(worst, float(envelope(t, p2[0], p2[1], 2.0)) / y)
    res["worst_envelope_over_data"] = worst
    res["implied_shot_factor"] = worst ** 2
    print(f"worst envelope/data = {worst:.3f}  -> the envelope can ask for "
          f"{worst**2:.2f}x too few shots")
    print(f"a lower bound therefore needs a factor <= {1/worst:.3f}")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    r = report()
    if a.json:
        p = pathlib.Path(__file__).parent / "data" / "signal_fit.json"
        p.parent.mkdir(exist_ok=True)
        p.write_text(json.dumps(r, indent=1))
        print(f"wrote {p}")

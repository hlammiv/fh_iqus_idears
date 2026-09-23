#!/usr/bin/env python3
"""How far does the Trotter calibration actually reach? (second-pass #2, tests 2/4/5)

    python3 calibration/domain_check.py [--json]

Three of the review's five tests for finding #2 can be answered from the raw
errors already in data/trotter_cal.json, without any new simulation. Two of them
come back badly, and the model's claims have been narrowed to match.

TEST 4 -- is err ~ r^-2 at the OPERATING step count, or only asymptotically?
  Clean at r >= 8 (slopes -2.00 to -2.11). It breaks at small r and long time:
  at tau = 2 the slope over r = 1-4 is between -0.58 and +0.63, because the
  error there is O(1) and saturated rather than asymptotic. The second-order arm
  operates at r = 11-400, comfortably inside the good window. The MULTIPRODUCT
  arm operates at r = 3.2 upward and is below r = 8 for m <= 13, so it sits in
  the window where the assumed power law is not reliable.

TEST 2 -- hold a time out of the fit and predict it.
  Drop tau = 1 and interpolate it from tau = 0.5 and 2: the prediction overshoots
  by +124% at n = 9 and +86% at n = 12. Interpolation INSIDE the measured range
  is good only to about a factor of two at long time, which bounds what the
  clamped extrapolation beyond tau = 2 can be worth.

TEST 5 -- the size-independence that licenses extrapolating in m.
  It is a short-time statement and it decays with time. Spread of W_eff across
  n >= 6:

      tau = 0.25   1.01x     tau = 1.0   1.26x
      tau = 0.5    1.27x     tau = 2.0   2.91x

  The adopted convention t = sqrt(m)/v puts EVERY plotted point at tau >= 1 and
  most far beyond tau = 2, which is precisely where the justification weakens.

WHAT THIS DOES NOT SETTLE
  The trajectory itself. Every number here comes from fixed-time size sweeps.
  Measuring W_eff at the pairs (m, tau = sqrt(m)/v) is review test 3 and needs
  new runs; see the STILL NEEDED block printed at the end.
"""
from __future__ import annotations
import argparse, json, math, pathlib
import numpy as np

HERE = pathlib.Path(__file__).parent
RAW = HERE / "data" / "trotter_cal.json"
V_B = 2.0
R_ASYMPTOTIC = 8          # the window the stored coefficients were fitted on


def load():
    return [r for r in json.loads(RAW.read_text())["rows"] if r["steps"] > 0]


def w_eff(rows, n, tau, rmin=R_ASYMPTOTIC, rmax=64):
    v = [r["abs_err"] * r["steps"] ** 2 / tau ** 3 for r in rows
         if r["n"] == n and r["tau"] == tau and rmin <= r["steps"] <= rmax
         and r["abs_err"] > 0]
    return float(np.median(v)) if v else None


def slope(rows, n, tau, rmin, rmax):
    v = [(r["steps"], r["abs_err"]) for r in rows
         if r["n"] == n and r["tau"] == tau and rmin <= r["steps"] <= rmax
         and r["abs_err"] > 0]
    if len(v) < 2:
        return None
    s, e = zip(*sorted(v))
    return float(np.polyfit(np.log(s), np.log(e), 1)[0])


def main(as_json=False):
    rows = load()
    ns = sorted({r["n"] for r in rows})
    taus = sorted({r["tau"] for r in rows})
    out = {"n": ns, "tau": taus, "r_asymptotic": R_ASYMPTOTIC}

    print("TEST 4  err ~ r^-2 at the operating step count?\n")
    print(f"  {'n':>3} {'tau':>5} {'slope r=1-4':>12} {'slope r>=8':>11} "
          f"{'W(r>=8)':>9} {'W(r=1-4)':>9}")
    t4 = []
    for n in ns:
        for tau in taus:
            a, b = slope(rows, n, tau, 1, 4), slope(rows, n, tau, 8, 64)
            wa, wb = w_eff(rows, n, tau, 1, 4), w_eff(rows, n, tau)
            if None in (a, b, wa, wb):
                continue
            t4.append({"n": n, "tau": tau, "slope_small_r": a, "slope_large_r": b,
                       "W_large_r": wb, "W_small_r": wa})
            print(f"  {n:>3} {tau:>5} {a:>12.2f} {b:>11.2f} {wb:>9.4f} {wa:>9.4f}")
    out["power_law"] = t4
    worst = max(abs(x["slope_small_r"] + 2.0) for x in t4)
    best = max(abs(x["slope_large_r"] + 2.0) for x in t4)
    print(f"\n  worst |slope + 2|: {worst:.2f} at r = 1-4, {best:.2f} at r >= 8")
    out["worst_dev_small_r"], out["worst_dev_large_r"] = worst, best

    print("\nTEST 5  size-independence, tau by tau (n >= 6)\n")
    print(f"  {'tau':>5} " + " ".join(f"{'n=' + str(n):>9}" for n in ns) + "   spread")
    t5 = {}
    for tau in taus:
        ws = [w_eff(rows, n, tau) for n in ns]
        big = [w for n, w in zip(ns, ws) if w and n >= 6]
        sp = max(big) / min(big)
        t5[str(tau)] = sp
        print(f"  {tau:>5} " + " ".join(f"{(f'{w:.4f}' if w else '-'):>9}" for w in ws)
              + f"   {sp:5.2f}x")
    out["size_spread"] = t5

    print("\nTEST 2  hold tau = 1 out and predict it from 0.5 and 2\n")
    t2 = []
    f = (math.log(1.0) - math.log(0.5)) / (math.log(2.0) - math.log(0.5))
    for n in ns:
        a, b, c = w_eff(rows, n, 0.5), w_eff(rows, n, 1.0), w_eff(rows, n, 2.0)
        if None in (a, b, c):
            continue
        pred = math.exp((1 - f) * math.log(a) + f * math.log(c))
        t2.append({"n": n, "measured": b, "predicted": pred, "rel_err": (pred - b) / b})
        print(f"  n={n:>2}: measured {b:.4f}  predicted {pred:.4f}  "
              f"{(pred - b) / b:+.0%}")
    out["time_holdout"] = t2
    out["worst_time_holdout"] = max(abs(x["rel_err"]) for x in t2)
    print(f"\n  worst time hold-out error: {out['worst_time_holdout']:+.0%}")

    print("\n" + "=" * 70)
    print("STILL NEEDED to close finding #2 (review test 3)\n")
    print("  Measure W_eff along the ADOPTED trajectory tau = sqrt(m)/v, not on a")
    print("  fixed-time size sweep. Everything above is a fixed-time sweep, and")
    print("  the whole m^(7/4) claim rests on the trajectory instead:\n")
    print(f"  {'patch':>14} {'sites':>6} {'tau':>6} {'sector dim':>13} {'GB/vec':>8}  where")
    for name, n in (("square2", 4), ("cross5", 5), ("rectangle2x3", 6),
                    ("square3", 9), ("rectangle3x4", 12), ("rectangle2x7", 14),
                    ("square4", 16)):
        up, dn = (n + 1) // 2, n // 2
        dim = math.comb(n, up) * math.comb(n, dn)
        gb = dim * 16 / 2 ** 30
        print(f"  {name:>14} {n:>6} {math.sqrt(n) / V_B:>6.2f} {dim:>13,} "
              f"{gb:>8.3f}  {'local' if gb < 0.5 else 'lenore'}")
    print("\n  n <= 14 fits locally (< 1 GB peak). n = 16 needs ~10 GB and belongs")
    print("  on lenore_remote. Also queued: the O11b multiproduct rerun (finer r,")
    print("  state infidelity alongside the observable error, more tau).")

    if as_json:
        p = HERE / "data" / "domain_check.json"
        p.write_text(json.dumps(out, indent=1))
        print(f"\nwrote {p}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    main(ap.parse_args().json)

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

    traj = HERE / "data" / "trotter_traj.json"
    if traj.exists():
        import sys as _sys
        _sys.path.insert(0, str(HERE.parent))
        from fhcost import hubbard as _H
        T = json.loads(traj.read_text())
        print("\nTEST 3  W_eff AT the operating points, tau = sqrt(m)/v_B\n")
        print("  The model does not assume W_eff is CONSTANT along the trajectory")
        print("  -- it assumes W depends on tau alone, and reads the fixed-time")
        print("  table at tau = sqrt(m)/v. That is the assumption under test.\n")
        print(f"  {'n':>3} {'tau':>6} {'measured':>9} {'table':>9} {'ratio':>7}"
              f"  substep conv")
        pts = []
        CONV_TOL = 1e-10     # the reference must be far better than the signal
        for mm in sorted(T["meta"], key=lambda x: x["n"]):
            rs = [r for r in T["rows"] if r["patch"] == mm["patch"]
                  and r["steps"] >= R_ASYMPTOTIC and r["W_eff"]]
            if not rs:
                continue
            if mm.get("substep_convergence", 0) > CONV_TOL:
                # the memory-capped Krylov basis was not compensated by enough
                # substepping, so the "exact" reference is wrong at a level that
                # swamps the Trotter error being measured. Report and skip.
                print(f"  {mm['n']:>3} {mm['tau']:>6.3f} {'EXCLUDED':>9} "
                      f"{'':>9} {'':>7}  reference conv "
                      f"{mm['substep_convergence']:.1e} > {CONV_TOL:.0e}")
                continue
            w = float(np.median([r["W_eff"] for r in rs]))
            tab = _H.w_measured(mm["tau"], 4.0)
            pts.append({"n": mm["n"], "tau": mm["tau"], "W_measured": w,
                        "W_table": tab, "ratio": w / tab,
                        "substep_convergence": mm["substep_convergence"]})
            print(f"  {mm['n']:>3} {mm['tau']:>6.3f} {w:>9.5f} {tab:>9.5f} "
                  f"{w / tab:>7.2f}  {mm['substep_convergence']:.1e}")
        out["trajectory"] = pts
        big = [x for x in pts if x["n"] >= 6]
        if len(big) >= 3:
            rr = [x["ratio"] for x in big]
            out["traj_ratio_range"] = [min(rr), max(rr)]
            print(f"\n  the fixed-time table predicts the trajectory to within "
                  f"{min(rr):.2f}-{max(rr):.2f}x (n >= 6)")
            # implied exponent: r = t^1.5 sqrt(W/eps), t = sqrt(m)/v, so
            # G = c_g m r ~ m^1.75 * sqrt(W(m)); fit W against m directly
            ns = np.array([x["n"] for x in big], float)
            ws = np.array([x["W_measured"] for x in big])
            gam = float(np.polyfit(np.log(ns), np.log(ws), 1)[0])
            alpha = 1.75 + gam / 2.0
            tabs = np.array([x["W_table"] for x in big])
            gam_t = float(np.polyfit(np.log(ns), np.log(tabs), 1)[0])
            out["alpha_measured"], out["alpha_table"] = alpha, 1.75 + gam_t / 2.0
            print(f"  W ~ m^{gam:+.3f} measured, m^{gam_t:+.3f} from the table")
            print(f"  implied alpha in G ~ m^alpha: {alpha:.3f} measured, "
                  f"{1.75 + gam_t / 2.0:.3f} from the table "
                  f"(Campbell's bound gives 2.25)")
            print(f"  at m = 400 the measured and table exponents differ by "
                  f"{400.0 ** (alpha - (1.75 + gam_t / 2.0)):.2f}x in gate count")

    clamp = HERE / "data" / "clamp_probe.json"
    if clamp.exists():
        import sys as _s2
        _s2.path.insert(0, str(HERE.parent))
        from fhcost import hubbard as _H2
        C = json.loads(clamp.read_text())
        print("\nTHE CLAMP  W_eff BEYOND tau = 2, where alpha = 1.75 comes from\n")
        print("  The model clamps W at its tau = 2 value and every plotted point")
        print("  sits in the clamped region, so the clamp -- not the measured")
        print("  tau-dependence -- is what sets the exponent.\n")
        print(f"  {'tau':>5} " + "".join(f"{'n=' + str(n):>11}"
                                         for n in sorted({m['n'] for m in C['meta']}))
              + f"{'clamp':>10} {'agree':>7}")
        ns = sorted({m["n"] for m in C["meta"]})
        rows_out, agree_all = [], True
        for tau in C["clamp_taus"]:
            got = {m["n"]: m["W_eff"] for m in C["meta"] if m["tau"] == tau}
            vals = [got[n] for n in ns if n in got]
            clamped = _H2.w_measured(tau, 4.0)     # the model's clamped value
            ag = (max(vals) / min(vals)) if len(vals) > 1 else float("nan")
            agree_all &= (ag < 1.5) if len(vals) > 1 else True
            rows_out.append({"tau": tau, "W": got, "clamped": clamped, "spread": ag})
            print(f"  {tau:>5} " + "".join(f"{got.get(n, float('nan')):>11.5f}"
                                           for n in ns)
                  + f"{clamped:>10.5f} {ag:>7.2f}x")
        out["clamp"] = rows_out
        big = [r for r in rows_out if len(r["W"]) > 1]
        if big:
            if agree_all:
                ratio = np.mean([min(r["W"].values()) / r["clamped"] for r in big])
                print(f"\n  The two sizes AGREE (spread < 1.5x at every tau), so the")
                print(f"  tau-trend past 2 is real and not a finite-size artefact.")
                print(f"  Measured W is {ratio:.2f}x the clamped value, so the clamp is")
                print(f"  {'CONSERVATIVE' if ratio < 1 else 'OPTIMISTIC'} and alpha = 1.75 is "
                      f"{'an upper bound' if ratio < 1 else 'too low'}.")
                # r ~ sqrt(W) so the step count moves as sqrt of this
                print(f"  Step count moves by {math.sqrt(ratio):.2f}x at fixed tau.")
            else:
                print(f"\n  The two sizes DISAGREE, so a lattice this small cannot be")
                print(f"  asked about tau > 2: at tau = 4 the cone spans 2 v tau = 16")
                print(f"  sites, wider than either patch. The clamp stands for want of")
                print(f"  evidence, not because it was tested.")

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

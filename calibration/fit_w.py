#!/usr/bin/env python3
"""Regenerate fhcost.hubbard.W_MEASURED from the raw calibration data.

The coefficients in hubbard.py were previously hard-coded with no way to check
them. This script derives them from calibration/data/*.json and prints the dict
verbatim, so the model's central empirical input is reproducible from the repo.

    python3 calibration/fit_w.py            # print the dict and the domain
    python3 calibration/fit_w.py --holdout  # leave-one-out prediction test

SELECTION, stated explicitly because it is a choice:
  * W_eff = err * r^2 / tau^3, the second-order form. Justified by W_eff being
    constant in r to four digits over r = 4..64 (printed below).
  * medians over patches with n >= 9 and step counts r >= 8. n = 4,5,6 are
    excluded because the lattice clips the observable's causal cone there; the
    exclusion is visible in the printed table.
"""
import json, math, pathlib, statistics as st, sys

HERE = pathlib.Path(__file__).resolve().parent
FILES = {0.0: "cal_u0.json", 4.0: "trotter_cal.json", 8.0: "cal_u8.json"}
TAUS = (0.25, 0.5, 1.0, 2.0)
N_MIN, R_MIN = 9, 8


def load():
    out = {}
    for u, fn in FILES.items():
        d = json.load(open(HERE / "data" / fn))
        assert abs(d["u_over_j"] - u) < 1e-9, (fn, d["u_over_j"])
        out[u] = [r for r in d["rows"] if r["steps"] > 0]
    return out


def w_table(rows, n_min=N_MIN, r_min=R_MIN, drop=()):
    w = {}
    for tau in TAUS:
        if (tau,) in drop:
            continue
        vals = [r["abs_err"] * r["steps"] ** 2 / tau ** 3 for r in rows
                if r["tau"] == tau and r["n"] >= n_min and r["steps"] >= r_min]
        if vals:
            w[tau] = st.median(vals)
    return w


def main():
    data = load()
    sizes = sorted({r["n"] for rs in data.values() for r in rs})
    print(f"domain: n in {sizes}, tau in {list(TAUS)}, U/J in {sorted(FILES)}\n")

    print("W_eff constancy in r (U=4, n=9) -- justifies the tau^3/r^2 form")
    for tau in TAUS:
        vs = [(r["steps"], r["abs_err"] * r["steps"] ** 2 / tau ** 3)
              for r in data[4.0] if r["tau"] == tau and r["n"] == 9]
        vs.sort()
        print(f"  tau={tau}: " + "  ".join(f"r={k}:{v:.4f}" for k, v in vs if k >= 4))

    print("\nW_MEASURED = {")
    for u in sorted(FILES):
        w = w_table(data[u])
        body = ", ".join(f"{t}: {w[t]:.4f}" for t in sorted(w))
        print(f"    {u}: {{{body}}},")
    print("}")

    if "--holdout" in sys.argv:
        print("\nleave-one-out: refit without a size, predict it")
        print(f"  {'U':>4}{'tau':>6}{'full':>10}{'without n=12':>14}{'shift':>8}")
        for u in sorted(FILES):
            full = w_table(data[u])
            held = w_table([r for r in data[u] if r["n"] != 12])
            for tau in TAUS:
                if tau in full and tau in held:
                    print(f"  {u:>4.0f}{tau:>6}{full[tau]:>10.4f}{held[tau]:>14.4f}"
                          f"{100*(held[tau]/full[tau]-1):>7.1f}%")


if __name__ == "__main__":
    main()


def mpf_table(src="data/trotter_cal.json", taus=(0.25, 0.5), r_min=4,
              report_spread=False):
    """Higher-order multiproduct coefficients W_2k from err = W_2k t^(2k+1)/r^(2k).

    Only r >= r_min is used: r = 1, 2 are not yet asymptotic. The default taus
    stop at 0.5 because at tau >= 1 the observable error passes through zero
    crossings and the extracted coefficient swings by 3-13x -- which is O11b,
    and why `src` exists: the finely-sampled rerun steps over the crossings.

    report_spread returns max/min across base step counts as well as the median,
    so a swinging coefficient is visible instead of being hidden by the median.
    """
    rows = json.loads((HERE / src).read_text())["rows"]
    out = {}
    for order in (4, 6, 8):
        out[order] = {}
        for tau in taus:
            vals = [r["abs_err"] * (-r["steps"]) ** order / tau ** (order + 1)
                    for r in rows
                    if r["steps"] < 0 and r.get("mp_order") == order
                    and r["tau"] == tau and r["n"] >= N_MIN
                    and -r["steps"] >= r_min and r["abs_err"] > 0]
            if vals:
                out[order][tau] = ((st.median(vals), max(vals) / min(vals), len(vals))
                                   if report_spread else st.median(vals))
    return out


if "--mpf" in sys.argv:
    t = mpf_table()
    print("\nW_MPF = {   # order: {tau: W_2k}, r >= 4, n >= 9, U/J = 4, tau <= 0.5 only")
    for o in sorted(t):
        body = ", ".join(f"{k}: {v:.5f}" for k, v in sorted(t[o].items()))
        print(f"    {o}: {{{body}}},")
    print("}")

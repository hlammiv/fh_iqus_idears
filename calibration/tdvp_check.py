#!/usr/bin/env python3
"""Is the published TDVP error truncation-limited? (O10)

    python3 calibration/tdvp_check.py [--json]

NO NEW TDVP CODE WAS WRITTEN FOR THIS. The experimental collaboration released
their full simulation output (Zenodo 17799843, `fermionic_dynamics.zip`, 112 MB,
CC-BY-4.0), which contains TDVP at chi = 256, 512, 1024, 2048 *and* the exact
fermionic-linear-optics reference at U = 0, for every observable including the
one this model costs. Their metadata says `Tenpy $\\chi=512$`: the runs are
TeNPy, which is also what any new run here would use.

Put the extracted `U_0/` and `U_4/` directories under refs/fermionic_dynamics/
(gitignored -- 1 GB). Needs `pip install tables uncertainties` to read their
pandas/pytables HDF5.

WHAT THE MODEL CLAIMED, AND WHAT THE FULL DATA SAYS
  classical.TDVP_FLOOR asserted a "chi-independent error floor" at t = 0.1, and
  used it to argue the published TDVP is not bond-dimension-limited, so its
  error cannot bound tensor-network capability. Those stored numbers reproduce
  EXACTLY from the full release -- 7.93e-3, 7.38e-3, 7.36e-3, 7.35e-3 on the
  dimer links -- so the extraction was right.

  The claim is right at t = 0.1 and wrong as a general statement. The chi
  dependence is NON-MONOTONIC in time:

      t = 0.1   slope -0.03   floored; chi = 256 is already ample
      t = 0.5   slope -0.45   genuinely truncation-limited, 2.5x over 8x in chi
      t = 1.0   slope -0.26   converging
      t = 1.5   slope -0.07   stalling
      t = 2.0   slope -0.04   stalled, and the error EQUALS the signal

  So "not truncation-limited" holds at the very short and the very long end, and
  fails in between. A floor at t = 0.1, where chi = 256 is plainly enough, is not
  truncation -- it points at integration/time-step error or the initial-state
  representation. Their release does not vary the time step, so it cannot settle
  which; that is the one thing a new run must do, and it is a dt scan at fixed
  chi in TeNPy rather than a new TDVP implementation.
"""
from __future__ import annotations
import argparse, json, pathlib
import numpy as np

HERE = pathlib.Path(__file__).parent
DATA = HERE.parent / "refs" / "fermionic_dynamics"
OBS = "spin_correlator_neighbours"
CHIS = (256, 512, 1024, 2048)
REPORT_T = (0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0)


def load(u="U_0", ref="Exact"):
    """The reference is `Exact`, NOT `FLO`.

    The deposit ships three rows under `exact`: `Exact`, `FLO` and
    `Majorana Propagation`. The first two disagree with each other by up to
    2.9e-2 at t = 2 -- the same size as the TDVP errors being measured -- so
    which one is used is not a detail. `calibration/free_fermion.py --deposit`
    settles it independently: our polynomial estimator reproduces `Exact` to
    3e-10 and `FLO` only to 2.9e-2, so `Exact` is the exact one.
    """
    import pandas as pd
    p = DATA / u / "exp_vals.h5"
    if not p.exists():
        raise SystemExit(f"missing {p}\n"
                         f"download Zenodo 17799843 and extract {u}/ there")
    t = pd.read_hdf(p, "TDVP")
    e = pd.read_hdf(p, "exact")
    return t[t.obs_type == OBS], e[(e.obs_type == OBS) & (e.method == ref)]


def dimer_links(exact):
    """The initial triplet covering: |C| = 1 at t = 0, everything else 0."""
    z = exact[exact.ev_time == 0.0]
    return {r.site_coords for r in z.itertuples() if abs(r.obs_val) > 0.5}


def errors(tdvp, exact, links=None):
    out = {}
    for tt in sorted(set(tdvp.ev_time) & set(exact.ev_time)):
        ex = {r.site_coords: r.obs_val for r in exact[exact.ev_time == tt].itertuples()
              if links is None or r.site_coords in links}
        if not ex:
            continue
        row = {}
        for c in CHIS:
            tv = tdvp[(tdvp.ev_time == tt) & (tdvp.method == f"TDVP chi={c}")]
            d = [abs(r.obs_val - ex[r.site_coords]) for r in tv.itertuples()
                 if r.site_coords in ex]
            if d:
                row[c] = float(np.mean(d))
        if len(row) == len(CHIS):
            row["signal"] = float(np.mean([abs(v) for v in ex.values()]))
            row["slope"] = float(np.polyfit(np.log(CHIS),
                                            np.log([row[c] for c in CHIS]), 1)[0])
            out[round(tt, 3)] = row
    return out


def main(as_json=False):
    tdvp, exact = load("U_0")
    links = dimer_links(exact)
    print(f"U = 0, observable {OBS}, reference method 'Exact', {len(links)} "
          f"dimer links of {len(exact[exact.ev_time == 0.0])} ordered pairs\n")
    out = {}
    # despite the name, `spin_correlator_neighbours` carries all 28*27 ORDERED
    # site pairs, not only neighbours -- 756 rows, of which 26 are the dimers
    for label, sel in (("DIMER LINKS (what the model stores)", links),
                       ("ALL 756 ORDERED SITE PAIRS (the name is misleading)",
                        None)):
        err = errors(tdvp, exact, sel)
        out["dimer" if sel else "all"] = err
        print(label)
        print(f"  {'t':>5}" + "".join(f"{c:>12}" for c in CHIS)
              + f"{'slope':>8}{'signal':>9}   verdict")
        for tt in REPORT_T:
            if tt not in err:
                continue
            r = err[tt]
            v = ("floored" if r["slope"] > -0.05 else
                 "stalling" if r["slope"] > -0.15 else "truncation-limited")
            print(f"  {tt:>5.2f}" + "".join(f"{r[c]:>12.3e}" for c in CHIS)
                  + f"{r['slope']:>8.3f}{r['signal']:>9.4f}   {v}")
        print()

    d = out["dimer"]
    print("The stored TDVP_FLOOR is reproduced from the full release:")
    print("  " + ", ".join(f"{c}: {d[0.1][c]:.2e}" for c in CHIS))
    late = [tt for tt in d if tt >= 1.5]
    print(f"\nAt t >= 1.5 the error is {np.mean([d[t][2048] for t in late]):.3f} "
          f"against a signal of {np.mean([d[t]['signal'] for t in late]):.3f} "
          f"-- comparable, so the published run says nothing about what a "
          f"converged tensor network could do there.")
    if as_json:
        p = HERE / "data" / "tdvp_check.json"
        p.write_text(json.dumps(out, indent=1))
        print(f"\nwrote {p}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    main(ap.parse_args().json)

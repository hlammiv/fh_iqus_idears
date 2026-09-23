#!/usr/bin/env python3
"""Render scenario bands and verify they are actually VISIBLE.

    python3 check_figure.py        # exit 1 if any claimed band draws nothing

Second-pass review #9 asked for rendered checks, not just logic checks: plot
synthetic bands with one zero edge and confirm the uncertainty stays visible,
cover the all-infeasible and partially feasible cases, and make sure every band
the figure claims has an identifiable visual encoding.

`selftest` checks the clipping arithmetic. This checks that matplotlib puts ink
on the page, by rendering to an in-memory buffer and measuring the filled area
of each PolyCollection. A band that has been NaN-masked away has zero area, and
that is precisely the failure the review found: on the default configuration the
NISQ+PEC scenario span had a zero lower edge at every plotted n, so nothing was
drawn where the model is most pessimistic.
"""
from __future__ import annotations
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fhcost import curves
from fhcost.budget import DEFAULT

FAILS: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not cond:
        FAILS.append(name)


def filled_area(ax) -> float:
    """Total area of every filled polygon on the axes, in axes coordinates."""
    tot = 0.0
    for coll in ax.collections:
        for path in coll.get_paths():
            v = path.vertices
            if len(v) < 3:
                continue
            v = ax.transData.transform(v)
            x, y = v[:, 0], v[:, 1]
            tot += abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))) / 2.0
    return tot


def render(lo, hi, floor=2.5, hatch=True):
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.set_xscale("log"); ax.set_yscale("log")
    n = np.logspace(2, 8, len(lo))
    ax.set_xlim(1e2, 1e8); ax.set_ylim(floor, 1e3)
    ld, hd, off = curves.band_for_plot(lo, hi, floor)
    ax.fill_between(n, ld, hd, alpha=0.3, lw=0)
    if hatch and off.any():
        ax.fill_between(n, floor, np.minimum(hd, curves.M_FLOOR), where=off,
                        alpha=0.3, lw=0, hatch="///")
    a = filled_area(ax)
    plt.close(fig)
    return a, off


def main() -> int:
    N = 40
    print("rendered band visibility\n")

    # 1. one zero edge -- the review's case
    lo = np.zeros(N)
    hi = np.linspace(10.0, 20.0, N)
    a, off = render(lo, hi)
    check("a band with a zero lower edge renders visible ink",
          a > 100.0 and off.all(),
          f"filled area {a:.0f} px^2, off-scale at all {N} points")

    # 2a. all infeasible -- nothing should be drawn, and that is correct
    a0, off0 = render(np.zeros(N), np.zeros(N))
    check("an all-infeasible band draws nothing (the figure annotates instead)",
          a0 == 0.0 and not off0.any(), f"filled area {a0:.0f}")

    # 2b. partially feasible
    lo2 = np.where(np.arange(N) < N // 2, 0.0, 6.0)
    a2, off2 = render(lo2, np.linspace(8.0, 40.0, N))
    check("a partially feasible band renders both halves",
          a2 > 100.0 and off2[: N // 2].all() and not off2[N // 2:].any(),
          f"filled area {a2:.0f} px^2, {off2.sum()} of {N} off-scale")

    # 2c. the regression: without clipping, case 1 disappears
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e2, 1e8); ax.set_ylim(2.5, 1e3)
    n = np.logspace(2, 8, N)
    ax.fill_between(n, np.where(lo > 0, lo, np.nan),
                    np.where(hi > 0, hi, np.nan), alpha=0.3, lw=0)
    a_old = filled_area(ax)
    plt.close(fig)
    check("and NaN-masking really did erase it, which is the bug",
          a_old == 0.0 and a > 0.0,
          f"masked {a_old:.0f} px^2 against clipped {a:.0f} px^2")

    # 3. every band on the real figure has ink
    ng = np.logspace(2, np.log10(3e8), 60)
    B = curves.evaluate_band(ng, DEFAULT)
    for k in ("nisq_pec", "star", "surface", "pinnacle"):
        lo_k, hi_k = B[k]
        ak, offk = render(lo_k, hi_k)
        check(f"scenario band {k!r} has visible extent",
              ak > 10.0,
              f"filled area {ak:.0f} px^2, {int(offk.sum())}/{len(ng)} off-scale")

    print()
    if FAILS:
        print(f"{len(FAILS)} FAILED: " + ", ".join(FAILS))
        return 1
    print("all bands render")
    return 0


if __name__ == "__main__":
    sys.exit(main())

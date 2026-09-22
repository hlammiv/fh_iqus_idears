# VENDORED COPY -- do not edit the original.
# source: /home/hlamm/Desktop/QC/nvidia/qlx_su3_estimate/paper/resource_model.py
# copied: 2026-09-22  (for fh_iqus; original is read-only)

#!/usr/bin/env python3
"""Self-contained fault-tolerant resource model for the SU(3) discrete-subgroup cost paper.

Every headline number is produced here from published, code-calibrated inputs, independent of any
compiler. The CUDA-Q QLX (alpha) tool is used ONLY to validate the Clifford+T gate counts (~2%).

CORRECTED IDLE MODEL (floor-respecting).
A distance-d code's sub-threshold logical error cannot decay faster than (p/p_th)^(d/2): >= d/2
faults already cause a logical error, so the combinatorial (d/2)-fault term is always present.
We therefore use exponent d_circ/2 (circuit distance), NOT d_circ, for bivariate-bicycle (BB) codes.

INPUTS AND SOURCES
  Storage / logical qubit:
    surface (rotated): 2 d^2 - 1.
    BB [[n,k,d]]: block = 2n physical (data + X- + Z-checks); per-logical = 2n/k.
    Hamiltonian sim holds the whole lattice: N_log = d_sp L^{d_sp} w  (w=8 Sigma36x3, 9 Sigma72x3);
    L=10 -> 24,000; L=2 -> 192 (216 for Sigma72x3).
  Idle logical error / logical qubit / cycle:
    surface: 0.1 (p/1e-2)^(d/2)                             [Fowler et al. 2012].
    BB: C (p/7e-3)^(d_circ/2), with C fixed so the gross code [[144,12,12]] (d_circ=10) matches the
        MEASURED memory rate ~2e-7 per 12-logical block => 1.7e-8 per logical at p=1e-3
        [Bravyi et al. 2024, arXiv:2308.07915]. d_circ: gross 10, two-gross [[288,12,18]] 18,
        big [[360,12,24]] 24. Pseudo-threshold ~0.7%.
    Decoder band: BP-OSD -> correlated BP shifts LER ~1-2 orders [Maan et al. arXiv:2510.14060] --
        a two-sided sensitivity (the extrapolation to unsimulated p is the larger uncertainty).
    [[288,50,8]] (arXiv:2606.02418): no circuit-level idle data; modeled optimistically surface-style d=8.
  Erasure conversion: p -> p/1.7 (BB circuit-level boost with erasure-aware BP-OSD
    [Pecorari & Pupillo arXiv:2502.20189]; the surface 4.4x [Wu et al. 2201.03540] does NOT transfer).
  Syndrome rounds: N_rounds = N_steps * r. r (rounds per Trotter step) = (per-link logical T-depth ~6e4)
    x (rounds per magic-state consume ~ few..d) x (link-serialization factor). SERIAL links -> r~7e7;
    fully PARALLEL links (all d_sp L^d_sp evolved concurrently, needing a ~L^d_sp-fold larger factory to
    feed them) -> r~1e6. The idle budget scales LINEARLY with r, and factory provisioning sets the
    achievable parallelism -- this is the least-pinned input. Full scale N_t=50 (serial) -> N_rounds~3.5e9.
  Adequacy budget: cumulative logical error over the whole run < 1e-2.
"""
from __future__ import annotations
import math

BUD, ERA = 1e-2, 1.7
PTH_BB, PTH_S = 7e-3, 1e-2
C = 1.7e-8/((1e-3/PTH_BB)**(10/2))          # anchor: gross d_circ=10 -> 1.7e-8/logical @1e-3
RPPS = 7e7                                   # rounds per Trotter step (factory-limited; floor ~1.6e6)

def idle_bb(dc, p):   return C*(p/PTH_BB)**(dc/2.0)
def idle_surf(d, p):  return 0.1*(p/PTH_S)**(d/2.0)
def surf_perlog(d):   return 2*d*d - 1
def bb_perlog(dc):    return 2*(144 + 18*(dc-12))/12.0     # BB [[n,12,d]]: n ~ 144+18(d-12)

CHAR_BB = {"gross [[144,12,12]]": (10, 24.0), "two-gross [[288,12,18]]": (18, 48.0),
           "big BB [[360,12,24]]": (24, 60.0)}            # name: (d_circ, per-logical)

def _fit(idle_fn, budget, p, lo, hi, step):
    d = lo
    while d <= hi:
        if idle_fn(d, p) < budget: return d
        d += step
    return None

# =============================================================== full scale
def full_scale(nlog=24000, nt=50):
    N = nt*RPPS
    imax = BUD/(nlog*N)
    print("="*74)
    print(f"FULL SCALE  {nlog:,} logical, N_rounds={N:.1e} (N_t={nt} x {RPPS:.0e})")
    print(f"  idle budget < {imax:.1e} per logical/cycle\n")
    ds = _fit(idle_surf, imax, 1e-3, 15, 60, 2)
    print(f"  surface: d={ds}, {nlog*surf_perlog(ds):,} qubits, idle_cum={nlog*N*idle_surf(ds,1e-3):.1e}")
    dc = _fit(idle_bb, imax, 1e-3, 12, 50, 2)
    print(f"  BB packing target: d_circ={dc}, per-logical {bb_perlog(dc):.0f}, "
          f"{nlog*bb_perlog(dc):,.0f} qubits -> {nlog*surf_perlog(ds)/(nlog*bb_perlog(dc)):.0f}x vs surface")
    print("  characterized BB codes at this scale:")
    for name,(d,_pl) in CHAR_BB.items():
        cum = nlog*N*idle_bb(d,1e-3)
        print(f"    {name:24} d_circ={d:>2}  idle_cum={cum:.1e}  [{'ADEQUATE' if cum<BUD else 'FAILS idle'}]")

# =============================================================== near-term spec sheet
def spec_sheet(steps=10, nlog=192, fac=3364):
    N = steps*RPPS
    imax = BUD/(nlog*N)
    wall_h = N*250e-6/3600.0
    print("="*74)
    print(f"NEAR-TERM SPEC SHEET: {steps} Trotter steps, {nlog} data qubits, error<{BUD:g}, "
          f"erasure x{ERA}")
    print(f"  required idle < {imax:.1e}; wall-clock (atoms 250us) ~ {wall_h:.0f} h\n")
    def req_p(idle_fn, d):
        lo, hi = 1e-6, 2e-2
        for _ in range(200):
            m = math.sqrt(lo*hi)
            if idle_fn(d, m) > imax: hi = m
            else: lo = m
        return math.sqrt(lo*hi)
    rows = [("surface", idle_surf, 31, nlog*surf_perlog(31)),
            ("big BB [[360,12,24]]", idle_bb, 24, math.ceil(nlog/12)*720),
            ("two-gross [[288,12,18]]", idle_bb, 18, math.ceil(nlog/12)*576)]
    print(f"  {'memory':24}{'d':>4}{'req p (raw)':>13}{'vs 2.5e-3':>11}{'phys q':>10}")
    for name, fn, d, stor in rows:
        praw = req_p(fn, d)*ERA
        print(f"  {name:24}{d:>4}{praw:>13.1e}{2.5e-3/praw:>10.1f}x{stor+fac:>10,}")
    print("  factory alternative (hold p=2.5e-3, cut idle-waiting rounds):")
    pe = 2.5e-3/ERA
    for name,(d,_pl) in list(CHAR_BB.items())[1:]:
        need = BUD/(nlog*steps*idle_bb(d,pe))
        print(f"    {name:24} rounds/step < {need:.1e}  ({RPPS/need:.0f}x more factory throughput)")

if __name__ == "__main__":
    full_scale(); print(); full_scale(nlog=27000)   # Sigma(72x3) w=9
    print(); spec_sheet()

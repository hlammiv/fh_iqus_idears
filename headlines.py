#!/usr/bin/env python3
"""Single source of truth for every headline number, so none goes stale.

README and METHODS quote from here. Regenerate with `python3 headlines.py`.
"""
import sys, pathlib, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from fhcost.budget import DEFAULT as D
from fhcost import nisq, ftqc, classical, curves, converged
from fhcost.hubbard import error_ledger

fow = D.but(pl_model="fowler")
ext = D.but(trotter="extensive")
lo, hi = classical.band(D)["band"]
s = curves.summary()

def decades(a, b):
    return f"{math.log10(b / a):.0f} decades"

L = []
L.append(f"- Classical (ED) frontier: **m = {lo}-{hi}**, an estimated capacity of")
L.append(f"  specified methods, not an impossibility boundary. The tensor-network arm is")
L.append(f"  unbounded by available data (`OPEN_ITEMS.md` O10).")
L.append(f"- Mitigated NISQ saturates: m = {nisq.max_m(1e3, D, 'pec'):.1f} at n = 10^3 to "
         f"{nisq.max_m(1e8, D, 'pec'):.1f} at n = 10^8 -- "
         f"**{decades(1e3, 1e8)} of qubits buy "
         f"{nisq.max_m(1e8, D, 'pec') / nisq.max_m(1e3, D, 'pec'):.2f}x in m.**")
L.append(f"- At n = 10^6: NISQ+PEC {nisq.max_m(1e6, D, 'pec'):.0f}, "
         f"STAR {ftqc.max_m_star(1e6, D):.0f}, "
         f"surface FT {ftqc.max_m_surface(1e6, fow):.0f}, "
         f"Pinnacle {ftqc.max_m_pinnacle(1e6, D):.0f}.")
L.append(f"- Clearing the ED frontier: Pinnacle {curves.fmt_crossing(s['n_pinnacle_clears_classical_hi'])}, "
         f"STAR {curves.fmt_crossing(s['n_star_clears_classical_hi'])}, "
         f"NISQ+PEC {curves.fmt_crossing(s['n_pec_clears_classical_hi'])}, "
         f"surface FT {curves.fmt_crossing(s['n_ft_clears_classical_hi'])}.")
L.append(f"- Continuous m is a capacity proxy. At n = 10^6 the feasible integer lattices are "
         f"NISQ {curves.max_integer_L(nisq.max_m(1e6, D, 'pec'))}x"
         f"{curves.max_integer_L(nisq.max_m(1e6, D, 'pec'))}, "
         f"FT {curves.max_integer_L(ftqc.max_m_surface(1e6, fow))}x"
         f"{curves.max_integer_L(ftqc.max_m_surface(1e6, fow))}.")
_pe = nisq.max_m(1e6, ext, "pec")
L.append("- Trotter: under the commutator bound mitigated NISQ does not reach even a"
         " 2x2 lattice; under the exact-diagonalisation calibration it reaches"
         f" {nisq.max_m(1e6, D, "pec"):.0f}."
         " Curves are drawn as a band between the two."
         if _pe <= 0 else
         f"- Trotter: the bound is {nisq.max_m(1e6, D, "pec") / _pe:.1f}x conservative in m.")
L.append(f"- Converged answers: the finite-size buffer needs {converged.m_required(0.0, D):.0f}"
         f" sites at t=0, so NOTHING classical converges (ED holds {classical.ed_frontier(D)}),"
         " and neither do NISQ or STAR. Only surface FT and Pinnacle do, by n = 1e8.")
L.append(f"- Error ledger sums to {error_ledger(D)['TOTAL']:.2f} of the tolerance, at "
         f"per-time two-sided 95%.")

txt = "\n".join(L) + "\n"
pathlib.Path("HEADLINES.md").write_text("# Headline numbers (generated)\n\n"
    "Produced by `headlines.py` from the current model. Do not hand-edit; do not\n"
    "quote figures elsewhere without regenerating.\n\n" + txt)
print(txt)

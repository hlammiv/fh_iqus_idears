# Headline numbers (generated)

Produced by `make_record.py` from `RESULTS.json`. Do not hand-edit; do not
quote figures elsewhere without regenerating.

- Model `9d5be6177a7a`. Every number below is read out of `RESULTS.json`; regenerate with `python3 make_record.py`.
- Classical (ED) frontier: **m = 24-26**, an estimated capacity of specified methods, not an impossibility boundary. The tensor-network arm is unbounded by available data (`OPEN_ITEMS.md` O10).
- Mitigated NISQ saturates: m = 10.1 at n = 10^3 to 14.9 at n = 10^8 -- **5 decades of qubits buy 1.47x in m.**
- At n = 10^6: NISQ+PEC 13, STAR 22, surface FT 13, Pinnacle 27.
- Clearing the ED frontier: Pinnacle 9.61e+05, STAR 5.39e+06, NISQ+PEC not reached below n = 1e+10, surface FT 8.07e+06.
- Continuous m is a capacity proxy. At n = 10^6 the feasible integer lattices are NISQ 3x3, FT 3x3.
- Trotter: under the commutator bound mitigated NISQ does not reach even a 2x2 lattice; under the exact-diagonalisation calibration it reaches 13. Curves are drawn as a band between the two.
- Converged answers: at t_max = sqrt(m)/v the Lieb-Robinson buffer is self-defeating -- the best arm at n = 10^6 reaches m = 27, whose own t_max = 2.59 would need m >= 1669 to CERTIFY a thermodynamic-limit answer. Nothing on the figure closes that loop, classical included. A bound failing is not proof of non-convergence (`OPEN_ITEMS.md` O8).
- Error ledger sums to 1.00 of the tolerance, at per-time two-sided 95%.

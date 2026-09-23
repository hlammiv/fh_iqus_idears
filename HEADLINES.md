# Headline numbers (generated)

Produced by `make_record.py` from `RESULTS.json`. Do not hand-edit; do not
quote figures elsewhere without regenerating.

- Model `5f032fcd3a4f`. Every number below is read out of `RESULTS.json`; regenerate with `python3 make_record.py`.
- Classical (ED) frontier: **m = 24-26**, an estimated capacity of specified methods, not an impossibility boundary. The tensor-network arm is unbounded by available data (`OPEN_ITEMS.md` O10).
- Mitigated NISQ saturates: m = 10.1 at n = 10^3 to 14.9 at n = 10^8 -- **5 decades of qubits buy 1.47x in m.**
- At n = 10^6: NISQ+PEC 13, STAR 22, surface FT 13, Pinnacle 27 -- but Pinnacle is costed on a chip with TWO COUPLER LAYERS, which its generalised bicycle codes require and the slide's nearest-neighbour grid does not provide. On that grid it is 0: not small, unavailable.
- Clearing the ED frontier: Pinnacle 9.61e+05, STAR 5.39e+06, NISQ+PEC not reached below n = 1e+10, surface FT 8.07e+06.
- Continuous m is a capacity proxy, and not every integer lattice can hold the state: half filling at S^z_tot = 0 with one holon, one doublon and a perfect triplet covering needs an EVEN site count, so floor(sqrt(m)) can name a lattice that does not exist. At n = 10^6 the largest admissible lattices are NISQ 3x4, STAR 4x5, FT 3x4, Pinnacle 4x6 (aspect capped at 2; beyond that a ribbon is quasi-1D).
- Trotter: under the commutator bound mitigated NISQ does not reach even a 2x2 lattice; under the exact-diagonalisation calibration it reaches 13. Curves are drawn as a band between the two.
- Converged answers: at t_max = sqrt(m)/v the Lieb-Robinson buffer is self-defeating -- the best arm at n = 10^6 reaches m = 27, whose own t_max = 2.59 would need m >= 1669 to CERTIFY a thermodynamic-limit answer. Nothing on the figure closes that loop, classical included. A bound failing is not proof of non-convergence (`OPEN_ITEMS.md` O8).
- Error ledger sums to 1.00 of the tolerance, at per-time two-sided 95%.

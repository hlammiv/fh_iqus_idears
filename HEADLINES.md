# Headline numbers (generated)

Produced by `headlines.py` from the current model. Do not hand-edit; do not
quote figures elsewhere without regenerating.

- Classical (ED) frontier: **m = 24-26**, an estimated capacity of
  specified methods, not an impossibility boundary. The tensor-network arm is
  unbounded by available data (`OPEN_ITEMS.md` O10).
- Mitigated NISQ saturates: m = 12.0 at n = 10^3 to 17.6 at n = 10^8 -- **5 decades of qubits buy 1.47x in m.**
- At n = 10^6: NISQ+PEC 16, STAR 27, surface FT 26, Pinnacle 48.
- Clearing the ED frontier: Pinnacle 1.57e+05, STAR 8.58e+05, NISQ+PEC not reached below n = 1e+10, surface FT 1.13e+06.
- Continuous m is a capacity proxy. At n = 10^6 the feasible integer lattices are NISQ 3x3, FT 5x5.
- Trotter: under the commutator bound mitigated NISQ does not reach even a 2x2 lattice; under the exact-diagonalisation calibration it reaches 16. Curves are drawn as a band between the two.
- Converged answers: the finite-size buffer needs 2 sites at t=0, so NOTHING classical converges (ED holds 26), and neither do NISQ or STAR. Only surface FT and Pinnacle do, by n = 1e8.
- Error ledger sums to 1.00 of the tolerance, at per-time two-sided 95%.

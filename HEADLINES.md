# Headline numbers (generated)

Produced by `headlines.py` from the current model. Do not hand-edit; do not
quote figures elsewhere without regenerating.

- Classical (ED) frontier: **m = 24-26**, an estimated capacity of
  specified methods, not an impossibility boundary. The tensor-network arm is
  unbounded by available data (`OPEN_ITEMS.md` O10).
- Mitigated NISQ saturates: m = 13.0 at n = 10^3 to 18.9 at n = 10^8 -- **5 decades of qubits buy 1.45x in m.**
- At n = 10^6: NISQ+PEC 17, STAR 29, surface FT 32, Pinnacle 57.
- Clearing the ED frontier: Pinnacle 1.10e+05, STAR 4.29e+05, NISQ+PEC not reached below n = 1e+10, surface FT 6.81e+05.
- Continuous m is a capacity proxy. At n = 10^6 the feasible integer lattices are NISQ 4x4, FT 5x5.
- Trotter: under the commutator bound mitigated NISQ does not reach even a 2x2 lattice; under the exact-diagonalisation calibration it reaches 17. Curves are drawn as a band between the two.
- Converged answers: the finite-size buffer needs 144 sites at t=0, so NOTHING classical converges (ED holds 26), and neither do NISQ or STAR. Only surface FT and Pinnacle do, by n = 1e8.
- Error ledger sums to 1.00 of the tolerance, at per-time two-sided 95%.

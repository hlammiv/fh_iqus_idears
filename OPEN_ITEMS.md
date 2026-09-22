# Open items

Deferred questions found while working through the scientific review. Separate
from the review's own list, which is tracked below.

## Deferred

### O1. Is the NISQ-vs-STAR sampling comparison apples-to-apples?
*Raised by finding #3. Partially resolved; the residual question is narrower than
first flagged.*

Correcting the NISQ PEC cost upward (1.875x in the exponent) moved STAR from
losing below `n ~ 5x10^4` to beating bare NISQ at every `n` it can run (1.04x at
10^4 rising to 1.94x at 10^8). The worry was that the two arms draw their
overheads from different sources.

**Checked, and the worry is mostly unfounded.** STAR's injected-rotation error is
a Z-type Pauli channel with RUS-total probability `q = 4p/15`. Its signed inverse
is `a = (1-q)/(1-2q)`, `b = -q/(1-2q)`, giving `gamma = 1/(1-2q)` and
`2 ln(gamma)/p = 1.0670` — which reproduces the STAR paper's
`gamma^2 = exp(8 P_Z,1 N)` exactly (`8 x 2/15 = 1.0667`). Verified against the
channel's Pauli transfer matrix. So **both arms are cancellation one-norms**, not
one of each.

What is left to check:
1. NISQ pays the one-norm over *all* two-qubit gates in the cone (`G_cone`);
   STAR pays it only over *rotations* (`n_rot_cone`), its Cliffords being
   error-corrected. That is a real architectural difference, but the margin is
   then set almost entirely by the assumed ratio `c_rot/c_g = 5/15`. That ratio
   has never been checked against a compiled circuit.
2. Confirm the STAR Clifford layer really contributes **zero** sampling overhead
   and only residual logical error (which is separately budgeted through `p_L`
   and the code distance).

### O2. Signal parameters are placeholders pending a fit
`s_res_min = 0.02` and `s_res_slope = 0.10` are estimates; the residual
antiferromagnetic correlation is not quoted numerically in arXiv:2510.26300.
A fit to their Zenodo deposit (record 17799843) is outstanding. `s_short = 1.0`
and `t_melt0 ~ 0.4-0.7` are sourced.

### O3. Pinnacle code family exhausts above n ~ 1e10
Only five published generalised bicycle codes, topping out at `d = 24`, so the
curve saturates artificially beyond that. Model limitation, not physics. The
plotted range (`n <= 1e8`) is inside the valid region.

### O4. The ZNE response model is unvalidated
Every ZNE conclusion rests on `s(lambda) = s(0) exp(-lambda Lambda)`. A response
with less curvature raises every bias ceiling. Measuring the real noise response
of a small compiled Hubbard circuit would settle it.

## Review status

| # | finding | status |
|---|---|---|
| 1 | initial state and connected correlator incompatible | **fixed** — matched to arXiv:2510.26300 |
| 2 | ZNE feasibility omits residual bias | **fixed** — bias-limited, not variance-limited |
| 3 | gatewise PEC overhead uses the wrong coefficient | **fixed** — attenuation separated from one-norm |
| 4 | mitigation theorem is overstated | open |
| 5 | multiproduct gains lack an error bound | open |
| 6 | Pinnacle calibration needs reconstruction | open |
| 7 | FT resource and error accounting incomplete | open |
| 8 | classical band is heuristic, not a ceiling | open |
| 9 | light-cone geometry inconsistent (1/3 vs 2/3) | open |
| 10 | error components do not combine to the tolerance | open |
| 11 | implementation and reporting issues | open |

# Trotter-error calibration

The model's central empirical input. Everything needed to regenerate it is here.

```
python3 calibration/trotter_cal.py [patch ...]   # produces data/*.json
python3 calibration/fit_w.py --holdout           # produces hubbard.W_MEASURED
```

`fit_w.py` reproduces the hard-coded `W_MEASURED` dict exactly. If it ever stops
doing so, the dict is stale.

## What was computed

Exact evolution of the dimerised S^z_tot = 0 triplet state, measuring the
equal-time connected nearest-neighbour correlator
`C^zz_ij = 4(<S^z_i S^z_j> - <S^z_i><S^z_j>)` on the first dimer link, against
second-order Trotter at step counts r = 1..64.

| convention | value |
|---|---|
| lattice / ordering | `rectangle(nx,ny)` boustrophedon snake, matching the collaborator's `small_dynamics.py` |
| patches | `square2` (4), `cross5` (5), `rectangle2x3` (6), `square3` (9), `rectangle3x4` (12); `square4` (16) still running |
| boundaries | **open** (no flux, no periodicity) -- the experiment uses a doubly periodic 7x4 torus with Phi = pi, so this differs |
| sector | fixed (N_up, N_down), half filling, one holon + one doublon, dimers on the rest |
| hopping split | two colour classes by `x%2` / `y%2`, as in `small_dynamics.py` |
| Trotter step | symmetric 5-stage `onsite(dt/2) pink(dt/2) gold(dt) pink(dt/2) onsite(dt/2)` |
| U/J | 0, 4, 8 |
| times | tau = 0.25, 0.5, 1.0, 2.0 |
| reference | Lanczos exponential, validated against `scipy.expm_multiply` to 2e-16 |
| sanity | `C^zz(0) = -1.000000` on every patch, matching the triplet algebra |

## Domain, and how far the model runs past it

**Calibrated:** 4-12 sites, tau 0.25-2.0, U/J 0-8.

**Every operating point on the figure is outside it.** At n = 10^8 the surface-FT
point sits at m = 236 and tau = 7.7; Pinnacle at m = 383 and tau = 9.8. The
m-extrapolation is argued from locality (W_eff is flat across n = 6, 9, 12 at
fixed tau); the tau-extrapolation is a clamp and has no such argument.

## Known weaknesses

1. **Boundaries differ from the experiment** (open here, doubly periodic torus
   there). Not corrected.
2. **Leave-one-out stability degrades with time.** Refitting without the n = 12
   patch moves the coefficient 0-11% at tau = 0.25-0.5 but up to 25% at tau = 1
   and 82% at tau = 2 -- worst precisely where the model extrapolates hardest.
3. **The U = 0, tau = 2 entry is non-monotonic** (below U = 4), an accidental
   zero-crossing in the error. Kept unsmoothed; not physics.
4. **Only the first dimer link** is probed, not an average over links.
5. **The multiproduct remainder is not calibrated.** The model reuses the
   second-order coefficient at higher orders. The convergence ORDER is validated
   (3.5-4.1 measured for order 4, 5.7-6.5 for order 6) but the coefficient is not.

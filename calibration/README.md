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

## free_fermion.py — the U = 0 easy limit (second-pass review #7)

`C^zz` on the triplet-covering state at `U = 0`, in polynomial time, checked
against exact many-body evolution.

    python3 calibration/free_fermion.py      # ~4 minutes, < 1 GB

**What it computes.** At `U = 0` the Hamiltonian is quadratic but the initial
state is not Gaussian: a triplet covering expands into `2^{m/2}` Fock branches.
`C^zz` is weight-4, and a four-fermion operator connects branches differing in
at most one triplet, so the branch sum collapses into a diagonal part (fixed by
the one- and two-mode occupation statistics) plus a coherent part local to one
triplet. Ordering each triplet's four modes contiguously makes the outside
Jordan–Wigner strings cancel, which is what keeps the coherent part local; the
4-mode matrix elements are evaluated on a 16-dimensional Fock space rather than
by hand, because that is where sign errors live.

**Result.** Worst `|C_poly − C_exact| = 1.7e-15` over `m = 4, 6, 8, 9` and
`t = 0, 0.25, 0.5, 1, 2`. The `O(m)` form (sparse Krylov propagator,
block-collapsed sums) reproduces the `O(M²)` form to 1e-16 and measures a
log–log exponent of 0.979 out to `m = 16384`.

### Memory

Two guards, both added after hitting them:

* `hopping()` refuses a dense `2m × 2m` above 8192 modes — it is 128 GiB at
  `m = 65536`. Use `hopping_sparse()`.
* `occupations()` refuses a dense `M × M` `nu2` above 8192 modes — 8.6 GB at
  `m = 16384`. The `O(m)` path uses `nu_only()`.

The exact reference builds a fixed-particle-number Fock sector and uses
`expm_multiply`, so `m = 9` (3×3, 18 modes) costs ~2 s and a few hundred MB.
Anything larger belongs on `lenore_remote`.

### Known weaknesses

1. The measured wall time is single-core unoptimised CPython. It is used as a
   *conservative* frontier precisely because it needs no extrapolation, but it
   is 2–3 orders of magnitude off what a vectorised implementation would do.
2. Exact at `U = 0` only. No accuracy claim is made at small non-zero `U`; a
   controlled expansion in `Ut` is the obvious next step and is not done.
3. The geometry is a plain `Lx × Ly` torus with a dimer covering along x, not
   the experiment's 7×4 double-periodic lattice at flux `Φ = π`. The validation
   compares two calculations on the *same* geometry, so this does not affect the
   correctness claim, only the specific numbers.

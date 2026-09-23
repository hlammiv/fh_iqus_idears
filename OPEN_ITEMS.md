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

### O2. Signal parameters — **fitted**, with residual caveats
Fitted to Zenodo 17799843 (dimer-link C^zz, TFLO+GPR, both U). The decay is
**Gaussian**, not exponential (free stretch exponent 2.15 +- 0.09 and
1.98 +- 0.05). `s_short = 1.000 +- 0.004` confirms the triplet algebra.
The placeholder `s_res_slope = 0.10` was **5x too large** (fitted 0.021).

Remaining caveats:
1. **Only two U values.** `s_res_slope` and `t_melt_slope` are a two-point
   interpolation with zero constraint on curvature. No basis for extrapolating
   beyond U/J = 4; the linear form is an assumption.
2. **`s_res` is not a real plateau.** |C| oscillates with period ~1.0-1.1 J^-1 and
   an amplitude comparable to the residual itself (U=0: |C| spans 0.021-0.077
   over t in [1,2], a +-50% swing about 0.043). A single scalar is a coarse
   summary; the honest form is Gaussian decay plus a damped revival.
3. **No per-link exact reference at U=4** in the deposit, so `s_res(U=4)` rests on
   mitigated data validated only globally (via n_triplets vs Majorana
   propagation, which agree to ~3%).
4. Do **not** use TDVP as ground truth for the residual: even chi=2048 gives
   0.125 at U=0 against 0.043 exact, a 3x overestimate.

### O3. Pinnacle code family exhausts above n ~ 1e10
Only five published generalised bicycle codes, topping out at `d = 24`, so the
curve saturates artificially beyond that. Model limitation, not physics. The
plotted range (`n <= 1e8`) is inside the valid region.

### O4. The ZNE response model is unvalidated
Every ZNE conclusion rests on `s(lambda) = s(0) exp(-lambda Lambda)`. A response
with less curvature raises every bias ceiling. Measuring the real noise response
of a small compiled Hubbard circuit would settle it.

### O5. Which gates damp the observable — measured, but the scaling is not
*Raised while working #4. Partially supersedes review #9.*

Measured on the Phasecraft/Quantinuum circuit (2415 two-qubit gates, their quoted
`p = 1e-3`): observed attenuation `Lambda = 0.20`, against **2.58** from the cone
model — a **12.9x overcharge**. The mechanism was checked on their own data by a
weight test: `Lambda(ZZ)/Lambda(Z) = 1.91`, where damping proportional to operator
weight predicts 2 and uniform-per-gate predicts 1. So a depolarizing error damps
the observable only where the Heisenberg-evolved operator has support.

`damping_model = "support"` implements this and reproduces the measurement to 8%
(0.184 vs 0.200). **It is not the default**, because the single measurement
constrains the *value* at one operating point but not the *scaling*:
`support_growth = 0` is a one-point fit and cannot be right at long times, when
the operator must eventually fill the lattice. Adopting it by default would be
extrapolating an unvalidated law in the direction that flatters the results.

Effect where it does apply: the discount goes as `w_obs/q ~ 1/m`, so it is large
for a small observable in a big register (their case, 13x) and modest at the m
this model actually reaches (1.21x at `n = 1e6`).

**Relation to review #9.** The reviewer showed the cone fraction should be 2/3,
not 1/3 — a 2x correction *unfavourable* to the results. The measurement says the
cone framing is the wrong quantity and overcharges by ~13x, *favourably*. Both
cannot stand. #9 remains open for the rest of its content (the integration
geometry, cluster-truncation matching, and the velocity taxonomy); only the
question of which gates damp is addressed here, and by measurement rather than
by geometry.

### O6. The demonstrated-mitigation arm reaches nothing at our step count
`strategy = "expcal"` costs TFLO+GPR at its measured effective overhead (0.08 —
*below* one, because GPR borrows statistics across correlated time points) and is
drawn only out to the largest `Lambda` actually demonstrated, 0.20. Under this
model's workload it reaches **nothing**: `Lambda(m=4) = 2.87`. Adopt the
experiment's fixed step density `r = 4 tau` and it immediately reaches `m = 20`
at `n = 1e6`.

So the whole distance between this model and a real 56-qubit experiment is the
**Trotter step count** (134 vs 4 at m=4), not the mitigation scheme, not the
channel conventions, and not the damping geometry. That makes review #5 — which
asks for an actual error bound behind the step count — the highest-value item
remaining.

### O7. The ideal p=0 curve should be capped by m_required(t)
*Queued behind #5.*

Once the lattice exceeds `m_required(t)` the finite answer already IS the
infinite answer, so extra sites buy nothing. At fixed `t = 1` that ceiling is
244 sites — about **730 physical qubits** — while the plotted p = 0 line reaches
333,333 at `n = 10^6`, a **1370x overshoot**. Past the ceiling the qubits should
buy time, not lattice.

It does not show up under the default convention because `t_max = sqrt(m)/v_B`
ties the ceiling to the lattice: `m_req(t_max) ~ 4m` at every size, so the target
recedes twice as fast as the lattice grows and nothing is ever converged. Within
that convention the slope-1 line is technically correct — which exposes the worse
problem: **if nothing is ever converged, the y axis is not measuring a physics
answer at all**, only the largest finite-lattice calculation that fits.

Fix: draw `m_required(t)` as an explicit ceiling (a BAND — it inherits the 11x
spread from O5, 244 to 2664 at t = 1) and cap the ideal curve by it in fixed-t
mode; relabel the axis in sqrt(m) mode to say what it measures.

This is the third distinct place the time-window convention has driven a wrong
conclusion, after finite-size extrapolation buying nothing and the classical band
being a hard wall.

### O8. ~~The Trotter calibration is U/J = 4 only~~ **CLOSED**
Calibrated at U/J = 0, 4 and 8. `W_eff` spans **68x** across that range at
tau = 0.25, almost all of it in the step from U = 0 (free hopping, where the two
colour groups nearly commute and the on-site term carries all the error) to
U = 4; beyond that it flattens, consistent with the dynamics slowing as 4J^2/U.

A consequence worth noting: because the circuit gets harder with U while the
signal gets *bigger* with U, **m(U) is non-monotonic and the conventional
U/J = 4 is the worst case** (m = 58 / 29 / 37 at U = 0 / 4 / 8, n = 10^6).

One entry is unphysical: U = 8, tau = 2 sits below U = 4. That is the same
accidental zero-crossing artefact seen in the raw error there. It is kept rather
than smoothed, and should not be read as physics.

### O9. The calibration is extrapolated well past its range
Measured on **4-12 sites, tau in [0.25, 2]**. The curves use it at m up to ~2000
and tau up to ~13.

* The **m-extrapolation** is supported by locality: `W_eff` is flat across
  n = 6, 9, 12 at fixed tau (1.461 / 1.448 / 1.453 at tau = 0.25, i.e. constant to
  1%), which is what a two-site observable must do once the lattice exceeds its
  cone. m = 16 is running to widen the range from 3x to 4x.
* The **tau-extrapolation is only a clamp.** `W_eff` falls ~15x over the measured
  range and is held fixed beyond tau = 2. That is conservative against the
  measured trend, but a tau^3 law stretched 6.6x past its data is not defensible
  on its own.
* **Every plotted operating point is outside the calibrated domain** -- worst case
  m = 383 and tau = 9.8 against m <= 12 and tau <= 2. `hubbard.calibration_status`
  reports this per point and the figure marks the calibrated m.
* **Leave-one-out stability degrades with time**: refitting without the n = 12
  patch moves the coefficient 0-11% at tau = 0.25-0.5 but up to 82% at tau = 2,
  worst precisely where the extrapolation is longest.
* The two Trotter models are **SCENARIOS, not an uncertainty interval**. Neither is
  a guaranteed edge and the truth is not required to lie between them.

Because of this every arm is now reported as a **band** between Campbell's bound
and the calibration, rather than as a point value. Neither edge is the answer.

### O10. The tensor-network question is OPEN, not resolved
*Raised by #8. This is the item most likely to reverse a headline.*

The classical band is now **ED-only, m = 24-26**, down from 24-62. That narrowing
is what flips "mitigated NISQ never clears classical" into "it clears at
n ~ 1.7e4". **The band narrowed because an unvalidated arm was removed, not
because tensor networks were shown to fail** -- and removing an arm is a change in
the direction that flatters the quantum curves, which deserves suspicion.

**Sharpened by a diagnostic that cost nothing to run.** Splitting the published
TDVP error by TIME rather than averaging it:

| t | exact | chi=256 | chi=2048 | ratio |
|---|---|---|---|---|
| 0.1 | 0.9418 | 7.93e-3 | 7.35e-3 | **1.1x** |
| 0.5 | 0.2429 | 6.36e-2 | 2.56e-2 | 2.5x |
| 2.0 | 0.0211 | 1.72e-1 | 1.57e-1 | 1.1x |

At t = 0.1 the state is barely entangled and chi = 256 is wild overkill; a
truncation-limited calculation would be at machine precision. Instead the error
is 7.9e-3, and an **eightfold** increase in bond dimension removes **8%** of it.
The deposit's `max_bond_dimension` column confirms every run saturated its cap,
so truncation was binding -- it simply was not what limited the accuracy.

So the published TDVP carries a **chi-independent error floor of ~7e-3 from the
earliest times** (plausibly two-site TDVP projection error on a snake MPS with
long-range Jordan-Wigner strings, or the time step, or the GPR smoothing). It is
**not a converged tensor-network calculation**, the chi-extrapolation to 2e12 is
meaningless, and this dataset **cannot bound what tensor networks can do here**.

This also corrects the reasoning given in METHODS 5. The conclusion there --
drop the entropy-derived arm -- stands, because the entropy model cannot
represent a floor of this kind. But the argument offered for it, "the error is
flat in chi, so tensor networks fail", was wrong: flatness in chi is evidence the
run was not chi-limited, not evidence that chi cannot help.

The floor sits at ~7e-3, essentially AT our absolute tolerance of ~6e-3. **A
clean implementation that removed it could plausibly reach this accuracy at
modest chi**, which would move the classical upper edge well above 26 and could
reverse the FT/STAR/Pinnacle crossings.

What is NOT established, and would move the band back up:
1. **That was not a best-effort classical attack.** They ran TDVP as a
   comparison, not as an optimised attempt to beat their own experiment.
2. **The plateau looks systematic, not truncation-limited.** If it is TDVP
   time-step or projection error, the chi-extrapolation to 2e12 is meaningless
   and a better-converged run could do far better.
3. **The geometry is close to worst case for MPS**: a snake on a 7x4 torus,
   periodic in both directions. A cylinder or PEPS would do better.
4. **No PEPS, neural-quantum-state or Pauli-path attempt exists** for this
   observable at this accuracy.
5. The ED edge itself (24-26) is a memory estimate, never benchmarked.

**Does leadership-class scaling change this?** Checked, and the two arms answer
differently.

*ED: no.* Out-of-core looks attractive -- a 700 PB filesystem holds an m = 30
vector where 10 PB of memory does not -- but time evolution touches the whole
vector once per Krylov matvec, so the binding cost is bandwidth. At m = 28 that
is ~1300 s per matvec and **50 weeks** for the evolution; m = 30 is 837 weeks.
The in-memory limit of 26 stands, and is now enforced in code.

*Tensor networks: massively, and this is the weak point.* The published run used
chi = 2048, which is **4.8e14 flops = 0.28 MILLISECONDS of exascale compute**. A
week of exascale is **2.1e9 times more**. Within the same budget the quantum arms
are given, chi could be 3e6 -- 1500x higher. **I used a sub-millisecond
calculation as the classical frontier.** That is a category error.

On the measured slope even chi = 3e6 gives 0.032 against a 0.0064 target, so the
conclusion may survive. But a 1500x increase in chi buying only 2.4x in error is
itself evidence the error is not truncation-limited, which means the slope cannot
be extrapolated in **either** direction. The tensor-network arm is therefore
**unbounded by the available data**, not bounded and small.

**Decided not to run it here.** chi = 1e5 needs 16.7 TB for the state alone
(~50 TB with TDVP working space) against lenore's 122 GB -- short by 400x. It
needs ~100-200 leadership nodes, which is an allocation request, not a background
job. lenore's ceiling is chi ~ 4800, only 2x the published value, which would
test the measured slope over one more factor of two and predict 0.077 -> 0.072:
too weak to settle anything.

**The decisive cheap test, if this is picked up again**, is not more chi at all.
Hold chi fixed and vary the TDVP time step. If the error moves with dt at fixed
chi, the 0.077 plateau is integration error, the chi-extrapolation is void in
both directions, and no bond dimension settles the question. If it does not move,
the plateau is real and the band is on firmer ground. Either outcome is decisive
and it runs at chi ~ 2048. It needs a correct fermionic 2D TDVP on a doubly
periodic torus -- use a library (TeNPy has Fermi-Hubbard and two-site TDVP built
in); hand-rolling one risks a confidently wrong number.

Until at least one purpose-built, leadership-scale classical attack is run, the
honest headline is "clears the **exact-diagonalisation** frontier", not "clears
classical".

### O11b. Multiproduct coefficients are calibrated only to tau <= 0.5
*Tier-1 done; tier-2 rerun queued. (Second-pass finding #2.)*

The model reused the SECOND-order coefficient at every multiproduct order.
Validating the convergence order does not justify that -- the order and the
coefficient are different things. Order-2k coefficients are now extracted from
the calibration data (`calibration/fit_w.py --mpf`), r >= 4, U/J = 4:

    W_MPF = {4: {0.25: 0.118, 0.5: 0.161},
             6: {0.25: 0.031, 0.5: 0.017},
             8: {0.25: 0.007, 0.5: 0.004}}

At tau = 0.5 order 4, W_4 = 0.161 against W_2 = 1.290, so the proper coefficient
moves the step count materially: MPF at n = 1e6 gives 22 / 28 / 28 at orders
4 / 6 / 8 against 23 / 23 / 21 before, putting the interior optimum at order 6.

**Only tau <= 0.5 is usable.** At tau >= 1 the observable error passes through
zero crossings and the extracted coefficient swings 3-13x across r. Every MPF
point on the figure is beyond even this reduced domain -- a stronger caveat than
the second-order calibration carries, and `selftest` asserts it.

**Queued rerun (~1-2 hours on lenore, all patches <= 12 sites):**
1. finer r sampling (r = 1..128, not powers of two) to locate the asymptotic
   plateau and step over the zero crossings;
2. record **state infidelity** alongside the observable error -- it is monotone
   and has no zero crossings, so it pins the coefficient where the observable
   cannot;
3. more tau values, for interpolation rather than clamping.

Nothing here fixes tau > 2 or m > 12; that is O9.

## First-pass review status

| # | finding | status |
|---|---|---|
| 1 | initial state and connected correlator incompatible | **fixed** — matched to arXiv:2510.26300 |
| 2 | ZNE feasibility omits residual bias | **fixed** — bias-limited, not variance-limited |
| 3 | gatewise PEC overhead uses the wrong coefficient | **fixed** — attenuation separated from one-norm |
| 4 | mitigation theorem is overstated | **fixed** — citation corrected, claims downgraded, demonstrated arm added |
| 5 | multiproduct gains lack an error bound | **partial** — Trotter error calibrated by exact diagonalisation (alpha 9/4 -> 7/4); branch-cost accounting still open |
| 6 | Pinnacle calibration needs reconstruction | **fixed** — footprint reproduced with no free parameter; single-engine T supply corrected; connectivity flagged as O11 |
| 7 | FT resource and error accounting incomplete | **fixed** — HWP workspace charged, magic-state error budgeted with per-point factory selection, failure-to-bias factor 2 |
| 8 | classical band is heuristic, not a ceiling | **partial** — relabelled as capacity, entropy bound fixed, TDVP measured; a validated TN estimate is still missing (O10) |
| 9 | light-cone geometry inconsistent (1/3 vs 2/3) | **partial** — which-gates-damp measured (O5); geometry, cluster matching and velocities still open |
| 10 | error components do not combine to the tolerance | open |
| 11 | implementation and reporting issues | **fixed** — all eleven; headlines now generated from one record |

## Second-pass review status

| # | finding | status |
|---|---|---|
| 1 | the convergence criterion fails an exact physical limit | **fixed** — factorial Lieb-Robinson form, exact t -> 0 limit |
| 2 | the default Trotter calibration is used beyond its evidence | **partial** — order-2k coefficients measured (O11b); domain still exceeded (O9) |
| 3 | Python, the explorer and the documents disagree | **fixed** — one model, one record; see below |
| 4 | factory cleanup neglects circuit-level errors | **fixed** — published operating points replace the cubic law; selection is p-aware and plant-level |
| 5 | Pinnacle lacks the common FT resource/error ledger | **fixed** — workspace charged, engine selected and certified, stalls and rejection scheduled |
| 6 | Hamming-weight workspace does not match the cited circuit | **fixed** — Campbell Thm 2; workspace, T-count, depth and synthesis all from one construction |
| 7 | the classical baseline misses the U = 0 easy limit | **fixed** — free-fermion estimator implemented and validated to 1.7e-15; band relabelled |
| 8 | the signal fit does not support per-time relative accuracy | **fixed** — envelope separated from bound, per-time tolerances, floor independent of the fit |
| 9 | zero uncertainty edges disappear from the plot | **fixed** — clipped and hatched, scenario span labelled, rendered check added |
| 10 | integer lattice reporting ignores initial-state constraints | **fixed** — admissibility from the state spec; rectangles with an aspect cap |
| extras | ledger enforcement, inert `s_sig`, cluster prose, fixed-t vs certified, narrative tests | **fixed** (last is acknowledged, not removed) |

All ten findings and the five additional checks have been worked. Nine are
closed; #2 is partial (order-2k coefficients measured, but the calibration
domain is still exceeded -- O9 and O11b). The review's "recommended next work"
list maps onto them as follows, with what remains:

| recommended | state |
|---|---|
| 1. exact-limit semantics: t = 0, U = 0, admissible states | done (#1, #7, #10) |
| 2. reproducible evidence: calibration code, raw data, held-out checks | done (`calibration/`, four programs + data + README) |
| 3. Python/JavaScript parity and a shared result record | done (#3; verified in a headless browser on every publish) |
| 4. a common FT resource ledger | done (#4, #5, #6); Pinnacle per-candidate block packing still open |
| 5. precision at each time | done (#8); `s_data_rel_unc` and `s_abs_floor` are now the leveraged assumptions |
| 6. recompute architecture comparisons | done, but every arm moved: see HEADLINES.md, and the classical baseline is an ED capacity, not a hardness result |

### What #10 closed, and what it did not

`floor(sqrt(m))` named lattices the specified state cannot occupy. Half filling
at `S^z_tot = 0` with one holon, one doublon and a perfect triplet covering needs
an EVEN site count (N_up = N/2) and a dimer-coverable remainder; a 5x5 fails
both, leaving 23 sites. `lattice_admissible` checks it and says why;
`best_lattice` enumerates rectangles (the experiment uses 7x4) with the most
square winning and an aspect cap of 2, beyond which a ribbon is quasi-1D and a
materially easier problem.

The correction runs both ways: floor(sqrt) named an impossible lattice at four of
five sizes checked AND understated the capacity -- at m = 13 it said 3x3 = 9,
while 3x4 = 12 is both admissible and larger.

**Still open from #10:** the review also asked for compiled resource feasibility
per admissible candidate, including code-block and replica packing. The model
still costs a continuous m and reports the admissible lattice afterwards; it does
not re-cost each candidate with its own qubit count rounded to code blocks. For
Pinnacle in particular, k = 14-16 logical qubits per block means the packing is
lumpy and the honest answer per candidate could differ by a block.

### The additional implementation checks

* **Ledger enforcement** was a side effect of `eps_absolute`. It did work --
  `frac_stat = 0.9` is rejected at all six public entry points -- but
  `hubbard.validate(cfg)` is now explicit and also checks eps, p against
  threshold, n_times and the absolute floor.
* **`s_sig` was inert** under `signal_regime = "curve"`, so two sensitivity rows
  reprinted the baseline. They now vary the fixed-signal scenario and the
  parameters the curve uses.
* **The cluster prose was wrong**: `cluster_t_reach(xi=0.2) = 0.368`, not zero.
  The claim is now confined to `xi >~ 0.5`.
* **Capacity at fixed t vs certified time** are different questions;
  `max_m_fixed_t_detail` carries the method and the non-claim.
* **Tests asserting narrative** -- acknowledged. Those touched in this pass now
  test a mechanism against a measured or published number; the rest are
  regression tests and are not evidence.

### What #9 closed, and what it did not

Worse than the review found: on the default configuration the NISQ+PEC scenario
band had a zero lower edge at **every one of the 60 plotted points**, so that arm
had no shading at all -- not a gap, a complete absence. The slide figure's Willow
band lost half its range the same way.

`curves.band_for_plot` clips a zero lower edge to the axis floor and returns the
clipped region as a mask; the figures hatch it under a dashed `m = 4` line, so
the band visibly runs off the bottom. All-infeasible columns return NaN on both
edges and are annotated rather than drawn. The shading is labelled at the foot of
the panel as a **scenario span, not a confidence interval**, since this style has
no legend to put that in.

`check_figure.py` is a rendered check, not a logic check: it measures the filled
area of each polygon. The NaN-masked version renders 0 px^2 where the clipped one
renders 13402.

**Still open from #9:** the hatched strip sits between the axis floor (2.5) and
`m = 4`, which is only about 13% of a decade -- legible but cramped, and four
overlapping hatched arms at the bottom of panel (a) are busy. Lowering the axis
floor would give it room at the cost of empty space everywhere else.

### What #8 closed, and what it did not

Three things were conflated and are now separate: the fitted **envelope**
(physics), a **lower bound** the costing may assume (data-based), and an
**absolute floor** (a specification, not the fitted residual). And one tolerance
was doing two jobs; `eps_absolute` now bounds systematics at the tightest time on
the grid while `eps_statistical` carries the per-time shot sum.

The envelope failed on its own terms: a leave-one-out refit
(`calibration/fit_signal.py`, which preserves the data, fit, residuals and
covariance) misses the tail by -52% to +55%, and at U/J = 4, t = 1.5 it sits 37%
above the data -- 1.89x too few shots, exactly as the review computed.

Everything tightened by 10-50%, the binding time moved off the endpoint (t = 1.6
at m = 256, not t_max = 8), and multiproduct NISQ no longer clears the ED
frontier at any n below 1e10.

Two prior results were overturned in passing. The surface arm's drop at n = 1e6
is a **factory-ladder edge**, not an architectural finding -- it sits exactly on
the cultivation rung's boundary, one step from a 68x larger plant. And the
model's agreement with `m ~ s^(2/9)` was a **floor artefact**: the s = 0.03 point
was clipped by the floor at the fitted residual, flattening the slope to 0.275;
unclipped it is 0.379.

**Still open from #8:**
* `s_data_rel_unc = 0.20` is an assumption. The source quotes no per-point error
  bar, and this number now sets every tolerance in the model.
* `s_abs_floor = 0.02` is a specification with no derivation. It is the single
  most leveraged number in the ledger and deserves an argument from what the
  experiment can actually resolve.
* Seven data points per U, two U values. The interpolation in U between 0 and 4
  is a two-point line, and beyond t = 2 there is no data at all -- the model
  falls back to the envelope with a 0.65 factor, which the held-out test says is
  about right but cannot confirm.
* The review also asked that smoothed data points not be treated as independent.
  They still are: the published values are GPR output, so neighbouring points
  share a kernel and the residual RMS understates the true uncertainty.

### What #7 closed, and what it did not

At `U = 0` this observable is polynomial and the ED band was simply the wrong
baseline. `calibration/free_fermion.py` implements the collapse of the
non-Gaussian branch sum (diagonal part from occupation statistics, coherent part
local to one triplet once its modes are ordered contiguously) and **validates it
against exact many-body evolution to 1.7e-15** at `m = 4, 6, 8, 9`. An `O(m)`
form reproduces it to 1e-16 and measures a log-log exponent of 0.979 out to
`m = 16384`.

Measured frontier: `m ~ 4e7` on ONE core of unoptimised CPython in a week,
against 26 for ED. No quantum arm is within six orders of magnitude at that
point. The band is relabelled "estimated ED capacity" everywhere, and the
explorer's U slider now reaches 0 and says so.

The TDVP peak-FLOP arithmetic is labelled as such, with an explicit
`assumed_fraction_of_peak = 1.0`.

**Still open from #7:** the review's third test -- vary the TDVP time step AND
bond dimension against an exact small reference, to identify what causes the
chi-independent plateau -- is not done. It needs a correct fermionic 2D TDVP on
a doubly periodic torus, which is O10, not something to hand-roll. Until then
the U > 0 band stays an ED capacity with the tensor-network question open.

Nothing here touches `U > 0`. The free-fermion claim is confined to `U = 0`
exactly; no accuracy claim is made at small non-zero U, where the natural next
step would be a controlled expansion in `Ut`.

### What #6 closed, and what it did not

The workspace formula read a synthesis T-count as a register size. It is now
Campbell Thm 2's `alpha = b - w(b)`, reproducing 63 / 255 / 428 at
`b = 64 / 256 / 432`, and the T-count, depth and synthesis allocation all come
from the same batching rather than from three different places. The synthesis
allowance now covers every synthesised rotation in the shot, branches included,
instead of dividing by the five rotation groups per step.

The old formula was wrong in both directions: too large at small `m` (a bogus
~60-qubit gradient register) and far too small at large `m` (log instead of
linear). Surface FT 4.9 -> 9.4 at `n = 1e5` and 251 -> 222 at `1e8`.

**Still open from #6 (O13):** the batch size is a free parameter, not an
optimised one. It is a genuine space-time trade -- at `m = 256`, `b = 1` costs
0 ancillas and 2.4e7 T gates while `b = m` costs 255 ancillas and 1.8e6 -- and
the model simply takes the full batch. The optimum depends on `d`, which depends
on the T-count, so choosing it needs an inner loop. `crossovers.md` shows what
the knob is worth; nothing chooses it.

### What #5 closed, and what it did not

Pinnacle was costed more loosely than the surface code in three ways, and the
paper supplied every missing piece.

* **Workspace.** It took the Hamming-weight T-count saving without paying the
  workspace the surface code pays. Both arms compile the same circuit; they may
  not keep different books on it.
* **Magic fidelity.** Nothing certified its T states. At `n = 1e8` the union
  bound was 9.5x over allowance -- the review's arithmetic, reproduced. The
  engine is now selected from the paper's four published specifications
  (Eqs. 7-10, 592 to 5430 qubits, reject rates 0.2% to 10%), and `n >= 1e7`
  forces the 5430-qubit 1e-11 engine.
* **Schedule.** "One state per logical cycle" was an assumption. Their Eq. (11)
  gives `t_me = max(2 d_a + 4r, t_r + 4r, d_a + t_r + 3r)` with `t_r = 10`,
  reproducing 14 / 18 / 23 / 26, and says it lower-bounds the processing unit's
  logical cycle. At `p = 1e-3` that is 23 cycles against the `d = 16` code's 18,
  so the processor STALLS -- live at `n = 1e5`. A T state now costs
  `max(dt, t_me) / (1 - p_reject)` code cycles.

Above `p = 1e-3` the arm refuses, as the surface code does. It previously
reached `m = 13.9` at `3e-3` purely because nothing was looking.

Effect: 20.3 -> 11.5 at `1e5`, 48.3 -> 38.6 at `1e6`, 144 -> 126 at `1e7`,
383 -> 350 at `1e8`. Pinnacle remains the strongest arm; it is now the strongest
arm on the same books. `crossovers.md` carries the field-by-field comparison.

**Still open from #5:** `pin_engines > 1` is charged its footprint but has no
dependency schedule -- the paper equips each processing unit with exactly one
engine, and more than one is our extrapolation. The GB code family also still
tops out at d = 24 (O3), which now matters more: the stall condition pushes
every `p = 1e-3` point onto that code.

### What #4 closed, and what it did not

The magic plant is now **thirteen published operating points from Litinski's
Table 1 plus two cultivation points**, not a formula. The cubic input law
`35 p_in^3` is gone: it described the suppression of input-state error and
ignored faults in the distillation circuitry, giving a two-level factory 3.2e-21
at 42.6 cycles where Litinski's simulation says 4.5e-20 at 128.

Selection now consults `cfg.p` (it did not — the same specification came out at
1e-5, 1e-3 and 3e-3) and minimises **units x footprint** inside the distance
loop, rather than picking the smallest unit and counting afterwards. The magic
budget carries the same factor-of-two failure-to-bias conversion as the logical
check. Rejection was verified to be already inside both sources' published time
costs: Litinski's cycles are `6 d_m / (1 - p_fail)`, cultivation's volume is
quoted including retries at a 99% discard rate.

Consequence at `n = 1e6`: surface FT goes 48.7 -> 227 at `p = 1e-5`, 48.7 -> 108
at `1e-4`, is unchanged at `1e-3`, and **refuses** above it, because neither
source is characterised there.

**Still open from #4:** cultivation's expected volume (~3e4 qubit-rounds at
`p = 1e-3`) is read off a log-log scatter plot in arXiv:2409.17595 Fig. 1 and is
good to about a factor of two. It is the weakest number in `ftqc.py`, and it
matters, because cultivation wins the plant comparison at every point below
`n ~ 1e7`. A digitised value, or the authors' tabulated volume, would settle it.

### What #3 closed, and what it did not

`fhcost/` is now the only model. `make_record.py` writes `RESULTS.json` and
stamps every consumer from it: the explorer's entire constants block, the
headline prose in `HEADLINES.md` and `README.md`, and a model fingerprint on the
figure, the slide, `crossovers.md` and the explorer header. `make_record.py
--check` fails if any of them is stale, and `selftest` asserts the fingerprints
agree.

The JavaScript is still a second implementation of the *functions*. It is no
longer a second set of *assumptions*, and it is no longer trusted: the page
recomputes eight configurations x six lattice sizes x sixteen intermediates,
plus eight qubit counts x seven architectures, plus the classical band, against
the Python answers on load, and `check_parity.py` reads the verdict out of a
headless Chrome. A 1e-7 nudge to one line of the port trips 209 probes.

**Still open from #3:** the port is verified, not generated. A future change to
`fhcost/` will make the parity check go red rather than fixing the JavaScript
automatically. Generating the port from a shared intermediate representation
would remove the remaining hand-maintenance; it is not worth it yet.

**Also fixed in passing:** one headline line asserted that nothing classical
converges, on the strength of the finite-size buffer at t = 0 (2 sites) -- which
actually shows the opposite. The headline now evaluates the criterion
self-consistently at each arm's own t_max, where nothing on the figure closes the
loop, exact diagonalisation included.
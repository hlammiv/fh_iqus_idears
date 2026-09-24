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

**Item 1 is now checked, and the assumed ratio was wrong.** NISQ pays the
one-norm over *all* two-qubit gates in the cone (`G_cone`); STAR pays it only
over *rotations* (`n_rot_cone`), its Cliffords being error-corrected. That is a
real architectural difference, and the margin is then set almost entirely by the
ratio `c_rot/c_g`, which was `5/15 = 0.333` by estimate.

The experimental paper gives its compiled swap-network gate count in **closed
form** (App. C 5, Eqs. C20–C27) and tabulates it. `fhcost/compiled.py`
implements it and reproduces their Table I exactly for the two `(Lx even, Ly
odd)` rows including the `4 × 7` they ran, and reproduces their quoted total:
`(588 × 4) + 24 + 39 = 2415`. It refuses on the other two rows, which come out
6 and 5 low — close, but not reproduced, so not used.

| quantity | measured on their circuit | model | verdict |
|---|---|---|---|
| two-qubit gates/site/step, **excluding** FSWAPs | 13.0 | `c_g = 15` | **+15%, conservative** |
| … including FSWAPs | 21.0 | — | not comparable: routing is charged separately through `routing_power` |
| arbitrary-angle rotations/site/step | **9.0** | `c_rot` was 5 | **−44%, optimistic** |
| `c_rot/c_g` | **0.69** | 0.333 | **2× optimistic** |

`c_g` is fine — better than fine, once the swap network (38% of their gates) is
separated out, since this model charges routing on its own axis and folding their
FSWAPs into `c_g` would double-count. **`c_rot` was not, and it was not even
close.**

The rotation count is not a compilation detail, it is **term counting, and it is
encoding-independent**. On a doubly periodic square lattice each site owns 2
bonds, so there are `2m` bonds × 2 spins = `4m` hopping terms; a second-order
step applies each *twice*; plus `m` on-site terms:

> rotations per step = `8m + m` = **`9m`, exactly**

and that is what their circuit gives — `2(n_h + n_v + n_boundary) + Lx Ly = 252`
on 28 sites. No compilation can use fewer arbitrary-angle rotations than there
are term exponentials, so **`c_rot = 9` is a floor, not a fit, and the `5` it
replaced was below the number of terms in the Hamiltonian.** Their quoted 4627
one-qubit gates put a loose ceiling of 41 per site per step on it, so 9 is the
bottom of the range, not the top.

**`c_rot = 9` is now the default**, and it moves things:

| | n = 10⁶ | n = 10⁸ |
|---|---|---|
| STAR | 20.8 → **15.7** (−24%) | 30.9 → **22.9** (−26%) |
| surface FT | 9.5 → **7.9** (−17%) | 112 → **73.3** (−35%) |
| NISQ+PEC | unchanged | unchanged |

NISQ is untouched because it pays over gates, not rotations — which is exactly
the asymmetry item 1 was about. The surface arm takes the largest hit because its
T-count is rotation-synthesis dominated, and **the surface/STAR crossover moves
out by a decade**, from below `n = 10⁷` to between `10⁷` and `10⁸`.

Item 2 (does the STAR Clifford layer really contribute zero sampling overhead)
is unchanged and still open.
2. Confirm the STAR Clifford layer really contributes **zero** sampling overhead
   and only residual logical error (which is separately budgeted through `p_L`
   and the code distance).

### O2. Signal parameters — **fitted**, with residual caveats
Fitted to Zenodo 17799843 (dimer-link C^zz, TFLO+GPR, both U). The decay is
**Gaussian**, not exponential (free stretch exponent 2.15 +- 0.09 and
1.98 +- 0.05). `s_short = 1.000 +- 0.004` confirms the triplet algebra.
The placeholder `s_res_slope = 0.10` was **5x too large** (fitted 0.021).

**Caveats 2-4 are now MEASURED against the exact answer at U = 0** — possible
only since `free_fermion.py` started reproducing the deposit, and run by
`calibration/tdvp_check.py`:

| t | 1.0 | 1.2 | 1.4 | **1.6** | 1.8 | 2.0 |
|---|---|---|---|---|---|---|
| exact mean \|C\| | 0.0753 | 0.0630 | 0.0421 | **0.0309** | 0.0343 | 0.0411 |
| TDVP χ = 2048 | 0.1005 | 0.1102 | 0.1229 | 0.1379 | 0.1582 | 0.1782 |

Remaining caveats:
1. **Only two U values.** `s_res_slope` and `t_melt_slope` are a two-point
   interpolation with zero constraint on curvature. No basis for extrapolating
   beyond U/J = 4; the linear form is an assumption. *(Unchanged — the deposit
   has no third U.)*
2. **`s_res` is not a real plateau — confirmed exactly.** The exact residual
   averages **0.0467** over `t ∈ [1, 2]` and swings **±47%** about it, decaying to
   a minimum of 0.0309 at `t = 1.6` and reviving to 0.0411. A single scalar is a
   coarse summary; the honest form is Gaussian decay plus a damped revival. The
   model's fitted `s_res_min = 0.043` sits **8% below** the exact mean, which is
   the conservative side: a smaller signal buys more shots.
3. **No per-link exact reference at U=4** — *verified, not assumed*: the `U_4`
   deposit's `exact` group carries only Majorana propagation, and only for `Z`,
   `doublon_sum` and `triplet_density`. There is no
   `spin_correlator_neighbours` row to check `s_res(U=4)` against, so it still
   rests on mitigated data validated globally (n_triplets vs Majorana
   propagation, ~3%).
4. **Do not use TDVP as ground truth for the residual**, and the reason is worse
   than an error bar. At χ = 2048 it averages **0.134, 2.9× the exact residual**,
   and it rises **monotonically** across the window where the exact answer decays
   to a minimum and revives. It does not merely have an error there — it
   manufactures signal. Three selftest checks hold this.

### O3. The Pinnacle arm has a hard ceiling — and it is the MAGIC, not the code
*This item said the wrong thing. Asking the model which branch refuses corrects it.*

The claim was: only five published generalised bicycle codes, topping out at
`d = 24`, so the curve saturates artificially past `n ~ 1e10`, and the plotted
range is comfortably inside the valid region. Two of those three are wrong.

**The ceiling is `m ≈ 391`, at ANY budget, and it is the magic engine.** Above it
`pinnacle_point` returns `None` — the arm does not flatten, it stops existing.
`ftqc.pinnacle_ceiling_cause` asks which branch refuses: just above the ceiling
the required T-state error is `9.7e-12`, and the cleanest engine
`PIN_ENGINE_TABLE` publishes is `1e-11`. `select_engine` returns `None` **before
any code is tried**. Raising `pin_engines` from 1 to 16 does not move the ceiling
by one part in 10⁴, which is the same statement from the other side: the problem
is T-state *quality*, not throughput. The surface arm has no such wall, because
the Litinski ladder cascades to far lower `p_out`.

**The code family is fully consumed from `n ~ 1e6`, not `1e10`.** The largest
code is the one selected at every plotted point from a million qubits up
(`family_exhausted` is now a field on every `pinnacle_point`). So the reassurance
that the plotted range sits inside the valid region was false — it sits on the
last row of the table. The *values* there are still sound, because `d = 24` does
meet the error budget; what does not exist is headroom. If anything tightened —
a longer evolution, a smaller tolerance, a worse `p` — there is no next code.

So O3 splits in two. The part that limits the curve is the engine table, which is
the same kind of gap as the factory ladder and would be closed by the architects
publishing a cleaner engine. The part about the code family is real but
non-binding, and is now a flag rather than a docstring.

### O4. The ZNE response model is unvalidated — and their data CANNOT validate it
Every ZNE conclusion rests on `s(lambda) = s(0) exp(-lambda Lambda)`. A response
with less curvature raises every bias ceiling.

The obvious external test is the published hardware data: twenty evolution times
on a real 56-qubit machine against an exact reference. `calibration/
noise_response.py` runs it, and the answer is that **the dataset is structurally
incapable of testing the response shape**, for a reason worth writing down.

**Their circuit does not get deeper with time.** Sec. C: *"We execute k = 4
second-order Trotter steps ... and time-evolve the initial state up to time
t = 2."* A fixed **count**, not a step density — the step size grows with `t`
and the gate count does not. So the twenty reported times are twenty *repeats of
one noise strength*, not a scan over it, and no amount of statistics constrains
the curvature of `s(lambda)` from them.

What the data do give is that one point, well measured. On the dimer links, in
the six times where the exact signal is above 0.1 and the ratio is
well-conditioned:

- best-conditioned point, `t = 0.1`: attenuation `0.8096 ± 0.0039`, **`Lambda = 0.211`**
- pooled over the six: `Lambda = 0.15`, attenuation spread `±0.12`

The spread is **7x** the quoted shot errors and the fitted trend in `t` is
*negative* — less damping at longer times, which no noise process does. That
residual is Trotter error: their circuit carries it, and no noiseless Trotterised
series was deposited to divide it out. The `t = 0.1` point is where the step size
is smallest and therefore where that contamination is smallest, which is why it
is the one to anchor on — and it lands on the `0.20` the model already uses.

**What would settle it** is unchanged in kind but now specific: noise-amplified
runs of the *same* circuit at several gains, which is what a ZNE experiment
produces and what this deposit does not contain.

### O5. Which gates damp the observable — measured, but the scaling is not
*Raised while working #4. Partially supersedes review #9.*

Measured on the Phasecraft/Quantinuum circuit (2415 two-qubit gates, their quoted
`p = 1e-3`): observed attenuation `Lambda = 0.20`, against **1.76** from the cone
model — an **8.8x overcharge**. (It read 2.58 and 12.9x until the decomposition
was generated rather than typed: that figure left the cone fraction out, i.e.
charged every gate in the circuit, which is not what the cone model does.)
The mechanism was checked on their own data by a
weight test: `Lambda(ZZ)/Lambda(Z) = 1.91`, where damping proportional to operator
weight predicts 2 and uniform-per-gate predicts 1. So a depolarizing error damps
the observable only where the Heisenberg-evolved operator has support.

`damping_model = "support"` implements this and reproduces the measurement to 8%
(0.184 vs 0.200). **It is not the default**, because the single measurement
constrained the *value* at one operating point but not the *scaling*.

**The scaling is now measured, and the obvious repair is excluded.** At `U = 0`
"the operator must eventually fill the lattice" is a calculation, not a worry,
and `calibration/support_growth.py` performs it on the same 7x4 instance the
attenuation was measured on -- the effective support being the inverse
participation ratio of the Heisenberg-evolved operator's single-particle weight,
which is exactly 4 at `t = 0` and so is the same quantity as `w_obs0`:

| t | 0.1 | 0.25 | 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|---|
| effective support, of 56 modes | 4.25 | 5.76 | 14.3 | 31.3 | 43.8 | 21.0 |

It reaches half the register by `t = 0.7` and averages **25.7 of 56** over their
window. So the support does fill the lattice -- and putting that into the damping
model makes it *worse*, not better:

| damping fraction from | fraction | predicted Lambda | vs measured 0.200 |
|---|---|---|---|
| cone, time-averaged | 0.683 | 1.76 | **8.8x** |
| measured Heisenberg support | 0.459 | 1.18 | **5.9x** |
| bare observable weight, `w = 4` | 0.071 | 0.184 | **0.9x** |

Damping does not track the support even though the support spreads. Their weight
test says the same from the other side: `Lambda(ZZ)/Lambda(Z) = 1.91` tracks the
*bare* weights 2 and 1, which an operator spread over 26 modes could not do.

So `support_growth = 0` is no longer an unconstrained fit -- it is the only one
of the three candidates that survives their measurement, and the stated reason
for distrusting it has been tested and does not hold. **The default still stays
`cone`**, so no headline moves; what this removes is the hole, not the caveat.

What is still unmeasured is the same statement at `U > 0`, where the operator
spreads faster and no polynomial reference exists. The `U = 0` support is a lower
bound on the spreading, and the measurement above says even the lower bound
already over-predicts the damping, so a faster-spreading operator does not rescue
the support law.

Effect where it does apply: the discount goes as `w_obs/q ~ 1/m`, so it is large
for a small observable in a big register (their case, 13x) and modest at the m
this model actually reaches (1.21x at `n = 1e6`).

**Relation to review #9.** The reviewer showed the cone fraction should be 2/3,
not 1/3 — a 2x correction *unfavourable* to the results. The measurement says the
cone framing is the wrong quantity and overcharges by ~13x, *favourably*. Both
cannot stand. Only the question of which gates damp is addressed here, and by
measurement rather than by geometry.

**The velocity taxonomy is now measured too, and it survives.** The model carries
three velocities for three jobs — `v = 2` for the light cone and `t_max`,
`v_corr = 4` for correlation spreading, `v_lr = 20` for the Lieb-Robinson
constant — and #9 called that inconsistent. At `U = 0` it is measurable:
`calibration/support_growth.py` puts the observable at the centre of a 29 × 29
open lattice and tracks the Chebyshev radius of the Heisenberg weight.

| what is tracked | measured velocity |
|---|---|
| bulk (mean radius) | **1.61** |
| front at weight 10⁻² | **2.12** |
| front at 10⁻⁴ | 2.35 |
| front at 10⁻⁸ | 3.20 |
| front at 10⁻¹² | **3.40** |

The velocity is **threshold-dependent**, which is precisely why one number cannot
serve: the exponentially small tail outruns the bulk. The `10⁻²` front — the
physical cone — moves at 2.12 against the exact max axial group velocity of
`2J = 2` for `-2J(cos kx + cos ky)`, and against the model's `v = 2`. The tail
keeps accelerating as the threshold drops, which is exactly what a
Lieb-Robinson constant has to cover. So the three numbers are three different
questions and they bracket the measurement in the right order,
`1.61 < 2.12 < 3.40 < 4 ≤ 20`. A selftest asserts that ordering.

**And the cluster buffer with it.** `classical.cluster_radius` charges
`v t + xi ln(1/eps)` so the tail outside the cone is below `eps`, with `xi = 1` a
Config default that nothing stood behind. The same scan measures it as
`front(eps) - v t`:

| eps | mean buffer (sites) | ln(1/eps) | implied xi |
|---|---|---|---|
| 10⁻² | 1.86 | 4.61 | **0.40** |
| 10⁻⁴ | 3.50 | 9.21 | 0.38 |
| 10⁻⁸ | 6.00 | 18.42 | 0.33 |
| 10⁻¹² | 8.25 | 27.63 | **0.30** |

The buffer is flat in `t` at fixed `eps`, so the `v t + xi ln(1/eps)` *shape* is
right. The implied `xi` is **0.30–0.40** against the default `1`, so the cluster
radius is charged **2.5–3.3x** more buffer than the measured tail needs —
conservative, in the direction that costs the classical arm. And `xi` drifts
*down* as `eps` falls, because a free-fermion tail is super-exponential rather
than exponential: a single `xi` is a bound, not a fit, and the model should keep
treating it as one.

Measured at `U = 0`; neither the butterfly velocity nor the tail shape at `U = 4`
is constrained by this.

**The integration geometry is the last piece, and it is now actually checked.**
`hubbard.mean_cone_fraction`'s docstring said its closed form
`2/3 + 1/(2√m) − 1/(6 m^{3/2})` "is checked against numerical integration in
selftest". It said that before the check existed. It does now, and the closed
form reproduces the integral to **1.8e-8** over `m = 4…4096`, tending to `2/3` at
large `m` and rising to `0.896` at `m = 4`. Review #9 was right that `1/3` — the
average of an *uncapped* cone — was the wrong constant; what makes it `2/3` is
the cone hitting the lattice and stopping.

**So #9 is closed.** Which gates damp: measured, and the support law that would
have repaired the cone model excluded (O5). Velocity taxonomy: measured, three
questions, correctly ordered. Cluster buffer: measured, `ξ = 1` conservative by
2.5–3.3×. Integration geometry: the promised check exists. The one thing not
settled is that all of it is at `U = 0`.

### O6. The demonstrated-mitigation arm reaches nothing at our step count
`strategy = "expcal"` costs TFLO+GPR at its measured effective overhead (0.08 —
*below* one, because GPR borrows statistics across correlated time points) and is
drawn only out to the largest `Lambda` actually demonstrated, 0.20. Under this
model's workload it reaches **nothing**: `Lambda(m=4) = 2.87`. Adopt the
experiment's fixed step density `r = 4 tau` and it immediately reaches `m = 20`
at `n = 1e6`.

So the whole distance between this model and a real 56-qubit experiment is the
**Trotter step count**, not the mitigation scheme, not the channel conventions,
and not the damping geometry. That makes review #5 — which asks for an actual
error bound behind the step count — the highest-value item remaining.

**And the experiment's step rule is worse than "r = 4 tau" — it is `k = 4`, full
stop.** The paper is explicit (Sec. C): four second-order steps for the whole
evolution to `t = 2`, so the step size grows with the time and the depth does
not. `trotter_mode = "fixed_count"` implements what they actually ran. Under it
the demonstrated arm does not reach `m = 20`, it reaches `m = 3.3e5` at
`n = 1e6` — which is the point. Both fixed conventions buy reach by abandoning
accuracy, and the amount is now computed rather than asserted:

| step rule | m reached at n = 1e6 | implied Trotter error vs its own budget |
|---|---|---|
| model's own, error-controlled | 0 | within budget by construction |
| slide's `r = 4 tau` | 20 | **25x** over |
| experiment's `k = 4` | 3.3e5 | **2.2e8x** over |

A fixed step count is a budget for one lattice at one time. It is not a rule, and
extrapolating it is what produces the absurd reach — so the `k = 4` row is
evidence about the convention, not a capability claim.

### O7. The ideal p=0 curve is capped by its convergence ceiling — **DONE**

Once the lattice exceeds `m_certified(t)` the finite answer already IS the
infinite answer, so extra sites buy nothing; past that point qubits should buy
**time**, not lattice. `nisq.converged_cap` is the largest `m` with
`m <= m_certified(t_max(m))`, and `max_m_ideal` now takes the minimum of that and
the qubit/clock limit. Ported to the explorer (which needed `lr_error`,
`m_certified` and `fs_speed` alongside) and verified in the browser.

**In fixed-t mode it bites.** At `t = 1` the ceiling is **484 sites** while the
uncapped line reached 333,333 at `n = 10⁶` — a **689×** overshoot. The line now
flattens there at every larger `n`, which is the correct shape: the answer stops
improving.

*(The numbers in the previous version of this item — 244 sites, 1370× — were
typed before the finite-size criterion was replaced with the factorial
Lieb-Robinson form, and never followed it. `check_docs.py` now holds these.)*

**Under the default convention it is inert, and that is the finding.** With
`t_max = sqrt(m)/v` the ceiling recedes with the lattice:

| m | 4 | 256 | 65,536 | 4.2 × 10⁶ |
|---|---|---|---|---|
| `m_certified(t_max(m)) / m` | 121× | 40× | 30× | 30× |

It falls towards an asymptote near **30** and never reaches 1. So the cap can
never bind there — not because the calculations are converged but because **none
of them is**, at any size. Within that convention the slope-1 line is technically
correct, which exposes the worse problem: the y axis is not measuring a physics
answer, only the largest finite lattice that fits. A selftest asserts the ratio
stays above 20, so the claim cannot quietly become false.

This was the third distinct place the time-window convention drove a wrong
conclusion, after finite-size extrapolation buying nothing and the classical band
being read as a hard wall. Still not fixed: the axis label in `sqrt(m)` mode
should say what it measures.

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

* The **m-extrapolation is supported by locality only at SHORT time**, and the
  earlier version of this item quoted the best case as if it were general. The
  spread of `W_eff` across n >= 6, measured (`calibration/domain_check.py`):

      tau = 0.25   1.01x        tau = 1.0   1.26x
      tau = 0.5    1.27x        tau = 2.0   2.91x

  Flat to 1% at tau = 0.25 is what a two-site observable must do once the lattice
  exceeds its cone. By tau = 2 it is a factor of three, and the adopted
  convention t = sqrt(m)/v puts EVERY plotted point at tau >= 1 with most far
  beyond 2. The locality argument therefore does not license the extrapolation
  where the model actually uses it.
* The **tau-extrapolation is only a clamp**, and interpolation INSIDE the range
  is already unreliable. Holding tau = 1 out and predicting it from tau = 0.5 and
  2 overshoots by **+124% at n = 9 and +86% at n = 12**. `W_eff` falls ~15x over
  the measured range and is held fixed beyond tau = 2; a tau^3 law stretched 6.6x
  past data that cannot interpolate to better than a factor of two is not
  defensible on its own.
* The **r^-2 power law is clean where the second-order arm operates** (slopes
  -2.00 to -2.11 for r >= 8, and the arm runs at r = 11-400) but **not where the
  multiproduct arm does**: MPF-4 operates at r = 3.2 upward and below r = 8 for
  m <= 13, and at small r with long tau the measured slope ranges from -0.58 to
  +0.63 because the error is saturated rather than asymptotic.
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

**Review test 3 is now MEASURED, and it answers half the question.** `W_eff` at
the pairs `(m, tau = sqrt(m)/v)` -- the trajectory the model actually walks, not
a fixed-time sweep -- for every patch up to `n = 14`, each with an exact
reference converged to `~1e-15`:

| n | tau | W measured | W from the table | ratio |
|---|---|---|---|---|
| 4 | 1.000 | 0.25608 | 0.18990 | 1.35 |
| 5 | 1.118 | 0.19386 | 0.17614 | 1.10 |
| 6 | 1.225 | 0.07099 | 0.16564 | **0.43** |
| 9 | 1.500 | 0.10231 | 0.14447 | 0.71 |
| 12 | 1.732 | 0.12250 | 0.13112 | 0.93 |
| 14 | 1.871 | 0.08103 | 0.12448 | 0.65 |
| **16** | **2.000** | **0.08026** | **0.11900** | **0.67** |

**The good half.** For every `n >= 6` the fixed-time table **over**-predicts `W`
at the operating point, by 1.1-2.3x. That is the conservative direction: a larger
`W` buys more Trotter steps and so more gates than the trajectory needs. The
central assumption under test -- that `W` depends on `tau` alone, so the table can
be read at `tau = sqrt(m)/v` -- holds to about a factor of two, in our favour.

**The half that is not settled, and cannot be at this scale.** Fitting `W` against
`m` gives `W ~ m^+0.31` and an implied `alpha = 1.75 + gamma/2 = 1.91`, against
`1.58` from the table -- a 7x difference in gate count by `m = 400`, in the
*anti*-conservative direction. But that is four points with real scatter:

> `alpha = 1.81 +- 0.16`, 95% CI **1.50 to 2.13** -- and including `n = 4, 5`
> flips it to **1.43**.

The adopted `1.75` sits inside. These lattices do not determine the exponent, and
saying they give `1.81` would be reading noise. `domain_check.py` prints the
interval next to the point estimate so the number cannot be quoted without it.

**`n = 16` landed, and it did what was predicted in advance and no more.** Before
it, the fit was `1.91 +- 0.20` over four points with a CI of `1.50-2.31`; the
forecast was that a fifth point would narrow the interval by about a third and
not close it. It narrowed the top from 2.31 to 2.13 and moved the estimate toward
the adopted value. One sub-claim did change: Campbell's `2.25` is now *outside*
the interval, which is what a loose upper bound should look like rather than a
contradiction — worst-case Trotter bounds are known to overestimate, which is why
`f_emp` exists at all. `domain_check.py --drop16` re-runs the fit without the
point and reaches the same headline.

**The inclusion rule changed to admit it, and that deserves saying plainly.** The
old rule was `||psi(sub) - psi(2 sub)|| < 1e-10`, with the stated purpose "the
reference must be far better than the signal" — a *ratio* statement behind an
absolute proxy that held only because substepping to machine precision is cheap
up to `n = 14`. At `n = 16` it is not: 165M amplitudes on a 12-vector Krylov
basis, 1.5-3 h per halving. The proxy is also the wrong quantity — `W_eff` comes
from `C^zz`, and at `n = 16` the two differ by three orders of magnitude
(`||dψ|| = 2.8e-7` against `|dC^zz| = 5.9e-10`).

So the rule is now the ratio it stood for: the reference's error must be a
negligible fraction of the smallest Trotter error it measures. It is applied
uniformly and is **not tuned to admit `n = 16`** — the six smaller patches pass by
ten orders of magnitude and `n = 16` by six:

| n | 4 | 5 | 6 | 9 | 12 | 14 | **16** |
|---|---|---|---|---|---|---|---|
| reference err / Trotter err | 8e-12 | 1e-11 | 6e-11 | 3e-11 | 2e-11 | 4e-11 | **9e-07** |

(The point's *as-run* convergence of 2.8e-6 was itself an artefact: it ran before
the inverted probe was fixed, so that number is the coarse comparison's own
error. Re-measured with the corrected probe it is 2.8e-7 on the state and 5.9e-10
on the observable. `--verify-ref` writes a sidecar and `domain_check` applies it
as an explicit override, keeping both values.)

The O11b multiproduct rerun is the remaining outstanding compute, plus the
`n = 16` clamp probe now running.

### O10. The tensor-network question is OPEN, not resolved
*Raised by #8. This is the item most likely to reverse a headline.*

The classical band is now **ED-only, m = 24-26**, down from 24-62. That narrowing
used to flip "mitigated NISQ never clears classical" into "it clears at
n ~ 1.7e4". It no longer does: since `c_rot` was measured, bare NISQ+PEC clears
neither band at any `n` up to 10¹³. What the narrowing flips now is the
MULTIPRODUCT arm — `n = 4.0e8` against the narrow band, never against the wide
one. **The band narrowed because an unvalidated arm was removed, not
because tensor networks were shown to fail** -- and removing an arm is a change in
the direction that flatters the quantum curves, which deserves suspicion.

**Sharpened by a diagnostic that cost nothing to run, then corrected by doing it
properly.** The first pass compared two bond dimensions at three times. Reading
*all four* chi at *every* time -- `calibration/tdvp_check.py`, against the full
Zenodo deposit rather than the summary that was stored -- shows the chi
dependence is **non-monotonic in time**, which the two-point comparison hid:

| t | mean \|C\| | chi=256 | chi=512 | chi=1024 | chi=2048 | slope | reading |
|---|---|---|---|---|---|---|---|
| 0.1 | 0.9418 | 7.94e-3 | 7.38e-3 | 7.37e-3 | 7.35e-3 | **-0.03** | floored |
| 0.3 | 0.5869 | 2.98e-2 | 2.12e-2 | 1.78e-2 | 1.60e-2 | -0.29 | truncation-limited |
| 0.5 | 0.2423 | 6.43e-2 | 4.84e-2 | 3.37e-2 | 2.62e-2 | **-0.44** | truncation-limited |
| 0.7 | 0.0931 | 6.93e-2 | 5.76e-2 | 4.15e-2 | 3.28e-2 | -0.37 | truncation-limited |
| 1.0 | 0.0753 | 5.03e-2 | 4.19e-2 | 3.36e-2 | 2.91e-2 | -0.27 | truncation-limited |
| 1.5 | 0.0347 | 1.19e-1 | 1.05e-1 | 1.04e-1 | 1.01e-1 | -0.07 | stalling |
| 2.0 | 0.0411 | 1.63e-1 | 1.38e-1 | 1.38e-1 | 1.51e-1 | **-0.04** | floored |

(slope = d log err / d log chi; dimer links, U = 0, against the deposit's own
`Exact` rows -- NOT its `FLO` rows, which differ from `Exact` by up to 2.9e-2 at
t = 2; `free_fermion.py --deposit` shows `Exact` is the exact one. The script
also prints the all-756-pairs version, which tells the same story.)

The numbers previously stored in `classical.py` **reproduce exactly** from the
full release, so the extraction was right. The *interpretation* was not. "A
chi-independent floor from the earliest times" holds at t = 0.1 and again from
t >= 1.5, and is **false through the middle of the window**, where the run is
genuinely bond-dimension-limited -- 8x the chi buys 2.5x the error at t = 0.5.
The pooled `TDVP_SLOPE = -0.12` averages a floor, a converging regime and a
second floor, and cannot be extrapolated in either direction.

Both ends still defeat the dataset as a bound, for different reasons:

- **t = 0.1.** The state is barely entangled, chi = 256 is wild overkill, and a
  truncation-limited calculation would be at machine precision. Instead the error
  is 7.9e-3 and an **eightfold** increase in bond dimension removes **8%** of it.
  The deposit's `max_bond_dimension` column confirms every run saturated its cap,
  so truncation was binding -- it simply was not what limited the accuracy.
- **t >= 1.5.** The chi = 2048 error (0.13) is about **four times** the mean
  |C^zz| there (0.034), six times the signed mean (0.021-0.031). A calculation
  whose error exceeds its answer bounds nothing.

So the published TDVP is **not a converged tensor-network calculation**, the
chi-extrapolation to 2e12 is meaningless, and this dataset **cannot bound what
tensor networks can do here**. What it is *not* is uniformly chi-insensitive, and
the repo no longer says that: `TDVP_SLOPE_BY_T` carries the time resolution and
five selftest checks tie it to the regenerated JSON.

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

**The dt axis has now been run, and it answers most of the question.**
`calibration/tdvp_dt.py` varies the one axis the deposit never varies, on the
same instance, in the same library. Three results, in order of how firmly they
hold:

**1. The integrator is not the problem — settled.** At 3 × 4 with χ = 128, where
the bond dimension is ample, the error falls **cleanly to zero** with the time
step: `2.79e-2 → 2.25e-3` over 16× in dt, with the halving ratio climbing to 2.
Two-site TDVP in TeNPy has no intrinsic dt floor on this problem, so no floor of
the published size can be blamed on the integrator.

**2. At 28 sites, at fixed χ, what is left is bond dimension — not dt.** Fitting
`err = a + b·dt^p` at χ = 256 gives a dt-converged residual of **`a = 1.11e-2`**
(p ≈ 0.9), and holding dt at its finest and raising χ *does* move it:

| χ | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|
| err at dt = 0.00625 | 1.81e-2 | 1.36e-2 | 1.17e-2 | **9.65e-3** |
| discarded weight | 6.6e-4 | 7.0e-5 | 1.5e-4 | 1.1e-4 |

slope **−0.25** from χ = 256 to 1024. So *our* run at t = 0.1 is
truncation-limited, and converging it in dt alone does not remove the error.

**3. Which is exactly why we cannot claim to have reproduced their floor.** Their
run at the same χ = 256 has a *smaller* error (7.94e-3 against our 1.36e-2) and
is already **flat** in χ (slope −0.03) where ours still falls at −0.25. We reach
their χ = 256 number only near χ = 1024. The likely difference is the **MPS
ordering**: ours is row-major on the torus, theirs is the Jordan-Wigner snake
matched to their swap network, and a better ordering means less entanglement
across the MPS cuts. That is a property of their compilation, not of TDVP.

So the narrow question — *is their* `t = 0.1` *floor the time step?* — is not
closed by this, and saying otherwise would be over-reading a run that is not a
faithful replica. What **is** closed is that the integrator cannot be the
explanation, which was the alternative that would have voided the
chi-extrapolation in both directions.

**And a sharper statement arrived from a different direction** (see O2): at
χ = 2048 their run's |C^zz| rises *monotonically* across `t ∈ [1, 2]`, 2.9× the
exact residual, where the exact answer decays to a minimum at `t = 1.6` and
revives. It does not merely carry an error at late times, it produces the wrong
shape. That settles "not a converged tensor-network calculation" without needing
the dt question at all.

*What exists already, checked before writing anything.* The collaboration's own
runs are **TeNPy** -- their metadata literally reads `Tenpy $\chi=512$` -- and
TeNPy ships Fermi-Hubbard and two-site TDVP. So the dt scan is a configuration of
their library at their parameters, not a new fermionic 2D TDVP on a doubly
periodic torus; hand-rolling one would risk a confidently wrong number for no
gain. TeNPy 1.1.0 and quimb 1.15.0 are installed locally; neither is on lenore
yet, which is a one-line install when the scan is ready for production sizes.

Until at least one purpose-built, leadership-scale classical attack is run, the
honest headline is "clears the **exact-diagonalisation** frontier", not "clears
classical".

### O15. The clamp beyond tau = 2 is NOT testable at n <= 14 -- measured
The clamp probe ran on lenore: W_eff beyond tau = 2 on two lattice sizes, with
the decision rule fixed in `calibration/apply_clamp.py` before the data existed.

| tau | n = 12 | n = 14 | clamped | spread | vs clamp |
|---|---|---|---|---|---|
| 2.5 | 0.0773 | 0.0821 | 0.1190 | **1.06x** | 0.69x |
| 3.0 | 0.0674 | 0.0332 | 0.1190 | 2.03x | 0.28x |
| 4.0 | 0.0481 | 0.0319 | 0.1190 | 1.51x | 0.27x |

**Verdict: the sizes DISAGREE and nothing was written.** The n = 12 series alone
was a clean monotone decline that looked like solid support for "the clamp is
conservative"; n = 14 breaks it. That is precisely why the rule was committed
first.

The one exception is **tau = 2.5**, where the two sizes agree to 1.06% and the
causal cone (2 v tau = 10 sites) fits inside both lattices with margin. At
tau = 3 the cone is exactly 12 sites and at tau = 4 it is 16, wider than either
patch -- so the cone argument explains which points are testable, though not
perfectly: tau = 3 nominally fits n = 12 and still disagrees by 2x.

So a tau = 2.5 extension is *available* on one validated point, worth 0.69x in W
and 0.83x in step count. It has NOT been taken, because taking it would mean
relaxing the pre-committed rule after seeing the data. Closing this properly
needs n = 16 or larger, where the cone fits at tau = 3-4.

### O14. Platform arms depend on numbers that are moving fast
The neutral-atom surface arm is **fidelity-limited**: at the measured 5e-3 it
does not exist, and at p <~ 2e-3 it does. That is a falsifiable prediction
against a number that has been improving steadily, and it should be rechecked
against the current best atom two-qubit fidelity rather than Evered et al. 2023.

Helios is **clock-limited** at 55 ms/layer, which is a first-generation figure
for a machine whose paper explicitly discusses clock speed as its main scaling
challenge. Both arms should be re-run when either number moves.

**O14 part two RESOLVED, and it changed the answer for Helios.** The model could
already express the combination -- `max_m_pinnacle` on an all-to-all platform is
high-rate qLDPC *with* transversal rounds -- it just was not being reported,
because the figure computes the Pinnacle line on `sc_long_range` only.

Crediting Helios with the codes it actually proposes -- and the answer is
narrower than first stated:

| n | Helios, surface | Helios, qLDPC |
|---|---|---|
| 1e9 | **9.9** | 6.5 |
| 1e10 | 7.9 | **9.5** |
| 1e12 | 31.5 | **68.6** |

qLDPC does NOT rescue an arm that did not exist -- transversal rounds already
gave Helios a surface arm from `n = 1.5e8`, where its own qLDPC arm only starts
at `6.3e8`. It overtakes only above `n = 4.8e9`, and is worth
~2.2x by `1e12`. Below the crossover the surface code is better,
because Pinnacle's single 4410-qubit magic engine serialises T supply while the
surface plant parallelises: the storage saving (101 physical per logical against
2500 for a d = 25 patch, 25x) does not pay until the lattice is large enough for
storage to dominate.

A first draft of this item asserted "0 at every n" for the surface arm and had to
be corrected against the model; the arm is on the platform figure now with the
crossover visible.

Neutral atoms still get nothing, and for the same reason as their surface arm:
no magic engine is characterised at their measured 5e-3, since the engine table
stops at 1e-3. Fidelity, not code choice and not clock.

### O11. Connectivity was a caveat, not a constraint -- **CLOSED**
`ftqc.pinnacle_point` and `METHODS.md` both pointed at an O11 that was never
written. It said generalised bicycle codes need connectivity the slide's grid
does not provide, and the model costed the arm anyway. `pin_nonlocal` was
supposed to mark this and was read by NOTHING -- setting it changed `model_id`
and no number.

Closed by `fhcost/platform.py`: connectivity now gates code admissibility, so
`max_m_pinnacle` returns 0 on a nearest-neighbour grid. Bravyi et al.
(arXiv:2308.07915) also make the requirement milder than the caveat implied --
degree 6, two edge-disjoint planar subgraphs, i.e. two coupler layers rather than
all-to-all. The arm is plotted on that chip and labelled with it.

**Closed:** `pin_nonlocal` is removed. It was documented as marking an arm
costed under a different hardware assumption and was read by nothing, so setting
it moved `model_id` and no number.

### O11b. Multiproduct coefficients -- **RESOLVED**, and the diagnosis was wrong
*Domain doubled from tau <= 0.5 to tau <= 1.0.*

The earlier note said the coefficient "swings 3-13x at tau >= 1 because the
observable error passes through zero crossings". The finely-sampled rerun on
lenore -- 8 multiproduct bases, 18 step counts, 6 taus -- shows that was wrong on
every count:

* **There are no zero crossings.** The extrapolated error falls smoothly over
  four decades in every case that is not at the precision floor.
* **Order 8 at small tau was measuring the double-precision floor.** The exact
  reference is good to ~1e-15, so a median across bases there averaged numerical
  noise and produced a 9291x "spread".
* **The observed convergence order drifts above nominal at tau >= 1** (4.3-4.6
  for the order-4 formula), so extracting W with a FIXED r^-2k law makes the
  answer drift with base by construction.
* **Lattice size and base were pooled**, so genuine finite-size dependence
  appeared as base-to-base scatter.

Adding bases made the reported spread WORSE (order 4 at tau = 1: 8.7x with two
bases, 15.2x with five), which is what forced looking at the raw error curves
instead of the summary statistic.

`fit_w.mpf_fit` now fits the exponent AND the coefficient, at fixed n, above the
precision floor. Fitted orders: **3.98-4.64, 6.01-6.45, 7.66-9.08** -- nominal,
within the drift subleading terms explain. Order 8 has no tau = 0.25 entry
because every point there is noise, and saying so is the correct output.

Consequences: MPF clears the ED frontier again, at n = 3.8e+08 rather than the
original 4e6; and the optimum is **arm-dependent** -- NISQ is shot-limited so
extra branches pay to order 8, FT is magic-limited so they stop paying after
order 4. The old test asserted a universal interior optimum and was hiding that.

**Still open:** tau > 1 entries exist but carry finite-size spread of 2.9-9.1x
between n = 9 and n = 12, so they are excluded. Widening needs n = 16, which is
the trajectory run's endpoint.

### O11b-old. Multiproduct coefficients are calibrated only to tau <= 0.5
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
| 1 | initial state and connected correlator incompatible | **fixed** — matched to arXiv:2510.26300. Worth recording that the *specification* stayed right while one *implementation* of it did not: `free_fermion.py` used singlet dimers and no π flux until the external check against the deposit caught it (see #7). The spec is not the code |
| 2 | ZNE feasibility omits residual bias | **fixed** — bias-limited, not variance-limited |
| 3 | gatewise PEC overhead uses the wrong coefficient | **fixed** — attenuation separated from one-norm |
| 4 | mitigation theorem is overstated | **fixed** — citation corrected, claims downgraded, demonstrated arm added |
| 5 | multiproduct gains lack an error bound | **partial** — branch-cost accounting is DONE and the status row saying otherwise was stale: each branch is charged as its own circuit at `k_i` × the gates on both the NISQ and the FT side, and charging only `‖c‖₁²` would understate the cost by 405× at order 4 and 654× at order 8. Convergence orders measured (3.98–4.64, 6.01–6.45, 7.66–9.08). What is still missing is a *bound* rather than a validated empirical law — the same status as the second-order calibration |
| 6 | Pinnacle calibration needs reconstruction | **fixed** — footprint reproduced with no free parameter; single-engine T supply corrected; connectivity flagged as O11 |
| 7 | FT resource and error accounting incomplete | **fixed** — HWP workspace charged, magic-state error budgeted with per-point factory selection, failure-to-bias factor 2 |
| 8 | classical band is heuristic, not a ceiling | **partial** — relabelled as capacity, entropy bound fixed, TDVP measured; a validated TN estimate is still missing (O10) |
| 9 | light-cone geometry inconsistent (1/3 vs 2/3) | **fixed** — which-gates-damp measured and the support law that would have repaired it excluded (O5); velocity taxonomy measured (1.61 bulk, 2.12 front, 3.40 tail, correctly bracketed by v = 2, v_corr = 4, v_lr = 20); cluster buffer measured, ξ = 0.30–0.40 against the default 1 so 2.5–3.3× conservative; cone integral now actually checked against numerical integration (1.8e-8). All at U = 0 |
| 10 | error components do not combine to the tolerance | **partial, and the residual is named** — the seven shares sum to exactly 1.00 of the tolerance and `validate` enforces that, but summing is only half an answer: a share allocated and never spent bounds nothing. Six of seven ARE spent at the point of use (Trotter in `trotter_steps`, synthesis in `surface_point`, logical in each arm's rejection test, magic in `select_factory`/`select_engine`, mitigation bias in `time_required`, statistics in the shot count). The seventh, **finite size, is allocated and never spent** — `frac_finite` appears only in `converged.m_certified`, which no plotted arm consults, and under `t_max = √m/v` it *cannot* be spent: certifying the thermodynamic limit needs ~30× the sites at every size (O7). `hubbard.ledger_enforcement` reports this per component and two selftest checks hold it |
| 11 | implementation and reporting issues | **fixed** — all eleven; headlines now generated from one record |

## Second-pass review status

| # | finding | status |
|---|---|---|
| 1 | the convergence criterion fails an exact physical limit | **fixed** — factorial Lieb-Robinson form, exact t -> 0 limit |
| 2 | the default Trotter calibration is used beyond its evidence | **partial** — order-2k coefficients measured (O11b); the trajectory is now measured at every n ≤ **16** and the table is conservative there by 1.1–2.3×, but the exponent is unresolved (α = 1.81 ± 0.16, CI 1.50–2.13) and the domain is still exceeded (O9) |
| 3 | Python, the explorer and the documents disagree | **fixed** — one model, one record; see below |
| 4 | factory cleanup neglects circuit-level errors | **fixed** — published operating points replace the cubic law; selection is p-aware and plant-level |
| 5 | Pinnacle lacks the common FT resource/error ledger | **fixed** — workspace charged, engine selected and certified, stalls and rejection scheduled |
| 6 | Hamming-weight workspace does not match the cited circuit | **fixed** — Campbell Thm 2; workspace, T-count, depth and synthesis all from one construction |
| 7 | the classical baseline misses the U = 0 easy limit | **fixed** — free-fermion estimator implemented, validated internally to 7.8e-16 and against the published exact data to 2.9e-10 (which caught a wrong dimer state and a missing pi flux); band relabelled |
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

**#10 remainder CLOSED, by checking the property that makes the shortcut safe.**
Every admissible integer candidate is now re-costed at its own site count, and
the answer matches `best_lattice(continuous m)` at every point tested -- for the
surface arm and for the block-quantised Pinnacle arm, where k = 14-16 logical
qubits per block makes the footprint lumpy in m.

The reason is worth more than the result: **feasibility is monotone in m**.
Every integer m from 4 to the maximum is feasible at six (arm, n) points, so
`best_lattice`, which takes the largest admissible site count BELOW the
continuous maximum, is feasible by construction. Block quantisation cannot break
that unless it makes a *smaller* lattice infeasible where a larger one is not,
and it does not. The monotonicity is asserted rather than assumed, so if a future
change breaks it the shortcut stops being safe and the test says so.

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

**Closed from #9:** the hatched strip between the axis floor (2.5) and `m = 4`
is only 0.20 of a decade, and up to four arms were hatched into it *on top of
each other* -- legible in principle, unreadable in practice. Lowering the axis
floor would have bought room at the cost of empty space everywhere else. Instead
the strip is now divided into **one lane per arm that runs off**: the vertical
extent down there carries no information beyond "below `m = 4`", so it is spent
on separation. Same information, four times the room, axis unchanged. The note
block says "one lane per arm" so the lanes are not read as values.

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

**#8 remainder CLOSED, and one of its claims was wrong.**

* `s_data_rel_unc` is now **measured, not assumed**: the RMS relative scatter of
  the published points about the fitted envelope, at the small-signal times where
  the bound binds, is **30%** -- the model had assumed 20%. For a LOWER bound on
  the signal, assuming less uncertainty than is observed is the wrong direction,
  so the default moves to 0.30. Effect at n = 1e6: surface FT 12.5 -> 9.5,
  NISQ 13.1 -> 12.4. The measure overstates their error (it contains our
  envelope's model error, and seven points cannot separate the two), which is
  again the conservative direction.
* `s_abs_floor` **was NOT "the single most leveraged number in the ledger"**.
  That claim, written in the #8 writeup, is false: the floor is INERT below
  ~0.038, because the data-based bound is 0.0383 at the binding time and wins the
  max(). Every value in the measured bracket -- 0.013 (envelope scatter) through
  0.020 (the model's) to 0.030 (smallest reported) -- gives identical answers.
  The concern is closed by showing it does not matter, and `s_data_rel_unc` is
  the number that does.

`hubbard.signal_uncertainty_anchors()` returns all four anchors with provenance:
raw shot noise 0.079 (160 shots, bounded +-1 estimator), envelope scatter 0.0126
absolute, 30% relative at small signal, smallest reported value 0.0295.
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
against exact many-body evolution to 7.8e-16** at `m = 4, 6, 8, 9`. An `O(m)`
form reproduces it to 1e-16 and measures a log-log exponent of 0.979 out to
`m = 16384`.

**And then against someone else's answer, which is what caught the physics.**
`--deposit` compares the estimator to the *published* exact `C^zz` on the same
7x4 instance. It did not agree, and chasing that down found two errors that an
internal check structurally cannot see, because it evolves our own state under
our own Hamiltonian:

| | was | is | difference it makes |
|---|---|---|---|
| dimer state | singlet, `(A - B)/sqrt2` | **S^z = 0 triplet**, `(A + B)/sqrt2` | up to **0.20** by t = 0.5 |
| Hamiltonian | no flux | **pi flux through the short cycle**, Peierls `pi/4` per y-bond | up to **0.13** |

Both give `C^zz = -1` at `t = 0`, which is exactly why neither showed up.
With both corrected the estimator reproduces the deposit to **2.9e-10**.

Neither error touches the cost -- same operation count, same memory, complex
arithmetic included -- so the frontier below and the classical band are
unchanged. The claim that moved is the one about fidelity to the experiment.

A third thing fell out: the deposit's own `Exact` and `FLO` rows **disagree**, by
up to 2.9e-2 at t = 2, which is the size of the TDVP errors O10 measures. Our
independent implementation reproduces `Exact` to 3e-10 and `FLO` only to 2.9e-2,
so `Exact` is the reference and `calibration/tdvp_check.py` now uses it.

Measured frontier: `m ~ 4e7` on ONE core of unoptimised CPython in a week,
against 26 for ED. No quantum arm is within six orders of magnitude at that
point. The band is relabelled "estimated ED capacity" everywhere, and the
explorer's U slider now reaches 0 and says so.

The TDVP peak-FLOP arithmetic is labelled as such, with an explicit
`assumed_fraction_of_peak = 1.0`.

**Still open from #7:** the review's third test -- vary the TDVP time step AND
bond dimension against an exact small reference, to identify what causes the
floor -- is not done. The bond-dimension half is now done from the published
deposit (O10): the chi dependence is resolved in time and the floor is isolated
to t = 0.1 and t >= 1.5. The dt half needs a run, but not a hand-rolled one --
the published runs are TeNPy, which has Fermi-Hubbard and two-site TDVP, so it is
a configuration of their library at their parameters. Until it happens the U > 0
band stays an ED capacity with the tensor-network question open.

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

**O13 CLOSED, with a negative result.** The batch is now optimised
(`hwp_batch = -1`), and the optimum is the full batch at all 20 (regime, n)
points tested -- baseline, Willow p_L, one day, one month, eps = 0.01, p = 1e-4,
and each magic family alone. The knob is real but inert here, which is now known
rather than assumed.

The instructive part is the first attempt. Minimising PHYSICAL QUBITS picks
`b = 1` and makes the arm *worse*: reach at `n = 1e8` falls from 125 to 30,
because `b = 1` has no workspace but seven times the T states, so the plant
delivers more slowly and the arm goes clock-bound. The objective has to be
whichever constraint binds, not the footprint. Both facts are asserted.

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
383 -> 350 at `1e8`. Pinnacle remains the strongest arm on the same ERROR ledger
-- but not on the same hardware: once connectivity is enforced (O11) the arm is 0
on the slide's grid, and the figure plots it on a two-coupler-layer chip nobody
has built. `crossovers.md` carries the field-by-field comparison.

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

**#4 remainder CLOSED, from the authors' released stats** (Zenodo
10.5281/zenodo.13777072, mirrored in `calibration/data/cultivation_stats.csv`).
Their `end2end-inplace-distillation` row at `p = 1e-3, d1 = 5, d2 = 15` carries a
complementary-gap histogram -- 117 kept-count bins, 113 error-count bins, 1e12
shots. Reconstructing the error/discard trade from it puts their headline 2e-9 at
gap cut 100: measured error **1.90e-9**, **73.0 attempts** per accepted state,
98.6% discard, against their quoted 99%.

Three of the four inputs are now exact rather than read off a figure:

| | |
|---|---|
| footprint | **463** qubits (was a guessed 450 = 2 x 15^2) |
| rounds per attempt | **20** |
| attempts | **73.0** |
| volume per state | 9.3e3 .. 6.8e5, their integrated ~3e4 inside it |

Only the volume stays uncertain, because discarded attempts terminate early and
the qubit count ramps during cultivation -- which their integration handles and a
flat `q*r*attempts` cannot.

**And it does not matter.** A 73x change in the cycle count moves surface FT at
`n = 1e7` by nothing at all. The plant-level optimisation added in #4 minimises
total plant qubits over ALL admissible sources, so degrading cultivation simply
hands the job to the Litinski ladder. The item called this "the weakest number in
ftqc.py, and it matters"; the first half was true and the second was not.

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

**But the verification itself had a hole, and closing it found a bug.** Parity
ran over fourteen hand-picked variants. A field the JavaScript silently ignores
passes all fourteen if none of them happens to touch it — which is not a
hypothetical: `check_parity.py --fields` now perturbs **every** `Config` field,
keeps the 52 that move a probe, and makes the browser recompute all of them in a
throwaway copy of the page (nothing shipped grows). On its first run it failed on
exactly one field, `simultaneous`, with 34 mismatches.

The JavaScript had `return 2.807033768343811;` commented *"Bonferroni at T = 20,
matching scipy in Python"*. It was neither: 2.807 is the `T = 10` value and
Python gives **3.0233** at `T = 20`, so the explorer understated the simultaneous
confidence factor by 7% and would not have tracked `n_times` in any case. Python
was right — `selftest` asserts `z > 3.0` there — and the port was never
exercised. Fixed by porting Wichura's AS241 inverse normal, which matches `scipy`
to `5e-16` over `T = 1…1000`, so the factor now tracks `n_times` instead of being
frozen at one value of it.

All 52 fields now agree. The sweep also reports the **29 that move no probe**,
which is its own finding rather than reassurance: several of them
(`xi`, `fs_speed`, `v_lr`, `obs_support_sites`) are inert only because the probe
set does not reach `fhcost/converged.py` at all. That is a stated coverage gap
now, not a hidden one.

**Also fixed in passing:** one headline line asserted that nothing classical
converges, on the strength of the finite-size buffer at t = 0 (2 sites) -- which
actually shows the opposite. The headline now evaluates the criterion
self-consistently at each arm's own t_max, where nothing on the figure closes the
loop, exact diagonalisation included.
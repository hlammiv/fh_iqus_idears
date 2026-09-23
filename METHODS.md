# Methods: quantitative resource estimate for 2D Fermi–Hubbard dynamics

Replaces the qualitative whiteboard sketch on slide 2 of `Q1_Fermi_Hubbard-1.pdf`.
Every number in `figures/fh_resource_estimate.pdf` and `crossovers.md` is produced
by `fhcost/`; nothing is hand-placed. Reproduce with

```
python3 -m fhcost.selftest     # ~135 checks, all arithmetic, <2 s, ~27 MB
python3 make_record.py         # the shared result record; stamps every consumer
python3 make_figure.py
python3 crossover_table.py
python3 check_parity.py        # Python vs JavaScript, in a headless browser
```

## 0. The question, and two things slide 1 leaves open

Given a √m × √m Fermi–Hubbard lattice and one week on a √n × √n grid of qubits
(p = 10⁻³, 10 ns gates, 100 ns measurement, plus a classical data centre), how
large an m can ⟨Z_i(t)Z_j(0)⟩ be resolved on, to **relative** error ε = 10⁻¹,
over T time points?

**(a) The time window.** Slide 1 says t ∈ [0, 1/m]. Taken literally this is not a
hard problem: at t = 1/m the causal cone of the correlator is smaller than one
lattice spacing, a cluster expansion is exact at a cost independent of m, and the
classical reach is unbounded at *every* m. Panel (c) of the figure shows exactly
this. We therefore adopt **t_max = √m / v_B** — evolve until the light cone
crosses the lattice — which is the regime where the question has content. Slide 1
should be corrected to match.

**(b) The state and the measurement — resolved by matching the experiment.**
An earlier version of this model assumed a Néel/product start together with the
*connected* correlator. Those are incompatible: a product state in the occupation
basis is a Z eigenstate, `Z_j|psi_0> = z_j|psi_0>`, so

    <Z_i(t) Z_j(0)>_c = <Z_i(t) Z_j(0)> - <Z_i(t)><Z_j(0)> = 0

identically, at every time and for every Hamiltonian. The signal scale was
therefore undefined, and it enters the cost twice.

We now match the state and observable of the Phasecraft/Quantinuum experiment,
[arXiv:2510.26300](https://arxiv.org/abs/2510.26300):

* **Initial state** — a dimerised half-filled configuration: a maximally entangled
  `S^z_tot = 0` TRIPLET on each link of a dimer covering, plus one holon and one
  doublon. The triplet is *not* a Z eigenstate, which is exactly what makes the
  connected correlator nonzero.
* **Observable** — the equal-time connected nearest-neighbour spin correlation
  `C^zz_ij(t) = 4[<S^z_i(t) S^z_j(t)> - <S^z_i(t)><S^z_j(t)>]`, with
  `S^z_i = (n_i,up - n_i,dn)/2`. Read out in the **real-space occupation basis**:
  no Hadamard test, no ancilla. The ancilla is charged only if `observable` is
  switched back to the two-time correlator.
* **Signal** — `|C^zz| = 1` exactly at `t = 0` (triplet algebra), melting over
  `t ~ 0.4-0.7` and leaving a residual antiferromagnetic correlation. The paper
  reports that the melt is **slower** and the residual **larger** at larger `U`.
  Modelled as `s(t) = s_res + (1 - s_res) exp(-t/t_melt)` with `t_melt` and
  `s_res` both rising with `U/J`. The residual magnitude is not quoted
  numerically in the paper, so `s_res_min` and `s_res_slope` are ESTIMATES.
* **Grid** — 7x4 = 28 sites, double-periodic, flux `Phi = pi`, `t in [0.1, 2]` in
  steps of 0.1 (20 points, units of 1/J), `U/J = 0` and 4, second-order Trotter,
  160 shots per point, TFLO + GPR mitigation (not PEC).

The `signal_regime` knob selects **short** (t before the melt, `s ~ 1`), **long**
(the residual), or **curve** (the decay, default). This matters because the two
regimes respond to `U` in *opposite* directions: raising `U` enlarges the
commutator norm, which hurts at short time, but also raises the residual, which
loosens the absolute tolerance and helps at long time. At `n = 10^6`, going from
`U/J = 4` to 8 moves mitigated NISQ **down** 14.3 -> 11.6 in the short regime and
surface FT **up** 48 -> 77 in the long regime.

Two checks worth noting: their lattice, `m = 28`, sits right at the top of the
exact-diagonalisation band estimated in §5 (24-26), an independent confirmation
of that edge; and they report that Trotter error "may be substantially lower than
worst-case commutator bounds would indicate", which supports treating the
bound-derived step count as an upper bound (§8, item 3).

### What the costing model may assume about the signal (review #8)

The tolerance was set from the fitted Gaussian envelope, evaluated once at
`t_max`. Three things were wrong with that, and they compound.

**The envelope is not a lower bound.** The published data is not monotone: at
`U/J = 4` it falls to 0.0466 at `t = 1.5` and comes back to 0.0532 at `t = 2.0`,
where the envelope reads 0.0640. That is 37% high, which is **1.89× too few
shots** at the point that binds. A leave-one-out refit
(`calibration/fit_signal.py`) is worse: the envelope is good to ~1% out to
`t = 0.7` and then misses by **−52% to +55%**.

**One tolerance was doing two jobs.** A Trotter or synthesis error is a *bias*
and must clear the tolerance at every requested time; shots are paid *per time*.
These are now separate:

| | bounded by | used for |
|---|---|---|
| `eps_absolute` | `signal_min` — the tightest time on the grid | Trotter, synthesis, logical, magic |
| `eps_statistical` | `signal_effective` — the harmonic-RMS, i.e. `Σ_t 1/s(t)²` as one number | shot counts |

Evaluating once at `t_max` charged the easy times at the hard time's price and
the hard times at nothing.

**The floor was the fit.** `eps_absolute` floored at `s_res_min`, a *fitted*
late-time residual, which let the fit set its own tolerance. The floor is now
`s_abs_floor`, a stated absolute error scale, independent of any fit — and it is
the single most leveraged number in the ledger.

So the costing path uses the **data**, reduced by `s_data_rel_unc` (an assumed
20%; the source quotes no per-point error bar), inside the measured window
`t ∈ [0.1, 2]`, and the envelope times `s_extrap_factor = 0.65` outside it. The
envelope survives as the physical description and as a diagnostic;
`signal_bound = "envelope"` restores the old behaviour for comparison only.

**Effect.** Everything tightens, and the binding time moves off the endpoint —
at `m = 256` it is `t = 1.6`, not `t_max = 8`.

| n | NISQ+PEC | STAR | surface FT | Pinnacle |
|---|---|---|---|---|
| 10⁵ | 14.4 → 12.1 | 19.3 → 15.5 | 9.4 → 9.3 | 14.1 → 10.6 |
| 10⁶ | 15.6 → **13.1** | 26.5 → **22.0** | 25.8 → **12.5** | 38.0 → **26.9** |
| 10⁷ | 16.6 → 14.0 | 33.1 → 27.6 | 50.4 → 30.6 | 114.6 → 73.4 |
| 10⁸ | 17.6 → 14.9 | 39.1 → 32.6 | 222 → 125 | 305 → 204 |

Multiproduct NISQ no longer clears the ED frontier at any `n` below 10¹⁰: it
reaches 19.7 at `n = 10⁶` against 26.

**Two things this overturned.** The surface arm's large drop at `n = 10⁶` is a
*factory-ladder edge*, not an architectural result — it sits at
`p_T = 2.0×10⁻⁹` against a target of exactly `2.0×10⁻⁹`, one rung from a plant
68× larger. And the model's previous agreement with `m ~ s^{2/9}` was an
artefact: with the floor at the fitted residual 0.043 the `s = 0.03` point was
*clipped*, flattening the slope to 0.275. With the floor set independently the
clip lifts and the true local slope is 0.379. The 2/9 derivation assumed
`α = 9/4` and an `s`-independent step count; neither holds under the measured
calibration.

## 1. Circuit cost — the exponent everything rides on

With W₂ = w·m (second-order Trotter commutator sum) and t_max = √m/v_B,

```
r(m)       ~ t_max^{3/2} sqrt(W2 / eps_T)   ~ m^{5/4}
G_total(m) = c_g · m · r(m)                 ~ m^{9/4}      <- alpha = 9/4
```

`fhcost/selftest.py` asserts α = 2.250000 by log–log fit. Every downstream
exponent is 1/α = 4/9.

**ε is relative.** The Trotter and synthesis budgets are absolute, so
ε_abs = ε_rel·s. Getting this wrong inflates r by √10 and costs 1.67× in m.

**The light-cone refinement is void here.** For a local observable only the causal
cone contributes, W_eff ∝ (v t)². But at t = √m/v_B the cone *is* the lattice, so
this changes a prefactor (1/3 of the space–time box in d = 2), not the exponent.
It only bites at fixed short t — which is panel (c).

**Encoding differs by architecture.** NISQ and STAR use Derby–Klassen (1.5
qubits/mode → 3m+1, O(1)-depth layers on a nearest-neighbour grid). Surface-code
FT uses Jordan–Wigner (2m+1): under lattice surgery the circuit is a sequence of
π/8 Pauli-*product* rotations whose cost is independent of Pauli weight, and
fermionic swap networks are Clifford, hence free in T-count.

## 2. NISQ: the mitigation ceiling under a one-week cap

```
Λ(m)      = (16/15) p G_cone(m)            expected faults (exact depolarizing factor)
γ²        = e^{2Λ}                          PEC overhead; ZNE order k pays e^{2(k+1)Λ}
N_req     = T ‖c‖₁² γ² / (s ε)²             shots
N_max     = budget · (n / q_per_copy) / T_circuit
```

Solving N_req ≤ N_max gives **m ∝ (ln n / p)^{4/9}** — slower than logarithmic.
Extra qubits buy only parallel copies, and copies enter through a logarithm.

This is the slide's "NISQ saturates" box made quantitative, and it is *stronger*
than the slide claims: Takagi, Endo, Minagawa & Gu (npj QI **8**, 114 (2022))
prove the sampling overhead is exponential in Λ for **any** mitigation strategy,
so this ceiling is not an artifact of PEC. The user's constraint — any mitigation
is allowed, but the run must fit in a week — is exactly what makes it bind.

Measured result: **m = 6.7 at n = 10³ → 8.5 at n = 10⁸.** Ten decades of qubits
buy ~1.3× in lattice size. The *unmitigated* device does not reach m = 4 at any n.

### What the lower-bound literature actually gives (review #4)

Takagi, Endo, **Minagawa & Gu**, npj QI **8**, 114 (2022) — an earlier version of
this document credited the wrong two authors — prove a worst-case
estimator-spread lower bound, exponential in circuit **depth**, for a defined
class of mitigation protocols under a layered local-depolarizing noise model.

That is not what the plotted curve is. The curve is the cost of **one specified
implementation** (gatewise PEC) as a function of **gate count**, for this state
and observable. Depth and gate count scale differently in m, so the theorem does
not license the curve's m-dependence; a worst-case bound does not say this
instance is hard; and the PEC-optimality result is for a particular dephasing
setting, so PEC is not "optimal mitigation" here in any proven sense. The curve
is now labelled "PEC (as implemented)".

### Measured: which gates actually damp the observable

The cone model charges every gate in the causal cone. Measured against the
Phasecraft/Quantinuum circuit — 2415 two-qubit gates, their own quoted
`p = 1e-3`, raw against exact (FLO is exact at U=0) — the observed attenuation is
`Lambda = 0.20`, against **2.58** from the cone model. A **12.9x overcharge**.

The mechanism was checked on their data rather than assumed:
`Lambda(ZZ)/Lambda(Z) = 1.91`, where damping proportional to operator weight
predicts 2 and uniform-per-gate predicts 1. A depolarizing error damps the
observable only where the Heisenberg-evolved operator has support.
`damping_model = "support"` implements `w_obs/q` and reproduces the measurement
to 8%. It is **not** the default: one measurement fixes the value at one
operating point, not the scaling, and `support_growth = 0` cannot hold at long
times. See `OPEN_ITEMS.md` O5, and note this partially supersedes review #9.

### The 160-shot tension, resolved

Their 56-qubit run reaches ~0.005 absolute on C^zz with 160 shots per point,
where this model prices mitigation at `exp(4 p G)`. Decomposing against their
actual circuit (m = 28, t = 2):

| factor | ratio |
|---|---|
| Trotter step count — 417 here vs **4** there | **104x** |
| per-gate attenuation — 1.067 vs measured 0.083 per pG | **12.9x** |
| per-step gate count — 420 here vs 604 there | 0.70x (ours optimistic) |
| cone fraction | 0.33x |
| **net overestimate of Lambda** | **311x** |

Two by-products. `c_g = 15` gets its first external check: their compiled circuit
is 21.6 two-qubit gates per site per step, so ours is 30% optimistic but the right
order. And a `strategy = "expcal"` arm, costed at TFLO+GPR's measured effective
overhead of 0.08 — *below* one, since GPR borrows statistics across correlated
time points — and drawn only out to the largest `Lambda` demonstrated (0.20),
reaches **nothing** under this model's workload, yet reaches `m = 20` at
`n = 10^6` the moment the experiment's step density is adopted.

**So the whole distance between this model and a real experiment is the Trotter
step count.** Not the mitigation scheme, not the channel conventions, not the
damping geometry.

### A correction: attenuation is not the cancellation one-norm

These are two different quantities and an earlier version of this model used one
number, 32/15, for both. Derived from the Pauli transfer matrix of the two-qubit
depolarizing channel `D(rho) = (1-p) rho + (p/15) sum_{P != I} P rho P` (the PTM
is diagonal; `D^-1 D = 1` was checked numerically):

* 7 of the 15 non-identity Paulis commute with a given non-identity observable and
  8 anticommute, so it is damped by `1 - 16p/15` per gate. The **attenuation**
  coefficient is `-ln(1 - 16p/15)/p = 1.067236`.
* The signed Pauli inverse is `D^-1 = a . + b sum_{P != I} P . P` with
  `a = (15-p)/(15-16p)`, `b = -p/(15-16p)`, giving a **cancellation one-norm**
  `gamma = |a| + 15|b| = (15+14p)/(15-16p)`. The PEC coefficient is
  `2 ln(gamma)/p = 4.000268` — a factor **3.748** larger than the attenuation.

Conflating them made PEC **1.875x too cheap in the exponent**. Correcting it moves
mitigated NISQ at `n = 10^6` from 8.76 to **6.89**. The attenuation — and with it
the unmitigated bias and every ZNE bias ceiling — is unchanged, as it must be: a
costing convention cannot move a physical damping rate. `selftest` now asserts
that invariance directly.

*Is STAR costed the same way?* Checked, and yes. Its injected-rotation error is a
Z-type Pauli channel with RUS-total probability `q = 4p/15`; the signed inverse
gives `gamma = 1/(1-2q)` and `2 ln(gamma)/p = 1.0670`, reproducing the STAR
paper's `gamma^2 = exp(8 P_Z,1 N)` exactly. **Both arms are cancellation
one-norms.** With NISQ corrected upward STAR now beats it at every `n` it can run
(1.04x at 10^4 rising to 1.94x at 10^8), where before it lost below `n ~ 5x10^4`.
The margin is set by NISQ paying over all two-qubit gates while STAR pays only
over rotations — a real architectural difference, but one whose size rests on the
assumed ratio `c_rot/c_g = 5/15`, which has not been checked against a compiled
circuit. See `OPEN_ITEMS.md` O1.

## 3. Surface-code FT

The correction that matters most: **an FT estimate is not a qubit count.** Each of
N = T/(sε)² ≈ 2×10⁵ shots is a full circuit of r ~ 10³–10⁴ Trotter steps at d
rounds per logical layer. Checking only qubits, or only a single shot against the
week, overstates m by more than 10× and wrongly makes the curve slope 1 forever.
With the shot budget imposed the curve is

```
m = min[ qubit-limited (slope 1),  shot/time-limited (slope 4/9) ]
```

and the bend is inside the plotted range (`selftest` asserts slope → 4/9).

Other choices: d is selected against the **per-shot** logical volume (a logical
fault biases the shot it lands in; it does not compound across independent shots —
summing over 10⁵ shots over-provisions d badly). Hamming-weight phasing collapses
the T-count ~30×, so the magic-state factory is not the bottleneck.

**p_L is the largest uncertainty, so it is drawn as a band:**

| model | expression | m at n = 10⁶ |
|---|---|---|
| Fowler (idealised) | 0.1 (p/10⁻²)^{(d+1)/2} | 38 |
| Willow (measured) | 1.43×10⁻³ · 2.14^{−(d−7)/2} | 12 |

## 4. STAR

Error-corrected Cliffords plus directly injected analog rotations. The injected
rotation error is **P_Z,1 = 2p/15, independent of d** — raising the code distance
cannot fix it — giving Λ_STAR = 0.533·p·N_rot(cone). So STAR saturates
logarithmically for the same reason NISQ does.

Λ_STAR/Λ_NISQ = 1/6: STAR's noise is six times gentler per unit of circuit. Against
that it pays a full surface-code footprint per logical qubit (~600× fewer parallel
copies) and a μs·d code cycle instead of a 10 ns gate clock (~2500× slower per
shot). Both penalties enter only through a logarithm — together they cost ~7 in
Λ_max — so the 6× gain is not wiped out:

| n | 10⁴ | 10⁵ | 10⁶ | 10⁷ | 10⁸ |
|---|---|---|---|---|---|
| NISQ + PEC | 7.1 | 7.5 | 7.9 | 8.2 | 8.6 |
| STAR | 5.1 | 9.2 | 10.6 | 11.7 | 12.8 |

**STAR overtakes mitigated NISQ at n ≈ 5×10⁴** and settles ~1.4× above it. Below
that its surface-code footprint starves it of parallel copies. It is still beaten
by full FT above n ≈ 2×10⁵, and it never clears the classical band.

### Which knob makes STAR effective

Only one, and it is not the obvious one. At n = 10⁶, changing a single input:

| knob | STAR | NISQ | ratio |
|---|---|---|---|
| baseline | 10.6 | 7.9 | 1.35 |
| p = 10⁻⁴ | 28.3 | 20.9 | **1.35** |
| p = 10⁻⁵ | 71.3 | 55.3 | 1.29 |
| injection error ÷2.5 (angle-dependent RUS) | 14.8 | 7.9 | **1.88** |
| injection error ÷10 | 24.9 | 7.9 | **3.17** |
| clock: 5 logical layers/step, not 20 | 11.4 | 7.9 | 1.44 |
| 240 ns code round, not 1 μs | 11.4 | 7.9 | 1.45 |
| footprint 2d² (1 tile, not 2) | 11.0 | 7.9 | 1.40 |

**Lowering p does not help STAR relative to NISQ at all.** Both Λ's are linear in
p, so better gates lift the two curves together; the d that STAR can drop to is
worth only ½ln of a footprint ratio. Likewise the clock and the footprint are
logarithmic levers, worth ~10% each.

**The injection error is the only knob with real leverage**, because it is the only
one that enters Λ_STAR itself rather than the shot budget. The angle-dependent RUS
variant (ε_RUS,θ ≈ α_RUS·θ·p with α_RUS/k ≈ 0.40) is the concrete route: rotations
in a Trotter step are mostly small-angle, so the average injected error falls with
the angle. That is where STAR's case should be argued, not on p.

## 3b. FT accounting: workspace, magic states, and failure-to-bias (review #7)

Three things the surface-code arm was getting for free.

**Hamming-weight phasing was taking its T-count saving without paying for it.**
Phasing k identical-angle rotations computes their Hamming weight into a
ceil(log2 k)-qubit register and applies one rotation per register bit; a
phase-gradient register of `n_syn` qubits is also needed (catalytic, so paid
once). That workspace is now charged: 52-72 logical qubits across the plotted
range.

**Magic-state infidelity was never budgeted.** The model consumed T states from a
fixed `(15-to-1)` spec at every parameter point, with no check that its 4.5e-8
output was clean enough. Union-bounding over the T count:

| n | T per shot | p_T required | fixed spec gave | shortfall |
|---|---|---|---|---|
| 10^5 | 2.5e4 | 2.6e-8 | 4.5e-8 | 2x |
| 10^6 | 1.4e5 | 4.4e-9 | 4.5e-8 | 10x |
| 10^8 | 4.6e6 | 1.4e-10 | 4.5e-8 | **320x** |

The factory is now **selected** per parameter point from a ladder (cultivation,
15-to-1, cultivation + 15-to-1 cleanup, two-level 15-to-1) as the cheapest whose
output meets the requirement, and `frac_magic = 0.05` of the tolerance is
reserved for it. The model switches from cultivation to cultivation + cleanup
around n = 10^7. If nothing on the ladder is clean enough the point returns None
rather than passing silently.

**A logical failure biases by twice its probability.** A failure flips a +-1
outcome, so the distance-selection rule carries a factor 2 it did not have.

Effect at n = 10^6: surface FT falls from 49.2 to **35.1**; at n = 10^8, from 531
to **297**. The curve also acquires genuine steps where the factory protocol
changes, so its slope is no longer a clean 4/9.

### The magic plant is published operating points, not a formula (review #4)

The earlier ladder built its two-level entries from the cubic input law,
`35 p_in^3`, giving `35 (4.5e-8)^3 = 3.2e-21` at 42.6 cycles. Both numbers are
optimistic, because the cubic law describes the suppression of **input-state**
error and says nothing about faults in the distillation circuitry — whose
footprint and cycle count do not shrink when the inputs get cleaner. Litinski
determines `p_out` numerically (a five-qubit density-matrix simulation including
storage errors and faulty T measurements) and reports **4.5e-20 at 128 cycles**
for the comparable two-level protocol: 14× worse in error, 3× slower.

The table is therefore the model. Every row is a published operating point at a
stated physical error rate — thirteen from Litinski's Table 1 at `p = 1e-4` and
`1e-3`, plus magic-state cultivation at `1e-3` and `5e-4`. Nothing is
interpolated between rows and nothing is interpolated in `p`.

**Cross-check.** A `(15-to-1)_{dX,dZ,dm}` block costs `2(dX + 4dZ)·3dX + 4dm`
physical qubits, which reproduces 810, 1150, 2070 and 4620 to the table's
rounding, and takes `6 dm / (1 - p_fail)` code cycles — so the published counts
exceeding `6 dm` by 0–2% **are** the rejection rate. Rejection is already paid.
Cultivation's cost is quoted by its authors as expected volume in qubit·rounds
*including retries*, which matters: its end-to-end discard rate at `p = 1e-3` is
**99%**.

**Selection now consults `p`, and optimises the plant.** Before, an output target
of 1e-10 returned the same specification at `p = 1e-5`, `1e-3` and `3e-3`, which
cannot be read as a hardware sensitivity; and the cheapest *unit* was chosen
before asking how many units the T rate demands. Both are fixed: rows tabulated
at any `p_phys ≥ p` are admissible (pessimistic off the tabulated points, exact
on them), and the source minimising **units × footprint** wins, chosen inside the
distance loop because the unit count depends on `d`.

**Above `p = 1e-3` the model refuses.** Neither source is characterised there, so
`surface_point` returns nothing rather than extrapolating a simulated
infidelity. That is why the FT arms vanish above `1e-3` in the explorer.

**The ledger is now consistent.** A faulty T state corrupts a ±1 measurement
exactly as a logical failure does, so the magic budget carries the same factor of
two the logical check already carried. It costs a factor of two in the per-state
target, which moves the plant one rung up the ladder.

**What this changes** (`n = 10^6`, surface FT):

| p | before | after |
|---|---|---|
| 10⁻⁵ | 48.7 | **227** |
| 10⁻⁴ | 48.7 | **108** |
| 10⁻³ | 25.7 | 25.7 |
| 3×10⁻³ | 5.0 | **0 — not tabulated** |

The two identical entries at 10⁻⁵ and 10⁻⁴ were the review's complaint visible in
the top-line answer. At `n = 10^7`–`10^8` the fabricated cleanup factory is
replaced by Litinski's `(15-to-1)^6_{11,5,5} × (15-to-1)_{25,11,11}`, and m rises
slightly (49.8 → 52.6, 236 → 251) because the plant-level optimisation needs 31
units where the old model needed 52.

**Not fixed here:** Pinnacle's magic engine is not held to this ledger at all, so
it survives at `p = 3×10⁻³` where the surface code does not. That is second-pass
finding #5.

## 4b. Pinnacle (QLDPC)

Webster *et al.* (Iceberg Quantum), arXiv:2602.11457v2 (2026);
`refs/2602.11457_pinnacle.pdf`. Generalised bicycle codes with generalised
lattice surgery and a **magic engine** that distils and injects inside a single
code block.

**Footprint, reproduced with no fitted parameter.** Their published Hubbard
formula is `n = 1620 ceil((L^2+1)/8) + 4410`, which decodes as: 1620 = n_pb for
the d = 24 code [[510,16,24]]; `ceil((L^2+1)/8) = ceil((2m+2)/16)` processing
**blocks** at k = 16 logical qubits each; and 4410 for **one** magic engine.
One processing unit, one engine, no memory. Our implementation reproduces all
seven rows of their Table IV to within 1.9%, worst case, with nothing fitted.

An earlier version of this model fitted a multiplicative engine factor to their
L = 8 point and got the structure wrong in two ways, both now corrected:

* the engine was **scaled with the block count**; there is exactly one, and
* so T states were assumed to arrive at one per logical cycle **per block**. They
  arrive at one per logical cycle for the **whole machine**. On a T-heavy
  dynamics workload the old model was therefore roughly m times too fast.

The earlier explanation that the residual L-dependence came from cheaper
idle-memory blocks was simply wrong; there is no memory in their calculation.

**Effect of the correction** (n = 10^6): Pinnacle falls from 195 to **61**.

| n | Pinnacle | surface FT |
|---|---|---|
| 10^5 | 24 | 12 |
| 10^6 | **61** | 49 |
| 10^7 | 168 | 166 |
| 10^8 | 434 | **531** |

So with the engine count as published, **Pinnacle is comparable to the surface
code on this workload, not far above it, and the surface code overtakes it by
n = 10^8.** The qLDPC advantage here is in *storage*, and serialised T supply
gives most of it back. Adding engines is the obvious lever and is exposed as
`pin_engines`, but it is **our** extrapolation, not theirs, and it has an optimum
near 16: each engine costs 4410 qubits, so 64 engines is worse than 16.

**Connectivity caveat, and it is not small.** Generalised bicycle codes require
non-local qLDPC connectivity. The slide specifies a **nearest-neighbour 2D
grid**, which does not provide it. This arm is therefore costed under a different
hardware assumption from every other curve on the figure, and the comparison is
not like-for-like. Either the grid assumption or this arm has to give. See
`OPEN_ITEMS.md` O11.

**Their own Fermi-Hubbard section is a different workload** -- ground-state energy
at 0.5% relative energy error, not dynamics -- so their Table IV is a footprint
check only. The logical-error model is a refit of their Table III (A = 5.84,
threshold B = 1.58%, exponent (d+1)/2, reproducing every entry to within 26%
across twelve orders of magnitude); it is labelled a refit, not an independently
demonstrated hardware threshold.

### Hamming-weight phasing: one construction (review #6)

The workspace was a log-sized weight register plus a "phase-gradient register"
whose size was **the rotation synthesis T-count**. A T-count is not a register
size, and this was not the workspace of the construction whose T-count the model
was already charging. It happened to land near the right answer at one batch
size and was wrong by 3.6× at 256 and 5.9× at 432.

Campbell's Theorem 2 (arXiv:2012.09238, Appendix E) is explicit. A batch of `b`
identical-angle phase gates `∏ⱼ exp(iθZⱼ)` costs

```
k     = floor(log2 b) + 1     arbitrary rotations
alpha = b - w(b)              Toffoli gates AND clean ancillas
```

with `w` the popcount. Everything now comes from that one statement — workspace,
T-count, synthesis count and depth:

| | |
|---|---|
| clean ancillas | `alpha`, reused across batches |
| T gates per group | `4·(m/b)·alpha + (m/b)·k·n_syn` |
| depth per group | `(m/b)·(2⌈log₂b⌉ + n_syn)` — batches **serialise**, they share the workspace |
| synthesised rotations per shot | `Σᵢ |cᵢ| kᵢ · r · c_rot · (m/b)·k` |

`alpha` reproduces the review's 63 / 255 / 428 at `b` = 64 / 256 / 432 exactly,
and is non-decreasing in `b`, so the bisection stays monotone. `hwp_group` also
reproduces Campbell's Eq. (E16) batching map, which is the closest thing to a
compiled count available without a compiler.

**Batch size is now a real space–time knob**, because both sides of it come from
the same theorem: `b = 1` is plain synthesis (no ancillas, `m` rotations),
`b = m` is minimum T-count at maximum workspace. At `m = 256` that is 0 ancillas
and 2.4×10⁷ T gates against 255 ancillas and 1.8×10⁶. The default is the full
batch, which is the construction whose T-count the model always used — the
change is that it now pays for it. `crossovers.md` carries the trade.

**Synthesis error is allocated over the whole shot.** It was divided by the
number of rotation *groups* per step (`c_rot` = 5), ignoring the step count, the
`k` weight-register rotations each group produces, and the multiproduct branches.
The union bound now runs over every synthesised rotation, weighted by `|cᵢ|kᵢ`
across branches — 3.0× the single-circuit count at order 4, 8.2× at order 6.

**Effect.** The old formula was wrong in *both* directions, and the correction
shows it:

| n | surface FT before | after | Pinnacle before | after |
|---|---:|---:|---:|---:|
| 10⁵ | 4.9 | **9.4** | 11.5 | **14.1** |
| 10⁶ | 25.7 | 25.8 | 38.6 | 38.0 |
| 10⁷ | 52.6 | **50.4** | 126 | **115** |
| 10⁸ | 251 | **222** | 350 | **305** |

It was too *large* at small `m` — the bogus ~60-qubit gradient register dominated
— and far too *small* at large `m`, where the true workspace is linear in `m`,
not logarithmic.

### Pinnacle on the common ledger (review #5)

This arm was costed more loosely than the surface-code arm in three ways, and
the paper itself supplies what was missing.

**It took the Hamming-weight T-count saving without paying the workspace.** The
same `hwp_workspace` the surface code is charged is now in Pinnacle's logical
qubit count. Both arms compile the same circuit; they may not keep different
books on it.

**Nothing was certifying its magic states.** The model held a single 4410-qubit
engine at every point and never checked its 10⁻⁹ output against the allowance.
At `n = 10⁸` the union bound was **9.5× over** — the review's own arithmetic,
reproduced. The engine is now *selected* from the paper's four published
specifications (their Eqs. 7–10):

| p | p_out | engine qubits | reject rate | d_a, r | t_me / t_c |
|---|---|---:|---:|---|---:|
| 10⁻⁴ | 10⁻⁹ | 592 | 0.2% | 1, 1 | 14 |
| 10⁻⁴ | 10⁻¹¹ | 1807 | 2% | 5, 1 | 18 |
| 10⁻³ | 10⁻⁹ | **4410** | 10% | 7, 2 | 23 |
| 10⁻³ | 10⁻¹¹ | 5430 | 10% | 9, 2 | 26 |

At `n ≥ 10⁷` the allowance forces the 5430-qubit 10⁻¹¹ engine. Above `p = 10⁻³`
no engine is characterised and the arm refuses — as the surface code now does.
Before, Pinnacle reached `m = 13.9` at `p = 3×10⁻³` purely because nothing was
looking.

**One state per logical cycle was an assumption, not a schedule.** Their Eq. (11)
gives the distillation time `t_me = max(2d_a + 4r, t_r + 4r, d_a + t_r + 3r)`
with a reaction time `t_r = 10`, reproducing 14, 18, 23 and 26 for the four rows
— and states that this "places a lower bound on the logical cycle time of the
associated processing unit". At `p = 10⁻³` that is 23 cycles, so a processor on
the `d = 16` code, whose logical cycle is 18, **stalls**. A T state therefore
costs `max(d_t, t_me)` code cycles, divided by the acceptance rate `1 − p_r`.
That stall is live at `n = 10⁵`.

Their own 15-to-1 engine error model, incidentally, has exactly the structure
review #4 asked the surface-code plant for:

```
p_out ≈ 35 p_rot³ + 6 p_rot p_m²
```

— an input-suppression term *and* a circuit-fault term. The second term is what
the cubic law alone was missing.

**Effect** (`m` at each `n`): 20.3 → 11.5 at `10⁵`, 48.3 → 38.6 at `10⁶`,
144 → 126 at `10⁷`, 383 → 350 at `10⁸`. Pinnacle is still the strongest arm on
the figure, and on the same *error* ledger as the surface code.

**But not on the same hardware** — see §4c. Generalised bicycle codes are not
embeddable in the slide's nearest-neighbour grid, and once that is enforced the
arm is **0** there. It is plotted on a chip with two coupler layers, which is
what its codes require and what nobody has built. "Strongest arm" is therefore a
statement about a machine that does not exist.

`crossovers.md` carries the field-by-field comparison at `m = 64`. The shared
rows agree by construction. The rows that differ — 1.04M magic qubits against
5430, and 0.56 s per shot against 6.66 s — are the actual trade: qLDPC buys
storage and gives it back in serialised T supply.

## 4c. Which machine? Connectivity, parallelism and clock

Every arm above was costed on **one implicit machine**: a nearest-neighbour 2D
superconducting grid, 10 ns two-qubit gates, unlimited gate parallelism. That
was prose in a docstring, never a parameter, and it is the most load-bearing
input in the model — three of four operating points are clock-limited, not
qubit-limited. `fhcost/platform.py` makes it data, in the same style as the
magic-state tables: published operating points, sourced line by line, with a
refusal outside the tabulated range. Every number is in `refs/README.md`.

| platform | connectivity | t₂q | layer | parallel 2q | p₂q | built |
|---|---|---:|---:|---:|---:|---:|
| superconducting grid | nearest-neighbour | 10 ns | derived | ∞ | 1e-3 | 105 |
| Quantinuum Helios | all-to-all | 70 µs | **55 ms** | 4 | **7.9e-4** | 98 |
| neutral atom | reconfigurable | 275 ns | derived | 60 | 5e-3 | 60 |
| SC, two coupler layers | degree-6 | 10 ns | derived | ∞ | 1e-3 | **0** |

### The layer time is not the gate time

The thing that changed how the module is written. Helios needs ~70 µs for a
two-qubit gate and **55 ms for a circuit layer** — their own "depth-1 time",
which they call their characteristic figure of merit for processor speed. The
breakdown is Rotate 18.2 + Global Shift 7.9 + Junction 4.5 + Other Shifts 4.3 +
Split/Combine 4.1 + Four-ion Shift 1.7 + Static 0.3 ms ≈ **41 ms of ion
transport**, against 70 µs of gating: transport dominates by **786×**.

No gate count predicts ion sorting, so `Platform.layer_seconds` uses a *measured*
layer time wherever a paper publishes one and only derives one otherwise. Helios
has one; neutral atoms do not for a *reconfiguring* circuit, so that row leaves
it unset rather than inventing it.

Per-shot time is now the binding of two constraints, with a measured layer time
overriding both:

```
t_circuit = max( depth × layer_time,  G_total / n_parallel_2q × t_2q ) + t_meas
```

With `n_parallel_2q = ∞` this is algebraically the old `depth·dt_gate + dt_meas`,
so the default is unchanged — asserted, not assumed.

### Connectivity now gates the codes

`pin_nonlocal` was documented as marking an arm "costed under a DIFFERENT
hardware assumption" and **was read by nothing**. Setting it changed the model
fingerprint and no number. It is superseded by `platform.ADMITS`:

* **nearest-neighbour grid** — surface code and STAR. Generalised bicycle codes
  are *not* embeddable, so `max_m_pinnacle` returns **0**, not a small number.
* **two coupler layers** — adds BB/GB codes. Bravyi *et al.* (arXiv:2308.07915)
  put the requirement at vertex degree six with **two edge-disjoint planar
  subgraphs**, which is far weaker than all-to-all. The Pinnacle caveat had
  implied more than the codes actually need.
* **all-to-all** — everything, and the Kivlichan swap network disappears: at
  m = 256, 12 layers per Trotter step instead of 44, with the routing penalty
  going to zero.

So Pinnacle is no longer plotted against grid-costed curves as though it ran on
the grid. It is computed on the two-coupler-layer chip it requires and labelled
with it, and `curves.summary` reports both — 26.9 on that chip at n = 10⁶, and
**0 on the slide's grid**.

### They do have a plan, and it is not lattice surgery

An earlier version of this section reported "no FT arm at any n" for Helios and
neutral atoms. That was **the lattice-surgery cost model applied to machines
nobody proposes running that way**, and it is worth being precise about what is
wrong with it.

Lattice surgery needs **O(d)** syndrome-extraction rounds per logical layer,
because a measurement error has to be caught by repetition. Zhou *et al.*
(arXiv:2406.17653, *Nature* 2025) show that transversal gates with correlated
decoding need only **Θ(1)** — the deviation from the ideal logical measurement
distribution is exponentially small in `d` with a *constant* number of rounds.
This requires exactly the flexible connectivity that mobile-qubit machines have
and a planar grid does not, and it has been demonstrated on both neutral atoms
and trapped ions.

The literature is blunt about how much this matters. Existing estimates put
RSA-2048 at 20M qubits / 8 hours on superconducting with a 1 µs cycle, and
**several years** on ions or atoms with a 1 ms cycle — that "several years" is
the same calculation as our "0". With the transversal architecture, Chen *et al.*
(arXiv:2505.15907) get **19M qubits in 5.6 days at the same 1 ms cycle**, "close
to 50× speed-up … with no increase in space footprint".

`platform.se_rounds` now returns 1 on a transversal-capable machine and `d`
otherwise. At `d = 21` that is a **19× reduction** in rounds, and it changes the
answer: neutral atoms get a STAR arm at `n = 10⁶` (m = 12.8, up from 7.0) and a
surface arm from `n ≈ 10⁷`.

### So why do they still lose here? Three different reasons

Not one shared reason, which is the useful part:

| platform | what actually stops it |
|---|---|
| **Helios** | the clock. A good gate — 7.9×10⁻⁴, *better* than the slide's 10⁻³ — on a 55 ms layer. |
| **neutral atoms** | **fidelity**, not clock. Their measured 5×10⁻³ is half the surface-code threshold, so `p_L(d=41)` is 4.8×10⁻⁸ against 1.0×10⁻²² at 10⁻³ and the footprint explodes. The arm appears at `p ≲ 2×10⁻³` — a falsifiable prediction, and atom fidelities have been improving. |
| **superconducting grid** | connectivity: no qLDPC codes, and it pays the swap network. |

### And the deeper reason transversal FT does not rescue this workload

**It fixes depth overhead, not shot count.** Shor's algorithm is *one* deep
circuit run once, so removing a factor of `d` from the logical cycle removes it
from the whole runtime. This observable needs **9.1×10⁶ shots**, and the clock is
paid once per shot. At `m = 12` on atoms the week is still overrun 157× even at
an optimistic 10⁻³ gate.

That is the honest statement of the trade, and it is more specific than "ions are
slow": a shot-dominated workload punishes per-shot latency linearly, and no
fault-tolerance scheme addresses per-shot latency.

### Is the trade worth taking? No, and not close

All-to-all deletes the swap network and Helios has a *better* gate than the slide
assumes. Against that, at m = 64 one shot costs 1.2×10⁻⁵ s on the grid and
**65.8 s** on Helios. Reach at n = 10⁶:

| arm | m |
|---|---:|
| SC grid, NISQ+PEC | 13.1 |
| SC ×2 layers, QLDPC | **24.8** |
| Helios, NISQ+PEC | 7.6 |
| neutral atom, NISQ+PEC | 4.8 |
| atoms, STAR (transversal) | 12.8 |
| Helios, STAR (transversal) | 5.5 |

Transversal fault tolerance is what puts the last two rows above zero. The
connectivity is real, Helios's gate is better than the slide's, and the FT scheme
is the right one — the shortfall is per-shot latency against a shot-dominated
workload. `figures/fh_platforms.pdf` panel (b) shows
which constraint binds for each, because that is what decides whether gate speed
matters at all.

**What this does not say.** These are today's published machines on *this*
workload. A shallower circuit, a smaller shot budget, or a machine with the same
connectivity and a faster layer would move the answer, and the model computes
rather than assumes it.

## 5. The classical frontier

**This is an estimated CAPACITY of specified methods under stated machine
assumptions, not a classical impossibility boundary.** Exceeding it does not
establish that every competitive classical method fails; sitting below it does
not establish that a point is easy in practice.

**Band: m = 24-26**, set by exact diagonalisation alone (dim = C(m,m/2)^2 at half
filling, 16 B/amplitude, 4 Krylov vectors, 1-100 PB).

### U = 0: this observable is polynomial, and ED is the wrong baseline (review #7)

The classical reference was exact diagonalisation at every parameter point. At
`U = 0` that is simply the wrong method. The Hamiltonian is quadratic, so
evolution is a `2m × 2m` single-particle matrix exponential — and although the
triplet-covering initial state is **not** Gaussian (expanding it gives `2^{m/2}`
Fock branches, which is why the experimental paper calls it a magic state),
`C^zz` is a **weight-4** observable, and a four-fermion operator can connect
branches differing in at most one triplet. The branch sum therefore collapses to

* a **diagonal** part, fixed by the one- and two-mode occupation statistics of
  the branches — Wick's theorem applies branch by branch, and the sum becomes an
  average over an independent per-triplet coin flip;
* a **coherent** part from pairs differing in exactly one triplet. Ordering each
  triplet's four modes contiguously makes the Jordan–Wigner strings from the rest
  of the lattice cancel, so this is a local 4-mode object.

Neither piece needs the branches. The argument is the experimental paper's
Appendix E (arXiv:2510.26300); `calibration/free_fermion.py` is an independent
implementation of it.

**Validated, not asserted.** Against exact many-body evolution at `m = 4, 6, 8, 9`
and `t = 0, 0.25, 0.5, 1, 2`, the worst disagreement is **1.7 × 10⁻¹⁵** — machine
precision. An `O(m)` form (sparse Krylov propagator, block-collapsed sums)
reproduces the `O(M²)` form to 10⁻¹⁶ and was timed out to `m = 16384`, giving a
log–log exponent of **0.979**: linear, as claimed.

**The frontier.** At 0.70 ms per site per correlator on **one core** of
unoptimised CPython, twenty time points fit `m ≈ 4 × 10⁷` into a week. Against
the ED band's 26. That is a measured number with no extrapolation of the
implementation — a vectorised or parallel version would go further, and memory
(`O(m)`) is nowhere near binding.

So at `U = 0` no quantum arm on this figure is within six orders of magnitude of
the classical capacity for this observable, and the ED band should never have
been drawn there. The claim is confined to `U = 0` exactly; the model makes no
accuracy claim at small non-zero `U`.

**Two labelling consequences.** The figure now says *estimated ED capacity*
rather than *classically easy*, because the band is a capacity of specified
methods at a specified parameter point, not an impossibility boundary. And
estimating a low-weight observable is a different task from preparing the full
state or sampling from it — only the first is what this figure asks about.

**The TDVP arithmetic is labelled.** Dividing a `χ³` operation count by a
machine's peak rate assumes perfect strong scaling and ignores memory traffic,
communication, and the SVD/QR and MPO applications that dominate a real sweep.
`TDVP_LEADERSHIP` now carries that caveat and an explicit
`assumed_fraction_of_peak = 1.0`; the affordable `χ` is an upper bound on what
the arithmetic permits, not a run anyone has done.

### Why the tensor-network arm was removed rather than widened

The old upper edge of 62 came from an entropy argument: assume an entanglement
density, convert to a bond dimension, convert to cost. Two things killed it.

First, a normalisation error the review caught: a balanced cut of m spinful sites
carries at most **m bits**, so `ent_rate <= 1`. The old range topped out at 1.5.

Second, and decisively, the model can be checked against real data. The TDVP
series published with the experiment (Zenodo 17799843) gives |C^zz_nn| on dimer
links at U = 0 against the exact free-fermion result:

| chi | 256 | 512 | 1024 | 2048 |
|---|---|---|---|---|
| mean abs. error | 0.100 | 0.085 | 0.078 | **0.077** |

Pooled over t in [0.5, 2] the error is **flat in chi** -- `err ~ chi^-0.12`, so
doubling the bond dimension buys 9%.

**But flatness in chi does not mean tensor networks fail, and the pooled number
hides what is going on.** Splitting the error by time rather than averaging it
(`calibration/tdvp_check.py`, all four chi at every time) gives a slope
`d log err / d log chi` that is **non-monotonic in time**:

| t | 0.1 | 0.3 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|---|---|
| slope | **-0.03** | -0.30 | **-0.45** | -0.37 | -0.26 | -0.07 | **-0.04** |
| mean abs. err at chi = 2048 | 7.4e-3 | 1.6e-2 | 2.6e-2 | 3.3e-2 | 3.0e-2 | 9.9e-2 | 1.6e-1 |
| exact mean \|C^zz\| | 0.942 | 0.587 | 0.243 | 0.093 | 0.075 | 0.036 | 0.034 |

Through the middle of the window the run *is* bond-dimension-limited: 8x the chi
buys 2.5x the error at t = 0.5. It is the two ends that are pathological, for
opposite reasons.

At **t = 0.1** the state is barely entangled and chi = 256 is wild overkill, yet
the error is 7.9e-3 and an *eightfold* increase in bond dimension removes only 8%
of it (7.35e-3 at chi = 2048). A truncation-limited calculation would be at
machine precision there. The deposit's `max_bond_dimension` column confirms every
run saturated its cap, so truncation was binding -- it just was not what limited
the accuracy. Plausibly two-site TDVP projection error on a snake MPS with
long-range Jordan-Wigner strings, or the time step, or the GPR smoothing. The
deposit varies chi and **never varies dt**, so it cannot distinguish them.

At **t >= 1.5** the chi = 2048 error is about four times the mean |C^zz| there.
A calculation whose error exceeds its answer bounds nothing.

The published TDVP is therefore **not a converged tensor-network calculation**,
any chi-extrapolation from it is meaningless, and it cannot bound what tensor
networks can do on this problem. The t = 0.1 floor sits at roughly our entire
absolute tolerance, so a clean implementation that removed it could plausibly
reach this accuracy at modest chi.

At the same operating point the entropy model predicts `chi ~ 6.6e3` suffices.
It is not mis-calibrated, it is **the wrong shape**: an entropy argument cannot
represent an error that saturates in chi. So it was removed rather than refitted.

### What this does and does not license

Removing it narrows the band from 24-62 to 24-26, and that narrowing is what
flips "mitigated NISQ never clears classical" into "clears at n ~ 1.7e4". **The
band narrowed because an unvalidated arm was deleted, not because tensor networks
were shown to fail** -- a change in the direction that flatters the quantum
curves, and one to be suspicious of accordingly.

That TDVP run was a comparison baseline, not a best-effort classical attack; the
error plateau looks systematic rather than truncation-limited, which would void
the chi-extrapolation entirely; a snake MPS on a doubly periodic 7x4 torus is
close to worst-case geometry; and no PEPS, neural-quantum-state or Pauli-path
attempt exists for this observable at this accuracy.

**So the defensible claim is "clears the exact-diagonalisation frontier", not
"clears classical".** See `OPEN_ITEMS.md` O10, which is the item most likely to
reverse a headline.

### Noise-induced classical simulability

Considered and **rejected** as a band: at p = 1e-3 and depth ~1e4 the required
Pauli-weight truncation is far larger than the register, so the published
quasi-polynomial results are asymptotically true but numerically vacuous here.

## 6. The three extrapolations

**(a) Noise extrapolation.** PEC and Richardson ZNE fail for opposite reasons and
must not be costed with the same formula.

*PEC* samples the inverse-noise quasiprobability decomposition: signed, weighted
samples, per-shot variance `gamma^2 = exp(2 Lambda)`, but an **unbiased**
estimator. Cost is exponential in `Lambda`, and shots genuinely buy accuracy.

*Richardson ZNE* runs the circuit at amplified noise `lambda_i` and averages
**bounded** outcomes. Its variance factor is only `(sum_i |c_i| sqrt(lambda_i))^2`
under the optimal allocation `N_i ~ |c_i| sqrt(v_i / tau_i)` (with `tau_i` the
longer runtime of the noise-amplified circuit) — polynomial, not exponential.
What it pays instead is **residual bias**: the extrapolant recovers the
noiseless value only to `|1 - sum_i c_i exp(-lambda_i Lambda)|`, and **no number
of shots removes that.** With half the relative tolerance reserved for bias:

| strategy | max Lambda before bias blows the allowance | variance factor there |
|---|---|---|
| unmitigated | 0.05 | 1 |
| ZNE order 1 | 0.21 | 5.6 |
| ZNE order 2 | 0.35 | 23.8 |
| ZNE order 3 | 0.46 | 94.8 |
| PEC | unbounded (unbiased) | `exp(2 Lambda)` |

**Higher ZNE order buys MORE noise headroom, not less** — an earlier version of
this model asserted the opposite, and even encoded it as a test. It was wrong:
that claim came from charging ZNE an exponential variance
`sum |c_i| exp(lambda_i Lambda)`, i.e. assuming attenuation had to be inverted at
every node, while never checking the bias at all. Those two errors together
reported ZNE points carrying **~0.9 relative bias against a 0.1 target** — points
that are infeasible at any shot count.

The consequence is blunt. At `p = 10^-3` the smallest real lattice, `m = 4`, already
needs `Lambda = 1.55`, which is **above every ZNE order's bias ceiling**. So under
this response model **Richardson ZNE cannot do even a 2x2 Hubbard lattice**, and
the ZNE curves are removed from the figure rather than drawn as feasible. ZNE
becomes usable again once the noise is low enough (it works at `p = 10^-5`).

*Caveat.* `s(lambda) = s(0) exp(-lambda Lambda)` is a test model, not an
established description of this observable under gate-local noise. A response
with less curvature would raise every bias ceiling above, and with it every ZNE
conclusion here. Measuring the actual noise response of a small compiled Hubbard
circuit is the check that would settle it.

**(b) Trotter-step extrapolation.** Order-2k multiproduct formulas, costed as
**classical extrapolation of expectation values**: branch i is its own circuit at
k_i times the base step count.

*The convergence is validated.* Against exact diagonalisation the order-4 formula
measures order **3.5-4.1** and order-6 measures **5.7-6.5**, with error gains up
to 5x10^6 over plain Trotter at the same base step count. The review's doubt
about the higher-order remainder was unfounded -- the formulas work.

*The costing was not.* Branch i has k_i times the gates, hence k_i times the PEC
exponent, and k_i times the runtime. Charging only `||c||_1^2` ignores all of
that, and the deepest branch dominates. Correct accounting uses the optimal
allocation over branches, minimum total time
`(sum_i |c_i| sqrt(v_i tau_i))^2 / delta^2` with `v_i = exp(k_i log Gamma^2)` and
`tau_i ~ k_i`. Effect at n = 10^6:

| order | NISQ + PEC | (as ||c||_1^2) | surface FT | (as ||c||_1^2) |
|---|---|---|---|---|
| 2 | 29.3 | 29.3 | 231 | 231 |
| 4 | **41.3** | 62.5 | **266** | 266 |
| 6 | 41.1 | 80.4 | 136 | 249 |
| 8 | 38.0 | 89.1 | 63 | 143 |

So multiproduct is worth about **1.4x for NISQ and 1.15x for FT**, not the ~3x
the old charge implied, and there is a genuine **optimum at order 4-6** rather
than a monotone gain.

The size of the undercharge itself depends on the Trotter model, because it
scales with `Lambda`: **11-14x** under the loose commutator bound but only
**2.7-6.7x** under the measured calibration, since shorter circuits mean a
smaller PEC exponent for the deepest branch to amplify. An earlier estimate of
13-280x was computed under the loose bound and does not survive the calibration. FT falls off faster than NISQ because it has no PEC
overhead to amortise, so the deepest branch's runtime is the whole cost.

**(c) Finite-size extrapolation.** For a local observable, once the lattice is
bigger than the boundary disturbance has travelled the finite answer *is* the
infinite answer. So m has a **ceiling you hit, not a runway**: above

    m_required(t) = [ 2 ( v t + xi ln(1/eps_fs) ) ]^2

extra sites buy literally nothing, and every further qubit should be spent on t.
Fitting a ladder of lattice sizes rather than waiting the exponential out shifts
the threshold down by xi ln(gain) -- 31% in m at t = 1 -- a constant shift in L,
not a change of scaling, paid for with K x the shots.

**But this criterion is far weaker than it looks, and it is conditional on the
initial state.** Three separate caveats, none of them small:

1. *It is a quench statement.* The exponentially-small-outside-the-cone tail is
   clean for evolution from a low-entanglement product/Neel state, where
   Calabrese-Cardy quasiparticle pairs are emitted at t = 0. From a **ground or
   thermal** state the lattice already carries equilibrium correlations of range
   xi_eq before any evolution; near half filling at low T the 2D Hubbard
   antiferromagnetic xi_eq grows, that term dominates, and the criterion
   understates m_required badly.
2. *The velocity is the butterfly speed, not the correlation front.* Boundary
   contamination is a single disturbance propagating inward at v_B; it does not
   need a quasiparticle pair. Charging 2 v_B (as an earlier version of this model
   did) makes the ballistic term twice too large.
3. *The rigorous Lieb-Robinson speed is ~10x v_B* and would inflate m_required by
   ~100x. The tight quasiparticle version used here is the **optimistic edge of a
   band, not a bound.**

| t | v_B, xi=1 (default) | 2v_B, xi=1 | v_B, xi=3 (thermal-like) | Lieb-Robinson |
|---|---|---|---|---|
| 0.5 | 185 | 244 | 1358 | 1000 |
| 1.0 | 244 | 385 | 1510 | 2664 |
| 2.0 | 385 | 763 | 1837 | 8394 |

**An 11x spread at t = 1.** The structure of the result -- a step rather than
growth -- is robust; its location is not. And the initial state is still
unspecified (see 0(b)), which is now the blocking ambiguity for this whole
argument. `fs_speed` and `xi` are knobs; m_required should be drawn as a band.

**Under slide 1's own convention the step is unreachable anyway.** With
t_max = sqrt(m)/v_B the boundary disturbance has always crossed the lattice by
construction, so no point on the m-vs-n plot is ever converged and finite-size
extrapolation buys exactly nothing. That is a property of the time window, not
of the method, and it is the strongest argument for changing it.

## 6b. The error ledger

Every contribution is a share of the **same** absolute tolerance, and the shares
must sum to at most 1. They did not: Trotter, synthesis and logical each took
0.30 while sampling took 0.50 (NISQ) or the **whole** budget (FT) -- 1.4x and
1.9x over respectively, and inconsistent between the two arms, which quietly
tilted every NISQ-vs-FT comparison toward FT.

| contribution | share |
|---|---|
| Trotter | 0.25 |
| gate synthesis | 0.10 |
| logical (per shot) | 0.10 |
| residual mitigation bias | 0.05 |
| statistical half-width | 0.50 |
| **total** | **1.00** |

Synthesis and logical get the small shares because both enter the resource only
logarithmically, so buying them down is cheap.

**Statistics are now a confidence half-width, not a 1-sigma spread.** The old
`1/(s eps)^2` was a ~68% statement. The guarantee is **per-time two-sided 95%**,
`z = 1.96`; `simultaneous = True` applies a Bonferroni correction across the time
grid (`z = 3.02` at 20 points, 2.38x the shots) if the claim is about the whole
curve rather than each point.

**The connected subtraction is free here.** `C^zz = 4(<S_i S_j> - <S_i><S_j>)`
estimates three correlated quantities from the same shots, but in the dimerised
triplet state `<S^z_i> = 0` by symmetry, so the disconnected term vanishes and
the estimator's variance collapses to the plain bounded-observable value. That is
a property of this state, not a general result.

Cost of getting this right, at n = 10^6:

| | before | after |
|---|---|---|
| NISQ + PEC | 29.3 | **26.8** |
| STAR | 58.5 | **42.6** |
| surface FT | 230.8 | **49.2** |
| Pinnacle | 707 | **195** |

FT and Pinnacle fall hardest because they were the ones spending the entire
tolerance on statistics. Two consequences worth naming: multiproduct **stops
helping FT altogether** (order 2/4/6/8 gives 49/45/23/8 -- the extra branches
never pay for themselves once FT is properly shot-limited), and ZNE lands on a
knife edge, `Lambda(m=4) = 0.222` against an order-3 ceiling of 0.217. **ZNE's
viability here is set by how much of the budget its residual bias is allocated,
not by the physics** -- at the old 0.5 share it runs, at 0.05 it does not.

### Zero band edges are a statement, not missing data (review #9)

Scenario bands were NaN-masked at any `n` where an edge was zero, so the shading
vanished exactly where the model is most pessimistic. On the default
configuration this was not an edge case: the **NISQ+PEC band had a zero lower
edge at every one of the 60 plotted points**, so no shading was drawn for that
arm at all. The Willow arm on the slide figure lost half its range the same way.

A zero edge says the commutator-bound scenario does not reach even a 2×2
lattice. That is the most informative thing the band has to say. `band_for_plot`
now clips the lower edge to the axis floor, returns the clipped region as a
mask, and the figure hatches it — so the band visibly runs off the bottom of the
axis under a dashed `m = 4` line. Where the *upper* edge is zero too, nothing is
feasible, both edges come back NaN, and the panel annotates instead of drawing a
band that does not exist.

The shading is also now labelled for what it is, at the foot of the panel rather
than in a legend this style does not have: **a scenario span between two Trotter
models, not a confidence interval.** Neither edge is a guaranteed bound and the
truth is not required to lie between them.

`check_figure.py` renders these cases and measures the filled area of each
polygon, because the logic test and the ink are different things: the
NaN-masked version renders **0 px²** where the clipped version renders 13402.

## 7. What the figure changes about the slide

1. The sketch's "logical advantage = where FT overtakes NISQ" marks the wrong
   crossing. NISQ+PEC and STAR both sit *below* the classical band at every n, so
   FT-vs-NISQ is a race between two losers. The meaningful crossing is
   **FT vs classical, at n ≈ 8×10⁵** (idealised p_L) or ≈ 1.5×10⁷ (measured p_L).
2. The unmitigated p = 10⁻³ curve does not exist — it never reaches m = 4.
3. STAR does not sit above NISQ and FT at p = 10⁻³; it is a wash with NISQ.
4. Surface-code FT starts later than sketched (factory bank dominates below
   n ≈ 10⁵) and then rises faster, but bends to slope 4/9 once shot-limited.
5. The time-window convention moves the answer more than any other input.

## 8. Inputs the literature does not pin down — parameterized, not guessed

Exposed in `fhcost/budget.py`; swept in `crossovers.md`.

1. **s = |⟨Z_i(t)Z_j(0)⟩_c|** at the light-cone front. Enters twice (shots ~ 1/s²,
   Trotter r ~ s^{−1/2}), net m ∝ s^{2/9}. Unpublished for 2D FH at t = L/v.
   *The single largest un-pinned number in the model.*
2. **α, entanglement density** at t_max. Sets the classical band, m = 12→58.
3. **κ, empirical-vs-worst-case Trotter ratio** for 2D FH, somewhere in [1, 100].
4. **p_L(d)** — idealised vs measured, a 3× spread in m_FT.
5. **v_B** — we use 2J (free-fermion max axial group velocity); the correlation
   front moves at ~2v_B = 4J. The rigorous Lieb–Robinson bound is ~10× looser and
   is deliberately not used.
6. **t_round** — 1 μs (decoder-latency-limited, as Google runs) vs the 240 ns the
   hardware spec alone implies.
7. **Whether PEC is calibratable at Λ ≈ 8.** This is an information-theoretic cap
   from the shot budget, not a demonstrated capability; sparse Pauli–Lindblad PEC
   has been shown at γ² ~ 10–10², not e¹⁶. The NISQ curve is an upper bound.

## 7b. Which integer lattices can actually hold the state (review #10)

`floor(sqrt(m))` is not enough, and the headline it produced was wrong. The
specified state is half-filled at `S^z_tot = 0`: one holon, one doublon, and
`S^z = 0` triplets covering every remaining site in pairs. Two constraints
follow, and a 5×5 fails both:

1. **The site count must be even.** Each triplet contributes one up and one down,
   and so does the doublon, so `N_up = (N−2)/2 + 1 = N/2`. On 25 sites that is
   not an integer.
2. **The remainder must admit a perfect dimer covering.** On 25 sites, 23 remain
   after the defects — odd, so no covering exists at all.

Both are satisfiable exactly when `N` is even. A grid is bipartite, `N` even
forces at least one side even, which makes the two colour classes equal; putting
the holon and doublon on **opposite sublattices** then leaves equal classes, and
Gomory's theorem gives a perfect matching of the rest. So the rule is *even site
count, defects on opposite sublattices*, and the defect placement is ours.

Squares are not required — the experiment itself uses 7×4 — so rectangles are
admitted, with the most square winning on a tie and an **aspect cap of 2**:
beyond that a ribbon is quasi-1D, which is a different and materially easier
problem.

The correction runs in both directions. `floor(sqrt(m))` named an *impossible*
lattice at four of five sizes checked, and it also **understated** the capacity:

| capacity m | floor(sqrt) | admissible |
|---:|---|---|
| 13 | 3×3 = 9 — impossible | **3×4 = 12** |
| 25 | 5×5 = 25 — impossible | **4×6 = 24** |
| 63 | 7×7 = 49 — impossible | **6×10 = 60** |
| 125 | 11×11 = 121 — impossible | **10×12 = 120** |

At `n = 10⁶` the largest admissible lattices are NISQ 3×4, STAR 4×5, surface FT
3×4, Pinnacle 4×6. `max_integer_L` is kept only so an old caller fails loudly.

## 7c. The additional implementation checks

Five smaller items from the same review, each verified against the code before
being acted on.

**The error ledger was enforced only as a side effect.** It is checked inside
`eps_absolute`, so every costing path did hit it — `frac_stat = 0.9` is rejected
at all six public entry points. But a precondition should not depend on being
reached through a tolerance calculation, so `hubbard.validate(cfg)` is now
explicit and also rejects `eps ∉ (0,1)`, `p` at or above threshold,
`n_times < 1`, and a zero absolute floor.

**`s_sig` was inert.** Under the default `signal_regime = "curve"` it is not
used, so the "weak signal" and "strong signal" sensitivity rows were reprinting
the baseline. Those rows now vary the fixed-signal *scenario* (where `s_sig` does
bite: 11.5 → 26.8) plus the parameters the curve actually uses —
`s_data_rel_unc`, `s_abs_floor`, and the superseded envelope bound.

**The cluster prose contradicted the implementation.** It claimed the buffered
cluster expansion is uncompetitive "at any xi in the plausible range (checked
from 0.2 to 1.0)". It is not: `cluster_t_reach` gives 0.368 at `xi = 0.2` and 0
at 0.5 and 1.0. The buffer is `xi ln(1/eps)`, so a short correlation length makes
the cluster cheap. The surviving claim is narrower and is now what the docstring
says.

**Two different questions were being read as one.** `max_m_fixed_t(0.02) ≈ 10⁵`
sites and `classical_t_reach ≈ 0.013` are both correct: the first is a *capacity
at fixed t* from an entropy-derived bond dimension with no convergence
guarantee, the second is the time to which a *bound certifies* a
thermodynamic-limit answer. `max_m_fixed_t_detail` now returns the method and
the non-claim alongside the number.

**Tests that assert narrative.** The review is right that a test deriving its
expectation from the same formula it checks catches regressions but is not
independent validation. Where such tests have been touched in this pass they
have been rewritten to test a mechanism against a measured or published number —
Campbell's Theorem 2, Litinski's Table 1, the Pinnacle Eq. (11) timings, the
`(ln n)^{4/9}` law, exact many-body evolution at `U = 0`, rendered polygon area.
The ones not yet touched remain regression tests and are not evidence.

## 8b. Hygiene (review #11)

Eleven implementation and reporting defects, all confirmed against the code
before fixing:

| # | defect | resolution |
|---|---|---|
| 1 | the two finite-size functions used different velocities and buffers (255 vs 196 at t=1, gain=1), and `n_sizes` was ignored | one law; `gain=1` now recovers `m_required` exactly, and the size ladder returns its shot multiplier |
| 2 | `classical.band` overrode the supplied RAM, so RAM sensitivity rows were decorative | `band` is documented as a FIXED 1-100 PB envelope; `ed_at_cfg_ram` is the single-machine sensitivity |
| 3 | STAR was documented as compact-encoded but coded as Jordan-Wigner | the **code** was right (lattice surgery makes Pauli weight free); the doc is fixed |
| 4 | Python ignored `lanes` for surface/STAR while the explorer divided by it | neither should: the factory bank is already sized for throughput, and STAR injects locally. The explorer now matches |
| 5 | fractional replicas | integer packing. Leftover qubits below one full copy cannot run anything, so the staircase is real |
| 6 | `crossing` reported "never" having searched only to 1e10 | reports "not reached below n = 1e10" |
| 7 | continuous m presented as an actual square lattice | `max_integer_L()`; at n = 10^6 the feasible lattices are 5x5, not 26.8 sites |
| 8 | README called 1e3-1e8 "ten decades" (it is five) and quoted stale crossovers | all headline figures generated by `headlines.py` into `HEADLINES.md` |
| 9 | the Willow p_L has no `p` in it, so its p-sensitivity column was fixed by construction | labelled a fixed experimental anchor |
| 10 | `m ~ s^(2/9)` claimed generally | **measured**: 0.25 for noise-limited NISQ (so 2/9 was right there) but **1.37** for shot-limited FT. The general claim was wrong |
| 11 | one-qubit, idle, preparation and measurement noise omitted | stated as a restricted channel model; `noise_channels = 1.89` adds them |

## 8c. One model, one record (second-pass review #3)

The Python package, the explorer and the prose had drifted into three models
reporting different numbers at the same nominal settings — at `n = 10⁶`, NISQ+PEC
was 16.8 in Python, 7.9 in the explorer and 27 in this repository's README. Not
rounding: different PEC coefficients, different Trotter models, a fixed signal
against a measured decay, and a superseded Pinnacle construction.

`fhcost/` is now the single model and everything else is generated from
`RESULTS.json`:

| artefact | what is generated |
|---|---|
| `explorer.html` | the whole constants block, straight out of `Config`, plus the probe grid it checks itself against |
| `HEADLINES.md`, `README.md` | the headline numbers, between generated markers |
| figure, slide, `crossovers.md`, explorer header | the model fingerprint |

**The fingerprint** is a hash of the model source *and* the default configuration,
so editing either invalidates every stamped artefact. Two artefacts carrying
different ids came from different models and must not be compared.

**The port is checked, not trusted.** Every JavaScript function keeps its Python
name — `eps_absolute` here is `fhcost.hubbard.eps_absolute` there — and the page
recomputes the probe grid on load: eight configurations (default, `p = 10⁻⁴`,
`ε = 0.01`, fixed `τ`, MPF order 6, commutator bound, `U/J = 8`, and the
independent model's full assumption set) × six lattice sizes × sixteen
intermediates, plus eight qubit counts × seven architectures, plus the classical
band. Tolerance 10⁻⁹ on intermediates, 10⁻⁷ on the maxima. The verdict goes in
the page header; `check_parity.py` reads it back out of a headless Chrome and
exits non-zero on a mismatch. Perturbing one line of the port by 10⁻⁷ trips 209
probes.

`python3 make_record.py --check` closes the loop: a model change that is not
followed by a regeneration fails instead of leaving a stale number in a document.

## 9. References

Fowler, Mariantoni, Martinis & Cleland, PRA **86**, 032324 (2012) ·
Google Quantum AI, Nature **638**, 920 (2025) ·
Litinski, Quantum **3**, 128 (2019) and **3**, 205 (2019) ·
Gidney, Quantum **2**, 74 (2018) ·
Ross & Selinger, QIC **16**, 901 (2016) ·
Campbell, QST **7**, 015007 (2021) ·
Childs, Su, Tran, Wiebe & Zhu, PRX **11**, 011020 (2021) ·
Kivlichan et al., PRL **120**, 110501 (2018) ·
Derby, Klassen, Bausch & Cubitt, PRB **104**, 035118 (2021) ·
Temme, Bravyi & Gambetta, PRL **119**, 180509 (2017) ·
van den Berg, Minev, Kandala & Temme, Nat. Phys. **19**, 1116 (2023) ·
Takagi, Endo, Minagawa & Gu, npj QI **8**, 114 (2022) ·
Quek, França, Khatri, Meyer & Eisert, Nat. Phys. **20**, 1648 (2024) ·
Akahoshi et al., PRX Quantum **5**, 010337 (2024) ·
Low, Kliuchnikov & Wiebe, arXiv:1907.11679 · Vazquez et al., Quantum **7**, 1067 (2023) ·
Schuster, Yin, Gao & Yao, PRX **15**, 041018 (2025) ·
Aharonov, Gao, Landau, Liu & Vazirani, STOC 2023 ·
Yoshioka et al., npj QI **10**, 45 (2024).

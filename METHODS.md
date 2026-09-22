# Methods: quantitative resource estimate for 2D Fermi–Hubbard dynamics

Replaces the qualitative whiteboard sketch on slide 2 of `Q1_Fermi_Hubbard-1.pdf`.
Every number in `figures/fh_resource_estimate.pdf` and `crossovers.md` is produced
by `fhcost/`; nothing is hand-placed. Reproduce with

```
python3 -m fhcost.selftest     # 40 checks, all arithmetic, <1 s, ~27 MB
python3 make_figure.py
python3 crossover_table.py
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

## 4b. Pinnacle (QLDPC)

Webster, Berent, Chandra, Hockings, Baspin, Thomsen, Smith & Cohen (Iceberg
Quantum), arXiv:2602.11457v2 (2026); local copy `refs/2602.11457_pinnacle.pdf`.
Generalised bicycle codes with generalised lattice surgery, one **magic engine**
per processing unit (distillation and injection inside a single code block,
one |T> per logical cycle, so there is no separate factory bank), and Clifford
frame cleaning for parallel memory access.

Their Table I gives the code family; `dt = d + 2` code cycles per logical cycle:

| [[n,k,d]] | n_pb | physical per logical | surface `4d²` | advantage |
|---|---|---|---|---|
| [[30,8,4]] | 140 | 17.5 | 64 | 3.7× |
| [[126,12,10]] | 452 | 37.7 | 400 | 10.6× |
| [[254,14,16]] | 860 | 61.4 | 1024 | 16.7× |
| [[510,16,24]] | 1620 | 101.2 | 2304 | 22.8× |

Their logical error ansatz is `p_L = (A/k)(p/B)^(d/2+C)`. Fitting their Table III
gives **C = 1/2** — the same exponent form as the surface code — with
**A = 5.84 and a threshold B = 1.58%**. The fit reproduces every Table III entry
to within 26% across twelve orders of magnitude (`selftest` asserts this). The
three sources of advantage are separable: a better threshold (1.58% vs 1%), the
`1/k` sharing across a block, and `n_pb/k` instead of `4d²` per logical qubit.

**Their own Fermi–Hubbard section is a different workload from ours** — ground
state energy at 0.5% relative energy error, not dynamics — so their Table IV is
not directly comparable. What is borrowed is the *architecture*. Their
logical-qubit count `N = 2L² + 2` is the same Jordan–Wigner count we already use.
Our footprint model reproduces their Table IV exactly at L = 8 (19 kq, the
calibration point for the one free parameter, the magic-engine factor) and runs
14–26% high at larger L, because they place idle data in cheaper memory blocks
that we do not model. Our Pinnacle numbers are therefore mildly conservative.

Running **our dynamics workload** on it:

| n | MASQ+PEC | STAR | surface FT | Pinnacle |
|---|---|---|---|---|
| 10⁵ | 7.5 | 9.2 | 14 | **52** |
| 10⁶ | 7.9 | 10.6 | 38 | **97** |
| 10⁷ | 8.2 | 11.7 | 99 | **264** |
| 10⁸ | 8.5 | 12.8 | 246 | **712** |

2.5–3.7× more lattice at fixed n, and — the number that matters — it clears the
classical band at **n ≈ 2.9 × 10⁵ instead of 2.3 × 10⁶, eight times earlier**.
In the converged framing it passes the classical light-cone limit (t ≈ 1) at
n ≈ 10⁶ where the surface code needs ≈ 10⁷.

The gain here is smaller than their headline 10× because that figure is for
RSA-2048, a far deeper and more T-heavy circuit where the rate advantage
compounds over a much longer schedule. On a shallow dynamics workload at modest
m, block granularity (k = 14–16 logical qubits per block) wastes some capacity.

## 5. The classical frontier

**This is an estimated CAPACITY of specified methods under stated machine
assumptions, not a classical impossibility boundary.** Exceeding it does not
establish that every competitive classical method fails; sitting below it does
not establish that a point is easy in practice.

**Band: m = 24-26**, set by exact diagonalisation alone (dim = C(m,m/2)^2 at half
filling, 16 B/amplitude, 4 Krylov vectors, 1-100 PB).

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

The error is **flat in chi** -- `err ~ chi^-0.12`, so doubling the bond dimension
buys 9% -- and it has plateaued at 0.077, which at late times is several times
larger than the signal itself (exact |C^zz| = 0.021-0.031 for t >= 1.5).
Extrapolating that slope to our tolerance would need `chi ~ 2e12`.

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

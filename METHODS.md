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

## 5. The classical band

Drawn horizontally because the x axis is qubits and classical compute does not
live on it. At t_max = √m/v_B the light cone has crossed the lattice, which kills
every structural shortcut at once:

| method | reach | why |
|---|---|---|
| Krylov / state vector | m = 24 (1 PB) – 26 (100 PB) | dim = C(m,m/2)², memory-bound |
| snake MPS / PEPS | m = 12–58 | S ≈ α·m ⇒ χ ~ 2^{αm}; χ³ runtime binds, not memory |
| light-cone cluster | m ≈ 24 | cone is ~4m sites > m; strictly worse than ED |

**Band: m = 24–56.** The spread is dominated by α, the entanglement density per
site — the least-pinned classical input. Published 2D Hubbard ED sits near m = 20.

**On noise-induced classical simulability:** Pauli-path methods give quasi-poly
simulation of noisy circuits at fixed p, which would swallow the unmitigated-NISQ
region. We considered drawing this as a second band and **rejected it**: at
γ = 10⁻³ and depth ~10⁴ the required Pauli-weight truncation is ℓ ≈ 7000 ≫ 3m, so
the theorem is asymptotically true but numerically vacuous at these sizes. The
simpler statement is also the stronger one — the entire mitigated-NISQ and STAR
region sits below the ED band anyway.

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

**(b) Trotter-step extrapolation — the largest single lever.** Order-2k
multiproduct formulas give r ∝ t^{1+1/2k}(W/ε)^{1/2k}. The shot-noise cost is the
coefficient 1-norm, and for well-conditioned MPFs ‖c‖₁ = O(log k) — polylogarithmic,
not exponential (our computed ‖c‖₁ = 5/3 at order 4 reproduces the published
value). Charged net of that cost:

| | order 2 | 4 | 6 | 8 |
|---|---|---|---|---|
| NISQ + PEC | 8 | 31 | 55 | 73 |
| surface FT | 38 | 140 | 163 | 107 |

**There is an optimum, and it differs by architecture.** NISQ is Λ-limited, so
depth reduction keeps paying. FT is shot-limited, so beyond order ~6 the ‖c‖₁²
shot cost overtakes the depth saving. This is the single biggest algorithmic win
available and it moves NISQ+PEC from "hopeless" to the edge of the classical band.
Caveat: MPFs at this ε and t have not been demonstrated for 2D FH — treat as an
upper bound on the gain.

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

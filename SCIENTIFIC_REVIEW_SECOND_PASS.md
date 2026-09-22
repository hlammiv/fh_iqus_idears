# Scientific review: second pass

Review date: 2026-09-22.

This document records a new examination of the current problem definition,
implementation, outputs, and selected primary references. The previous review
document was not used as a checklist or opened during this pass. The earlier
conversation remained in context, so this was not a fully blinded review.

All existing arithmetic self-tests passed. Additional checks included direct
evaluation of the initial dimer correlation, architecture operating points,
factory choices, calibration domains, and the JavaScript explorer's calculation
functions. No model code or existing results were changed during the review.
This file records findings and proposed work; it does not implement the fixes.

## Assessment and priorities

The current crossover results are not yet scientifically reliable. The main
problems are a convergence criterion used as a necessary condition, extrapolation
of empirical calibration beyond its evidence, and unequal resource and error
accounting between fault-tolerant architectures.

| Priority | Finding | Consequence |
|---|---|---|
| Critical | Finite-size criterion fails the exact initial-state limit | The headline that nothing classical converges is false |
| High | Default Trotter calibration is extrapolated in both size and time | The new resource scaling and operating points are assumptions, not calibrated predictions |
| High | Python and HTML implement different models | The user receives incompatible results depending on the interface |
| High | Factory cleanup neglects circuit faults | FT fidelity guarantees and hardware sensitivities are unsupported |
| High | Pinnacle omits resources charged to surface FT | Architecture rankings are not based on a common accuracy/resource ledger |
| High | Hamming-weight workspace does not match the cited construction | Qubit count, T-count, and depth do not describe one demonstrated circuit |
| High | Classical baseline misses a tractable limit | Parameter comparisons can imply difficulty where a specialized classical method applies |
| High | Signal fit does not support per-time relative accuracy | Shot budgets can be too small near signal minima |
| Medium | Plot masks zero-valued uncertainty edges | The figure hides the scenario in which NISQ cannot run the smallest lattice |
| Medium | Integer lattice reporting ignores state constraints | Reported feasible sizes can be incompatible with the specified initial state |

## 1. The convergence criterion fails an exact physical limit

**Locations:** [converged.py](fhcost/converged.py), `m_required_extrapolated`,
`classical_t_reach`; [headlines.py](headlines.py); [HEADLINES.md](HEADLINES.md).

The code predicts approximately **143.59 sites at t = 0** as the required size
for a converged answer. Since ED capacity is estimated as 26 sites, it concludes
that classical methods cannot produce converged answers.

The specified initial state is a product of triplet dimers, apart from the
holon/doublon defects. On an intact dimer,

```text
|psi_dimer> = (|up,down> + |down,up>)/sqrt(2)
C^zz_ij(0) = 4[<S^z_i S^z_j> - <S^z_i><S^z_j>] = -1.
```

This is exact on two sites. Adding distant dimers does not change the answer.
A direct two-site calculation during this review returned -1 to floating-point
precision, with zero variance up to rounding.

The conceptual error is treating a **sufficient system size from a conservative
error estimate as a necessary minimum size**. Failure to satisfy a bound does not
prove failure of convergence. The bound also does not reproduce the vanishing
boundary-evolution error at t = 0 for compatible initial reduced states.

### Suggested fixes

- Withdraw the claim that nothing classical converges at this accuracy.
- Separate “certified by this bound” from “actually converged” and “impossible.”
- Derive or calibrate an error estimate with the correct short-time limit, initial
  state, observable support, boundary conditions, and prefactor.
- Distinguish the initial correlation length from the spatial decay parameter of
  a dynamical truncation bound.
- Give finite-size error its own share of the error ledger. Currently this
  calculation reuses `frac_trotter`, which is also spent on simulation error.
- Apply the same physical accuracy criterion to classical and quantum methods,
  without interpreting an upper error bound as a lower complexity bound.

### New tests

1. **Exact t = 0 test:** evaluate a dimer-link correlation on an isolated dimer and
   on larger compatible product-of-dimers states. The result must agree exactly.
2. **Short-time cluster test:** compare increasing clusters at t = 0, 0.01, 0.05,
   and 0.1 using exact, memory-bounded dynamics.
3. **Certification semantics test:** if an upper bound exceeds the tolerance but
   the measured error does not, report “not certified by this bound,” not failure.
4. **Budget test:** include both finite-size and Trotter error in the total ledger
   for a thermodynamic-limit calculation.

## 2. The default Trotter calibration is used beyond its evidence

**Locations:** [hubbard.py](fhcost/hubbard.py), `W_MEASURED`, `w_measured`,
`trotter_steps`; [curves.py](fhcost/curves.py), `evaluate_band`.

The hard-coded calibration is described as measured on 4–12 sites at times
0.25–2. Outside its time interval, the coefficient is clamped. The current
default model uses it at these operating points:

| Architecture, n = 10^8 | Reachable sites m | Final time sqrt(m)/2 |
|---|---:|---:|
| NISQ + PEC | 18.9062 | 2.1741 |
| STAR | 42.2790 | 3.2511 |
| Surface FT | 290.5537 | 8.5228 |
| Pinnacle | 433.9176 | 10.4153 |

Independence of system size at **fixed short time** does not establish
independence along t proportional to sqrt(m), where the explored region grows
with the system. The new m^(7/4) gate scaling follows from this extrapolation
assumption; it is not independently established by the calibration.

The experimental paper qualifies its locality argument by fixed time and warns
about increasing Trotter effects around t = 1.5–2. A conservative commutator
upper bound is not invalid merely because an observable-specific error is much
smaller than the bound.

No calibration programs or raw error tables were found in this folder. The
stored coefficients therefore cannot be independently regenerated here.

### Suggested fixes

- Preserve the calibration program, raw errors, Hamiltonian conventions, initial
  states, link selections, boundaries, fluxes, and Trotter ordering.
- State the calibration domain in both lattice size and time on plots and outputs.
- Separate interpolation from extrapolation and report when a point leaves the
  calibrated domain.
- Treat the bound and empirical extrapolation as two scenarios, not automatically
  as guaranteed lower/upper uncertainty edges or a confidence interval.
- Validate the error coefficient along the actual size–time trajectory.
- For multiproduct formulas, calibrate or bound the higher-order remainder;
  validating an asymptotic convergence order does not determine its coefficient
  or justify reusing the second-order coefficient.

### New tests

1. Regenerate the stored coefficients from raw calculations with explicit metadata.
2. Hold out at least one size and time from fitting and test prediction accuracy.
3. Test pairs of increasing size and time satisfying the adopted t = sqrt(m)/v
   convention, rather than only fixed-time size sweeps.
4. Compare error-versus-step-count curves with the assumed power law at the
   proposed operating point, not just at very large step count.
5. Verify that all out-of-domain operating points are identified in generated
   tables and figures.

Source: [Experimental paper, Appendix C.3.4](https://arxiv.org/html/2510.26300).

## 3. Python, the explorer, and the documents disagree

**Locations:** [explorer.html](explorer.html), [README.md](README.md),
[HEADLINES.md](HEADLINES.md), [METHODS.md](METHODS.md).

Direct evaluation of the explorer's calculation functions with its baseline
second-order settings gives:

| At n = 10^6 | Current Python | HTML explorer |
|---|---:|---:|
| NISQ + PEC | 16.7982 | 7.8801 |
| STAR | 29.0004 | 10.6111 |
| Surface FT | 31.6542 | 37.8560 |
| Pinnacle | 57.2890 | 95.5000 |
| Classical band | 24–26 | 30–30 |

The explorer retains the old PEC coefficient, extensive Trotter model, fixed
signal, old Pinnacle construction, and older statistical costing. This is not
rounding or numerical solver precision: the models are different.

README reports NISQ 27 and STAR 43 at n = 10^6, whereas HEADLINES reports 17 and
29. METHODS contains duplicated sections with incompatible descriptions and
numbers. A generated headline file does not make copied prose self-updating.

### Suggested fixes

- Use one executable model, or export a versioned parameter/result record for
  the explorer and documentation.
- If a JavaScript implementation is retained, give every physical parameter and
  costing rule an explicit mapping to Python.
- Mark scenario differences explicitly instead of presenting them as equivalent
  implementations.
- Remove superseded duplicate sections and generate headline text from a shared
  result record.
- Attach a model/configuration identifier to tables, figures, and explorer output.

### New tests

1. Numerical parity tests between Python and JavaScript for defaults and several
   nondefault configurations, including precision, p, time convention, and MPF order.
2. Compare intermediate outputs as well as maxima: steps, signal, gates, shot
   cost, distance, workspace, factories, physical footprint, and time per shot.
3. Check that every displayed headline matches the shared result record.
4. Include zero-capacity cases and transitions between factory/code choices.

## 4. Factory cleanup neglects circuit-level errors

**Location:** [ftqc.py](fhcost/ftqc.py), `FACTORIES` and `select_factory`.

The factory table assigns cultivation plus cleanup an output error of

```text
35 * (2e-9)^3 = 2.8e-25.
```

This accounts for ideal suppression of input-state errors but not faults
introduced by the cleanup circuit. The cleanup retains a fixed physical
footprint and cycle count. Cleaner input states do not eliminate faults in the
distillation circuitry itself.

`select_factory()` also does not use `cfg.p`. At an output target of 1e-10, the
review obtained identical factory specifications at physical error rates 1e-5,
1e-3, and 3e-3. Those cannot be interpreted as a calibrated hardware sensitivity.

### Suggested fixes

- Include both input-state and circuit-fault contributions to output infidelity.
- Select code distances, footprint, fidelity, rejection probability, and throughput
  jointly for the requested physical error and output target.
- Use a consistent failure-probability-to-observable-bias conversion. The surface
  logical-error check uses a factor of two, while its magic-error allocation does
  not; either justify the distinction or use compatible conservative conversions.
- Optimize complete machine throughput/cost, rather than choosing a factory only
  by its individual footprint before determining how many copies are needed.

### New tests

1. Reproduce published factory operating points at their stated physical errors.
2. Reduce input infidelity while keeping factory protection fixed: the output
   should approach a circuit-error floor, not zero without limit.
3. Sweep physical error and required output error; verify valid changes in
   protection and resource choice.
4. Check that rejection/stall costs and accepted-output rate are included.

Source: [Litinski, Magic State Distillation: Not as Costly as You Think](https://arxiv.org/pdf/1905.06903).

## 5. Pinnacle lacks the common FT resource/error ledger

**Location:** [ftqc.py](fhcost/ftqc.py), `pinnacle_point`.

Pinnacle uses the shared Hamming-weight-phased T-count, but does not include its
workspace in the logical-qubit count. It also does not test magic-state fidelity.

At the current n = 10^8 Pinnacle point:

| Quantity | Value |
|---|---:|
| T states per shot | 3.314e6 |
| Magic-error allowance | 0.00032 |
| N_T p_T using nominal p_T = 1e-9 | 0.003314 |

The nominal engine specification does not certify the requested allocation: the
union-bound contribution is over ten times the allowance, even before applying
a worst-case failure-to-bias factor. A cleaner engine could address this, but
must be selected and costed.

The 4,410-qubit engine also has a physical production time and approximately 10%
rejection rate. The model assumes one successful state every processing-code
cycle. At n = 10^5 it selects distance 16, hence an 18-round logical cycle, while
the cited engine requires at least 22 rounds for its distillation measurements
alone, before reaction-time constraints.

### Suggested fixes

- Apply the same data, synthesis, magic-state, and measurement error ledger to
  Pinnacle and surface FT.
- Add the workspace required by the chosen compiled algorithm.
- Select engine fidelity and protection for the workload instead of retaining a
  fixed 4,410-qubit engine at all operating points.
- Set the supply period from both processing and engine timing; include rejection.
- Charge any additional engine/lane parallelism with a feasible dependency schedule
  and corresponding footprint.
- Continue distinguishing Pinnacle's nonlocal connectivity assumption from the
  original nearest-neighbour physical grid.

### New tests

1. Verify the per-shot magic-state error allocation at every plotted point.
2. Check that T consumption never exceeds successful engine production.
3. Include a case where the processor is faster than the engine and must stall.
4. Compare architecture ledgers field-by-field for an identical compiled circuit.

Source: [Pinnacle architecture paper, Section V.2](https://arxiv.org/html/2602.11457v2).

## 6. Hamming-weight workspace does not match the cited circuit

**Location:** [ftqc.py](fhcost/ftqc.py), `hwp_workspace` and `t_counts`.

The workspace formula charges a logarithmic weight register plus a purported
phase-gradient register whose size comes from the rotation-synthesis T-count.
That is not the workspace formula for the cited construction.

For a batch of m rotations, Campbell's Appendix E gives a construction requiring
`m - popcount(m)` clean ancillas. Examples:

| Batch size | Model workspace | Cited construction's clean ancillas |
|---|---:|---:|
| 64 | 64 | 63 |
| 256 | 70 | 255 |
| 432 | 73 | 428 |

Agreement near one small batch size does not validate the scaling. A different
low-workspace implementation may be possible, but its T-count and depth must be
derived consistently. The synthesis allocation also divides by the number of
rotation groups rather than all the synthesized weight-register rotations.

### Suggested fixes

- Choose one explicit circuit construction, including batching and workspace reuse.
- Derive qubits, T-count, synthesis count, and depth from that same construction.
- Do not interpret a synthesis T-count as a register size without a specified
  phase-gradient implementation and its preparation cost.
- Allocate synthesis error across all synthesized operations, including the
  weight-register rotations and any multiproduct branches.

### New tests

1. Compile batch sizes 8, 16, 64, and 256 and count peak live ancillas.
2. Verify that the logical-depth estimate respects circuit dependencies.
3. Compare the hand-derived resource formula with compiled counts.
4. Check the accumulated synthesis-error allowance for the entire shot.

Source: [Campbell, Early fault-tolerant simulations of the Hubbard model, Appendix E](https://arxiv.org/pdf/2012.09238).

## 7. The classical baseline misses the U = 0 easy limit

**Locations:** [classical.py](fhcost/classical.py), `band`;
[make_figure.py](make_figure.py); [crossovers.md](crossovers.md).

The classical reference remains the same ED-only band at U = 0. However, the
selected low-weight density/spin correlations are classically tractable in that
limit, even for the non-Gaussian triplet initial state. The experimental paper
develops the applicable calculation in Appendix E. Difficulty of sampling the
full distribution is a different question from computing these observables.

For U = 4, the 24–26-site band is an estimated ED capacity. The code acknowledges
that tensor-network capability is unresolved, but the figure still labels the
region “classically easy” and visually emphasizes crossings as though this were
the decisive classical comparison.

### Suggested fixes

- Include specialized classical methods applicable at each parameter point.
- Distinguish low-weight observable estimation from full-state or sampling tasks.
- Label the interacting band explicitly as “estimated ED capacity.”
- Do not infer a classical impossibility boundary from a nonconverged TDVP run.
- Treat peak-FLOP scaling of a TDVP workload as an optimistic arithmetic estimate,
  not demonstrated exascale wall time; memory, communication, SVD/MPO costs, and
  achievable parallelism need separate accounting.

### New tests

1. Implement the appropriate U = 0 correlation calculation and validate it against
   exact small-system evolution.
2. Ensure the U = 0 comparison does not present the ED memory threshold as the
   applicable classical frontier for the chosen observable.
3. For U = 4, vary both TDVP time step and bond dimension against an exact small
   reference. A plateau in bond dimension alone does not identify its cause.
4. Record achieved observable error and measured wall time for classical benchmarks.

Source: [Experimental paper, Appendix E](https://arxiv.org/html/2510.26300).

## 8. The signal fit does not support per-time relative accuracy

**Locations:** [hubbard.py](fhcost/hubbard.py), `signal`, `signal_at`,
`eps_absolute`; [selftest.py](fhcost/selftest.py), stored `PUB` data.

The Gaussian envelope sets the absolute tolerance, but omits the oscillations
acknowledged in OPEN_ITEMS. At U = 4, t = 1.5:

| Quantity | Value |
|---|---:|
| Model signal | 0.064012 |
| Stored comparison datum | 0.0466 |
| Ratio of squared signal scales | 1.8869 |

Using the envelope permits a tolerance about 37% larger at this point; the
statistical shot cost alone differs by approximately 1.89 times. The stored
datum is not being asserted as exact ground truth, but the discrepancy shows
that the fit check does not establish the precision needed by the costing model.

Zero initial derivative does not imply Gaussian decay at later times. The test
labelled “the decay is Gaussian” checks the configured exponent and a comparison
with an exponential; it is not an independent validation of the functional form.

### Suggested fixes

- Distinguish an empirical envelope from a lower confidence bound on signal size.
- Evaluate tolerances and shot requirements at each requested time.
- Include oscillatory minima, zeros, and uncertainty in the inferred signal.
- Specify an absolute-error floor independently of a fitted late-time residual.
- Preserve fit code, data selection, covariance assumptions, and residuals;
  smoothed data points should not automatically be treated as independent.

### New tests

1. Held-out prediction tests for times not used in fitting.
2. Per-time coverage tests using known signals with oscillations and zeros.
3. Check that a monotone envelope does not silently loosen the tolerance near
   an actual signal minimum.
4. Test short-time derivatives as local constraints without treating them as proof
   of a global Gaussian law.

## 9. Zero uncertainty edges disappear from the plot

**Location:** [make_figure.py](make_figure.py), `msk` and the uncertainty
`fill_between` loop.

Zero capacities are replaced with NaN before plotting. For NISQ, the
extensive-bound edge is zero at the checked budgets, while the empirical edge
is positive. The shading therefore disappears instead of showing the possibility
that the smallest lattice cannot run.

At n = 10^3, 10^6, and 10^8 the checked NISQ band endpoints were respectively
`(0, 13.00)`, `(0, 16.80)`, and `(0, 18.91)`.

### Suggested fixes

- Represent “below minimum feasible lattice” explicitly.
- Clip the plotted edge to the axis floor only with an annotation explaining the
  clipping, or use hatching/separate feasibility markers.
- Distinguish scenario envelopes from statistical uncertainty in the legend.

### New tests

1. Plot synthetic bands with one zero edge and verify that their uncertainty
   remains visible.
2. Check rendered figures for all-infeasible and partially feasible scenarios.
3. Ensure every claimed band has an identifiable visual encoding.

## 10. Integer lattice reporting ignores initial-state constraints

**Locations:** [curves.py](fhcost/curves.py), `max_integer_L`;
[headlines.py](headlines.py).

The current headline reports a feasible 5×5 surface-FT lattice. But the stated
half-filled, zero-total-spin initial state contains one holon, one doublon, and
triplet pairs on all remaining sites. On 25 sites, 23 remain after placing the
defects, so that state cannot be formed. Equal up/down populations at half
filling also require an even total number of sites.

Taking `floor(sqrt(m))` is therefore insufficient. This is separate from whether
continuous m is a useful capacity proxy.

### Suggested fixes

- Define the admissible lattice family, boundaries, flux pattern, filling, and
  defect/dimer placement.
- Enumerate valid integer lattices and evaluate their resource feasibility directly.
- If odd-size states are desired, specify their changed filling/spin/defect
  configuration and recalibrate the corresponding workload.

### New tests

1. Construct the initial state specification for each reported integer size and
   verify particle number, spin sector, and complete dimer covering.
2. Reject incompatible sizes instead of rounding into them.
3. Check compiled resource feasibility for each admissible candidate, including
   code-block and replica packing.

## Additional implementation checks

- `error_ledger()` rejects over-allocation only when called explicitly. Resource
  entry points accept, for example, `frac_stat=0.9` with the other default shares
  unchanged. Validate configurations before costing.
- The weak/strong `s_sig` sensitivity rows are identical because the default
  `signal_regime="curve"` bypasses `s_sig`. Either switch to the fixed-signal
  scenario or vary the parameters actually used by the curve.
- The claim that the buffered cluster approach fails for every xi from 0.2 to 1.0
  disagrees with the implementation: `cluster_t_reach(xi=0.2)` returned about
  0.3902. This does not validate the cluster model, but it contradicts the prose.
- `max_m_fixed_t(0.02)` returns about 99,998 sites from an entropy-based MPS model,
  while the converged-classical function reports zero time reach. These are not
  interchangeable questions; outputs must identify the method and convergence
  assumptions supporting each result.
- Several tests assert narrative conclusions from the same formulas that produce
  those conclusions. Such tests can detect regressions but are not independent
  scientific validation.

## Recommended next work

1. **Fix exact-limit semantics first:** t = 0, U = 0, and admissible initial states.
2. **Make evidence reproducible:** include calibration code, raw data, parameter
   metadata, and held-out checks.
3. **Unify the implementations:** establish Python/JavaScript parity and a shared
   result record before publishing another figure.
4. **Build a common FT resource ledger:** explicit circuits, workspace, synthesis,
   engines/factories, error allocations, and feasible schedules.
5. **Validate precision at each time:** use actual signal variation and uncertainty,
   rather than relying on a smooth envelope.
6. **Recompute architecture comparisons:** only after the above, with extrapolated
   scenarios and unresolved classical capability clearly identified.

Start new dynamics tests with small, explicitly memory-bounded systems. The
repository's existing warning about large ED memory requirements still applies.
No large simulations, model edits, or figure regeneration were performed as
part of this review.

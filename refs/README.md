# References — what each number comes from

The PDFs in this directory are **gitignored** (they are other people's papers).
This file is committed, so the provenance survives even though the files do not.
Every constant in `fhcost/` that came from outside this repository is listed
here with the paper, the quantity, and where in the paper it is stated.

Re-download anything with `curl -sL -o refs/<id>.pdf https://arxiv.org/pdf/<id>`.

| file | paper |
|---|---|
| `2012.09238_campbell_hubbard.pdf` | Campbell, *Early fault-tolerant simulations of the Hubbard model* |
| `1905.06903_litinski_distillation.pdf` | Litinski, *Magic State Distillation: Not as Costly as You Think* |
| `2409.17595_cultivation.pdf` | Gidney, Shutty & Jones, *Magic state cultivation* |
| `2602.11457_pinnacle.pdf` | Webster *et al.* (Iceberg Quantum), the Pinnacle architecture |
| `2510.26300_phasecraft_quantinuum.pdf` | the Fermi–Hubbard dynamics experiment this model targets |
| `2511.05465.pdf` | *Helios: A 98-qubit trapped-ion quantum computer* |
| `2308.07915.pdf` | Bravyi *et al.*, *High-threshold and low-overhead fault-tolerant quantum memory* (BB codes) |
| `2304.05420.pdf` | Evered *et al.*, *High-fidelity parallel entangling gates on a neutral-atom quantum computer* |
| `1711.04789.pdf` | Kivlichan *et al.*, linear-depth fermionic swap networks |
| `2408.13687.pdf` | Google Quantum AI, *Quantum error correction below the surface code threshold* (Willow) |
| `quantinuum_helios_datasheet_v1.2.pdf` | Quantinuum Helios Product Data Sheet v1.2, 2026-06-09 |

---

## `fhcost/platform.py`

### Helios — arXiv:2511.05465, and the data sheet

| constant | value | where |
|---|---|---|
| `t_layer` | **55 ms** | §II: *"resulting in an average of 55 ms per layer"*. Their own **"depth-1 time"**, defined as the duration of a depth-10 random-pairing program divided by 10, and described as *"our characteristic figure of merit for processor speed"*. |
| `t_2q` | ~70 µs | §II: *"The 2Q gate operation itself requires approximately ∼70 µs to execute."* |
| `n_parallel_2q` | 4 | *"2Q beams are applied in the four 2Q operation zones"*; data sheet, *Two qubit quantum logic zones: 4*. |
| `p_2q` | 7.9(2)×10⁻⁴ | Abstract. **Better than the slide's 10⁻³.** |
| `n_demonstrated` | 98 | Title. |
| connectivity | all-to-all | Via ion transport through an X-junction, not wiring. |

**The number that matters.** The layer time is **not** set by the gate. The
paper's per-layer breakdown is Rotate 18.2 ms + Global Shift 7.9 ms + Junction
4.5 ms + Other Shifts 4.3 ms + Split/Combine 4.1 ms + Four-ion Shift 1.7 ms +
Static 0.3 ms ≈ **41 ms of transport**, against 70 µs of gate — transport
dominates by ~600×. This is why `Platform.layer_seconds` uses a *measured* layer
time where one exists instead of building one up from gate counts: no gate count
can predict ion sorting.

### Cultivation — arXiv:2409.17595, and its released stats

The paper's figure gives the expected volume; the **released simulation data**
(Zenodo 10.5281/zenodo.13777072, `stats.csv`, mirrored at
`calibration/data/cultivation_stats.csv`) gives the underlying counts. The
`end2end-inplace-distillation` row at `p = 1e-3, d1 = 5, d2 = 15` carries
`q = 463`, `r = 20`, 10¹² shots, and a complementary-gap histogram of 117
kept-count bins plus 113 error-count bins.

Reconstructing the error/discard trade from that histogram reproduces the
paper's headline: at gap cut 100 the measured error is **1.90×10⁻⁹** against
their quoted 2×10⁻⁹, with **73.0 attempts** per accepted state — a 98.6%
discard rate against their quoted 99%.

| constant | value | status |
|---|---|---|
| `CULT_FOOTPRINT` | 463 | exact (`q`) |
| `CULT_ROUNDS_PER_ATTEMPT` | 20 | exact (`r`) |
| `CULT_ATTEMPTS_1E3` | 73.0 | exact, reconstructed |
| `CULT_VOLUME_1E3` | 3×10⁴ | their integrated value, inside the bracket 9.3×10³–6.8×10⁵ |

### Neutral atom — arXiv:2304.05420

| constant | value | where |
|---|---|---|
| `t_2q` | 275 ns | Extended data, time-optimal single-pulse gate, `T = 275 ns`. |
| `p_2q` | 5×10⁻³ | Abstract: 99.5% fidelity. |
| `n_parallel_2q` | 60 | Abstract: *"on up to 60 atoms in parallel"*. |

`t_layer` is deliberately **left unset**. Atom rearrangement costs 10²–10⁴ µs per
move against a 275 ns gate, so a layer needing reconfiguration is transport-bound
exactly as Helios is — but no layer time for a *reconfiguring* circuit is
published, so the model declines to invent one.

### Superconducting — the slide spec, plus Willow (arXiv:2408.13687)

| constant | value | where |
|---|---|---|
| `t_2q`, `t_meas` | 10 ns, 100 ns | The slide's own hardware spec. Not measured here. |
| `t_round` | 1.1 µs | Willow abstract: *"a cycle time of 1.1 µs"*. **Measured**, and slightly slower than the 1 µs the model had assumed. |
| `n_demonstrated` | 105 | Willow. |

### `sc_long_range` — hypothetical, and labelled so

Bravyi *et al.* arXiv:2308.07915: BB-code Tanner graphs have **vertex degree
six** and *"consist of two edge-disjoint planar subgraphs"* (§1). So hosting a
BB/GB code needs **two coupler layers**, not all-to-all — a much weaker
requirement than the Pinnacle caveat implied. The gross code is `[[144,12,12]]`.

`n_demonstrated = 0` because no such device exists. This row is an upper bound on
what connectivity alone could buy, not a machine.

---

## The published simulation deposit — Zenodo 17799843

`calibration/tdvp_check.py` (O10) needs more than the summary numbers: it needs
every bond dimension at every time, and the exact reference to compare them to.
Both are in the data release that accompanies arXiv:2510.26300.

| item | value |
|---|---|
| DOI / record | Zenodo **17799843** |
| archive | `fermionic_dynamics.zip`, 112 MB |
| licence | **CC-BY-4.0** |
| extracted to | `refs/fermionic_dynamics/{U_0,U_4}/exp_vals.h5` (566 + 489 MB) |
| in git? | **no** — 1 GB; `refs/fermionic_dynamics/` is gitignored, this row is the record |

```
curl -sL -o fermionic_dynamics.zip https://zenodo.org/records/17799843/files/fermionic_dynamics.zip
unzip -q fermionic_dynamics.zip -d refs/fermionic_dynamics
pip install tables uncertainties        # pandas-in-HDF5 with pickled error bars
python3 calibration/tdvp_check.py --json
```

`exp_vals.h5` holds pandas frames keyed by method: raw hardware, TFLO/GPR
mitigated, **TDVP at chi = 256, 512, 1024, 2048**, Majorana propagation, and the
**exact FLO** reference at U = 0. `obs_type = spin_correlator_neighbours` is the
observable this model costs.

Two traps worth recording, since both silently return the wrong thing:

- the bond dimension lives in the `method` string (`TDVP chi=512`), **not** in
  `max_bond_dimension`, for this observable;
- the `spins` column is `None` throughout, which makes any pandas merge that
  includes it drop every row without error.

**And the tooling answer.** The deposit's own metadata records the simulator as
`Tenpy $\chi=512$`: their TDVP is **TeNPy**. That settles what the missing dt
scan should be written in — `calibration/tdvp_dt.py` configures the same library
rather than hand-rolling a fermionic 2D TDVP on a doubly periodic torus. TeNPy
1.1.0 and quimb 1.15.0 are installed locally; neither is on lenore yet.

---

## Constants sourced earlier (unchanged)

| constant | source |
|---|---|
| `W_MEASURED`, α | measured here — `calibration/trotter_cal.py`, not from a paper |
| `w_commutator` | Campbell arXiv:2012.09238, PLAQ bound |
| `hwp_group`, `alpha = b - w(b)` | Campbell arXiv:2012.09238 App. E, Thm 2 |
| `MAGIC_SOURCES` (13 rows) | Litinski arXiv:1905.06903 Table 1 |
| cultivation rows | Gidney *et al.* arXiv:2409.17595 Figs 1–2 |
| `PIN_ENGINE_TABLE`, `engine_cycles` | Webster *et al.* arXiv:2602.11457 Eqs 7–11 |
| `SIGNAL_DATA` | Zenodo 17799843, accompanying arXiv:2510.26300 |
| free-fermion collapse | arXiv:2510.26300 App. E; reimplemented and validated in `calibration/free_fermion.py` |
| `step_depth` JW swap term | Kivlichan *et al.* arXiv:1711.04789 |
| `p_logical` fowler/willow | Fowler *et al.* PRA 86 032324; Google arXiv:2408.13687 |

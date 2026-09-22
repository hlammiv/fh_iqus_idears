# Q1 · 2D Fermi–Hubbard resource estimate

Quantitative replacement for the whiteboard sketch on slide 2 of
`Q1_Fermi_Hubbard-1.pdf`. Self-contained: nothing outside this folder is
imported, written, or executed.

| | |
|---|---|
| `fhcost/` | the cost models (`hubbard` → `nisq` / `ftqc` / `classical` → `curves`) |
| `make_figure.py` | 3-panel paper figure → `figures/fh_resource_estimate.pdf` |
| `make_slide.py` | 1-panel slide version → `figures/*_slide.png` |
| `crossover_table.py` | sensitivity sweep → `crossovers.md` |
| `check_palette.py` | CVD validation of the figure palette (exit 1 on fail) |
| `explorer.html` | interactive version, published as an Artifact |
| `METHODS.md` | derivations, constants, citations, open parameters |
| `vendor/` | copies of three read-only helpers from elsewhere in `~/Desktop/QC` |

```
python3 -m fhcost.selftest    # 40 checks, pure arithmetic, <1 s, ~27 MB peak
python3 make_figure.py && python3 make_slide.py && python3 crossover_table.py
python3 check_palette.py
```

Nothing here runs a simulation: the classical entries are closed-form memory and
flop counts, not state vectors or tensors. If the models are ever validated
against a real small-lattice ED run, note that m = 16 already needs ~2.6 GB per
state vector and m = 18 needs ~38 GB — run that on `lenore_remote` (port 60022),
not on this machine, and put a hard memory guard on it.

## Headline

- Mitigated NISQ saturates: **m = 6.7 → 8.5 from n = 10³ to 10⁸.** Ten decades of
  qubits buy ~1.3× in lattice size, and by Takagi *et al.* this holds for *any*
  mitigation strategy once the run has to fit in a week.
- The unmitigated p = 10⁻³ device never reaches m = 4 at any n.
- **Classical band m = 24–56.** NISQ+PEC and STAR both sit below it at every n,
  so the slide's "logical advantage = FT overtakes NISQ" marks the wrong crossing.
- Surface-code FT clears the band at **n ≈ 2×10⁶** (idealised p_L) or later
  (measured p_L). STAR overtakes NISQ at n ≈ 5×10⁴ and settles ~1.4× above it.
- The only knob that improves STAR *relative to* NISQ is the rotation injection
  error. Lowering p lifts both curves equally; footprint and clock are
  logarithmic levers. See `METHODS.md` §4.
- **Pinnacle (QLDPC, arXiv:2602.11457)** beats the surface code by 2.5–3.7× on
  this workload and clears the classical band **8× earlier in n** (2.9e5 vs 2.3e6).
- Trotter-step (multiproduct) extrapolation is the largest single lever — and it
  has an optimum, because the shot cost ‖c‖₁² eventually overtakes the depth
  saving for the shot-limited FT curve.

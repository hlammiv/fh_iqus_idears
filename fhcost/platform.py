"""Hardware platforms: connectivity, gate parallelism and clock speed.

WHY THIS EXISTS
  Every arm in this model was costed on ONE implicit machine: a nearest-neighbour
  2D superconducting grid, 10 ns two-qubit gates, and unlimited gate parallelism.
  That assumption was never a field, only prose in budget.py's docstring, and it
  is doing more work than any other input -- three of four operating points are
  clock-limited, not qubit-limited.

  It is also already violated. Generalised bicycle codes need non-local
  connectivity that a nearest-neighbour grid does not provide; the model flagged
  this in a docstring and costed the arm anyway. `pin_nonlocal` existed as a
  Config field and NOTHING read it.

THE SAME DISCIPLINE AS THE MAGIC-STATE TABLES
  Published operating points, sourced line by line, with a refusal outside the
  tabulated range -- as in ftqc.MAGIC_SOURCES and ftqc.PIN_ENGINE_TABLE. Nothing
  here is fitted or invented. Every number below is quoted in refs/README.md
  with its paper, and the PDFs are in refs/.

THE ONE THING THAT SURPRISED US
  For trapped ions the layer time is NOT set by the gate. Helios needs ~70 us for
  a two-qubit gate but ~55 ms for a circuit layer, because ion transport
  dominates by a factor of ~600: Rotate 18.2 + Global Shift 7.9 + Other Shifts
  4.3 + Junction 4.5 + Split/Combine 4.1 + Four-ion Shift 1.7 + Static 0.3 ms.
  So a platform that quotes a measured layer time gets used directly, and only
  platforms without one are built up from gate counts.
"""
from __future__ import annotations
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Platform:
    """A machine, as published. See refs/README.md for every number's source."""
    name: str
    connectivity: str        # "nn_grid" | "thickness_2" | "all_to_all"
    t_2q: float              # seconds, one two-qubit gate
    t_meas: float            # seconds, readout
    t_round: float           # seconds, one syndrome-extraction cycle
    n_parallel_2q: float     # simultaneous two-qubit gates; inf = unlimited
    t_layer: float | None    # MEASURED seconds per circuit layer, if published
    p_2q: float              # two-qubit infidelity
    n_demonstrated: int      # largest physical qubit count actually built
    # Can this connectivity run TRANSVERSAL logical gates? If so the machine
    # needs O(1) syndrome-extraction rounds per logical layer instead of O(d),
    # which is the whole reason a slow-clock machine is not simply hopeless.
    transversal: bool
    source: str

    def layer_seconds(self, gates_in_layer: float = 1.0) -> float:
        """Wall-clock for one layer of `gates_in_layer` two-qubit gates.

        A measured layer time wins outright -- it already contains the transport,
        cooling and sorting that dominate on a mobile-qubit machine and that no
        gate count can predict. Otherwise the layer is the slower of its gate
        duration and its throughput limit.
        """
        if self.t_layer is not None:
            return self.t_layer
        serial = gates_in_layer / self.n_parallel_2q * self.t_2q
        return max(self.t_2q, serial)


# --- the table. Every row is a citation; see refs/README.md -----------------
SUPERCONDUCTING = Platform(
    name="superconducting grid",
    connectivity="nn_grid",
    t_2q=10e-9,              # the slide's own spec
    t_meas=100e-9,           # the slide's own spec
    t_round=1.1e-6,          # MEASURED on Willow, arXiv:2408.13687
    n_parallel_2q=math.inf,  # a planar grid gates every disjoint pair at once
    t_layer=None,
    p_2q=1e-3,               # the slide's own spec
    n_demonstrated=105,     # Willow
    transversal=False,      # planar grid: lattice surgery, O(d) rounds
    source="slide spec for the clock; Google Quantum AI, arXiv:2408.13687 "
           "(1.1 us cycle, 105 qubits) for the syndrome round",
)

HELIOS = Platform(
    name="Quantinuum Helios",
    connectivity="all_to_all",
    t_2q=70e-6,              # "the 2Q gate operation itself requires ~70 us"
    t_meas=1e-3,             # ground-state cooling ~3 ms dominates; conservative
    t_round=55e-3,           # no surface-code demo; a round is at least a layer
    n_parallel_2q=4.0,       # four two-qubit operation zones
    t_layer=55e-3,           # MEASURED "average of 55 ms per layer" (depth-1 time)
    p_2q=7.9e-4,             # 7.9(2)e-4, better than the slide's 1e-3
    n_demonstrated=98,
    transversal=True,       # transversal gates demonstrated on trapped ions
    source="Helios: A 98-qubit trapped-ion quantum computer, arXiv:2511.05465; "
           "55 ms/layer is their own 'depth-1 time' figure of merit",
)

NEUTRAL_ATOM = Platform(
    name="neutral atom (reconfigurable)",
    connectivity="all_to_all",   # by moving atoms, not by wiring
    t_2q=275e-9,             # time-optimal single-pulse gate, T = 275 ns
    t_meas=1e-3,
    t_round=1e-3,            # rearrangement-dominated; see the caveat below
    n_parallel_2q=60.0,      # "99.5% fidelity on up to 60 atoms in parallel"
    # Rearrangement is 100 us - 10 ms per move against a 275 ns gate, so a layer
    # that needs reconfiguration is transport-bound exactly as Helios is. There
    # is no published layer time for a reconfiguring circuit, so this platform
    # is costed WITHOUT one and the model says so rather than inventing it.
    t_layer=None,
    p_2q=5e-3,               # 99.5%
    n_demonstrated=60,
    transversal=True,       # Zhou et al. arXiv:2406.17653; atom arrays are
                            # the platform that motivated the scheme
    source="Evered et al., arXiv:2304.05420 (275 ns gate, 99.5%, 60 in parallel)",
)

SC_LONG_RANGE = Platform(
    name="superconducting, two coupler layers",
    connectivity="thickness_2",
    t_2q=10e-9,
    t_meas=100e-9,
    t_round=1.1e-6,
    n_parallel_2q=math.inf,
    t_layer=None,
    p_2q=1e-3,
    n_demonstrated=0,       # NOT BUILT. A hypothetical, and labelled as one.
    transversal=False,      # two coupler layers is still not the flexible
                            # connectivity transversal gates need
    source="Bravyi et al., arXiv:2308.07915: BB-code Tanner graphs have vertex "
           "degree six and decompose into TWO edge-disjoint planar subgraphs, so "
           "two coupler layers suffice -- this is far weaker than all-to-all. "
           "No such device exists; n_demonstrated = 0 marks it as hypothetical.",
)

PLATFORMS = {
    "superconducting": SUPERCONDUCTING,
    "helios": HELIOS,
    "neutral_atom": NEUTRAL_ATOM,
    "sc_long_range": SC_LONG_RANGE,
}

# Which codes each connectivity class can host. This is what pin_nonlocal was
# documented to mean and never did.
#   nn_grid      a planar degree-4 lattice: the surface code, and STAR on top of
#                it. Bivariate/generalised bicycle codes are NOT embeddable.
#   thickness_2  degree-6, two edge-disjoint planar subgraphs (Bravyi et al.).
#                Enough for BB/GB codes; still not all-to-all.
#   all_to_all   anything.
ADMITS = {
    "nn_grid":     {"surface", "star"},
    "thickness_2": {"surface", "star", "gb"},
    "all_to_all":  {"surface", "star", "gb"},
}


def get(name: str) -> Platform:
    if name not in PLATFORMS:
        raise ValueError(f"unknown platform {name!r}; have {sorted(PLATFORMS)}")
    return PLATFORMS[name]


def admits(code: str, platform_name: str) -> bool:
    """Can this platform's connectivity host this code family?"""
    return code in ADMITS[get(platform_name).connectivity]


# Syndrome-extraction rounds per logical layer.
#
#   LATTICE SURGERY needs O(d) rounds, because a measurement error has to be
#   caught by repetition. That is the cost model every arm in this repository
#   used, and it is the right one for a planar grid.
#
#   TRANSVERSAL gates need O(1). Zhou et al. (arXiv:2406.17653, Nature 2025)
#   prove that with transversal operations and correlated decoding the deviation
#   from the ideal logical measurement distribution is exponentially small in d
#   with only a CONSTANT number of rounds -- "transversal algorithmic fault
#   tolerance". It requires the flexible connectivity that mobile-qubit machines
#   have and a planar grid does not.
#
# This matters enormously here. Charging a trapped-ion or atom machine O(d)
# rounds at its slow clock is the calculation that gives "several years" for
# RSA-2048 in the literature; the transversal architecture gives 5.6 days for
# the same problem at the same 1 ms cycle (Chen et al., arXiv:2505.15907),
# "close to 50x speed-up ... with no increase in space footprint". Costing these
# platforms with lattice surgery is costing them at something nobody proposes.
SE_ROUNDS_TRANSVERSAL = 1.0      # Theta(1); they insert one round after the gate


def se_rounds(d: int, platform_name: str) -> float:
    """Syndrome-extraction rounds per logical layer on this machine."""
    return SE_ROUNDS_TRANSVERSAL if get(platform_name).transversal else float(d)


def swap_network_free(platform_name: str) -> bool:
    """All-to-all removes the Kivlichan swap network and the routing penalty.

    This is where the connectivity advantage is actually paid back, and it has to
    be computed against the clock cost rather than asserted.
    """
    return get(platform_name).connectivity == "all_to_all"

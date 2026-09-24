"""Their compiled circuit, in closed form -- the external check on c_g and c_rot.

O1 left one thing hanging: STAR pays its PEC one-norm only over ROTATIONS while
NISQ pays it over every two-qubit gate in the cone, so the margin between the two
arms "is set almost entirely by the assumed ratio c_rot/c_g = 5/15. That ratio
has never been checked against a compiled circuit."

It can be. The experimental paper gives the swap-network gate count in closed
form (arXiv:2510.26300 App. C 5, Eqs. C20-C27) and tabulates it for four
lattices, so the count is not an estimate and not ours.

    n_h = 2 Lx Ly                                       horizontal hops
    n_v = 2 Lx Ly            (Ly even)                  vertical hops
        = 2 Lx (Ly - 1)      (Ly odd)                   -- the wrap is missed
    n_f = (2 Lx^2 - 5 Lx + 4) Ly   (Lx even)            FSWAPs
        = (2 Lx^2 - 4 Lx + 2) Ly   (Lx odd)
    C_nw = n_f + 2 n_v + n_h                            one forward pass
    C_b  = 0                 (Ly even)                  boundary-hop circuit
         = Lx^2 + 3 Lx       (Ly odd)
    C_c  = Lx Ly                                        Coulomb
    C_merge = 4 Ly           (Ly even)
            = 4 (Ly - 1)     (Ly odd)
    C_tot = N (2 C_nw + 2 C_b + C_c - C_merge) + C_merge

REPRODUCTION, and where it stops. This reproduces their Table I exactly for the
two (Lx even, Ly odd) rows -- 4x5 -> 428N + 16 and 4x7 -> 588N + 24, the latter
being the lattice they ran. It comes out 6 low on 4x6 and 5 low on 5x5, i.e. by
Ly and by Lx in the two cases the text's piecewise C_nw covers differently.
Rather than reverse-engineer the difference, `two_qubit_gates` refuses outside
the validated case. The instance this model cares about is inside it.
"""
from __future__ import annotations

# (Lx, Ly): (per-step coefficient, constant) from their Table I
TABLE_I = {(4, 5): (428, 16), (4, 6): (486, 24), (5, 5): (674, 16),
           (4, 7): (588, 24)}
VALIDATED = ((4, 5), (4, 7))          # the rows the formulas reproduce exactly

# Their own instance: Lx x Ly = 4 x 7 in the appendix's convention (28 sites,
# 56 qubits), k = 4 second-order steps, plus 39 two-qubit gates of state
# preparation -> "(588 x 4) + 24 + 39 = 2415 two-qubit gates overall".
THEIRS = {"Lx": 4, "Ly": 7, "steps": 4, "prep_2q": 39, "total_2q": 2415,
          "one_qubit_gates": 4627}


def _parts(Lx: int, Ly: int) -> dict:
    n_h = 2 * Lx * Ly
    n_v = 2 * Lx * Ly if Ly % 2 == 0 else 2 * Lx * (Ly - 1)
    n_f = ((2 * Lx ** 2 - 5 * Lx + 4) * Ly if Lx % 2 == 0
           else (2 * Lx ** 2 - 4 * Lx + 2) * Ly)
    return {"n_h": n_h, "n_v": n_v, "n_f": n_f,
            "C_nw": n_f + 2 * n_v + n_h,
            "C_b": 0 if Ly % 2 == 0 else Lx ** 2 + 3 * Lx,
            "C_c": Lx * Ly,
            "C_merge": 4 * Ly if Ly % 2 == 0 else 4 * (Ly - 1),
            # hops the swap network misses, both spin sectors
            "n_boundary": 0 if Ly % 2 == 0 else 2 * Lx}


def two_qubit_gates(Lx: int, Ly: int, steps: int) -> int:
    """Eq. (C27). Refuses outside the cases that reproduce their Table I."""
    if (Lx, Ly) not in VALIDATED and (Lx % 2, Ly % 2) != (0, 1):
        raise ValueError(
            f"{Lx}x{Ly}: the closed form reproduces their Table I only for "
            f"Lx even and Ly odd; it is {TABLE_I.get((Lx, Ly), ('?',))[0]} vs "
            f"ours on the other rows. Refusing rather than guessing.")
    p = _parts(Lx, Ly)
    per = 2 * p["C_nw"] + 2 * p["C_b"] + p["C_c"] - p["C_merge"]
    return steps * per + p["C_merge"]


def rotations(Lx: int, Ly: int, steps: int) -> int:
    """Arbitrary-angle rotations: one per Hamiltonian term per application.

    A second-order step applies every hopping term TWICE (forward then reverse)
    and the on-site term once. This is the count STAR pays its one-norm over,
    its Cliffords being error-corrected, and it is a LOWER bound on the
    one-qubit gate count, which also carries basis changes.

    IT IS ALSO ENCODING-INDEPENDENT, which is what makes it safe to adopt as
    c_rot. On a doubly periodic square lattice each site owns 2 bonds, so there
    are 2m bonds x 2 spins = 4m hopping terms, applied twice by a second-order
    step, plus m on-site terms:

        rotations per step = 8m + m = 9m   exactly

    and that is what their circuit gives, 252 on 28 sites. No compilation can
    use fewer arbitrary-angle rotations than there are term exponentials, so
    c_rot = 9 is a floor, not a fit -- and the c_rot = 5 it replaced was BELOW
    the number of terms in the Hamiltonian.
    """
    p = _parts(Lx, Ly)
    per = 2 * (p["n_h"] + p["n_v"] + p["n_boundary"]) + p["C_c"]
    return steps * per


def per_site_per_step(Lx: int, Ly: int) -> dict:
    """c_g and c_rot as this model defines them, measured on their circuit."""
    m = Lx * Ly
    g = two_qubit_gates(Lx, Ly, 2) - two_qubit_gates(Lx, Ly, 1)   # the marginal step
    r = rotations(Lx, Ly, 1)
    return {"m": m, "c_g": g / m, "c_rot": r / m, "ratio": r / g,
            "gates_per_step": g, "rotations_per_step": r}


def check_table() -> list[tuple]:
    """(lattice, ours, theirs, agrees) for every row of their Table I."""
    out = []
    for (Lx, Ly), (coef, const) in sorted(TABLE_I.items()):
        p = _parts(Lx, Ly)
        ours = (2 * p["C_nw"] + 2 * p["C_b"] + p["C_c"] - p["C_merge"],
                p["C_merge"])
        out.append(((Lx, Ly), ours, (coef, const), ours == (coef, const)))
    return out

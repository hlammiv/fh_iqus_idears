#!/usr/bin/env python3
"""U = 0: C^zz on the triplet-covering state, in POLYNOMIAL time -- and checked.

Second-pass review #7. The classical baseline was exact diagonalisation at every
parameter point, but at U = 0 the Hamiltonian is quadratic and the observable is
weight-4, so this correlator is classically cheap even though the initial state
is NOT Gaussian. The experimental paper develops the argument in its Appendix E
(arXiv:2510.26300); this is an independent implementation, validated against
exact many-body evolution.

THE STATE.  |psi> = prod_T (|A_T> - |B_T>)/sqrt(2) x (lonely sites), where each
triplet T spans two sites (i, j) with A_T = {i up, j down} and B_T = {i down,
j up}. Expanding gives 2^{N_T} Fock branches -- exponentially many, which is why
the paper calls it non-Gaussian.

WHY IT COLLAPSES.  Under free evolution n_i(t) = sum_ab R*_ia R_ib c†_a c_b with
R = exp(-i h t), so <n_i n_j>(t) is a sum of FOUR-fermion expectations in the
initial state. A four-fermion operator changes occupations on at most four
modes, and two branches differ by whole triplets (four modes each). So only two
kinds of term survive:

  DIAGONAL   b = b'. Signs cancel and the branch sum becomes an average over an
             independent per-triplet coin flip, so it is fixed by the one- and
             two-mode occupation statistics nu_a and nu2_ab. Wick's theorem then
             applies branch by branch. Cost O(M^2) in the mode count M = 2m.

  COHERENT   b, b' differing in exactly ONE triplet, always with sign -1. All
             four operators must lie inside that triplet, so with the modes of
             each triplet ordered CONTIGUOUSLY the Jordan-Wigner strings from
             outside cancel (the operator has an even number of factors) and
             the matrix element is a local 4-mode object -- computed here on a
             16-dimensional Fock space rather than by hand, because that is
             where sign errors live. Cost O(m).

Neither piece needs the 2^{N_T} branches. Total cost is one sparse matrix
exponential applied to four basis vectors, then O(M^2) per correlator.

    python3 calibration/free_fermion.py            # validate + time
    python3 calibration/free_fermion.py --sizes 4 6 8
"""
from __future__ import annotations
import argparse, itertools, json, math, pathlib, time
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

# --------------------------------------------------------------------------
# geometry: an Lx x Ly lattice, dimer covering along x, two lonely sites


def lattice(Lx: int, Ly: int, periodic: bool = True) -> list[tuple[int, int]]:
    """Nearest-neighbour bonds between site indices s = x + Lx*y."""
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            s = x + Lx * y
            if periodic or x + 1 < Lx:
                bonds.append((s, (x + 1) % Lx + Lx * y))
            if periodic or y + 1 < Ly:
                bonds.append((s, x + Lx * ((y + 1) % Ly)))
    # a 2-site direction with periodic wrap doubles the bond; keep it once
    return sorted({(min(a, b), max(a, b)) for a, b in bonds if a != b})


def covering(Lx: int, Ly: int) -> tuple[list[tuple[int, int]], int, int]:
    """(dimer list, doublon site, holon site). Dimers along x, last pair lonely."""
    sites = list(range(Lx * Ly))
    doublon, holon = sites[-2], sites[-1]
    rest = sites[:-2]
    dimers = [(rest[i], rest[i + 1]) for i in range(0, len(rest) - 1, 2)]
    spare = rest[2 * len(dimers):]          # an odd lattice leaves one site over
    return dimers, doublon, holon, spare


def mode_order(Lx, Ly):
    """Modes ordered so each triplet's four modes are CONTIGUOUS.

    This is the whole reason the coherent term is local: with the four modes of
    a triplet adjacent in the ordering, an even-parity operator supported on
    them carries no Jordan-Wigner string from the rest of the lattice.
    """
    dimers, doublon, holon, spare = covering(Lx, Ly)
    order = []
    for (a, b) in dimers:
        order += [(a, 0), (a, 1), (b, 0), (b, 1)]     # 0 = up, 1 = down
    for s in [doublon, holon] + list(spare):
        order += [(s, 0), (s, 1)]
    idx = {sp: k for k, sp in enumerate(order)}
    return order, idx, dimers, doublon, holon, list(spare)


def hopping_sparse(Lx, Ly, idx, J: float = 1.0, periodic: bool = True):
    """Always build this one. A dense 2m x 2m is 128 GiB at m = 65536."""
    M = 2 * Lx * Ly
    r, c, v = [], [], []
    for (a, b) in lattice(Lx, Ly, periodic):
        for s in (0, 1):
            i, j = idx[(a, s)], idx[(b, s)]
            r += [i, j]; c += [j, i]; v += [-J, -J]
    return csr_matrix((v, (r, c)), shape=(M, M))


def hopping(Lx, Ly, idx, J: float = 1.0, periodic: bool = True) -> np.ndarray:
    """Dense form, for the small exact reference only."""
    M = 2 * Lx * Ly
    if M > 8192:
        raise MemoryError(f"refusing a dense {M}x{M}; use hopping_sparse")
    return hopping_sparse(Lx, Ly, idx, J, periodic).toarray()


# --------------------------------------------------------------------------
# initial-state occupation statistics (closed form, no branch enumeration)


def nu_only(idx, dimers, doublon, holon, spare=()):
    """Just the one-mode occupations. O(m) memory.

    occupations() below also returns the full M x M nu2, which is 8.6 GB at
    m = 16384 -- fine for the small exact cross-checks, fatal at scale. The O(m)
    path must never call it, and occupations() now refuses rather than swapping.
    """
    M = 2 * (len(dimers) * 2 + 2 + len(spare))
    nu = np.zeros(M)
    for (a, b) in dimers:
        for m_ in (idx[(a, 0)], idx[(a, 1)], idx[(b, 0)], idx[(b, 1)]):
            nu[m_] = 0.5
    nu[idx[(doublon, 0)]] = 1.0
    nu[idx[(doublon, 1)]] = 1.0
    return nu


def occupations(idx, dimers, doublon, holon, spare=(), _max_modes=8192):
    """(nu_a, nu2_ab): one- and two-mode occupation statistics of the branches.

    Within a triplet the two branches are {i up, j down} and {i down, j up}, each
    with probability 1/2, so nu = 1/2 on all four modes and nu2 is 1/2 on the two
    co-occupied pairs and 0 on the rest. Different factors are independent.
    """
    M = 2 * (len(dimers) * 2 + 2 + len(spare))
    if M > _max_modes:
        raise MemoryError(f"refusing a dense {M}x{M} nu2 "
                          f"({M * M * 8 / 2**30:.1f} GiB); use nu_only")
    nu = np.zeros(M)
    nu2 = np.zeros((M, M))
    groups = []
    for (a, b) in dimers:
        A = (idx[(a, 0)], idx[(b, 1)])                 # branch A occupations
        B = (idx[(a, 1)], idx[(b, 0)])                 # branch B occupations
        groups.append(([A, B], 0.5))
    groups.append(([(idx[(doublon, 0)], idx[(doublon, 1)])], 1.0))
    groups.append(([()], 1.0))                          # holon: nothing occupied
    for s in spare:                                     # leftover sites: empty
        groups.append(([()], 1.0))
    mode_of_group = []
    for confs, w in groups:
        ms = sorted({m for c in confs for m in c})
        mode_of_group.append(ms)
        for c in confs:
            for m in c:
                nu[m] += w
            for m1 in c:
                for m2 in c:
                    nu2[m1, m2] += w
    # across independent factors nu2 factorises
    g_of_mode = {}
    for gi, ms in enumerate(mode_of_group):
        for m in ms:
            g_of_mode[m] = gi
    for a in range(M):
        for b in range(M):
            if g_of_mode.get(a, -1) != g_of_mode.get(b, -2):
                nu2[a, b] = nu[a] * nu[b]
    return nu, nu2


# --------------------------------------------------------------------------
# a tiny exact Fock engine, used only on FOUR modes for the coherent term


def _sign(mask: int, a: int) -> int:
    return -1 if bin(mask & ((1 << a) - 1)).count("1") % 2 else 1


def _apply(op: str, a: int, mask: int):
    """(new mask, sign) for c_a or c^dag_a on a bitmask; None if annihilated."""
    occ = (mask >> a) & 1
    if (op == "c" and not occ) or (op == "d" and occ):
        return None
    return mask ^ (1 << a), _sign(mask, a)


def four_mode_element(bra: int, ket: int, ops) -> complex:
    """<bra| c^dag_A c_B c^dag_G c_H |ket> on a 4-mode Fock space, bitmask basis.

    Computed, not hand-derived: this is exactly where fermionic sign errors hide.
    """
    A, B, G, H = ops
    st, s = ket, 1
    for op, a in (("c", H), ("d", G), ("c", B), ("d", A)):
        r = _apply(op, a, st)
        if r is None:
            return 0.0
        st, sg = r
        s *= sg
    return float(s) if st == bra else 0.0


# --------------------------------------------------------------------------
# the polynomial estimator


def czz_free(t: float, h: np.ndarray, nu, nu2, idx, dimers, site_i, site_j):
    """C^zz_{ij}(t) = <(n_iu - n_id)(n_ju - n_jd)> - <n_iu - n_id><n_ju - n_jd>."""
    R = expm(-1j * h * t)
    iu, idn = idx[(site_i, 0)], idx[(site_i, 1)]
    ju, jdn = idx[(site_j, 0)], idx[(site_j, 1)]

    def n1(p):                       # <n_p(t)>
        return float(np.real((np.abs(R[p]) ** 2) @ nu))

    tri_modes = [tuple(idx[(a, s)] for a in (x, y) for s in (0, 1))
                 for (x, y) in dimers]
    # relabel so each triplet's modes are 0..3 locally; A = {a up, b down}
    local = []
    for (x, y) in dimers:
        au, ad, bu, bd = (idx[(x, 0)], idx[(x, 1)], idx[(y, 0)], idx[(y, 1)])
        g = [au, ad, bu, bd]
        loc = {m: k for k, m in enumerate(g)}
        maskA = (1 << loc[au]) | (1 << loc[bd])
        maskB = (1 << loc[ad]) | (1 << loc[bu])
        local.append((g, loc, maskA, maskB))

    def n2(p, q):                    # <n_p(t) n_q(t)>
        rp, rq = R[p], R[q]
        ap, aq = np.abs(rp) ** 2, np.abs(rq) ** 2
        d1 = ap @ nu2 @ aq                                   # <n_a><n_g> term
        Mx = np.outer(np.conj(rp), rp) * np.outer(rq, np.conj(rq))
        d2 = np.sum(Mx * (nu[:, None] - nu2))                # exchange term
        x = 0.0 + 0j
        for (g, loc, mA, mB) in local:
            rp_g, rq_g = rp[g], rq[g]
            for a, b, c_, e in itertools.product(range(4), repeat=4):
                K = np.conj(rp_g[a]) * rp_g[b] * np.conj(rq_g[c_]) * rq_g[e]
                if K == 0:
                    continue
                el = (four_mode_element(mA, mB, (a, b, c_, e))
                      + four_mode_element(mB, mA, (a, b, c_, e)))
                if el:
                    x -= 0.5 * K * el          # the (-1)^{|b|+|b'|} = -1 sign
        return float(np.real(d1 + d2 + x))

    sz_i = n1(iu) - n1(idn)
    sz_j = n1(ju) - n1(jdn)
    nn = n2(iu, ju) - n2(iu, jdn) - n2(idn, ju) + n2(idn, jdn)
    return nn - sz_i * sz_j


def blocks_of(idx, dimers, doublon, holon, spare=(), nu=None):
    """[(modes, nu2_block - outer(nu, nu))] -- the correction to independence.

    nu2 factorises across the independent factors of the initial product state,
    so the only departure from nu_a nu_b lives inside a triplet (4 modes) or a
    lonely site (2 modes). That is what makes the O(m) form possible.
    """
    out = []
    for (a, b) in dimers:
        g = [idx[(a, 0)], idx[(a, 1)], idx[(b, 0)], idx[(b, 1)]]
        n2 = np.zeros((4, 4))
        for conf in ((0, 3), (1, 2)):                 # branch A, branch B
            for m1 in conf:
                for m2 in conf:
                    n2[m1, m2] += 0.5
        nn = np.array([nu[m] for m in g])
        out.append((np.array(g), n2 - np.outer(nn, nn)))
    for site, filled in ((doublon, True), (holon, False)):
        g = [idx[(site, 0)], idx[(site, 1)]]
        n2 = np.ones((2, 2)) if filled else np.zeros((2, 2))
        nn = np.array([nu[m] for m in g])
        out.append((np.array(g), n2 - np.outer(nn, nn)))
    return out


def czz_free_sparse(t, hsp, nu, blocks, idx, dimers, site_i, site_j):
    """Same estimator in O(m) time and memory -- what a frontier claim needs.

    Two changes from czz_free. The four rows of R = exp(-i h t) come from a
    Krylov exponential on the SPARSE hopping matrix instead of a dense expm, and
    the O(M^2) double sums are collapsed using the block structure of nu2:

        sum_ab f_a g_b nu2_ab = (f.nu)(g.nu)
                              + sum_blocks sum_{a,b in block} f_a g_b C_ab

    with C the correction returned by blocks_of. Both pieces are O(M), so the
    whole correlator is linear in the lattice size.
    """
    M = hsp.shape[0]
    rows = {}
    for p_ in (idx[(site_i, 0)], idx[(site_i, 1)],
               idx[(site_j, 0)], idx[(site_j, 1)]):
        e = np.zeros(M, dtype=complex)
        e[p_] = 1.0
        # h is real symmetric, so row p of exp(-iht) is its column p
        rows[p_] = expm_multiply(-1j * hsp * t, e)

    local = []
    for (x, y) in dimers:
        au, ad, bu, bd = (idx[(x, 0)], idx[(x, 1)], idx[(y, 0)], idx[(y, 1)])
        g = [au, ad, bu, bd]
        loc = {m_: k for k, m_ in enumerate(g)}
        local.append((np.array(g), (1 << loc[au]) | (1 << loc[bd]),
                      (1 << loc[ad]) | (1 << loc[bu])))

    def n1(p_):
        return float(np.real((np.abs(rows[p_]) ** 2) @ nu))

    def n2(p_, q_):
        rp, rq = rows[p_], rows[q_]
        ap, aq = np.abs(rp) ** 2, np.abs(rq) ** 2
        f, g_ = np.conj(rp) * rq, rp * np.conj(rq)
        d1 = (ap @ nu) * (aq @ nu)
        d2 = (f @ nu) * g_.sum() - (f @ nu) * (g_ @ nu)
        for (gg, C) in blocks:
            d1 += ap[gg] @ C @ aq[gg]
            d2 -= f[gg] @ C @ g_[gg]
        x = 0.0 + 0j
        for (gg, mA, mB) in local:
            rp_g, rq_g = rp[gg], rq[gg]
            for a, b, c_, e in itertools.product(range(4), repeat=4):
                K = np.conj(rp_g[a]) * rp_g[b] * np.conj(rq_g[c_]) * rq_g[e]
                if K == 0:
                    continue
                el = (four_mode_element(mA, mB, (a, b, c_, e))
                      + four_mode_element(mB, mA, (a, b, c_, e)))
                if el:
                    x -= 0.5 * K * el
        return float(np.real(d1 + d2 + x))

    iu, idn = idx[(site_i, 0)], idx[(site_i, 1)]
    ju, jdn = idx[(site_j, 0)], idx[(site_j, 1)]
    sz_i, sz_j = n1(iu) - n1(idn), n1(ju) - n1(jdn)
    nn = n2(iu, ju) - n2(iu, jdn) - n2(idn, ju) + n2(idn, jdn)
    return nn - sz_i * sz_j


# --------------------------------------------------------------------------
# exact many-body reference


def fock_basis(M: int, N: int):
    states = [sum(1 << p for p in c) for c in itertools.combinations(range(M), N)]
    return states, {s: k for k, s in enumerate(states)}


def many_body_H(h: np.ndarray, states, index):
    M = h.shape[0]
    rows, cols, vals = [], [], []
    for k, s in enumerate(states):
        for b in range(M):
            if not (s >> b) & 1:
                continue
            s1, g1 = _apply("c", b, s)
            for a in range(M):
                if h[a, b] == 0 or (s1 >> a) & 1:
                    continue
                s2, g2 = _apply("d", a, s1)
                rows.append(index[s2]); cols.append(k)
                vals.append(h[a, b] * g1 * g2)
    return csr_matrix((vals, (rows, cols)), shape=(len(states), len(states)),
                      dtype=complex)


def czz_exact(t, h, idx, dimers, doublon, site_i, site_j):
    M = h.shape[0]
    occ = []
    for (x, y) in dimers:
        occ.append(((idx[(x, 0)], idx[(y, 1)]), (idx[(x, 1)], idx[(y, 0)])))
    occ.append(((idx[(doublon, 0)], idx[(doublon, 1)]),))
    N = sum(len(o[0]) for o in occ)
    states, index = fock_basis(M, N)
    psi = np.zeros(len(states), dtype=complex)
    for choice in itertools.product(*[range(len(o)) for o in occ]):
        mask, sgn, amp = 0, 1, 1.0
        for gi, ci in enumerate(choice):
            for m in occ[gi][ci]:
                mask |= 1 << m
            if len(occ[gi]) == 2:
                amp *= (1.0 if ci == 0 else -1.0) / math.sqrt(2.0)
        # canonical ordering: build the sign by creating in increasing index
        st, s = 0, 1
        for m in sorted(p for p in range(M) if (mask >> p) & 1):
            st, g = _apply("d", m, st)
            s *= g
        psi[index[st]] += amp * s
    psi /= np.linalg.norm(psi)
    H = many_body_H(h, states, index)
    pt = expm_multiply(-1j * H * t, psi)

    def nop(p):
        return np.array([1.0 if (s >> p) & 1 else 0.0 for s in states])
    iu, idn = idx[(site_i, 0)], idx[(site_i, 1)]
    ju, jdn = idx[(site_j, 0)], idx[(site_j, 1)]
    w = np.abs(pt) ** 2
    szi = (nop(iu) - nop(idn))
    szj = (nop(ju) - nop(jdn))
    return float(w @ (szi * szj) - (w @ szi) * (w @ szj))


# --------------------------------------------------------------------------


def validate(sizes, times, periodic=True):
    out = []
    for (Lx, Ly) in sizes:
        order, idx, dimers, doublon, holon, spare = mode_order(Lx, Ly)
        h = hopping(Lx, Ly, idx, periodic=periodic)
        nu, nu2 = occupations(idx, dimers, doublon, holon, spare)
        si, sj = dimers[0]
        for t in times:
            t0 = time.time()
            a = czz_free(t, h, nu, nu2, idx, dimers, si, sj)
            t_poly = time.time() - t0
            t0 = time.time()
            b = czz_exact(t, h, idx, dimers, doublon, si, sj)
            t_exact = time.time() - t0
            out.append({"Lx": Lx, "Ly": Ly, "m": Lx * Ly, "t": t,
                        "czz_poly": a, "czz_exact": b, "abs_err": abs(a - b),
                        "s_poly": t_poly, "s_exact": t_exact})
            print(f"  m={Lx * Ly:3d} ({Lx}x{Ly})  t={t:4.2f}  "
                  f"poly {a:+.10f}  exact {b:+.10f}  |d| {abs(a - b):.2e}  "
                  f"[{t_poly * 1e3:7.1f} ms vs {t_exact * 1e3:9.1f} ms]")
    return out


def scaling(sizes, t=1.0, periodic=True, cross_check_to=64):
    """Measured wall time of the O(m) estimator, where exact is impossible."""
    out = []
    for (Lx, Ly) in sizes:
        m = Lx * Ly
        order, idx, dimers, doublon, holon, spare = mode_order(Lx, Ly)
        nu = nu_only(idx, dimers, doublon, holon, spare)
        blocks = blocks_of(idx, dimers, doublon, holon, spare, nu)
        hsp = hopping_sparse(Lx, Ly, idx, periodic=periodic)
        si, sj = dimers[0]
        t0 = time.time()
        v = czz_free_sparse(t, hsp, nu, blocks, idx, dimers, si, sj)
        el = time.time() - t0
        row = {"m": m, "t": t, "czz": v, "seconds": el}
        if m <= cross_check_to:      # the O(m) form must reproduce the O(M^2) one
            nu_d, nu2_d = occupations(idx, dimers, doublon, holon, spare)
            row["czz_dense"] = czz_free(t, hopping(Lx, Ly, idx, periodic=periodic),
                                        nu_d, nu2_d, idx, dimers, si, sj)
            row["dense_err"] = abs(row["czz_dense"] - v)
        out.append(row)
        tag = f"  (dense {row['czz_dense']:+.8f}, |d| {row['dense_err']:.1e})" \
            if "czz_dense" in row else ""
        print(f"  m={m:6d} ({Lx}x{Ly})  C^zz = {v:+.8f}   {el:8.3f} s{tag}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--open", action="store_true", help="open boundaries")
    args = ap.parse_args()
    per = not args.open
    print("U = 0 validation: polynomial estimator vs exact many-body evolution")
    print(f"(periodic={per}, J=1, dimer covering along x, one doublon + one holon)\n")
    v = validate([(2, 2), (3, 2), (4, 2), (3, 3)], [0.0, 0.25, 0.5, 1.0, 2.0], per)
    print("\nwall-time scaling of the polynomial estimator (no exact reference):")
    s = scaling([(4, 4), (8, 8), (16, 16), (32, 32), (64, 64), (128, 128),
                 ], 1.0, per)
    worst = max(x["abs_err"] for x in v)
    print(f"\nworst |C_poly - C_exact| = {worst:.3e}")
    p = pathlib.Path(__file__).parent / "data" / "free_fermion.json"
    p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps({"validation": v, "scaling": s,
                             "worst_abs_err": worst, "periodic": per}, indent=1))
    print(f"wrote {p}")

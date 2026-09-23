#!/usr/bin/env python3
"""Exact Trotter-error calibration for 2D Fermi-Hubbard, review finding #5.

WHAT THIS MEASURES
  The resource model derives the Trotter step count r from Campbell's worst-case
  commutator bound. Three independent lines (a collaborator's model, the
  Phasecraft/Quantinuum experiment, and this review) say that bound is far too
  loose. This script measures the ACTUAL error, for the actual observable and
  the actual initial state, by exact evolution on small patches, and reports
  kappa = (bound-implied error) / (measured error).

CONVENTIONS MATCHED TO THE COLLABORATOR'S small_dynamics.py
  - rectangle(nx,ny) snake (boustrophedon) site ordering
  - hopping edges two-coloured by x%2 / y%2
  - second-order Trotter as the symmetric 5-stage sweep
        onsite(dt/2) pink(dt/2) gold(dt) pink(dt/2) onsite(dt/2)
  - fixed (N_up, N_down) sector
  Their patches cross5 / square2 / square3 / rectangle2x3 are included so the
  two implementations can be cross-checked.

WHAT DIFFERS, DELIBERATELY
  Their calibration uses a Neel start and the local doublon density. Ours must
  use the experiment's state and observable (arXiv:2510.26300): a dimerised
  S^z_tot = 0 triplet covering with one holon and one doublon, and the
  equal-time connected correlator C^zz_ij = 4(<S_i S_j> - <S_i><S_j>).
  A Neel start is a Z eigenstate and gives C^zz_c == 0 identically.

MULTIPRODUCT
  Costed as CLASSICAL EXTRAPOLATION OF EXPECTATION VALUES: run separate circuits
  at step counts k*r, combine the resulting numbers with Lagrange weights. Each
  branch is its own circuit with its own depth -- which is the cost the resource
  model was missing.
"""
from __future__ import annotations
import json, math, pathlib, sys, time
from itertools import combinations
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

U_OVER_J = 4.0


def rectangle(nx, ny):                      # their snake ordering, verbatim
    return [(x, y) for y in range(ny) for x in (range(nx) if y % 2 == 0 else range(nx - 1, -1, -1))]


def sector_masks(n, k):
    return np.array([sum(1 << i for i in occ) for occ in combinations(range(n), k)], dtype=np.int64)


def spin_hopping(masks, edges):
    lookup = {int(m): i for i, m in enumerate(masks)}
    rows, cols, vals = [], [], []
    for col, mask in enumerate(masks):
        mask = int(mask)
        for i, j in edges:
            if ((mask >> i) ^ (mask >> j)) & 1:
                lo, hi = sorted((i, j))
                between = ((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1)
                sign = (-1) ** (bin(mask & between).count("1"))
                rows.append(lookup[mask ^ (1 << i) ^ (1 << j)]); cols.append(col); vals.append(-sign)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(len(masks),) * 2, dtype=float)


class Patch:
    """Fixed-sector Hubbard patch, factorised so nothing of size dim^2 is formed."""

    def __init__(self, sites, u=U_OVER_J):
        self.sites = list(sites); self.n = len(sites)
        idx = {xy: i for i, xy in enumerate(self.sites)}
        self.dimers, self.holon, self.doublon, self.singles = self._cover()
        n_up = 1 + len(self.dimers) + sum(1 for s in self.singles if s[1] == 0)
        n_dn = 1 + len(self.dimers) + sum(1 for s in self.singles if s[1] == 1)
        self.up, self.dn = sector_masks(self.n, n_up), sector_masks(self.n, n_dn)
        self.nu, self.nd = len(self.up), len(self.dn)
        self.dim = self.nu * self.nd
        edges = [[], []]
        for i, (x, y) in enumerate(self.sites):
            for xy, colour in (((x + 1, y), x % 2), ((x, y + 1), y % 2)):
                if xy in idx:
                    edges[colour].append((i, idx[xy]))
        self.edges = edges
        self.hop = [(spin_hopping(self.up, g), spin_hopping(self.dn, g)) for g in edges]
        occ_u = ((self.up[:, None] >> np.arange(self.n)) & 1).astype(np.int8)
        occ_d = ((self.dn[:, None] >> np.arange(self.n)) & 1).astype(np.int8)
        self.occ_u, self.occ_d = occ_u, occ_d
        # V = u * sum_i (n_iu - 1/2)(n_id - 1/2), as a (nu, nd) outer sum
        self.V = u * ((occ_u - .5).astype(np.float64) @ (occ_d - .5).astype(np.float64).T)
        self.psi0 = self._initial()

    def _cover(self):
        """Holon, doublon, then dimers along the snake; odd leftovers stay single."""
        order = list(range(self.n))
        holon, doublon = order[0], order[1]
        rest = order[2:]
        dimers = [(rest[i], rest[i + 1]) for i in range(0, len(rest) - 1, 2)]
        singles = [(rest[-1], 0)] if len(rest) % 2 else []
        return dimers, holon, doublon, singles

    def _initial(self):
        """Dimerised S^z_tot=0 triplet covering + one holon + one doublon."""
        psi = np.zeros((self.nu, self.nd), dtype=complex)
        up_i = {int(m): i for i, m in enumerate(self.up)}
        dn_i = {int(m): i for i, m in enumerate(self.dn)}
        base_u = 1 << self.doublon
        base_d = 1 << self.doublon
        for s, sp in self.singles:
            (base_u := base_u | (1 << s)) if sp == 0 else (base_d := base_d | (1 << s))
        amp = 2.0 ** (-len(self.dimers) / 2)
        for bits in range(1 << len(self.dimers)):
            mu, md = base_u, base_d
            for k, (a, b) in enumerate(self.dimers):
                if (bits >> k) & 1:   mu |= 1 << a; md |= 1 << b     # up on a, down on b
                else:                 mu |= 1 << b; md |= 1 << a     # the other triplet term
            psi[up_i[mu], dn_i[md]] += amp                            # + sign => triplet
        nrm = np.linalg.norm(psi)
        assert abs(nrm - 1) < 1e-10, nrm
        return psi

    # ---- observables ----
    def czz(self, P, i, j):
        w = np.abs(P) ** 2
        su_i, sd_i = self.occ_u[:, i], self.occ_d[:, i]
        su_j, sd_j = self.occ_u[:, j], self.occ_d[:, j]
        Si = 0.5 * (su_i[:, None] - sd_i[None, :])
        Sj = 0.5 * (su_j[:, None] - sd_j[None, :])
        e_ij = float((w * Si * Sj).sum()); e_i = float((w * Si).sum()); e_j = float((w * Sj).sum())
        return 4.0 * (e_ij - e_i * e_j)

    def doublon_density(self, P, site_index):
        w = np.abs(P) ** 2
        return float(w[self.occ_u[:, site_index] == 1][:, self.occ_d[:, site_index] == 1].sum())

    # ---- dynamics (all matvec-based; nothing of size dim x dim is built) ----
    def _H(self, P):
        out = self.V * P
        for hu, hd in self.hop:
            out = out + hu @ P + P @ hd.T
        return out

    def exact(self, tau, krylov=None, sub=None):
        """Lanczos/Arnoldi exponential, substepped. Validated against scipy below.

        The Krylov basis holds `krylov` full state vectors. At 14 sites that is
        188 MB each, so a fixed 40 would want 7.7 GB. The basis is capped by
        KRYLOV_BUDGET_GB and the substep count raised to compensate -- accuracy
        comes from substepping, memory from the cap.
        """
        vec_gb = self.dim * 16 / 2 ** 30
        if krylov is None:
            krylov = int(max(8, min(40, KRYLOV_BUDGET_GB / max(vec_gb, 1e-9))))
        # A Krylov space of dimension k resolves exp(-iHt) only while
        # ||H|| dt <~ k/3. Capping k for memory therefore REQUIRES more substeps,
        # and the first version of this did not: at n = 14 the cap took k to 8,
        # left sub at 10, and the "exact" reference came out wrong by 2e-2 while
        # every smaller patch was at 1e-15. Scale the substep by the spectral
        # radius, not by tau alone.
        if sub is None:
            nrm = self._hnorm()
            sub = max(1, int(np.ceil(tau * 4)),
                      int(np.ceil(3.0 * tau * nrm / max(krylov, 1))))
        P = self.psi0.copy()
        for _ in range(sub):
            P = self._expv(P, tau / sub, krylov)
        return P

    def _hnorm(self):
        """Cheap upper bound on ||H||: hopping bandwidth plus the on-site term."""
        if getattr(self, "_hn", None) is None:
            hop = sum(max(abs(A).sum(axis=1).max(), abs(B).sum(axis=1).max())
                      for A, B in self.hop)
            self._hn = float(hop + np.abs(self.V).max())
        return self._hn

    def _expv(self, P, dt, m):
        beta = np.linalg.norm(P)
        Vk = [P / beta]; H = np.zeros((m + 1, m + 1), complex)
        for j in range(m):
            w = self._H(Vk[j])
            for i in range(max(0, j - 1), j + 1):
                H[i, j] = np.vdot(Vk[i], w); w = w - H[i, j] * Vk[i]
            H[j + 1, j] = np.linalg.norm(w)
            if H[j + 1, j].real < 1e-13:
                m = j + 1; break
            Vk.append(w / H[j + 1, j])
        from scipy.linalg import expm
        E = expm(-1j * dt * H[:m, :m])
        out = np.zeros_like(P)
        for i in range(m):
            out += beta * E[i, 0] * Vk[i]
        return out

    def trotter(self, tau, steps):
        dt = tau / steps
        P = self.psi0.copy()
        ons = np.exp(-0.5j * dt * self.V)
        (pu, pd), (gu, gd) = self.hop
        def ev(A, B, Q, c):
            Q = expm_multiply((c * 1j) * A, Q) if A.nnz else Q
            Q = expm_multiply((c * 1j) * B, Q.T).T if B.nnz else Q
            return Q
        for _ in range(steps):
            P = ons * P
            P = ev(pu, pd, P, -0.5 * dt)
            P = ev(gu, gd, P, -1.0 * dt)
            P = ev(pu, pd, P, -0.5 * dt)
            P = ons * P
        return P


# ======================================================================= driver
PATCHES = {
    "square2":      rectangle(2, 2),      # 4  -- theirs
    "cross5":       [(1, 0), (0, 1), (1, 1), (2, 1), (1, 2)],   # 5 -- theirs
    "rectangle2x3": rectangle(2, 3),      # 6  -- theirs
    "square3":      rectangle(3, 3),      # 9  -- theirs
    "rectangle3x4": rectangle(3, 4),      # 12
    "rectangle2x7": rectangle(2, 7),      # 14 -- fits locally, 176 MB/vector
    "square4":      rectangle(4, 4),      # 16 -- 2.5 GB/vector, lenore only
}

# The ADOPTED trajectory: the model evolves to t = sqrt(m)/v_B, so the
# calibration has to be measured there too. Everything else in this file is a
# fixed-time size sweep, which is a different -- and weaker -- statement
# (second-pass review #2, test 3).
V_B = 2.0
TRAJECTORY = ("square2", "cross5", "rectangle2x3", "square3", "rectangle3x4",
              "rectangle2x7", "square4")
TRAJ_STEPS = (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64)
KRYLOV_BUDGET_GB = 1.5      # cap the Lanczos basis; substep to make up accuracy
TAUS = (0.25, 0.5, 1.0, 2.0)
STEPS = (1, 2, 4, 8, 16, 32, 64)
MEM_CAP_GB = 3.0        # this machine has ~10 GB free; n = 16 belongs on lenore
CAMPBELL_W = 9.5            # commutator norm per site at U/J = 4


def mp_weights(k):
    """Lagrange weights extrapolating h^2 -> 0 from step counts 1..k."""
    x = [(1.0 / i) ** 2 for i in range(1, k + 1)]
    out = []
    for i, xi in enumerate(x):
        c = 1.0
        for j, xj in enumerate(x):
            if i != j: c *= xj / (xj - xi)
        out.append(c)
    return out


def trajectory(only=None, out="data/trotter_traj.json"):
    """Measure W_eff AT the operating points, tau = sqrt(m)/v_B (review #2, test 3).

    The whole m^(7/4) exponent rests on W_eff being flat along this trajectory.
    Every other measurement in this file is a fixed-time size sweep, which shows
    something weaker: locality at SHORT time, decaying to a 2.9x spread by
    tau = 2 (calibration/domain_check.py).

    Each patch is also checked for convergence of the exact reference by halving
    the substep, because a Krylov basis capped for memory is only as good as the
    substepping that compensates for it.
    """
    names = only or list(TRAJECTORY)
    rows, meta = [], []
    for name in names:
        sites = PATCHES[name]
        P = Patch(sites)
        gb = P.dim * 16 / 2 ** 30
        if gb * 8 > MEM_CAP_GB:
            print(f"SKIP {name}: {gb:.2f} GB/vector x 8 Krylov exceeds the "
                  f"{MEM_CAP_GB} GB cap -- run this one on lenore", flush=True)
            continue
        tau = math.sqrt(P.n) / V_B
        i, j = P.dimers[0]
        t0 = time.monotonic()
        ex = P.exact(tau)
        ex2 = P.exact(tau, sub=max(2, int(np.ceil(tau * 8))))
        conv = float(np.linalg.norm(ex - ex2))
        cz_ex = P.czz(ex, i, j)
        print(f"[{name}] n={P.n} tau={tau:.3f} dim={P.dim:,} "
              f"({gb * 1024:.1f} MB/vec)  C^zz_exact = {cz_ex:+.8f}  "
              f"substep convergence {conv:.2e}", flush=True)
        vals = {}
        for r in TRAJ_STEPS:
            ap = P.trotter(tau, r)
            vals[r] = P.czz(ap, i, j)
            err = abs(vals[r] - cz_ex)
            rows.append({"patch": name, "n": P.n, "tau": tau, "steps": r,
                         "czz_exact": cz_ex, "czz_trotter": vals[r],
                         "abs_err": err,
                         "W_eff": err * r ** 2 / tau ** 3 if err > 0 else None,
                         "infidelity": float(max(0.0, 1 - abs(np.vdot(ex, ap)) ** 2))})
        meta.append({"patch": name, "n": P.n, "tau": tau, "dim": P.dim,
                     "gb_per_vec": gb, "czz_exact": cz_ex,
                     "substep_convergence": conv,
                     "seconds": time.monotonic() - t0})
        good = [x["W_eff"] for x in rows
                if x["patch"] == name and x["steps"] >= 8 and x["W_eff"]]
        if good:
            print(f"   W_eff(r>=8) = {float(np.median(good)):.5f}", flush=True)
    path = pathlib.Path(__file__).parent / out
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps({"rows": rows, "meta": meta, "v_b": V_B,
                                "u_over_j": U_OVER_J}, indent=1))
    print(f"\nwrote {path}")
    return rows, meta


CLAMP_TAUS = (2.5, 3.0, 4.0)


def clamp_probe(patches=("rectangle3x4", "rectangle2x7"), out="data/clamp_probe.json"):
    """W_eff BEYOND the calibrated range, on two lattice sizes (review #2).

    alpha = 1.75 comes from CLAMPING W at its tau = 2 value, and every plotted
    point sits in the clamped region -- so the clamp, not the measured
    tau-dependence, is what sets the exponent. This is the only affordable test
    of it.

    Two sizes, because a small lattice cannot honestly be asked about long times:
    at tau = 4 the cone spans 2 v tau = 16 sites, wider than a 12-site patch, so
    the answer wraps. If n = 12 and n = 14 AGREE the tau-trend is real and the
    clamp is conservative; if they DIVERGE, small lattices cannot test it at all
    and the clamp stands for want of evidence. Either outcome is decisive.
    """
    rows, meta = [], []
    for name in patches:
        P = Patch(PATCHES[name])
        gb = P.dim * 16 / 2 ** 30
        if gb * 8 > MEM_CAP_GB:
            print(f"SKIP {name}: {gb:.2f} GB/vector -- lenore", flush=True)
            continue
        i, j = P.dimers[0]
        for tau in CLAMP_TAUS:
            t0 = time.monotonic()
            ex = P.exact(tau)
            conv = float(np.linalg.norm(ex - P.exact(tau, sub=None if False else
                                                     int(np.ceil(6.0 * tau * P._hnorm() / 8)))))
            cz = P.czz(ex, i, j)
            best = []
            for r in (8, 12, 16, 24, 32, 48, 64):
                err = abs(P.czz(P.trotter(tau, r), i, j) - cz)
                if err > 0:
                    best.append(err * r ** 2 / tau ** 3)
                rows.append({"patch": name, "n": P.n, "tau": tau, "steps": r,
                             "czz_exact": cz, "abs_err": err,
                             "W_eff": err * r ** 2 / tau ** 3 if err > 0 else None})
            W = float(np.median(best)) if best else float("nan")
            meta.append({"patch": name, "n": P.n, "tau": tau, "W_eff": W,
                         "czz_exact": cz, "substep_convergence": conv,
                         "seconds": time.monotonic() - t0})
            print(f"[{name}] n={P.n} tau={tau:.2f}  C^zz={cz:+.8f}  "
                  f"W_eff={W:.5f}  conv={conv:.1e}  "
                  f"({time.monotonic()-t0:.0f}s)", flush=True)
    path = pathlib.Path(__file__).parent / out
    path.write_text(json.dumps({"rows": rows, "meta": meta,
                                "clamp_taus": list(CLAMP_TAUS)}, indent=1))
    print(f"\nwrote {path}")
    return meta


def main():
    if "--clamp" in sys.argv:
        args = [a for a in sys.argv[1:] if not a.startswith("-")]
        return clamp_probe(tuple(args) if args else
                           ("rectangle3x4", "rectangle2x7"))
    if "--trajectory" in sys.argv:
        args = [a for a in sys.argv[1:] if not a.startswith("-")]
        return trajectory(args or None)
    only = sys.argv[1:] or list(PATCHES)
    rows, meta = [], []
    for name in only:
        sites = PATCHES[name]
        t0 = time.monotonic()
        P = Patch(sites)
        gb = P.dim * 16 / 2**30
        if gb > MEM_CAP_GB:
            print(f"SKIP {name}: {gb:.1f} GB/vector exceeds cap", flush=True); continue
        i, j = P.dimers[0]
        meta.append({"patch": name, "n": P.n, "dim": P.dim, "gb_per_vec": gb,
                     "n_dimers": len(P.dimers), "probe_link": [i, j],
                     "build_s": time.monotonic() - t0})
        print(f"[{name}] n={P.n} dim={P.dim:,} ({gb:.3f} GB/vec) dimers={len(P.dimers)}", flush=True)
        c0 = P.czz(P.psi0, i, j)
        print(f"   C^zz(t=0) = {c0:+.6f}   (triplet algebra predicts -1)", flush=True)
        for tau in TAUS:
            te = time.monotonic()
            ex = P.exact(tau)
            cz_ex = P.czz(ex, i, j)
            vals = {}
            for r in STEPS:
                ap = P.trotter(tau, r)
                vals[r] = P.czz(ap, i, j)
                rows.append({"patch": name, "n": P.n, "tau": tau, "steps": r,
                             "czz_exact": cz_ex, "czz_trotter": vals[r],
                             "abs_err": abs(vals[r] - cz_ex),
                             "infidelity": float(max(0.0, 1 - abs(np.vdot(ex, ap)) ** 2))})
            # classical multiproduct extrapolation of the EXPECTATION VALUES
            for k in (2, 3, 4):
                w = mp_weights(k)
                for base in (1, 2, 4, 8):
                    ks = [base * i for i in range(1, k + 1)]
                    for ki in ks:                       # branches are their own circuits
                        if ki not in vals:
                            vals[ki] = P.czz(P.trotter(tau, ki), i, j)
                    est = sum(wi * vals[ki] for wi, ki in zip(w, ks))
                    rows.append({"patch": name, "n": P.n, "tau": tau,
                                 "steps": -base, "mp_order": 2 * k,
                                 "mp_branches": ks, "mp_l1": sum(abs(x) for x in w),
                                 "czz_exact": cz_ex, "czz_trotter": est,
                                 "abs_err": abs(est - cz_ex),
                                 "deepest_branch_steps": max(ks)})
            print(f"   tau={tau}: exact {cz_ex:+.6f}  err(r=8) {abs(vals[8]-cz_ex):.3e}"
                  f"  err(r=64) {abs(vals[64]-cz_ex):.3e}  [{time.monotonic()-te:.1f}s]", flush=True)
        json.dump({"rows": rows, "meta": meta, "campbell_w": CAMPBELL_W,
                   "u_over_j": U_OVER_J}, open("trotter_cal.json", "w"), indent=1)
    print("wrote trotter_cal.json", flush=True)


if __name__ == "__main__":
    main()

"""v2 of axiom212.py: same mathematics, restructured for speed.

* exact rationals via FLINT fmpq (C), scalar arithmetic only (no FLINT polynomial multiplication, so the
  fmpz_poly_mul_SS bug fixed by FLINT PR 2790 cannot affect us);
* per stratum: cells, slabs, moment vectors computed once;
  omega(P,Q,e,f) = sum_slabs  E_{P,f} . Itilde_{Q,e}, with Itilde = moment-convolved y-integral (cached);
* Omega(cell, Lam, e, f) = sum_{P,Q} r_Lam(P,Q) omega(...) cached;
* strata are independent -> multiprocessing.
Validated against v1 (axiom212.py) entry by entry, see test_v2().
"""
from flint import fmpq
from functools import lru_cache
from math import comb, factorial
from itertools import combinations
from collections import defaultdict
import os

from simplex_sieve import p_to_m, n_monomials, _canon

Z = fmpq(0)
ONE = fmpq(1)


def Fr(a, b=1):
    return fmpq(a, b)


class Params:
    def __init__(self, K=45, L=Fr(53, 50), g=Fr(41, 625), tau=Fr(249, 250), caps=None):
        self.K, self.L, self.g, self.tau = K, L, g, tau
        B1 = [Fr(b, 5000) for b in (777, 794, 875, 917, 953, 983, 1016, 1042, 1063, 1081)]
        self.caps = caps if caps is not None else ([Z] + [4 * b for b in B1] + [4 * Fr(1081, 5000)] * 40)

    def C(self, m):
        return Z if m <= 0 else self.caps[min(m, len(self.caps) - 1)]


@lru_cache(maxsize=None)
def D(lam, k):
    if k == 0:
        return 1 if not lam else 0
    tot = 0
    for mu, c in p_to_m(lam, k).items():
        if len(mu) > k:
            continue
        w = c * n_monomials(mu, k)
        for e in mu:
            w *= factorial(e)
        tot += w
    return tot


def expand_shift(Lam, s, t, nplain, g):
    res = {((), (), ()): ONE}
    for j in Lam:
        opts = [(None, 0, (s + t) * g ** j)]
        for i in range(1, j + 1):
            cf = comb(j, i) * g ** (j - i)
            if s:
                opts.append(("U", i, cf))
            if t:
                opts.append(("Z", i, cf))
        if nplain:
            opts.append(("W", j, ONE))
        nxt = defaultdict(lambda: Z)
        for (lu, lz, lw), c in res.items():
            for grp, i, cf in opts:
                if cf == 0:
                    continue
                if grp is None:
                    key = (lu, lz, lw)
                elif grp == "U":
                    key = (_canon(lu + (i,)), lz, lw)
                elif grp == "Z":
                    key = (lu, _canon(lz + (i,)), lw)
                else:
                    key = (lu, lz, _canon(lw + (i,)))
                nxt[key] += c * cf
        res = dict(nxt)
    return res


def radial(Lam, s, t, n, g):
    out = defaultdict(lambda: Z)
    for (lu, lz, lw), c in expand_shift(Lam, s, t, n - t, g).items():
        du, dzw = sum(lu), sum(lz) + sum(lw)
        if s == 0:
            if lu:
                continue
            fa, P = ONE, -1
        else:
            d = D(lu, s)
            if not d:
                continue
            fa, P = Fr(d, factorial(s - 1 + du)), s - 1 + du
        db = D(lz, t) * D(lw, n - t)
        if not db:
            continue
        out[(P, n - 1 + dzw)] += c * fa * Fr(db, factorial(n - 1 + dzw))
    return {k: v for k, v in out.items() if v != 0}


# ---------------- I ----------------
def inc_beta(P, N, M, R):
    full = R ** (P + N + 1) * Fr(factorial(P) * factorial(N), factorial(P + N + 1))
    if M == R:
        return full
    d = R - M
    return full - sum((comb(P, k) * (-1) ** k * R ** (P - k) * d ** (N + k + 1) / (N + k + 1) for k in range(P + 1)), Z)


def I_strata(pr):
    out = []
    for s in range(0, pr.K + 1):
        if s and pr.C(s) - s * pr.g <= 0:
            break
        n = pr.K - s
        for t in range(0, n + 1):
            R1 = pr.L - (s + t) * pr.g
            if R1 <= 0:
                break
            out.append((s, t))
    return out


def I_stratum_values(pr, s, t, keys):
    """{(A, Lam): contribution of stratum (s,t)} for all keys."""
    n = pr.K - s
    R1 = pr.L - (s + t) * pr.g
    mult = comb(pr.K, s) * comb(n, t) * (-1) ** t
    M = min(pr.C(s) - s * pr.g, R1) if s else None
    out = {}
    rcache = {}
    for A, Lam in keys:
        if Lam not in rcache:
            rcache[Lam] = radial(Lam, s, t, n, pr.g)
        acc = Z
        for (P, Q), r in rcache[Lam].items():
            N = Q + A + 1
            fq = Fr(factorial(Q) * factorial(A), factorial(N))
            acc += r * fq * (R1 ** N if s == 0 else inc_beta(P, N, M, R1))
        out[(A, Lam)] = mult * acc
    return out


# ---------------- J ----------------
def phi_terms(alpha, m, cellU):
    if cellU == "V":
        return [(alpha + m + 1, 0, Fr(factorial(alpha) * factorial(m), factorial(alpha + m + 1)))]
    return [(alpha - k, m + k + 1, Fr(comb(alpha, k) * (-1) ** k, m + k + 1)) for k in range(alpha + 1)]


def inner_terms(b, cellU):
    alpha, lam = b
    out = defaultdict(lambda: Z)
    idx = range(len(lam))
    for r in range(len(lam) + 1):
        for A in combinations(idx, r):
            m = sum(lam[i] for i in A)
            nu = _canon([lam[i] for i in idx if i not in A])
            for eV, fX, c in phi_terms(alpha, m, cellU):
                out[(nu, eV, fX)] += c
    return [(k, v) for k, v in out.items() if v != 0]


def _affpow(c0, c1, e):
    return [comb(e, i) * c0 ** (e - i) * c1 ** i for i in range(e + 1)]


def _polymul(p, q):
    r = [Z] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x != 0:
            for j, y in enumerate(q):
                r[i + j] += x * y
    return r


class Stratum:
    def __init__(self, pr, s, t):
        self.pr = pr
        g, L, tau = pr.g, pr.L, pr.tau
        self.s, self.t, self.n = s, t, (pr.K - 1) - s
        self.R1 = L - (s + t) * g
        self.Tp = tau - (s + t) * g
        self.Cap = pr.C(s) - s * g if s else Z
        self.Wc = pr.C(s + 1) - s * g
        self.aW = self.Wc - g
        self.bV = L - pr.C(s + 1) - t * g
        self.cg = self.R1 - g
        self.ok = self.Tp > 0 and (s == 0 or self.Cap > 0)
        self.cells = {}
        if self.ok:
            base = [(-1, 0, Z), (0, -1, Z), (1, 1, self.Tp)] + ([(1, 0, self.Cap)] if s else [(1, 0, Z)])
            defs = {1: ("V", [(1, 0, self.aW), (0, -1, -self.bV)]),
                    2: ("W", [(1, 0, self.aW), (0, 1, self.bV)]),
                    3: ("V", [(-1, 0, -self.aW), (-1, -1, -self.cg)]),
                    4: ("g", [(-1, 0, -self.aW), (1, 1, self.cg)])}
            for k, (u, cons) in defs.items():
                cons = base + cons
                sl = self._slabs(cons) if s else self._range1d(cons)
                if sl:
                    self.cells[k] = (u, sl)
        self._E, self._It, self._om, self._Om, self._mom = {}, {}, {}, {}, {}

    @staticmethod
    def _vertices(cons):
        pts = set()
        for i in range(len(cons)):
            for j in range(i + 1, len(cons)):
                a1, b1, c1 = cons[i]; a2, b2, c2 = cons[j]
                det = a1 * b2 - a2 * b1
                if det == 0:
                    continue
                x = (c1 * b2 - c2 * b1) / Fr(det)
                y = (a1 * c2 - a2 * c1) / Fr(det)
                if all(al * x + be * y <= ga for al, be, ga in cons):
                    pts.add((x, y))
        return pts

    def _slabs(self, cons):
        V = self._vertices(cons)
        if len(V) < 3:
            return []
        xs = sorted({p[0] for p in V})
        out = []
        for a0, a1 in zip(xs, xs[1:]):
            am = (a0 + a1) / 2
            lo = hi = None
            for al, be, ga in cons:
                if be == 1:
                    val = ga - al * am
                    if hi is None or val < hi[0]:
                        hi = (val, al, ga)
                elif be == -1:
                    val = al * am - ga
                    if lo is None or val > lo[0]:
                        lo = (val, al, ga)
            if hi[0] <= lo[0]:
                continue
            out.append((a0, a1, self._yb(lo, -1), self._yb(hi, 1)))
        return out

    def _yb(self, bnd, be):
        _, al, ga = bnd
        k0, k1 = (ga, -al) if be == 1 else (-ga, al)
        if 1 + k1 == 0:
            return ("y", self.R1 - k0)
        assert 1 + k1 == 1
        return ("c", k0)

    def _range1d(self, cons):
        lo, hi = Z, None
        for al, be, ga in cons:
            if be == 1:
                hi = ga if hi is None else min(hi, ga)
            elif be == -1:
                lo = max(lo, -ga)
            elif 0 > ga:
                return []
        if hi is None or hi <= lo:
            return []
        return [("1d", self.R1 - hi, self.R1 - lo)]  # y-range

    @staticmethod
    def _Fy(Q, e, bound):
        kind, v = bound
        if kind == "y":
            out = [Z] * (Q + 1)
            for k in range(Q + 1):
                out[Q - k] = comb(Q, k) * (-1) ** k * v ** (k + e + 1) / (k + e + 1)
            return out
        out = [Z] * (Q + e + 2)
        out[Q + e + 1] = Fr(factorial(Q) * factorial(e), factorial(Q + e + 1))
        if v != 0:
            for j in range(e + 1):
                out[e - j] -= comb(e, j) * (-1) ** j * v ** (Q + j + 1) / (Q + j + 1)
        return out

    def _moments(self, cell, si, upto):
        key = (cell, si)
        a0, a1 = self.cells[cell][1][si][:2]
        c0, c1 = self.R1 - a1, self.R1 - a0
        m = self._mom.get(key)
        if m is None or len(m) < upto:
            p0, p1 = [ONE], [ONE]
            m = []
            for k in range(upto):
                p0.append(p0[-1] * c0); p1.append(p1[-1] * c1)
                m.append((p1[-1] - p0[-1]) / (k + 1))
            self._mom[key] = m
        return m

    def _E_(self, P, f, Xkind):
        key = (P, f, Xkind)
        if key not in self._E:
            pa = _affpow(self.R1, -ONE, P)
            px = _affpow(self.Wc - self.R1, ONE, f) if Xkind == "W" else [self.pr.g ** f]
            self._E[key] = _polymul(pa, px)
        return self._E[key]

    def _Itilde(self, cell, si, Q, e, need):
        key = (cell, si, Q, e)
        it = self._It.get(key)
        if it is None or len(it) < need:
            lo, hi = self.cells[cell][1][si][2:]
            pl, ph = self._Fy(Q, e, lo), self._Fy(Q, e, hi)
            n = max(len(pl), len(ph))
            Ip = [(pl[j] if j < len(pl) else Z) - (ph[j] if j < len(ph) else Z) for j in range(n)]
            mom = self._moments(cell, si, n + need + 1)
            it = [sum((Ip[j] * mom[i + j] for j in range(n) if Ip[j] != 0), Z) for i in range(need)]
            self._It[key] = it
        return it

    def omega(self, cell, P, Q, e, f):
        key = (cell, P, Q, e, f)
        v = self._om.get(key)
        if v is not None:
            return v
        Ukind, slabs = self.cells[cell]
        Xkind = "W" if Ukind == "W" else "g"
        if P == -1:
            _, yl, yh = slabs[0]
            R1 = self.R1
            F = lambda y: sum((comb(Q, k) * (-1) ** k * R1 ** (Q - k) * y ** (k + e + 1) / (k + e + 1) for k in range(Q + 1)), Z)
            X = self.Wc if Xkind == "W" else self.pr.g
            v = X ** f * (F(yh) - F(yl))
        else:
            E = self._E_(P, f, Xkind)
            v = Z
            for si in range(len(slabs)):
                it = self._Itilde(cell, si, Q, e, len(E))
                v += sum((E[i] * it[i] for i in range(len(E)) if E[i] != 0), Z)
        self._om[key] = v
        return v

    def Omega(self, cell, Lam, e, f):
        key = (cell, Lam, e, f)
        v = self._Om.get(key)
        if v is None:
            v = Z
            for (P, Q), r in radial(Lam, self.s, self.t, self.n, self.pr.g).items():
                v += r * self.omega(cell, P, Q, e, f)
            self._Om[key] = v
        return v


def J_strata(pr):
    out = []
    for s in range(0, pr.K):
        if s and pr.C(s) - s * pr.g <= 0:
            break
        n = (pr.K - 1) - s
        for t in range(0, n + 1):
            if pr.tau - (s + t) * pr.g <= 0:
                break
            out.append((s, t))
    return out


def J_stratum_matrix(pr, s, t, basis):
    st = Stratum(pr, s, t)
    nb = len(basis)
    Jm = [[Z] * nb for _ in range(nb)]
    if not st.ok:
        return Jm
    n = (pr.K - 1) - s
    mult = comb(pr.K - 1, s) * comb(n, t) * (-1) ** t * pr.K
    for cell, (Ukind, _) in st.cells.items():
        terms = [inner_terms(b, Ukind) for b in basis]
        for i in range(nb):
            for j in range(i, nb):
                acc = Z
                for (nu1, e1, f1), c1 in terms[i]:
                    for (nu2, e2, f2), c2 in terms[j]:
                        acc += c1 * c2 * st.Omega(cell, _canon(nu1 + nu2), e1 + e2, f1 + f2)
                if acc != 0:
                    Jm[i][j] += mult * acc
    for i in range(nb):
        for j in range(i):
            Jm[i][j] = Jm[j][i]
    return Jm


# ---------------- driver with multiprocessing ----------------
def _ser(x):
    return (int(x.p), int(x.q))


def _de(t):
    return fmpq(t[0], t[1])


def _job_I(args):
    pr, s, t, keys = args
    return {k: _ser(v) for k, v in I_stratum_values(pr, s, t, keys).items()}


def _job_J(args):
    pr, s, t, basis = args
    return [[_ser(x) for x in row] for row in J_stratum_matrix(pr, s, t, basis)]


def gram(basis, pr=None, workers=None, verbose=False):
    import time
    from multiprocessing import get_context
    pr = pr or Params()
    t0 = time.time()
    nb = len(basis)
    keys = sorted({(basis[i][0] + basis[j][0], _canon(basis[i][1] + basis[j][1])) for i in range(nb) for j in range(i, nb)})
    workers = workers or max(1, (os.cpu_count() or 2) - 2)
    ctx = get_context("fork")
    with ctx.Pool(workers) as pool:
        Ivals = defaultdict(lambda: Z)
        for d in pool.imap_unordered(_job_I, [(pr, s, t, keys) for s, t in I_strata(pr)]):
            for k, v in d.items():
                Ivals[k] += _de(v)
        if verbose:
            print(f"  I: {len(keys)} keys, {len(I_strata(pr))} strata, {time.time()-t0:.1f}s", flush=True)
        I = [[Ivals[(basis[i][0] + basis[j][0], _canon(basis[i][1] + basis[j][1]))] for j in range(nb)] for i in range(nb)]
        J = [[Z] * nb for _ in range(nb)]
        jobs = [(pr, s, t, basis) for s, t in J_strata(pr)]
        jobs.sort(key=lambda a: -(a[1] * 3 + (44 - a[2])))  # rough heavy first
        for k, M in enumerate(pool.imap_unordered(_job_J, jobs)):
            for i in range(nb):
                for j in range(nb):
                    J[i][j] += _de(M[i][j])
        if verbose:
            print(f"  J: {len(jobs)} strata, {time.time()-t0:.1f}s", flush=True)
    return I, J


def test_v2():
    """Compare with v1 on a small basis (exact equality)."""
    from fractions import Fraction
    import axiom212 as v1
    basis = [(0, ()), (1, ()), (0, (2,)), (2, ()), (1, (2,))]
    I1, J1 = v1.gram(basis)
    I2, J2 = gram(basis)
    eq = lambda a, b: Fraction(int(b.p), int(b.q)) == a
    ok = all(eq(I1[i][j], I2[i][j]) and eq(J1[i][j], J2[i][j]) for i in range(5) for j in range(5))
    print("v2 == v1 exactly:", ok)
    return ok


if __name__ == "__main__":
    test_v2()

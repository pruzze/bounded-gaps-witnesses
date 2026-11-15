"""Gram matrices of I_T and J_T on a capped support, for the even monomial basis B_d, in ball arithmetic.

Basis: b = (a, alpha) with alpha a partition into even parts, a + |alpha| <= d, F_b = (R - |x|)^a m_alpha(x).
(B_21 has dimension 846 and spans the same space as Axiom's B_21 of even power sums.)

Idea: every quadratic form used so far is linear in "keys":
    I_T(F)  = sum over (nu, A)        of  [F^2]_{nu,A} * WI(nu, A),         WI = int_T (R - |x|)^A m_nu
    J_T(F)  = sum over cells and keys of  [G^2]_{nu,key} * Psi_cell(nu, key)
so we tabulate WI and Psi once for all keys that B_d can produce and assemble the Gram matrices from them. The
per-stratum integrals are the same as in capped_mono.py / kintali_J.py / axiom212_v2.py (validated exact code);
here they are reorganised so that the fibre weights factor through two-way splits nu = lambda_U + mu and all
tables are computed with Arb balls (rigorous enclosures) instead of exact rationals.

Fibre weights. With S_e(x) = sum_j C(e,j) g^(e-j) j! x^j and prodS[lam] = prod_{e in lam} S_e:
    GU_s[lam]  = nmono(lam, s) prodS[lam]                      (rough group, s coordinates, degree index jU)
    GY[mu]     = sum_{lamZ + lamW = mu} nmono(lamZ, t) nmono(lamW, m) prod(lamW!) x^|lamW| prodS[lamZ]
    w_nu(P, Q) = sum_{lam + mu = nu} GU[lam]_jU GY[mu]_q / (P! Q!),   P = s - 1 + jU,  Q = n - 1 + q.
"""
import os
import pickle
import sys
import time
from collections import Counter, defaultdict
from itertools import product as iproduct
from math import comb, factorial

from flint import arb, arb_mat, ctx, fmpq

import axiom212_v2 as v2
import capped_mono as cm
from mono_engine import canon

PREC = int(os.environ.get("PREC", "1536"))
ctx.prec = PREC
ZERO = arb(0)
NOMAT = bool(int(os.environ.get("NOMAT", "0")))
TAG = ("_nomat" if NOMAT else "") + (f"_k{os.environ['KVARS']}" if "KVARS" in os.environ else "")


def partitions(n, largest=None):
    if n == 0:
        yield ()
        return
    largest = n if largest is None else min(largest, n)
    for p in range(largest, 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def polymul(p, q):
    r = [ZERO] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r


class Universe:
    """Combinatorial data for B_d that does not depend on the support."""

    def __init__(self, d, k):
        self.d, self.k = d, k
        self.dE = d - d % 2
        self.Wmax = 2 * self.dE
        self.parts = [tuple(2 * x for x in p) for w in range(self.Wmax // 2 + 1) for p in partitions(w)]
        self.idx = {p: i for i, p in enumerate(self.parts)}
        self.wt = [sum(p) for p in self.parts]
        self.alphas = [p for p in self.parts if sum(p) <= self.dE]
        self.basis = [(a, al) for al in self.alphas for a in range(d - sum(al) + 1)]
        self.EmaxV = 2 * d + 2          # largest exponent sum in the marginal keys
        # ordered two-way splits nu = lam + mu
        self.pairs = []
        for inu, nu in enumerate(self.parts):
            items = sorted(Counter(nu).items(), reverse=True)
            for chosen in iproduct(*(range(m + 1) for _, m in items)):
                lam, mu = [], []
                for (v, m), c in zip(items, chosen):
                    lam += [v] * c
                    mu += [v] * (m - c)
                self.pairs.append((self.idx[tuple(lam)], self.idx[tuple(mu)], inu))
        self.byclass = defaultdict(list)
        for i, w in enumerate(self.wt):
            self.byclass[w].append(i)


class Support:
    def __init__(self, K, L, g, tau, caps):
        self.K, self.L, self.g, self.tau = K, fmpq(L), fmpq(g), fmpq(tau)
        self.caps = [fmpq(c) for c in caps]
        self.pr = v2.Params(K=K, L=self.L, g=self.g, tau=self.tau, caps=list(self.caps))
        assert self.L - self.tau <= self.g, "cell-2 rectangle shortcut needs L - tau <= g"
        for m in range(1, len(self.caps) - 1):
            assert self.caps[m] <= self.caps[m + 1] <= self.caps[m] + self.g, "support shape"

    def C(self, m):
        return fmpq(0) if m <= 0 else self.caps[min(m, len(self.caps) - 1)]


def axiom_support(K=None):
    p = v2.Params()
    return Support(K or int(os.environ.get("KVARS", p.K)), p.L, p.g, p.tau, p.caps)


_PRODS = {}


def prodS(U, g):
    key = (id(U), g.p, g.q)
    if key not in _PRODS:
        ga = arb(g)
        gp = [ga ** j for j in range(U.Wmax + 1)]
        S = {e: [comb(e, j) * factorial(j) * gp[e - j] for j in range(e + 1)] for e in range(2, U.Wmax + 1, 2)}
        out = {}
        for i, lam in enumerate(U.parts):
            out[i] = [arb(1)] if not lam else polymul(S[lam[0]], out[U.idx[lam[1:]]])
        _PRODS[key] = out
    return _PRODS[key]


def build_GU(U, pS, s):
    return {i: [cm.nmono(lam, s) * c for c in pS[i]] for i, lam in enumerate(U.parts) if len(lam) <= s}


def build_GY(U, pS, t, m):
    """GY[mu] as coefficient lists (index q = jZ + |lamW|) for t excess and m plain coordinates."""
    GY = {}
    for ilz, ilw, inu in U.pairs:
        lz, lw = U.parts[ilz], U.parts[ilw]
        if len(lz) > t or len(lw) > m:
            continue
        c = cm.nmono(lz, t) * cm.nmono(lw, m)
        for e in lw:
            c *= factorial(e)
        row = GY.get(inu)
        if row is None:
            row = GY[inu] = [ZERO] * (U.wt[inu] + 1)
        sh = U.wt[ilw]
        for j, x in enumerate(pS[ilz]):
            row[j + sh] += c * x
    return GY


# ---------------------------------------------------------------- I
def I_strata(sup):
    return cm.bad_strata(sup.K, sup.L, sup.g, tuple(sup.caps))


def stratum_I(U, sup, s, t):
    """{inu: [value for A = 0 .. 2d - |nu|]}: contribution of stratum (s, t) to int_bad (R-|x|)^A m_nu."""
    ctx.prec = PREC
    k, d = sup.K, U.d
    n = k - s
    R1q = sup.L - (s + t) * sup.g
    a0q = max(fmpq(0), sup.C(s) - s * sup.g)
    R1, b1 = arb(R1q), arb(R1q - a0q)
    mult = comb(k, s) * comb(n, t) * (-1) ** t
    pS = prodS(U, sup.g)
    GU, GY = build_GU(U, pS, s), build_GY(U, pS, t, n - t)
    Amax = 2 * d
    Mmax = U.Wmax + Amax
    R1p = [R1 ** j for j in range(s + U.Wmax + 1)]
    b1p = [b1 ** j for j in range(n + Mmax + s + U.Wmax + 3)]
    # TT[jU][m] = int_0^{b1} (R1 - b)^P b^N db / (P! N!),  P = s - 1 + jU,  N = n + m
    TT = []
    for jU in range(U.Wmax + 1):
        P = s - 1 + jU
        row = []
        for m in range(Mmax + 1):
            N = n + m
            v = ZERO
            for kk in range(P + 1):
                v += (comb(P, kk) * (-1) ** kk) * R1p[P - kk] * b1p[N + kk + 1] / (N + kk + 1)
            row.append(v / (arb(factorial(P)) * arb(factorial(N))))
        TT.append(row)
    # step 1: y[imu][A][jU] = sum_q GY[mu]_q TT[jU][q + A]
    y = {}
    for imu, row in GY.items():
        w = U.wt[imu]
        nJ = U.Wmax - w + 1
        ya = []
        for A in range(Amax - w + 1):
            ya.append([sum((row[q] * TT[jU][q + A] for q in range(w + 1)), ZERO) for jU in range(nJ)])
        y[imu] = ya
    out = {}
    for il, imu, inu in U.pairs:
        gu = GU.get(il)
        ym = y.get(imu)
        if gu is None or ym is None:
            continue
        o = out.get(inu)
        nA = Amax - U.wt[inu] + 1
        if o is None:
            o = out[inu] = [ZERO] * nA
        for A in range(nA):
            yA = ym[A]
            acc = ZERO
            for jU, x in enumerate(gu):
                acc += x * yA[jU]
            o[A] += acc
    for inu, o in out.items():
        for A in range(len(o)):
            o[A] = o[A] * (mult * factorial(A))
    return out


def WI_free(U, sup, inu, A):
    nu, k = U.parts[inu], sup.K
    c = cm.nmono(nu, k)
    for e in nu:
        c *= factorial(e)
    D = k + A + sum(nu)
    return arb(sup.L) ** D * (arb(c * factorial(A)) / arb(factorial(D)))


# ---------------------------------------------------------------- J
def J_strata(sup):
    return v2.J_strata(sup.pr)


def _Fy(Q, e, bound, vp):
    """Arb version of axiom212_v2.Stratum._Fy; vp(j) = v^j."""
    kind, v = bound
    if kind == "y":
        out = [ZERO] * (Q + 1)
        for kk in range(Q + 1):
            out[Q - kk] = (comb(Q, kk) * (-1) ** kk) * vp(kk + e + 1) / (kk + e + 1)
        return out
    out = [ZERO] * (Q + e + 2)
    out[Q + e + 1] = arb(factorial(Q) * factorial(e)) / arb(factorial(Q + e + 1))
    if v != 0:
        for j in range(e + 1):
            out[e - j] -= (comb(e, j) * (-1) ** j) * vp(Q + j + 1) / (Q + j + 1)
    return out


class _Rows:
    """Minimal matrix stand-in (list of rows) with the two operations stratum_J needs."""

    def __init__(self, rows):
        self.rows = rows

    def __add__(self, other):
        return _Rows([[a + b for a, b in zip(r1, r2)] for r1, r2 in zip(self.rows, other.rows)])

    def __getitem__(self, ij):
        return self.rows[ij[0]][ij[1]]


class _Pow:
    def __init__(self, x):
        self.x, self.c = arb(x), [arb(1)]

    def __call__(self, j):
        while len(self.c) <= j:
            self.c.append(self.c[-1] * self.x)
        return self.c[j]


def stratum_J(U, sup, s, t):
    """Tables for stratum (s, t) of the marginal:
    PsiV, Psig: {inu: [E = 0 .. Emax - |nu|]};  X: {inu: {(m1, m2): value}} (cell 2, s > 0);
    Phi0: {inu: [E]} (cell 2, s = 0, without the factor Wc^F)."""
    ctx.prec = PREC
    k, d = sup.K, U.d
    st = v2.Stratum(sup.pr, s, t)
    res = dict(V={}, g={}, X={}, Phi0={})
    if not st.ok:
        return res
    n = (k - 1) - s
    mult = comb(k - 1, s) * comb(n, t) * (-1) ** t
    pS = prodS(U, sup.g)
    GY = build_GY(U, pS, t, n - t)
    Emax = U.EmaxV
    Wm = U.Wmax
    R1 = arb(st.R1)
    Qs = [n - 1 + q for q in range(Wm + 1)]
    invQf = [arb(1) / arb(factorial(Q)) for Q in Qs]

    def addto(tab, inu, vals):
        o = tab.get(inu)
        if o is None:
            tab[inu] = list(vals)
        else:
            for i, v in enumerate(vals):
                o[i] += v

    if s == 0:
        R1p = _Pow(st.R1)
        for cid, (u, slabs) in st.cells.items():
            _, yl, yh = slabs[0]
            pl, ph = _Pow(yl), _Pow(yh)
            phi = []       # phi[q][E]
            for Q in Qs:
                rowE = []
                for E in range(Emax + 1):
                    v = ZERO
                    for kk in range(Q + 1):
                        v += (comb(Q, kk) * (-1) ** kk) * R1p(Q - kk) * (ph(kk + E + 1) - pl(kk + E + 1)) / (kk + E + 1)
                    rowE.append(v)
                phi.append(rowE)
            tab = res["V"] if u == "V" else res["g"] if u == "g" else res["Phi0"]
            for inu, row in GY.items():
                w = U.wt[inu]
                vals = [mult * sum((row[q] * invQf[q] * phi[q][E] for q in range(w + 1)), ZERO)
                        for E in range(Emax - w + 1)]
                addto(tab, inu, vals)
        return res

    GU = build_GU(U, pS, s)
    Ps = [s - 1 + j for j in range(Wm + 1)]
    invPf = [arb(1) / arb(factorial(P)) for P in Ps]
    Imax = Ps[-1]
    Jmax = Qs[-1] + Emax + 1
    R1p = _Pow(st.R1)
    # ---- cells with U = V or U = g: omega(P, Q, E) = sum_slabs sum_i Ec[P][i] It[i]
    Ec = arb_mat(Wm + 1, Imax + 1, [(comb(P, i) * (-1) ** i) * R1p(P - i) if i <= P else ZERO
                                    for P in Ps for i in range(Imax + 1)])
    ncol = (Wm + 1) * (Emax + 1)
    om = {"V": None, "g": None}
    for cid, (u, slabs) in st.cells.items():
        if u == "W":
            continue
        for a0, a1, lo, hi in slabs:
            c0, c1 = _Pow(st.R1 - a1), _Pow(st.R1 - a0)
            mom = [(c1(j + 1) - c0(j + 1)) / (j + 1) for j in range(Imax + Jmax + 2)]
            Hank = arb_mat(Imax + 1, Jmax + 1, [mom[i + j] for i in range(Imax + 1) for j in range(Jmax + 1)])
            plo, phi_ = _Pow(lo[1]), _Pow(hi[1])
            Ip = [[ZERO] * ncol for _ in range(Jmax + 1)]
            for qi, Q in enumerate(Qs):
                for E in range(Emax + 1):
                    col = qi * (Emax + 1) + E
                    for j, v in enumerate(_Fy(Q, E, lo, plo)):
                        Ip[j][col] += v
                    for j, v in enumerate(_Fy(Q, E, hi, phi_)):
                        Ip[j][col] -= v
            if NOMAT:
                # Scalar Arb operations only. arb_mat products go through FLINT's integer matrix/FFT code, and the
                # stock FLINT 3.6.0 wheel has a known signed-FFT bug (fixed in 7ad753d), so the certificate run
                # avoids them and the fast run is compared against it.
                HI = [[sum((Hank[i, j] * Ip[j][c] for j in range(Jmax + 1) if Ip[j][c] != 0), ZERO)
                       for c in range(ncol)] for i in range(Imax + 1)]
                rows = [[sum((Ec[p, i] * HI[i][c] for i in range(p + s)), ZERO) for c in range(ncol)]
                        for p in range(Wm + 1)]
                part = _Rows(rows)
            else:
                IpM = arb_mat(Jmax + 1, ncol, [v for row in Ip for v in row])
                part = Ec * (Hank * IpM)
            om[u] = part if om[u] is None else om[u] + part
    for u in ("V", "g"):
        M = om[u]
        if M is None:
            continue
        # omt[jU][q][E] = omega / (P! Q!)
        omt = [[[M[jU, qi * (Emax + 1) + E] * invPf[jU] * invQf[qi] for E in range(Emax + 1)]
                for qi in range(Wm + 1)] for jU in range(Wm + 1)]
        z = {}
        for imu, row in GY.items():
            w = U.wt[imu]
            nJ = Wm - w + 1
            z[imu] = [[sum((row[q] * omt[jU][q][E] for q in range(w + 1)), ZERO) for jU in range(nJ)]
                      for E in range(Emax - w + 1)]
        tab = res[u]
        for il, imu, inu in U.pairs:
            gu, zm = GU.get(il), z.get(imu)
            if gu is None or zm is None:
                continue
            nE = Emax - U.wt[inu] + 1
            vals = []
            for E in range(nE):
                zE = zm[E]
                acc = ZERO
                for jU, x in enumerate(gu):
                    acc += x * zE[jU]
                vals.append(mult * acc)
            addto(tab, inu, vals)
    # ---- cell 2 (U = W): omega2 = sum_i C(E,i) Aint(P, F+i) Bint(Q, E-i); store X(m1, m2) = sum A w B
    if 2 in st.cells:
        aWp, Wcp, bVp = _Pow(st.aW), _Pow(st.Wc), _Pow(st.bV)
        Mx = Emax
        Aint = [[sum(((comb(m, j) * (-1) ** j) * Wcp(m - j) * aWp(P + j + 1) / (P + j + 1) for j in range(m + 1)), ZERO)
                 * invPf[jU] for m in range(Mx + 1)] for jU, P in enumerate(Ps)]
        Bint = [[bVp(Q + m + 1) * (arb(factorial(m)) / arb(factorial(Q + m + 1))) for m in range(Mx + 1)]
                for Q in Qs]   # the Q! cancels against 1/Q!
        aU = {il: [sum((x * Aint[jU][m1] for jU, x in enumerate(gu)), ZERO) for m1 in range(Mx - U.wt[il] + 1)]
              for il, gu in GU.items()}
        bY = {imu: [sum((row[q] * Bint[q][m2] for q in range(U.wt[imu] + 1)), ZERO)
                    for m2 in range(Mx - U.wt[imu] + 1)] for imu, row in GY.items()}
        X = res["X"]
        for il, imu, inu in U.pairs:
            a, b = aU.get(il), bY.get(imu)
            if a is None or b is None:
                continue
            tot = Mx - U.wt[inu]
            xo = X.get(inu)
            if xo is None:
                xo = X[inu] = {}
            for m1 in range(tot + 1):
                am = mult * a[m1]
                for m2 in range(tot - m1 + 1):
                    key = (m1, m2)
                    v = am * b[m2]
                    xo[key] = xo[key] + v if key in xo else v
    return res


# ---------------------------------------------------------------- serialisation of balls
def ser(x):
    m, r = x.mid(), x.rad()
    return (tuple(map(int, m.man_exp())), tuple(map(int, r.man_exp())))


def de(t):
    (mm, me), (rm, re) = t
    mid = arb(mm) * arb(2) ** me
    if rm == 0:
        return mid
    rad = arb(rm) * arb(2) ** re
    return mid + arb(0, 1) * rad


def ser_I(out):
    return {inu: [ser(v) for v in o] for inu, o in out.items()}


def ser_J(res):
    return dict(V={i: [ser(v) for v in o] for i, o in res["V"].items()},
                g={i: [ser(v) for v in o] for i, o in res["g"].items()},
                Phi0={i: [ser(v) for v in o] for i, o in res["Phi0"].items()},
                X={i: {key: ser(v) for key, v in o.items()} for i, o in res["X"].items()})


class Tables:
    """Accumulated sums over strata."""

    def __init__(self):
        self.I, self.V, self.g, self.Phi0, self.X = {}, {}, {}, {}, {}
        self.doneI, self.doneJ = set(), set()

    @staticmethod
    def _addlist(tab, i, vals):
        o = tab.get(i)
        if o is None:
            tab[i] = [de(v) for v in vals]
        else:
            for j, v in enumerate(vals):
                o[j] += de(v)

    def add_I(self, st, out):
        for i, vals in out.items():
            self._addlist(self.I, i, vals)
        self.doneI.add(st)

    def add_J(self, st, res):
        for name in ("V", "g", "Phi0"):
            tab = getattr(self, name)
            for i, vals in res[name].items():
                self._addlist(tab, i, vals)
        for i, o in res["X"].items():
            xo = self.X.setdefault(i, {})
            for key, v in o.items():
                xo[key] = xo[key] + de(v) if key in xo else de(v)
        self.doneJ.add(st)

    def dump(self, path):
        data = dict(I={i: [ser(v) for v in o] for i, o in self.I.items()},
                    V={i: [ser(v) for v in o] for i, o in self.V.items()},
                    g={i: [ser(v) for v in o] for i, o in self.g.items()},
                    Phi0={i: [ser(v) for v in o] for i, o in self.Phi0.items()},
                    X={i: {key: ser(v) for key, v in o.items()} for i, o in self.X.items()},
                    doneI=self.doneI, doneJ=self.doneJ, prec=PREC)
        tmp = path + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, path)

    @classmethod
    def load(cls, path):
        T = cls()
        data = pickle.load(open(path, "rb"))
        for name in ("I", "V", "g", "Phi0"):
            setattr(T, name, {i: [de(v) for v in o] for i, o in data[name].items()})
        T.X = {i: {key: de(v) for key, v in o.items()} for i, o in data["X"].items()}
        T.doneI, T.doneJ = data["doneI"], data["doneJ"]
        return T


# ---------------------------------------------------------------- Gram assembly
def assemble(U, sup, T):
    """Gram matrices (lists of lists of arb) of I_T and J_T in the basis U.basis."""
    ctx.prec = PREC
    d, k, g = U.d, sup.K, arb(sup.g)
    nb = len(U.basis)
    # WI = free - bad
    WI = {}
    for inu in range(len(U.parts)):
        if U.wt[inu] > U.Wmax:
            continue
        bad = T.I.get(inu)
        WI[inu] = [WI_free(U, sup, inu, A) - (bad[A] if bad else ZERO) for A in range(2 * d - U.wt[inu] + 1)]
    Wc0 = _Pow(sup.C(1))

    def PsiW(inu, E, F):
        v = ZERO
        xo = T.X.get(inu)
        if xo:
            for i in range(E + 1):
                x = xo.get((F + i, E - i))
                if x is not None:
                    v += comb(E, i) * x
        p0 = T.Phi0.get(inu)
        if p0:
            v += Wc0(F) * p0[E]
        return v

    PV, Pg = T.V, T.g
    cV = lambda a, e: arb(factorial(a) * factorial(e)) / arb(factorial(a + e + 1))
    gp = _Pow(sup.g)
    gam = lambda a, e, j: (comb(a, j) * (-1) ** j) * gp(e + j + 1) / (e + j + 1)
    dlt = lambda a, e, j: arb(comb(a, j) * (-1) ** j) / (e + j + 1)
    cache = {}

    def lam(inu, a, e, a2, e2):
        key = (inu, a, e, a2, e2) if (a, e) <= (a2, e2) else (inu, a2, e2, a, e)
        v = cache.get(key)
        if v is None:
            v = ZERO
            pv = PV.get(inu)
            if pv:
                v += cV(a, e) * cV(a2, e2) * pv[a + e + a2 + e2 + 2]
            pg = Pg.get(inu)
            if pg:
                for j in range(a + 1):
                    gj = gam(a, e, j)
                    for j2 in range(a2 + 1):
                        v += gj * gam(a2, e2, j2) * pg[a - j + a2 - j2]
            for j in range(a + 1):
                dj = dlt(a, e, j)
                for j2 in range(a2 + 1):
                    v += dj * dlt(a2, e2, j2) * PsiW(inu, a - j + a2 - j2, e + j + e2 + j2 + 2)
            cache[key] = v
        return v

    pos = {b: i for i, b in enumerate(U.basis)}
    GI = [[ZERO] * nb for _ in range(nb)]
    GJ = [[ZERO] * nb for _ in range(nb)]
    al = U.alphas
    for i1, a1 in enumerate(al):
        n1 = d - sum(a1)
        for a2 in al[: i1 + 1]:
            n2 = d - sum(a2)
            for nu, c in cm.mmul(a1, a2).items():
                wi = WI[U.idx[nu]]
                for x in range(n1 + 1):
                    r = pos[(x, a1)]
                    for y in range(n2 + 1):
                        GI[r][pos[(y, a2)]] += c * wi[x + y]
            for e1 in [0] + sorted(set(a1)):
                b1 = cm.remove_one(a1, e1) if e1 else a1
                for e2 in [0] + sorted(set(a2)):
                    b2 = cm.remove_one(a2, e2) if e2 else a2
                    for nu, c in cm.mmul(b1, b2).items():
                        inu = U.idx[nu]
                        for x in range(n1 + 1):
                            r = pos[(x, a1)]
                            for y in range(n2 + 1):
                                GJ[r][pos[(y, a2)]] += c * lam(inu, x, e1, y, e2)
    for i1, a1 in enumerate(al):          # fill the other triangle of the alpha blocks
        for a2 in al[:i1]:
            for x in range(d - sum(a1) + 1):
                r = pos[(x, a1)]
                for y in range(d - sum(a2) + 1):
                    c = pos[(y, a2)]
                    GI[c][r] = GI[r][c]
                    GJ[c][r] = GJ[r][c]
    return GI, GJ


def top_ratio(GI, GJ, k, squarings=None, verbose=False):
    """Largest generalized eigenvalue of k J x = lambda I x, as a rigorous lower bound.

    C = I^-1 (k J) after diagonal scaling; the dominant eigenvector is extracted by repeated squaring of C
    (all eigenvalues are real and nonnegative). The returned ball is the Rayleigh quotient of a vector with
    exact dyadic entries, so its lower end is a valid lower bound whatever the convergence."""
    nb = len(GI)
    sc = [(arb(1) / GI[i][i].sqrt()).mid() for i in range(nb)]
    A = arb_mat(nb, nb, [GI[i][j] * sc[i] * sc[j] for i in range(nb) for j in range(nb)])
    B = arb_mat(nb, nb, [k * GJ[i][j] * sc[i] * sc[j] for i in range(nb) for j in range(nb)])
    C = A.solve(B)
    x = arb_mat(nb, 1, [arb(1)] * nb)
    last, q = None, None
    for it in range(squarings or 60):
        y = C * x
        nrm = max(abs(y[i, 0].mid()) for i in range(nb))
        x = arb_mat(nb, 1, [(y[i, 0] / nrm).mid() for i in range(nb)])
        q = (x.transpose() * B * x)[0, 0] / (x.transpose() * A * x)[0, 0]
        if verbose:
            print(f"   squaring {it}: ratio = {q.str(30)}", flush=True)
        qm = q.mid()
        if last is not None and abs(qm - last) < abs(qm) * arb(10) ** -40:
            break
        last = qm
        C = C * C
        top = max(abs(C[i, i].mid()) for i in range(nb))
        C = C * (arb(1) / top).mid()
        C = arb_mat(nb, nb, [C[i, j].mid() for i in range(nb) for j in range(nb)])   # any matrix will do
    xs = arb_mat(nb, 1, [x[i, 0] * sc[i] for i in range(nb)])
    return q, xs


# ---------------------------------------------------------------- driver
_CTX = {}


def _job(args):
    kind, s, t = args
    U, sup = _CTX["U"], _CTX["sup"]
    t0 = time.time()
    if kind == "I":
        return kind, s, t, ser_I(stratum_I(U, sup, s, t)), time.time() - t0
    return kind, s, t, ser_J(stratum_J(U, sup, s, t)), time.time() - t0


def run(d, sup=None, workers=2, path=None, log=print, save_every=600):
    from multiprocessing import get_context
    sup = sup or axiom_support()
    U = Universe(d, sup.K)
    path = path or f"gramtab_d{d}{TAG}.pkl"
    T = Tables.load(path) if os.path.exists(path) else Tables()
    jobs = [("I", s, t) for s, t in I_strata(sup) if (s, t) not in T.doneI]
    jobs += [("J", s, t) for s, t in J_strata(sup) if (s, t) not in T.doneJ]
    log(f"d={d}: dim {len(U.basis)}, {len(U.parts)} nu, {len(U.pairs)} splits, {len(jobs)} strata to do "
        f"({len(T.doneI)} I and {len(T.doneJ)} J done), prec {PREC}")
    _CTX.update(U=U, sup=sup)
    prodS(U, sup.g)
    t0 = last = time.time()
    if jobs:
        # No maxtasksperchild: re-forking a worker while the main thread is inside print() copies a locked
        # stdout buffer into the child, which then hangs at its first flush (seen at d = 12).
        with get_context("fork").Pool(workers) as pool:
            for n_done, (kind, s, t, data, dt) in enumerate(pool.imap_unordered(_job, jobs), 1):
                (T.add_I if kind == "I" else T.add_J)((s, t), data)
                log(f"  {kind} stratum s={s} t={t}: {dt:.0f}s   [{n_done}/{len(jobs)}, {time.time()-t0:.0f}s]")
                if time.time() - last > save_every:
                    T.dump(path); last = time.time()
        T.dump(path)
    return U, sup, T


if __name__ == "__main__":
    d = int(sys.argv[1])
    workers = int(os.environ.get("WORKERS", "2"))
    t0 = time.time()
    lg = lambda *a: os.write(1, (f"[{time.time()-t0:7.0f}s] " + " ".join(map(str, a)) + "\n").encode())
    U, sup, T = run(d, workers=workers, log=lg)
    GI, GJ = assemble(U, sup, T)
    lg("Gram matrices assembled")
    with open(f"gramtab_d{d}{TAG}_gram.pkl", "wb") as f:
        pickle.dump(dict(basis=U.basis, GI=[[ser(v) for v in r] for r in GI], GJ=[[ser(v) for v in r] for r in GJ],
                         K=sup.K, prec=PREC), f, protocol=pickle.HIGHEST_PROTOCOL)
    q, x = top_ratio(GI, GJ, sup.K, verbose=True)
    lg(f"d={d} dim={len(U.basis)}: K*J/I = {q.str(40)}  "
       f"> 4: {bool(q > 4)}")
    with open(f"gramtab_d{d}{TAG}_vector.pkl", "wb") as f:
        pickle.dump(dict(basis=U.basis, x=[ser(x[i, 0]) for i in range(len(U.basis))], ratio=ser(q)), f)

"""Kintali capped J_T, cells with U = V (cells 1, 3) and U = g (cell 4), exact; cell 2 (U = B - rho1) separately.

J_T = int_{x^ in T^} ( int_0^{U(x^)} F(x^,u) du )^2,  U = min{V, max{g, B - rho1}},  V = R - |x^|.
G on a cell = sum_beta m_beta(x^) * sum_h gamma[beta][h] V^h ; square first -> sum_nu m_nu * poly_nu(V);
per x^-stratum: fibre weights w_nu(P,Q) from 3-way splits (capped_mono), 2-D cell integrals from axiom212_v2.Stratum.
"""
import sys, time, pickle, os
from math import comb, factorial
from collections import defaultdict
from flint import fmpq
import capped_mono as cm, axiom212_v2 as v2
from mono_engine import remove_one
from kintali_capped import W, k, R, g, B, r, c as zc


def G_terms(kind):
    out = defaultdict(lambda: defaultdict(lambda: fmpq(0)))
    for a, al, cf in W:
        for e in [0] + sorted(set(al)):
            beta = remove_one(al, e) if e else al
            if len(beta) > k - 1:
                continue
            if kind == "V":
                out[beta][a + e + 1] += cf * fmpq(factorial(a) * factorial(e), factorial(a + e + 1))
            else:  # U = g:  sum_j C(a,j)(-1)^j V^(a-j) g^(e+j+1)/(e+j+1)
                for j in range(a + 1):
                    out[beta][a - j] += cf * comb(a, j) * (-1) ** j * g ** (e + j + 1) / (e + j + 1)
    return {b: {h: v for h, v in d.items() if v != 0} for b, d in out.items()}


def square_G(G):
    bs = sorted(G)
    sq = defaultdict(lambda: defaultdict(lambda: fmpq(0)))
    for i, b1 in enumerate(bs):
        for b2 in bs[: i + 1]:
            conv = defaultdict(lambda: fmpq(0))
            for p, x in G[b1].items():
                for q, y in G[b2].items():
                    conv[p + q] += x * y
            sym = 1 if b1 == b2 else 2
            for nu, cn in cm.mmul(b1, b2).items():
                if len(nu) > k - 1:
                    continue
                d = sq[nu]
                for E, v in conv.items():
                    d[E] += sym * cn * v
    return {nu: {E: v for E, v in d.items() if v != 0} for nu, d in sq.items()}


PR = v2.Params(K=k, L=R, g=g, tau=r, caps=[fmpq(0)] + [B] * 80)


def fibre_weights(nu, s, t, n):
    Wt = defaultdict(lambda: fmpq(0))
    for nU, nZ, nW in cm.splits3(nu, s, t, n - t):
        if (not s and nU) or (not t and nZ):
            continue
        cnt = cm.nmono(nU, s) * cm.nmono(nZ, t) * cm.nmono(nW, n - t) if s else cm.nmono(nZ, t) * cm.nmono(nW, n - t)
        if not cnt:
            continue
        fw = 1
        for e in nW:
            fw *= factorial(e)
        pu = cm.shift_profile(nU, g) if s else (fmpq(1),)
        pz = cm.shift_profile(nZ, g) if t else (fmpq(1),)
        dW = sum(nW)
        for jU, xu in enumerate(pu):
            if xu == 0:
                continue
            for jZ, xz in enumerate(pz):
                if xz == 0:
                    continue
                P = (s - 1 + jU) if s else -1
                Q = n - 1 + jZ + dW
                Wt[(P, Q)] += cnt * fw * xu * xz / (factorial(P) if P >= 0 else 1) / factorial(Q)
    return Wt


def J_stratum(args):
    sqV, sqg, s, t = args
    st = v2.Stratum(PR, s, t)
    if not st.ok:
        return s, t, (0, 1), {}
    n = (k - 1) - s
    mult = comb(k - 1, s) * comb(n, t) * (-1) ** t
    tot = fmpq(0)
    percell = {}
    cells = {cid: u for cid, (u, _) in st.cells.items()}
    for nu in set(sqV) | set(sqg):
        wts = None
        for cid, u in cells.items():
            if u == "W":
                continue
            poly = (sqV if u == "V" else sqg).get(nu)
            if not poly:
                continue
            if wts is None:
                wts = fibre_weights(nu, s, t, n)
            for (P, Q), w in wts.items():
                for E, cE in poly.items():
                    val = w * cE * st.omega(cid, P, Q, E, 0)
                    tot += val
                    percell[cid] = percell.get(cid, fmpq(0)) + val
    tot *= mult
    return s, t, (int(tot.p), int(tot.q)), {c: (int((mult * v).p), int((mult * v).q)) for c, v in percell.items()}


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] in ("run134",):
    from multiprocessing import get_context
    t0 = time.time()
    if os.path.exists("kintali_Gsq.pkl"):
        sqV, sqg = pickle.load(open("kintali_Gsq.pkl", "rb"))
    else:
        sqV = square_G(G_terms("V")); print(f"square V-type: {len(sqV)} nu [{time.time()-t0:.0f}s]", flush=True)
        sqg = square_G(G_terms("g")); print(f"square g-type: {len(sqg)} nu [{time.time()-t0:.0f}s]", flush=True)
        pickle.dump((sqV, sqg), open("kintali_Gsq.pkl", "wb"))
    strata = v2.J_strata(PR)
    print(f"{len(strata)} marginal strata", flush=True)
    if len(sys.argv) > 1 and sys.argv[1] == "time1":
        tt = time.time(); print(J_stratum((sqV, sqg, *strata[int(sys.argv[2])]))[:2], f"{time.time()-tt:.1f}s"); sys.exit()
    J = fmpq(0); cellsum = defaultdict(lambda: fmpq(0))
    with get_context("fork").Pool(14) as pool:
        for s, t, (p, q), pc in pool.imap_unordered(J_stratum, [(sqV, sqg, s, t) for s, t in strata]):
            J += fmpq(p, q)
            for cid, (pp, qq) in pc.items():
                cellsum[cid] += fmpq(pp, qq)
            print(f"   J stratum s={s} t={t} [{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(dict(J_cells134=J, percell={c: v for c, v in cellsum.items()}), open("kintali_J134.pkl", "wb"))
    print("J cells 1,3,4 =", float(J), {c: float(v) for c, v in cellsum.items()})


# ---------------- cell 2 (U = W = B - rho1): bivariate in (V, W) ----------------
def G_terms_W():
    out = defaultdict(lambda: defaultdict(lambda: fmpq(0)))
    for a, al, cf in W:
        for e in [0] + sorted(set(al)):
            beta = remove_one(al, e) if e else al
            if len(beta) > k - 1:
                continue
            for j in range(a + 1):   # V^(a-j) W^(e+j+1)
                out[beta][(a - j, e + j + 1)] += fmpq(cf * comb(a, j) * (-1) ** j, e + j + 1)
    return {b: {h: v for h, v in d.items() if v != 0} for b, d in out.items()}


def square_GW(G):
    bs = sorted(G)
    sq = defaultdict(lambda: defaultdict(lambda: fmpq(0)))
    for i, b1 in enumerate(bs):
        for b2 in bs[: i + 1]:
            conv = defaultdict(lambda: fmpq(0))
            for (p1, q1), x in G[b1].items():
                for (p2, q2), y in G[b2].items():
                    conv[(p1 + p2, q1 + q2)] += x * y
            sym = 1 if b1 == b2 else 2
            for nu, cn in cm.mmul(b1, b2).items():
                if len(nu) > k - 1:
                    continue
                d = sq[nu]
                for EF, v in conv.items():
                    d[EF] += sym * cn * v
    return {nu: {EF: v for EF, v in d.items() if v != 0} for nu, d in sq.items()}


def J_stratum_cell2(args):
    sqW, s, t = args
    st = v2.Stratum(PR, s, t)
    if not st.ok or 2 not in st.cells:
        return s, t, (0, 1)
    n = (k - 1) - s
    mult = comb(k - 1, s) * comb(n, t) * (-1) ** t
    tot = fmpq(0)
    for nu, poly in sqW.items():
        wts = fibre_weights(nu, s, t, n)
        for (P, Q), w in wts.items():
            for (E, Fx), cE in poly.items():
                tot += w * cE * st.omega(2, P, Q, E, Fx)
    tot *= mult
    return s, t, (int(tot.p), int(tot.q))


# ---------------- optimised full stratum (all four cells) ----------------
def _incl(P, m, X):
    """int_0^X a^P (X' - a)^m da is not needed; we need A[P][m] = int_0^{aW} a^P (Wc - a)^m da."""
    raise NotImplementedError


class Cell2Tables:
    def __init__(self, st):
        self.st = st
        self.A, self.Bt = {}, {}

    def Aint(self, P, m):  # int_0^{aW} a^P (Wc - a)^m da
        key = (P, m)
        if key not in self.A:
            st = self.st
            aW, Wc = st.aW, st.Wc
            # expand (Wc - a)^m
            self.A[key] = sum((comb(m, j) * Wc ** (m - j) * (-1) ** j * aW ** (P + j + 1) / (P + j + 1) for j in range(m + 1)), fmpq(0))
        return self.A[key]

    def Bint(self, Q, m):  # int_0^{bV} b^Q (bV - b)^m db = bV^(Q+m+1) Q! m!/(Q+m+1)!
        key = (Q, m)
        if key not in self.Bt:
            bV = self.st.bV
            self.Bt[key] = bV ** (Q + m + 1) * fmpq(factorial(Q) * factorial(m), factorial(Q + m + 1))
        return self.Bt[key]

    def omega2(self, P, Q, E, F):
        # V = W + Delta, Delta = bV - b ; V^E W^F = sum_i C(E,i) W^(F+i) Delta^(E-i)
        return sum((comb(E, i) * self.Aint(P, F + i) * self.Bint(Q, E - i) for i in range(E + 1)), fmpq(0))


def J_stratum_all(args):
    sqV, sqg, sqW, s, t = args
    st = v2.Stratum(PR, s, t)
    if not st.ok:
        return s, t, (0, 1), {}
    n = (k - 1) - s
    mult = comb(k - 1, s) * comb(n, t) * (-1) ** t
    cells = {cid: u for cid, (u, _) in st.cells.items()}
    T = {cid: defaultdict(lambda: fmpq(0)) for cid in cells}
    for nu in set(sqV) | set(sqg) | set(sqW):
        wts = None
        for cid, u in cells.items():
            poly = {"V": sqV, "g": sqg, "W": sqW}[u].get(nu)
            if not poly:
                continue
            if wts is None:
                wts = fibre_weights(nu, s, t, n)
            Tc = T[cid]
            for (P, Q), w in wts.items():
                for key, cE in poly.items():
                    Tc[(P, Q, key)] += w * cE
    percell = {}
    if 2 in cells and s > 0:
        assert R - r <= g, 'cell-2 rectangle shortcut needs R - r <= g'
    c2 = Cell2Tables(st) if 2 in cells and s > 0 else None
    for cid, Tc in T.items():
        acc = fmpq(0)
        for (P, Q, key), v in Tc.items():
            if v == 0:
                continue
            if cells[cid] == "W":
                E, F = key
                om = c2.omega2(P, Q, E, F) if c2 is not None else st.omega(cid, P, Q, E, F)
            else:
                om = st.omega(cid, P, Q, key, 0)
            acc += v * om
        percell[cid] = mult * acc
    tot = sum(percell.values(), fmpq(0))
    return s, t, (int(tot.p), int(tot.q)), {c: (int(v.p), int(v.q)) for c, v in percell.items()}


def main_all():
    from multiprocessing import get_context
    t0 = time.time()
    if os.path.exists("kintali_Gsq_all.pkl"):
        sqV, sqg, sqW = pickle.load(open("kintali_Gsq_all.pkl", "rb"))
    else:
        sqV = square_G(G_terms("V")); print(f"square V: {len(sqV)} nu [{time.time()-t0:.0f}s]", flush=True)
        sqg = square_G(G_terms("g")); print(f"square g: {len(sqg)} nu [{time.time()-t0:.0f}s]", flush=True)
        sqW = square_GW(G_terms_W()); print(f"square W: {len(sqW)} nu, {sum(len(v) for v in sqW.values())} terms [{time.time()-t0:.0f}s]", flush=True)
        pickle.dump((sqV, sqg, sqW), open("kintali_Gsq_all.pkl", "wb"))
    strata = v2.J_strata(PR)
    print(f"{len(strata)} marginal strata", flush=True)
    if len(sys.argv) > 2 and sys.argv[1] == "time1":
        tt = time.time(); res = J_stratum_all((sqV, sqg, sqW, *strata[int(sys.argv[2])]))
        print(res[:2], {c: float(fmpq(*v)) for c, v in res[3].items()}, f"{time.time()-tt:.1f}s"); return
    J = fmpq(0); cellsum = defaultdict(lambda: fmpq(0))
    with get_context("fork").Pool(int(os.environ.get("WORKERS", "14"))) as pool:
        for s, t, (p, q), pc in pool.imap_unordered(J_stratum_all, [(sqV, sqg, sqW, s, t) for s, t in strata]):
            J += fmpq(p, q)
            for cid, (pp, qq) in pc.items():
                cellsum[cid] += fmpq(pp, qq)
            print(f"   J stratum s={s} t={t} [{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(dict(J_T=J, percell=dict(cellsum)), open("kintali_J_all.pkl", "wb"))
    print("J_T =", float(J), {c: float(v) for c, v in cellsum.items()})


if __name__ == "__main__" and (len(sys.argv) == 1 or sys.argv[1] in ("time1", "all")):
    main_all()

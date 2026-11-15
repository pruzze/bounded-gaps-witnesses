"""Exact I, J for witnesses written in the monomial-symmetric basis, on an (enlarged) simplex.

F(x) = sum_{(a, alpha)} c_{a,alpha} (R - |x|_1)^a m_alpha(x),   supported on R * simplex_k.
  I(F) = int_{|x| <= R} F^2
  J(F) = int_{|x^| <= r} ( int_0^{R - |x^|} F(x^, u) du )^2 dx^       (x^ in k-1 variables)
Polymath8b's M_{k,eps}:  R = 1 + eps, r = 1 - eps, criterion k J / I > 4.

Key quantity (same idea as the 'factorial moment' in Axiom's PrimeGapsLib, implemented independently):
  FM_n(alpha, beta) = sum over monomials x^e of m_alpha * m_beta in n variables (with multiplicity) of prod e_i!
so that  int_{|x|<=R} (R-|x|)^A m_alpha m_beta = R^{n+A+d} A! FM_n(alpha,beta) / (n+A+d)!,  d = |alpha|+|beta|.
FM_n = sum_r C(n, r) N_r(alpha, beta), N_r = ordered fillings of r active coordinates (memoised recursion).
"""
from functools import lru_cache
from math import comb, factorial
from collections import defaultdict
from flint import fmpq, fmpq_mat


def canon(t):
    return tuple(sorted(t, reverse=True))


def remove_one(t, v):
    l = list(t); l.remove(v); return tuple(l)


@lru_cache(maxsize=None)
def nil(alpha, beta):
    """{r: N_r(alpha, beta)} as a tuple of (r, count)."""
    if not alpha and not beta:
        return ((0, 1),)
    out = defaultdict(int)
    xs = [0] + sorted(set(alpha)); ys = [0] + sorted(set(beta))
    for x in xs:
        a2 = remove_one(alpha, x) if x else alpha
        for y in ys:
            if x == 0 and y == 0:
                continue
            b2 = remove_one(beta, y) if y else beta
            w = factorial(x + y)
            for r, c in (nil(a2, b2) if a2 <= b2 else nil(b2, a2)):
                out[r + 1] += w * c
    return tuple(sorted(out.items()))


def FM(n, alpha, beta):
    a, b = (alpha, beta) if alpha <= beta else (beta, alpha)
    return sum(comb(n, r) * c for r, c in nil(a, b))


def load_witness(path):
    """Terms [a, alpha, c] of F = sum c (R - |x|_1)^a m_alpha(x). Accepts a bare list of terms or a witness file
    of this repository (a JSON object whose "terms" carry the integer coefficients as strings)."""
    import json
    data = json.load(open(path))
    if isinstance(data, dict):
        data = data["terms"]
    return [(a, canon(al), c if isinstance(c, int) else int(c)) for a, al, c in data]


def _blocks(vecs):
    """vecs: {key alpha: {p: coef}} -> grouped by degree |alpha| -> (list of alphas, list of p's, fmpq_mat coef[p_idx][alpha_idx])."""
    bydeg = defaultdict(list)
    for al in vecs:
        bydeg[sum(al)].append(al)
    out = {}
    for d, als in bydeg.items():
        ps = sorted({p for al in als for p in vecs[al]})
        M = fmpq_mat(len(ps), len(als))
        for j, al in enumerate(als):
            for i, p in enumerate(ps):
                c = vecs[al].get(p)
                if c:
                    M[i, j] = c
        out[d] = (als, ps, M)
    return out


def quadratic(vecs, n, kernel, verbose=False):
    """sum_{al, al'} FM_n(al, al') * sum_{p,p'} v_al[p] v_al'[p'] kernel(p, p', |al|+|al'|)."""
    import time
    t0 = time.time()
    B = _blocks(vecs)
    degs = sorted(B)
    total = fmpq(0)
    for i1, d1 in enumerate(degs):
        als1, ps1, M1 = B[d1]
        for d2 in degs[i1:]:
            als2, ps2, M2 = B[d2]
            d = d1 + d2
            K = fmpq_mat(len(ps1), len(ps2))
            for i, p in enumerate(ps1):
                for j, q in enumerate(ps2):
                    K[i, j] = kernel(p, q, d)
            H = M1.transpose() * K * M2  # H[al1, al2]
            sub = fmpq(0)
            for a, al1 in enumerate(als1):
                for b, al2 in enumerate(als2):
                    h = H[a, b]
                    if h != 0:
                        sub += h * FM(n, al1, al2)
            total += sub if d1 == d2 else 2 * sub
        if verbose:
            print(f"   deg {d1}: {time.time()-t0:.1f}s, nil cache {nil.cache_info().currsize}", flush=True)
    return total


def IJ(witness, k, R, r, verbose=False):
    # I: vectors over a for each alpha
    vI = defaultdict(dict)
    for a, al, c in witness:
        vI[al][a] = vI[al].get(a, 0) + c
    kerI = lambda a, b, d: R ** (k + a + b + d) * fmpq(factorial(a + b), factorial(k + a + b + d))
    I = quadratic(vI, k, kerI, verbose)
    # J: G(x^) = sum g (R - S)^p m_beta(x^),  p = a + e + 1, beta = alpha minus one part e (or e = 0)
    vJ = defaultdict(lambda: defaultdict(lambda: fmpq(0)))
    for a, al, c in witness:
        for e in [0] + sorted(set(al)):
            beta = remove_one(al, e) if e else al
            if len(beta) > k - 1:
                continue
            vJ[beta][a + e + 1] += c * fmpq(factorial(a) * factorial(e), factorial(a + e + 1))
    vJ = {b: {p: v for p, v in d.items() if v != 0} for b, d in vJ.items()}
    n = k - 1
    Rr = R - r

    @lru_cache(maxsize=None)
    def kerJ(p, q, d):
        P = p + q
        return sum((comb(P, i) * Rr ** (P - i) * r ** (n + i + d) * fmpq(factorial(i), factorial(n + i + d)) for i in range(P + 1)), fmpq(0))
    J = quadratic(vJ, n, kerJ, verbose)
    return I, J

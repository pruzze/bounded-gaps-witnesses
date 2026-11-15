"""Exact variational quantities for the Maynard-Tao sieve on the simplex.

Support: R_k = {t in [0,inf)^k : t_1 + ... + t_k <= 1}.
Basis functions: b = (1 - P1)^a * P_lam, P_lam = prod_j P_{lam_j}, P_r = sum_i t_i^r.

  I(F)  = int_{R_k} F^2
  J(F)  = int_{R_{k-1}} ( int_0^{1 - P1(t')} F(t', u) du )^2 dt'
  M_k  >= k * J(F) / I(F)      (Maynard's criterion: need > 4 for DHL[k,2] via Bombieri-Vinogradov)

Everything reduces to G(k, A, Lam) = int_{R_k} (1 - P1)^A P_Lam, evaluated exactly by
expanding P_Lam into monomial symmetric functions m_mu and the Dirichlet identity
  int_{R_k} (1-P1)^A prod t_i^{e_i} = A! prod e_i! / (k + A + sum e_i)!.
"""
from fractions import Fraction
from functools import lru_cache
from math import factorial
from itertools import combinations
from collections import Counter


def _canon(parts):
    return tuple(sorted((p for p in parts if p), reverse=True))


@lru_cache(maxsize=None)
def p_to_m(lam, k):
    """Power-sum product P_lam in the monomial symmetric basis (k variables): {mu: coeff}."""
    res = {(): 1}
    for r in lam:
        nxt = Counter()
        for mu, c in res.items():
            mu_l = list(mu)
            # add r to an existing part (each distinct value once; coefficient = multiplicity in result)
            for v in set(mu_l):
                i = mu_l.index(v)
                nu = _canon(mu_l[:i] + [v + r] + mu_l[i + 1:])
                nxt[nu] += c * Counter(nu)[v + r]
            if len(mu_l) < k:  # new part
                nu = _canon(mu_l + [r])
                nxt[nu] += c * Counter(nu)[r]
        res = dict(nxt)
    return res


@lru_cache(maxsize=None)
def n_monomials(mu, k):
    out = factorial(k) // factorial(k - len(mu))
    for m in Counter(mu).values():
        out //= factorial(m)
    return out


@lru_cache(maxsize=None)
def G(k, A, lam):
    lam = _canon(lam)
    total = Fraction(0)
    for mu, c in p_to_m(lam, k).items():
        num = factorial(A)
        for e in mu:
            num *= factorial(e)
        total += Fraction(c * n_monomials(mu, k) * num, factorial(k + A + sum(mu)))
    return total


def inner_integral(a, lam):
    """int_0^{T} (T-u)^a prod_j (Q_{lam_j} + u^{lam_j}) du, T = 1 - S, as {(e, nu): coeff} meaning coeff*(1-S)^e*Q_nu."""
    out = Counter()
    idx = range(len(lam))
    for r in range(len(lam) + 1):
        for A in combinations(idx, r):
            m = sum(lam[i] for i in A)
            rest = _canon([lam[i] for i in idx if i not in A])
            out[(a + m + 1, rest)] += Fraction(factorial(a) * factorial(m), factorial(a + m + 1))
    return out


def gram_matrices(k, basis):
    """basis: list of (a, lam). Returns exact I, J (lists of lists of Fractions)."""
    n = len(basis)
    I = [[G(k, basis[i][0] + basis[j][0], basis[i][1] + basis[j][1]) for j in range(n)] for i in range(n)]
    inner = [inner_integral(a, lam) for a, lam in basis]
    J = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            s = Fraction(0)
            for (e1, n1), c1 in inner[i].items():
                for (e2, n2), c2 in inner[j].items():
                    s += c1 * c2 * G(k - 1, e1 + e2, n1 + n2)
            J[i][j] = J[j][i] = s
    return I, J


def top_ratio(k, I, J, dps=60):
    """Largest generalized eigenvalue of k*J v = lambda I v (high precision), its eigenvector,
    and an exact rational lower bound from a rationalised eigenvector."""
    import mpmath as mp
    mp.mp.dps = dps
    n = len(I)
    # diagonal scaling D = diag(I_ii^-1/2): entries are ~1e-200, and mpmath's Cholesky uses an absolute tolerance
    d = [1 / mp.sqrt(mp.mpf(I[t][t].numerator) / I[t][t].denominator) for t in range(n)]
    Im = mp.matrix([[d[a] * d[b] * mp.mpf(I[a][b].numerator) / I[a][b].denominator for b in range(n)] for a in range(n)])
    Jm = mp.matrix([[k * d[a] * d[b] * mp.mpf(J[a][b].numerator) / J[a][b].denominator for b in range(n)] for a in range(n)])
    L = mp.cholesky(Im)
    Li = mp.inverse(L)
    S = Li * Jm * Li.T
    E, Q = mp.eigsy(S)
    i = max(range(n), key=lambda t: E[t])
    w = Li.T * Q[:, i]
    v = [d[t] * w[t] for t in range(n)]
    vmax = max(abs(v[t]) for t in range(n))
    vq = [Fraction(int(mp.nint(v[t] / vmax * 10**30)), 10**30) for t in range(n)]
    num = sum(vq[a] * k * J[a][b] * vq[b] for a in range(n) for b in range(n))
    den = sum(vq[a] * I[a][b] * vq[b] for a in range(n) for b in range(n))
    return E[i], num / den, vq


if __name__ == "__main__":
    import sys, time
    # Validation 1: Maynard (2013), k = 105, basis (1-P1)^b P2^c with b + 2c <= 11 -> lambda ~ 4.0020697
    t0 = time.time()
    k = 105
    basis = [(b, (2,) * c) for b in range(12) for c in range(6) if b + 2 * c <= 11]
    print("basis size", len(basis))
    I, J = gram_matrices(k, basis)
    lam, exact, _ = top_ratio(k, I, J)
    print(f"k={k}: top eigenvalue {mp_str(lam) if False else lam}")
    print(f"exact rational lower bound: {float(exact):.10f}  (> 4: {exact > 4})  [{time.time()-t0:.1f}s]")

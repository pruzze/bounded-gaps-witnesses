"""Capped-support integrals for monomial-basis witnesses (Stadlmann-type one-band support, uniform cap).

Support (in the witness's own coordinates):  T = { x in [0,inf)^k : |x|_1 <= R,  sum_{x_i > g} x_i <= B  (if any x_i > g) }.
I_T(F) = I_free(F) - I_bad(F),   bad = { rough mass > B }.

Method:
  1. square first:  F^2 = sum_nu sum_A C[nu][A] (R - |x|)^A m_nu(x)   (exact integers; products m_a m_b by an
     ordered-filling recursion, c_nu = F_nu * prod(mult!) / len(nu)!);
  2. stratify the bad set by the exact rough set (s coordinates, shifted by g) and, for the box x_i <= g on the
     others, inclusion-exclusion over t 'excess' coordinates (also shifted by g); fibre integrals per group are
     Dirichlet; a single m_nu splits over the three groups multiplicatively part by part;
  3. the remaining 2-D integral over {a > B - s g, a + b <= R - (s+t) g} is a closed-form incomplete beta.
"""
from functools import lru_cache
from math import comb, factorial
from collections import defaultdict, Counter
from itertools import product as iproduct
from flint import fmpq
from mono_engine import canon, remove_one


# ---------- products of monomial symmetric functions ----------
@lru_cache(maxsize=None)
def _fill(alpha, beta):
    """ordered fillings of active coordinates: {nu (sorted multiset of x+y): count}."""
    if not alpha and not beta:
        return {(): 1}
    out = defaultdict(int)
    for x in [0] + sorted(set(alpha)):
        a2 = remove_one(alpha, x) if x else alpha
        for y in [0] + sorted(set(beta)):
            if x == 0 and y == 0:
                continue
            b2 = remove_one(beta, y) if y else beta
            for nu, c in _fill(a2, b2).items():
                out[canon(nu + (x + y,))] += c
    return dict(out)


@lru_cache(maxsize=None)
def mmul(alpha, beta):
    """m_alpha * m_beta = sum c_nu m_nu (independent of the number of variables, as long as it is >= len(nu))."""
    a, b = (alpha, beta) if alpha <= beta else (beta, alpha)
    out = {}
    for nu, f in _fill(a, b).items():
        c = f
        for m in Counter(nu).values():
            c *= factorial(m)
        q, r = divmod(c, factorial(len(nu)))
        assert r == 0
        out[nu] = q
    return out


def square(witness, kmax):
    """F^2 grouped: {nu: {A: coeff}} with F = sum c (R-|x|)^a m_alpha."""
    by = defaultdict(dict)
    for a, al, c in witness:
        by[al][a] = by[al].get(a, 0) + c
    als = sorted(by)
    sq = defaultdict(lambda: defaultdict(int))
    for i, a1 in enumerate(als):
        v1 = by[a1]
        for a2 in als[: i + 1]:
            v2 = by[a2]
            conv = defaultdict(int)
            for p, x in v1.items():
                for q, y in v2.items():
                    conv[p + q] += x * y
            sym = 1 if a1 == a2 else 2
            for nu, c in mmul(a1, a2).items():
                if len(nu) > kmax:
                    continue
                d = sq[nu]
                for A, v in conv.items():
                    d[A] += sym * c * v
    return {nu: {A: v for A, v in d.items() if v} for nu, d in sq.items()}


# ---------- fibre weights ----------
def nmono(mu, n):
    if len(mu) > n:
        return 0
    out = factorial(n) // factorial(n - len(mu))
    for m in Counter(mu).values():
        out //= factorial(m)
    return out


@lru_cache(maxsize=None)
def shift_part(e, g):
    """(g + v)^e fibre weight: list over j of C(e,j) g^(e-j) j!  (degree j in the group variable)."""
    return tuple(comb(e, j) * g ** (e - j) * factorial(j) for j in range(e + 1))


@lru_cache(maxsize=None)
def shift_profile(mu, g):
    prof = [fmpq(1)]
    for e in mu:
        sp = shift_part(e, g)
        new = [fmpq(0)] * (len(prof) + e)
        for i, x in enumerate(prof):
            if x != 0:
                for j, y in enumerate(sp):
                    new[i + j] += x * y
        prof = new
    return tuple(prof)


def splits3(nu, maxU, maxZ, maxW):
    """all (nuU, nuZ, nuW) sub-multiset splits with length limits."""
    items = sorted(Counter(nu).items(), reverse=True)
    out = []
    def rec(i, U, Z, W):
        if len(U) > maxU or len(Z) > maxZ or len(W) > maxW:
            return
        if i == len(items):
            out.append((tuple(U), tuple(Z), tuple(W)))
            return
        v, m = items[i]
        for u in range(m + 1):
            for z in range(m - u + 1):
                w = m - u - z
                rec(i + 1, U + [v] * u, Z + [v] * z, W + [v] * w)
    rec(0, [], [], [])
    return out


# ---------- the bad-set integral ----------
def inc_beta_lower(P, N, a0, R1):
    """int_{a0}^{R1} a^P (R1 - a)^N da, exact (0 <= a0 < R1)."""
    full = R1 ** (P + N + 1) * fmpq(factorial(P) * factorial(N), factorial(P + N + 1))
    if a0 == 0:
        return full
    # int_0^{a0} a^P (R1-a)^N da = sum_k C(N,k) R1^(N-k) (-1)^k a0^(P+k+1)/(P+k+1)
    low = sum((comb(N, k) * R1 ** (N - k) * (-1) ** k * a0 ** (P + k + 1) / (P + k + 1) for k in range(N + 1)), fmpq(0))
    return full - low


def bad_strata(k, R, g, B):
    Cf = (lambda m, B=B: B[min(m, len(B) - 1)]) if isinstance(B, tuple) else (lambda m, B=B: B)
    out = []
    for s in range(1, k + 1):
        if s * g >= R:
            break
        n = k - s
        for t in range(0, n + 1):
            R1 = R - (s + t) * g
            a0 = max(fmpq(0), Cf(s) - s * g)
            if R1 <= a0:
                break
            out.append((s, t))
    return out


def I_bad_stratum(sq, k, R, g, B, s, t):
    Cf = (lambda m, B=B: B[min(m, len(B) - 1)]) if isinstance(B, tuple) else (lambda m, B=B: B)
    n = k - s
    R1 = R - (s + t) * g
    a0 = max(fmpq(0), Cf(s) - s * g)
    mult = comb(k, s) * comb(n, t) * (-1) ** t
    geom = {}
    total = fmpq(0)
    for nu, Apoly in sq.items():
        W = defaultdict(lambda: fmpq(0))   # (P, Q) -> weight
        for nU, nZ, nW in splits3(nu, s, t, n - t):
            cnt = nmono(nU, s) * nmono(nZ, t) * nmono(nW, n - t)
            if not cnt:
                continue
            fw = 1
            for e in nW:
                fw *= factorial(e)
            pu = shift_profile(nU, g)
            pz = shift_profile(nZ, g) if t else (fmpq(1),)
            if not t and nZ:
                continue
            dW = sum(nW)
            for jU, xu in enumerate(pu):
                if xu == 0:
                    continue
                for jZ, xz in enumerate(pz):
                    if xz == 0:
                        continue
                    W[(s - 1 + jU, (n - 1 + jZ + dW) if n else -1)] += cnt * fw * xu * xz
        for (P, Q), w in W.items():
            fP = fmpq(1, factorial(P)); fQ = fmpq(1, factorial(Q)) if Q >= 0 else fmpq(1)
            for A, c in Apoly.items():
                key = (P, Q, A)
                if key not in geom:
                    if Q < 0:   # no smooth coordinates: only the rough group, int_{a0}^{R1} a^P (R1-a)^A da
                        geom[key] = inc_beta_lower(P, A, a0, R1)
                    else:
                        N = Q + A + 1
                        geom[key] = fmpq(factorial(Q) * factorial(A), factorial(N)) * inc_beta_lower(P, N, a0, R1)
                total += c * w * fP * fQ * geom[key]
    return mult * total


def _job(args):
    sq, k, R, g, B, s, t = args
    v = I_bad_stratum(sq, k, R, g, B, s, t)
    return s, t, (int(v.p), int(v.q))


def I_bad(sq, k, R, g, B, workers=14, verbose=False):
    from multiprocessing import get_context
    import time
    t0 = time.time()
    st = bad_strata(k, R, g, B)
    tot = fmpq(0)
    with get_context("fork").Pool(workers) as pool:
        for s, t, (p, q) in pool.imap_unordered(_job, [(sq, k, R, g, B, s, t) for s, t in st]):
            tot += fmpq(p, q)
            if verbose:
                print(f"   stratum s={s} t={t} done [{time.time()-t0:.0f}s]", flush=True)
    return tot

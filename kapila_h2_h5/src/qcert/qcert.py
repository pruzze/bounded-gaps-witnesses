"""Rigorous lower bound for the quotient Q_{k,eps}(F_0) of a Kapila trial (hexagon:2610.00419v2, eq. (4.4)),
in Arb ball arithmetic via python-flint. Own code; from the bundle only the trial JSON is read.

  Q = (k/m2) N / D,
  D = E[P(S_k, Z)^2]                       = sum_W sum_cells a_{W,cell} mu^{(k)}_W(cell),
  N = E[1{S_{k-1} <= tau} Pm(S_{k-1}, Z)^2] = sum_W int_0^tau K_W(s) mu^{(k-1)}_W(ds),
  mu^{(n)}_W(ds) = E[Z_W ; S_n in ds],  K_W = sum over (B, B') with B + B' = W of kappa_B kappa_B',
  kappa_B(s) = sum_e sum_m coef[B][e][m] Wc_e(e_m - s),  Wc_e(x) = int_0^{clamp(x, 0, T)} (t/T)^e g(t) dt.

Distribution functions and densities of mu_W come from the Fourier data Phi_W(xi_j), xi_j = 2 pi j / P0,
j <= J (transforms.py), with three separately bounded errors:
  * truncation j > J: |Phi_W(xi)| <= m_W |phi(xi)|^{n-4} (|W| <= 4), and sum_{j>J} |phi(xi_j)|^{n-4} is bounded by
    Lipschitz blocks (|phi'| <= E X) and, for large xi, by |phi(xi)| <= 2 f(0)/xi (f decreasing);
  * periodisation: 0 <= C_per(y) - C(y) <= E[Z_W 1{S_n >= P0}] <= e^{-lam P0} E[Z_W e^{lam S_n}];
  * the kink points e_m - T are balls; the segment ends use their midpoints (error term below).
The numerator integral is done segment by segment (between consecutive break points e_m and e_m - T) with
acb.integral (Petras); on a segment kappa_B is analytic: a polynomial in s plus logarithms log(c + b(e_m - s)),
where the logarithms of far edges are expanded in a power series with explicit remainder and the near ones are
evaluated directly with the analytic flag forwarded (non-finite ball on the branch cut).

Usage: python qcert.py <trial.json> [workers]
"""
import itertools
import math
import multiprocessing as mp
import sys
import time

from flint import acb, arb, ctx

import transforms as tf
from transforms import Trial, log, set_partitions, weighted

PREC = 192
ctx.prec = PREC
P0 = arb(int(sys.argv[4]) if len(sys.argv) > 4 else 8)
LAMS = [arb(x) for x in (2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 40)]
NSER = 80
NTAY = 28
NAN = acb(arb("nan"))
INSPECT = None
TEST = len(sys.argv) > 3 and sys.argv[3] == "test"

path = sys.argv[1]
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
tr = Trial(path)
k = tr.k
b, c, T, ag, m2 = tr.b, tr.c, tr.T, tr.ag, tr.m2
edges, psi, basis = tr.edges, tr.psi, tr.basis
ncell = len(edges) - 1
log(f"{path}: k = {k}, basis {basis}, {ncell} cells, tau = {tr.tau.str(15)}, prec {PREC}")
log(f"c = {c.str(15)}, T = {T.str(15)}, a_g = {ag.str(15)}, m2 = {m2.str(20)}")

key = lambda A, B: tuple(sorted(A + B))

# ---------------- marginal recipe ----------------
recipe = {}
for Ai, A in enumerate(basis):
    for keep in itertools.product((0, 1), repeat=len(A)):
        Bp = tuple(sorted(p for p, kp in zip(A, keep) if kp))
        e = sum(p for p, kp in zip(A, keep) if not kp)
        recipe.setdefault(Bp, []).append((Ai, e))
keys = sorted(recipe, key=lambda B: (len(B), B))
evals = sorted({e for lst in recipe.values() for _, e in lst})
# dpsi[A][m] = psi[A][m-1] - psi[A][m], m = 0..ncell (psi[-1] = psi[ncell] = 0)
dpsi = [[(row[m - 1] if m >= 1 else arb(0)) - (row[m] if m < ncell else arb(0)) for m in range(ncell + 1)]
        for row in psi]
coef = {B: {e: [arb(0)] * (ncell + 1) for e in evals} for B in keys}
for B in keys:
    for Ai, e in recipe[B]:
        coef[B][e] = [x + y for x, y in zip(coef[B][e], dpsi[Ai])]
pairs_num = []                                     # (B, B', multiplicity, W)
for i, B in enumerate(keys):
    for B2 in keys[i:]:
        pairs_num.append((B, B2, 1 if B == B2 else 2, key(B, B2)))
W_num = sorted({p[3] for p in pairs_num}, key=lambda W: (len(W), W))
W_den = sorted({key(A, B) for A in basis for B in basis}, key=lambda W: (len(W), W))
assert max(len(W) for W in W_num + W_den) <= 4
PMAX = max(sum(W) for W in W_num + W_den)
log(f"marginal keys {keys}, e-values {evals}; {len(W_num)} numerator and {len(W_den)} denominator weights")

# W_e(x) for 0 <= x <= T:  Poly_e(x) + lam_e log(1 + b x / c)
negcb = -c / b
lam_e = {e: negcb ** e / (b * T ** e) for e in evals}
polyc = {}                                          # polyc[e][i] = coefficient of x^i in Poly_e
for e in evals:
    pc = [arb(0)] * (e + 2)
    for i in range(1, e + 1):
        pc[i] = negcb ** (e - i) / (b * i) / T ** e
    pc[e + 1] += ag / (e + 1) / T ** e
    polyc[e] = pc


def W_full(e, x):
    v = sum((pc * x ** i for i, pc in enumerate(polyc[e]) if i), arb(0))
    return v + lam_e[e] * (1 + b * x / c).log()


WT = {e: W_full(e, T) for e in evals}
for e in evals:   # cross-check the closed form against direct integration at x = T/3
    x = T / 3
    direct = acb.integral(lambda t, _: (t / T) ** e * (1 / (c + b * t) + ag), 0, x, rel_tol=arb(2) ** -150).real
    assert direct.overlaps(W_full(e, x)), (e, direct, W_full(e, x))
supkappa = {B: sum((abs(coef[B][e][m]) * WT[e] for e in evals for m in range(ncell + 1)), arb(0)) for B in keys}

# ---------------- Fourier data ----------------
EX = (T * tr.moment(1, arb(0))).real
f0 = (1 / c + ag) ** 2 / m2
log(f"E X = {EX.str(12)}, f(0) = {f0.str(12)}")


def xi_of(j):
    return 2 * arb.pi() * j / P0


# choose J: |phi(xi_J)|^(k-5) < 1e-45
J = 200
while True:
    a = tr.phi_closed(xi_of(J)).abs_upper()
    if bool((k - 5) * a.log() < -104):
        break
    J = int(J * 1.1)
log(f"J = {J} frequencies (|phi(xi_J)| = {tr.phi_closed(xi_of(J)).abs_upper().str(15)})")

phis, mus = [], []
for j in range(1, J + 1):
    xi = xi_of(j)
    ph = tr.phi_closed(xi)
    phis.append(ph)
    mus.append({p: tr.moment(p, xi) for p in range(2, PMAX + 1)})
    if j in (1, 10, 100):
        chk = tr.moment(0, xi)
        assert chk.overlaps(ph), (j, chk, ph)
    if j % 500 == 0:
        log(f"  transforms {j}/{J}")
mu0 = {p: tr.moment(p, arb(0)) for p in range(2, PMAX + 1)}
mom = {}                                            # m_W^{(n)} = E Z_W
for n in (k - 1, k):
    for W in set(W_num) | set(W_den):
        mom[n, W] = weighted(n, W, acb(1), mu0).real
Phi = {}
for n, Ws in ((k - 1, W_num), (k, W_den)):
    for W in Ws:
        Phi[n, W] = [weighted(n, W, phis[j], mus[j]) for j in range(J)]
log("weighted transforms done")

# truncation: U_n >= sum_{j > J} |phi(xi_j)|^(n-4)
def tail_sum(n):
    U = arb(0)
    ja = J + 1
    xi_big = 4 * f0                                  # beyond: |phi| <= 2 f0 / xi <= 1/2
    nblk = 0
    rho = 0.01
    while bool(xi_of(ja) < xi_big):
        w = max(1, int(ja * rho))
        while True:
            q = tr.phi_closed(xi_of(ja)).abs_upper() + xi_of(w) * EX
            if bool(q < 1) and bool((n - 4) * q.log() < -230) or w == 1:
                break
            w = max(1, w // 2)
        if not bool(q < 1):
            raise SystemExit(f"tail block at j = {ja} not certified")
        U += w * (q.log() * (n - 4)).exp()
        ja += w
        nblk += 1
        rho = min(0.05, rho * 1.05)
    q0 = 2 * f0 / xi_of(ja)
    U += (q0.log() * (n - 4)).exp() * (1 + arb(ja) / (n - 5))
    return U, nblk


tails = {}
for n in (k - 1, k):
    tails[n], nb = tail_sum(n)
    log(f"truncation: sum_(j>J) |phi|^(n-4) <= {tails[n].upper().str(5)}  (n = {n}, {nb} blocks)")

# periodisation: eps_per[n, W] >= E[Z_W 1{S_n >= P0}]
expmom = []
for LAM in LAMS:
    expmom.append((LAM, tr.moment(0, acb(0, -LAM)), {p: tr.moment(p, acb(0, -LAM)) for p in range(2, PMAX + 1)}))
eps_per = {}
for n, Ws in ((k - 1, W_num), (k, W_den)):
    for W in Ws:
        eps_per[n, W] = min((((-LAM * P0).exp() * weighted(n, W, mphiL, mmuL).real).upper()
                             for LAM, mphiL, mmuL in expmom), key=lambda x: float(x))
log(f"periodisation (best lambda per weight): largest bound "
    f"{max((eps_per[x] / mom[x]).upper() for x in eps_per).str(5)} (relative to m_W)")

# ---------------- denominator ----------------
TWO_PI_I = 2 * arb.pi() * acb(0, 1)


def unit_powers(y):
    """[e^{-i xi_j y} for j = 1..J]. Repeated multiplication by a rotation inflates rectangular complex balls by up to
    sqrt2 per step (wrapping effect), so the power is recomputed directly by exp every 32 steps."""
    w1 = acb(0, -2 * arb.pi() * y / P0).exp()
    out = []
    wp = acb(1)
    for j in range(1, J + 1):
        if j % 32 == 0:
            wp = acb(0, -2 * arb.pi() * j * y / P0).exp()
        else:
            wp = wp * w1
        out.append(wp)
    return out


def horner(co, x):
    """sum co[i] x^i by Horner in scalar ball arithmetic (no polynomial routines of FLINT are used)."""
    v = acb(0)
    for cc in reversed(co):
        v = v * x + cc
    return v


def cdf_data(n, W):
    beta_c = [Phi[n, W][j] / (TWO_PI_I * (j + 1)) for j in range(J)]
    S = sum(beta_c, acb(0))
    poly = beta_c
    err_tr = 2 / arb.pi() * mom[n, W] * tails[n] / (J + 1)
    return S, poly, err_tr


def cdf_raw(n, W, data, y):
    """C_per,J(y) without error terms."""
    S, poly, err_tr = data
    pw = unit_powers(y) if not isinstance(y, list) else y
    return mom[n, W] * y / P0 + 2 * (S - sum((bc * p for bc, p in zip(poly, pw)), acb(0))).real


def cdf(n, W, data, y):
    """Ball for the true C_W(y) = E[Z_W 1{S_n <= y}], 0 <= y < P0."""
    S, poly, err_tr = data
    pw = UPW[y] if isinstance(y, int) else unit_powers(y)
    yy = edges[y] if isinstance(y, int) else y
    val = mom[n, W] * yy / P0 + 2 * (S - sum((bc * p for bc, p in zip(poly, pw)), acb(0))).real
    lo = val - err_tr - eps_per[n, W]
    hi = val + err_tr
    return arb.union(lo, hi)


UPW = {m: unit_powers(edges[m]) for m in range(1, ncell + 1)}
D = arb(0)
for W in W_den:
    data = cdf_data(k, W)
    Cs = [arb(0)] + [cdf(k, W, data, m) for m in range(1, ncell + 1)]
    for m in range(ncell):
        a = sum((psi[Ai][m] * psi[Bi][m] for Ai, A in enumerate(basis) for Bi, B in enumerate(basis)
                 if key(A, B) == W), arb(0))
        if not a.is_zero():
            D += a * (Cs[m + 1] - Cs[m])
log(f"D = E[P^2] = {D.str(20, radius=True)}   (bundle den = {tr.spec.get('den')})")

# ---------------- numerator ----------------
# density of mu^{(k-1)}_W (periodised, truncated at J): rho(s) = (1/P0)[m_W + 2 Re sum_j Phi_j e^{-i xi_j s}];
# on a segment it is replaced by its Taylor polynomial at the centre s_c, with real coefficients
#   a_n = (2/P0) (2 pi/P0)^n / n! Re[(-i)^n sum_j Phi_j e^{-i xi_j s_c} j^n]   (+ m_W/P0 for n = 0)
# and remainder <= (2/P0) A_W t^(N+1)/(N+1)! e^t,  t = xi_J |z - s_c|,  A_W = sum_j |Phi_j|.
cdf_num = {W: cdf_data(k - 1, W) for W in W_num}
AW = {W: sum((x.abs_upper() for x in Phi[k - 1, W]), arb(0)) for W in W_num}
XIJ = xi_of(J)
FACT = [arb(1)]
for i in range(1, NTAY + 2):
    FACT.append(FACT[-1] * i)
NEG_I_POW = [acb(1), acb(0, -1), acb(-1), acb(0, 1)]


def taylor_density(sc):
    pw = unit_powers(sc)
    out = {}
    for W in W_num:
        v = [Phi[k - 1, W][j] * pw[j] for j in range(J)]
        co = []
        jp = [arb(j + 1) for j in range(J)]
        for n in range(NTAY + 1):
            if n == 0:
                ssum = sum(v, acb(0))
            else:
                ssum = acb(0)
                for j in range(J):
                    v[j] = v[j] * jp[j]
                    ssum += v[j]
            a = 2 / P0 * (2 * arb.pi() / P0) ** n / FACT[n] * (NEG_I_POW[n % 4] * ssum).real
            if n == 0:
                a += mom[k - 1, W] / P0
            co.append(acb(a))
        out[W] = co
    return out

# break points
tau = tr.tau
bps = [(edges[m].mid(), arb(0)) for m in range(1, ncell + 1) if bool(edges[m] < tau)]
kink_rad = arb(0)
for m in range(1, ncell + 1):
    x = edges[m] - T
    if bool(x > 0) and bool(x < tau):
        bps.append((x.mid(), x.rad()))
        kink_rad = arb.max(kink_rad, x.rad())
    elif x.overlaps(arb(0)) or x.overlaps(tau):
        raise SystemExit("kink point not separated from 0 or tau")
pts = sorted({p for p, _ in bps}, key=lambda a: float(a))
pts = [arb(0)] + [arb(p) for p in pts] + [tau]
for a1, a2 in zip(pts[:-1], pts[1:]):
    assert bool(a1 < a2)
nseg = len(pts) - 1
log(f"{nseg} segments on [0, tau]; kink radius {kink_rad.str(3)}")
edge_mid = [arb(e.mid()) for e in edges]
T_mid = arb(T.mid())


def segment_integral(si):
    lo, hi = pts[si], pts[si + 1]
    h = hi - lo
    s0 = hi
    smid = (lo + hi) / 2
    act, sat = [], []
    for m in range(1, ncell + 1):
        em = edges[m]
        if bool(em <= smid) or em.overlaps(smid):
            if not bool(em <= lo) and not em.overlaps(lo):
                raise RuntimeError("edge inside segment")
            continue
        km = em - T
        if bool(km >= smid):
            sat.append(m)
        else:
            act.append(m)
    # polynomial part and log constants
    cpoly = {}
    near, far = [], []
    for m in act:
        x = c + b * (edges[m] - s0)
        (far if bool(x >= 8 * b * h) else near).append(m)
    rmax = arb(0)
    for m in far:
        rmax = arb.max(rmax, b * h / (c + b * (edges[m] - s0)))
    logc = c.log()
    for B in keys:
        const = sum((coef[B][e][m] * WT[e] for e in evals for m in sat), arb(0))
        # moments of active edges
        pc = [arb(0)] * (max(evals) + 2)
        for e in evals:
            mo = [arb(0)] * (e + 2)
            for m in act:
                d = edges[m] - s0
                cm = coef[B][e][m]
                if cm.is_zero():
                    continue
                dp = arb(1)
                for l in range(e + 2):
                    mo[l] += cm * dp
                    dp *= d
            for i in range(1, e + 2):
                ci = polyc[e][i]
                if ci.is_zero():
                    continue
                for l in range(i + 1):
                    pc[l] += ci * math.comb(i, l) * mo[i - l] * h ** l
        L = {m: sum((coef[B][e][m] * lam_e[e] for e in evals), arb(0)) for m in act}
        for m in far:
            const += L[m] * ((c + b * (edges[m] - s0)).log() - logc)
        for m in near:
            const -= L[m] * logc
        pc[0] += const
        # far series coefficients: sum_far L r_m^n (-1)^(n+1)/n Y^n, r_m = b h / x_m
        ser = [arb(0)] * (NSER + 1)
        Lfar_abs = arb(0)
        for m in far:
            r = b * h / (c + b * (edges[m] - s0))
            rn = arb(1)
            Lm = L[m]
            Lfar_abs += abs(Lm)
            for n in range(1, NSER + 1):
                rn *= r
                ser[n] += Lm * rn
        for n in range(1, NSER + 1):
            ser[n] = ser[n] * (1 if n % 2 else -1) / n
        for i in range(len(pc)):
            ser[i] += pc[i]
        cpoly[B] = ([acb(x) for x in ser], [(L[m], c + b * (edges[m] - s0)) for m in near], Lfar_abs)
    if not act:
        # kappa constant on the segment: integral = sum_W K_W (C_W(hi) - C_W(lo)), periodised and truncated
        kapc = {B: cpoly[B][0][0] for B in keys}
        tot = acb(0)
        for W in W_num:
            KW = sum((mult * kapc[B] * kapc[B2] for B, B2, mult, WW in pairs_num if WW == W), acb(0))
            tot += KW * (cdf_raw(k - 1, W, cdf_num[W], hi) - cdf_raw(k - 1, W, cdf_num[W], lo))
        return si, tot, 0, 0
    # the Taylor density is used within |z - s_c| <= 1/xi_J, so long segments are cut into pieces
    npiece = max(1, int(float((h * XIJ).upper())) + 1)
    if npiece > 1:
        tot = acb(0)
        for ip in range(npiece):
            a1 = lo + h * ip / npiece if ip else lo
            a2 = lo + h * (ip + 1) / npiece if ip + 1 < npiece else hi
            a1, a2 = arb(a1.mid()) if ip else a1, arb(a2.mid()) if ip + 1 < npiece else a2
            tot += piece_integral(a1, a2, (a1 + a2) / 2, cpoly, far, rmax, s0, h)
        return si, tot, len(act), len(near)
    return si, piece_integral(lo, hi, smid, cpoly, far, rmax, s0, h), len(act), len(near)


def piece_integral(lo, hi, sc, cpoly, far, rmax, s0, h):
    tay = taylor_density(sc)

    def parts(z, analytic):
        """(kappa_B for all B, rho_W for all W) at z, or None if not certified analytic there."""
        Y = (s0 - z) / h
        q = Y.abs_upper() * rmax
        if far and not bool(q < 0.5):
            return None
        rem_q = (q ** (NSER + 1) / ((NSER + 1) * (1 - q))) if far else arb(0)
        kap = {}
        for B in keys:
            poly, nr, Lfar_abs = cpoly[B]
            v = horner(poly, Y)
            for Lm, x in nr:
                arg = x + b * h * Y
                lg = arg.log(analytic=analytic)
                if not lg.is_finite():
                    return None
                v += Lm * lg
            if far:
                v += acb(arb(0, (Lfar_abs * rem_q).upper()), arb(0, (Lfar_abs * rem_q).upper()))
            kap[B] = v
        x = z - sc
        t = XIJ * x.abs_upper()
        if not bool(t < 3):
            return None
        tn = t ** (NTAY + 1) / FACT[NTAY + 1] * t.exp() * 2 / P0
        rho = {}
        for W in W_num:
            rr = (tn * AW[W]).upper()
            rho[W] = horner(tay[W], x) + acb(arb(0, rr), arb(0, rr))
        return kap, rho

    def integrand(z, analytic):
        pr = parts(z, analytic)
        if pr is None:
            return NAN
        kap, rho = pr
        tot = acb(0)
        for W in W_num:
            rhoW = rho[W]
            KW = sum((mult * kap[B] * kap[B2] for B, B2, mult, WW in pairs_num if WW == W), acb(0))
            tot += KW * rhoW
        return tot

    if INSPECT is not None:
        return INSPECT(lo, hi, parts)
    return acb.integral(integrand, lo, hi, abs_tol=arb(2) ** -125, rel_tol=arb(2) ** -70,
                        eval_limit=10 ** 7, depth_limit=10 ** 5)


def pack(x):
    """Exact, picklable form of an arb ball (mid and radius as mantissa-exponent integers)."""
    m, e = x.mid().man_exp()
    r, f = x.rad().man_exp()
    return int(m), int(e), int(r), int(f)


def unpack(t):
    m, e, r, f = t
    return arb(m) * arb(2) ** e + arb(0, arb(r) * arb(2) ** f)


def seg_packed(si):
    si, val, na, nn = segment_integral(si)
    val = acb(val)
    return si, pack(val.real), pack(val.imag), na, nn


def run_segments():
    ctxm = mp.get_context("fork")
    out = [None] * nseg
    t0 = time.time()
    done = 0
    with ctxm.Pool(WORKERS) as pool:
        for si, vr, vi, na, nn in pool.imap_unordered(seg_packed, range(nseg), chunksize=4):
            out[si] = (acb(unpack(vr), unpack(vi)), na, nn)
            done += 1
            if done % 100 == 0 or done == nseg:
                log(f"  segments {done}/{nseg}")
    return out


if __name__ == "__main__":
    if TEST:
        import random
        random.seed(1)
        for si in [0, nseg - 1] + random.sample(range(nseg), 4):
            t0 = time.time()
            _, val, na, nn = segment_integral(si)
            log(f"segment {si} [{pts[si].str(8)}, {pts[si+1].str(8)}]: {na} active, {nn} near, "
                f"{val.real.str(15, radius=True)}  ({time.time()-t0:.2f}s)")
        sys.exit(0)
    res = run_segments()
    Nsum = acb(0)
    for val, na, nn in res:
        Nsum += val
    assert Nsum.imag.contains(0) or bool(abs(Nsum.imag).upper() < arb(2) ** -100)
    Nraw = Nsum.real
    # error terms
    supK = {W: sum((mult * supkappa[B] * supkappa[B2] for B, B2, mult, WW in pairs_num if WW == W), arb(0))
            for W in W_num}
    err_trunc = sum((supK[W] * tau * 2 / P0 * mom[k - 1, W] * tails[k - 1] for W in W_num), arb(0))
    err_alias = sum((supK[W] * eps_per[k - 1, W] for W in W_num), arb(0))
    suprho = {W: (mom[k - 1, W] + 2 * sum((x.abs_upper() for x in Phi[k - 1, W]), arb(0))
                  + 2 * mom[k - 1, W] * tails[k - 1]) / P0 for W in W_num}
    nkink = sum(1 for _, r in bps if not r.is_zero())
    err_kink = 4 * kink_rad * nkink * sum((supK[W] * suprho[W] for W in W_num), arb(0))
    err = err_trunc + err_alias + err_kink
    N = Nraw + arb(0, err.upper())
    Q = arb(k) / m2 * N / D
    print("RESULT")
    print(f"  trial {path}")
    print(f"  J = {J}, P0 = {P0}, prec = {PREC}, segments = {nseg}")
    print(f"  D = {D.str(20, radius=True)}")
    print(f"  N (integrals) = {Nraw.str(20, radius=True)}")
    print(f"  error terms: truncation {err_trunc.upper().str(3)}, periodisation {err_alias.upper().str(3)}, "
          f"kinks {err_kink.upper().str(3)}")
    print(f"  N = {N.str(20, radius=True)}")
    print(f"  Q = {Q.str(20, radius=True)}")
    print(f"  Q lower endpoint = {Q.lower().str(20)}")
    log("done")

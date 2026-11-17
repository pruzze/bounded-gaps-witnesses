"""Independent floating-point evaluation of the quotient Q_{k,eps}(F_0) for the trials of Kapila's paper
(hexagon:2610.00419v2, eqs. (4.1)-(4.5)), by Fourier inversion. No lattice in the single coordinate, so it works
for every k of the paper (2305, 110070, 5090000, 232500000).

Written from the formulas of the paper. From the bundle only the trial file is read (k, theta, beta,
relative_shift, eps_cells, delta, cell edges E, basis, profiles). None of the bundle's code is imported or run.
Not rigorous: double precision, Gauss-Legendre quadrature, truncated Fourier series.

  F_0(t) = P(S, Z) prod g(t_i),  P = sum_A psi_A(S) Z_A,  Z_j = sum (t_i/T)^j,
  g(t) = 1/(c + (k-1) t) + a_g on [0, T],  c = theta/log k,  T = beta/log k,  a_g = r_g/(c + (k-1) T),
  Q = (k/m2) E[ 1{S_{k-1} <= tau} Pm(S_{k-1}, Z)^2 ] / E[ P(S_k, Z)^2 ],   X_i i.i.d. with density g^2/m2,
  Pm(s, Z) = int_0^T P(s + t, Z + (t/T)^.) g(t) dt.

Method: phi(xi) = E e^{i xi X} and mu_j(xi) = E[(X/T)^j e^{i xi X}] by quadrature at the frequencies 2 pi m / PER;
the measures E[Z_A ; S_n in ds] have transform sum over set partitions (eq. (4.5)); their distribution functions
at the cell edges give E[P^2]; their densities on a graded quadrature grid give the numerator.

Usage: python q_fourier.py <trial.json> [frequencies M = 3000] [period = 8]"""
import itertools
import json
import sys
import time

import numpy as np

spec = json.load(open(sys.argv[1]))
M = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
PER = float(sys.argv[3]) if len(sys.argv) > 3 else 8.0
t00 = time.time()
log = lambda s: print(f"[{time.time()-t00:5.0f}s] {s}", flush=True)
k = spec["k"]; theta = spec["theta"]; beta = spec["beta"]; rg = spec["relative_shift"]
basis = [tuple(a) for a in spec["basis"]]
edges = np.array(spec["E"], dtype=float) * spec["delta"]
psi = np.array(spec["profiles"], dtype=float)                    # psi[A][cell], cell m = [edges[m], edges[m+1])
assert psi.shape == (len(basis), len(edges) - 1)
eps = spec["eps_cells"] * spec["delta"]; sigma, tau = 1 + eps, 1 - eps
assert abs(edges[-1] - sigma) < 1e-15, (edges[-1], sigma)
b = k - 1; c = theta / np.log(k); T = beta / np.log(k); ag = rg / (c + b * T)
xg, wg = np.polynomial.legendre.leggauss(16)


def panels(lo, hi, first, ratio, maxlen):
    """Break [lo, hi] into panels: the first of length `first`, then growing by `ratio`, capped at `maxlen`."""
    out = [lo]; length = first
    while out[-1] < hi:
        out.append(min(hi, out[-1] + length)); length = min(maxlen, length * ratio)
    return np.array(out)


def gl(breaks, x=xg, w=wg):
    a, bb = breaks[:-1], breaks[1:]
    nodes = (a[:, None] + bb[:, None]) / 2 + (bb - a)[:, None] / 2 * x[None, :]
    weights = (bb - a)[:, None] / 2 * w[None, :]
    return nodes.ravel(), weights.ravel()


# one-coordinate transforms at xi_m = 2 pi m / PER, m = 0..M
xi = 2 * np.pi * np.arange(M + 1) / PER
tq, wq = gl(panels(0.0, T, c / b / 8, 1.5, PER / (2 * M)))
g2 = (1 / (c + b * tq) + ag) ** 2
m2 = float(np.sum(wq * g2))
m2_exact = (1 / c - 1 / (c + b * T)) / b + (2 * ag / b) * np.log1p(b * T / c) + ag * ag * T
dens = wq * g2 / m2
ang = np.outer(xi, tq)
phim1 = (-2 * np.sin(ang / 2) ** 2 + 1j * np.sin(ang)) @ dens       # phi - 1, accurately
logphi = np.log1p(phim1)
ex = np.exp(1j * ang)
JMAX = 4 * max([max(a) for a in basis if a] + [0])
mu = {j: ex @ (dens * (tq / T) ** j) for j in range(1, JMAX + 1)}
del ang, ex
log(f"k = {k}, c = {c:.6g}, T = {T:.6g}, a_g = {ag:.6g}; m2 = {m2:.14e} (closed form {m2_exact:.14e}); "
    f"{len(tq)} quadrature nodes, {M} frequencies, period {PER}")


def set_partitions(items):
    if not items:
        yield []; return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        yield [[first]] + part
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]


def transform(n, wgt):
    """Fourier transform at xi_m of the measure E[ Z_A ; S_n in ds ] for the labelled list wgt (eq. (4.5))."""
    tr = np.zeros(M + 1, complex)
    for part in set_partitions(list(range(len(wgt)))):
        r = len(part); fall = 1.0
        for i in range(r): fall *= n - i
        term = fall * np.exp((n - r) * logphi)
        for blk in part: term = term * mu[sum(wgt[i] for i in blk)]
        tr += term
    return tr


def cdf_at(tr, e):
    """Distribution function at the points e of the measure with transform tr (periodised, period PER)."""
    e = np.asarray(e, float)
    z = (1 - np.exp(-1j * np.outer(e, xi[1:]))) / (1j * xi[1:])
    return (tr[0].real * e + 2 * (z @ tr[1:]).real) / PER


def density_at(trs, s):
    """Densities at the points s of the measures with transforms trs (array: weights x frequencies)."""
    out = np.empty((len(s), trs.shape[0]))
    for a in range(0, len(s), 8000):
        z = np.exp(-1j * np.outer(s[a:a + 8000], xi[1:]))
        out[a:a + 8000] = (trs[:, 0].real[None, :] + 2 * (z @ trs[:, 1:].T).real) / PER
    return out


key = lambda A, B: tuple(sorted(A + B))

# denominator E[P(S_k, Z)^2]
den = 0.0; cache = {}
for A, B in itertools.product(range(len(basis)), repeat=2):
    kk = key(basis[A], basis[B])
    if kk not in cache:
        F = cdf_at(transform(k, kk), edges); cache[kk] = np.diff(F)
    den += float(np.sum(psi[A] * psi[B] * cache[kk]))
tail = max(abs(transform(k, ())[-1]), abs(transform(k, key(basis[-1], basis[-1]))[-1]))
log(f"E[P^2] = {den:.16e}   (bundle: den = {spec.get('den')});  |transform| at the last frequency: {tail:.1e}")

# marginal profiles: Pm(s, Z) = sum_{A'} Psi_{A'}(s) Z_{A'},
# Psi_{A'} = sum over basis A and sub-lists kept = A' of  a_{A, e}(s),  e = sum of the removed exponents,
# a_{A, e}(s) = int_0^T psi_A(s + t) (t/T)^e g(t) dt = sum_m (psi_{A,m-1} - psi_{A,m}) W_e(edges[m] - s)
recipe = {}                                                        # A' -> list of (A index, e)
for Ai, A in enumerate(basis):
    for keep in itertools.product((0, 1), repeat=len(A)):
        Ap = tuple(sorted(p for p, kp in zip(A, keep) if kp)); e = sum(p for p, kp in zip(A, keep) if not kp)
        recipe.setdefault(Ap, []).append((Ai, e))
zero = np.zeros((len(basis), 1))
dpsi = np.concatenate([zero, psi], axis=1) - np.concatenate([psi, zero], axis=1)      # index m = 0..len(edges)-1


def H(j, z):
    """int_0^z v^j / (1 + v) dv."""
    z = np.asarray(z, float); out = np.empty_like(z); small = z < 0.5
    zs = z[small]; acc = np.zeros_like(zs); term = zs ** (j + 1)
    for n in range(0, 70):
        acc += term / (j + n + 1); term = -term * zs
    out[small] = acc
    zl = z[~small]; h = np.log1p(zl)
    for i in range(1, j + 1): h = zl ** i / i - h
    out[~small] = h
    return out


def W(e, x):
    """int_0^x (t/T)^e g(t) dt for 0 <= x, constant for x >= T."""
    x = np.clip(x, 0.0, T)
    return (c ** e / b ** (e + 1) * H(e, b * x / c) + ag * x ** (e + 1) / (e + 1)) / T ** e


# quadrature grid in s on [0, tau]: on every segment between consecutive break points the variable is
# x = (next edge above) - s on a graded grid, because Psi varies on the scale c/b just below each edge
kinks = edges - T                                                   # W is constant beyond x = T: a kink of Psi at edge - T
brk = np.unique(np.concatenate([edges[edges < tau], [tau], kinks[(kinks > 0) & (kinks < tau)]]))
x8, w8 = np.polynomial.legendre.leggauss(8)
S_nodes, S_w = [], []
xmin = c / b * 1e-5
for lo_, hi_ in zip(brk[:-1], brk[1:]):
    ref = edges[np.searchsorted(edges, hi_ - 1e-15)]                # smallest edge >= hi_
    x0, x1 = ref - hi_, ref - lo_
    xb = [x0] if x0 >= xmin else [x0, xmin]                        # the short first panel [0, xmin] is kept
    while xb[-1] < x1:
        step = min(0.002, max(xb[-1] * 0.6, xmin))
        xb.append(min(x1, xb[-1] + step))
    xn, xw = gl(np.array(xb), x8, w8)
    S_nodes.append(ref - xn); S_w.append(xw)
s_nodes = np.concatenate(S_nodes); s_w = np.concatenate(S_w)
log(f"numerator grid: {len(s_nodes)} nodes on [0, tau], {len(brk) - 1} segments")

need_e = sorted({e for lst in recipe.values() for _, e in lst})
Psi = {Ap: np.zeros(len(s_nodes)) for Ap in recipe}
for a in range(0, len(s_nodes), 20000):
    X = edges[None, :] - s_nodes[a:a + 20000, None]
    for e in need_e:
        We = W(e, X)
        for Ap, lst in recipe.items():
            for Ai, ee in lst:
                if ee == e: Psi[Ap][a:a + 20000] += We @ dpsi[Ai]
log("marginal profiles done")

keys = sorted(recipe)
wlist = sorted({key(A, B) for A in keys for B in keys})
trs = np.array([transform(k - 1, wk) for wk in wlist])
rho = density_at(trs, s_nodes)
col = {wk: i for i, wk in enumerate(wlist)}
# self-check of the s-grid and the densities: the grid integral of each density against the distribution function
worst = 0.0
for wk in wlist:
    direct = float(cdf_at(trs[col[wk]], [tau])[0]); grid = float(np.sum(s_w * rho[:, col[wk]]))
    worst = max(worst, abs(grid / direct - 1))
log(f"self-check: grid integral of the densities over [0, tau] against their distribution functions: largest relative difference {worst:.1e}")
integrand = np.zeros(len(s_nodes))
for A in keys:
    for B in keys:
        integrand += Psi[A] * Psi[B] * rho[:, col[key(A, B)]]
num = float(np.sum(s_w * integrand))
Q = k / m2 * num / den
log(f"numerator = {num:.16e}")
log(f"Q = {Q:.12f}")
print(f"bundle: Q_exploratory = {spec.get('Q_exploratory')}  (m = {spec.get('m')}, k = {k})")

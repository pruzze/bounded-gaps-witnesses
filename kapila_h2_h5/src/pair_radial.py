"""Independent floating-point check of the radial pair cost used for H_4 and H_5 in Kapila's paper
(hexagon:2610.00419v2, eqs. (8.7), (8.10), (10.3), (10.4)):  C = 2 zeta^2 J_K(K eta) E_nu[ R(S/sigma)^2 ],
R(x) = (1 + (eta - 1) x)/(1 - zeta x),  zeta = 2 rho sigma (limit),  nu = F_0^2 / I(F_0) for the k-dimensional source.
The expectation is a one-dimensional integral against the measures E[Z_A ; S_k in ds] (Fourier inversion as in
q_fourier.py). J_K by quadrature of the one-dimensional form of (8.7). Only the trial file is read.
Usage: python pair_radial.py <trial.json> <K> <eta numerator> <eta denominator> <omega>"""
import itertools, json, sys
import numpy as np
from scipy import integrate

spec = json.load(open(sys.argv[1])); K = int(sys.argv[2]); eta = int(sys.argv[3]) / int(sys.argv[4]); omega = float(sys.argv[5])
M, PER = 3000, 8.0
k = spec["k"]; theta = spec["theta"]; beta = spec["beta"]; rg = spec["relative_shift"]
basis = [tuple(a) for a in spec["basis"]]; edges = np.array(spec["E"], dtype=float) * spec["delta"]
psi = np.array(spec["profiles"], dtype=float); eps = spec["eps_cells"] * spec["delta"]; sigma = 1 + eps
b = k - 1; c = theta / np.log(k); T = beta / np.log(k); ag = rg / (c + b * T)
rho = 0.25 + omega; zeta = 2 * rho * sigma
xg, wg = np.polynomial.legendre.leggauss(16)
def gl(breaks, x=xg, w=wg):
    a, bb = breaks[:-1], breaks[1:]
    return ((a[:, None] + bb[:, None]) / 2 + (bb - a)[:, None] / 2 * x[None, :]).ravel(), ((bb - a)[:, None] / 2 * w[None, :]).ravel()
br = [0.0]; length = c / b / 8
while br[-1] < T:
    br.append(min(T, br[-1] + length)); length = min(PER / (2 * M), length * 1.5)
tq, wq = gl(np.array(br)); g2 = (1 / (c + b * tq) + ag) ** 2; m2 = float(np.sum(wq * g2)); dens = wq * g2 / m2
xi = 2 * np.pi * np.arange(M + 1) / PER; ang = np.outer(xi, tq)
logphi = np.log1p((-2 * np.sin(ang / 2) ** 2 + 1j * np.sin(ang)) @ dens); ex = np.exp(1j * ang)
mu = {j: ex @ (dens * (tq / T) ** j) for j in range(1, 4 * max(max(a) for a in basis if a) + 1)}
def set_partitions(items):
    if not items:
        yield []; return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        yield [[first]] + part
        for i in range(len(part)): yield part[:i] + [[first] + part[i]] + part[i + 1:]
def transform(n, wgt):
    tr = np.zeros(M + 1, complex)
    for part in set_partitions(list(range(len(wgt)))):
        r = len(part); fall = 1.0
        for i in range(r): fall *= n - i
        term = fall * np.exp((n - r) * logphi)
        for blk in part: term = term * mu[sum(wgt[i] for i in blk)]
        tr += term
    return tr
x4, w4 = np.polynomial.legendre.leggauss(6)
s, sw = gl(edges, x4, w4); cell = np.repeat(np.arange(len(edges) - 1), 6)
R2 = ((1 + (eta - 1) * s / sigma) / (1 - zeta * s / sigma)) ** 2
num = den = numS = 0.0; cache = {}
for A, B in itertools.product(range(len(basis)), repeat=2):
    kk = tuple(sorted(basis[A] + basis[B]))
    if kk not in cache:
        tr = transform(k, kk); d = np.empty(len(s))
        for a in range(0, len(s), 4000):
            z = np.exp(-1j * np.outer(s[a:a + 4000], xi[1:])); d[a:a + 4000] = (tr[0].real + 2 * (z @ tr[1:]).real) / PER
        cache[kk] = d
    w = psi[A][cell] * psi[B][cell] * cache[kk] * sw
    den += float(w.sum()); num += float((w * R2).sum()); numS += float((w * s / sigma).sum())
Mrad = num / den
a = K * eta; gam = np.sqrt(1 - 2 / K)
f1 = lambda t: 4 * np.arctanh(a * gam * t / (2 + (a - 2) * t)) / (a * gam * (2 + (a - 2) * t))
J = sum(integrate.quad(f1, lo_, hi_, limit=400, epsabs=0, epsrel=1e-12)[0] for lo_, hi_ in zip([0, 1e-6, 1e-4, 1e-2, 0.5, 0.99, 0.9999], [1e-6, 1e-4, 1e-2, 0.5, 0.99, 0.9999, 1]))
print(f"m = {spec.get('m')}, k = {k}, K = {K}, eta = {eta}, zeta = {zeta:.9f}")
print(f"E[P^2] = {den:.13e};  E_nu[R(S/sigma)^2] = {Mrad:.12f};  E_nu[S/sigma] = {numS / den:.12f};   R(1)^2 = {((eta)/(1-zeta))**2:.6f}")
print(f"J_K(K eta) = {J:.13e}")
print(f"C = 2 zeta^2 J E_nu[R^2] = {2 * zeta**2 * J * Mrad:.12e}")

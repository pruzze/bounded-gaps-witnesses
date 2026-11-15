"""Brute-force check of the table engine against the definition of the support and of I, J as stated in
Gap212Challenge/Basic.lean (AxiomMath/BGP212), one band, in dilated coordinates, for a small number of variables.

  T = {x in [0, inf)^K : sum x < L, sum_{x_i > g} x_i <= C(#{i : x_i > g})}
  I(F) = int_T F^2
  J(F) = int over (q, t, t') with sum q <= tau of F(q, t) F(q, t'),   F = 0 outside T

No cell decomposition is used here: the indicator of T is evaluated pointwise at uniform random points
(Monte Carlo, N million samples), and every Gram entry is compared with its standard error.
Usage: python bruteforce_def.py [K=4] [d=4] [N=200 (million samples)]
"""
import itertools
import sys
import tempfile

import numpy as np
from flint import ctx, fmpq

import gramtab as gt
from gram_report import rayleigh

K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
d = int(sys.argv[2]) if len(sys.argv) > 2 else 4
N = int(sys.argv[3]) if len(sys.argv) > 3 else 200
# at most two coordinates can exceed g (three would need a sum above 9/10 = L), so K = 4 leaves a smooth one
L, g, tau = fmpq(9, 10), fmpq(3, 10), fmpq(4, 5)
caps = [fmpq(0), fmpq(1, 2), fmpq(7, 10)] + [fmpq(7, 10)] * 20
ctx.prec = gt.PREC
sup = gt.Support(K, L, g, tau, caps)
with tempfile.TemporaryDirectory() as tmp:
    U, sup, T = gt.run(d, sup=sup, workers=1, path=f"{tmp}/t.pkl", log=lambda *a: None)
GI, GJ = gt.assemble(U, sup, T)
Lf, gf, tf = float(L.p) / float(L.q), float(g.p) / float(g.q), float(tau.p) / float(tau.q)
capf = np.array([float(c.p) / float(c.q) for c in caps])


def in_T(x):                       # x: array (..., K)
    big = x > gf
    return (x.sum(-1) < Lf) & ((x * big).sum(-1) <= capf[big.sum(-1)])


def basis_fn(a, alpha, x):
    exps = set(itertools.permutations(tuple(alpha) + (0,) * (K - len(alpha))))
    m = sum(np.prod(x ** np.array(e), axis=-1) for e in exps)
    return (Lf - x.sum(-1)) ** a * m


rng = np.random.default_rng(20261115)
dim = len(U.basis)
iu = np.triu_indices(dim)
sI = np.zeros(len(iu[0])); sI2 = np.zeros_like(sI); sJ = np.zeros_like(sI); sJ2 = np.zeros_like(sI)
B = 1_000_000
for _ in range(N):
    x = rng.random((B, K)) * Lf
    ind = in_T(x)
    v = np.stack([basis_fn(a, al, x) * ind for a, al in U.basis])          # (dim, B)
    p = v[iu[0]] * v[iu[1]]
    sI += p.sum(1); sI2 += (p * p).sum(1)
    q = rng.random((B, K - 1)) * Lf
    t1, t2 = rng.random((B, 1)) * Lf, rng.random((B, 1)) * Lf
    x1, x2 = np.hstack([q, t1]), np.hstack([q, t2])
    ok = q.sum(1) <= tf
    v1 = np.stack([basis_fn(a, al, x1) for a, al in U.basis]) * (in_T(x1) & ok)
    v2 = np.stack([basis_fn(a, al, x2) for a, al in U.basis]) * in_T(x2)
    p = 0.5 * (v1[iu[0]] * v2[iu[1]] + v1[iu[1]] * v2[iu[0]])
    sJ += p.sum(1); sJ2 += (p * p).sum(1)
n = N * B
volI, volJ = Lf ** K, Lf ** (K + 1)
print(f"K={K} d={d} dim={dim}, {N} million samples; caps {capf[1:3]}, g={gf}, tau={tf}, L={Lf}")
worst = 0.0
for idx, (i, j) in enumerate(zip(*iu)):
    mI, mJ = sI[idx] / n, sJ[idx] / n
    eI = np.sqrt(max(sI2[idx] / n - mI * mI, 0) / n) * volI
    eJ = np.sqrt(max(sJ2[idx] / n - mJ * mJ, 0) / n) * volJ
    Ie, Je = float(GI[i][j].mid()), float(GJ[i][j].mid())
    zI, zJ = (mI * volI - Ie) / eI, (mJ * volJ - Je) / eJ
    worst = max(worst, abs(zI), abs(zJ))
    if i == j or max(abs(zI), abs(zJ)) > 3:
        print(f"  {U.basis[i]} x {U.basis[j]}: I engine {Ie:.6e} MC {mI*volI:.6e} (rel err {eI/abs(Ie):.1e}, z={zI:+.1f});  "
              f"J engine {Je:.6e} MC {mJ*volJ:.6e} (rel err {eJ/abs(Je):.1e}, z={zJ:+.1f})")
print(f"largest |z| over all {len(iu[0])} I entries and {len(iu[0])} J entries: {worst:.2f}")

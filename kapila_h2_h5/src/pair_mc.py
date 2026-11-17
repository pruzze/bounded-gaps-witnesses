"""Independent plausibility check of Kapila's pair cost for H_2 (hexagon:2610.00419v2, eqs. (8.7)-(8.9)):
C <= 2 zeta^2 J_K(K eta) E_nu[ av_{i<j} W_ij / (1 - zeta x + zeta (y_i + y_j))^2 ],
y_i = t_i/sigma, x = sum y_i, nu = density F^2 / I(F), F = P(S, Z) prod g(t_i) on S <= sigma.
J_K by quadrature of (8.7); the expectation by Monte Carlo (importance sampling from the product law g^2/m2).
The tail restriction of the support is ignored (its mass is of order 1e-7). Only the data file m2.json is read.
Per-sample rows (weight, pair statistic, radial statistic, S) are saved for the error analysis.
Usage: python pair_mc.py <m2.json> <samples> <seed>"""
import json, sys, time
import numpy as np
from scipy import integrate

spec = json.load(open(sys.argv[1])); nsamp = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
k = spec["k"]; theta = spec["theta"]; beta = spec["beta"]; rg = spec["relative_shift"]
E = np.array(spec["E"], dtype=float) * spec["delta"]; psi = np.array(spec["profiles"], dtype=float)
b = k - 1; c = theta / np.log(k); T = beta / np.log(k); ag = rg / (c + b * T)
eps = spec["eps_cells"] * spec["delta"]; sigma = 1 + eps
omega = spec["omega"]; rho = 0.25 + omega
eta = 46843 / 40000                                  # verified_final/m2_final.json
zeta = 2 * rho * sigma                               # limit of s/b for lambda -> rho, b -> 1/2
K = k; a = K * eta; gam = np.sqrt(1 - 2 / K)
# J_K(a): definition (8.7) as a double integral, and the one-dimensional form
f1 = lambda t: 4 * np.arctanh(a * gam * t / (2 + (a - 2) * t)) / (a * gam * (2 + (a - 2) * t))
J1, err1 = integrate.quad(f1, 0, 1, limit=400, epsabs=0, epsrel=1e-12)
def inner(u):
    g = lambda w: 1.0 / ((1 - u - w + a * u) * (1 - u - w + a * w) + a * a / (2 * K) * (u - w) ** 2)
    return integrate.quad(g, 0, 1 - u, limit=200, epsrel=1e-10)[0]
J2 = integrate.quad(inner, 0, 1, limit=400, epsrel=1e-9, points=[1e-6, 1e-4, 1e-2])[0]
print(f"J_K(K eta): one-dimensional form {J1:.10e}, double integral {J2:.10e}   (bundle: J in [8.479277052736967e-06])", flush=True)

# sampling from g^2/m2 by inverse CDF on a fine table
tt = np.concatenate([[0.0], np.geomspace(1e-9, T, 400000)])
prim2 = lambda t: -1 / (b * (c + b * t)) + (2 * ag / b) * np.log(c + b * t) + ag * ag * t
cdf = prim2(tt) - prim2(0.0); cdf /= cdf[-1]
t0 = time.time(); sw = swf = swf2 = sw2 = 0.0; used = 0
rows = []
for n in range(nsamp):
    X = np.interp(rng.random(K), cdf, tt)
    S = X.sum()
    if S >= sigma: continue
    z = X / T
    Z = [1.0, (z ** 2).sum(), (z ** 3).sum(), (z ** 4).sum(), (z ** 2).sum() ** 2]
    cell = np.searchsorted(E, S, side="right") - 1
    P = sum(psi[A, cell] * Z[A] for A in range(5))
    w = P * P
    y = X / sigma; x = y.sum()
    A1 = 1 - x + a * y
    Wm = np.outer(A1, A1) + a * a / (2 * K) * (y[:, None] - y[None, :]) ** 2
    Dm = (1 - zeta * x + zeta * (y[:, None] + y[None, :])) ** 2
    Rm = Wm / Dm; f = (Rm.sum() - np.trace(Rm)) / (K * (K - 1)); rad = ((1 + (eta - 1) * x) / (1 - zeta * x)) ** 2; rows.append((w, f, rad, S))
    sw += w; swf += w * f; swf2 += (w * f) ** 2; sw2 += w * w; used += 1
    if (n + 1) % 5000 == 0 or n + 1 == nsamp:
        np.save(f"pair_mc_rows_{sys.argv[3]}.npy", np.array(rows))
        m = swf / sw
        # delta-method standard error of a ratio estimator (covariance term neglected -> rough)
        se = np.sqrt(max(swf2 / used - (swf / used) ** 2, 0) / used) / (sw / used)
        print(f"[{time.time()-t0:5.0f}s] {n+1} samples ({used} inside): E_nu[...] = {m:.5f} +- ~{se:.5f};  "
              f"C = 2 zeta^2 J E = {2 * zeta**2 * J1 * m:.6e}   (bundle: moment <= 4.20512265610959, C <= 1.9869304490977e-05)", flush=True)

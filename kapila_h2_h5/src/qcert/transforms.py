"""One-coordinate transforms and the weighted Fourier data for Kapila's trials, in Arb ball arithmetic.

Written from the formulas of hexagon:2610.00419v2 (Sections 4 and 6). From the bundle only the trial JSON is read.

  b = k-1, c = theta/log k, T = beta/log k, u = c + bT, a_g = r_g/u,  g(t) = 1/(c+bt) + a_g on [0, T],
  f = g^2/m2 (density of X),  phi(xi) = E e^{i xi X},  mu_p(xi) = E[(X/T)^p e^{i xi X}].
  Phi_W(xi) = E[Z_W e^{i xi S_n}] = sum over set partitions pi of the labelled list W of
              (n)_{|pi|} phi^{n-|pi|} prod_{B in pi} mu_{sum B}   (eq. (4.5) of the paper).
"""
import itertools
import json
import time

from flint import acb, arb, ctx, fmpq

t00 = time.time()


def log(s):
    print(f"[{time.time() - t00:7.1f}s] {s}", flush=True)


def exact(x):
    """The binary64 value x as an exact ball."""
    return arb(fmpq(*float(x).as_integer_ratio()))


class Trial:
    def __init__(self, path):
        spec = json.load(open(path))
        self.spec = spec
        self.k = spec["k"]
        k = self.k
        self.theta, self.beta, self.rg = exact(spec["theta"]), exact(spec["beta"]), exact(spec["relative_shift"])
        self.delta = exact(spec["delta"])
        self.edges = [arb(e) * self.delta for e in spec["E"]]        # dyadic, exact
        self.psi = [[exact(h) for h in row] for row in spec["profiles"]]
        self.basis = [tuple(a) for a in spec["basis"]]
        self.eps = arb(spec["eps_cells"]) * self.delta
        self.sigma, self.tau = 1 + self.eps, 1 - self.eps
        assert self.edges[-1].overlaps(self.sigma)
        assert all(len(r) == len(self.edges) - 1 for r in self.psi)
        lk = arb(k).log()
        self.b = arb(k - 1)
        self.c = self.theta / lk
        self.T = self.beta / lk
        self.u = self.c + self.b * self.T
        self.ag = self.rg / self.u
        b, c, u, ag, T = self.b, self.c, self.u, self.ag, self.T
        self.m2 = (1 / c - 1 / u) / b + 2 * ag / b * (u / c).log() + ag * ag * T

    # ---------- one coordinate ----------
    def phi_closed(self, xi):
        """phi(xi) for real xi != 0 (arb), closed form with E_1 (Section 6.1 of the paper)."""
        b, c, u, ag, T = self.b, self.c, self.u, self.ag, self.T
        a = xi / b
        I = acb(0, 1)
        E = acb(0, -a * c).expint(1) - acb(0, -a * u).expint(1)
        eac = acb(0, -a * c).exp()
        B0 = (acb(0, xi * T).exp() - 1) / (I * xi)
        B1 = eac * E / b
        B2 = eac / b * (acb(0, a * c).exp() / c - acb(0, a * u).exp() / u + I * a * E)
        return (B2 + 2 * ag * B1 + ag * ag * B0) / self.m2

    def moment(self, p, xi, rel_bits=150):
        """mu_p(xi) = int_0^T f(t) (t/T)^p e^{i xi t} dt by Petras integration (integrand meromorphic)."""
        b, c, ag, T = self.b, self.c, self.ag, self.T
        xi = acb(xi)

        def fn(t, _analytic):
            g = 1 / (c + b * t) + ag
            v = g * g * (I * xi * t).exp()
            if p:
                v = v * (t / T) ** p
            return v
        I = acb(0, 1)
        # split near 0, where the integrand varies on the scale c/b
        pts = [arb(0)]
        x = c / b
        while bool(x < T / 4):
            pts.append(x)
            x = x * 8
        pts.append(T)
        tot = acb(0)
        for lo, hi in zip(pts[:-1], pts[1:]):
            tot += acb.integral(fn, lo, hi, rel_tol=arb(2) ** -rel_bits, eval_limit=10 ** 6)
        return tot / self.m2


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        yield [[first]] + part
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]


def weighted(n, W, phi, mu):
    """Phi_W at one frequency from phi and mu (dict p -> value), labelled list W."""
    tot = acb(0)
    if phi.is_zero():
        raise ValueError
    lphi = phi.log()
    for part in set_partitions(list(range(len(W)))):
        r = len(part)
        fall = arb(1)
        for i in range(r):
            fall *= n - i
        term = fall * ((n - r) * lphi).exp()
        for blk in part:
            term *= mu[sum(W[i] for i in blk)]
        tot += term
    return tot

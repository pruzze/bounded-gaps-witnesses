"""Independent check of the segment representation in qcert.py at random real points of random segments:
  kappa_B(s) against the direct sum  sum_e sum_m coef[B][e][m] Wc_e(e_m - s), with Wc_e by Petras integration
  of (t/T)^e g(t) (no closed form, no series), and rho_W(s) against direct evaluation of the truncated Fourier
  series (1/P0)[m_W + 2 Re sum_j Phi_j e^{-i xi_j s}] with every phase from exp.
Also checks the total mass: int_0^tau rho_() ds (Petras over all segments) against C_(){tau} from the CDF formula.
Usage: python check_segments.py <trial.json>"""
import random
import sys

trial = sys.argv[1]
sys.argv = ["qcert.py", trial, "1"]
import qcert as q
from flint import acb, arb

T, b, c, ag = q.T, q.b, q.c, q.ag


def Wc_direct(e, x):
    if bool(x <= 0):
        return arb(0)
    x = arb.min(x, T)
    return acb.integral(lambda t, _: (t / T) ** e * (1 / (c + b * t) + ag), 0, x, rel_tol=arb(2) ** -120).real


def rho_direct(W, s):
    tot = acb(0)
    for j in range(q.J):
        tot += q.Phi[q.k - 1, W][j] * acb(0, -2 * arb.pi() * (j + 1) * s / q.P0).exp()
    return (q.mom[q.k - 1, W] + 2 * tot.real) / q.P0


worst_k, worst_r = 0.0, 0.0
results = {}


def inspect(lo, hi, parts):
    global worst_k, worst_r
    for _ in range(2):
        u = random.random()
        s = lo + (hi - lo) * arb(u)
        s = arb(s.mid())
        kap, rho = parts(acb(s), False)
        for B in q.keys:
            direct = sum((q.coef[B][e][m] * Wc_direct(e, q.edges[m] - s) for e in q.evals
                          for m in range(1, q.ncell + 1) if not q.coef[B][e][m].is_zero()), arb(0))
            if not kap[B].real.overlaps(direct):
                raise SystemExit(f"kappa mismatch B {B}: {kap[B]} vs {direct}")
            rel = float((abs(kap[B].real - direct) / (abs(direct) + 1e-300)).upper())
            worst_k = max(worst_k, rel)
        for W in q.W_num[:4] + q.W_num[-2:]:
            d = rho_direct(W, s)
            if not rho[W].real.overlaps(d):
                raise SystemExit(f"rho mismatch W {W}: {rho[W]} vs {d}")
            worst_r = max(worst_r, float((abs(rho[W].real - d) / (abs(d) + 1e-300)).upper()))
    return acb(0)


random.seed(7)
q.INSPECT = inspect
segs = random.sample(range(1, q.nseg), 12)
for si in segs:
    q.segment_integral(si)
q.log(f"kappa and rho agree on {len(segs)} random segments (2 points each); largest relative difference of "
      f"midpoints: kappa {worst_k:.1e}, rho {worst_r:.1e}")

# mass check with the density alone
q.INSPECT = None
W0 = ()
tot = q.cdf_raw(q.k - 1, W0, q.cdf_num[W0], q.pts[1])          # segment 0 by the CDF


def mass_inspect(lo, hi, parts):
    def f(z, analytic):
        pr = parts(z, analytic)
        return q.NAN if pr is None else pr[1][W0]
    return acb.integral(f, lo, hi, abs_tol=arb(2) ** -130, rel_tol=arb(2) ** -80)


q.INSPECT = mass_inspect
for si in range(1, q.nseg):
    tot += q.segment_integral(si)[1].real
direct = q.cdf_raw(q.k - 1, W0, q.cdf_num[W0], q.tr.tau)
q.log(f"mass check: sum of segment integrals of rho_() = {tot.str(20, radius=True)}; "
      f"C_()(tau) = {direct.str(20, radius=True)}; overlap: {tot.overlaps(direct)}")

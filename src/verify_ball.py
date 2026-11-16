"""Re-check a witness file with ball arithmetic.

Usage: [WORKERS=n] [NOMAT=1] python verify_ball.py <witness.json>

Reads the support (K, L, tau, g, caps C(m)) and the terms [a, alpha, c] of
    F(x) = sum c * (L - |x|_1)^a * m_alpha(x)
from the witness file, recomputes the Gram matrices of I and J on that support for the space B_d (d the degree of
the witness) with Arb balls (gramtab.py), and encloses K * J(F) / I(F) for the integer coefficient vector using
scalar ball operations only. Nothing but the support, the scale and the coefficients is taken from the file: the
ratio stored in it is printed for comparison and not used.

Normalisation: a file without "scale" is in the coordinates x = 4t of the 212 datum and is compared with 4. A file
with "scale" c is in coordinates x = c*t, and K*J/(c*I), the value in the physical scale, is compared with 1.

With NOMAT=1 the tables are built without FLINT matrix products as well (slower; see the README on FLINT 3.6.0).
"""
import json
import os
import sys
import tempfile
import time

from flint import ctx, fmpq

import gramtab as gt
from gram_report import rayleigh

ctx.prec = gt.PREC
fr = lambda s: fmpq(*map(int, s.split("/"))) if "/" in s else fmpq(int(s))
W = json.load(open(sys.argv[1]))
K, d = W["K"], W["degree"]
caps = [fr(c) for c in W["caps_C"]]
sup = gt.Support(K, fr(W["L"]), fr(W["g"]), fr(W["tau"]), [fmpq(0)] + caps + [caps[-1]] * 80)
t0 = time.time()
log = lambda *a: os.write(1, (f"[{time.time()-t0:7.0f}s] " + " ".join(map(str, a)) + "\n").encode())
with tempfile.TemporaryDirectory() as tmp:
    U, sup, T = gt.run(d, sup=sup, workers=int(os.environ.get("WORKERS", "2")), path=f"{tmp}/tables.pkl", log=log)
GI, GJ = gt.assemble(U, sup, T)
log("Gram matrices assembled")
pos = {(a, tuple(al)): i for i, (a, al) in enumerate(U.basis)}
vec = [0] * len(U.basis)
for a, al, c in W["terms"]:
    vec[pos[(a, tuple(al))]] += int(c)
q, I, J = rayleigh(vec, GI, GJ, K)
log(f"K = {K}, degree {d}, {len(W['terms'])} terms, precision {gt.PREC} bits")
if "scale" in W:
    q = q / fr(W["scale"])
    log(f"K*J/(scale*I) recomputed: {q.str(60)}")
    log(f"stored in the file:       {W.get('ratio_enclosure', 'n/a')}")
    log(f"I > 0: {bool(I > 0)}   K*J/(scale*I) > 1: {bool(q > 1)}")
else:
    log(f"K*J/I recomputed:   {q.str(60)}")
    log(f"stored in the file: {W.get('ratio_enclosure', 'n/a')}")
    log(f"I > 0: {bool(I > 0)}   K*J/I > 4: {bool(q > 4)}")
    log(f"in the physical scale (x = 4t) the same quantity is K*J/I = {(q / 4).str(40)}, to be compared with 1")

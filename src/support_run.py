"""Gram tables (gramtab.py) for one of the supports in profiles.py, the optimum over B_d for nested degrees, and,
when d >= 25 and WITNESS1295 names the file, the ratio of the public 1295-term witness (AxiomMath/PrimeGapsLib,
Gap246) on that support.

Usage: python support_run.py <profile> <d> [K]        (K overrides the number of variables, for scouting)
The ratio reported is K * J_T / (c * I_T), c the dilation of the profile (same normalisation as run_profile.py).
"""
import os
import pickle
import sys
import time

from flint import arb, ctx, fmpq

import gramtab as gt
import mono_engine as me
from gram_report import integer_vector, rayleigh
from profiles import PROFILES, witness_coords

name, d = sys.argv[1], int(sys.argv[2])
k, R, g, Cz, r, c = witness_coords(name)
K = int(sys.argv[3]) if len(sys.argv) > 3 else k
tag = f"{name}_k{K}_d{d}"
ctx.prec = gt.PREC
t0 = time.time()
log = lambda *a: os.write(1, (f"[{tag} {time.time()-t0:6.0f}s] " + " ".join(map(str, a)) + "\n").encode())
sup = gt.Support(K, R, g, r, [fmpq(0)] + [Cz(m) for m in range(1, 80)])
U, sup, T = gt.run(d, sup=sup, workers=int(os.environ.get("WORKERS", "2")), path=f"gramtab_{tag}.pkl", log=log)
GI, GJ = gt.assemble(U, sup, T)
log("Gram matrices assembled")
with open(f"gramtab_{tag}_gram.pkl", "wb") as f:
    pickle.dump(dict(basis=U.basis, GI=[[gt.ser(v) for v in row] for row in GI],
                     GJ=[[gt.ser(v) for v in row] for row in GJ], K=K, c=(int(c.p), int(c.q)), prec=gt.PREC), f,
                protocol=pickle.HIGHEST_PROTOCOL)
ca = arb(c)
if d >= 25 and os.environ.get("WITNESS1295"):
    # path to PrimeGapsCert/Gap246/k50e25d25n1295.json of AxiomMath/PrimeGapsLib
    W = me.load_witness(os.environ["WITNESS1295"])
    pos = {b: i for i, b in enumerate(U.basis)}
    vec = [0] * len(U.basis)
    for a, al, cf in W:
        vec[pos[(a, al)]] += cf
    q, I, J = rayleigh(vec, GI, GJ, K)
    log(f"public 1295-term witness: K*J_T/(c*I_T) = {(q / ca).str(45)}")
    log(f"target from the paper: {PROFILES[name]['target']}")
for dd in sorted({dd for dd in (d - 4, d - 2, d) if dd >= 2}):
    idx = [i for i, (a, al) in enumerate(U.basis) if a + sum(al) <= dd]
    sI = [[GI[i][j] for j in idx] for i in idx]
    sJ = [[GJ[i][j] for j in idx] for i in idx]
    q, x = gt.top_ratio(sI, sJ, K)
    cvec = integer_vector([x[i, 0].mid() for i in range(len(idx))])
    qc, _, _ = rayleigh(cvec, sI, sJ, K)
    log(f"optimum over B_{dd} (dim {len(idx)}): K*J_T/(c*I_T) = {(qc / ca).str(35)}")

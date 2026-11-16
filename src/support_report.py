"""Nested optima from a Gram pickle written by support_run.py, and witness files for that support.

Usage: python support_report.py <profile> <gram.pkl> <dmin> <dmax> [witness degrees ...]

For every d' in [dmin, dmax] the best ratio K*J_T/(c*I_T) over B_{d'} is reported for an integer vector (about 50
digits), enclosed with scalar Arb operations. For the listed witness degrees (default: the smallest d' whose ratio
exceeds 1) the vector is written to witness_<profile>_k<K>_d<d'>.json in the witness coordinates of profiles.py.
"""
import json
import sys

from flint import arb, ctx, fmpq

import gramtab as gt
from gram_report import integer_vector, load, rayleigh
from profiles import PROFILES, witness_coords

ctx.prec = gt.PREC
name, path, dmin, dmax = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
want = [int(v) for v in sys.argv[5:]]
basis, GI, GJ, K = load(path)
k, R, g, Cz, r, c = witness_coords(name)
ca = arb(c)
res = {}
for dd in range(dmin, dmax + 1):
    idx = [i for i, (a, al) in enumerate(basis) if a + sum(al) <= dd]
    sI = [[GI[i][j] for j in idx] for i in idx]
    sJ = [[GJ[i][j] for j in idx] for i in idx]
    q, x = gt.top_ratio(sI, sJ, K)
    cvec = integer_vector([x[i, 0].mid() for i in range(len(idx))])
    qc = rayleigh(cvec, sI, sJ, K)[0] / ca
    res[dd] = (idx, cvec, qc)
    print(f"{name} K={K} d={dd:2d} dim={len(idx):4d}: K*J_T/(c*I_T) = {qc.str(35)}   > 1: {bool(qc > 1)}", flush=True)
first = next((dd for dd in sorted(res) if res[dd][2] > 1), None)
for dd in sorted(set(want) | ({first} if first else set())):
    if dd not in res:
        continue
    idx, cvec, qc = res[dd]
    terms = [[basis[i][0], list(basis[i][1]), str(ci)] for i, ci in zip(idx, cvec) if ci != 0]
    caps = [Cz(m) for m in range(1, 13)]
    while len(caps) > 1 and caps[-1] == caps[-2]:
        caps.pop()
    out = dict(
        description="F = sum c * (L - |x|_1)^a * m_alpha(x), m_alpha the monomial symmetric polynomial in K variables, "
                    "on the support T~ = {|x|_1 <= L, sum of coordinates > g is <= C(number of such coordinates)}, "
                    "in the dilated coordinates x = (26/(25U)) t of the profile (profiles.py); marginal with moving "
                    "endpoint as in kintali_J.py. Criterion: K*J_T/(scale*I_T) > 1.",
        profile=name, source=PROFILES[name].get("target", ""), K=K, L=str(fmpq(R)), tau=str(fmpq(r)), g=str(fmpq(g)),
        caps_C=[str(v) for v in caps], scale=str(c), degree=dd, dimension=len(idx), terms=terms,
        ratio_enclosure=qc.str(60), ratio_lower_bound=qc.lower().str(60, radius=False),
        exceeds_1=bool(qc > 1), arb_precision_bits=gt.PREC)
    out_name = f"witness_{name}_k{K}_d{dd}.json"
    json.dump(out, open(out_name, "w"), indent=0)
    print(f"wrote {out_name}: {len(terms)} terms, ratio {qc.str(30)}")

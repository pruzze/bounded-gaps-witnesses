"""Read a Gram pickle from gramtab.py, report the best ratio K*J_T/I_T for every nested degree d' <= d, and write a
witness (integer coefficients) with a rigorously enclosed ratio.

Usage: python gram_report.py gramtab_d21_gram.pkl [second_gram.pkl to cross-check]

All sums that enter the certified ratio are done with scalar Arb operations (no arb_mat products).
"""
import json
import pickle
import sys

from flint import arb, arb_mat, ctx, fmpz

import gramtab as gt

ctx.prec = gt.PREC


def load(path):
    G = pickle.load(open(path, "rb"))
    GI = [[gt.de(v) for v in r] for r in G["GI"]]
    GJ = [[gt.de(v) for v in r] for r in G["GJ"]]
    return G["basis"], GI, GJ, G["K"]


def rayleigh(c, GI, GJ, K):
    """K * c GJ c / c GI c with scalar ball operations; c are exact integers."""
    n = len(c)
    I = J = arb(0)
    for i in range(n):
        if c[i] == 0:
            continue
        si = sj = arb(0)
        for j in range(n):
            if c[j] != 0:
                si += GI[i][j] * c[j]
                sj += GJ[i][j] * c[j]
        I += si * c[i]
        J += sj * c[i]
    return K * J / I, I, J


def integer_vector(x, digits=50):
    """Scale the dyadic vector to integers with about `digits` significant digits in the largest entry."""
    mx = max(abs(v) for v in x)
    scale = arb(10) ** digits / mx
    out = []
    for v in x:
        w = (v * scale).mid()
        out.append(int((w + arb(1) / 2).floor().unique_fmpz()) if True else 0)
    return out


def main():
    basis, GI, GJ, K = load(sys.argv[1])
    d = max(a + sum(al) for a, al in basis)
    if len(sys.argv) > 2:
        b2, GI2, GJ2, _ = load(sys.argv[2])
        assert b2 == basis
        n = len(basis)
        okI = all(GI[i][j].overlaps(GI2[i][j]) for i in range(n) for j in range(n))
        okJ = all(GJ[i][j].overlaps(GJ2[i][j]) for i in range(n) for j in range(n))
        worst = max(float(((GJ[i][j] - GJ2[i][j]) / GJ[i][j]).abs_upper()) for i in range(n) for j in range(n)
                    if not GJ[i][j].contains(0))
        print(f"cross-check against {sys.argv[2]}: I balls overlap: {okI}, J balls overlap: {okJ}, "
              f"largest relative difference in J: {worst:.1e}")
    relrad = max(float((GJ[i][i].rad() / GJ[i][i].abs_lower())) for i in range(len(basis)))
    print(f"dimension {len(basis)}, largest relative radius on the J diagonal: {relrad:.1e}")
    results = {}
    for dd in range(min(14, d), d + 1):
        idx = [i for i, (a, al) in enumerate(basis) if a + sum(al) <= dd]
        sI = [[GI[i][j] for j in idx] for i in idx]
        sJ = [[GJ[i][j] for j in idx] for i in idx]
        q, x = gt.top_ratio(sI, sJ, K)
        xs = [x[i, 0].mid() for i in range(len(idx))]
        c = integer_vector(xs)
        qc, I, J = rayleigh(c, sI, sJ, K)
        results[dd] = (idx, c, qc)
        print(f"d={dd:2d} dim={len(idx):3d}: K*J/I = {qc.str(35)}   > 4: {bool(qc > 4)}", flush=True)
    first = next((dd for dd in sorted(results) if results[dd][2] > 4), None)
    for dd in sorted({d} | ({first} if first else set())):
        idx, c, qc = results[dd]
        terms = [[basis[i][0], list(basis[i][1]), str(ci)] for i, ci in zip(idx, c) if ci != 0]
        out = dict(
            description="Witness for the variational step of Axiom Math's H_1 <= 212 draft (3 Sep 2026): "
                        "F = sum c * (L - |x|_1)^a * m_alpha(x), m_alpha the monomial symmetric polynomial in K variables, "
                        "on the support T~ = {|x|_1 <= L, sum of coordinates > g is <= C(number of such coordinates)}; "
                        "marginal with moving endpoint as in axiom212_v2.py / kintali_J.py.",
            K=K, L="53/50", tau="249/250", g="41/625",
            caps_C=[str(v) for v in gt.axiom_support().caps[1:11]],
            degree=dd, dimension=len(idx), terms=terms,
            ratio_enclosure=qc.str(60), ratio_lower_bound=qc.lower().str(60, radius=False),
            exceeds_4=bool(qc > 4), arb_precision_bits=gt.PREC)
        name = f"witness_axiom212_d{dd}.json"
        json.dump(out, open(name, "w"), indent=0)
        print(f"wrote {name}: {len(terms)} terms, ratio {qc.str(30)}")


if __name__ == "__main__":
    main()

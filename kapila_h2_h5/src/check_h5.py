"""Independent finite check of a modified H_5 certificate (Kapila, hexagon:2610.00419v2, section 12.3) and
selection of the best window for a given pair budget.

Base data from the bundle: h5_residues.u32.gz and the three patch files. On top of that, our own changes:
  <extra patches csv>: prime, missing_residue  (a deleted pair +-p is restored, the missing class of p moves)
  <removals txt>     : signed parent elements that are removed (one per line), may be empty or "-"
Everything is recomputed here: primes up to 2432107459, the parent on both signs, the check that no element lies
in the missing class of any certificate prime (class 0 for the deleted primes: only +-p could lie there), the
windows (n smallest negative, K - n smallest positive entries) and their pair counts.
Usage: python check_h5.py <bundle root> <budget> <extra patches csv or -> <removals txt or -> [out prefix]"""
import sys, time, gzip, csv, hashlib
import numpy as np

t0 = time.time(); log = lambda s: print(f"[{time.time()-t0:6.0f}s] {s}", flush=True)
root, budget, EXTRA, REMOVE = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
OUT = sys.argv[5] if len(sys.argv) > 5 else None
CUT, UP, L, K = 10_000_000, 234_360_000, 2_432_107_459, 232_465_000
odd = np.ones((L + 1) // 2, dtype=bool); odd[0] = False
i = 1
while (2 * i + 1) ** 2 <= L:
    if odd[i]: p = 2 * i + 1; odd[(p * p) // 2::p] = False
    i += 1
log(f"sieve done: pi({L}) = {int(odd.sum()) + 1}")
cp = 2 * np.flatnonzero(odd[: (UP + 1) // 2]) + 1; cp = cp[(cp > CUT) & (cp <= UP)].astype(np.int64)
res = np.frombuffer(gzip.decompress(open(f"{root}/baseline/inherited/tuples/h5_residues.u32.gz", "rb").read()), dtype="<u4").astype(np.int64).copy()
assert len(cp) == len(res) == 12_204_094
for fn in ("baseline/h5/h5_new_patches.csv", "new/h5_dependent_patches.csv", "new/h5_dependent2_patches.csv"):
    for row in csv.DictReader(open(f"{root}/{fn}")):
        j = int(np.searchsorted(cp, int(row["prime"]))); assert res[j] == 0; res[j] = int(row["missing_residue"])
pairs = np.empty((len(cp), 2), dtype="<u4"); pairs[:, 0] = cp; pairs[:, 1] = res
digest = hashlib.sha256(pairs.tobytes()).hexdigest(); del pairs
log(f"paper's certificate: zero residues {int((res == 0).sum())}; sha256 matches section 12.3: {digest == 'f223f9823aef305e4eeeec4f64c933ba20e1dd42872408c557e5a3011cd4c191'}")
nex = 0
if EXTRA != "-":
    for row in csv.DictReader(open(EXTRA)):
        p, r = int(row["prime"]), int(row["missing_residue"]); j = int(np.searchsorted(cp, p))
        assert cp[j] == p and res[j] == 0 and 0 < r < p, f"bad extra patch {p}"
        res[j] = r; nex += 1
log(f"extra patches: {nex} deleted pairs restored; zero residues now {int((res == 0).sum())}")
odd[: (CUT + 1) // 2] = False; odd[cp[res == 0] >> 1] = False
pos = odd; neg = odd.copy()                                       # pos[i]: +(2i+1) in the parent; neg[i]: -(2i+1)
nrm = 0
if REMOVE != "-":
    for line in open(REMOVE):
        x = int(line.split()[0])
        arr = pos if x > 0 else neg
        assert arr[abs(x) >> 1], f"removed element {x} is not in the parent"
        arr[abs(x) >> 1] = False; nrm += 1
log(f"removed elements: {nrm}; parent: {int(pos.sum())} positive, {int(neg.sum())} negative entries")

# every certificate prime with a nonzero class r: no +q = r (mod p) and no -q = r (mod p), i.e. q = p - r
nz = np.flatnonzero(res != 0); fails = 0; looked = 0; CH = 40000
for a in range(0, len(nz), CH):
    idx = nz[a:a + CH]; p = cp[idx]; r = res[idx]
    J = np.arange(L // int(p[0]) + 1, dtype=np.int64)
    for t, arr in ((r, pos), (p - r, neg)):
        q = t[:, None] + p[:, None] * J[None, :]
        ok = (q <= L) & (q & 1 == 1)
        hit = arr[q[ok] >> 1]; looked += hit.size; fails += int(hit.sum())
# deleted primes (class 0): the only possible elements are +-p
zp = cp[res == 0]; fails0 = int(pos[zp >> 1].sum() + neg[zp >> 1].sum())
log(f"nonzero classes checked: {len(nz)} primes, {looked} lookups, failures {fails}; deleted primes present: {fails0}")
log(f"ADMISSIBLE (class 0 for p <= {CUT}, recorded class for {CUT} < p <= {UP}): {fails == 0 and fails0 == 0}")

qp = (2 * np.flatnonzero(pos) + 1).astype(np.int64); qn = (2 * np.flatnonzero(neg) + 1).astype(np.int64)
def long_pairs(nneg, npos, D):
    bb = qp[:npos]; aa = qn[int(np.searchsorted(qn[:nneg], D - bb[-1], side="right")):nneg]
    return int(np.sum(npos - np.searchsorted(bb, D - aa, side="right"), dtype=np.int64))
def best_D(nneg, npos):
    hi_ = int(qn[nneg - 1] + qp[npos - 1]); lo_ = hi_ - 40_000_000
    assert long_pairs(nneg, npos, lo_) > budget
    while lo_ < hi_:
        mid = (lo_ + hi_) // 2
        if long_pairs(nneg, npos, mid) <= budget: hi_ = mid
        else: lo_ = mid + 1
    return lo_ + (lo_ % 2)
nmin, nmax = max(1, K - len(qp)), min(len(qn), K - 1)
log(f"budget {budget}: {nmax - nmin + 1} windows (negative count from {nmin} to {nmax})")
best = None; step = max(1, (nmax - nmin) // 80); cands = list(range(nmin, nmax + 1, step)) + [nmax]
while True:
    for n in cands:
        v = (best_D(n, K - n), n)
        if best is None or v < best: best = v
    if step == 1: break
    step = max(1, step // 20)
    cands = [n for n in range(best[1] - 25 * step, best[1] + 25 * step + 1, step) if nmin <= n <= nmax]
Dn, n = best
log(f"RESULT: {n} negative and {K - n} positive entries: first {-int(qn[n - 1])}, last {int(qp[K - n - 1])}, "
    f"diameter {int(qn[n - 1] + qp[K - n - 1])}; pairs further apart than {Dn}: {long_pairs(n, K - n, Dn)} <= {budget}; "
    f"beyond {Dn - 2}: {long_pairs(n, K - n, Dn - 2)}")
if OUT:
    open(OUT + "_window.txt", "w").write(f"negative_count {n}\npositive_count {K - n}\nfirst {-int(qn[n - 1])}\nlast {int(qp[K - n - 1])}\nD {Dn}\nlong_pairs {long_pairs(n, K - n, Dn)}\nbudget {budget}\n")

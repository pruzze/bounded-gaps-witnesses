"""Exact rational value of k*J_T/I_T for a witness on a one-band support profile (profiles.py).
Usage: WITNESS=<witness file> [WORKERS=n] python run_profile.py <profile>
Writes result_<profile>.txt (the exact fraction) and checkpoints progress_<profile>.jsonl in the working directory."""
import sys, os, time, pickle
from collections import defaultdict
from multiprocessing import get_context
from flint import fmpq
import mono_engine as me, capped_mono as cm, axiom212_v2 as v2
import kintali_J as kj
from profiles import PROFILES, witness_coords

name = sys.argv[1]
k, R, g, Cz, r, c = witness_coords(name)
caps = tuple([fmpq(0)] + [Cz(m) for m in range(1, 80)])
W = me.load_witness(os.environ["WITNESS"])
# patch kintali_J module globals (used by fork-inherited workers)
kj.k, kj.R, kj.g, kj.r, kj.W = k, R, g, r, W
kj.PR = v2.Params(K=k, L=R, g=g, tau=r, caps=list(caps))
t0 = time.time()
# os.write, not print: the pool re-forks workers (maxtasksperchild), and a fork taken while the main thread
# holds the buffered-stdout lock leaves the child hanging at its first flush (seen 2026-11-15).
log = lambda *a: os.write(1, (f"[{name} {time.time()-t0:6.0f}s] " + " ".join(map(str, a)) + "\n").encode())
assert R - r <= g, "rectangle shortcut for cell 2 needs 2*eps <= delta"
for m in range(1, 40):
    assert caps[m] <= caps[m + 1] <= caps[m] + g, "support shape: C nondecreasing, C(m+1) <= C(m) + delta"

I_free, J_free = me.IJ(W, k, R, r)
log(f"cap-free k*J/I = {float(k * J_free / (c * I_free)):.15f}")

cache = f"cache_{name}_sq.pkl"
if os.path.exists(cache):
    sqI, sqV, sqg, sqW = pickle.load(open(cache, "rb"))
else:
    sqI = cm.square(W, k); log("square F^2 done")
    sqV = kj.square_G(kj.G_terms("V")); log("square G(V) done")
    sqg = kj.square_G(kj.G_terms("g")); log("square G(g) done")
    sqW = kj.square_GW(kj.G_terms_W()); log("square G(W) done")
    pickle.dump((sqI, sqV, sqg, sqW), open(cache, "wb"))

# Memory: the big tables are module globals inherited by fork (copy-on-write); tasks carry only (s, t).
SQ = dict(I=sqI, V=sqV, g=sqg, W=sqW)

def jobI(st):
    s, t = st
    v = cm.I_bad_stratum(SQ["I"], k, R, g, caps, s, t)
    return s, t, (int(v.p), int(v.q))

def jobJ(st):
    s, t = st
    return kj.J_stratum_all((SQ["V"], SQ["g"], SQ["W"], s, t))

import json
ckpt = f"progress_{name}.jsonl"
done = {}
if os.path.exists(ckpt):
    for line in open(ckpt):
        d = json.loads(line); done[(d["kind"], d["s"], d["t"])] = d
def save(kind, s, t, p, q, pc=None):
    with open(ckpt, "a") as f:
        f.write(json.dumps(dict(kind=kind, s=s, t=t, p=str(p), q=str(q), pc={c: [str(a), str(b)] for c, (a, b) in (pc or {}).items()})) + "\n")
workers = int(os.environ.get("WORKERS", "8"))
with get_context("fork").Pool(workers, maxtasksperchild=4) as pool:
    Ib = fmpq(0)
    todoI = []
    for st in cm.bad_strata(k, R, g, caps):
        d = done.get(("I", *st))
        if d: Ib += fmpq(int(d["p"]), int(d["q"]))
        else: todoI.append(st)
    log(f"I strata: {len(todoI)} to do, {len(done)} checkpointed")
    for s, t, (p, q) in pool.imap_unordered(jobI, todoI):
        Ib += fmpq(p, q); save("I", s, t, p, q)
        log(f"  I stratum s={s} t={t}")
    log(f"I_bad/I_free = {float(Ib / I_free):.6e}")
    J = fmpq(0); percell = defaultdict(lambda: fmpq(0))
    todoJ = []
    for st in v2.J_strata(kj.PR):
        d = done.get(("J", *st))
        if d:
            J += fmpq(int(d["p"]), int(d["q"]))
            for cid, (a, b) in d["pc"].items(): percell[int(cid)] += fmpq(int(a), int(b))
        else: todoJ.append(st)
    for s, t, (p, q), pc in pool.imap_unordered(jobJ, todoJ):
        save("J", s, t, p, q, pc)
        log(f"  J stratum s={s} t={t}")
        J += fmpq(p, q)
        for cid, (pp, qq) in pc.items():
            percell[cid] += fmpq(pp, qq)
    log("J_T done")
I_T = I_free - Ib
q = k * J / (c * I_T)
digits = q.p * 10 ** 40 // q.q
thr = 4 if PROFILES[name].get("direct") else 1      # dilated profiles are stated with threshold 4
log(f"capped k*J_T/I_T = {str(digits)[0]}.{str(digits)[1:]}   (> {thr}: {q > thr})")
log(f"target: {PROFILES[name]['target']}")
pickle.dump(dict(I_free=I_free, J_free=J_free, I_bad=Ib, J_T=J, percell=dict(percell), ratio=q), open(f"result_{name}.pkl", "wb"))
open(f"result_{name}.txt", "w").write(f"{q}\n")

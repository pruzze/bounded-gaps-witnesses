"""Kintali's capped one-band D182 support (k=48): exact I_T, and later J_T, for the public 1295-term witness."""
import sys, time, pickle, os
from flint import fmpq
import mono_engine as me, capped_mono as cm

# Default witness: the public 1295-term witness of AxiomMath/PrimeGapsLib (PrimeGapsCert/Gap246), if a path to it
# is given in WITNESS. run_profile.py replaces W, k, R, g, r and the caps by those of the chosen profile.
_w = os.environ.get("WITNESS", "")
W = me.load_witness(_w) if os.path.exists(_w) else []
U, Lb, delta, Bt = fmpq(521, 2000), fmpq(491, 2000), fmpq(91, 5000), fmpq(1899, 10000)
c = fmpq(26, 25) / U                     # z = c t  (witness coordinates)
k, R, g, B, r = 48, fmpq(26, 25), c * delta, c * Bt, c * Lb
if __name__ == "__main__":
    t = time.time()
    if os.path.exists("kintali_square.pkl"):
        sq = pickle.load(open("kintali_square.pkl", "rb"))
    else:
        sq = cm.square(W, k)
        pickle.dump(sq, open("kintali_square.pkl", "wb"))
    nterms = sum(len(v) for v in sq.values())
    print(f"square: {len(sq)} signatures nu, {nterms} (nu,A) terms [{time.time()-t:.1f}s]", flush=True)
    st = cm.bad_strata(k, R, g, B)
    print(f"bad strata: {len(st)}: {st[:6]} ... {st[-3:]}", flush=True)
    if len(sys.argv) > 1 and sys.argv[1] == "time1":
        t = time.time(); s0, t0 = st[0]
        v = cm.I_bad_stratum(sq, k, R, g, B, s0, t0); print(f"one stratum {st[0]}: {time.time()-t:.1f}s  value {float(v):.4e}")
        sys.exit()
    I_free, J_free = me.IJ(W, k, R, r)
    Ib = cm.I_bad(sq, k, R, g, B, verbose=True)
    I_T = I_free - Ib
    pickle.dump(dict(I_free=I_free, J_free=J_free, I_bad=Ib, I_T=I_T), open("kintali_I.pkl", "wb"))
    print(f"I_free={float(I_free):.15e}  I_bad={float(Ib):.6e}  rel={float(Ib/I_free):.6e}")
    print(f"ratio using J_free and I_T: {float(48*J_free/(c*I_T)):.15f}  (cap-free q0 = {float(48*J_free/(c*I_free)):.15f})")

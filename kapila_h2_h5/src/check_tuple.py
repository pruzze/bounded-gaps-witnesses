"""Check an admissible K-tuple and its number of long pairs.

The tuple file holds one integer per line (plain text, or gzip when the name ends in .gz). The script checks
  - that the entries are distinct (it sorts them),
  - admissibility: for every prime p <= K some residue class mod p contains no entry. Without a witness file this
    is checked by reducing every entry mod every prime. With a witness file (CSV "prime,missing_residue", plain or
    gzip) the stated class of every prime p <= K is checked to be empty, by a bitmap lookup along the progression;
  - the number of pairs i < j with h_j - h_i > D.

Usage: python check_tuple.py <tuple file> <D> [witness csv]
Requires numpy.
"""
import gzip
import sys
import time

import numpy as np


def read_ints(path, skip=0):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        lines = f.read().split("\n")[skip:]
    return lines


t0 = time.time()
path, D = sys.argv[1], int(sys.argv[2])
wpath = sys.argv[3] if len(sys.argv) > 3 else None
T = np.array(sorted(int(x) for x in read_ints(path) if x.strip()), dtype=np.int64)
K = len(T)
assert np.all(np.diff(T) > 0), "repeated entries"
lo, hi = int(T[0]), int(T[-1])
print(f"K = {K}, first {lo}, last {hi}, diameter {hi - lo}")

sieve = np.ones(K + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(K ** 0.5) + 1):
    if sieve[i]: sieve[i * i::i] = False
primes = np.flatnonzero(sieve).tolist()

if wpath is None:
    bad = [p for p in primes if len(np.unique(T % p)) == p]
    print(f"primes <= K: {len(primes)}; primes with every class occupied: {len(bad)} {bad[:10]}")
    ok = not bad
else:
    w = {}
    for line in read_ints(wpath, skip=1):
        if line.strip():
            p, r = line.split(","); w[int(p)] = int(r)
    member = np.zeros(hi - lo + 1, dtype=bool); member[T - lo] = True
    missing = [p for p in primes if p not in w]
    bad = [p for p in primes if p in w and member[(w[p] % p - lo) % p::p].any()]
    print(f"primes <= K: {len(primes)}; without a witness: {len(missing)}; witness class occupied: {len(bad)}")
    ok = not missing and not bad
print(f"admissible: {ok}")
n = int(np.sum(K - np.searchsorted(T, T + D, side="right")))
print(f"pairs further apart than {D}: {n}   [{time.time() - t0:.0f} s]")

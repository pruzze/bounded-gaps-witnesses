"""Physical support parameters (one band [0,U), marginal base |x| <= Lb, rough threshold delta, caps C(m))."""
from flint import fmpq

PROFILES = {
    "smoketest": dict(k=8, U=fmpq(1), Lb=fmpq(19, 20), delta=fmpq(3, 20), caps=lambda m: fmpq(0) if m <= 0 else (fmpq(2, 5) if m <= 2 else fmpq(1, 2)), target="n/a"),
    # Kintali, "Improved bounded gaps between primes" (1 Sep 2026), D182 support: k=48, uniform cap 0.1899
    "kintali": dict(k=48, U=fmpq(521, 2000), Lb=fmpq(491, 2000), delta=fmpq(91, 5000),
                    caps=lambda m: fmpq(0) if m <= 0 else fmpq(1899, 10000),
                    target="1.000774701187468444763475111497864"),
    # Stadlmann arXiv 2608.31126 v1 displayed choice (as quoted by Kintali): k=49, eps=0.0075, A=0.253,
    # delta=0.028, caps 0.15 for one or two large coordinates and 0.17 thereafter
    "stadlmann": dict(k=49, U=fmpq(2605, 10000), Lb=fmpq(2455, 10000), delta=fmpq(28, 1000),
                      caps=lambda m: fmpq(0) if m <= 0 else (fmpq(15, 100) if m <= 2 else fmpq(17, 100)),
                      target="> 1 (paper states 49J/I > 1, no value given)"),
    # Stadlmann's displayed support with one variable fewer (k = 48, H(48) = 236); witness of degree 22
    "stadlmann48": dict(k=48, U=fmpq(2605, 10000), Lb=fmpq(2455, 10000), delta=fmpq(28, 1000),
                        caps=lambda m: fmpq(0) if m <= 0 else (fmpq(15, 100) if m <= 2 else fmpq(17, 100)),
                        target="1.0001616319193665101960411879539314036267 for witnesses/stadlmann_k48_d22.json "
                               "(exact rational value, also enclosed by Arb)"),
    # Axiom Math draft (3 Sep 2026), Table 3, already in the dilated coordinates of its Section 9 (factor 4):
    # K = 45, L = 53/50, marginal base 249/250, g = 41/625, ten caps C(m). For our own witnesses in B_d.
    "axiom": dict(direct=True, k=45, R=fmpq(53, 50), r=fmpq(249, 250), g=fmpq(41, 625),
                  caps=lambda m: fmpq(0) if m <= 0 else 4 * fmpq(
                      (777, 794, 875, 917, 953, 983, 1016, 1042, 1063, 1081)[min(m, 10) - 1], 5000),
                  target="Axiom Theorem 11.1: 4.00438409833460131937... for their P*; ours (Arb): "
                         "d=21 4.00438409833460132024933958, d=20 4.00240137702709155691567799"),
}


def witness_coords(name):
    """Witness coordinates z = (26/(25U)) t: returns k, R, g(m-independent), Cz(m), r, scale c."""
    p = PROFILES[name]
    if p.get("direct"):
        return p["k"], p["R"], p["g"], p["caps"], p["r"], fmpq(1)
    c = fmpq(26, 25) / p["U"]
    Cz = lambda m, p=p, c=c: c * p["caps"](m)
    return p["k"], fmpq(26, 25), c * p["delta"], Cz, c * p["Lb"], c

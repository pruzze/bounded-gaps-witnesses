# Smaller H_2, ..., H_5 from the certificates of hexagon:2610.00419v2

A. Kapila, "Finite prime gaps from weighted Selberg sieves and pair penalties"
([hexagon:2610.00419v2](https://hexagonmath.org/2610.00419v2), 6 November 2026, marked by its author as primarily
AI-generated) proves H_2 <= 19 844, H_3 <= 1 357 780, H_4 <= 82 234 788 and H_5 <= 4 830 355 416, where
H_m = liminf (p_{n+m} - p_n). The final step (Section 11 of the paper) is the inequality

    rho * (K/k) * q_-  -  m (1 - beta_-)  -  B_D(H) * C_+  >  0,

where q_-, beta_- and C_+ come from numerical certificates and B_D(H) is the number of pairs of the admissible
K-tuple H that are further apart than D. The resulting bound is H_m <= D.

The certificates are not changed here. The tuples are. With the paper's own q_-, beta_-, C_+ and pair budgets:

| | paper | here | data |
|---|---|---|---|
| H_2 | 19 844 | **19 842** | `data/h2_k2305_D19842.txt` |
| H_3 | 1 357 780 | **1 357 700** | `data/h3_k110070_D1357700.txt` |
| H_4 | 82 234 788 | **82 205 490** | `data/h4_k5089150_D82205490.txt.gz`, `data/h4_witness.csv.gz` |
| H_5 | 4 830 355 416 | **4 827 554 952** | `data/h5_restored_pairs.csv`, `data/h5_window_D4827554952.txt` |

**Everything here is conditional on the paper.** These are the paper's theorems applied to other tuples. If the
paper's analysis or its certified values fail, so do these numbers. Section "What has been checked" below says
which parts of the paper were reproduced independently and which were not.

## The tuples

| | K | long pairs | budget | where the budget comes from |
|---|---|---|---|---|
| H_2 | 2305 | 7 beyond 19 842 | 7 | Section 11.2 (room for 7.04 pairs) |
| H_3 | 110 070 | 3 637 beyond 1 357 700 | 3 639 | Section 11.1 table (the printed values allow 3 640) |
| H_4 | 5 089 150 | 4 076 272 beyond 82 205 490 | 4 076 326 | `verified_final/m4_final.json`, `pair_budget: 4076327` |
| H_5 | 232 465 000 | 4 672 314 934 beyond 4 827 554 952 | 4 672 314 960 | `verified_final/m5_final.json`, `pair_budget: 4672314961` |

The budgets are floor((surplus)/C_+) from the printed rational endpoints of Section 11.1; for H_4 and H_5 the
certificate files give the same number plus one, and the smaller value is used here.

**H_2.** An admissible 2305-tuple from -14 to 19 860 with exactly 7 pairs further apart than 19 842:
(-14, 19840), (-14, 19846), (-14, 19860), (0, 19846), (0, 19860), (6, 19860), (16, 19860). It was found by an
integer program over residue-class switches with the pair count as a constraint.

**H_3.** The paper's 110 070-tuple is not maximal inside its own window: 4 positions are blocked by no prime.
Filling them and a local search over residue-class switches gives 110 083 points, and a run of 110 070 of them has
3 637 pairs further apart than 1 357 700.

**H_4.** A window position outside the paper's tuple can be added whenever every prime whose stated missing class
contains it has a second empty class (the stated class of that prime then moves). With the primes up to 250 000
kept fixed, 1 311 positions can be added; switches of the empty class of tight primes in (250 000, 800 000) and
repeated filling add some 470 more. A run of 5 089 150 of the resulting points is in the data file, with a
missing class for every prime up to K in `h4_witness.csv.gz`.

**H_5.** Two gains. (a) The selected window of the paper has 3 267 409 565 long pairs, but its final certificate
allows 4 672 314 961; `radial_pair/m5_d35000.json` of the paper's bundle still carries an older budget. The full
budget alone gives D = 4 830 013 236. (b) The paper's parent set {+-q : q prime, 10^7 < q <= 2 432 107 459}
deletes the pair +-p for 664 149 primes p whose missing class is 0. For 56 841 of the deleted primes in
[2*10^7, 2.25*10^7) some nonzero class contains no element of the parent; moving the missing class of p there
puts +-p back. `data/h5_restored_pairs.csv` lists these primes with their new missing class. The window
`data/h5_window_D4827554952.txt` takes the 116 856 576 smallest negative and the 115 608 424 smallest positive
entries of the modified parent.

## Checking

Requirements: Python 3 with numpy; scipy for the pair-cost checks. The H_5 check needs about 4 GB of memory.

    cd src
    python check_tuple.py ../data/h2_k2305_D19842.txt 19842
    python check_tuple.py ../data/h3_k110070_D1357700.txt 1357700
    python check_tuple.py ../data/h4_k5089150_D82205490.txt.gz 82205490 ../data/h4_witness.csv.gz

`check_tuple.py` checks distinctness, admissibility for every prime up to K (by reducing every entry, or by
checking each stated missing class against a bitmap) and counts the long pairs.

H_5 is too large to store as a list. It is defined by the paper's certificate plus `h5_restored_pairs.csv`.
Download the paper's bundle `prime_gap_verified_certificates.zip` from
[a-kapila/maths-proofs](https://github.com/a-kapila/maths-proofs/tree/main/prime-gaps)
(sha256 `ea5fd3c99270f1d9f93b2a99c6ffc5e97fc5da7b581ae00e6f6743aa317564b4`), unpack it, and run

    python check_h5.py <unpacked>/prime_gap_verified/prime_gap_boundary_bundle 4672314960 ../data/h5_restored_pairs.csv -

It sieves the primes up to 2 432 107 459 itself, rebuilds the paper's certificate from its residue file and patch
files (and checks the sha256 printed in Section 12.3 of the paper), applies the restored pairs, checks that no
element of the parent lies in the missing class of any certificate prime, and selects the best window for the
given budget. Only data files of the bundle are read; none of its code is run. Runtime about 4 minutes.

## What has been checked

All with programs written from the paper's formulas; from the paper's bundle only data files were read.

- **Finite data of the paper**, all four cases: the four tuples are admissible and have the stated pair counts;
  for H_5 the parent rebuilt with an own sieve has 233 950 978 elements, none in any of the 11 539 945 stated
  nonzero classes, and the sha256 of the prime/residue sequence is the one printed in Section 12.3.
- **The quotients**, all four cases (`q_fourier.py`, double precision, Fourier inversion):

  | m | k | here | bundle, `Q_exploratory` | certified lower endpoint |
  |---|---|---|---|---|
  | 2 | 2305 | 7.620007352530 | 7.62000735252986 | 7.6199897109 |
  | 3 | 110070 | 11.428912385 | (not stored) | 11.4288832314 |
  | 4 | 5090000 | 15.238803150 | 15.2388031515 | 15.2385787059 |
  | 5 | 232500000 | 19.048129063 | 19.0481293149 | 19.0480476275 |

- **The pair costs.** H_2 (`pair_mc.py`): J_K(K eta) = 8.4792770527e-06 by quadrature (bundle 8.479277052737e-06);
  the expectation in (8.9) by Monte Carlo over 146 057 samples is 4.2030 +- 0.0012 against the certified upper
  endpoint 4.20512. H_4, H_5 (`pair_radial.py`, the radial bound (10.3)): E[R(S/sigma)^2] = 4.44273 and 4.57946
  against the certified 4.44466 and 4.58080; J_K agrees with the bundle to 10 digits; beta_- from (10.2) equals
  the printed values.
- **The terminal inequalities** with exact fractions from the printed endpoints of Section 11.1.
- **The distribution input.** Theorem 2.3 of the paper (from the paper's reference [2]) is proved in the Lean
  development [openai/PrimeGaps186](https://github.com/openai/PrimeGaps186) (commit `61340d0`, theorem
  `source_mpz3_primeIndicator_coherent_divisor_weight_log_saving_of_deligne`) in the generality used, given two
  Deligne-type bounds that are axioms there. The statements were compared, not every definition.

**Not checked:** the rigorous interval enclosures of the paper (its Appendix B; the values above are floating point
and only show that the certified endpoints are plausible), the support-loss bounds of Section 5, the four-moment
majorants of Section 9, the slicing argument of Section 10, and the pair costs for H_3. Sections 2, 3, 7 and 8 of
the paper were read and no error was found, which is a reading, not a verification.

## Files

- `data/`: the tuples, the H_4 witness file, the H_5 patches and window.
- `src/check_tuple.py`, `src/check_h5.py`: the finite checks.
- `src/q_fourier.py`, `src/pair_mc.py`, `src/pair_radial.py`: the floating-point checks of the certified values;
  each reads a trial file (`mpfr_runs/m2.json` etc.) from the paper's bundle.

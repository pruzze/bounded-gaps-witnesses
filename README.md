# Witnesses for the variational step in bounded gaps between primes

Recent bounds for H_1 = liminf (p_{n+1} - p_n) in the framework of Stadlmann (arXiv:2608.31126) end with a
finite computation: one needs a symmetric function F on an explicit support T_k with

    k * J(F) / I(F) > 1.

This repository contains two such functions that were computed independently of the papers, the exact rational
value of the ratio for each, and the programs that produced and checked them.

| witness | for | k | degree | terms | k J / I (physical scale) |
|---|---|---|---|---|---|
| `witnesses/axiom212_k45_d21.json` | the certificate of Axiom Math's H_1 <= 212 | 45 | 21 | 846 | 1.00109602458365033006... (exact) |
| `witnesses/axiom212_k45_d20.json` | the same | 45 | 20 | 707 | 1.00060034425677288922... (ball arithmetic) |
| `witnesses/stadlmann_k48_d22.json` | Stadlmann's displayed support with k = 48 | 48 | 22 | 1041 | 1.00016163191936651019... (exact) |

Neither is a new bound. The first fills in a hypothesis of a published formalization; the second is an
observation about a published parameter choice. Both rest entirely on the analytic results of the papers, which
are not examined here.

A third, separate part is in [`kapila_h2_h5/`](kapila_h2_h5/): admissible tuples that, with the unchanged
certificates of A. Kapila's hexagon:2610.00419v2, give H_2 <= 19 842, H_3 <= 1 357 700, H_4 <= 82 205 490 and
H_5 <= 4 827 554 952 (the paper: 19 844, 1 357 780, 82 234 788, 4 830 355 416). With the paper's quotient
bounds recomputed rigorously in ball arithmetic (they come out larger than the certified ones), the same
certificates give H_3 <= 1 357 652, H_4 <= 82 190 538 and H_5 <= 4 827 255 418. All of these are conditional on
that paper; its support-loss bounds, pair costs and parts of the analysis were not checked. See the README in
that folder.

## 1. The numerical certificate of H_1 <= 212

Axiom Math's Lean formalization of H_1 <= 212
([AxiomMath/BGP212](https://github.com/AxiomMath/BGP212), read at commit `6641780`) is conditional on a
numerical certificate, `Gap212.Gap212Certificate`: the existence of a symmetric, square-integrable function F
supported on the set T_45 of the parameter datum `gap212Params` with

    0 < I(F) < 45 * J(F).

At the time of writing, no such F is published in that repository or in the draft paper
("A New Bound for Small Gaps Between Primes", Axiom Math, preliminary draft of 3 September 2026), whose
Theorem 11.1 states the value of 45 J / I for a vector P* that is not printed. The F given here is not the
authors' F.

In the dilated coordinates x = 4t of Section 9.3 of the draft (support parameters L = 53/50, tau = 249/250,
g = 41/625, caps C(m) = 4 B_m),

    F~(x) = sum over (a, alpha) of c * (L - x_1 - ... - x_45)^a * m_alpha(x)     on the support, 0 outside,

where m_alpha is the monomial symmetric polynomial of the partition alpha in 45 variables, alpha has even parts
only, and a + |alpha| <= d. These functions span the space B_d of the draft. In the physical scale of the Lean
file, F(t) = F~(4t); then I(F) = 4^-45 I(F~) and J(F) = 4^-46 J(F~), so `Iint < 45 * Jint` is the same statement
as 45 J(F~) / I(F~) > 4.

| file | degree d | terms | 45 J(F~)/I(F~) | 45 J(F)/I(F) | how it is known |
|---|---|---|---|---|---|
| `witnesses/axiom212_k45_d21.json` | 21 | 846 | 4.00438409833460132024933958040890879577... | 1.00109602458365033006233489510222719894... | exact rational arithmetic, and ball arithmetic |
| `witnesses/axiom212_k45_d20.json` | 20 | 707 | 4.00240137702709155691567799325574511554... | 1.00060034425677288922891949831393627888... | ball arithmetic only |

The coefficients are integers of at most 51 digits. `witnesses/axiom212_k45_d21_exact_ratio.txt` holds the
dilated value 45 J(F~)/I(F~) = 4.00438... for d = 21 as a fraction (numerator 387 digits, denominator 386 digits);
the physical value 45 J(F)/I(F) is that fraction divided by 4.

The draft states 4.00438409833460131937... for its P*. The value above agrees with it to 17 significant digits and
is slightly larger, as expected if the vector here is the maximiser over B_21 and P* a rounded maximiser.

The definitions were read against `Gap212Challenge/Basic.lean` (`T`, `Iint`, `Jint`, `gap212Params`,
`Certificate`).

## 2. Stadlmann's displayed support already gives k = 48

Stadlmann (arXiv:2608.31126v1) proves H_1 <= 240 with k = 49 variables on the one-band support with
eps = 0.0075, A_1 = 0.253, delta = 0.028 and caps B_{1,1} = B_{1,2} = 0.15, B_{1,m} = 0.17 for m >= 3.
On our reading of the paper, k enters the argument only through the final variational step: the conditions on
the support do not involve k.

With these parameters left unchanged, the variational step also goes through for k = 48:

    48 * J(F) / I(F) = 1.00016163191936651019604118795393140362670...

for the function F of degree 22 with 1041 terms in `witnesses/stadlmann_k48_d22.json`. The value is an explicit
fraction (`witnesses/stadlmann_k48_d22_exact_ratio.txt`, numerator and denominator of 420 digits each; this file
is in the physical scale, unlike the one for the 212 witness), computed
in exact rational arithmetic and enclosed by ball arithmetic; the two agree on all 60 digits compared. Since the
narrowest admissible 48-tuple has diameter H(48) = 236, her theorem with her own parameters yields H_1 <= 236.

The margin is small. The best ratio over B_d for this support is 0.99953 at d = 21, 1.00016 at d = 22,
1.00062 at d = 23 and 1.00124 at d = 25 (ball arithmetic), so degree 22 is the lowest that works.

The bound 236 itself is not new: Kintali ("Improved bounded gaps between primes", 1 September 2026) obtained it
first, with modified parameters (delta = 0.0182 and a single cap 0.1899), and smaller bounds have been announced
since. The observation here is only that the parameters displayed in the paper need no change for k = 48.

The witness is written in the coordinates x = c*t with c = (26/25)/(A_1 + eps) = 2080/521, in which the support
is {|x|_1 <= 26/25, ...}. The file records c as `scale`; in the physical scale the ratio is
48 J(F~) / (c * I(F~)).

## File format

A witness file is a JSON object. `K`, `L`, `tau`, `g` and `caps_C` describe the support in the coordinates of the
file: T~ = {x >= 0, |x|_1 <= L, the coordinates exceeding g sum to at most C(their number)}, with
`caps_C[m-1]` = C(m) and the last value repeating; the marginal J integrates over the first K-1 coordinates with
sum at most `tau`. `terms` is a list of `[a, alpha, c]` with the integer c written as a string, for
F~ = sum c (L - |x|_1)^a m_alpha(x). If `scale` is present the coordinates are x = scale * t and the criterion is
K J/(scale * I) > 1; if it is absent the coordinates are x = 4t and the criterion is K J/I > 4.
`ratio_enclosure` is the Arb ball that was obtained when the file was written.

## Checking a witness

Requirements: Python 3 with `python-flint` (and `numpy` for the Monte Carlo test). See `requirements.txt`.

**Ball arithmetic** (about 50 minutes for the 212 witness, about 15 minutes for the k = 48 witness, on two processes of
a laptop):

    cd src
    WORKERS=2 python verify_ball.py ../witnesses/axiom212_k45_d21.json
    WORKERS=2 python verify_ball.py ../witnesses/stadlmann_k48_d22.json

This recomputes the Gram matrices of I and J on the support for B_d with Arb balls at 1536 bits and encloses the
ratio for the integer coefficient vector. Only the support, the scale and the coefficients are read from the file.

**Exact rational arithmetic** (about 12 hours each on two processes; a different program):

    cd src
    WITNESS=../witnesses/axiom212_k45_d21.json WORKERS=2 python run_profile.py axiom
    WITNESS=../witnesses/stadlmann_k48_d22.json WORKERS=2 python run_profile.py stadlmann48

This squares F first, integrates stratum by stratum in exact rationals, and writes the fraction to
`result_<profile>.txt`. It checkpoints after every stratum and resumes when restarted.

**How the witnesses were found.** For the 212 datum: `python gramtab.py 21` builds the tables and Gram matrices,
and `python gram_report.py gramtab_d21_gram.pkl` extracts the dominant generalized eigenvector by repeated
squaring, rounds it to integers and encloses its Rayleigh quotient with scalar ball operations. For a support in
`profiles.py`: `python support_run.py stadlmann48 25`, then
`python support_report.py stadlmann48 gramtab_stadlmann48_k48_d25_gram.pkl 21 25` does the same for every degree
from 21 to 25 and writes the witness of the lowest sufficient degree.

## What has been checked

- The two programs compute differently: one tabulates integrals of basis products in ball arithmetic and
  assembles Gram matrices, the other squares F first and integrates in exact rationals. They are not fully
  independent: they share the description of the strata and cells of the marginal (`axiom212_v2`), the
  multiplication of monomial symmetric polynomials (`capped_mono.mmul`) and small helpers, and the ball
  program's cell integrals are a port of the exact formulas. For the 212 witness they agree on all 40 digits
  compared, for the k = 48 witness on all 60.
- Independent recomputation by another person: Carlos Toledo (`carlostoledo1891`) reports that a program written
  from scratch, which shares no code with this repository and uses only the witness files, the coordinates above
  and the Lean definitions, gives for the d = 21 witness exactly the fraction in
  `witnesses/axiom212_k45_d21_exact_ratio.txt`, and for the d = 20 witness an exact value inside the ball given
  here ([comment of 16 November 2026 in AxiomMath/BGP212 issue 2](https://github.com/AxiomMath/BGP212/issues/2#issuecomment-6918984119),
  [his code](https://github.com/carlostoledo1891/cert-machine/tree/main/instruments/bgp212-witness)). I have not
  run his program myself.
- 282 random table entries of the ball program agree with the exact program to 1e-60. At degree 4 the Gram
  matrices agree entry by entry with an earlier implementation in the power-sum basis (not included here).
- A run of the ball program that avoids FLINT matrix products (`NOMAT=1`) gives the same balls for every Gram
  entry of the 212 datum at d = 21. See the note on FLINT below.
- As a test of the reading of the definitions, `src/bruteforce_def.py` compares the ball program with a Monte
  Carlo evaluation that uses no cell decomposition, only the pointwise membership test for T and the condition of
  `Jregion` in the Lean file (4 variables, 10 basis functions, 200 million samples): all 110 Gram entries of I
  and J agree within two standard errors (relative errors 0.04 % to 0.4 %). This tests the set-up, not the
  values in 45 or 48 variables.
- On another support both programs reproduce a published number: the 34 digits printed in Kintali's paper,
  1.000774701187468444763475111497864, for the 1295-term witness of `AxiomMath/PrimeGapsLib`
  (`PrimeGapsCert/Gap246/k50e25d25n1295.json`, read at commit `1faa7b1`, not included here) on his k = 48
  support. The exact program gives 48 J/I = 1.00077470118746844476347511149786443813176... as the fraction in
  `witnesses/kintali_k48_exact_ratio.txt` (numerator and denominator of 482 digits each; 46 hours on two
  processes, `WITNESS=<that file> WORKERS=2 python run_profile.py kintali`).
- The agreement with Theorem 11.1 of the Axiom draft to 17 digits is evidence from outside these programs that
  the support and the marginal of the 212 datum are the ones the authors used.
- H(48) = 236 was confirmed by exhaustive search (no admissible 48-tuple of diameter at most 235; an explicit
  one of diameter 236). That search is not part of this repository.

## What this is not

- These are computer-algebra computations, not formal proofs. In particular the 212 witness is not a Lean proof
  of `Gap212Certificate`, and it does not replace the authors' own F or the matrices behind Theorem 11.1.
- Nothing here examines the analytic parts of the papers: the equidistribution estimates, the Harman reduction,
  bilinear Bombieri-Vinogradov, Sections 5-8 of the Axiom draft, or the proofs in arXiv:2608.31126. The k = 48
  observation is conditional on Stadlmann's theorem exactly as her own k = 49 result is.

## A note on FLINT 3.6.0

`python-flint` 0.9.0 ships FLINT 3.6.0, which has a known bug in an FFT-based multiplication routine (fixed
upstream in FLINT pull request 2790). The exact program uses scalar `fmpq` arithmetic only, no polynomial or
matrix products. The ball program uses `arb_mat` products while building tables unless `NOMAT=1` is set; the
final enclosure of the ratio always uses scalar ball operations. Both modes gave identical Gram matrices for the
212 datum.

## Files

- `src/gramtab.py`, `src/gram_report.py`, `src/verify_ball.py`, `src/support_run.py`, `src/support_report.py`:
  the ball program.
- `src/run_profile.py`, `src/kintali_J.py`, `src/kintali_capped.py`, `src/capped_mono.py`,
  `src/mono_engine.py`: the exact program. `src/axiom212_v2.py` holds the cell geometry of the marginal, used by
  both. `src/profiles.py` lists the supports.
- `src/bruteforce_def.py`: the Monte Carlo test against the definitions.
- `kapila_h2_h5/`: tuples, checks and notes for the H_2 to H_5 bounds of hexagon:2610.00419v2 (own README).

## Authorship and licence

Fredrik Prüzelius (ORCID [0009-0006-9202-9610](https://orcid.org/0009-0006-9202-9610)). The programs and this
text were written with Claude (Anthropic's AI model); the runs were made on my own machine.

Code: MIT licence (see `LICENSE`). Witness files and this text: CC BY 4.0.

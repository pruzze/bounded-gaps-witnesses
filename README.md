# An independent witness for the numerical certificate of H_1 <= 212

Axiom Math's Lean formalization of the bound H_1 <= 212 for gaps between consecutive primes
([AxiomMath/BGP212](https://github.com/AxiomMath/BGP212), read at commit `6641780`) is conditional on a
numerical certificate, `Gap212.Gap212Certificate`: the existence of a symmetric, square-integrable function F
supported on the set T_45 of the parameter datum `gap212Params` with

    0 < I(F) < 45 * J(F).

At the time of writing, no such F is published in that repository or in the draft paper
("A New Bound for Small Gaps Between Primes", Axiom Math, preliminary draft of 3 September 2026), whose
Theorem 11.1 states the value of 45 J / I for a vector P* that is not printed.

This repository contains an F that was computed independently, the exact value of 45 J(F) / I(F) for it, and the
programs that produced and checked it. It is not the authors' F.

## The witness

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

The coefficients are integers of at most 51 digits. `witnesses/axiom212_k45_d21_exact_ratio.txt` holds the value
for d = 21 as a fraction (numerator 387 digits, denominator 386 digits).

The draft states 4.00438409833460131937... for its P*. The value above agrees with it to 17 significant digits and
is slightly larger, as expected if the vector here is the maximiser over B_21 and P* a rounded maximiser.

### File format

A witness file is a JSON object. `K`, `L`, `tau`, `g` and `caps_C` describe the support in dilated coordinates
(`caps_C[m-1]` = C(m); the last value repeats). `terms` is a list of `[a, alpha, c]` with the integer c written
as a string. `ratio_enclosure` is the Arb ball for 45 J/I that was obtained when the file was written.

## Checking it

Requirements: Python 3 with `python-flint` (and `numpy` for the Monte Carlo test). See `requirements.txt`.

**Ball arithmetic** (about 50 minutes on two processes of a laptop):

    cd src
    WORKERS=2 python verify_ball.py ../witnesses/axiom212_k45_d21.json

This recomputes the Gram matrices of I and J on the support for B_21 with Arb balls at 1536 bits and encloses
45 J/I for the integer coefficient vector. Only the support and the coefficients are read from the file.

**Exact rational arithmetic** (about 12 hours on two processes; a different program):

    cd src
    WITNESS=../witnesses/axiom212_k45_d21.json WORKERS=2 python run_profile.py axiom

This squares F first, integrates stratum by stratum in exact rationals, and writes the fraction to
`result_axiom.txt`. It checkpoints after every stratum and resumes when restarted.

**How the witness was found:** `python gramtab.py 21` builds the tables and Gram matrices;
`python gram_report.py gramtab_d21_gram.pkl` extracts the dominant generalized eigenvector by repeated squaring,
rounds it to integers and encloses its Rayleigh quotient with scalar ball operations.

## What has been checked

- The two programs compute differently: one tabulates integrals of basis products in ball arithmetic and
  assembles Gram matrices, the other squares F first and integrates in exact rationals. They are not fully
  independent: they share the description of the strata and cells of the marginal (`axiom212_v2`), the
  multiplication of monomial symmetric polynomials (`capped_mono.mmul`) and small helpers, and the ball
  program's cell integrals are a port of the exact formulas. For the d = 21 witness they agree on all 40 digits
  compared.
- 282 random table entries of the ball program agree with the exact program to 1e-60. At degree 4 the Gram
  matrices agree entry by entry with an earlier implementation in the power-sum basis (not included here).
- A run of the ball program that avoids FLINT matrix products (`NOMAT=1`) gives the same balls for every Gram
  entry at d = 21. See the note on FLINT below.
- The definitions were read against `Gap212Challenge/Basic.lean` (`T`, `Iint`, `Jint`, `gap212Params`,
  `Certificate`). As a test of that reading, `src/bruteforce_def.py` compares the ball program with a Monte Carlo
  evaluation that uses no cell decomposition, only the pointwise membership test for T and the condition of
  `Jregion` (4 variables, 10 basis functions, 200 million samples): all 110 Gram entries of I and J agree within
  two standard errors (relative errors 0.04 % to 0.4 %). This tests the set-up, not the 45-variable value.
- On another support the ball program reproduces a published number: the 34 digits printed in Kintali's paper
  (1 September 2026) for the 1295-term witness of `AxiomMath/PrimeGapsLib` on his k = 48 support.
- The agreement with Theorem 11.1 of the draft to 17 digits (above) is evidence from outside these programs that
  the support and the marginal are the ones the authors used.

## What this is not

- It is a computer-algebra computation, not a Lean proof of `Gap212Certificate`.
- It does not replace the authors' own F or the matrices behind Theorem 11.1.
- It says nothing about the analytic hypotheses of the formalization (the equidistribution estimates, the Harman
  reduction, bilinear Bombieri-Vinogradov) or about Sections 5-8 of the draft.

## A note on FLINT 3.6.0

`python-flint` 0.9.0 ships FLINT 3.6.0, which has a known bug in an FFT-based multiplication routine (fixed
upstream in FLINT pull request 2790). The exact program uses scalar `fmpq` arithmetic only, no polynomial or
matrix products. The ball program uses `arb_mat` products while building tables unless `NOMAT=1` is set; the
final enclosure of the ratio always uses scalar ball operations. Both modes gave identical Gram matrices here.

## Files

- `src/gramtab.py`, `src/gram_report.py`, `src/verify_ball.py`: the ball program.
- `src/run_profile.py`, `src/kintali_J.py`, `src/kintali_capped.py`, `src/capped_mono.py`,
  `src/mono_engine.py`: the exact program. `src/axiom212_v2.py` holds the cell geometry of the marginal, used by
  both. `src/profiles.py` lists the supports.
- `src/bruteforce_def.py`: the Monte Carlo test against the definitions.

## Authorship and licence

Fredrik Prüzelius (ORCID [0009-0006-9202-9610](https://orcid.org/0009-0006-9202-9610)). The programs and this
text were written with Claude (Anthropic's AI model); the runs were made on my own machine.

Code: MIT licence (see `LICENSE`). Witness files and this text: CC BY 4.0.

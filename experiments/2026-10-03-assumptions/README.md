# VERDICT: the wrong premise is `41, 47, 53 | d` (the D0 convention). Both optima that exist in the literature or in this repo as *proved/published* facts lie outside the all-in-d family that our engine searches: the proved-minimal n=35 AP pays 29^2 on a single hit term, and Pfoertner's published 37-term sums-of-two-squares record pays 31^2 on a single hit term. Two for two. The square option is not a 5x garnish; it is what real optima look like. Assumptions 1 and (c) are closed negatively, with proofs.

H4, laptop only, 2026-10-03. Code here; every numeric claim below was produced by one of these
scripts, and every Loeschian test routes through `src/loeschian.py` (the C bitmap in `joint.c` was
cross-checked against `is_loeschian` and against A003136 before use).

| file | what it does |
|---|---|
| `audit.py` | structure audit of every progression in `records.json` |
| `joint.c` | exact Loeschian bitmap + **exact** counts of n-term APs with last term <= X |
| `count.py` | per-prime local-density count model, calibrated against those exact counts |
| `bottom.py` | the naive independence model (kept because it **fails**, see §3) |
| `joint.py` | superseded by `joint.c` (no numpy on this laptop); ignore it |

Reproduce: `python3 audit.py`; `cc -O2 -o joint joint.c && ./joint 2000000 12 18 && ./joint 19000000 22 28`; `python3 count.py`.

---

## 1. DEMONSTRATED: the two independently-established optima both use the square option

`audit.py` output, all terms re-verified Loeschian:

| progression | how established | bad q in (n/2, n] left out of d | where the q^2 lands |
|---|---|---|---|
| n=35, a=219830911, d=2709630 | **exhaustive sweep, proved minimal** (`records.json`) | **29** | 29^2 \|\| t_15 |
| n=36 / 43 / 47 / 55 (ours) | engine search | none | — |

d for the proved-minimal 35 is `2 * 3^2 * 5 * 7 * 11 * 17 * 23`: the forced primes 2,5,11,17 plus
23, but **not** 29. The one progression in this repo whose optimality is proved does not live in
the all-in-d family. Every engine-found record does, because the engine can only emit that family.

Independent second data point, from the only published long AP in an analogous set
(bad primes = 3 mod 4), Pfoertner's 37-term sums-of-two-squares AP a=2215647809, d=25438644:

```
d = 2^2 * 3^3 * 7^2 * 11 * 19 * 23        bad (3 mod 4) primes <= 37: 3,7,11,19,23,31
all 37 terms are sums of two squares      31 is NOT in d: 31^2 || t_12
```

He forced the primes up to 23 into d and took the square option on exactly one prime, 31. Same
shape as our proved-minimal 35. For n=58 the corresponding choice is: 2,5,11,17,23,29 | d (forced),
and 41/47/53 each either in d or paying q^2 on a single hit index.

**Quantified.** The local density for a bad q with n/2 < q <= n, over (a, d mod q^inf), is

```
in d:      (q-1)/q
square:    (1/q)(q-1)/q + ((q-1)/q)(s/q)(1/(q+1)),  s = 2q-n single-hit classes
           (classes with two hits are impossible: the two hit terms differ by q*d)
```
using `P(v_q even and > 0 | q | t) = 1/(q+1)`. Ratio square-family / all-in-d per prime:
**1.585 (q=41, s=24), 1.766 (q=47, s=36), 1.906 (q=53, s=48) — product 5.33x.** This reproduces
the PROGRESS.md figure (19% / 81%, ~5x) from an independent derivation, so that number is now
confirmed twice. The engine searches the 19%.

## 2. PROVED: assumption 1 has no escape, including the non-primitive one offered in the register

Let pi(n) = {bad q : v_q(n) odd}; n is Loeschian iff pi(n) = {}. pi(mn) = pi(m) XOR pi(n).

*Claim.* If a, a+d, ..., a+(n-1)d are all Loeschian with g = gcd(a,d) and a = g a', d = g d',
then g is Loeschian and a', a'+d', ... are all Loeschian. So every non-primitive solution is
(Loeschian constant) x (primitive solution) and nothing new can live there.

*Proof.* All terms Loeschian means pi(t'_k) = pi(g) for every k. Suppose q in pi(g). Then
q | t'_k for all k, so q | d' and q | a', contradicting gcd(a',d')=1. Hence pi(g) = {}, g is
Loeschian, and pi(t'_k) = {} for all k. QED

The same XOR argument kills the "variable multiplier" variant: a progression inside a scaled set
c*L is c times a progression in L. Register item 1 is closed; do not reopen it.

## 3. MEASURED: the independence model overpredicts long runs, by a factor that compounds

`joint.c` gives **exact** counts in the forced family
(d = base(n)*m, base(n) = 3*prod(bad q <= n/2), a = 1 mod 3, a != 0 mod the forced bad q,
last term <= X):

| n | X | admissible pairs | per-term rate rho1 | EXACT #APs | rho_eff = (count/pairs)^(1/n) | rho_eff/rho1 |
|---|---|---|---|---|---|---|
| 12 | 2e6 | 8.080e8 | 0.6153 | 434143 | 0.534 | 0.868 |
| 14 | 2e6 | 6.836e8 | 0.6153 | 95725 | 0.531 | 0.853 |
| 16 | 2e6 | 5.925e8 | 0.6153 | 16975 | 0.520 | 0.836 |
| 18 | 2e6 | 5.227e8 | 0.6153 | 2851 | 0.510 | 0.820 |
| 22 | 1.9e7 | 3.156e9 | 0.6218 | 1819 | 0.521 | 0.837 |
| 24 | 1.9e7 | 2.881e9 | 0.6218 | 285 | 0.511 | 0.821 |
| 26 | 1.9e7 | 2.651e9 | 0.6218 | 48 | 0.504 | 0.810 |
| 28 | 1.9e7 | 2.454e9 | 0.6218 | 5 | 0.489 | 0.787 |

A directly measured per-term rate is **not** the right exponent base: using it overpredicts the
count by (1/0.82)^n, which is 50x at n=22 and ~1e5 at n=58. `bottom.py` shows the damage: fed the
measured rho1 it predicts 280 APs of length 35 with last term <= 3.1e8 where exactly **one**
exists (the proved minimum) -- a 280x overprediction, consistent with (1/0.82)^35 = 900 to within
the sampling error. The correction is negative correlation between terms at each prime
(for one prime q, true `1 - n/q + n/q^2` vs independence `(1-1/q)^n = 1 - n/q + C(n,2)/q^2`).

*Caveat, labelled:* PROGRESS.md's quoted `observed/expected = 0.82 +- 0.09` is probably the same
0.82 appearing as a single factor because its rho was fitted to observed run statistics (which
absorbs the correlation) rather than measured per-term. If so the 0.4/hour rate stands. **If any
planning number was ever computed as pool x (measured per-term rate)^58, it is too high by ~1e5**
and should be recomputed. That is worth 10 minutes of someone's time to check.

## 4. MEASURED: there is no unsearched "bottom". Term size 1e16-1e19 is the right place.

Only 3 | d and the bad q <= n/2 are truly forced, so base(58) = 3*2*5*11*17*23*29 = **3741870**,
not D0 = 382160924970 -- a factor 102131 smaller, so at a given X the square-option family holds
~1e5 times more (a,d) pairs. That looked like a huge unsearched region below the engine's floor.
It is not, and here is why: each of those pairs is 5.24e-5 as likely to work
(the three per-prime ratios of §1 inverted), so the net gain is the same 5.33x, and **both
families' counts scale as X^2**, so the X at which the first 58 appears differs by only sqrt(5.33)
= 2.3x between them.

Plugging the calibrated rho_eff (rho1 from PROGRESS.md's Q=53 Monte Carlo, times the 0.78-0.82
correlation factor measured in §3) into pairs(58,X) = 2.47e-10 X^2:

| X | E[#58-term APs with all terms <= X] |
|---|---|
| 1e12 | ~1e-6 |
| 1e14 | ~4e-4 |
| 1e16 | ~0.1 |
| 1e18 | ~7 |
| 1e20 | ~1e4 |

Sanity check at a length where we have data: the same formula gives ~33 expected 55-term APs with
terms <= 2e16, and two were found in a partial sweep of that region. Consistent.

**Conclusion: a 58-term AP with all terms below ~1e16 very likely does not exist** (E ~ 0.1), the
first ones appear around 1e17-1e18, and the engine at 1e16-1e20 is in the right band. H1's premise
(that a large win hides below the term-size floor) is therefore worth much less than it looks:
there are no solutions down there to find. The gain from going smaller is bounded by the ~6x/decade
yield slope only over the 1-2 decades between 1e17 and the engine's actual working size.

## 5. CLOSED: (c) no map from the sums-of-two-squares problem

Pfoertner's 37 terms: **0 of 37 are Loeschian**, and not by luck -- 3 | d and a = 2 mod 3, so all
37 terms are = 2 mod 3, and no Loeschian number is = 2 mod 3. Their bad-prime (2 mod 3) odd-
valuation vectors take **33 distinct values** over the 37 terms, so by the XOR lemma of §2 no
multiplier, and no scaling of any kind, can transport that progression into L. The two problems'
prime structures are independent (a prime splits in Z[i] and in Z[omega] according to independent
conditions mod 4 and mod 3); there is no ring map and no congruence map. What transfers is only
the *shape* of his answer, which is §1 -- and that is the valuable part.

## 6. Speculation (clearly labelled, not demonstrated)

- The structural wall behind "no search-free construction": fixing the factorisation shape of m
  terms of an AP imposes m-2 quadratic conditions (any three terms satisfy one linear relation, so
  three forced shapes give a conic, four give a genus-1 curve, five give a surface with finitely
  many rational points). That matches, and explains, the measured outcome of `exp_auto` (3 terms ->
  1-parameter families) and `exp_ec` (4 terms -> rank 0/1 curves, 72-digit points). It predicts
  that **no construction can control more than ~4 of the 58 terms**, so the 57-term record was not
  a construction: it was a small clever search in the right family. Untested but checkable.
- Under that reading, "exploit the prime structure" + "two of {large d, square option, automatic
  terms}" most plausibly means: forced primes in d, square option on 41/47/53, and a hand-sized
  search at term size ~1e17. §1 says the square option is the half of the space we are not
  searching; that is where I would spend the next engine-week.

## 7. What I attacked and could NOT break

- Register 1 (scaling / primitivity): proved closed, §2.
- Register 2 (linear function identically a norm): unchanged; the dimension count in §6 explains
  why the recurrence/conic/elliptic escapes cap out where `exp_ec` measured them.
- (c) cross-problem map: proved closed, §5.
- (a) long APs in multiplicative semigroups: nothing found beyond the above; the only usable
  consequence of multiplicativity is the XOR lemma, and it is an obstruction, not a tool.

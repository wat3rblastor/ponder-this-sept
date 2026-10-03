# VERDICT: LOSS — square mode is CORRECT and VALIDATES (it finds the proved-optimal 35-term AP in 0.5 s, which the all-in-d family provably cannot contain), but its long-run yield per unit of work at n=58 is BELOW 1, not above 3: see the table below. And the term-size claim does not survive either: the 1e5 drop in d is cancelled by a 2.5e5 rise in the a-floor, so the smallest reachable last term is 2.05e12 (all-in-d) vs 3.12e12 (all-square) — slightly WORSE, not 1e5 better. The "5.33x" is real but it is the ratio of the UNION of all 8 modes to all-in-d; the pure all-square family is 0.406x of all-in-d, and the per-prime factors 1.585/1.766/1.906 decompose as 1 (in-d) + (2q-n)/q (square). Searching 19% is right; the other 81% is spread over 7 modes and costs proportional work to reach.

Date: 2026-10-03. Laptop (8 cores, 8 GB), CPU engine only.

--------------------------------------------------------------------------

## 1. What was implemented

`apsearch_sq.c` = a copy of `src/c/apsearch.c` (untouched) plus per-unit MODES.
A unit is now `(K, shift, mode)`:

* `--ind q,..`     bracket prime IN d:  `d = K * D0 * prod(ind)`
* `--square q,..`  bracket prime in SQUARE mode: `q !| d`, and stage 1 gains a
  component `(modulus q^2, start -(max(0,n-q))*d, step -d,
  count min(n-1,q-1)-max(0,n-q)+1)` — i.e. `2q-n` classes when `q < n`.

The mode is part of the unit identity: `"ind"` and `"sq"` appear in every
`--out` record (unit lines and hit lines) and on every `*** n=` line.
Square-mode primes are explicitly excluded from tier B **and from tier C**
(tier C's mask says "q misses the whole window", the exact complement of the
square family; applying it would reject every candidate the mode generates).
The original CRT self-check now also covers the `q^2` components, and a new
per-unit semantic self-check confirms that every residue of a `q^2` component
makes `q^2` divide exactly one in-window term and `q` divide no other.

See `BUILD.md` for flags and exact command lines. Raw logs in `raw/`.

## 2. Validation — the decisive n=35 test: PASS

```
./build/apsearch_sq --nterms 35 --D0 5610 --ind 23 --square 29 \
    --kmin 21 --kmax 21 --shifts 1 --report 33 --histmin 33 --modcap 200000000000
...
*** n=35 a=219830911 d=2709630 (K=21 shift=0 ind=[23] sq=[29])
K=21 shift=0 ... R=2861568 surv=1972 conf=1946 ge35=1 minlast=311958331 best=35 0.5s
```

It finds the proved minimal-last-term 35-term Loeschian AP
(a=219830911, d=2709630, last term 311958331) in **0.5 s**, after 2.86e6
residues. Verified independently:

* `python3 src/verify.py 219830911 2709630 35` → OVERALL: PASS
* `python3 src/crosscheck.py 219830911 2709630 35` → PASS (35 terms, each with
  an explicit x^2+xy+y^2 representation)

This AP is *not in the all-in-d family at all* (29 !| d, 29^2 || t_15), so no
amount of CPU time in that family can reach it. Control run in the all-in-d
family (`--D0 162690`, i.e. 29 | d, all 56 d-values with last term below the
optimum, full a-range to 3.9e8): best run **31**, `raw/n35_alld_control.*`.

So the square path is implemented correctly and does exactly what it is for.

## 3. Regressions

* n=47 record admissible in all-in-d mode:
  `--nterms 47 --D0 3741870 --ind 41,47,53 --kmin 205 --atest 2646171143023357`
  → every stage-1 component OK, no stage-2 rejection, run length 47.
  (`--atest` is the cheap form of the regression: a full unit at that modulus
  is 4.7e9 residues ≈ 14 min.)
* all-in-d path bit-identical to production: `build/apsearch` and
  `build/apsearch_sq --D0 3741870 --ind 41,47,53` with the same args
  (n=58, K=7, shift 0, modcap 1e13) both report
  `MOD=93760753290 R=1732900 surv=19 conf=19 best=18 a=22225747809637
  d=2675126474790`.
* engine `--selftest` still OK.

## 4. Where the 5.33x went (the most important result here)

`density.py` recomputes, from scratch, the count of 58-term candidates per mode
at a fixed last-term bound X. For each bad prime r <= 2000 the local density of
admissible a is

* `r | d`            : `1 - 1/r`        (r must not divide a)
* `r` square mode    : `(2r - n)/r^2`   (r^2 hits exactly one term)
* `r !| d`, `r > n`  : `1 - n/r`        (r misses the window)

and the number of available d below D is `D / dmin(mode)`, with
`dmin = 3741870 * prod(in-d bracket primes)`. Then
`N(X) = (X/((n-1)*dmin)) * X * rho`.

```
mode                       dmin          rho      1/rho   (n-1)*dmin    N(X)/N_all-in-d
all-in-d           382160924970    5.172e-12  1.933e+11 21783172723290     1.0000
sq=41                9320998170    7.569e-14  1.321e+13   531296895690     0.6000
sq=47                8131083510    8.612e-14  1.161e+13   463471760070     0.7826
sq=53                7210583490    9.008e-14  1.110e+13   411003258930     0.9231
sq=41,47              198319110    1.260e-15  7.934e+14    11304189270     0.4696
sq=41,53              175867890    1.318e-15  7.586e+14    10024469730     0.5538
sq=47,53              153416670    1.500e-15  6.667e+14     8744750190     0.7224
sq=41,47,53             3741870    2.195e-17  4.556e+16      213286590     0.4334
                                              sum over all 8 modes = 5.485
```

The per-prime factors quoted in the brief, 1.585 / 1.766 / 1.906, are exactly
`1 + (2q-n)/q` for q = 41/47/53, and they factor as

    1          <- the in-d option (what the engine already searches)
  + (2q-n)/q   <- the square option  (0.585, 0.766, 0.906)

So the 5.33x is the ratio of the **union of all 8 modes** to the all-in-d mode
alone. It is NOT "the square families hold 5.33x more than the all-in-d
family": the pure all-square family holds `(24/41)(36/47)(48/53) = 0.406x`,
i.e. **less than half** of what the engine already searches. The statement
"the engine currently searches the 19%" is literally true (1/5.33 = 0.1875) but
the remaining 81% is split over seven other modes, each of which costs work in
proportion to its size — it is extra *space*, not a shortcut.

Term size, same table: the modulus budget spent on `q^2` raises the smallest a
that is admissible at all (column `1/rho`, the expected first admissible a) by
almost exactly the factor by which it lowers d. Solving `N(X) = 1`:

```
smallest last term with N(X) = 1
  all-in-d          X = 2.05e12        sq=41,47       X = 3.00e12
  sq=41             X = 2.65e12        sq=41,53       X = 2.76e12
  sq=47             X = 2.32e12        sq=47,53       X = 2.42e12
  sq=53             X = 2.14e12        sq=41,47,53    X = 3.12e12
```

Every mode bottoms out at the same 2-3e12. The brief's "dropping 41,47,53 from
d drops the step base from 3.82e11 to 3.74e6, i.e. 1e5 in reachable term size"
counts the d side of the ledger only; the a side moves 2.5e5 the other way
(1/rho goes 1.93e11 -> 4.56e16). Net, all-square reaches term sizes 1.5x
*larger*, not 1e5 smaller.

## 5. Measurement at n = 58, matched wall clock (900 s each, same laptop)

Same binary, same `--nterms 58 --D0 3741870 --kmin 1 --kmax 999 --shifts 8
--modcap 1.4e15 --mintierb 1 --b2 2000 --report 30 --histmin 30 --tlimit 900`,
run concurrently on separate cores; only the mode differs. Raw: `raw/n58_*`.

| | all-in-d (production family) | sq=53 (one square) | sq=41,47,53 (all square) |
|---|---|---|---|
| d | K * 3.8216e11 | K * 7.2106e9 | K * 3.74187e6 |
| MOD | 1.1337e15 | 2.6337e14 | 1.3108e15 |
| tier-B primes pinned | 59,71,83,89,101,107,113 | 59,71,83,89,101 | 59,71 |
| residues examined | 5.006e9 | 9.898e9 | 1.329e10 |
| residues/s | 5.5e6 | 1.10e7 | 1.48e7 |
| **survivors (stage-3 tests)** | **133 407** | **60 030** | **4 554** |
| survivors/s | 148 | 66.7 | 5.06 |
| runs >= 20 | 390 | 156 | 5 |
| runs >= 25 | 46 | 19 | 1 |
| runs >= 30 | 5 | 2 | 0 |
| runs >= 35 | 0 | 2 (n=37, same AP twice) | 0 |
| longest run | 34 | **37** | 29 |
| smallest last term among runs >= 30 | 2.182e14 | 2.131e14 | — |

**Long-run yield per unit of work, square vs all-in-d:**

| threshold | sq=53 / all-in-d | sq=41,47,53 / all-in-d |
|---|---|---|
| >= 20 | **0.40** | 0.013 |
| >= 25 | **0.41** | 0.022 |
| >= 30 | **0.40** | < 0.2 (0 of 5) |

The ratio is flat at 0.40 across thresholds, which is a sign it is real and not
small-number noise at the tail. Required for a WIN was >= 3.0.

Per *survivor* the three modes are equivalent (2.92e-3, 2.60e-3, 1.10e-3 runs
>= 20 per survivor) — exactly as the structure predicts, since a survivor in
any mode has been cleared against every bad prime <= 2000. The entire
difference is **survivors per second**, i.e. enumeration efficiency:

* A square component costs `q^2` of modulus (1681 / 2209 / 2809) and returns
  nothing in enumeration efficiency — it selects a *different, smaller* family.
* A pinned tier-B prime costs `q` of modulus and returns `q/(q-n)` in
  enumeration efficiency (3.32x at 83, 2.87x at 89, 2.35x at 101, 2.18x at
  107, 2.05x at 113).

The modulus budget is hard-capped (MOD must stay below ~3e16 for the 51-bit
Barrett in `loesch_core.h`, and below 2^64/64 regardless). Spending 1.04e10 of
it on `41^2*47^2*53^2` leaves room for only 59 and 71, so 83,89,101,107,113
fall back to tier C, where their window conditions are paid by *constructing
and rejecting* candidates: 4.5x fewer survivors per second for sq=53 (measured
2.2x, the rest recovered because that mode's MOD is 4.3x smaller and its terms
correspondingly easier), and ~100x fewer for all-square.

Note the comparison is if anything generous to square mode: sq=53 searches
terms 4.3x smaller than all-in-d at the same shift, which raises its per-term
pass probability. Normalised for term size it would look worse than 0.40.

A 37-term AP was found in the sq=53 family and verified independently:
`a=212569971883993, d=14421166980, n=37` —
`src/verify.py` PASS, `src/crosscheck.py` PASS (37 explicit x^2+xy+y^2 forms).
53 does not divide d there; 53^2 divides exactly one term. So the new path
produces real, verifiable progressions that the production engine's family does
not contain.

## 6. What to do with this

* **Do not switch the engine to square mode.** At equal wall clock it finds
  0.40x (one square prime) to 0.013x (three) as many long runs.
* **Keep the code.** It is correct, it is 60 lines, it costs nothing when
  `--square` is empty (bit-identical to production), and it is the only way to
  reach 81% of the solution space. The n=35 case is the proof that this matters
  in principle: the proved-optimal 35-term AP is in a square family, so if the
  58-term answer happens to be too, the production engine will never see it no
  matter how long it runs.
* **If the space ever has to be covered, the order is by yield:** all-in-d
  (1.00), then sq=53 (0.40), then sq=47 (expect ~0.3), sq=41, then the
  two-square modes, and all-square last (0.013). Covering all 8 modes costs
  ~1/0.19 = 5.3x the time of covering all-in-d alone and finds 5.3x the
  progressions, which is the cost/gain cancellation already measured in the
  pin-set experiment, rediscovered here in a different guise.
* **Correct the two prior results.** (1) "Square families hold 5.33x more" is a
  misreading of `prod (1 + (2q-n)/q)`: that is the union of all 8 modes over
  all-in-d; the square families themselves hold 0.41x. (2) "Square mode reaches
  term sizes 1e5 smaller" double-counts: the q^2 congruences raise the
  smallest admissible a by 2.5e5 while d falls by 1e5.

Exact command line for the best square configuration, should it be wanted:

```sh
./build/apsearch_sq --nterms 58 --D0 3741870 --ind 41,47 --square 53 \
    --kmin 1 --kmax 999 --shifts 8 --modcap 1400000000000000 --mintierb 1 \
    --report 50 --out out.jsonl
```

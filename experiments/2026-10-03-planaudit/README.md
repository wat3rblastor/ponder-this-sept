# LOSS: the corrected ranking delivers 1.00x the measured yield of the current one (top 10% of residues; 1.08x at 20%, 1.00x at 50%) -- the plan ordering is not where the problem is.

Audit of the unit score in `tools/plan_units.py`, 2026-10-03, laptop only.
Nothing in `tools/` or `src/c/` was modified; everything here is new.

Measured numbers up front:

| quantity | value |
|---|---|
| Spearman(measured log yield, **old** score), n=100 units | **+0.881** |
| Spearman(measured log yield, **new** score), n=100 units | **+0.892** |
| OLS slope of measured on old score (1.0 = correct scale) | **1.007 +- 0.049** |
| residual scatter sd(M - score) | **0.274-0.283** vs **0.300** measurement noise |
| yield ratio new/old ranking, best 10% of the residue budget | **1.000** |
| same, best 20% / 30% / 50% | 1.108 / 1.084 / 1.003 |
| perfect-oracle ranking vs old, best 10% | 1.080 (and that number still contains noise-fitting) |
| Spearman(old score, new score) over 3000 plan units | 0.99928 |
| plan-wide absolute E58: `plan_units.py` **underestimates** by | **2.28x** (+- 0.06 sem) |

The score's *relative* ordering is right to within measurement error: the
residual scatter around it is statistically indistinguishable from the noise of
the measurement itself, so the unexplained per-unit modelling error is bounded
at well under 1.3x (1 sigma). None of the three errors in the project's
analytic formula is present in the planner's score.

## 1. What the score actually computes

For a unit `(K, s)`, `v = ly + lgood + st` is `log` of the expected number of
58-term Loeschian APs **per stage-1 residue**, up to a unit-independent
constant; `E58 = CONST * exp(v) * res`, `res = 4 * prod_{q pinned} (q-58)`.
Structurally:

* baseline `rho(t)^58` with `rho(t) = 0.71992*(ln t/39.144)^-0.485` stands in
  for the product over *all* bad primes of the generic AP-local factor. This is
  legitimate: the generic factor for a bad prime `q > n` in an `n`-term AP is
  `(q-n)/q + n/q^2 = 1 - n/q + O(1/q^2)`, and the single-term Loeschian factor
  is `(q/(q+1))^n = 1 - n/q + O(1/q^2)`. They agree to first order, so
  `rho^58` is a valid baseline and the remaining per-prime work is only to
  *correct* the primes whose situation differs from generic.
* `ly` applies exactly those corrections: `q/(q-58)` for pinned (tier-B) q
  (the odometer never materialises a term divisible by q), and
  `((q-1)/q)/((q-58)/q)` for q | d (not filtered anywhere in the engine; the
  whole window dies iff q | a).
* `st` is the size term, `lgood` the correction for good primes p | K.

**The three errors of the quoted analytic model are all absent:**

(a) `corr_q = 1 - n/(q+1)` applied at q < n. Cannot happen here: `BAD` is
`[q prime : q = 2 mod 3, q > 58]`, and *every* bad prime below 58
(2,5,11,17,23,29,41,47,53) divides `D0`, so the planner never needs a local
factor for `q < n`, `n/2 < q < n`, or `q = 59..` other than the three cases it
actually implements (pinned / `q | d` / generic). Checked: `D0 = 3 * 2*5*11*17*23*29*41*47*53`
is exactly 3 times the product of all bad primes below 58.

(b) the missing `prod_{p | d} (1 - 1/p)`. Present for the primes that matter
for ranking (bad `q | K`: `plan_units.py:52-54`; good `p | K`:
`plan_units.py:62-66`). Missing only for the tier-A primes `p | D0`, which are
the same for every unit -- a constant, so it cannot reorder anything.

(c) the triangular 1/2 in the pair count: belongs to the "expected number of
58-APs in a range" formula, which this score does not use.

**No double counting of the pinned primes**: `ly` removes `log((q-58)/q)` for
each pinned q exactly once, and the `(q-58)` factors appear once, in `res`.

## 2. Discrepancies found (with measured size at n = 58)

1. **`wbar` is a CUDA-only quantity applied to a score that also drives the CPU
   engine.** `apsearch_cuda.cu:419-422` walks *unreduced* class representatives
   `A = Rred + w*MOD`, so a CUDA unit covers window multiples
   `[64s+w, 64s+64+w)` with `w < c1+c2`; `plan_units.py:55-61` models the mean
   of that, `wbar`, which reaches 57 in my sample. `apsearch.c:352` keeps
   `R < MOD`, so for the CPU engine `w = 0` and the offset is only
   `E[R]/MOD = 0.5`. Effect of setting `wbar = 0.5`: sd 0.059 in log over 3000
   plan units, sd 0.127 over shift-0 units, **worst single unit 0.87 log =
   2.39x** (small K, shift 0). Ranking effect: Spearman(old,new) = 0.99928.
2. **Term indexing off by one**: `size_term` evaluates `rho` at
   `a + k*d + d` (`plan_units.py:78`), i.e. the terms `a+d .. a+58d` instead of
   `a .. a+57d`.
3. **Coarse grids**: `ln ln t` is averaged over 6 values of k and 8 of the 64
   window positions.
4. **Hard-coded `ln T = 40.5`** in the good-prime correction
   (`plan_units.py:66`) instead of the unit's own term size; `GOOD` truncated at
   1000 and `BAD` at 10000, so a bad/good prime factor of K above that bound
   gets no correction (K runs to 600000 in the deployed plan).
   Items 2-4 together are worth sd **0.011** in log (1.1%); they are swamped by
   item 1.
5. **`MODCAP` default mismatch** (planner 2e15, `apsearch.c` 6e15,
   `apsearch_cuda.cu` 2e15, `tools/run_plan.sh` 2e16, `GOAL.md` 2e16). With the
   canonical 2e16 everywhere the pinned set agrees with the engine; with the
   bare defaults it does not (the planner would drop 131 from the pin set of
   every K divisible by one of 59..113 -- exactly the top-ranked units). Audit
   and measurement below use 2e16 on both sides.
6. **Absolute calibration is off by 2.28x** (measured, see below), and in the
   opposite direction from the quoted model's error (b). Two separate
   constants cancel most of the way: `BASE_UNPINNED = 6.4e-14` double counts
   the bad primes that `rho` (which is a *post-sieve* conditional rate, not the
   density of Loeschian integers -- measured below) already contains, while the
   tier-A factor `prod_{p|D0}(1-1/p)` is missing. Only the absolute E58 is
   affected, not the ranking.

Everything in items 1-5 is implemented in `score2.py` (`v_new`). Where I was
unsure, both readings are implemented: `v_new(..., wbar=...)` can be given the
CUDA `wbar` (reported as `v_new_cudaw`); I trust `wbar = 0.5` for the CPU
engine and the CUDA `wbar` for the GPU engine, and the measurement here was
done on the CPU enumeration, where the corrected reading wins by a hair
(Spearman +0.892 vs +0.884).

## 3. Measurement

`probe.c` copies the stage-1 odometer and stage-2 tier-C bitmask of
`apsearch.c` verbatim (same pin set, same CRT walk, same masks) and replaces
stage 3 with a statistics pass: for every surviving window it tests all 58
terms and records per-position pass counts. Observable:

    M = log(win_ok / res) + sum_{k=0..57} log p_k

the log yield per stage-1 residue, which is exactly what `v` models.
`win_ok` excludes windows killed by a bad prime `q | d` dividing `a` (all 58
terms are then divisible by q; that is a `prod (1-1/q)` factor, not a per-term
rate, and must not be folded into `p_k`).

Why not count runs >= 30 directly: measured, a full-strength unit
(`apsearch --nterms 58 --kmin 32841 --shifts 1 --b2 10000 --report 28`) produced
**59 survivors and a longest run of 16 in 25 s** -- 2.4 windows/s. With the
measured per-term rate 0.709, P(run >= 30 in a window) = 29*0.709^30*(1-0.709)^2
= 8e-5 and P(run >= 34) = 1.8e-5, so one run >= 30 per ~1.4 core-hours and one
run >= 34 per ~6 core-hours. Counting those to +-30% would cost ~10 core-hours
*per unit*; there is no laptop budget in which run counting discriminates
between units. Lowering b2 to 200 instead raises the survivor rate from
2.8e-7 to 1.2-4.2 windows per stage-1 residue and makes the per-term rate
measurable to 0.3% in 10 s, at 58x the leverage (the quantity being ranked is a
product of 58 per-term factors, so a 0.3% error in the rate is 18% in the
58-term yield).

* 100 units, `maxres = 400` stage-1 residues each (deterministic work, the
  machine was loaded at load-average 80 by other jobs): K from 2015 to 59893,
  shifts {0,1,2,4,8,16,32}, 40% of them with a bad prime factor in 59..113
  (the units the score ranks highest). 31k-78k term tests per unit.
* Statistical power: sigma(M) = 0.300 (median) per unit, against a score spread
  of 2.92 in log (18x) across the sample -- signal/noise ~10 per unit, and the
  100-unit regression resolves the slope to +-0.049.
* Independence caveat: `M` multiplies the 58 measured per-position rates, which
  is the model's own assumption. A positive correlation between terms of one
  window would make the true 58-yield larger than `M`; the measurement cannot
  rule that out and it is the one place where the absolute 2.28x could move.

Checks that the observable is sound:
* b2 invariance: unit (K=10050, s=0) measured at b2=200 gives M = -35.289
  +- 0.305, at b2=10000 M = -35.522 +- 0.286 (diff 0.23, inside noise).
* `rho` amplitude: at b2=10000 the measured per-term rate is 0.7090 +- 0.0034
  at mean `ln t` = 39.36 (K=10050) and 0.6948 +- 0.0034 at `ln t` = 41.02
  (K=33611), against `rho` = 0.7180 / 0.7037. Same ratio 0.9875 / 0.9873 at
  both sizes: the `(ln t)^-0.485` *shape* is right, the amplitude is 1.26% per
  term too high, i.e. 0.48x over 58 terms. (This also proves `rho` is the
  conditional rate for sieved windows, not the density of Loeschian integers:
  a random integer of that size is Loeschian with probability ~0.11.)

## 4. Verdict and what it means

The ranking is effectively right. Both scores correlate with measured yield at
Spearman ~0.88-0.89, the old score's scale is 1.007 +- 0.049, and the top 10%
of the residue budget selects **the same 22 units** under either score
(ratio 1.000). Even an oracle that ranks by the *measured* yield gains only
1.08x at a 10% budget on this sample, and part of that is noise-fitting; so
there is at most a few percent to be had from re-ranking, not 2x. The correct
reading of the observed 16-29x yield falloff across the plan is that the score
is successfully sorting a real 29x spread, not that it is mis-sorting one.

Two things worth carrying forward anyway:
* use `wbar = 0.5` when the plan is to be run on the CPU engine
  (`score2.py --wbar` default); it is worth up to 2.4x on individual small-K
  shift-0 units even though it is worth nothing to the ordering as a whole.
* the project's quoted expected yields are low by 2.28x (sem 0.06): the
  `E58 = 1.777` at 2e16 residues printed by `plan_units.py` corresponds to a
  measured `E58 = 4.06`. Pin `MODCAP`/`--modcap` to 2e16 on both sides so the
  pin sets agree.

## Files

* `probe.c` -- per-unit yield measurement (build: `cc -O3 -march=native -std=gnu11 -o build/probe experiments/2026-10-03-planaudit/probe.c -lm`)
* `score2.py` -- corrected scorer (`v_new`) + `v_old` wrapper around `plan_units.py`
* `plan_compare.py` / `plan_compare.txt` -- old vs new ranking over 3000 plan units
* `analyse.py` / `analysis.txt` -- correlation of both scores with the measurement, DECIDES statistic
* `measured.jsonl` (100 units, b2=200), `measured_b2_10000.jsonl` (2 units, full-strength b2)
* `abs_rho_check.txt`, `abs_calib.txt`, `b2_consistency.txt`, `joined.json`

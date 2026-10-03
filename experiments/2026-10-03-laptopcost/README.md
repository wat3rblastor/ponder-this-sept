# LOSS (not catastrophic): this 8-core laptop needs ~1.7e4 hours (~2 years) of continuous apsearch to expect one 58-term Loeschian AP; honest 1-sigma band 9e3 - 8e4 hours. This AGREES with the 8-GPU "0.4 per hour" figure (that figure implies 7.8e4 laptop-hours; measurement says 1.7e4, a factor 4.6 apart), so the heuristic formula's sign errors did NOT cost orders of magnitude.

Measured 2026-10-03 on the laptop (8 cores, 8 GB), `build/apsearch`, 762 s wall /
6092 core-seconds, n = 58, `--report 30`, `--b2 2000`, K = 20001..27010 in eight
disjoint blocks. Raw logs: `w0..w7.log`, `w0..w7.jsonl`. Analysis: `analyze.py`
-> `fit.txt`; histogram `histogram.tsv`; per-term measurements `q_scan.txt`,
`size_scan.txt`; measurement programs `mc_q.c`, `mc_win.c` (new files, nothing
in `src/` touched).

## The numbers

| quantity | measured |
|---|---|
| residues/s/core | 2.98e6 |
| residues/s, 8 cores | 2.38e7 |
| admissible windows per residue | 3.267e-5 |
| admissible windows/s, 8 cores | 778 |
| per-term pass probability inside a window, `q` | **0.64233 +- 0.00073** |
| pass probability of a term just outside the window, `p_out` | 0.4724 |
| P(engine reports a run >= 58) per window | 2.12e-11 |
| 58-term APs per laptop-core-hour | 7.4e-6 |
| **expected wall-clock, this laptop** | **1.69e4 h = 1.93 y** |

## What was actually measured, and what the engine's histogram means

Read `src/c/apsearch.c` stage 3 before trusting any histogram:

* `--report R` sets `want = R` and the early abort `if (cur + (nterms-k) < want) break;`.
  That abort can never fire before a run of length >= R has been found, so the
  logged hit set is **complete for run lengths >= R** and meaningless below R.
  With `R = 30` the complete range is L >= 30.
* The logged `n` is **not** the in-window run: after finding the best run inside
  the 58-term window, stage 3 extends outward in both directions
  (`while (a >= d && is_loeschian(a-d))` and the forward loop), so `n` is the
  length of a **maximal** run of the full AP. Counting is therefore of maximal
  runs, not of nested events. The terms outside the window are *not* protected
  by the tier-C sieve, which is why `p_out` (0.47) differs from `q` (0.64) and
  why the extension matters (it is worth a factor 2.99 at L = 58).
* 14 hit lines were produced; all 14 are distinct `(a, d, n)` triples, so no
  maximal run was double-counted. (The same maximal run *can* be reached from
  several admissible residues, so this was checked, not assumed.)
* The sieve only accepts windows in which **no** term is divisible by any bad
  prime <= `b2`. Solutions in which some term is divisible by a bad prime <= 2000
  to an even power are outside the engine's reach; the cost measured here is the
  cost *of this engine*, which is the right thing for comparison with the GPU run.

Empirical histogram (`histogram.tsv`), 592 070 admissible windows:

```
   L  obs N(>=L)  pred N(>=L)  obs/pred
  30          14        13.16      1.06
  31          11         8.22      1.34
  32           8         5.13      1.56
  33           5         3.20      1.56
  34           0         1.99      0.00
  35           0         1.24      0.00
```

Best run of the whole laptop run: n = 33, e.g. `a=131884807318232731`,
`d=8790083435234970` (K=23001), verified with
`python3 src/verify.py 131884807318232731 8790083435234970 33` -> PASS. Also
re-verified two older reported runs, n=37 and n=35 -> PASS, so the engine's
`is_loeschian` is not inflating run lengths.

## The decay: it is essentially constant, and here is why

The brief expected `r` to fall as L crosses bad primes. **It does not**, and the
reason is structural rather than empirical: inside an admissible window every
term is already free of *every* bad prime <= 2000, and every bad prime q > 2000
satisfies q > 58, so q can divide at most one of the 58 indices. The 58
Loeschian events are therefore independent Bernoulli(q) with the *same* q. The
only residual L-dependence is the slow growth of `ln t` across the window
(<1 %) and the window-edge/extension effect:

```
   L= 30  r=0.6247      L= 46  r=0.6108
   L= 34  r=0.6225      L= 50  r=0.6031
   L= 38  r=0.6197      L= 54  r=0.5903
   L= 42  r=0.6160      L= 57  r=0.5737
```

The decline from 0.625 to 0.574 is purely the shrinking number of window
alignments as L approaches 58, not arithmetic.

Because `r = q` is a near-identity here, `q` was measured directly instead of
being read off a 14-event tail: `mc_win.c` samples random `a` with the tier-A
conditions for the real `d = K*D0`, keeps the terms that the tier-C sieve would
have let through, and tests them with the project's own `is_loeschian`. Over
426 542 real AP terms at four K values, q = 0.64233 +- 0.00073 - fifty times more
precise than anything the tail can give, and it predicts N(>=30) = 13.2 against
14 observed.

### A factor-2 trap worth recording

A control that sampled *random integers* free of bad primes <= 2000 gave
q = 0.325, half the AP value (`size_scan.txt`). The missing factor is exactly 2
and is arithmetic, not statistical: `3 | d` and Loeschian numbers are never
2 mod 3, so apsearch pins `a = 1 mod 3` and **all 58 terms are 1 mod 3**. Half of
the integers coprime to 3 (the 2 mod 3 ones) can never be Loeschian, so
conditioning on 1 mod 3 doubles the density exactly. Any cost model for this
search that estimates the per-term density from unconditioned integers is wrong
by 2^58.

Scaling with term size, from `size_scan.txt` (unconditioned-mod-3 control,
multiply by 2 for the AP value): q ~ (ln t)^-0.548 over ln t = 35.5 .. 42.4.

## Uncertainty band, and why it is this wide

| anchor | laptop-hours |
|---|---|
| direct `q` measurement, +-1 sigma on q | 1.59e4 .. 1.80e4 |
| engine tail N(>=30) = 14, +-1 sigma Poisson (27 % -> dq/q = 0.9 %) | 9.2e3 .. 2.4e4 |
| no outward extension at all (pure q^58, a hard lower bound on yield) | 5.0e4 |
| 8-GPU 0.4/h figure scaled by measured throughput ratio 3.2e-5 | 7.8e4 |
| **adopted** | **1.7e4, band 9e3 .. 8e4** |

The band is only ~1 decade, not the several decades the brief anticipated,
because the extrapolation is *not* a 25-term blind geometric fit. It is a
one-parameter model whose parameter is measured to 0.11 % on 4.3e5 samples and
whose functional form (independence across the 58 indices) is forced by
"every bad prime > 2000 exceeds 58". The 25-term distance costs
58 * 0.0011 = 6.6 % in the final rate, not decades.

The genuinely weak alternative is to ignore the model and fit the tail alone.
The observed cumulative decay 30 -> 33 is r_obs = (5/14)^(1/3) = 0.709 against the
model's 0.624 (a ~1.2 sigma fluctuation on 14 events, driven by 5 hits at exactly
L = 33 and none above). Propagating r_obs^25 from L = 33 instead gives **225
hours**, and r = 0.66 gives 1382 hours. So a tail-only reading cannot rule out a
sub-1000-hour answer; it is just a far less constrained reading of far less data.
**The cheap decisive test is more of the same run**: at L = 34..38 the two
hypotheses differ by 3x to 30x in count, so ~8 more laptop-hours (about 25x the
present statistics) would separate them outright. That is the single highest-value
follow-up, and it needs no new code:

```
build/apsearch --nterms 58 --kmin <K> --kmax <K+9> --shifts 3 --b2 2000 \
               --report 30 --out w.jsonl    # 8 disjoint K blocks, K >= 30001
```

## Cross-check against the 8-GPU datapoint: agrees

8x RTX PRO 6000 measured 7.45e11 residues/s. Using the laptop's *measured*
3.267e-5 admissible windows per residue and 2.12e-11 per window gives
**1.85 expected 58-term APs per GPU-hour** against the **0.4** the project
quotes - a factor 4.6, the same direction and the same order. Since the two
engines differ in `modcap`, K range and residue bookkeeping, a factor ~5 is
within the inter-engine systematic. Conclusion: the 0.4/GPU-hour planning
number survives measurement. The sign errors in `corr(n)` evidently cancelled
to within an order of magnitude at n = 58.

Laptop/GPU throughput ratio: 2.38e7 / 7.45e11 = 3.2e-5, i.e. the 8-GPU rig is
~31 000 laptops. That ratio, not the probability model, is what makes the laptop
hopeless on its own.

## One model number that is NOT confirmed

`experiments/2026-10-03-structure/mathres/plan_units.py` uses a measured
per-term pass probability `rho(t) = 0.71992*(ln t/39.144)^-0.485`, i.e. 0.709 at
the sizes searched here. Measurement says 0.642. The ratio is 0.906 per term,
hence **0.0033 at n = 58** - that planner over-states the 58-term yield by ~300x.
This is not a close call: with rho = 0.709 the run above should have logged
N(>=30) = 202; it logged 14. Whatever `rho` in that script is measuring, it is
not the per-term pass probability inside an apsearch admissible window, and unit
*rankings* derived from it are fine but its absolute yields are not.

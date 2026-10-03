# LOSS: M(58) ≈ 3e17 (honest range 2e16 – 2e18), i.e. ~0.3x the heuristic's 1e18 — the heuristic is CONFIRMED, not 4x pessimistic; measured growth is x1.78–1.86 per extra term inside a plateau, rising to ~x2.7/term by n≈58.

The hoped-for outcome (M(58) ≤ 1e15) is **not supported**. The naive reading of the
data does give ~1e15, and that reading is wrong for a reason this experiment can
demonstrate from the data itself (§5). Re-aiming the project at ~1e15-sized terms
would be a mistake.

## 1. Decides

| quantity | value |
|---|---|
| **M(58) central estimate** | **3e17** (18 digits) |
| honest range (≈2σ, 3 independent routes + fit scatter) | 2e16 – 2e18 |
| ratio to the project's heuristic 1e18 | **0.3x** (range 0.02x – 2x) |
| measured growth factor per extra term, n=26..35 | **x1.788** (0.2523 dex/term) |
| ditto n=28..35 | x1.784 (0.2513 dex/term) |
| ditto n=22..31 (base-330 plateau only) | x1.72 (0.2355 dex/term) |
| modelled growth factor per extra term at n≈56..58 | **x2.6 – x2.7** (0.41–0.43 dex/term) |
| is the heuristic bias constant or n-dependent? | **constant** (slope of log10(M/heur) = −0.016 ± 0.018 dex/term, not significant) |

The single most useful artefact: **the growth rate is not constant.** It is
~0.25 dex/term at n=30 and ~0.42 dex/term at n=55. Any planning that assumes a
fixed factor per term calibrated at n≈35 will under-shoot n=58 by ~3 orders of
magnitude.

## 2. What was measured (exhaustive, provable)

A sweep over `d = base*m, m = 1..mmax`, all starts `a`, all terms `≤ N`, is a
**complete** search of every AP with last term `≤ min(N, mmax*base*(n-1))`.
So `M(n) = L` is PROVED iff a sweep found `L` and that sweep had
`N ≥ L` **and** `mmax*base*(n-1) ≥ L` **and** `base` divides the forced base of `n`.

Forced base (3 | d, and every bad prime p with 2p ≤ n must divide d):

| n range | forced base |
|---|---|
| 10–21 | 30 = 3·2·5 |
| 22–33 | 330 = 3·2·5·11 |
| 34–45 | 5610 = 3·2·5·11·17 |
| 46–57 | 129030 = 3·2·5·11·17·23 |
| 58+ | 3741870 = 3·2·5·11·17·23·29 |

Sweeps run (all finished; `mmax` is the *contiguous* prefix actually scanned):

| tag | base | N | mmax | proves |
|---|---|---|---|---|
| `nobase` | **1** | 1.6e5 | 7620 | n=14..24 with **no base assumption at all** |
| `b30` | 30 | 4.0e6 | 14814 | n=14..27 |
| `b330a` | 330 | 5.0e7 | 7215 | n=18..31 |
| `b330b` | 330 | 1.3e8 | 12707 | n≤31; **M(32) > 1.3e8** |
| `b5610a` | 5610 | 3.3e8 | 1782 | n=34,35; **M(36) > 3.3e8** |
| `b5610b` | 5610 | 6.0e8 | 1179 | (partial, stopped at budget) |

Two independent validations:
* The `nobase` sweep (base = 1, **no** forcing assumed) reproduces M(14)…M(24)
  exactly as the base-30/330 sweeps do. The forced-base argument is therefore
  verified empirically, not just asserted.
* `b5610a` independently re-derives the known record **M(35) = 311958331**
  at a = 219830911, d = 2709630.

Every reported progression passes `python3 src/verify.py a d n`.

## 3. The curve

`heur` = this experiment's re-derived model (§4); `heur_st` = the formula as
stated in the task brief, implemented verbatim. Ratio = measured / heur.

| n | base | status | M(n) | a | d | heur | ratio | heur_st | ratio_st |
|---|---|---|---|---|---|---|---|---|---|
| 14 | 30 | PROVED | 1117 | 727 | 30 | 1.0e2 | 11.2 | 1.0e2 | 11.2 |
| 15 | 30 | PROVED | 1147 | 727 | 30 | 1.0e2 | 11.5 | 1.0e2 | 11.5 |
| 16 | 30 | PROVED | 2383 | 1033 | 90 | 1.0e2 | 23.8 | 1.0e2 | 23.8 |
| 17 | 30 | PROVED | 2473 | 1033 | 90 | 3.74e2 | 6.62 | 1.71e2 | 14.5 |
| 18 | 30 | PROVED | 13051 | 2341 | 630 | 6.54e2 | 20.0 | 8.64e2 | 15.1 |
| 19 | 30 | PROVED | 13681 | 2341 | 630 | 1.26e3 | 10.9 | 1.33e3 | 10.3 |
| 20 | 30 | PROVED | 47119 | 21469 | 1350 | 2.64e3 | 17.9 | 2.26e3 | 20.8 |
| 21 | 30 | PROVED | 48469 | 21469 | 1350 | 5.17e3 | 9.37 | 3.67e3 | 13.2 |
| 22 | 330 | PROVED | 157831 | 137041 | 990 | 1.32e4 | 12.0 | 7.73e3 | 20.4 |
| 23 | 330 | PROVED | 158821 | 137041 | 990 | 3.49e4 | 4.55 | 1.93e4 | 8.23 |
| 24 | 330 | PROVED | 159811 | 137041 | 990 | 5.99e4 | 2.67 | 6.64e4 | 2.41 |
| 25 | 330 | PROVED | 160801 | 137041 | 990 | 1.09e5 | 1.47 | 1.12e5 | 1.44 |
| 26 | 330 | PROVED | 1673461 | 270961 | 56100 | 2.11e5 | 7.92 | 2.00e5 | 8.36 |
| 27 | 330 | PROVED | 1729561 | 270961 | 56100 | 4.41e5 | 3.92 | 3.84e5 | 4.50 |
| 28 | 330 | PROVED | 5719081 | 2235271 | 129030 | 1.01e6 | 5.64 | 8.09e5 | 7.07 |
| 29 | 330 | PROVED | 18530251 | 466051 | 645150 | 2.69e6 | 6.90 | 1.96e6 | 9.45 |
| 30 | 330 | PROVED | 19175401 | 466051 | 645150 | 4.51e6 | 4.25 | 5.84e6 | 3.28 |
| 31 | 330 | PROVED | 19820551 | 466051 | 645150 | 8.44e6 | 2.35 | 9.95e6 | 1.99 |
| 32 | 330 | **bracket** | (1.30e8, 2.796e8] | 95618329 | 5935380 | 1.63e7 | — | 1.74e7 | — |
| 33 | 330 | **bracket** | (1.30e8, 2.856e8] | 95618329 | 5935380 | 3.28e7 | — | 3.16e7 | — |
| 34 | 5610 | PROVED | 309248701 | 219830911 | 2709630 | 6.37e7 | 4.85 | 5.52e7 | 5.60 |
| 35 | 5610 | PROVED | 311958331 | 219830911 | 2709630 | 1.29e8 | 2.41 | 1.07e8 | 2.93 |
| 36 | 5610 | **lower bd** | > 3.30e8 | — | — | 2.73e8 | — | 2.14e8 | — |

n=32,33: the upper end comes from a base-5610 sweep (an incomplete sub-family for
those n, hence only an upper bound); the lower end `> 1.3e8` is a *proved*
exclusion from the complete base-330 sweep `b330b`. **Never read 2.796e8 as M(32).**

## 4. The heuristic, re-derived

The brief's formula is **structurally correct but wrong in three places**; both
versions are tabulated above.

For a set `D` of bad primes dividing `d` (so `d = 3·∏D·m`):

```
E_D(T) = T² · ∏_{p∈D}(1-1/p) / (18·(n-1)·∏D) · ρ_D(T)^n · ∏_{q bad, q∉D} corr_q(n)
ρ_D(T) = ρ₁(T) / ∏_{p∈D} p/(p+1)
ρ₁(T)  = 2·0.638909/√(ln T) · (1 + 0.36/ln T)
E(T)   = Σ over subsets D ⊇ forced(n)       ;  solve E(T) = 1
```

Corrections relative to the brief:
1. **Triangular factor ½.** `Σ_d (T − (n−1)d) ≈ T²/2`, not `T²`. The brief drops it.
2. **`∏_{p|d}(1−1/p)`.** When `p | d` one must also have `p ∤ a` (else every term
   needs `p²|`). That costs `(1−1/p)` per forced prime — a factor **3.2 at n=58**
   (p = 2,5,11,17,23,29). The brief drops it.
3. **`corr_q` for n/2 < q < n is wrong in the brief.** `1 − n/(q+1)` is the
   probability no term is divisible by q, valid only for `q ≥ n`. For `q < n`
   *some* term is always divisible by q, and the right value is
   `(1 − (n−q)/q)/(q+1)`. At n=58 the brief's expression is **negative** for
   q = 41,47,53 — it cannot be evaluated as written. (Derivation:
   `P(v_q even | q|x) = 1/(q+1)`, and `1−n/q + n/(q(q+1)) = 1−n/(q+1)` recovers
   the brief's form in the `q ≥ n` regime, so that part is exactly right.)

`ρ₁` was calibrated, not assumed: `density.c` sieves the Loeschian numbers and
counts those ≡1 mod 3. Measured density 0.35299 (1e6), 0.32532 (1e7), 0.30333
(1e8) vs asymptotic `2K/√(ln T)` = 0.34377, 0.31828, 0.29774 — the asymptotic is
2.0–2.7% low, fitted by the `(1+0.36/ln T)` term. (The `2K` rather than `K`
comes from conditioning on `a ≡ 1 mod 3`: `n ↦ 3n` bijects Loeschian onto
Loeschian-divisible-by-3, so 2/3 of Loeschian numbers are ≡1 mod 3 while only
1/3 of integers are.)

**I could not reproduce the project's claim that the heuristic predicts 1.3e9 at
n=35.** The brief's formula, implemented verbatim, gives **8.6e7** at n=35 (and
1.07e8 after a later refinement of which optional primes to consider) — i.e. it
is ~3x *optimistic*, not 4x pessimistic. The project's 1.3e9 is ~10x more
pessimistic than the brief's own formula; whatever produced it is a different
calculation. Its *shape* however matches mine almost exactly: 1.3e9→1e18 over
n=35→58 is 0.386 dex/term, mine is 0.394 dex/term. So the project heuristic ≈
(my model) × 10, and the n-dependence is agreed.

## 5. Why the ratio is constant, and why the naive fit's 1e15 is wrong

Bias of the model: `log10(M_measured / heur)` over the 12 proved points n≥22 has
slope **−0.0164 ± 0.0180 dex/term** — indistinguishable from zero. Mean ratio
**4.20x**, scatter 0.26 dex. So **the heuristic is biased by a roughly constant
factor, ~4x optimistic, with no detectable drift in n.** The scatter is not
noise: it is a sawtooth, because one lucky `(a,d)` with a long run serves many n
at once (d = 645150 supplies n = 29, 30 and 31 at nearly the same last term), so
the ratio decays across each plateau and jumps at the next one.

Three structure-aware routes to M(58):

| route | M(58) |
|---|---|
| (B) this model × constant measured bias 4.20x | 6.3e17  [1.9e17, 2.0e18] |
| (C) this model × bias trend extrapolated to n=58 | 2.0e17 |
| (D) project heuristic 1e18 × measured/predicted at n=35 (0.240) | 2.4e17 |
| geometric mean → **central estimate** | **3.1e17** |

And the structure-blind route, for contrast:

| route | M(58) |
|---|---|
| (A) plain `log10 M = a + b·n` on the proved points n≥26 | 6.0e14 [2.9e13, 1.2e16] |

**(A) is rejected, and the data say why.** It assumes b is constant. It is not:

```
modelled d(log10 M)/dn:  0.24 (n=24)  0.29 (n=30)  0.31 (n=35)
                         0.34 (n=38)  0.38 (n=44)  0.38 (n=50)  0.43 (n=58)
```

and the measured slope tracks the modelled slope to within the sawtooth
(measured 0.2523 dex/term over n=26..35 vs modelled 0.2875..0.3076 there).
The slope rises for two reasons that are *not* in the data range:
* every time n crosses a bad prime q, `corr_q` falls off a cliff (q goes from the
  "no term divisible by q" regime to the "some term always is" regime). Between
  n=35 and n=58 this happens at q = 41, 47, 53 — visible as the spikes at
  n=41,47,53 in the growth table, and *entirely outside* the measured range.
* the forced base jumps ×23 at n=46 and ×29 at n=58.

Trap (a) in the brief is real, but note the discontinuities are *mild*: the ×17
jump at n=34 only costs ×1.94 (heur(33)→heur(34)), because 17 was already being
paid for voluntarily in the optimal d below n=34 (d=56100 at n=26, d=645150 at
n=29 both contain 17). Measured: M(21)→M(22) across the ×11 jump is ×3.26, vs
×1.72 within-plateau — one extra factor ~1.9. The corr-cliffs, not the base
jumps, are what bend the curve.

## 6. Implication for the project

The searchable region is **not** much denser than modelled. 58-term progressions
should first appear near **3e17** — an 18-digit last term, ~0.5 digits *below*
the assumed 1e18, which is noise at this precision. The project's target size
stands. The one piece of genuinely good news is small and quantified: the
heuristic is ~4x pessimistic relative to a correctly-specified model of the same
shape (routes B–D), so the search bound can be relaxed by ~0.5 dex, not by 3.

What would actually move this: a proved point at n ≥ 40. The cost of proving
M(n) is `M(n)²/((n−1)·base)` sieve-cells ≈ 2.5 ns each, so n=38 (M ≈ 1.5e9)
costs ~8 core-hours and n=40 (M ≈ 5e9) ~100 core-hours — the first n ≥ 40 point
is the cheapest way to test the corr-cliff behaviour that drives everything above.

## 7. Files

* `mcurve.c` / `mcurve` — `src/c/loeschsearch.c` patched to record, in the *same*
  O(N) pass, the smallest start `a` for **every** run length, not just one
  `--nmin`. This is what made the whole curve affordable: one sweep yields all n.
  Emits `"mina":{k: a}` per d.
* `density.c` / `density` — empirical Loeschian density calibration.
* `sweep.sh`, `sweep2.sh` — parallel m-range splitting (8 cores).
* `agg.py` — per-sweep aggregation with the completeness test.
* `heuristic.py` — the model of §4 (both variants).
* `analyze.py` — the table, the fits, the growth decomposition.
* `analysis.txt` — full output of `analyze.py`.
* `table.json` — machine-readable M(n).
* `raw/*.jsonl` — raw per-d sweep output (~30k d values); `raw/*.log` — run logs.

Reproduce: `make bin && cc -O3 -march=native -o mcurve mcurve.c && ./sweep.sh … && python3 analyze.py`

# Forensics on the long Loeschian APs — verdict

**CLEAN LOSS: long Loeschian progressions are statistically indistinguishable from random
admissible ones. Every per-term feature tested came out within 1.06x of a properly matched
control (3486 real terms vs 7320 control terms, no p-value below 0.3), and the one large
deviation that does exist — the ~20x advantage of steps `d` divisible by the bad primes
`q <~ n` — is fully accounted for by the local-density product the repo already uses, with
no residual bias (model-corrected excess 0.90x for `23 | d`, 0.54x for `29 | d`).**

Nothing here is a new search bias to exploit. The strongest *statistically significant*
deviation from the matched control anywhere in this study is the hit-prime count,
obs/exp = 0.90x (137 observed vs 152.8 expected, z = -1.3, not significant), and the largest
point estimate on any per-term feature is `smooth1e4` at 1.06x (p = 0.37).

## What was examined

Five verified records (all re-verified with `src/verify.py` → PASS today):

| name | a | d | n |
|---|---|---|---|
| 55a | 11687581876345393 | 2·3³·5·11·17·23·29·41·47·53·59 | 55 |
| 55b | 296246969176050787 | 2²·3·5²·11·17·23·29·41·47·53²·59 | 55 |
| 47 | 2646171143023357 | 2·3·5²·11·17·23·29·41²·47·53 | 47 |
| 43 | 5413537078288507 | 2⁸·3·5·11·17·23·29·41·47·53 | 43 |
| 35opt | 219830911 | 2·3²·5·7·11·17·23 | 35 |

Plus a **new exhaustive corpus** built with `build/loeschsearch` (`corpus_5610.jsonl`,
`corpus_129030.jsonl`, merged `corpus_all.jsonl`): 174 steps `d`, each scanned exhaustively
over *every* start `a` with all terms ≤ 2·10⁹, giving the provably-longest run for that `d`
plus the exact count of starts whose run reaches each length k ≥ 25. 122 of those optimal
runs have n ≥ 27; 5 were spot-re-verified with `src/verify.py` → PASS. The corpus reproduces
the known 35-optimum exactly (d = 2709630 → a = 219830911, n = 35), which validates the
pipeline.

## The control

For a record `(a, d, n)`: sample `a' ≡ a (mod M)` with `M = Π_{bad p | d} p^(v_p(d)+2)`,
same magnitude as `a`, take a random index k, keep `t' = a' + k·d` **if it is Loeschian**.
This matches: the size, the arithmetic class mod `d`, and the local conditions at every bad
prime dividing `d` (which is what admissibility means). Null hypothesis: *the terms of a long
AP look like independent random Loeschian numbers of that size in that class.* An unmatched
control (e.g. random Loeschian integers of the same size) would have shown a huge fake
"squarefull part" effect, since bad primes are forced to even powers by Loeschianness itself.

## Per-term features — corpus (122 optimal runs, n ≥ 27)

| feature | real | control | ratio | p (2-sided) |
|---|---|---|---|---|
| `is_auto` (3·c·d·t square, c ≤ 60) | 0/3486 | 0/7320 | — | see below |
| `is_auto2` (c·d·t square, c ≤ 60) | 0/3486 | 0/7320 | — | see below |
| `is_csquare` (t = c·square, c ≤ 60) | 8/3486 = 0.0023 | 21/7320 = 0.0029 | 0.80x | 0.66 |
| `good_sq_big` (square part from good primes > 1) | 223/3486 = 0.0640 | 451/7320 = 0.0616 | 1.04x | 0.58 |
| `smooth1e4` (largest prime ≤ 10⁴) | 269/3486 = 0.0772 | 535/7320 = 0.0731 | 1.06x | 0.37 |
| `smooth1e6` (largest prime ≤ 10⁶) | 1227/3486 = 0.3520 | 2545/7320 = 0.3477 | 1.01x | 0.61 |
| `log P(t)/log t` (largest prime factor) | 0.7868 | 0.7911 | −0.02 sd | z = −1.30 |
| `log rad(t)/log t` (squarefree kernel) | 0.9812 | 0.9809 | +0.01 sd | z = +0.36 |
| distinct prime factors | 2.1024 | 2.0870 | +0.02 sd | z = +1.03 |
| prime factors with multiplicity | 2.3368 | 2.3202 | +0.01 sd | z = +0.85 |
| distinct bad primes | 0.1629 | 0.1643 | −0.00 sd | z = −0.21 |

## Per-term features — the five records (235 terms)

| feature | real | control | ratio | p |
|---|---|---|---|---|
| `is_auto` | 0/235 | 0/20000 | — | — |
| `is_csquare` | 0/235 | 0/20000 | — | — |
| `good_sq_big` | 6/235 = 0.0255 | 528/20000 = 0.0264 | 0.97x | 1.0 |
| `smooth1e4` | 4/235 = 0.0170 | 439/20000 = 0.0220 | 0.78x | 0.82 |
| `smooth1e6` | 11/235 = 0.0468 | 1259/20000 = 0.0630 | 0.74x | 0.38 |
| `log P(t)/log t` | 0.7955 | 0.7946 | +0.005 sd | z = 0.07 |
| `log rad(t)/log t` | 0.9972 | 0.9969 | +0.018 sd | z = 0.27 |
| distinct primes | 2.0255 | 2.0062 | +0.02 sd | z = 0.34 |

**Power.** With 3486 real corpus terms and a control rate of 0.06 we detect a 1.25x shift at
p = 0.05; for the 0.003-rate features, 2.0x. The 235 record terms alone only detect 1.9x /
8x respectively — so the records on their own are genuinely underpowered, which is exactly
why the exhaustive corpus was built. Nothing exceeded 1.06x in the well-powered sample.

## Automatic families: dead, and provably so (not just underpowered)

`3·c·d·t` is a perfect square iff `t = S_c · u²` where `S_c = squarefree kernel of 3·c·d`.
Because `d` is large and nearly squarefree, `S_c` is itself enormous, so the automatic-family
values are extraordinarily sparse near the record terms. Exact expected number of
automatic-family terms per record under the null (`autofamily.py`, c ≤ 60):

| record | min S_c | E[auto terms in window] | observed |
|---|---|---|---|
| 55a | 127386974990 | 9.9·10⁻¹² | 0 |
| 55b | 240352783 | 2.5·10⁻¹¹ | 0 |
| 47 | 10713791 | 1.3·10⁻⁹ | 0 |
| 43 | 1158063409 | 1.0·10⁻¹⁰ | 0 |
| 35opt | 16422 | 1.6·10⁻⁴ | 0 |

So the Monte-Carlo version of this test could never have had power — but it does not matter:
the analytic null says the expected count is ~10⁻¹¹ and the observed count is 0. The
automatic-family identity can only bite for terms of order `S_c·u²` with small `u`, i.e.
terms *smaller* than `d`, which no long AP has. **Do not revisit this** (consistent with the
already-recorded pentagonal/automatic-terms negative).

## Hit indices — exactly as forced, no extra structure

For a bad prime `q ∤ d`, the hit index is `k0 = −a·d⁻¹ mod q`. Two or more hits in the window
is fatal (it would force `q | d`); one hit survives iff `q² | t_{k0}`. The **correct**
conditional null, over `k0` uniform mod `q`, is
`P(hit | AP all-Loeschian) = B / (qA + B)` with `A = #{k0 : 0 hits}`, `B = #{k0 : 1 hit}`.

| sample | primes `q < 1000`, `q ∤ d`, with a hit | expected (unconditional) | expected (conditional) | ratio |
|---|---|---|---|---|
| 122 corpus optima | 137 | 1294.2 | 152.8 | **0.90x** (z = −1.3) |
| 5 records | 1 | 62.7 | ~0.2 | — (the single hit is 29²·… in 35opt, forced: 29 < n and 29 ∤ d) |

Relative position of the hit index inside the window: mean 0.482 (n = 137) against a uniform
null of 0.500 — **no edge clustering**. Beware: against the *unconditional* null the deficit
looks like a 60x effect with z = −9.4; that is an artifact of the wrong null, and is the one
trap in this dataset. It is entirely explained by "one hit requires `q² | t`".

## Which `d` are good — the only large effect, and it is already known

Exhaustive counts of starts reaching length k, regressed on `log d` plus an indicator per bad
prime dividing `d` (174 steps, `corpus_analysis.py` Part 1, k = 26):

| extra bad prime in `d` | measured gain | local-density model prediction | residual excess over model |
|---|---|---|---|
| 23 | **x18.2** (±0.14 in log) | x25.3 | **0.90x** (not significant) |
| 29 | x3.56 (±0.26) | x7.2 | 0.54x |
| 41 | x2.48 (±0.33) | x2.6 | — |
| 47 | x2.71 (±0.33) | x2.1 | — |
| 53 | x1.60 (±0.33) | x1.9 | — |
| `log d` | −0.016 ± 0.045 (null) | — | — |

The model (`density_model.py`): `P_q = (A + B/q)/q` for `q ∤ d`, `≈ 1 − 1/q` for `q | d`,
multiplied over bad `q ≤ 300` and by the number of available starts. `log(obs/pred)` has a
scatter of only ±0.66 (1 sd, i.e. x1.9, and the counts themselves are Poisson with means
1–100, so most of that is shot noise), with no significant residual for any prime. **The
choice of `d` is fully explained by the per-prime local densities; there is nothing extra.**

This also answers the "both 55s have `59 | d`, which is not forced" question: it is not a
coincidence and not a hidden structure, it is the density model. For n = 55, `P_59(59 ∤ d)`
= 0.567 vs 0.983 with `59 | d` — a x1.7 gain per the model, and the same reasoning makes
`29, 41, 47, 53 | d` worth x14.2, x2.9, x2.3, x2.0. The pattern "`d` is divisible by every bad
prime up to the first bad prime exceeding n" holds for 55a, 55b, 47 and 43 (and the 43 even
has the next one, 47, as well); the only violation is 35opt, which is a *last-term-minimising*
optimum where the extra factor 29 in `d` would have cost more in term size than it bought.
The `7 | d` in 35opt is a good prime and is, as expected, irrelevant — it buys nothing in any
of the statistics above; it is simply what the minimal last term happened to need.

## Relations between records — none

All pairwise `gcd(a_i, a_j) = 1` except `gcd(a_43, d_35) = 7`; `gcd(a, d) = 1` for every
record; the `gcd(d_i, d_j)` are just the shared forced primes. `a/d` (the "centre index")
is 57.6, 24.8, 33.8, 110.7, 81.1 — no pattern. 55a's `a` factors as `43 · 907813 · 299405527`,
55b's and 47's are prime, 43's is `7 · 773362439755501`. No ratio, no shared structure, no
rescaling relation. The two 55s are genuinely independent solutions.

## Everything tested (so nobody repeats it)

1. largest prime factor (absolute and `log P/log t`) — null
2. smoothness at 10⁴ and 10⁶ — null
3. number of prime factors, distinct and with multiplicity — null
4. squarefree-kernel ratio `log rad(t)/log t` — null
5. square part contributed by good primes — null
6. `t = c·square` for c ≤ 60 — null
7. `3·c·d·t` a perfect square (automatic families), c ≤ 60 — null, and analytically impossible
8. `c·d·t` a perfect square, c ≤ 60 — same
9. number of distinct bad primes per term — null
10. position (index k) of the few structured terms — no clustering
11. hit index `k0 = −a/d mod q` presence, against the correct conditional null — 0.90x
12. hit index relative position in window — 0.482 vs 0.500
13. which bad primes divide `d`, and to what power — fully explained by local density
14. good primes dividing `d` (e.g. `7 | d` in 35opt) — irrelevant, as theory says
15. `log d` effect after conditioning on prime structure — null
16. gcds / ratios / shared factors between all five records — nothing

## Files

- `forensics.py` — per-term features + matched control for the five records
- `compare.py` — real-vs-control tables with binomial p-values (`python3 compare.py REAL CTRL`)
- `structure.py` — `d` factorisation, forced vs extra bad primes, hit indices, record relations
- `autofamily.py` — exact analytic null for the automatic-family feature
- `corpus_analysis.py` — corpus: which `d` are good (Part 1), per-term forensics on the 122
  exhaustive optima (Part 2), hit-index statistics (Part 3)
- `density_model.py` — local-density model vs the exhaustive counts
- `corpus_5610.jsonl`, `corpus_129030.jsonl`, `corpus_all.jsonl` — the exhaustive corpus
- `real_terms.json`, `ctrl_terms.json`, `corpus_real_terms.json`, `corpus_ctrl_terms.json`
- `results.txt` — full captured output of every script

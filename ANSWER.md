# Longest Loeschian AP found so far: **n = 36 terms**, a = 73415383, d = 37418700

**This does NOT yet clear the mandatory bar of n ≥ 58** (published record: 57). Search is in
progress against the 2026-10-03 05:27 CDT deadline; this file is rewritten whenever the record
changes, so the headline above is always the current best *verified* result.

| goal | target | status |
|---|---|---|
| **G3 (mandatory)** | n ≥ 58 | **not yet** — best verified n = **36** |
| G2 (bonus `*`) | n ≥ 42 | not yet |
| G1 (secondary) | 35 terms, minimal last term | 35-term AP found, last term `1345651183`; **not optimized** (deferred by GOAL.md §1 until G3 is secured) |

## The progression

```
a = 73415383
d = 37418700
n = 36 terms  (indices k = 0 .. 35)
last term a + 35d = 1383069883
```

All 36 terms:

```
  73415383   110834083   148252783   185671483   223090183   260508883
 297927583   335346283   372764983   410183683   447602383   485021083
 522439783   559858483   597277183   634695883   672114583   709533283
 746951983   784370683   821789383   859208083   896626783   934045483
 971464183  1008882883  1046301583  1083720283  1121138983  1158557683
1195976383  1233395083  1270813783  1308232483  1345651183  1383069883
```

Regenerate with: `python3 -c "a,d=73415383,37418700; print([a+k*d for k in range(36)])"`

## Verification

```
$ python3 src/verify.py 73415383 37418700 36 --maximal
...
  maximality: a-d = 35996683 is NOT Loeschian
  maximality: a+36*d = 1420488583 is NOT Loeschian
  run is maximal in both directions

OVERALL: PASS  (n=36, a=73415383, d=37418700)
```

`src/verify.py` factors every term from scratch in Python ints (its own Miller–Rabin +
Pollard–Brent, no sieve, no fixed-width arithmetic) and checks that every prime `p ≡ 2 (mod 3)`
occurs to an even power. The run is **maximal at both ends**, so n = 36 is the true length here,
not a truncation. An independent cross-check by a second implementation is still outstanding
(required by GOAL.md §8 only for the final `n ≥ 58` claim); three implementations already agree
on the underlying Loeschian test (`make check`: factor-based test vs. brute-force
`x²+xy+y²` enumeration vs. enumeration sieve, over all n ≤ 20000).

## Why it works

`d = 37418700 = 2² · 3 · 5² · 11 · 17 · 23 · 29`. Two structural facts force this shape:

- A Loeschian number is never `≡ 2 (mod 3)` (for `n` coprime to 3,
  `n ≡ (−1)^(#bad prime factors) mod 3`, and Loeschian means every bad exponent is even). So
  `3 | d` is forced and every term is `≡ 1 (mod 3)` — here `a ≡ 1 (mod 3)`. With `3 ∤ d` the
  longest run found anywhere is 2.
- Each "bad" prime `p ≡ 2 (mod 3)` with `2p ≤ n` must divide `d`: two terms divisible by `p`
  would each need `p² |`, but their difference is `j·p·d` with `j < p`. For n = 36 that forces
  `2, 5, 11, 17 | d`; 23 and 29 are in `d` too, which removes their cost as well.

So no term is divisible by any of 2, 5, 11, 17, 23, 29, and every term sits in the residue class
`a mod d`, where the per-term probability of being Loeschian rises from ~0.14 to ~0.6.

## How it was found

`experiments/2026-10-02-calibration`, ~50 seconds of one core. `src/c/loeschsearch.c` sieves a
Loeschian bitmap up to 2·10⁹ by enumerating `x²+xy+y² ≤ N`, then for each step `d` computes the
run length for **every** start `a` in one descending pass (`run(a) = L[a] ? run(a+d)+1 : 0`,
O(N) time and O(d) memory), so a single pass covers all starts at once. The key pruning idea is
the paragraph above: only steps divisible by `3 · 2·5·11·17·23·29` were tried, which is a factor
of ~3.7·10⁶ fewer `d` values than a blind sweep.

This result is a by-product of calibrating the cost model, not a serious attempt at the record —
the serious engine (`src/c/apsearch.c`, a Wróblewski-style two-stage CRT/bitmask search) is what
is being run at the n = 58 target.

## Where to watch progress

- `records.json` — machine-readable current records.
- `PROGRESS.md` — dated log, one entry per campaign, including what was ruled out.
- `experiments/<date>-<slug>/*.jsonl` — live per-work-unit output from running campaigns; lines
  with `"hit":true` are candidate progressions, each of which must still pass `src/verify.py`
  before it counts.

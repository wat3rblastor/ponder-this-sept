# Longest Loeschian AP found so far: **n = 43 terms**, a = 5413537078288507, d = 48916598396160

**This does NOT yet clear the mandatory bar of n ≥ 58** (published record: 57). The search is
still running against the 2026-10-03 05:27 CDT deadline; this file is rewritten whenever the
record improves, so the headline is always the current best *verified* result.

| goal | target | status |
|---|---|---|
| **G3 (mandatory)** | n ≥ 58 | **not yet** — best verified n = **43** |
| **G2 (bonus `*`)** | n ≥ 42 | **CLEARED** — n = 43, verified |
| G1 (secondary) | 35 terms, minimal last term | best known last term `1345651183` (a=73415383, d=37418700); **not optimized** — deferred by GOAL.md §1 until G3 is secured |

## The progression

```
a    = 5413537078288507
d    = 48916598396160          ( = 128 * 3*2*5*11*17*23*29*41*47*53 )
n    = 43 terms  (indices k = 0 .. 42)
last = a + 42d = 7468034210927227      (16 digits)
```

All 43 terms:

```
5413537078288507  5462453676684667  5511370275080827
5560286873476987  5609203471873147  5658120070269307
5707036668665467  5755953267061627  5804869865457787
5853786463853947  5902703062250107  5951619660646267
6000536259042427  6049452857438587  6098369455834747
6147286054230907  6196202652627067  6245119251023227
6294035849419387  6342952447815547  6391869046211707
6440785644607867  6489702243004027  6538618841400187
6587535439796347  6636452038192507  6685368636588667
6734285234984827  6783201833380987  6832118431777147
6881035030173307  6929951628569467  6978868226965627
7027784825361787  7076701423757947  7125618022154107
7174534620550267  7223451218946427  7272367817342587
7321284415738747  7370201014134907  7419117612531067
7468034210927227
```

Regenerate with:
`python3 -c "a,d=5413537078288507,48916598396160; print([a+k*d for k in range(43)])"`

## Verification

```
$ python3 src/verify.py 5413537078288507 48916598396160 43 --maximal
...
  maximality: a-d = 5364620479892347 is NOT Loeschian
  maximality: a+43*d = 7516950809323387 is NOT Loeschian
  run is maximal in both directions

OVERALL: PASS  (n=43, a=5413537078288507, d=48916598396160)
```

`src/verify.py` factors all 43 terms from scratch in Python ints (its own Miller–Rabin +
Pollard–Brent, no sieve, no fixed-width arithmetic) and checks that every prime `p ≡ 2 (mod 3)`
occurs to an even power. The run is **maximal at both ends**, so 43 is its true length. Three
independent implementations agree on the underlying Loeschian test for every `n ≤ 20000`
(`make check`), and the search engine was separately validated by having it rediscover a
28-term AP found by other means.

## Why it works

`d = 48916598396160 = 2^8 · 3 · 5 · 11 · 17 · 23 · 29 · 41 · 47 · 53`. Two facts force the shape:

- **A Loeschian number is never `≡ 2 (mod 3)`.** For `t` coprime to 3,
  `t ≡ (−1)^(number of bad prime factors, with multiplicity) (mod 3)`, and "Loeschian" means
  every such exponent is even, forcing that count even. So `3 | d` is forced and every term is
  `≡ 1 (mod 3)` — measured: with `3 ∤ d` the longest run found anywhere is **2**. This one fact
  raises the per-term probability from ~0.14 to ~0.6.
- **Bad primes below `n` want to divide `d`.** For `2p ≤ n` it is a theorem: two terms divisible
  by `p` would each need `p² |`, but their difference is `j·p·d` with `j < p`. For `n/2 < p < n`
  it is not forced but costs a factor `~1/(p+1)`, so `d` carries 41, 47 and 53 too.

So no term is divisible by any of 2, 5, 11, 17, 23, 29, 41, 47, 53, and `a ≡ 1 (mod 3)`. On top
of that, `a` was chosen so that none of 59, 71, 83, 89, 101, 107, 113 divides *any* of the 43
terms — each of those primes can hit at most one term in the window, and `a`'s residue pushes
the hit outside it.

## How it was found

`experiments/2026-10-03-ap58`, ~20 minutes on 8 cores, work unit `K = 128`.
`src/c/apsearch.c` is a two-stage search after Wróblewski's AP26 engine, adapted from primes to
Loeschian numbers. The key pruning idea: the admissible residues of `a` for the primes
59…113 form a product of arithmetic progressions, so they are **enumerated directly by nested
additive loops** and an inadmissible `a` is never even constructed — worth about `3·10⁴` over
sieving `a`. A second stage then settles 64 candidates per modulo using precomputed 64-bit
masks for 148 further primes. Combined, roughly `10⁶` less work per candidate than scanning.

## Honest status on n ≥ 58

Measured throughput is 3.65·10¹⁴ candidate `(a,d)` pairs per second on 8 cores, and the
calibrated model (which reproduces the observed `n=31` and `n=36` frequencies within a factor
of 2) puts the probability of a random candidate yielding 58 terms at ~1·10⁻²¹ at this term
size. That is ~24 CPU-days. **n = 58 is therefore out of reach in a 12-hour budget on this
machine**; expected reach is n ≈ 54–55. See PROGRESS.md for the derivation and for the
approaches already ruled out.

## Where to watch progress

- `records.json` — machine-readable current records, with full history.
- `PROGRESS.md` — dated log: what was tried, what was found, what was ruled out.
- `experiments/2026-10-03-ap58/w*.jsonl` — live per-work-unit output; `"hit":true` lines are
  candidates, each of which must still pass `src/verify.py` before it counts.

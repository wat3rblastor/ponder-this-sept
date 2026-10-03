**WIN — 1.02 work units / 1.36 ns per admissible `a` emitted below X = 2.18e13 for a pin set of selectivity 1/delta = 118150 (pins 59..107), and 2.09 units / 6.7 ns at selectivity 242745 (pins 59..113). That is ~100x better than the "decisive WIN" threshold of 100 operations per output, and 1e4x better than the naive baselines. 45,000,000 emitted values were verified against an independent verifier with 0 failures, plus exact set equality against full brute-force scans in 30 configurations.**

**But read §5 before acting on it.** The algorithm wins decisively; the *200x project inference does not follow*, because the thing that is scarce at small term size is not enumeration speed but **supply**. The honest payoff is (a) a free ~1.2-1.5x for today's engine, and (b) exhaustive access to the small-term region, worth ~1.6x integrated over a 34-core-hour sweep, not 200x. The decisive algorithmic fact and the disappointing project consequence are both below.

---

## 1. The algorithm: wheel-under-X + table filter

Problem, as posed: given pin set Q of bad primes q > 58 with q∤d, enumerate every
`a < X` with `a ≡ j·d (mod q)`, `1 ≤ j ≤ q-58`, for all q ∈ Q — in the hard regime
X ≪ M = ∏Q.

Split Q = Q1 ∪ Q2 with **M1 = ∏Q1 ≤ X**:

- **Q1 — wheel.** Enumerate the ∏(q-58) admissible residues mod M1 by a nested
  additive CRT walk: substituting `a = d·u (mod M1)` makes the admissible set the
  box `u mod q ∈ [1, q-58]`, so `a = Σ_q e_q·j_q·d (mod M1)` is built with one add
  and one conditional subtract per level. **O(1) amortised per residue, no
  multiply, no modulo.**
- **Lift.** Each residue r gives `⌊(X-r)/M1⌋+1` values `r + t·M1 < X`.
- **Q2 — filter.** Reject with one byte-table lookup per leftover prime; `a mod q`
  is maintained additively across the lift loop.

Cost per emitted `a` = `∏_{q∈Q2} q/(q-58)` cheap ops — **O(1) when Q2 is small, and
exactly 1 when Q2 is empty**. This is approach (a) (meet-in-the-middle) with the
split chosen so one side fits under X, which makes the sort/binary-search
unnecessary; it subsumes approach (c) (the pruned CRT walk is the `Q2 = ∅`,
`M1 > X` special case, measured below as the loser).

I did not need a lattice (approach (b)). The reason the problem looked hard is
that both obvious baselines are bad in *different* directions, and the split
avoids both at once.

## 2. Measured cost per output (bench.out, single core, Apple M-series, `cc -O3`)

`work_per_out` = (wheel leaves + lift iterations) / outputs. It is the stable
metric; `ns_per_out` moves ±20% run to run at the low end.

| pins | X (term floor) | 1/delta (pins only) | 1/delta (with 2,3,5) | work/out | ns/out |
|---|---|---|---|---|---|
| 59..101 | 2.18e13 | 7214 | 54106 | **1.00** | 1.22 |
| 59..107 | 2.18e13 | 15753 | **118150** | **1.02** | **1.36** |
| 59..113 | 2.18e13 | 32366 | **242745** | **2.09** | 6.7 |
| 59..131 | 2.18e13 | 58082 | 435612 | 3.74 | 20.0 |
| 59..173 | 2.18e13 | 380115 | 2.85e6 | 24.50 | 382 |
| 59..113 | 1e14 | 32366 | 242745 | 1.38 | 3.20 |
| 59..113 | 1e15 | 32366 | 242745 | 1.04 | 1.56 |
| 59..113 | 7.26e16 | 32366 | 242745 | 1.00 | 1.26 |

The full uncapped headline run: pins 59..113, X = 57·D0 = 21783172723290 →
**673,025,293 outputs in 5.86 s = 8.7 ns each** (matching the predicted
delta·X = 6.73025e8 to 5 significant figures).

**Baselines, for contrast.**

- *Naive full-residue enumeration* (= pruned CRT walk, `--split 7`): cost per
  output measured at 37732 / 3779 / 379 for X = 1e9 / 1e10 / 1e11, i.e. exactly
  the predicted **M/X** (37800 / 3780 / 378). The pruning does *not* bite: a
  partial CRT value gives no information about whether the final value lands
  under X, so the walk pays for the whole residue class count ∏(q-58).
- *Sieving [0,X)*: 1/delta per output = 32366 (pins only) by construction.
- Both are >1e4 per output, which is the briefing's "decisive LOSS" band. The
  split beats them by ~1e4x.

So: the briefing's premise — that small-relative-to-modulus admissible values
cannot be produced at O(1) each — is **false**. The answer to the question as
asked is yes, at 1-2 operations each.

## 3. The structural rule that falls out

Pinning one extra prime q costs a factor `q/(q-58)` in enumeration **and** gains
a factor `((q+1)/q)^58` in the 58-term pass probability (q can no longer divide a
term, so it leaves the F product). Both equal `1 + 58/q + O(1/q²)`: **they
cancel**. Measured confirmation at X = 7.26e16, pins 59..173 vs 59..113: cost
6.99x, probability gain 6.90x — cancellation to 1.3%.

At second order the cost term is `58²/q² = 3364/q²` and the gain term is
`C(58,2)/q² = 1653/q²`, so **cost wins and extra filtered pins are a net loss**.
The measured hits/core-s column confirms it monotonically (59..113 → 59..131 →
59..173 at X = 2.18e13: 1.79e-5 → 9.3e-6 → 2.2e-6).

**Actionable rule, and it reverses today's engine's design order:** pick the term
budget X *first*, then pin exactly the primes that fit free in the wheel
(∏Q ≤ X), and **filter nothing**. Today's engine picks the pin set first and then
inflates the term size to 64·MOD to recover supply — which is the wrong way round.

## 4. Where that gives a real, free win today

At the engine's own term floor X = 64·MOD = 7.26e16, the prime 131 *fits free* in
the wheel (M = 4.95e15 ≤ X): measured cost 1.068 work/out (vs 1.0005), probability
gain 1.555x → **net 1.46x by work count (1.22x by wall-clock, 1.61 vs 1.26
ns/out, where the wheel's extra depth shows up)** — free, by adding 131 to the wheel and shortening
the lift range to keep the term ceiling where it is. 137 does not fit
(M·137 = 6.8e17 > X) and filtering it is a net loss (measured 1.85x cost vs 1.52x
gain). So: **one extra pinned prime, ~1.45x, no new code beyond a longer wheel.**

## 5. Why this is NOT 200x for the project — the part that matters

Re-derived from scratch (`derive.py` reproduces the briefing's selectivity column
exactly: 442, 23035, 54106, 118150, 242745), and combined with the briefing's
(inherited, *not* verified here) P58 model:

- Per unit of work, the best feasible configuration — pins 59..107, terms at the
  natural floor 2.18e13, 1.36 ns/output — beats today's structure (pins 59..113,
  terms at 7.26e16, 1.26 ns/output) by **497x** in expected 58s per core-second.
  The briefing's ~200x is real and if anything understated.
- **It is capped by supply.** There are only `delta_total · X` = **1.84e8**
  admissible `a` below 2.18e13 for that pin set — 0.25 core-seconds of work, with
  expected yield 1.3e-5 of a 58. The region where the 497x applies is exhausted
  almost instantly.
- Sweeping K (d = K·D0, term size 57·K·D0, wheel modulus held at 3.34e11 so the
  cost stays ~1.4 ns/candidate) gives the integrated picture:

| Kmax | term size | candidates | core-hours | E[58s], small-first | E[58s], same work at engine floor |
|---|---|---|---|---|---|
| 1 | 2.18e13 | 1.8e8 | 0.00008 | 1.3e-5 | 1.2e-8 |
| 100 | 2.18e15 | 1.8e10 | 0.0084 | 8.5e-5 | 1.2e-6 |
| 3331 | 7.26e16 | 6.1e11 | 0.28 | 2.0e-4 | 4.0e-5 |
| 10000 | 2.18e17 | 1.8e12 | 0.84 | 2.7e-4 | 9.1e-5 |
| 400000 | 8.71e18 | 7.4e13 | 33.6 | 7.1e-4 | 4.5e-4 |

The ratio collapses from ~1000x at K=1 to **1.6x at 34 core-hours**, because
candidates at term size ≤ T grow like T² while the per-candidate advantage grows
only like a power of log T. Most of any realistic budget is spent at large T
where the two structures are identical.

**So the correct reading of the result:** the enumeration floor is a genuine
inefficiency with a genuine O(1) fix, but it is **not the project's binding
constraint**. It buys (i) the free ~1.45x of §4, (ii) exhaustive coverage of the
small-term region the current engine cannot reach at all, and (iii) the freedom
to choose any term budget independently of the pin set — which is the right thing
to have, but is worth single-digit multiples, not 200x, over a realistic sweep.
I am flagging this as "the experiment answers the question asked, but the
question is not the 200x lever it was believed to be."

Caveats I will not paper over: the P58 model (`2·0.638/(√ln T · F)`, `F = ∏ p/(p+1)`)
is inherited from the briefing and **not verified** in this experiment — every
cost number here is measured, every probability number is modelled. The K-sweep
also assumes K's own prime factors add no new constraints.

## 6. Correctness

- `brute.py`: full scan of every `a < X` and exact set comparison (no missing,
  no extra) for 6 (d, X, Q) cases × 5 split choices = **30 configurations**,
  including d = D0, 7·D0, 1234567·D0 and Q up to 9 primes. All MATCH.
- `enum --verify`: recomputes `j = (a mod q)·d⁻¹ mod q` and asserts
  `1 ≤ j ≤ q-58` for every q, using a **brute-force linear search** for the
  inverse so it shares no table or code path with the search. Verified
  **45,000,000** emitted values (20M at pins 59..107, 20M at 59..113, 5M at
  59..173) — `badver=0` in all 29 runs of `bench.sh`.
- The verifier earned its keep: it caught a real buffer overflow (filter table
  sized 160, primes 167/173 overflowing it) that `brute.py` missed because its
  cases stopped at 137. Fixed (table 512 + explicit guard) and re-measured; the
  numbers above are post-fix.

## 7. Reproduce

```
cd experiments/2026-10-03-enumfloor
cc -O3 -march=native -o enum enum.c
python3 derive.py          # re-derives D0, forced primes, selectivity, P58 table -> derive.out
python3 brute.py           # brute-force correctness, ~9 s          -> brute.out
./bench.sh > bench.out     # all cost measurements, ~60 s
python3 analyze.py         # cost x P58 -> hits/core-s, K sweep     -> analyze.out
```

Files: `enum.c` (enumerator + independent verifier), `brute.py`, `bench.sh`,
`derive.py`, `analyze.py`; raw output in `derive.out`, `brute.out`, `bench.out`,
`analyze.out`.

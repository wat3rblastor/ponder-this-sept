# PROGRESS

Append-only. Newest entry at the bottom.

---

## 2026-10-02 (session 1, 17:00–20:00 CDT) — tooling, calibration, and the real engine

**Deadline set.** User set a hard 12-hour deadline (2026-10-03 05:27 CDT), banned blind brute
force, and required literature research to run concurrently with search. GOAL.md §0 records this.

### Built

- `src/loeschian.py` — bad-prime (`p ≡ 2 mod 3`) Loeschian test, own Miller–Rabin +
  Pollard–Brent, Python ints throughout. Plus an enumeration sieve and an independent
  brute-force `x²+xy+y²` test.
- `src/verify.py` — authoritative verifier. Factors every term from scratch, per-term PASS/FAIL,
  `--maximal` end test. Only this mints a record.
- `tests/` — 19 tests pinning the §2b invariants. Three independent implementations agree on
  every `n ≤ 20000`. `make check` passes.
- `src/c/loeschsearch.c` — flat bitmap sieve to `N`, then for each step `d` the longest run for
  **every** start `a` in one O(N) descending pass. Used for calibration.
- `src/c/loeschclass.c` — class-compressed searcher for steps too large to sieve flat.
- `src/c/apsearch.c` — **the real engine**: Wróblewski-style two-stage search (see below).

### Mathematical findings (the ones that mattered)

1. **Loeschian numbers are never `2 (mod 3)`.** For `t` coprime to 3,
   `t ≡ (−1)^(#bad prime factors with multiplicity) (mod 3)`, and Loeschian forces that count
   even. Hence `3 | d` is forced and every term is `≡ 1 (mod 3)`. *Measured:* with `3 ∤ d` the
   longest run found anywhere is **2**. This single fact raises the per-term density from ~0.14
   to ~0.6 and was the difference between "hopeless" and "tractable".
2. **A factoring-free exact test along a progression.** For `t ≡ 1 (mod 3)`, `t` is Loeschian iff
   every bad prime `p ≤ √t` has even valuation — because `t` has at most one prime factor above
   `√t` and the mod-3 parity then forces it to be good. So an exact Loeschian bitmap for
   `{r + jD}` costs one parity sieve (~1.6·J ops), no factoring. Used in `loeschclass.c`.
3. **The bad-prime cliff is at `p < n`, not `p ≤ n/2`** — but only as a *cost*, not a theorem.
   Measured with `29 ∤ d`: per-term survival drops 0.55 → 0.36 and hits exactly zero at `k = 32`.
   GOAL.md §2-6 originally overstated this as necessary; **corrected** — a verified 25-term AP
   exists with `d = 990`, omitting 17 and 23.

### Calibrated cost model (reproduces measurement within a factor of 2)

`P(all n terms Loeschian) ≈ ρ^n · corr(n)` with `ρ = 2·0.638/(√(ln T)·F_d)`,
`F_d = ∏_{bad p | d} p/(p+1)`, and `corr(n) = ∏_{bad q ∤ d} [(1−n/(q+1)) / e^{−n/(q+1)}]`.
Measured vs predicted: `P(31) = 3.9e−8` vs `5e−9`…`1.2e−8`; `P(36) = 4.8e−10` vs `8.1e−10`.
Note `corr` decays fast in `n`: `0.11` at n=31, `0.02` at n=36, `1.75e−3` at n=58. A constant
correction factor (as one research agent assumed) overstates feasibility by ~3 orders.

### Campaigns

- `experiments/2026-10-02-calibration` — flat sieve scan, ~80 s total.
  **Best: n = 36 at a = 73415383, d = 37418700**, verified by `src/verify.py` and maximal at both
  ends. Also established the geometric decay (ratio ≈ 0.555/term) and the 29-cliff.
- `experiments/2026-10-03-ap58` — the `n = 58` campaign, 8 workers on disjoint `K` windows.
  `d = K·D0`, `D0 = 382160924970`; tier B = {59,71,83,89,101,107,113} (MOD = 1.134e14,
  1.168e9 admissible residues per unit); tier C = 148 bitmask primes to 2000.

### Engine design (`src/c/apsearch.c`), after Wróblewski's AP26 search

Both research agents independently identified this as the only technique in the record-hunting
literature that transfers. Three tiers of bad primes:
- **A** (`p | d`): no term is ever divisible by `p`; only need `p ∤ a`.
- **B** (smallest bad `q > n`): good residues are `a ≡ −k·d (mod q)` for `k = n..q−1`, an
  arithmetic progression in `a`. The CRT-admissible set is enumerated **directly by nested
  additive loops**, so an inadmissible `a` is never materialized. Worth `∏ q/(q−n) ≈ 3e4`.
- **C** (bad primes up to 2000): 64-bit bitmask words, one modulo + AND deciding 64 candidates
  at once, short-circuiting on zero.
Stage 3 confirms with an exact test; `src/verify.py` re-derives independently.

**Validated end-to-end**: the engine independently rediscovered `a = 2235271, d = 129030`
(a 28-term AP found separately), and its own hits pass `src/verify.py`.

Two bugs caught before they cost compute: a Barrett magic that truncated a 128-bit quotient to
64 bits (replaced with plain `%` — on Apple silicon Lemire's fastmod measures only ~3% faster),
and a cofactor test that checked the parity of the *total* bad-factor count rather than each
exponent, which wrongly accepts `p·q` for distinct bad primes. Both are now in `--selftest`.

### Measured throughput and the honest projection

7.36e6 residues/s/core ⇒ **3.65e14 raw candidates/s on 8 cores**. With
`P(58) ≈ 1.1e−21` per raw candidate at the minimal term size, `n = 58` needs ~7.5e20 candidates
= **~24 days on this machine**. A research agent benchmarked a Metal GPU port of this exact
kernel at 1.83 G residues/s (7× the whole CPU), which brings it to ~3.4 days. We have 0.4 days.

**So `n = 58` is roughly 20–50× out of reach on this hardware in this budget.** Expected reach is
n ≈ 54–55; `n = 58` is a few percent. Recorded here so the next session does not re-derive it.

### Ruled out (do not redo)

- **Meet-in-the-middle is impossible here**, and two agents confirmed it appears nowhere in the
  AP-record literature. For fixed `d` every constraint is a congruence on the single unknown `a`,
  so the admissible set is a product of per-prime residue sets and direct CRT enumeration is
  already optimal. Splitting indices `0..28 / 29..57` gives two halves over the *same* parameter
  space — the join is free but generating either half costs the whole search.
- **No published `(a,d)` exists for any long Loeschian AP.** The IBM September 2026 solution page
  is 404 as of today; the blog lists only lengths (57 La Vallee; 55 ×5; 51, 50, 48…). No solver
  published code or numbers; no OEIS sequence for Loeschian APs exists even in draft. So there is
  nothing to extend or seed from.
- **The `q²` relaxation for `q = 59`** (allowing the hit index inside the window when
  `59² | t_k`) doubles the admissible density but requires MOD × 59, which inflates terms ~59×
  and costs `(ρ'/ρ)^58 ≈ 0.08`. **Net loss ~12×.** Skip.
- **The powerful-number trick** (`terms = m·k` for `k = 1..n`, which produces the records for APs
  of powerful numbers) cannot work: it needs `m·k` Loeschian for every `k ≤ n`, impossible for
  `k ≡ 2 (mod 3)`.
- **Extra multiplicity in `d`** (`3²`, `11²`, as seen in Pfoertner's sum-of-two-squares records)
  is pure cost for G3. Those records minimize the *last term*, a different objective; the
  multiplicity is an artifact of that. Possibly relevant to G1 later.

### Next

1. Keep the 8-worker `n = 58` campaign on the minimal-`T` band (`K ≤ 333`), then widen.
2. Metal GPU port of the stage-1/2 kernel (7× measured; skeleton already compiled and
   benchmarked on this machine). **Critical:** keep the index arithmetic 32-bit — 64-bit `%` on
   the Apple GPU is 15× slower and would make the GPU *slower* than the CPU.
3. Update `records.json` + `ANSWER.md` on every improvement; final write-up before the deadline.

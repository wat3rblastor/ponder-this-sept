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

---

## 2026-10-02/03 (session 1 continued, 20:00–00:00 CDT) — the real engine, and G1 opened

### Record improved: n = 36 → **n = 43**

`a = 5413537078288507`, `d = 48916598396160 = 128·D0`, last term `7468034210927227`.
Verified by `src/verify.py` (maximal at both ends) **and** by the new constructive
cross-check. This clears the G2 bonus (`n ≥ 42`). Found in `experiments/2026-10-03-ap58`
at work unit `K = 128`, ~20 minutes on 8 cores.

### Constructive independent verification (`src/crosscheck.py`)

`verify.py` applies the bad-prime criterion — the same fact the search is built on, so using
it alone is nearly circular. `crosscheck.py` instead *constructs* `x, y` with
`x² + xy + y² = t` for every term, via Eisenstein-integer arithmetic: factor `t`, represent
each prime power as a norm (`3 = N(2+ω)`; `p ≡ 1 mod 3` by Cornacchia on `p = u²+3v²`;
`p ≡ 2 mod 3` only ever through `p²`), multiply, and check the identity by direct arithmetic.
A representation is positive proof independent of any bad-prime theory. **PASS on n=43.**

### Throughput: 2.3× from three fixes (3.5 → 5.8 M residues/s/core)

Profiling said 277 ns/residue with ~240 of it in stage 3. Causes and fixes:
- `pollard()` took **one gcd per iteration** (~1 ms per factorization). Brent's variant with
  the gcd batched every 128 steps → the Loeschian test is ~20 µs and stage 3 is now ~8 ns per
  residue, i.e. no longer the bottleneck.
- Tier C was iterated in **ascending prime order**, so after 2 and 5 moved into stage 1 the AND
  chain began with the feeblest killers (11, 17, 23 cut only `1/p`) and almost never
  short-circuited. Now ordered by kill fraction.
- Modulo by tier-C primes folds three 17-bit limbs of `R` with precomputed `2^17`, `2^34`
  residues, so the divide is 32-bit instead of 64-bit.
Also `--isl` was added to pipe values through the searcher's own test: 300 random values in
`[1e16,1e17]` agree exactly with the Python implementation.

### Stage 1 now pins 2 and 5 as well (≈1.8×)

`a ≡ 1 mod 2` and `a ≢ 0 mod 5` belong in the CRT enumeration, not the bitmask stage:
candidates with `a` even or divisible by 5 are now never constructed. Stage-1 components are
generic `(modulus, start, step, count)` tuples. `MOD = 1133661268029390` over
`{3,2,5,59,71,83,89,101,107,113}`, 4.67e9 admissible residues per unit — the design an
optimization agent derived independently as the exact knapsack optimum for `MOD ≤ 1e16`.

A **CRT self-check** was added after a silent-failure scare: it verifies `R0` hits the first
good residue of every component and that each additive step moves only its own component.
(The scare itself was an observation error — piping the searcher through `head -3` SIGPIPE'd it
before it finished the unit. Don't diagnose a searcher through a truncating pipe.)

### Why n ≥ 58 is not reachable here, with numbers

Measured 5.8e6 residues/s/core ⇒ **7.2e14 candidate `(a,d)` pairs/s on 8 cores**.
Calibrated `P(58) ≈ 2.35e-22` per raw candidate at the minimal term size, times a 0.4415
penalty because avoid-only pinning cannot find APs where `q² |` a term. Expected 58-hits over
the whole remaining budget: **~1.6e-3**. Expected best run: **n ≈ 47–48**.
Independent agreement: a research agent reproduced `P(58)` to within 14% by a different route,
and a second put the total cost at 12–30 CPU-days (its figure was optimistic only because it
held the correlation correction constant in `n`; measured, that factor falls 0.11 → 0.02 →
1.75e-3 at n = 31, 36, 58).
The barrier is *existence plus coverage*, not cleverness: 58-term APs first become abundant
around last term `~1e18`, where the space holds `~1e22` candidates, and we can examine `~1e19`.

### Ruled out this session (do not redo)

- **Extending or repairing a found AP by scaling is impossible, probability exactly 0.**
  Define the bad-parity vector `π(n) = {q bad : v_q(n) odd}`; then `n ∈ L ⟺ π(n) = ∅` and
  `π(cn) = π(c) △ π(n)`. So `c·t_j ∈ L` for all `j` iff `π(t_j) = π(c)` for all `j` — scaling
  XORs every term's vector identically. In a primitive AP each bad prime divides at most one
  term, so `m ≥ 2` terms cannot share a nonempty vector; hence `π(c) = ∅`, `c` is Loeschian,
  and length is preserved. Corollary: **a 57-of-58 near miss is never repairable**, so
  collecting near misses is pointless. Verified on the 28-term AP: exhausting all 2⁷
  multipliers built from every bad prime near the window leaves the maximum at 28.
- **Sub/super-progressions**: a step-`d/2` super-progression needs 115 consecutive Loeschian
  terms — strictly harder.
- **Pool/meet-in-the-middle**: no combining operation exists (an AP is fixed by two terms, and
  the only structure-preserving maps are length-preserving Loeschian scalings).
- **The `q²` relaxation at `q = 59`**: admits 2× more candidates but finds only 1.9657× more
  solutions, and costs `log 59` of extra modulus. Net ≈ −1.6×.
- **More or fewer tier-B primes**: adding 131 makes `MOD = 1.49e17`, forcing terms ~130× larger
  (net −22×); dropping to 6 primes is −5×. The 7-prime set is the optimum because `64·MOD` sets
  the minimum term size and density falls as `1/√(log T)`.
- **Dropping 41/47/53 from `D0`** and handling them mod `p²`: −0.018× to −0.47×. Keep them.
- **`3²` or `3³` in `d`**: −0.33×, −0.11×. 3 is unconstrained; extra factors only inflate terms.
- **Ordering `K` by singular series**: worth ×1.00. Requiring `q | d` is arithmetically the same
  object as pinning `q` in `MOD` (gains agree to 0.01%), and that is already collected.
- **No published `(a,d)` exists to extend or seed from.** The IBM September 2026 solution page
  is still 404; the blog publishes only lengths. No solver published numbers or code; OEIS has
  no Loeschian-AP sequence even in draft. (Final standings, for the record: 57 Jackson La
  Vallee; 55 ×5; 51 ×2; 50; 48 ×3; 47 ×4; 46 ×3; 45 ×4; 44.)
- Apple M2 **does** expose a GPU to Metal and a benchmarked port of this kernel runs at 1.83 G
  residues/s, 7× the whole CPU — but only if the index arithmetic stays 32-bit (64-bit `%` on
  that GPU is 15× slower and would make it *lose* to the CPU). Not implemented: 7× buys ~1.4
  terms, and the CPU-side fixes above were cheaper per unit of gain.

### G1 opened (GOAL.md §1 sanctions this once G3 is blocked)

G1 — a 35-term AP minimising the last term — is the puzzle's *actual* main challenge, and
unlike G3 it is **exhaustively solvable**: for `n = 35` the bad primes with `2p ≤ 35` are
`2, 5, 11, 17`, and `p | d` is a theorem for those, so with the forced `3` **every** valid step
is a multiple of `5610`. Scanning `d = 5610·m` for all `m` with `34d ≤ L`, over all `a`, is
therefore a *complete* search for last term `≤ L`, and yields a provable minimum.

`loeschsearch.c` gained `--nmin`, reporting the smallest `a` whose run reaches `nmin` for each
step (the descending scan means the last qualifying `a` seen is the smallest).

- `experiments/2026-10-03-g1`, `L = 2e8`, `m = 1..1049`: **no 35-term AP exists with last term
  ≤ 2·10⁸** — longest runs in the family top out at 27–29. Exhaustive negative result.
- `L = 4e8`, `m = 1..2100` in progress. The productive sub-family is `m ≡ 0 mod 667`
  (`d` divisible by `23·29`), where removing those two cliffs is worth ~130×.

---

## 2026-10-03 (session 1 continued, ~19:00–) — GPU port, record n = 47

### Record improved: n = 43 → **n = 47**

`a = 2646171143023357`, `d = 78342989618850 = 2·3·5²·11·17·23·29·41²·47·53` (`K = 205`),
last term `6249948665490457`. Verified by `src/verify.py` (maximal at both ends) and by the
constructive cross-check. Found within minutes of the GPU engine going live.

### Metal GPU port (`src/c/apsearch_gpu.m`): 4.1e8 residues/s, ~9× the whole CPU

Profiling had put ~164 of the CPU engine's 172 ns per residue in stage 2, which is
embarrassingly parallel over ~5·10⁹ independent residues per work unit. One GPU thread per
prefix of the additive loop nest, each walking the two innermost levels.

- **The trap that decides it:** 64-bit `%` on the Apple GPU is ~15× slower than 32-bit and
  would make the GPU *lose* to the CPU. `R < MOD < 2^51` is therefore split into three 17-bit
  limbs and folded with precomputed `2^17`, `2^34` residues; the fold is `< 2^32` for
  `b2 ≤ 16000`, asserted at startup.
- Only 64-bit add/compare/subtract are needed on device (MSL has no 128-bit type), so the
  per-thread start residues are computed host-side where `mulmod` exists.
- Deepening tier C to 10000 drops survivors to ~480 per unit, so the CPU-side exact test costs
  1.6 s against the GPU's 11.3 s and is no longer a bottleneck.
- Validated by reproducing the n=43 record at its own work unit, K=128.
- `loesch_core.h` now holds the exact Loeschian test, shared by both engines rather than
  duplicated — divergence between two copies of that function is exactly the bug that would
  invalidate a record silently.

**Bug worth remembering:** the Metal buffers were allocated inside the work-unit loop but
released only when `main()` returned (a single outer `@autoreleasepool`), leaking ~100 MB per
unit. The machine went to 9.9 GB of swap and units took 420 s instead of 12. Also: `pkill -f`
with narrow patterns had left 14 stray workers alive earlier, costing 4.9 GB — check
`ps | grep -c` after killing, not just the pattern you think you used.

### Tooling

`tools/promote.py` is now the only way a record is recorded: it refuses any candidate that does
not pass **both** verifiers, appends the superseded record to history, and regenerates
`ANSWER.md` from `records.json` so the answer file cannot drift from the data.

### G1 SOLVED OPTIMALLY (proof, not best effort): last term 311958331

`a = 219830911`, `d = 2709630 = 2·3²·5·7·11·17·23` (`d/5610 = 483`), 35 terms,
last term **311958331** — 4.3× smaller than the 1345651183 recorded earlier (which was merely
the 35-prefix of the n=36 run found during calibration).

**Why it is optimal.** For `n = 35`, every bad prime `p` with `2p ≤ 35` (so 2, 5, 11, 17) must
divide `d`: among 35 consecutive indices at least two terms are divisible by `p`, each needs
`p² |`, and their difference is `p·j·d` with `j < p`, forcing `p | d`. Together with the forced
3, `d` is necessarily a multiple of `lcm(3,2,5,11,17) = 5610`. Any 35-term AP with last term
`≤ 311958331` therefore has `34d ≤ 311958331`, i.e. `m = d/5610 ≤ 1635`.
`experiments/2026-10-03-g1opt` swept **all** `m = 1..1638` against **all** `a`, with every term
`≤ 311958331`, and found **exactly one** 35-term AP — this one. Hence no 35-term Loeschian AP
has a smaller last term. Both verifiers pass and the run is maximal at both ends.

Method note worth keeping: the first attempt at this proof set the bound at the then-record
1345651183, which needed `m ≤ 7062` and a 169 MB bitmap per worker — ~3 h and it drove the
machine into swap. Once the sweep *found* 311958331, re-running with that as the bound cut both
the bitmap (4.3×) and the number of steps (4.3×) — 18× less work, and it finished in 10 minutes.
**Tighten the bound as soon as the search improves it.**

Also note this falsifies, in the useful direction, the calibrated prediction that the n=35
optimum would be near `1.3e9`: the truth is `3.1e8`. Model estimates of *where* a record lies
are worth far less than an exhaustive sweep when the space is small enough to exhaust.

---

## 2026-10-02 ~22:00 CDT — paused by user; clean resumable state

All search processes stopped at the user's request. State at pause:

- **G3 record: n = 47** (`a = 2646171143023357`, `d = 78342989618850`), verified both ways,
  maximal at both ends. G1 is **proved optimal** (35 terms, last term 311958331).
- Shift-0 exhaustion of `K = 1..3365` (d = K·D0, a < 64·MOD ≈ 7.3e16, tier B
  {3,2,5,59,71,83,89,101,107,113}): **complete through K = 784** (GPU logs
  `experiments/2026-10-03-gpu/g2/g4/g5.log`), plus scattered earlier CPU coverage of
  K ≤ 672 at the pre-fix engine (see 2026-10-03-ap58*). CPU band K = 3366..9365 partially
  covered (`2026-10-03-cpu2/c*.jsonl`, interrupted).
- **Resume with:**
  `./build/apsearch_gpu --nterms 58 --kmin 785 --kmax 3365 --shifts 1 --modcap 2000000000000000 --b2 10000 --report 48 --out experiments/2026-10-03-gpu/g6.jsonl`
  plus CPU workers on disjoint `--kmin/--kmax` bands, plus `tools/autopromote.sh 60 &`.
- Ready on the shelf: the 19-bit-limb GPU build (`build/apsearch_gpu2`, validated) for
  deep-pinned tier B including 131 (MOD up to 2^57) — the right tool only for a multi-day
  campaign at term sizes ~1e18, not for short runs.
- Expected-value note for the next session: exhausting the remaining shift-0 region is
  [0.04, 0.8] expected 57-finds depending on which calibration branch is right (n=47 anchor
  vs fixed-T decay model); after that, shift planes 1-2 repeat the same economics on fresh
  space. User directive on record: keep going until n >= 57 (GOAL.md §0), M2 on-chip GPU
  allowed, no rented hardware.

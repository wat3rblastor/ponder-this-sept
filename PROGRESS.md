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

---

## 2026-10-03 03:00–03:55 UTC — NVIDIA GB10: CUDA engine, record n = 55

Environment changed: rented Vast.ai box, NVIDIA GB10 (48 SMs, cc 12.1, CUDA 13.2), 20 cores.
User directives this session: target is **n >= 58**, wanted **by 2026-10-03 15:50 UTC**.

### Record: n = 47 -> **n = 55**

`a = 11687581876345393`, `d = 202927451159070` (`K = 531 = 9·59`, shift 8), verified by
`src/verify.py` (maximal both ends) and `src/crosscheck.py`, promoted by `tools/autopromote.sh`.

### Engine (`src/c/apsearch_cuda.cu`, `make cuda`)

Throughput ladder on one full unit (4.67e9 residues), every step validated against the CPU
engine's hit sets or by set containment:

| kernel | s/unit | note |
|---|---|---|
| Metal M2 (previous) | 11–12 | |
| CUDA 1 chain/thread | 1.31 | straight port, magic-multiply modulo |
| 4 chains/thread | 0.49 | kernel was latency-bound on dependent L1 gathers |
| strided groups, 7 residues × 4 primes per exit test | 0.25 | below |

**Strided groups** (the big one): walk the innermost loop *without* reducing mod MOD, so
`A mod r` advances by a fixed stride; store each tier-C table pre-rotated by that stride so 7
consecutive residues read 7 adjacent words (one fold + ~one cache line per prime per group).
An unreduced `A = R + w·MOD` is the same class with its 64-candidate window moved up by `w`,
so a unit `(K, shift)` now covers, per class, the multiples `[64·shift + w, 64·shift + 64 + w)`.

Tried and rejected (kept behind `--kernel`): per-lane flat chains (1, N chains, and flat
strided) — 2× slower, the constant cache beats per-lane row loads; tables in shared memory —
slower (occupancy). Pipeline bug worth remembering: a plain `cudaMemcpy` uses the legacy
default stream and silently waits for the other stream's kernel (cost 1/3 of the GPU).

### Calibration (full-window histograms, 205 units, `c1.jsonl` hist lines)

Runs ≥ n per unit fall by ~0.737/term for n = 25..47 after removing the position factor;
`E58 ≈ W·ρ^58 ≈ 1e-5` per generic shift-0 unit (W ≈ 546 survivor windows). That is ~3× more
optimistic than the previous session's pessimistic branch.

### K is not uniform: value-ordered planning (`tools/plan_units.py`)

Corrects "ordering K by singular series is worth ×1.00" (2026-10-02). For a bad prime
`q > 58` with `q | K`: no term is ever divisible by `q`, so the tier-C filter passes
`(q-1)/q` of candidates instead of `(q-58)/q` — ×1.79 for q = 131 at the same GPU cost — and
if `q` was a pinned prime its slot goes to the next one (131, 137, …). Measured: K = 131 gives
836 survivors/unit vs ~450; K = 71·83 gives 3× the survivors per residue. The planner scores
each `(K, shift)` by yield per residue, `∏` of those factors times `(ln T)^-29`, and the engine
consumes the list (`--units`). With the modulus cap lifted to 2e16 (the strided kernel only
needs the unreduced walk to fit 64 bits) the model gives E58 ≈ 0.69 per 12 GPU-hours versus
0.24 for plain shift planes on K ≤ 3365. The n = 55 came from a 59 | K unit minutes later.

### Coverage

- Reduced-kernel planes (multiples `[0,64)`): shift 0, K ≤ 1949 (`2026-10-03-gpu`, `c1.jsonl`).
- Strided planes: `c2.jsonl` (shift 0, K 1950..~2700), `s*.jsonl` (queue, brief),
  `p1.jsonl`, `p2.jsonl` (plan order; each unit is a line, `--resume` skips them and
  `plan_units.py` excludes them when replanning).

### Running

`tools/run_plan.sh experiments/2026-10-03-cuda/plan2.txt experiments/2026-10-03-cuda/p2.jsonl 58`
plus `tools/autopromote.sh 60`. Stops itself when records.json reaches n >= 58.

### Research verdicts and handoff (04:45 UTC)

- Mathematical research agent: nothing worth >= 1.3x. Checked d structure and pin sets (MODCAP
  2e16 is the optimum: E = 0.25/0.45/0.50/0.64/0.59 at 1.2e13/1.2e14/1.2e15/2e16/2e17), the 59^2
  family (~1.03x, not worth a kernel mode), raising the 0.737 per-term probability (impossible:
  it is the Poisson tail of bad primes > 1e4), literature (Choudhry gives 9-11 terms only; IBM
  page lists lengths only), longer windows and symmetry (re-enumerate the same APs). A slightly
  better planner score (per-term ln ln average instead of the last term) is worth ~1.035x.
- Stage 3 moved to a worker thread behind a queue (hit sets identical on the validation units).
- Multi-GPU tooling: `tools/build_here.sh`, `tools/multi_gpu.sh`, `tools/deploy_remote.sh`,
  `RENTED_GPU.md`, engine `--slice i n`. GOAL.md opens with a START HERE block.
- Local search stopped 2026-10-03 ~05:00 UTC at the user's request; all finished units are committed. The search continues on a rented multi-GPU
  machine from this repo alone. Units finished here after the final push may be repeated there.

### GPU research agent: warp-compaction kernel measured, not adopted (04:55 UTC)

Candidate `--kernel 30` (fixed unrolled prefix of T0 rows, ballot/shuffle compaction of surviving
groups, optional dual-copy tables with 128-bit loads) is kept, unmerged, in
`experiments/2026-10-03-cuda/kernel30-candidate/`. It is **correct** (hit sets identical to the
reference on the validation units in three configurations) but measured only **1.04–1.08×**:
0.23–0.24 s vs 0.25 s on a standard unit, 4.97 s vs 5.17 s on K=6319 — not the modelled
1.3–1.5×. The vector-load-only variant was slower (0.30 s). Not adopted: the gain does not justify
swapping the production kernel mid-campaign. The agent's CPU model says a group of 7 needs ~25
rows while its warp runs ~61, and infers the kernel is L1-bandwidth-bound; a two-stage
compaction is the remaining idea if someone wants to push further.
Large units (131/137 pinned, e.g. K=6319) run at 1.29e10 residues/s solo against 1.87e10 for
standard units: their early rows are weaker, so chains run longer.
Agent's estimate for other cards, from SM count × clock: RTX 5090 ≈ 3.5–4× the GB10,
RTX 4090 ≈ 2.5–3×. Unmeasured; `tools/build_here.sh` prints the real figure.

---

## 2026-10-03 04:20–06:50 UTC — 8x RTX PRO 6000 (Vast.ai): 2x throughput, structure hunt, handoff to Azure

Machine: 8x NVIDIA RTX PRO 6000 Blackwell (sm_120, 188 SMs), 208 logical cores, 1 TB RAM, CUDA 12.8.
Record unchanged at **n = 55**; a second, primitive 55 was found
(`a = 296246969176050787`, `d = 11950172123811900`, K = 31270) and is listed in ANSWER.md §5.
Coverage on this box: ~4.4e15 residues (`experiments/remote/*.jsonl`, tags r1..r7), no run >= 56.

### Throughput: 3.67e11 -> 7.45e11 residues/s (all committed)

| change | aggregate res/s | note |
|---|---|---|
| launch, 1 engine/GPU, strided kernel 20 | 3.67e11 | GPUs idle during serial per-K host setup |
| 3 engines/GPU (time-sliced) | ~4.2e11 | device saturated |
| kernel 31 (`--kernel 31 --t0 24`, warp compaction of words) | 5.11e11 | 42/42 old units identical |
| host prep pool, sparse table build, pinned buffers (`--prep 8`), 2 engines/GPU | 5.96e11 | 60/60 identical |
| `--report 55` (stage 3 stops when 55 is impossible) | 6.11e11 | runs < 55 no longer logged |
| CUDA MPS (`MPS=1` in multi_gpu.sh) | ~6.2-6.6e11 | +8.5% per GPU by logs |
| per-word load guards in kernel 31 | 7.45e11 | 60/60 identical, +14-20% |

Lessons: (1) do not cap OpenMP threads: stage 3 runs in a std::thread that ignores `--threads`
and needs ~25 cores per GPU; `OMP_NUM_THREADS=8` was 3.8x slower. (2) the `gpu=` log field is
the inter-launch interval, not kernel time. (3) host pipeline now loses <2%; the cards are
thermally capped (85 C, ~2200 of 2430 MHz, ~460 W): hardware. (4) compressed/shrunk tail tables
do nothing: all GPU time is in the first 24 rows. (5) never `pkill -f <pattern>` from a shell
whose own command line contains the pattern. Parked, unvalidated: queued-launch + blocking-sync
host patch (`experiments/2026-10-03-structure/botl/queued_launch.patch`, expected 0-3%);
Montgomery stage-3 arithmetic (2.2x less CPU per exact test, not needed while GPU-bound).

### Planner (tools/plan_units.py, corrected)

Rescaled copies: units with a good prime p | K re-find (p a, p d) images of earlier progressions;
every "new" n=55 on this box except K=31270 is the record times 13/19/31/37/43. The planner now
down-weights good primes in K (0.68/0.80/0.84/0.89 for 7/13/19/31), accounts for the strided
window offset and averages term size per term: ~1.15x. Exact removal of copies is worth only
1.3% of compute; pinning them out is a net loss. Calibration: per-term pass rate after the 1e4
sieve rho(T) = 0.72 (ln T / 39.14)^-0.485 (Monte Carlo; validated by exact counts up to 43-digit
terms); observed/expected long runs 0.82 +- 0.09 (2 sigma, hits are clustered). Rate: ~0.32-0.40
expected 58s per hour at 6e11 res/s. MODCAP 2e16 remains optimal; the 64-bit cap costs nothing.

### Hints from the 57-term record holder, and what was measured

Relayed by the user: his a, d are "orders of magnitude" larger; "you can do it by hand if you
think about it the right way"; "exploit the prime structure"; he used two of {large d, square
option, automatic terms}; "every 57 contains three 55s". Results (sources and logs in
`experiments/2026-10-03-structure/`):

- **Large d alone: no gain (measured).** d = 3 * prod(bad primes <= y), y = 53..200: per-term
  pass rate stays 0.44-0.50 unsieved and observed run counts match the model within 2% up to
  43-digit terms. Each extra prime costs ~7-8x per word beyond our pin set.
- **Automatic terms: real, but too thin.** General rule: if t_{k0} is Loeschian and
  3*c*d*t_{k0} is a perfect square (c >= 1 integer), then t_{k0 +- c j^2} is Loeschian for all j.
  Shapes: a = x^2, d = 3 m^2 (positions 0,1,4,...,49 free: 8 of 58; a bad prime can hit index k
  only if k is a non-residue); a second square in the progression frees 8 more (conic,
  m = 2AB, x = |6B^2 - A^2|); a = 3x^2 + m^2, d = 24 m^2 gives 13 free (generalized pentagonal
  k), 23 with one extra condition, 32/40 with two/three (elliptic curve / finitely many).
  Measured: 20-49x more runs of 10-20 terms than a generic AP of the same size, free positions
  never fail. But the families are 1-dimensional (~T/1e13 members below T vs T^2/4e13): after an
  equal-depth sieve the single-square family is 0.1-0.8x of baseline per walked value and has
  ~1e7x too few candidates; a 58 needs terms near 1e40. The elliptic-curve (3+ squares) version
  was being tested by an agent at handoff (`exp_ec/`): untested question is whether small points
  exist once 5, 11, 17, 23, 29 must divide m, or whether those primes can be confined to
  automatic positions.
- **Square option: modest gain, the one actionable lead.** Per candidate at equal pool size for
  n = 58: square classes for 59 (a ≡ -j d mod 59^2, j = 0..57) 1.5x; 53 left out of d with
  53^2 on its single hit (j = 5..52) 1.10x; 47 out 1.01x; 41 out 0.92x; 29 out (n = 57 only)
  0.30x. More important: 58-window yield falls ~6x per decade of term size (X^-0.79), and these
  families are fresh pools of small-term candidates (union ~12.5x the baseline density at
  n = 58, estimated ~5x fewer candidates to a first hit). An agent was adding them to the engine
  as per-unit "modes" with a brute-force oracle at handoff (`/workspace/sqopt` on the Vast box;
  not merged, not in this repo unless a later commit says so). If it is not here: the engine's
  stage-1 component (modulus q^2, start, step = -d mod q^2, count) already has the right shape;
  the unit key must include the mode; tier C must not also avoid q.
- **No search-free construction exists by any route tried**: a linear function of k cannot be
  identically a norm; scaling cannot repair a term; CRT-forced factors cover at most 2/58 of
  each term's digits (lattice check); polynomial identities reduce to an existing AP.
  Exhaustive small cases (n = 8..28) show optimal APs look like random admissible ones.
- Web: IBM has not published the September solution (2026-10-03); no public (a, d) for any 55/57.
  Choudhry (Integers 25 #A14) gives parametric 9-term APs and states no general method is known.

### Handoff

Vast.ai box is being dropped; the search continues on Azure from this repo (GOAL.md START HERE).
Units finished on the Vast box after the last coverage commit may be repeated on Azure.

Addendum (small-case agent, `exp_small/`): by counting, the all-in-d family (41, 47, 53 | d)
holds only (41/65)(47/83)(53/101) ~ 19% of the 58-term progressions at a given term size; the
square-option sub-families for 41/47/53 hold the other ~81% (~5x). Per-term pass rate (unsieved
beyond Q, Monte Carlo): 0.571 at 1e12, 0.468 at 1e18, 0.414 at 1e24 for Q = 53. Smaller terms
matter far more than anything else, which is why the square-option pools are the next step.

Addendum (elliptic-curve agent, `exp_ec/`, first numbers 06:55 UTC): the 3-square version works
mechanically but looks dead for 57/58. Shape a = 3X^2 + m^2-type with three bases gives at most
31 automatic positions of 57 (32 of 58); four bases 39/40 (greedy). Of 174 index-difference pairs
with >= 30 coverage, 104 curves have rank 0, 66 rank 1, 4 rank 2. The small bad primes cannot be
confined to automatic positions, so 5, 11, 17, 23, 29 must divide m; the smallest point with
623645 | m found has a 34-digit m (72-digit terms). Exact test of that point: automatic positions
never fail, but only 6 of 35 resolved non-automatic terms are Loeschian (~17%); longest run 6.

---

## 2026-10-03 — GOAL2 cycle 1 (orchestrator)

Method switched to `GOAL2.md`: the premise is that `n >= 58` is reachable with this laptop, the
orchestrator plans and judges, agents implement, and throughput is never a plan. Four independent
hypotheses dispatched in parallel, each with a threshold fixed in advance.

| # | hypothesis | the number that decides it | dir |
|---|---|---|---|
| H1 | The term-size floor `T >= 64*MOD` is an artefact of our enumerator, not the problem. If CRT-admissible `a` **small relative to the modulus** can be enumerated at ~O(1) each, the pin table says up to **200x** at the same term size. | cost per admissible `a < X` emitted, for a pin set with `1/delta >= 5e4`. WIN if < ~100 ops each; LOSS if > ~1e4 | `2026-10-03-enumfloor` |
| H2 | The automatic-terms covering problem has never been posed. Free indices are `k0 + c*j^2`; a *union* over several `(c,k0)` could cover far more than the 8 and 13 seen so far. | smallest term size `T` admitting >= 30 of 58 automatic indices. WIN if `T <= 1e20`, LOSS if `T > 1e30`; the coverage-vs-min-T frontier is the artefact either way | `2026-10-03-covering` |
| H3 | Our own long progressions are the only ground truth about what solutions look like, and have never been dissected against a matched control. | strongest deviation from matched control. WIN at >= 5x, LOSS if everything is within ~2x | `2026-10-03-forensics` |
| H4 | Every ruled-out idea rests on a premise; one of them is wrong or escapable. | names the assumption and exhibits the construction/counterexample, or confirms them | `2026-10-03-assumptions` |

Each brief carries the full mathematics restated, a falsifiable threshold, the laptop-only budget,
and the rule that a clean negative goes in the first line.

### Cycle 2 hypotheses, drafted while cycle 1 runs (queue must never be empty)

- **H5 — interleaving / genuine MITM.** A 58-term AP with step `d` is *exactly* a pair of 29-term
  APs with step `2d` whose starts differ by `d`. That is two coordinates, not one, so the old
  "MITM does not apply" argument (which assumed `d` fixed and `a` the only unknown) does not cover
  it. Counting check to do first: scanning `a` for a 29-run costs the same as for a 58-run
  (expected ~2 term tests either way), so the naive version is neutral — the question is whether
  the *sieve* costs differ, since a 29-window needs much weaker pinning (`q-29` good residues
  instead of `q-58`). If the sieve is cheaper per 29-run than per 58-run by more than the join
  costs, this wins.
- **H6 — cross the two failed ideas.** The square option (bad `p < 58` left out of `d`, paying
  `p^2` on its single hit) failed on cost; automatic terms failed on size. But a hit that lands on
  an *automatic* index is free regardless of what divides it. Do the forced small primes'
  hit-classes intersect the automatic set for any shape? Previously checked only for the
  pentagonal shape at `p = 5`.
- **H7 — the kernel formulation.** `t` Loeschian iff its squarefree kernel avoids all bad primes.
  Pose the problem as: 58 numbers in AP whose squarefree kernels avoid a fixed prime set. Fixing
  kernels for a subset of indices gives an overdetermined system on `(a,d)` — how many indices can
  be fixed before it has no solutions, and what does the boundary case look like?

### GOAL2 cycle 1 — verdicts

| # | hypothesis | verdict | what it changed |
|---|---|---|---|
| H1 | enumerator term-size floor | **WIN on threshold, small in practice** | 1.02 work units per admissible `a` emitted (threshold was <100), verified against brute force on 30 configurations, 45M values independently re-checked. But supply-capped: the 497x per-work edge collapses to 1.6x once enough candidates are needed. Free 1.46x available today by pinning 131 in the wheel. |
| H2 | automatic-terms covering | **LOSS, with the best structure of the cycle** | max automatic coverage is **17**, at any size |
| H3 | forensics on known solutions | **clean LOSS, well powered** | long progressions are indistinguishable from random admissible ones (16 features, all within 1.06x of a matched control, 3486 real terms vs 7320) |
| H4 | attack the register | (running) | |

**H1's structural result, which reverses the engine's design order.** Pinning an extra prime `q`
costs `q/(q-58)` in enumeration and gains `((q+1)/q)^58` in pass probability. Both are
`1 + 58/q + O(1/q²)` — they cancel (measured 6.99x cost vs 6.90x gain, 1.3% apart); at second
order cost wins, so filtered pins are a **net loss**. The pin set is therefore nearly irrelevant
and **the only lever is term size**. Rule: choose the term budget first, pin exactly what fits
free in the wheel (`prod Q <= X`), filter nothing.

**H2's classification.** A term is free only if `a+kd = L1² + 3L2²` identically, so `a` and `d`
are binary quadratic forms and the whole automatic set is governed by ONE quadratic
`P(k) = -disc(a+kd)/12`, with `A = {k in [0,57] : P(k) a perfect square}`. Regimes:
`α = 0` (P linear) is *all* the known theory and caps at 13; `α > 0` (Pell conic) reaches 15;
`α < 0` (ellipse) reaches **17** — both new and never exploited. Coverage frontier (min T):
10-12 → 7.1e7, 13-14 → 2.9e9, 15-16 → 3.8e9, 17 → 7.3e12, 18+ → nothing at any size.

**H2's other result: bad primes cannot be parked on automatic indices.** A bad prime divides
`X²+3Y²` only if it divides both `X` and `Y`, so parking `p`'s hit class needs `p | y_k` at every
hit index, which is unsatisfiable. Hence `2*5*11*17*23*29 = 1247290 | d` **always** — this is a
strengthening of the forced-prime theorem, and it is what keeps the square-shape families large.

**H3's methodological catch**, worth keeping: against the *unconditional* null the hit indices
`k0 = -a/d mod q` look like a 60x deficit (z = -9.4). The correct conditional null is
`B/(qA+B)`, under which obs/exp = 0.90x. The 60x was an artefact of the wrong null.

### Open items carried forward (flagged by the agents themselves)

- H2's "18+ impossible" is exhaustive only inside finite boxes, and the new Pell (15) and
  elliptic (17) mechanisms were screened **at the discriminant level only — no integral pencil
  was constructed**. A 17-hit family starting at `T >= 7.3e12` rather than `5.3e14` has never
  been costed. This is the strongest surviving lead from cycle 1.
- H1's wheel enumerator is implemented and verified but not merged into the production engine.

### Cycle 2, dispatched

- **H5 — is the counting heuristic right about *where* 58s first exist?** Everything above is
  conditioned on it, and it was **4x pessimistic** at n=35 (predicted 1.3e9, truth 3.12e8).
  Measure M(n) exhaustively for n = 20..40+, fit, extrapolate to 58 with honest error bars.
  WIN if M(58) <= ~1e15 (the search region is far denser than modelled and the project re-aims);
  LOSS if >= 1e17. `experiments/2026-10-03-mcurve`.

### GOAL2 cycle 2 — verdicts (all losses, all corrective)

| # | hypothesis | verdict |
|---|---|---|
| H4 | attack the register | **partial win**: the wrong premise is the `41,47,53 \| d` convention — but see H7 |
| H5 | is the counting heuristic right about where 58s first exist? | **LOSS** — `M(58) ≈ 3e17` (range 2e16–2e18), heuristic confirmed |
| H6 | construct the Pell/elliptic pencil (15/17 free terms) | **LOSS** — best family needs 3.6e18 members at `T ~ 5e24` |
| H7 | implement + measure the square option | **LOSS** — validates perfectly, yields 0.40x |
| — | has anyone published a long progression? | **LOSS** — nothing anywhere |

**H4 proved register item 1 shut.** If all terms of a non-primitive AP are Loeschian with
`g = gcd(a,d)`, then `pi(g)` is empty (a prime in it would divide every term, hence both `a'` and
`d'`), so `g` is Loeschian and the reduced AP is itself a solution. Every non-primitive solution is
(Loeschian constant) x (primitive solution). The same XOR argument kills the variable-multiplier
and scaled-set variants. Do not reopen.

**H7 corrected H4 and me, and this is the important one.** The 5.33x is the ratio of the *union of
all 8 modes* to all-in-d; the pure square family is **0.406x** of all-in-d. The other 81% is spread
over seven modes and costs work proportional to its size — extra space, not a speedup. And the
term-size claim double-counted: dropping 41,47,53 lowers `d` by 1e5 but the `q^2` congruences raise
the smallest admissible `a` by 2.5e5, so **all eight modes have the same existence threshold,
2.0–3.1e12**. Per *survivor* the modes are equivalent; the whole difference is survivors/second,
because a `q^2` component costs `q^2` of modulus and buys no enumeration efficiency — H1's
cancellation law, applied correctly.
Validation was clean: the patched engine finds the proved-optimal 35-term progression in 0.5 s
(the all-in-d engine cannot — control best run 31), and produced a new verified 37-term
progression in the square family, `a = 212569971883993, d = 14421166980` (verified, maximal).

**H5's by-products are worth more than its verdict.** New proved exhaustive points:
M(22)=157831, M(25)=160801, M(26)=1673461, M(28)=5719081, M(29..31)=1.853e7/1.918e7/1.982e7,
M(34)=309248701, M(35)=311958331, M(36)>3.30e8, plus n=14..21. A base=1 sweep assuming **no**
forcing reproduces M(14..24) exactly, which empirically validates the forced-base theorem. It also
found three errors in the analytic model this project quotes: a missing triangular 1/2; a missing
`prod_{p|d}(1-1/p)` worth 3.2x at n=58; and `corr_q = 1-n/(q+1)` is valid only for `q >= n` — at
n=58 it is applied to q=41,47,53 where it is **negative**, the correct factor being
`(1-(n-q)/q)/(q+1)`.

**A structural fact worth recording**: n=58 sits exactly on a cliff. `2*29 = 58`, so a 58-index
window covers every residue mod 29 exactly twice, forcing `29 | d`. At n=57 it is optional. Hence
**no 57-term progression with `29 ∤ d` can ever extend to 58** — the cheap 57s are dead ends, and
58 is genuinely harder than 57 by more than one term's worth.

### Cycle 3, dispatched — audit our own numbers

Since H5 found sign errors in the model at exactly the primes that dominate n=58, the headline
cost figures may be wrong too, and they are what say the laptop is ~1e8 short.

- **H8 — measure the laptop cost empirically** rather than from the model: run the engine, build
  the run-length histogram, fit the decay, extrapolate to 58, and cross-check against the one hard
  external point (8 GPUs, 7.45e11 res/s, ~0.4/hour). WIN if < ~1000 laptop-hours.
  `experiments/2026-10-03-laptopcost`.
- **H9 — audit the planner's ranking.** Yield per unit falls ~16x across the plan, so ordering is
  worth ~2.5x; if the score shares the model's errors we are searching in the wrong order. Measure
  old vs corrected ranking against actual engine counts. WIN at >= 2x.
  `experiments/2026-10-03-planaudit`.

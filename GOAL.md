# GOAL — IBM Ponder This, September 2026: Loeschian Arithmetic Progressions

Source: https://research.ibm.com/blog/ponder-this-september-2026
(puzzle suggested by Hugo Pfoertner; author Gadi Aleksandrowicz; posted 01 Sep 2026)

This file is the durable contract for this project. Any agent (me, a fresh session, or a
subagent) should be able to read **only this file** plus `records.json` and immediately know
what to do next. Keep it current; it outranks memory and chat history.

---

## START HERE if you are a fresh session on a new GPU machine (rewritten 2026-10-03 06:50 UTC)

State: best verified is **n = 55** (two distinct progressions, `records.json` / `ANSWER.md` §5);
target is a verified **n >= 58** (a 57 is a waypoint: promote it, keep going). The search moved
Mac -> GB10 -> 8x RTX PRO 6000 (Vast.ai) -> **Azure next**. This repo is the complete state:
every finished work unit is a line in `experiments/*/*.jsonl` (committed), and the planner skips
them. **Do not redesign or re-derive anything: build, plan, launch, watch, verify.** Read the last
two sections of `PROGRESS.md`, then run exactly:

```
make check                                   # python unit tests + CPU sieve build
tools/build_here.sh                          # builds for this GPU; must print best=47 for K=205
MODCAP=2e16 python3 tools/plan_units.py --budget-res 2e16 --kmax 600000 --smax 3000 \
    --out experiments/remote_plan.txt        # ranked (K, shift) units; finished ones are excluded (~3 min)
tools/multi_gpu.sh experiments/remote_plan.txt az1   # 2 engines per GPU under CUDA MPS, disjoint slices
setsid nohup tools/autopromote.sh 60 > /dev/null 2>&1 < /dev/null &
```

- `tools/multi_gpu.sh` defaults are the validated production settings: `PPG=2` engines per GPU,
  CUDA MPS on (`MPS=0` to disable if the daemon cannot start), flags
  `--kernel 31 --t0 24 --prep 8 --report 55`. Do **not** cap OpenMP threads (stage 3 needs ~25
  cores per GPU; a cap of 8 was 3.8x slower). Use a NEW tag (`az1`) on a new machine.
- The kernel was tuned on RTX PRO 6000 Blackwell (sm_120): 7.45e11 residues/s on 8 GPUs
  (~9e10 per GPU). On a different GPU model `build_here.sh` prints the real figure; if kernel 31
  misbehaves there, `EXTRA="--kernel 20 --prep 8 --report 55"` is the older strided kernel.
- Progress: `grep -h '\*\*\* n=' experiments/remote/*.log | sort -t= -k2 -n | tail`. With
  `--report 55` only runs >= 55 are logged. Most n=55 lines are the record times 13, 19, 31, 37,
  43 (rescaled copies); `tools/list55.py` reduces to primitive form and lists distinct ones.
- `tools/autopromote.sh` verifies (both checkers), records, regenerates `ANSWER.md`, commits and
  **pushes** any new best, keeps the 55+ list current, and commits+pushes the coverage logs every
  ~10 minutes. **Push after every commit** (user instruction, supersedes §7 "never push").
- The engines stop themselves when any log shows a run >= 58. Then: `tools/promote.py <a> <d> <n>`
  (must pass `src/verify.py` and `src/crosscheck.py`), confirm `records.json` and `ANSWER.md`,
  append a PROGRESS.md entry, commit, push, report to the user. That is the finish line (§8).
- If an engine dies, rerun the same `tools/multi_gpu.sh` line: finished units are skipped.
- Calibrated expectation (PROGRESS.md 2026-10-03 Vast.ai entry): ~0.4 expected 58s per hour at
  7e11 residues/s, ~0.7 per hour for "57 or better"; yield per residue falls slowly down the plan.
- Open leads and everything ruled out on 2026-10-03 (large d, automatic terms, square option,
  elliptic-curve multi-square family, hints from the 57-term record holder) are in the last
  PROGRESS.md section; sources are in `experiments/2026-10-03-structure/`.

---

## 0. Hard deadline and method constraints (added 2026-10-02 17:27 CDT)

- **DEADLINE (user, 2026-10-03 ~03:50 UTC, verbatim): "I need it 58 within 12 hours."** i.e. by
  **2026-10-03 15:50 UTC**. Measured projection on the single GB10: E[58] ~ 0.44 in 12 h
  (~35%), E[57] ~ 1.2. Throughput is the only lever left; units are independent, so extra GPUs
  scale it linearly.
- **AMENDED 2026-10-03 ~03:38 UTC (user, verbatim): "It'll be great if you find 57, but what I
  really need is you to find 58."** The stop condition is a verified **n >= 58**. A 57 is a
  waypoint: promote it, keep searching. Hardware is now a rented NVIDIA GB10 box (CUDA engine,
  `experiments/2026-10-03-cuda`), which supersedes the "no rented hardware" note below.

- **AMENDED 2026-10-02 ~21:45 CDT: the user's directive is now "keep going until you find 57".**
  The 12-hour deadline below is superseded as a stop condition; it remains the planning
  horizon for the current campaign phase. The stop condition is a verified n >= 57
  (n >= 58 remains the aspiration if 57 falls).
- Original deadline: 2026-10-03 05:27 CDT (12 hours from 2026-10-02 17:27 CDT). The `n >= 58`
  result must be verified, recorded in `records.json`, and written up in `ANSWER.md` before
  then. Budget backwards from the deadline: reserve the last 45 minutes for verification,
  independent cross-check, `ANSWER.md`, and commits. If the deadline is close and no `n >= 58`
  exists, record the best verified `n` reached and say plainly that the bar was not cleared —
  do not keep running past the deadline.
- **No blind brute force.** Undirected enumeration of `(a, d)` is forbidden as a strategy.
  Every campaign must be justified in writing by a structural argument that prunes the space
  first (see §2 and §2b) — which bad primes are *forced* into `d`, which residues of `a` are
  admissible, and why the remaining space is small enough to be worth the CPU. A search that
  cannot state its pruning argument does not get launched. Searching a space that theory has
  already cut down to size is not brute force; sweeping `(a, d)` and hoping is.
- **Mathematical reconnaissance runs in parallel with search, not before it.** Reaching
  `n = 58` on this budget plausibly needs a construction we do not currently have, so
  literature/web research for better ideas is a first-class, continuously-running activity.
  Keep research agents and search agents in flight at the same time (see §6): research agents
  hunt for constructions (analogous records for sums of two squares, admissible-tuple /
  prime-AP techniques, anything by Hugo Pfoertner, Green–Tao-style explicit constructions);
  search agents grind the disjoint families the current theory says are best. Fold any idea
  that survives scrutiny straight into the live campaign.

---

## 1. The problem (verbatim essentials)

A **Loeschian number** is an integer of the form `x² + x·y + y²` with `x, y ∈ ℤ`.
First few: `0, 1, 3, 4, 7, 9, 12, 13, 16, 19, 21, …` (OEIS A003136).

An arithmetic progression (AP) of Loeschian numbers is described by a start `a` and step `d`;
the `n`-term AP is `a, a+d, a+2d, …, a+(n-1)d`, all of which must be Loeschian.
Example: `1, 7, 13, 19` — `n=4`, `a=1`, `d=6`.

### Objectives

- **G3 (MANDATORY — the point of this project):** find a verified AP of Loeschian numbers with
  **n ≥ 58 terms**, beating the published record of 57. This is not optional and not "best
  effort": the project is not done until an `n ≥ 58` AP is verified and recorded. Everything
  else is secondary and may be deferred.
- **G2 (bonus `*`, waypoint):** an AP with **≥ 42 terms**. Last term need not be small. Useful
  only as a pipeline smoke test on the way to G3.
- **G1 (secondary):** a **35-term** AP **minimizing the last term** `a + 34d`. Report `(a, d)`.
  Work on this only once G3 is secured, or when a G3 search is blocked and waiting on something.

### Known bar (from the published solvers list — the contest is closed, these are the targets)

Best published lengths: **n = 57** (Jackson La Vallee), then 55, 51, 50, 48, 47, 46, 45, 44.
So **n ≥ 58 is the mandatory bar**. G1's optimum is *not* published on the page — treat our
best verified `a+34d` as the record and drive it down as a side quest.

`a = 0` is a legal Loeschian number, but note `0, d, 2d, …` with `d` Loeschian requires every
`k·d` Loeschian — do not accidentally "win" G1 with a degenerate reading. Progressions must be
strictly increasing with `d ≥ 1` and all terms `≥ 0`. If a trivially-scaled family (multiply a
known AP by a square) beats a record, that is legal and worth reporting, but also record the
primitive (gcd-reduced) form.

---

## 2. The mathematical lever (read this before writing any search code)

`x² + x·y + y²` is the norm form of the Eisenstein integers `ℤ[ω]` (discriminant −3, class
number 1). Hence:

> **`n ≥ 1` is Loeschian ⟺ every prime `p ≡ 2 (mod 3)` divides `n` to an even power.**

Primes `3` and `p ≡ 1 (mod 3)` are unconstrained. Call `p ≡ 2 (mod 3)` the **bad primes**:
`2, 5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83, 89, 101, …`

Consequences that shape the entire search:

1. **Killing a bad prime is cheap.** If `p | d` and `p ∤ a`, then `p ∤ a + kd` for every `k`,
   so `p` imposes *no* constraint on the AP. Therefore a long AP wants
   `d` divisible by all bad primes up to some bound `B`, with `gcd(a, d)` avoiding them.
2. **Bad primes `p ∤ d` are expensive.** They hit about `n/p` terms, and each hit must have
   `v_p ≥ 2` — i.e. the term must land in a sparse residue class mod `p²`. Small surviving bad
   primes are the dominant obstruction.
3. **Squares are free.** `m²·(Loeschian)` is Loeschian, and `(a, d) → (m²a, m²d)` preserves AP
   validity. So useful normal form: factor out square common factors and search primitive APs.
4. **Loeschian numbers are multiplicatively closed** (norms multiply), density `~C·N/√(log N)` —
   high enough for Green–Tao-style existence, which is why long APs exist at all.
5. **Loeschian numbers are never `2 (mod 3)`** — and this is the single most constraining fact
   in the problem. For `n` coprime to 3, `n = ∏ p_i^{e_i}` gives
   `n ≡ (-1)^(Σ_{p_i ≡ 2 (3)} e_i) (mod 3)`, so Loeschian (all those `e_i` even) forces
   `n ≡ 1 (mod 3)`. Consequences: `3 | d` is forced for any `n ≥ 3`, and the primitive form of
   every AP has `3 ∤ a`, hence **`a ≡ 1 (mod 3)` and every term `≡ 1 (mod 3)`**. (If `3 | a`
   and `3 | d` then `9 | d` and `(a/3, d/3)` is an equally long AP, so reduce.) Measured: with
   `3 ∤ d` the longest run found anywhere is **2**.
6. **Bad primes below `n` in `d`: forced for `p ≤ n/2`, merely very expensive above it.**
   `p | d` is a *theorem* only when `2p ≤ n` (§2b), i.e. `2,5,11,17,23,29 | d` for `n = 58`
   — this is the real hard constraint, and with the forced `3` it gives
   `d ≡ 0 mod 3,741,870`. For `n/2 < p < n` (41, 47, 53 at `n = 58`) omitting `p` from `d` is
   *legal*: one or two terms are divisible by `p` and each such term needs `v_p ≥ 2`, which
   costs a factor `~1/(p+1)`. **Do not state this as a necessary condition** — a verified
   25-term AP exists with `d = 990 = 2·3²·5·11`, omitting both 17 and 23.
   Taking all of 41, 47, 53 into `d` gives the working step base
   `D0 = 3·2·5·11·17·23·29·41·47·53 = 382,160,924,970` (last term `≥ 57·D0 ≈ 2.18e13`).
   Measured trade: including them vs. paying the `1/(p+1)` is close to neutral, because a
   larger `d` inflates every term and the per-term density falls as `1/√(log T)`. Treat the
   choice as a search dimension, not a constraint.
7. **Chinese remainder structure:** choosing `d = 3^e · ∏_{p ≡ 2 (3), p ≤ B} p · (stuff)` and then
   searching `a` over residues is the standard productive shape. Larger `B` buys length but
   inflates `a + 34d`, so **G1 and G3 pull in opposite directions** — expect different
   `d`-families for each.

A fast primality-free test for "is `n` Loeschian": trial-divide/factor `n`, check parity of
`v_p` for bad `p`. For bulk scanning, prefer a **sieve**: mark Loeschian numbers up to `N` by
sieving out `n` with odd `v_p` for some bad `p` (or directly enumerate `x² + xy + y² ≤ N`).

---

## 2b. Gotchas and invariants (read before writing any search code)

Every one of these has a plausible wrong version that silently produces garbage. Encode them as
assertions or unit tests, not as comments.

- **It is `p ≡ 2 (mod 3)`, not `3 (mod 4)`.** The sum-of-two-squares rule is the one most likely
  to come out of muscle memory. Wrong modulus ⇒ a "Loeschian test" that is wrong for
  2, 5, 11, 17, … Unit-test against the known prefix `0, 1, 3, 4, 7, 9, 12, 13, 16, 19, 21`
  *and* against brute-force `x² + xy + y²` enumeration, every time the test changes.
- **Integer width.** Terms can exceed 2^63. numpy `int64` overflows **silently** and will report
  a wrong Loeschian verdict with no error. Any array work must either stay provably below 2^62
  with an assertion on `a + (n-1)d`, or use Python ints / `dtype=object` / gmpy2. The verifier
  must use Python ints, always.
- **Sieve bound.** A sieve built to `N` says nothing about `m > N`; a lookup past the end that
  returns `False` reads as "not Loeschian" and will make a real record look dead. Assert
  `a + (n-1)d ≤ N` before any sieve lookup.
- **Both endpoints and the count.** An `n`-term AP has `n-1` steps; the last term is
  `a + (n-1)d`, not `a + nd`. Off-by-one here inflates or deflates reported lengths. The
  verifier should print `n` and the explicit term list so the count is auditable.
- **`0` and `1` are Loeschian** (`x=y=0`; `x=1,y=0`). Tests that start at 2 or treat 0 as invalid
  will mis-evaluate APs starting low.
- **Negative `x, y` are allowed** but add nothing: `x² + xy + y²` is symmetric and takes the same
  value set over ℤ as over ℕ with the right pairs, so brute-force enumeration must still sweep
  enough of the lattice. Only use enumeration as a cross-check on the prime-factorization test.
- **`d ≥ 1`, strictly increasing, all terms `≥ 0`.** A search that allows `d = 0` finds an
  infinitely long "progression" of one repeated Loeschian number. Assert `d ≥ 1`.
- **The bad-prime cliff is at `p < n`, not `p <= n/2`.** The easy argument ("two terms
  divisible by `p` both need `p^2 |`, their difference is `j·p·d` with `j < p`, so `p | d`")
  only proves `p | d` is *forced* for `n >= 2p`. But for `p < n < 2p` there is still at least
  one term divisible by `p`, and it needs `v_p >= 2`, which costs a factor `~1/(p+1)` in
  probability — three such primes (41, 47, 53 at `n = 58`) cost `~1/10^5` between them.
  Measured (experiments/2026-10-02-calibration): with `29 ∤ d`, run counts fall off a cliff
  immediately past `k = 29` (per-term survival drops 0.55 → 0.36 and hits zero at `k = 32`),
  while with `29 | d` the decay stays a clean 0.53/term. So treat `p | d` as the strongly
  preferred choice for every bad `p < n` — but **it is a cost heuristic, not a necessary
  condition**, and saying otherwise wrongly shrinks the space (see §2-6). Bad primes `p >= n`
  are different and much cheaper: at most one term is divisible, and `a` can be chosen mod `p`
  so that the hit index lands outside `[0, n-1]`.
- **Maximality.** When reporting length `n`, confirm `a - d` (if `≥ 0`) and `a + nd` are *not*
  Loeschian, or say explicitly that the run may extend. Otherwise the record understates itself
  and later sessions re-find the same thing.
- **Never trust the searcher's own verdict.** Searchers use fast filters that are allowed to be
  approximate in one direction. Only `src/verify.py`, factoring from scratch, mints a record.

---

## 2c. Solve this as efficiently as possible

Efficiency is a first-class requirement, not a nicety. The search space is unbounded; the only
way to reach `n ≥ 58` is to spend compute where it pays.

- **Think before you brute-force.** Every hour of CPU should be justified by an argument from §2.
  A pruning insight that shrinks the space is worth more than a faster inner loop. Derive the
  necessary conditions on `d` *first*, then search only what survives them.
- **Never re-search covered ground.** The campaign README + PROGRESS.md record exactly which
  `(a, d)` space has been exhausted. Check it before launching anything.
- **Reuse, don't rebuild.** One sieve, one verifier, one scoring function in `src/`. Campaign
  scripts import them; they do not reimplement the Loeschian test.
- **Cheap test first.** Filter candidates with small-bad-prime residue conditions (O(1) per
  candidate) before any factoring. Factor only what survives. Verify fully only the final hit.
- **Right tool for the inner loop.** Prototype in Python; if a sweep would take more than ~10
  minutes of pure-Python work, move the hot loop to a C extension, numpy vectorization, or a
  small C/Rust program. Measure before optimizing — profile one representative run.
- **Bounded campaigns.** Every search gets an explicit budget (wall clock or iterations) and
  reports partial results on exhaustion. No open-ended runs that produce nothing.
- **Prefer structured construction to blind enumeration.** CRT/backtracking over bad primes and
  greedy `d` refinement have found long APs for others; raw scanning of `(a, d)` will not reach
  58.
- **Checkpoint every long run.** See §3c — an uncheckpointed multi-hour search is wasted the
  first time it is interrupted, and makes "already covered" unknowable.

## 3. Repository contract

```
GOAL.md            this file — objectives, state of the art, protocol
ANSWER.md          THE DELIVERABLE — the human-readable answer (see §9)
records.json       canonical machine-readable record of best-known results (see schema below)
PROGRESS.md        append-only human log: one dated entry per work session
src/               reusable code (sieve, verifier, searchers)
  verify.py        authoritative verifier — single source of truth for "is this AP valid"
experiments/       one subdirectory per search campaign: README.md + script + raw output
results/           raw candidate dumps too big/ugly for records.json
```

Rules:

- **Nothing becomes a record until `src/verify.py` says so.** The verifier must independently
  factor every term (no reliance on the searcher's own sieve) and print PASS/FAIL per term.
- `records.json` schema:
  ```json
  {
    "g1_min_last_term_35": {"a": 0, "d": 0, "last": 0, "n": 35,
                            "verified": true, "method": "experiments/xxx", "date": "YYYY-MM-DD"},
    "g3_longest":          {"a": 0, "d": 0, "n": 0,
                            "verified": true, "method": "...", "date": "YYYY-MM-DD"},
    "history": [ /* every superseded record, same shape — never delete */ ]
  }
  ```
- Numbers may be huge. Use Python ints / strings in JSON when beyond 2^53, and say which.
- Negative results matter: if a campaign proves "no 35-term AP with `a+34d < X` for
  `d` in family F", write it in PROGRESS.md. That is progress toward optimality for G1.

---

## 3b. Environment (pinned — do not re-litigate per session)

- **Python 3.11+**, run as `python3`. Dependencies live in `requirements.txt`; install with
  `python3 -m pip install -r requirements.txt` into a venv at `.venv/` (git-ignored).
- **Allowed and encouraged:** `gmpy2` (fast factoring, `is_square`, big-int arithmetic — the
  single biggest easy win), `numpy` (vectorized sieving, subject to the int64 warning in §2c),
  `sympy` (only as an *independent* cross-check implementation, not in hot loops — it is slow),
  `ortools` if a CP-SAT encoding is tried. A C or Rust helper for an inner loop is fine; commit
  the source plus a one-line build command, never just a binary.
- **`make check` must exist and must pass** before any campaign is launched or any record is
  minted. It runs the unit tests for the Loeschian test, the sieve, and the verifier, including
  the §2c invariants. First command of every session after reading state.
- If a dependency is added, update `requirements.txt` and say so in the PROGRESS.md entry.
- Record machine facts that affect budgeting (core count, available RAM) in the campaign README,
  since sweep sizes are chosen against them.

---

## 3c. Checkpointing (required for any run over ~5 minutes)

Long searches must be resumable, otherwise an interruption loses hours and the
"never re-search covered ground" rule in §2c cannot be honored.

- Every campaign script writes `experiments/<date>-<slug>/checkpoint.json` at least every 60
  seconds: the search-space cursor (e.g. last `d` completed, `a`-range position), best-so-far
  `(a, d, n)`, and a count of candidates examined.
- Every campaign script accepts `--resume` and continues from that checkpoint.
- On normal completion, write `done.json` stating the space **fully exhausted**, in a form a
  later session can compare against (explicit family definition + ranges, not prose).
- Flush on signal: handle SIGINT/SIGTERM by checkpointing before exit, so Ctrl-C is cheap.
- Prefer many bounded shards over one unbounded run — a shard that finishes is a fact; a shard
  that was killed at 80% is only a fact if it checkpointed.

---

## 4. Work loop (what to do when a session starts)

1. Read `records.json` + the last 2–3 entries of `PROGRESS.md`.
2. Pick the next item from §5 "Backlog" (or add one). Priority order is fixed by §1:
   (a) tooling that G3 needs (sieve, verifier, scorer) — only until it is good enough,
   (b) **G3 length pushes toward n ≥ 58** — the default activity,
   (c) G1 sweeps / lower bounds — only when G3 is done or stalled.
3. Run the campaign in `experiments/<date>-<slug>/`.
4. Verify any candidate with `src/verify.py`. Update `records.json` only on PASS.
5. If a record changed, rewrite `ANSWER.md` (§9) in the same session — it must never lag
   `records.json`.
6. Append a PROGRESS.md entry: what was tried, what was found, what was ruled out, next step.
7. Commit (see §7).

A session should always end with the repo in a state where step 1 is enough to resume.

### Stall rule — switch strategy class rather than grinding

Open-ended searches fail by over-investing in the first approach that half-worked. So:

- **If three consecutive campaigns in the same strategy class produce no improvement in best `n`,
  that class is parked.** Write "PARKED: <class>, after <campaigns>" in PROGRESS.md with what it
  topped out at, and move to the next class.
- Strategy classes, in rough order of expected payoff:
  1. CRT / backtracking over bad primes (choose `d`'s bad-prime set, solve for `a` mod `p²`),
  2. greedy / beam refinement of a known good `d` (`d' = k·d`, extend the run),
  3. SAT / CP-SAT encoding of the residue constraints,
  4. literature and OEIS reconnaissance (§5) for a construction we have not thought of,
  5. randomized restart over `d`-families with many `p ≡ 1 (mod 3)` factors.
- A parked class may be revived, but only with a *stated new idea*, not just more compute.
- Conversely: if a class is improving `n` every campaign, stay in it and do not diversify.

---

## 5. Backlog (keep ordered; edit freely)

- [ ] **Recon first (cheap, do before burning CPU):** check OEIS A003136 and the related
      sequences/b-files for arithmetic progressions in Loeschian numbers, and look for anything
      published by **Hugo Pfoertner** (the puzzle's proposer, an active OEIS contributor) on this
      exact question. Note that the IBM page publishes only the record *lengths*, never the
      `(a, d)` values, so there is nothing to copy — what we want is the construction style and
      any known `d`-families. Record findings (including "nothing relevant") in PROGRESS.md and
      cite URLs. One session's worth, max.
- [ ] `src/loeschian.py`: bad-prime Loeschian test + segmented sieve up to `N`. Unit-test
      against the known prefix `0,1,3,4,7,9,12,13,16,19,21,…` and against brute-force
      `x²+xy+y²` enumeration.
- [ ] `src/verify.py`: CLI `verify.py A D N` → factors each of the `N` terms, prints per-term
      verdict and the overall PASS/FAIL. Must be the only thing that can mint a record.
- [ ] Baseline G2: reproduce *some* `n ≥ 42` AP by the §2-1 construction
      (`d = ∏ bad primes ≤ B`, CRT-search `a`). Establishes the pipeline end-to-end.
- [ ] **G3 push to `n ≥ 58` (mandatory)**: hill-climb / beam search over `d` (choice of bad-prime set and
      exponents, plus `3^e` and `p ≡ 1 (3)` factors), scoring `d` by the longest run of
      Loeschian terms achievable. Consider:
      - greedy extension: take a good `(a,d)` and try `d' = k·d`;
      - CRT/backtracking over bad primes `p ≤ n` that divide neither `a` nor `d`
        (each must be handled by `v_p ≥ 2`);
      - CP-SAT / SAT encoding of the residue constraints;
      - search over `d` with many `p ≡ 1 (mod 3)` factors (free multiplicative slack).
- [ ] Prune first: derive necessary conditions on `d` for a given length `n` (which bad primes
      *must* divide `d`), write the argument down, and use it to cut the search space before
      spending compute.
- [ ] G1 sweep (after G3): for 35 terms, enumerate `d` in structured families and scan `a`,
      minimizing `a + 34d`. Record the best and exactly which `(a,d)` space was exhausted.
- [ ] Sanity cross-check: independently re-verify current records with a second implementation
      (e.g. sympy factorint vs. our own) before claiming anything externally.

---

## 6. Subagent protocol

Use subagents where they genuinely buy something, and not otherwise. They are the right tool for
exactly two things here:

1. **Disjoint parallel search.** The `d`-family space partitions cleanly; N agents on N disjoint
   families is a real N× speedup. This is the main use.
2. **Independent cross-checks.** A second agent re-deriving a record with its own implementation
   catches the bug that would otherwise make us claim a wrong answer.

3. **Mathematical / literature reconnaissance, running concurrently with search** (§0). At
   least one research agent should be in flight whenever a search campaign is running, because
   the search alone may not close the gap to `n = 58` inside the deadline. Research agents own
   questions, not `(a, d)` ranges: known records for APs in other norm-form sets (sums of two
   squares especially), admissible-tuple and prime-AP search engineering, Green–Tao-style
   explicit constructions, anything by Hugo Pfoertner, and any paper on APs in multiplicatively
   defined sets. They must return concrete, actionable constructions with citations — or an
   explicit "nothing usable found", which is also a result. A research agent must never be
   trusted on a numeric claim without our own verifier confirming it.

Do **not** spawn a subagent for: a single sweep the orchestrator could run inline, reading one
file, "exploring" without an owned partition, or anything whose result you'd have to redo
yourself to trust. One agent per partition — no redundant duplicates of the same range.

Sizing: a handful of concurrent agents on well-separated families, not dozens on slivers. If a
partition is small enough that agent setup dominates its runtime, merge it into a neighbour.

When spawning:

- Give each subagent: the §2 math summary, the exact `d`-family / `a`-range it owns, the
  scoring objective (G1 vs G3), a wall-clock or iteration budget, and the output path
  `experiments/<date>-<slug>/<agent-id>.json`.
- Partition so that **no two agents search the same `(a, d)` space**; state the partition in
  the campaign README so later sessions know what is already covered.
- Required return value: best `(a, d, n)` found, the space actually exhausted, and runtime.
  "Found nothing" plus an exhausted range is a useful result — ask for it explicitly.
- Subagents **must not** edit `records.json`, `GOAL.md`, or commit. The orchestrator verifies,
  updates records, and commits. This keeps the record store single-writer.
- Run independent subagents in one batch so they execute concurrently.

**Model policy: every agent on this project runs Opus 5 — orchestrator, search agents, and
cross-check agents alike.** Do not downgrade search or cross-check agents to a cheaper tier to
save tokens. The real budget here is CPU-seconds, not tokens: a sharper pruning argument or a
tighter inner loop is worth far more than the token difference, and a cross-check is only
worth running if the checker is as strong as the thing it is checking.

---

## 7. Commit protocol

- Commit after each meaningful unit: new/changed tooling, a completed campaign, a new record,
  or a GOAL/backlog revision. Don't batch unrelated work into one commit.
- **No Claude/AI attribution in commits.** No `Co-Authored-By: Claude`, no "Generated with"
  trailer, no 🤖 — commit messages describe the change only.
- Message style: `<area>: <what changed>`, imperative, with the concrete result in the body.
  - `verify: add authoritative AP verifier`
  - `search: sweep d = 2·5·11·17·k for 35-term APs` (body: best `a+34d`, space covered)
  - `record: G3 n=58 at a=…, d=…` (body: verifier output summary)
- A record-changing commit must touch `records.json` **and** `PROGRESS.md` in the same commit.
- Keep large raw dumps out of git if they exceed a few MB; commit a summary + the script that
  regenerates them instead.
- Work on `main` is fine for solo exploration; use a branch for anything that might break the
  verifier. Never push without being asked.

---

## 8. Definition of done — STOP CONDITION

**The project is done the moment a verified AP with `n ≥ 58` exists.** That is the finish line,
not a milestone. Do not keep searching for longer progressions afterwards, and do not treat this
as an open-ended record chase.

Finishing means, in this order:

1. `src/verify.py` passes on all `n ≥ 58` terms.
2. An independent cross-check (separate implementation, Opus 5) agrees.
3. `records.json` updated, `ANSWER.md` written (§9), PROGRESS.md entry appended, committed.
4. Report to the user. Then stop.

Secondary goals exist only as fallbacks *before* the stop condition is met:

- G2 (`n ≥ 42`) is a smoke test and should fall out of G3 work; just note it when it happens.
- G1 (minimal 35-term endpoint) is something to work on only while a G3 search is blocked or
  stalled. It is not required for done-ness and should not delay it.

---

## 9. The deliverable: `ANSWER.md`

**This is where the user looks for the answer.** It lives at the repo root, is written in plain
prose + numbers, and is self-contained: readable without opening any other file, any script, or
this one. Create it as soon as there is any verified result worth reporting (even `n = 42`), and
rewrite it whenever a record changes — never let it lag `records.json`.

Required contents, in this order:

1. **Headline, one line, at the very top.** The answer and nothing else, e.g.
   `Longest Loeschian AP found: n = 58 terms, a = …, d = …` — and whether it clears the
   mandatory bar of 58.
2. **The progression**, written out: `a`, `d`, the last term `a + (n-1)d`, and all `n` terms
   listed explicitly (or, if they are too large to be readable, the first few / last few plus a
   one-line command that regenerates the full list).
3. **Verification**: the `src/verify.py` invocation used and its PASS output, plus a note that
   an independent cross-check agreed. A reader must be able to re-run one command and confirm it.
4. **Why it works**, in a short paragraph: the structure of `d` (which bad primes divide it,
   what residue class `a` sits in) so the result is understandable, not just assertable.
5. **G1 status** if any 35-term work was done: best `a + 34d` found and the space searched.
6. **How it was found**: which campaign directory, roughly how much compute, and the key pruning
   idea. Two or three sentences.

Keep it short — a page. Detail belongs in PROGRESS.md and the campaign READMEs; `ANSWER.md` is
the answer.

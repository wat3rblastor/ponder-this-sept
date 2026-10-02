# GOAL — IBM Ponder This, September 2026: Loeschian Arithmetic Progressions

Source: https://research.ibm.com/blog/ponder-this-september-2026
(puzzle suggested by Hugo Pfoertner; author Gadi Aleksandrowicz; posted 01 Sep 2026)

This file is the durable contract for this project. Any agent (me, a fresh session, or a
subagent) should be able to read **only this file** plus `records.json` and immediately know
what to do next. Keep it current; it outranks memory and chat history.

---

## 1. The problem (verbatim essentials)

A **Loeschian number** is an integer of the form `x² + x·y + y²` with `x, y ∈ ℤ`.
First few: `0, 1, 3, 4, 7, 9, 12, 13, 16, 19, 21, …` (OEIS A003136).

An arithmetic progression (AP) of Loeschian numbers is described by a start `a` and step `d`;
the `n`-term AP is `a, a+d, a+2d, …, a+(n-1)d`, all of which must be Loeschian.
Example: `1, 7, 13, 19` — `n=4`, `a=1`, `d=6`.

### Objectives

- **G1 (main):** find a **35-term** AP of Loeschian numbers **minimizing the last term**
  `a + 34d`. Report `(a, d)`.
- **G2 (bonus `*`):** find an AP with **≥ 42 terms**. Last term need not be small.
- **G3 (bonus `**`):** find the AP with the **most** terms (> 42).

### Known bar (from the published solvers list — the contest is closed, these are the targets)

Best published lengths: **n = 57** (Jackson La Vallee), then 55, 51, 50, 48, 47, 46, 45, 44.
So: **G2 is table stakes, G3 means n ≥ 58.** G1's optimum is *not* published on the page —
treat our best verified `a+34d` as the record and keep driving it down.

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
5. **Chinese remainder structure:** choosing `d = 3^e · ∏_{p ≡ 2 (3), p ≤ B} p · (stuff)` and then
   searching `a` over residues is the standard productive shape. Larger `B` buys length but
   inflates `a + 34d`, so **G1 and G3 pull in opposite directions** — expect different
   `d`-families for each.

A fast primality-free test for "is `n` Loeschian": trial-divide/factor `n`, check parity of
`v_p` for bad `p`. For bulk scanning, prefer a **sieve**: mark Loeschian numbers up to `N` by
sieving out `n` with odd `v_p` for some bad `p` (or directly enumerate `x² + xy + y² ≤ N`).

---

## 3. Repository contract

```
GOAL.md            this file — objectives, state of the art, protocol
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

## 4. Work loop (what to do when a session starts)

1. Read `records.json` + the last 2–3 entries of `PROGRESS.md`.
2. Pick the next item from §5 "Backlog" (or add one). Prefer: (a) anything that makes the
   verifier/sieve faster or more trustworthy, then (b) G3 length pushes, then (c) G1 lower
   bounds / exhaustive sweeps.
3. Run the campaign in `experiments/<date>-<slug>/`.
4. Verify any candidate with `src/verify.py`. Update `records.json` only on PASS.
5. Append a PROGRESS.md entry: what was tried, what was found, what was ruled out, next step.
6. Commit (see §7).

A session should always end with the repo in a state where step 1 is enough to resume.

---

## 5. Backlog (keep ordered; edit freely)

- [ ] `src/loeschian.py`: bad-prime Loeschian test + segmented sieve up to `N`. Unit-test
      against the known prefix `0,1,3,4,7,9,12,13,16,19,21,…` and against brute-force
      `x²+xy+y²` enumeration.
- [ ] `src/verify.py`: CLI `verify.py A D N` → factors each of the `N` terms, prints per-term
      verdict and the overall PASS/FAIL. Must be the only thing that can mint a record.
- [ ] Baseline G2: reproduce *some* `n ≥ 42` AP by the §2-1 construction
      (`d = ∏ bad primes ≤ B`, CRT-search `a`). Establishes the pipeline end-to-end.
- [ ] G1 sweep: for 35 terms, enumerate `d` in structured families and scan `a`, minimizing
      `a + 34d`. Record the best and the search bound actually covered (for a provable-ish
      lower bound, state exactly which `(a,d)` space was exhausted).
- [ ] G3 push to `n ≥ 58`: hill-climb / beam search over `d` (choice of bad-prime set and
      exponents, plus `3^e` and `p ≡ 1 (3)` factors), scoring `d` by the longest run of
      Loeschian terms achievable. Consider:
      - greedy extension: take a good `(a,d)` and try `d' = k·d`;
      - CRT/backtracking over bad primes `p ≤ n` that divide neither `a` nor `d`
        (each must be handled by `v_p ≥ 2`);
      - CP-SAT / SAT encoding of the residue constraints;
      - search over `d` with many `p ≡ 1 (mod 3)` factors (free multiplicative slack).
- [ ] Optimality pressure for G1: derive necessary conditions on `d` (which bad primes *must*
      divide `d` for 35 terms) to prune the sweep, and write the argument down.
- [ ] Sanity cross-check: independently re-verify current records with a second implementation
      (e.g. sympy factorint vs. our own) before claiming anything externally.

---

## 6. Subagent protocol

Subagents are encouraged — this problem parallelizes cleanly over disjoint `d`-families.

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

---

## 7. Commit protocol

- Commit after each meaningful unit: new/changed tooling, a completed campaign, a new record,
  or a GOAL/backlog revision. Don't batch unrelated work into one commit.
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

## 8. Definition of done

- G1: a verified 35-term AP whose last term we believe is minimal, with the searched space
  documented well enough that the claim is auditable.
- G2: verified `n ≥ 42`. (Should fall out of G3 work.)
- G3: verified `n ≥ 58`, beating the published record of 57.

Since G3 has no ceiling, this project is open-ended by design: each session should leave either
a better record, a larger exhausted search space, or better tooling.

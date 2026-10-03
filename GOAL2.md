# GOAL2 — the method contract: how to reach n ≥ 58 without immense computation

`GOAL.md` defines **what** (the problem, the mathematics, the verifier contract, the stop
condition). This file defines **how we work**. Where the two disagree about method, this file
wins. Where they disagree about the problem or what counts as a verified result, `GOAL.md` wins.

---

## 1. The premise — not a hypothesis, an axiom

> **A 58-term arithmetic progression of Loeschian numbers can be found with the computation
> available on this laptop. Someone has already found a 57-term progression without immense
> computation. Therefore a route exists that we have not yet seen.**

This is the starting assumption of every session, and it is **not up for re-litigation**. It is
not weakened by any amount of accumulated negative evidence. In particular:

- "The calibrated model says this needs N GPU-hours" is a statement about **the search we
  happen to have built**, never about the problem. It is evidence that the current method is
  wrong, not that the target is out of reach.
- Every impossibility result in this repo is an impossibility result *about a specific
  construction under specific assumptions*. Each one narrows the space; none of them closes it.
- "Rent more GPUs" is **not an answer** and must never be proposed as the plan. If the only
  thing a plan adds is throughput, it is not a plan.

The orchestrator's standing belief is that the gap between what we can compute and what we need
is a **gap in understanding**, and the job is to close it by thinking, not by spending.

**Compute budget: this laptop.** 8 cores, 8 GB, the integrated GPU. Any plan whose success
depends on hardware we do not have is out of scope by definition. This is a constraint on
*method*, and it is the point: it forces every plan to be an idea rather than a purchase.

---

## 2. Stop condition

Exactly `GOAL.md` §8: a **verified `n ≥ 58`** — passing `src/verify.py` *and*
`src/crosscheck.py`, recorded in `records.json`, written up in `ANSWER.md`, committed.

Nothing else ends the loop. Not a good idea, not a promising measurement, not a partial result,
not running out of plans. If the plan queue empties, the orchestrator's job is to refill it.

---

## 3. Division of labour — the rule that makes this work

**The orchestrator does not implement. Ever.**

| | orchestrator | agents |
|---|---|---|
| generates hypotheses | ✅ | only when asked for a survey |
| designs the falsifying experiment | ✅ | ✖ |
| writes code | ✖ | ✅ |
| runs searches / measurements | ✖ | ✅ |
| reads raw logs | ✖ (reads *verdicts*) | ✅ |
| decides what the result means | ✅ | ✖ (reports, does not conclude the project) |
| decides what happens next | ✅ | ✖ |

The reason for the rule: an orchestrator that implements gets absorbed into the implementation.
It starts optimizing the thing it is building instead of asking whether the thing should exist.
It mistakes "my program now runs 2× faster" for progress. Keeping the orchestrator's hands empty
keeps its attention on the only question that matters — *is this the right idea?* — and means a
failed plan costs an agent's time, not the orchestrator's conviction.

A corollary: **the orchestrator never runs out of energy, because it never does the work.** When
an agent returns "this does not work", the orchestrator has lost nothing but an assumption, and
assumptions are what we are trying to spend.

---

## 4. The loop

```
1. STATE a hypothesis: a specific claim about structure that, if true, makes n = 58 cheap.
2. DESIGN the cheapest experiment that could prove it FALSE.
      - it must be falsifiable, measurable on this laptop, and bounded (< ~1 h of agent time)
      - state in advance the number that decides it, and which way
3. DISPATCH one agent per independent hypothesis, in parallel. Give each the brief in §5.
4. READ the verdicts. For each:
      - CONFIRMED  -> design the follow-up that turns it into terms on the ground
      - REFUTED    -> record WHICH ASSUMPTION died (§6), then go to 5
      - AMBIGUOUS  -> redesign the measurement; ambiguity is a bug in the experiment
5. RE-PLAN. Generate new hypotheses (§7 lists where to look). Return to 1.
```

Run 2–4 hypotheses concurrently; they are independent and agents are cheap. Never serialize on
one idea. Never wait on a running experiment without starting another.

**Every cycle must end with the plan queue non-empty.** If the orchestrator cannot think of a
next hypothesis, that is itself the problem to attack: dispatch a survey agent whose brief is
"here is everything we have ruled out and the exact assumption each rested on; find the
assumption we have not questioned."

---

## 5. Agent brief template

Every agent gets a self-contained brief. Agents do not share the orchestrator's context and must
not need it.

```
OBJECTIVE   one sentence: the single question this agent answers.
BACKGROUND  the minimum mathematics needed, restated in full. Point at GOAL.md §2 and §2b.
            Never say "as discussed".
TASK        exactly what to build or measure.
DECIDES     the number that settles it, and the threshold. e.g. "if the matched measurement
            gives < 5x per unit of work, the idea is dead; report that plainly."
BUDGET      wall-clock and the fact that the only hardware is this laptop.
OUTPUT      experiments/<date>-<slug>/ : code, raw data, and a README stating the verdict,
            the number, and how to reproduce it.
RULES       - a clean negative is a success; say it in the first line
            - never report a numeric claim this repo's verifiers have not confirmed
            - do not edit records.json, GOAL.md, GOAL2.md, or ANSWER.md; do not commit
            - if the experiment turns out to measure the wrong thing, say so and stop
```

Verification is non-negotiable and belongs to the orchestrator: a candidate progression is real
only when `src/verify.py` and `src/crosscheck.py` both pass, via `tools/promote.py`.

---

## 6. The register — and how to use it

`PROGRESS.md` holds what has been ruled out. **Its purpose is not to stop us; it is to show us
which assumptions we have already paid for.** Every entry must name the assumption it rested on,
because the next idea usually lives inside one of them.

Current entries, each with the assumption to attack:

| ruled out | rested on the assumption that… |
|---|---|
| scaling cannot extend or repair a run | the AP is primitive and `c` is a single global constant |
| no linear `f(k)` is identically a norm | the representation `z_k` moves linearly in `k` |
| CRT forces ≤ 2/58 of each term's digits | the forced factors are distinct primes, one per term |
| meet-in-the-middle does not apply | `d` is fixed, so `a` is the only unknown |
| automatic terms lose by ~10⁵ | one square condition, shapes `a = x², d = 3m²` and `a = 3x²+m², d = 24m²` |
| the square option at 41/47/53 loses | the `p²` classes are pinned in stage 1, costing modulus |
| the pin set 59..113 is optimal | term size is floored by `64·MOD`, i.e. by *this* enumerator |
| large `d` alone gives nothing | "large" meant `d = 3·∏(bad p ≤ y)` |

**Re-running a dead idea unchanged is forbidden. Attacking the assumption underneath it is
exactly the job.**

---

## 7. Where to look for the next hypothesis

Not a backlog to work through — a set of directions that have *not* been exhausted. The
orchestrator should add to this list every cycle.

- **The automatic-terms identity, generalized.** `3·c·d·t_{k₀}` a perfect square makes
  `t_{k₀ ± cj²}` free. Only two shapes and `c = 1` have been tried. The free index set is
  `k₀ + c·j²`; a *union* over several `(c, k₀)` is a covering problem over `[0, 57]` that has
  never been posed properly, let alone optimized. Ask: what is the maximum coverage of 58
  consecutive indices by unions of quadratic progressions whose square conditions are
  simultaneously satisfiable at small parameters?
- **The second mechanism.** The pentagonal shape makes terms free by `t_k = 3x² + (ms)²`, i.e. by
  *being represented*, not by the `c·j²` rule. These are two different sources of free terms. Are
  there others? Enumerate systematically rather than by example.
- **Attack the enumerator, not the search.** Every "optimal configuration" result above is
  conditioned on `T ≥ 64·MOD` — a property of *our* stage-1 enumerator. An enumerator that
  produced admissible `a` **small relative to its modulus** would break that floor and, by the
  pin-set table, is worth up to 200× at the same term size. That is a lattice/CVP question
  (enumerate CRT-admissible residues below a bound), not a number-theory question.
- **Our own 55s.** Two verified distinct 55-term progressions exist. Both have `59 | d`. Look for
  more structure in them: factorizations, residues, relations. The record holder's remark that
  "every 57 contains three 55s" is a statement about *our own data* we have not mined.
- **Run the problem backwards.** Instead of choosing `(a, d)` and testing terms, choose the 58
  *factorizations* and solve for `(a, d)`. The terms must be 58 numbers in AP each of whose
  squarefree kernel avoids all primes ≡ 2 mod 3. Posed as a constraint system over the kernels,
  what does it force?
- **Small cases, hard.** Exhaustive optimal APs are known for `n ≤ 28` and we proved the `n = 35`
  minimum. These are the only ground truth about what long Loeschian APs *look like*. They have
  been glanced at ("they look random"), never dissected.

---

## 8. Standing rules

1. **Never propose throughput as a plan.** More cores, more GPUs, more hours, a faster kernel —
   none of these is a hypothesis. They are permitted only as the *execution* of a confirmed
   structural idea, never as the idea.
2. **Never conclude the target is unreachable.** That conclusion is excluded by §1. A negative
   result refutes a method; report it as such and re-plan.
3. **One agent per independent question**, dispatched in parallel, each with its own
   `experiments/` directory.
4. **Measure, don't argue.** When the orchestrator and a model disagree, the matched experiment
   decides. Prefer a 20-minute measurement over an hour of algebra.
5. **A negative result must name its assumption** (§6), or it is not finished.
6. **Record every cycle in PROGRESS.md**: hypothesis, experiment, number, verdict, next move.
7. **The loop does not end until `GOAL.md` §8 is satisfied.** If every plan has failed, the
   correct next action is to generate more plans.

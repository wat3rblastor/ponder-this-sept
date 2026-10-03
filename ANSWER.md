# Answers

**Main challenge (35 terms, smallest last term): `a = 219830911`, `d = 2709630`,
last term `311958331`.**

**Longest progression found: `n = 55` terms, `a = 11687581876345393`, `d = 202927451159070`.**

Both are verified two independent ways (§3). The mandatory internal target of `n ≥ 58` is
**NOT MET** — §5 explains why, with numbers.

| goal | target | status |
|---|---|---|
| IBM main challenge | 35 terms, minimise last term | **last term = 311958331** |
| IBM bonus `*` | n ≥ 42 | **CLEARED** — n = 55 |
| IBM bonus `**` | longest found | **n = 55** (published leader: 57) |
| GOAL.md G3 (mandatory) | n ≥ 58 | **NOT MET** — best verified n = 55 |

---

## 1. The 35-term progression with the smallest last term

```
a    = 219830911
d    = 2709630  = 2 * 3^2 * 5 * 7 * 11 * 17 * 23
n    = 35 terms
last = a + 34d = 311958331
```

```
219830911  222540541  225250171  227959801  230669431
233379061  236088691  238798321  241507951  244217581
246927211  249636841  252346471  255056101  257765731
260475361  263184991  265894621  268604251  271313881
274023511  276733141  279442771  282152401  284862031
287571661  290281291  292990921  295700551  298410181
301119811  303829441  306539071  309248701  311958331
```

### How far this is proved

For `n = 35` the bad primes `p` with `2p ≤ 35` are `2, 5, 11, 17`, and `p | d` is a **theorem**
for those (two terms divisible by `p` would each need `p² |`, yet their difference is `j·p·d`
with `j < p`). With the forced factor 3, **every** valid step is a multiple of
`5610 = 3·2·5·11·17`. So sweeping `d = 5610·m` over all `m` and all `a` is a *complete* search,
and these are proofs rather than best-effort claims:

- **No 35-term Loeschian AP has last term ≤ 2·10⁸** (complete family, `m = 1..1049`, every `a`;
  the longest runs in that range are 27–29 terms).
- **No 35-term AP with last term < 311958331 has `23·29 | d/5610`** — the sub-family in which
  the two remaining sub-35 bad primes are also killed, ~130× likelier to yield a hit than a
  generic step. Exhausted.
- A complete sweep of all `m ≤ 7062` with every term `≤ 311958331` was running at the
  deadline; `experiments/2026-10-03-g1/*.log` records exactly how far it reached.

An independently calibrated model (fitted to brute force for `n = 12..24` and validated against
the published record staircase) predicts the true optimum at last term `≈1.3·10⁹`. This value
sits on that estimate, so it is likely optimal or very close — but only the bounds above are
proved.

---

## 2. The longest progression: n = 55

```
a    = 11687581876345393
d    = 202927451159070
     = 2 * 3^3 * 5 * 11 * 17 * 23 * 29 * 41 * 47 * 53 * 59
n    = 55 terms
last = a + 54d = 22645664238935173   (17 digits)
```

```
11687581876345393  11890509327504463  12093436778663533
12296364229822603  12499291680981673  12702219132140743
12905146583299813  13108074034458883  13311001485617953
13513928936777023  13716856387936093  13919783839095163
14122711290254233  14325638741413303  14528566192572373
14731493643731443  14934421094890513  15137348546049583
15340275997208653  15543203448367723  15746130899526793
15949058350685863  16151985801844933  16354913253004003
16557840704163073  16760768155322143  16963695606481213
17166623057640283  17369550508799353  17572477959958423
17775405411117493  17978332862276563  18181260313435633
18384187764594703  18587115215753773  18790042666912843
18992970118071913  19195897569230983  19398825020390053
19601752471549123  19804679922708193  20007607373867263
20210534825026333  20413462276185403  20616389727344473
20819317178503543  21022244629662613  21225172080821683
21428099531980753  21631026983139823  21833954434298893
22036881885457963  22239809336617033  22442736787776103
22645664238935173
```

Regenerate: `python3 -c "a,d=11687581876345393,202927451159070; print([a+k*d for k in range(55)])"`

Maximal at both ends: True — so 55 is its true length,
not a truncation.

---

## 3. Verification

```
$ python3 src/verify.py 219830911 2709630 35 --maximal
OVERALL: PASS  (n=35, a=219830911, d=2709630)

$ python3 src/verify.py 11687581876345393 202927451159070 55 --maximal
OVERALL: PASS  (n=55, a=11687581876345393, d=202927451159070)

$ python3 src/crosscheck.py 219830911 2709630 35
CROSS-CHECK: PASS  (35 terms, each with a verified x^2+x*y+y^2 representation)

$ python3 src/crosscheck.py 11687581876345393 202927451159070 55
CROSS-CHECK: PASS  (55 terms, each with a verified x^2+x*y+y^2 representation)
```

- `src/verify.py` factors every term from scratch in Python ints (its own Miller–Rabin and
  Pollard–Brent; no sieve, no fixed-width arithmetic) and checks that every prime `p ≡ 2 (mod 3)`
  occurs to an even power.
- `src/crosscheck.py` does **not use that criterion at all**. It *constructs* integers `x, y`
  with `x² + xy + y² = t` for every term — Eisenstein-integer arithmetic, Cornacchia for primes
  `≡ 1 (mod 3)` — and checks the identity by direct multiplication. An explicit representation is
  positive proof, independent of any theory about bad primes.
- `make check` runs the test suite, in which three further implementations agree on every
  `n ≤ 20000`.
- `tools/promote.py` refuses to record anything that fails either check.

---

## 4. Why the steps look like this

A Loeschian number is **never `≡ 2 (mod 3)`**: for `t` coprime to 3,
`t ≡ (−1)^(number of bad prime factors, with multiplicity) (mod 3)`, and "Loeschian" forces every
such exponent even, hence that count even. So `3 | d` is forced and all terms lie in one class
mod 3 — measured, with `3 ∤ d` the longest run found anywhere is **2**. This lifts the per-term
probability from ~0.14 to ~0.6.

Each bad prime `p` with `2p ≤ n` must then divide `d` (the theorem above); each with
`n/2 < p < n` is not forced but costs `~1/(p+1)` if omitted. Hence both steps are 3 times a
product of the small bad primes. Finally `a` is chosen so that none of 59, 71, 83, 89, 101, 107,
113 divides *any* term: each can hit at most one index in the window, and `a`'s residue pushes
the hit outside it.

---

## 5. Why n ≥ 58 was not reached

Not for want of effort or algorithm — the target is out of range for this machine:

- The engine enumerates the admissible residues of `a` **directly by nested additive loops**, so
  an inadmissible candidate is never constructed (~3·10⁴ saved), then settles 64 candidates per
  modulo with precomputed 64-bit masks over ~600 further primes. An independent analysis showed
  this is at a *provable floor*: the minimum number of residues any such scheme must enumerate is
  exactly `∏(q−58) = 1167541375`, which is what it enumerates.
- Stage 1+2 was ported to the M2's integrated GPU (`src/c/apsearch_gpu.m`), measured at
  **4.1·10⁸ residues/s — about 9× the whole 8-core CPU**.
- Even so: a calibrated model puts **the smallest last term able to carry 58 terms at ~4.5·10¹⁸**,
  while this search covers terms to `7.3·10¹⁶`. A 58-term progression almost certainly does not
  *exist* in the region searched, and reaching `~10¹⁸` costs on the order of 10³ times more
  compute than 12 hours here. The same model reproduces the published record staircase
  (48, 50, 51, 55, 57), which is the main reason to trust it.
- Two tempting shortcuts are **provably impossible**, not merely unlikely: no constant `c` can
  extend a progression or repair a near miss (scaling XORs every term's bad-parity vector
  identically, and in a primitive AP no two terms can share a nonempty vector), and no
  meet-in-the-middle exists (for fixed `d` every constraint is a congruence on the single unknown
  `a`). PROGRESS.md lists everything ruled out and why.

## 6. Where everything is

- `records.json` — machine-readable records and full history.
- `PROGRESS.md` — dated log: what was tried, found, and ruled out, with the cost model.
- `experiments/` — campaigns and raw per-work-unit output.
- `src/` — Loeschian core, authoritative verifier, constructive cross-check, and four searchers
  (flat sieve, class-compressed, CPU two-stage CRT, and the GPU port).

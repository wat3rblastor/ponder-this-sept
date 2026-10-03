# Answers

**Main challenge (35 terms, smallest last term): `a = 219830911`, `d = 2709630`,
last term `311958331`.**

**Longest progression found: `n = 47` terms, `a = 2646171143023357`, `d = 78342989618850`.**

Both are verified two independent ways (§3). The mandatory internal target of `n ≥ 58` is
**NOT MET** — §5 explains why, with numbers.

| goal | target | status |
|---|---|---|
| IBM main challenge | 35 terms, minimise last term | **last term = 311958331** |
| IBM bonus `*` | n ≥ 42 | **CLEARED** — n = 47 |
| IBM bonus `**` | longest found | **n = 47** (published leader: 57) |
| GOAL.md G3 (mandatory) | n ≥ 58 | **NOT MET** — best verified n = 47 |

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

## 2. The longest progression: n = 47

```
a    = 2646171143023357
d    = 78342989618850
     = 2 * 3 * 5^2 * 11 * 17 * 23 * 29 * 41^2 * 47 * 53
n    = 47 terms
last = a + 46d = 6249948665490457   (16 digits)
```

```
2646171143023357  2724514132642207  2802857122261057
2881200111879907  2959543101498757  3037886091117607
3116229080736457  3194572070355307  3272915059974157
3351258049593007  3429601039211857  3507944028830707
3586287018449557  3664630008068407  3742972997687257
3821315987306107  3899658976924957  3978001966543807
4056344956162657  4134687945781507  4213030935400357
4291373925019207  4369716914638057  4448059904256907
4526402893875757  4604745883494607  4683088873113457
4761431862732307  4839774852351157  4918117841970007
4996460831588857  5074803821207707  5153146810826557
5231489800445407  5309832790064257  5388175779683107
5466518769301957  5544861758920807  5623204748539657
5701547738158507  5779890727777357  5858233717396207
5936576707015057  6014919696633907  6093262686252757
6171605675871607  6249948665490457
```

Regenerate: `python3 -c "a,d=2646171143023357,78342989618850; print([a+k*d for k in range(47)])"`

Maximal at both ends: True — so 47 is its true length,
not a truncation.

---

## 3. Verification

```
$ python3 src/verify.py 219830911 2709630 35 --maximal
OVERALL: PASS  (n=35, a=219830911, d=2709630)

$ python3 src/verify.py 2646171143023357 78342989618850 47 --maximal
OVERALL: PASS  (n=47, a=2646171143023357, d=78342989618850)

$ python3 src/crosscheck.py 219830911 2709630 35
CROSS-CHECK: PASS  (35 terms, each with a verified x^2+x*y+y^2 representation)

$ python3 src/crosscheck.py 2646171143023357 78342989618850 47
CROSS-CHECK: PASS  (47 terms, each with a verified x^2+x*y+y^2 representation)
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

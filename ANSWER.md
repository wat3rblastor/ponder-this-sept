# Answers

**Main challenge (35 terms, smallest last term): `a = 219830911`, `d = 2709630`,
last term `311958331`.**

**Longest progression found: `n = 56` terms, `a = 240533583825657151`, `d = 13967981807653500`.**

Both are verified two independent ways (§3). The mandatory internal target of `n ≥ 58` is
**NOT MET** — see §5b for the status of the search.

| goal | target | status |
|---|---|---|
| IBM main challenge | 35 terms, minimise last term | **last term = 311958331** |
| IBM bonus `*` | n ≥ 42 | **CLEARED** — n = 56 |
| IBM bonus `**` | longest found | **n = 56** (published leader: 57) |
| GOAL.md G3 (mandatory) | n ≥ 58 | **NOT MET** — best verified n = 56 |

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

## 2. The longest progression: n = 56

```
a    = 240533583825657151
d    = 13967981807653500
     = 2^2 * 3 * 5^3 * 11 * 17^2 * 23 * 29 * 41 * 43 * 47 * 53
n    = 56 terms
last = a + 55d = 1008772583246599651   (19 digits)
```

```
 240533583825657151   254501565633310651   268469547440964151
 282437529248617651   296405511056271151   310373492863924651
 324341474671578151   338309456479231651   352277438286885151
 366245420094538651   380213401902192151   394181383709845651
 408149365517499151   422117347325152651   436085329132806151
 450053310940459651   464021292748113151   477989274555766651
 491957256363420151   505925238171073651   519893219978727151
 533861201786380651   547829183594034151   561797165401687651
 575765147209341151   589733129016994651   603701110824648151
 617669092632301651   631637074439955151   645605056247608651
 659573038055262151   673541019862915651   687509001670569151
 701476983478222651   715444965285876151   729412947093529651
 743380928901183151   757348910708836651   771316892516490151
 785284874324143651   799252856131797151   813220837939450651
 827188819747104151   841156801554757651   855124783362411151
 869092765170064651   883060746977718151   897028728785371651
 910996710593025151   924964692400678651   938932674208332151
 952900656015985651   966868637823639151   980836619631292651
 994804601438946151  1008772583246599651
```

Regenerate: `python3 -c "a,d=240533583825657151,13967981807653500; print([a+k*d for k in range(56)])"`

Maximal at both ends: True — so 56 is its true length,
not a truncation.

---

## 3. Verification

```
$ python3 src/verify.py 219830911 2709630 35 --maximal
OVERALL: PASS  (n=35, a=219830911, d=2709630)

$ python3 src/verify.py 240533583825657151 13967981807653500 56 --maximal
OVERALL: PASS  (n=56, a=240533583825657151, d=13967981807653500)

$ python3 src/crosscheck.py 219830911 2709630 35
CROSS-CHECK: PASS  (35 terms, each with a verified x^2+x*y+y^2 representation)

$ python3 src/crosscheck.py 240533583825657151 13967981807653500 56
CROSS-CHECK: PASS  (56 terms, each with a verified x^2+x*y+y^2 representation)
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

## 5. Every distinct progression of 55 or more terms found

Each row is a different progression: rescaled copies (`m·a`, `m·d` for a Loeschian multiplier
`m`, which the search re-finds routinely) are reduced to the primitive form and listed once.
Every row passes both `src/verify.py` and `src/crosscheck.py` at the stated length.

| n | a | d | last term | d factored |
|---|---|---|---|---|
| 56 | `240533583825657151` | `13967981807653500` | `1008772583246599651` | 2^2 * 3 * 5^3 * 11 * 17^2 * 23 * 29 * 41 * 43 * 47 * 53 |
| 55 | `11687581876345393` | `202927451159070` | `22645664238935173` | 2 * 3^3 * 5 * 11 * 17 * 23 * 29 * 41 * 47 * 53 * 59 |
| 55 | `296246969176050787` | `11950172123811900` | `941556263861893387` | 2^2 * 3 * 5^2 * 11 * 17 * 23 * 29 * 41 * 47 * 53^2 * 59 |
| 55 | `248561417402960821` | `54046344492032310` | `3167064019972705561` | 2 * 3^2 * 5 * 11 * 17^2 * 23 * 29 * 41 * 47^2 * 53 * 59 |

Re-check any row: `python3 src/verify.py <a> <d> <n> --maximal`

## 5b. Status of the n ≥ 58 target

**NOT MET**. The search runs on 8 NVIDIA RTX PRO 6000 Blackwell GPUs at roughly 6·10¹¹ residues per
second (`src/c/apsearch_cuda.cu`, launched by `tools/multi_gpu.sh`). The engine enumerates the
admissible residues of `a` directly by nested additive loops, so an inadmissible candidate is
never constructed, then settles 64 candidates per residue with precomputed bitmasks over every
bad prime up to 10⁴, and confirms survivors with an exact test. A calibrated model puts the
expected yield at a few tenths of a 58-term progression per hour at that rate; PROGRESS.md has
the model, the measurements, and everything that was ruled out. Two shortcuts are provably
impossible: no constant `c` can extend a progression or repair a near miss, and no
meet-in-the-middle exists (for fixed `d` every constraint is a congruence on the single unknown
`a`).

## 6. Where everything is

- `records.json` — machine-readable records and full history.
- `PROGRESS.md` — dated log: what was tried, found, and ruled out, with the cost model.
- `experiments/` — campaigns and raw per-work-unit output.
- `src/` — Loeschian core, authoritative verifier, constructive cross-check, and four searchers
  (flat sieve, class-compressed, CPU two-stage CRT, and the Metal and CUDA GPU engines).

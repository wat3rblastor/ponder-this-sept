# Answers

**Main challenge (35 terms, smallest last term): `a = 219830911`, `d = 2709630`,
last term `311958331`.**

**Longest progression found: `n = 57` terms, `a = 416580650729355319`, `d = 42208909841086560`.**

Both are verified two independent ways (§3). The mandatory internal target of `n ≥ 58` is
**NOT MET** — see §5b for the status of the search.

| goal | target | status |
|---|---|---|
| IBM main challenge | 35 terms, minimise last term | **last term = 311958331** |
| IBM bonus `*` | n ≥ 42 | **CLEARED** — n = 57 |
| IBM bonus `**` | longest found | **n = 57** (published leader: 57) |
| GOAL.md G3 (mandatory) | n ≥ 58 | **NOT MET** — best verified n = 57 |

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

## 2. The longest progression: n = 57

```
a    = 416580650729355319
d    = 42208909841086560
     = 2^5 * 3^3 * 5 * 11 * 13 * 17 * 23 * 29 * 41 * 47 * 53 * 59
n    = 57 terms
last = a + 56d = 2780279601830202679   (19 digits)
```

```
 416580650729355319   458789560570441879   500998470411528439
 543207380252614999   585416290093701559   627625199934788119
 669834109775874679   712043019616961239   754251929458047799
 796460839299134359   838669749140220919   880878658981307479
 923087568822394039   965296478663480599  1007505388504567159
1049714298345653719  1091923208186740279  1134132118027826839
1176341027868913399  1218549937709999959  1260758847551086519
1302967757392173079  1345176667233259639  1387385577074346199
1429594486915432759  1471803396756519319  1514012306597605879
1556221216438692439  1598430126279778999  1640639036120865559
1682847945961952119  1725056855803038679  1767265765644125239
1809474675485211799  1851683585326298359  1893892495167384919
1936101405008471479  1978310314849558039  2020519224690644599
2062728134531731159  2104937044372817719  2147145954213904279
2189354864054990839  2231563773896077399  2273772683737163959
2315981593578250519  2358190503419337079  2400399413260423639
2442608323101510199  2484817232942596759  2527026142783683319
2569235052624769879  2611443962465856439  2653652872306942999
2695861782148029559  2738070691989116119  2780279601830202679
```

Regenerate: `python3 -c "a,d=416580650729355319,42208909841086560; print([a+k*d for k in range(57)])"`

Maximal at both ends: True — so 57 is its true length,
not a truncation.

---

## 3. Verification

```
$ python3 src/verify.py 219830911 2709630 35 --maximal
OVERALL: PASS  (n=35, a=219830911, d=2709630)

$ python3 src/verify.py 416580650729355319 42208909841086560 57 --maximal
OVERALL: PASS  (n=57, a=416580650729355319, d=42208909841086560)

$ python3 src/crosscheck.py 219830911 2709630 35
CROSS-CHECK: PASS  (35 terms, each with a verified x^2+x*y+y^2 representation)

$ python3 src/crosscheck.py 416580650729355319 42208909841086560 57
CROSS-CHECK: PASS  (57 terms, each with a verified x^2+x*y+y^2 representation)
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
| 57 | `416580650729355319` | `42208909841086560` | `2780279601830202679` | 2^5 * 3^3 * 5 * 11 * 13 * 17 * 23 * 29 * 41 * 47 * 53 * 59 |
| 57 | `605952547032888577` | `39345378030286350` | `2809293716728924177` | 2 * 3 * 5^2 * 11 * 17 * 23 * 29 * 41 * 47 * 53 * 59 * 349 |
| 56 | `240533583825657151` | `13967981807653500` | `1008772583246599651` | 2^2 * 3 * 5^3 * 11 * 17^2 * 23 * 29 * 41 * 43 * 47 * 53 |
| 56 | `410226150109637989` | `45523391543351370` | `2914012684993963339` | 2 * 3^2 * 5 * 11 * 17 * 23 * 29 * 41 * 47 * 53 * 59 * 673 |
| 56 | `358492265928380449` | `113517463314013770` | `6601952748199137799` | 2 * 3 * 5 * 11 * 17^2 * 23 * 29 * 41 * 47 * 53 * 101 * 173 |
| 55 | `11687581876345393` | `202927451159070` | `22645664238935173` | 2 * 3^3 * 5 * 11 * 17 * 23 * 29 * 41 * 47 * 53 * 59 |
| 55 | `296246969176050787` | `11950172123811900` | `941556263861893387` | 2^2 * 3 * 5^2 * 11 * 17 * 23 * 29 * 41 * 47 * 53^2 * 59 |
| 55 | `248561417402960821` | `54046344492032310` | `3167064019972705561` | 2 * 3^2 * 5 * 11 * 17^2 * 23 * 29 * 41 * 47^2 * 53 * 59 |

Re-check any row: `python3 src/verify.py <a> <d> <n> --maximal`

### Rescaled copies of the rows above, as the search found them

Each of these is a row of the table above with `a` and `d` both multiplied by the stated
Loeschian multiplier. They are valid progressions of the stated length (each verified by both
checkers) but not new ones.

| n | a | d | last term | multiplier | primitive a |
|---|---|---|---|---|---|
| 55 | `151938564392490109` | `2638056865067910` | `294393635106157249` | 13 | `11687581876345393` |
| 55 | `222064055650562467` | `3855621572022330` | `430267620539768287` | 19 | `11687581876345393` |
| 55 | `362315038166707183` | `6290750985931170` | `702015591406990363` | 31 | `11687581876345393` |
| 55 | `432440529424779541` | `7508315692885590` | `837889576840601401` | 37 | `11687581876345393` |
| 55 | `502566020682851899` | `8725880399840010` | `973763562274212439` | 43 | `11687581876345393` |
| 55 | `572691511940924257` | `9943445106794430` | `1109637547707823477` | 49 | `11687581876345393` |
| 55 | `712942494457068973` | `12378574520703270` | `1381385518575045553` | 61 | `11687581876345393` |
| 55 | `783067985715141331` | `13596139227657690` | `1517259504008656591` | 67 | `11687581876345393` |
| 55 | `853193476973213689` | `14813703934612110` | `1653133489442267629` | 73 | `11687581876345393` |
| 55 | `923318968231286047` | `16031268641566530` | `1789007474875878667` | 79 | `11687581876345393` |
| 55 | `1133695442005503121` | `19683962762429790` | `2196629431176711781` | 97 | `11687581876345393` |
| 55 | `1203820933263575479` | `20901527469384210` | `2332503416610322819` | 103 | `11687581876345393` |
| 55 | `1273946424521647837` | `22119092176338630` | `2468377402043933857` | 109 | `11687581876345393` |
| 55 | `1484322898295864911` | `25771786297201890` | `2875999358344766971` | 127 | `11687581876345393` |

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

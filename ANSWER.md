# Answers

**Main challenge (35 terms, smallest last term): `a = 73415383`, `d = 37418700`, last term `1345651183`.**
**Longest progression found: `n = 43` terms, `a = 5413537078288507`, `d = 48916598396160`.**

Both are verified two independent ways (see Verification). The mandatory internal target of
`n ≥ 58` was **not reached** — §"Why n ≥ 58 was not reached" explains why, with numbers.

| goal | target | status |
|---|---|---|
| IBM main challenge | 35 terms, minimise last term | **last term = 1345651183** (verified; exhaustive bounds below) |
| IBM bonus `*` | n ≥ 42 | **CLEARED** — n = 43, verified |
| IBM bonus `**` | longest found | **n = 43** (published leader: 57) |
| GOAL.md G3 (mandatory) | n ≥ 58 | **NOT MET** — best verified n = 43 |

---

## 1. The 35-term progression with the smallest last term

```
a    = 73415383
d    = 37418700  = 2^2 * 3 * 5^2 * 11 * 17 * 23 * 29
n    = 35 terms
last = a + 34d = 1345651183
```

All 35 terms:

```
  73415383   110834083   148252783   185671483   223090183   260508883
 297927583   335346283   372764983   410183683   447602383   485021083
 522439783   559858483   597277183   634695883   672114583   709533283
 746951983   784370683   821789383   859208083   896626783   934045483
 971464183  1008882883  1046301583  1083720283  1121138983  1158557683
1195976383  1233395083  1270813783  1308232483  1345651183
```

(This run in fact continues one term further, to `1383069883`, giving 36 terms — but the
35-term prefix is what minimises the last term.)

### How far this is proved

For `n = 35` the bad primes `p` with `2p ≤ 35` are `2, 5, 11, 17`, and for those `p | d` is a
**theorem** (two terms divisible by `p` would each need `p² |`, yet their difference is
`j·p·d` with `j < p`). With the forced factor 3, **every** valid step is therefore a multiple
of `5610 = 3·2·5·11·17`. So scanning `d = 5610·m` over all `m` and all `a` is a *complete*
search, and the following are proofs rather than best-effort statements:

- **No 35-term Loeschian AP has last term ≤ 2·10⁸.** Complete family, `m = 1..1049`, every `a`
  (longest runs in that range top out at 27–29 terms).
- **No 35-term AP with last term < 1345651183 has `23·29 | d/5610`** — the sub-family where the
  two remaining sub-35 bad primes are also killed, which is ~130× more likely to produce a hit
  than a generic step. Exhausted (`m = 667j`, `j = 1..10`).
- Sweeps of the `m = 23j` and `m = 29j` sub-families and of the full `m` range were still
  running at the deadline; `experiments/2026-10-03-g1/*.jsonl` records exactly how far each got.

An independently calibrated model (fitted to brute force for `n = 12..24` and validated against
the published record staircase) predicts the true optimum for `n = 35` at last term `≈1.3·10⁹`.
Our `1.345·10⁹` sits right on that estimate, so this is likely optimal or very close — but we
only *proved* the bounds above.

---

## 2. The longest progression: n = 43

```
a    = 5413537078288507
d    = 48916598396160  = 2^8 * 3 * 5 * 11 * 17 * 23 * 29 * 41 * 47 * 53   (= 128 * D0)
n    = 43 terms
last = a + 42d = 7468034210927227      (16 digits)
```

All 43 terms:

```
5413537078288507  5462453676684667  5511370275080827
5560286873476987  5609203471873147  5658120070269307
5707036668665467  5755953267061627  5804869865457787
5853786463853947  5902703062250107  5951619660646267
6000536259042427  6049452857438587  6098369455834747
6147286054230907  6196202652627067  6245119251023227
6294035849419387  6342952447815547  6391869046211707
6440785644607867  6489702243004027  6538618841400187
6587535439796347  6636452038192507  6685368636588667
6734285234984827  6783201833380987  6832118431777147
6881035030173307  6929951628569467  6978868226965627
7027784825361787  7076701423757947  7125618022154107
7174534620550267  7223451218946427  7272367817342587
7321284415738747  7370201014134907  7419117612531067
7468034210927227
```

Regenerate: `python3 -c "a,d=5413537078288507,48916598396160; print([a+k*d for k in range(43)])"`

---

## 3. Verification

Every claim above passes **two independent checks**.

```
$ python3 src/verify.py 73415383 37418700 35
OVERALL: PASS  (n=35, a=73415383, d=37418700)

$ python3 src/verify.py 5413537078288507 48916598396160 43 --maximal
  maximality: a-d = 5364620479892347 is NOT Loeschian
  maximality: a+43*d = 7516950809323387 is NOT Loeschian
  run is maximal in both directions
OVERALL: PASS  (n=43, a=5413537078288507, d=48916598396160)

$ python3 src/crosscheck.py 73415383 37418700 35
CROSS-CHECK: PASS  (35 terms, each with a verified x^2+x*y+y^2 representation)

$ python3 src/crosscheck.py 5413537078288507 48916598396160 43
CROSS-CHECK: PASS  (43 terms, each with a verified x^2+x*y+y^2 representation)
```

- `src/verify.py` factors every term from scratch in Python ints (its own Miller–Rabin and
  Pollard–Brent; no sieve, no fixed-width arithmetic) and checks that every prime `p ≡ 2 (mod 3)`
  occurs to an even power.
- `src/crosscheck.py` does **not use that criterion at all**. It *constructs* integers `x, y`
  with `x² + xy + y² = t` for every term — via Eisenstein-integer arithmetic, Cornacchia for
  primes `≡ 1 (mod 3)` — and checks the identity by direct multiplication. An explicit
  representation is positive proof, independent of any theory about bad primes.
- `make check` runs 21 tests in which three further implementations agree on every `n ≤ 20000`.

---

## 4. Why these steps look the way they do

A Loeschian number is **never `≡ 2 (mod 3)`**: for `t` coprime to 3,
`t ≡ (−1)^(number of bad prime factors, with multiplicity) (mod 3)`, and "Loeschian" forces
every one of those exponents to be even, hence that count even. So `3 | d` is forced and every
term sits in one class mod 3 — measured, with `3 ∤ d` the longest run found anywhere is **2**.
This single constraint lifts the per-term probability from ~0.14 to ~0.6.

Then each bad prime `p` with `2p ≤ n` must divide `d` (the theorem above), and each with
`n/2 < p < n` is not forced but costs a factor `~1/(p+1)` if omitted. That is exactly the shape
of both steps: `d` is 3 times the product of the small bad primes.

Finally, `a` is chosen so that none of 59, 71, 83, 89, 101, 107, 113 divides *any* term: each of
those can hit at most one index in the window, and `a`'s residue pushes the hit outside it.

---

## 5. Why n ≥ 58 was not reached

Not a shortfall of effort or of algorithm — the target is out of range for this hardware:

- The engine (`src/c/apsearch.c`) is a Wróblewski-style two-stage search: the admissible
  residues of `a` are enumerated **directly by nested additive loops**, so an inadmissible
  candidate is never constructed (~3·10⁴), and a bitmask stage then settles 64 candidates per
  modulo over 148 further primes. Measured: **7.2·10¹⁴ candidate `(a,d)` pairs per second** on
  8 cores. An independent analysis confirmed this configuration is at a *provable* floor — the
  minimum number of residues any such scheme must enumerate is exactly `∏(q−58) = 1167541375`,
  which is what it enumerates.
- A calibrated model puts the **smallest last term that can carry 58 terms at ~4.5·10¹⁸**, while
  this search covers terms up to `7.3·10¹⁵`. A 58-term progression almost certainly does not
  *exist* in the region searched, and reaching `~10¹⁸` costs about `10³` times more compute than
  12 hours on 8 cores — roughly 7500 eight-core-hours. The same model reproduces the published
  record staircase (48, 50, 51, 55, 57), which is the main reason to trust it.
- Two tempting shortcuts were **proved impossible**, not merely unlikely: no constant `c`
  can ever extend a progression or repair a near miss (scaling XORs every term's bad-parity
  vector identically, and in a primitive AP no two terms can share a nonempty vector), and no
  meet-in-the-middle exists (for fixed `d` every constraint is a congruence on the single
  unknown `a`). See PROGRESS.md for the full list of what was ruled out and why.

## 6. Where everything is

- `records.json` — machine-readable records and history.
- `PROGRESS.md` — dated log: what was tried, found, and ruled out, with the cost model.
- `experiments/2026-10-02-calibration`, `2026-10-03-ap58*`, `2026-10-03-g1` — campaigns and raw
  per-work-unit output.
- `src/` — the Loeschian core, the authoritative verifier, the constructive cross-check, and
  three searchers (flat sieve, class-compressed, and the two-stage CRT engine).

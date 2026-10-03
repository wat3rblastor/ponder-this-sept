# LOSS — no configuration reaches 30 automatic indices at any size T; the maximum automatic coverage is 17 (best found), 13 was the previous record, and "smallest T with >= 30 automatic" does not exist.

The covering problem is now posed properly and the answer is structural, not a
size question: the per-condition independence that the covering formulation
assumes is false. All simultaneous square conditions on one progression are
consequences of a **single quadratic polynomial**, and the automatic set is
always `{k : P(k) is a perfect square}` for that one polynomial. That set has
at most ~17 elements inside a window of 58, no matter how large the numbers are.

Second result, independent of the first and the reason the old families were so
large: **every bad prime p < 58 must divide d. This is not avoidable by parking
the forced hits on automatic indices** — a bad prime cannot divide an automatic
term without dividing the parameters. So `2*5*11*17*23*29 = 1247290 | d` always,
giving the universal floor `T >= 57*1247290 = 71095530` and, with the pentagonal
shape `d = 24m^2`, the familiar `T >= 5.3e14`.

---

## 1. The covering problem as literally posed (`cover.py`)

Conditions `(c_i, k0_i)`, automatic indices `union {k0_i + c_i j^2} ∩ [0,57]`.
Exact max-coverage (branch and bound over all 773 maximal sets):

| conditions r | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| max coverage | 8 | 16 | 23 | **30** | 36 | 41 | 45 | 49 | 53 | 56 | **58** |

So 4 conditions suffice for 30 and 11 for all 58 — e.g.
`(c,k0) = (1,0),(1,2),(1,4),(1,6)` covers 30. Every set of size 8 is a `c = 1`
translate; no `(c,k0)` covers more than 8.

**This table is useless**, and that is the main finding. The conditions cannot
be chosen independently; see §2.

## 2. Classification of all "free term" mechanisms (the real problem)

A term is Loeschian iff `t = X^2 + 3Y^2` is solvable. For a term to be free —
Loeschian for *every* member of a parametric family — the representation must be
an identity in the parameters. Writing the family's parameters as `(p,q)`:

* A family of integers given by a polynomial in the parameters whose values are
  all of the form `X^2+3Y^2` must itself be `L1^2 + 3 L2^2` for linear forms
  `L1, L2` (degree 2 is forced; rank > 2 is excluded because `L1^2+3L2^2` has
  rank <= 2, and a pencil of psd rank-<=2 forms with no common kernel does not
  occur). So **`a` and `d` are binary quadratic forms in two parameters**, and
  `a + k d` is a pencil of binary quadratic forms.
* `a + kd = L1^2 + 3L2^2 = M^T diag(1,3) M` forces
  `disc(a+kd) = -12 det(M)^2`. Define

      P(k) := -disc(a + k d)/12 ,

  a **quadratic polynomial in k** (`disc` is quadratic in `k`). Then

      automatic set  A  =  { k in [0,57] : P(k) is a perfect square }.

* Leading coefficient `alpha = -disc(d)/12`. Three regimes:
  - **`alpha = 0`** (`d` a perfect-square form): `P` linear → *squares in a
    58-term arithmetic progression*. **This is the whole of the known theory**:
    - `a = x^2, d = 3m^2` → `P(k) = m^4 k` → `A = {j^2}`, 8 indices.
    - `a = 3x^2+m^2, d = 24m^2` → `P(k) = m^2(24k+1)` → `A` = generalized
      pentagonal, 13 indices.
    - mechanism 1 (`3 c d t_{k0}` square) is the same thing: it gives
      `P(k) = c D^2 (k-k0)`, whose square values are exactly `k0 + c j^2`.
    Exhaustive search (`qsearch` with `alpha=0`, `|beta| <= 40000`, all `y <= 4000`):
    **max 13**, attained only by `P = e^2 (24k+1)` and its reversal
    `P = e^2 (24(57-k)+1)` (checked: all 13-hit configs in the log are one of these two). This matches the heuristic
    `r * sqrt(58/beta)` maximized at `beta = 24` (8 square roots of 1 mod 24).
  - **`alpha > 0` non-square `P`** (`d` positive definite): `y^2 = P(k)` is a
    Pell conic. New, unexploited, and slightly better: **15** hits found, e.g.
    `P(k) = 3(k-29)^2 - 143` gives 14, `P(k) = 15(k-26)^2 + 1309` gives 15.
  - **`alpha < 0`** (`d` *indefinite* — allowed: `d` only has to be a positive
    *value* of the form): `y^2 = P(k)` is an **ellipse**, hits are
    representations of one integer by `X^2 + 4|alpha| Y^2`. Best of all three:
    **17** hits, e.g. `P(k) = -(2k-67)^2 + 2210` gives 16 (2210 = 2*5*13*17 is a
    sum of two squares in many ways), and `alpha=-12, beta=672, gamma=784`
    gives 17.
  - **`P` a perfect-square polynomial** (`beta^2 = 4 alpha gamma`): then *all*
    58 indices are automatic — but this happens **only** when `d = lambda a`,
    i.e. `t_k = a*(1+k lambda)`; a rescaling of an actual 58-term solution.
    Proved by substitution: the tangency condition reduces to
    `(C-3A)^2 + 3B^2 = 0`. Circular, no information.

**Why coverage saturates.** In the elliptic case the count is a representation
number of `M` by `X^2+4nY^2`, which grows only like a divisor function, while the
`k`-window forces `|X| <= 114n` and hence `M <= (57n)^2`; growing `M` to buy more
representations spreads them over a wider `X`-range, and the congruence
`X = -beta (mod 2n)` costs another factor `~n/2^omega(n)`. In the Pell case the
solutions in one class grow geometrically (ratio >= 2+sqrt3), so a 58-wide window
holds 2-3 per class and extra classes need a larger discriminant. Both effects
push against each other and the observed maximum is 15-17.

## 3. The frontier: coverage vs minimum term size (`frontier.py`)

Two lower bounds on `T = a + 57d`:
1. `T >= 2 max_{k in A} y_k` (the minimum of a positive definite binary form of
   discriminant `-12 y^2` is `>= 2y`, and `T` is the largest term).
2. `T >= 57 * prod {p : p bad, p < 58, p not parkable}` (§4).

Exhaustive regions searched (each one *complete* for the stated box, i.e. every
quadratic in the box, hence every configuration whose hits all have `y <= YMAX`):

| region | box | max hits |
|---|---|---|
| `alpha = 0` (linear, the known mechanism) | `|beta| <= 40000`, `y <= 4000` | **13** |
| `alpha in [1,24]` (Pell) | `|beta| <= 3000`, `y <= 1500` | **15** |
| `alpha in [-60,-1]` (elliptic) | `|beta| <= 800`, `y <= 800` | **17** |
| vertex-centred `alpha in [1,2000]` | `y <= 1500` | **15** |
| vertex-centred `alpha in [-2500,-1]` | `y <= 1500` | **12** (large &#124;alpha&#124; is worse) |

Nothing above 17 anywhere; the maximum is attained at small `|alpha|` and the
hit count *decreases* as `|alpha|` grows, which is the behaviour the heuristic in
§2 predicts.

| coverage | T lower bound | 2*max y | 57*prod(forced p) | example P (alpha, beta, gamma) | primes forced into d |
|---|---|---|---|---|---|
| 10 | 7.11e7 | 438 | 71095530 | (-32, 664, 44521) | 2,5,11,17,23,29 |
| 11 | 7.11e7 | 282 | 71095530 | (-35, 560, 17956) | 2,5,11,17,23,29 |
| 12 | 7.11e7 | 80 | 71095530 | (2, -114, 1600) | 2,5,11,17,23,29 |
| 13 | 2.91e9 | 362 | 2914916730 | (-20, 540, 29241) | +41 |
| 14 | 2.91e9 | 212 | 2914916730 | (-15, 750, 2401) | +41 |
| 15 | 3.77e9 | 188 | 3768063090 | (-16, 688, 1444) | +53 |
| 16 | 3.77e9 | 188 | 3768063090 | (-16, 752, 4) | +53 |
| 17 | 7.26e12 | 200 | 7261057574430 | (-12, 672, 784) | +41,47,53 |
| 18..58 | **no configuration found at any size** | | | | |

These are *lower* bounds: they ignore whether the pencil `a+kd` is integrally
realizable with `d` divisible by the forced primes. The one fully realized
family below has `T = 5.3e14` at coverage 13, five orders of magnitude above its
own bound of 7.1e7 — so the true frontier is much higher than the table.
Note the frontier is monotone the wrong way for us: higher coverage forces
*more* primes into `d`, because parking 41/47/53 needs `p | y_k`, which fails
once `A` is large.

## 4. Bad primes are not avoidable (`badprimes.py`, and the obstruction)

Posed question: can a bad prime `p < 58` with `p` not dividing `d` have its one
forced residue class land inside the automatic set, where the term is Loeschian
regardless? **No.** If `t_k = X^2 + 3Y^2` and `p = 2 (mod 3)` divides `t_k`, then
`-3` is a non-residue mod `p`, so `p | X` and `p | Y`. In the family, `X, Y` are
the linear forms `L1, L2` with `det M_k = y_k`; two independent congruences
`L1 = L2 = 0 (mod p)` have a nonzero solution only when `p | y_k`. Hence:

> `p` may be omitted from `d` **iff** some class `r mod p` has every one of its
> indices in `A` **and** `p | y_k` for each of them.

Checked directly in the pentagonal family (`y_k = m*s_k`, `s_k^2 = 24k+1`):
`p | s_k` needs `s_k = p`, i.e. `k = (p^2-1)/24`, which is `5, 12, 22, 35` for
`p = 11, 17, 23, 29` — and in every case the rest of the residue class
(`{5,16,27,38,49}`, `{12,29,46}`, `{22,45}`, `{6,35}`) is not automatic. Naive
droppability (ignoring the `p | y_k` condition) wrongly suggests 23, 29, 41, 47,
53 can be dropped, which would have cut `T` from 5.3e14 to 1.2e9; the
corresponding congruence `3x^2 + (s m)^2 = 0 (mod p)` is unsatisfiable exactly
because `-3` is a non-residue. Verified: `family.py`'s `setup()` returns no
solution for any `m`, for all five primes.

Consequences: `1247290 = 2*5*11*17*23*29` divides `d` in every configuration
(2 needs 29 automatic indices in one class, 5 needs 12, 11 needs 5 — all beyond
reach), so `T > 7.1e7` universally, and with the pentagonal shape
`d = 24 m^2, 623645 | m`, `T >= 1369 * 623645^2 = 5.32e14`.

## 5. The best realized configuration, verified

Pentagonal family, 13 automatic indices, `A = {0,1,2,5,7,12,15,22,26,35,40,51,57}`:

    m = 623645 = 5*11*17*23*29,  d = 24 m^2 = 9334394064600,  a = 3x^2 + m^2

200 members (x even, coprime to 5,11,17,23,29) were tested: **2600/2600 of the
automatic terms are Loeschian, 0 failures**; a further 300 members gave a
measured pass rate 0.4742 on the 45 non-automatic indices. Best member found:
`a = 388933884793, d = 9334394064600`, 32 of the 45 non-automatic terms
Loeschian, `T = 532449395566993` (15 digits). Repo verifiers on the automatic
prefix `k = 0,1,2`:

    python3 src/verify.py     388933884793 9334394064600 3   -> OVERALL: PASS
    python3 src/crosscheck.py 388933884793 9334394064600 3   -> CROSS-CHECK: PASS

## 6. Why this cannot win, numerically

Measured here on 300 members of the realized family (13500 terms): the 45
non-automatic terms pass at **0.4742** each, so `0.4742^45 = 2.6e-15` per member.
The family has `9.8e11` members below `T = 1e20` (`x` even, `w <= 433`), so the
expected number of 58-term solutions below `1e20` in this family is
`9.8e11 * 2.6e-15 = 2.5e-3`. Member count grows linearly in `T`, so break-even
is `T ~ 4e23` and you would have to sieve `~1e15` members to get there — worse
than the unconstrained search, which has `~T^2/D0` candidates instead of `~T`.

At coverage 17 (the best structurally available) the remaining 41 terms give
`0.47^41 = 3.6e-14`, a 14x improvement per candidate, but the forced-prime
product jumps to `2*5*11*17*23*29*41*47*53 = 1.27e11`, so `d >= 1.27e11` and the
family below any fixed `T` shrinks by far more than 14x. The automatic-terms
idea therefore cannot be made to pay at `n = 58`.

## 7. Open items (what would overturn this)

* The exhaustive boxes above are complete but finite. The structural argument in
  §2 says the hit count should fall off outside them, and both scans confirm the
  fall-off, but there is no *proof* that 18+ hits is impossible for a quadratic.
  A proof would need a bound on the number of divisors of `N` in an interval of
  length `~AT` around `sqrt(N)` (square-`alpha` case) or on representation
  numbers of `M` by `X^2+4nY^2` with `X` in one class mod `2n` (elliptic case).
* The Pell (`alpha > 0` non-square) and elliptic (`alpha < 0`) mechanisms give
  14-17 automatic indices and have **never been used**; this work did not
  construct an integral `(a,d)` pencil realizing one (only the discriminant
  condition `disc(a+kd) = -12P(k)`, which is necessary, was imposed). If one is
  realizable with `1247290 | d`, it would be the first family beating 13 — but
  §6 shows it still cannot pay for itself at `n = 58`.

## Files

* `cover.py` — exact max-coverage for the covering problem as posed (§1).
* `qsearch.c` — exhaustive search over quadratics `P`, any sign of `alpha`.
  `./qsearch ALO AHI BMAX YMAX REPORT`.
* `sym.c` / `sym2.c` — deep scan of the vertex-centred case
  `P(k) = alpha(2k-h)^2 + C`. `./sym2 ALO AHI YMAX REPORT`.
* `frontier.py` — turns a search log into the coverage/T frontier, including the
  bad-prime parkability test. `python3 frontier.py exh_*.txt`
* `badprimes.py` — which bad primes a given automatic set could free (naive
  version; §4 explains why the `p | y_k` condition kills all of them).
* `family.py` — the realized pentagonal family, verification of the automatic
  indices, pass-rate measurement.
* `exh_*.txt`, `sym_*.txt`, `frontier.txt` — search logs / frontier output.

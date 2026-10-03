# LOSS — the best family constructible from automatic terms is **F17/B** below. It needs **1.4e18 members tested** to yield one 58-term progression, at term size **T ~ 2.5e28**. The laptop threshold was 1e13 members, so this misses by 5 orders of magnitude; at equal term size and equal work per candidate it captures **6e-15** of the unconstrained search's yield at `T = 1e18` and **2e-23** at its own break-even.

Three results came out of closing this gap. Two are permanent.

1. **THEOREM 1 kills the `alpha > 0` ("Pell conic") mechanism outright.** A pencil
   with a positive definite member forces `disc_k(S) = beta^2 - 4*alpha*gamma >= 0`.
   **All 32 of the previously reported 15-hit `alpha > 0` configurations violate it.**
   Sweeping every configuration with >= 12 discriminant-level hits and constructing
   the integral pencil explicitly, the maximum `alpha > 0` coverage actually
   realizable is **4**, not 15 (`sweep.txt`). That lead never existed.
2. **The `alpha < 0` ("elliptic") 17-hit mechanism is real and is constructed here** —
   explicit forms, explicit `L1/L2` certificates, **494 700 automatic terms
   verified with zero failures** (plus 4.8M more in the tail runs, where the
   all-automatic-pass assertion never fired). It is the first automatic family in this project beating 13,
   and it is **3.2e7 times better than the classical 13-hit pentagonal family**
   (whose corrected cost is 4.5e25 members at `T ~ 4.5e33`). It still loses by 15
   orders of magnitude.
3. **A costing error that applies to the old pentagonal numbers too.** The
   per-member success probability is **not** `rho^41` when a bad prime `p < 58` is
   left out of `d`. The primes 41, 47, 53 cannot be parked on automatic indices
   (their forced residue class always contains a non-automatic index), so each
   member carries ~3.3 indices whose pass rate is 0.0115, not 0.45, and the
   correct per-member probability is `rc^41 * E[prod_p (1/p)^m_p] = rc^41 * 1.0e-4`
   — a **1.0e4x** correction (`poison.log`, `poisondist.log`; the tail test in
   `tail.log` sees 0 members with >= 33/41 passes where the naive binomial predicts
   11.7). The expectation must be taken over the distribution of the forced
   classes, which is strongly non-uniform; using the mean exponent is 14x too
   pessimistic and assuming uniform classes is 74x too optimistic. The published
   pentagonal figure `0.4742^45 = 2.6e-15` per member misses this entirely.

---

## 1. Correct framing (one normalisation fix)

`t` is Loeschian iff `t = x^2 + x y + y^2`. A term of a parametric family is
automatic iff that representation is a polynomial identity, i.e. iff the binary
quadratic form `a + k d` (in the family parameters `p, q`) satisfies

    a + k d = Q o M_k  for an integer matrix M_k,   Q := x^2 + x y + y^2,

equivalently `a+kd = L1^2 + L1 L2 + L2^2` for integral linear forms `L1, L2`.
Since `disc(Q) = -3`, the governing quadratic is

    S(k) := -disc(a + k d)/3 = (4(A+kA')(C+kC') - (B+kB')^2)/3,
    a = (A,B,C),  d = (A',B',C'),   alpha = -disc(d)/3,  gamma = -disc(a)/3.

`experiments/2026-10-03-covering` used `P = -disc/12`, the norm form `x^2+3y^2` of
the index-2 suborder; that is `S/4` and silently assumes `det M_k` even. Its
*search* ran over all integer `(alpha,beta,gamma)`, so the hit counts carry over
verbatim; only the translation to `(a,d)` changes (`disc(d) = -3 alpha`,
`disc(a) = -3 gamma`). Everything below uses `S`.

`S(k)` being a perfect square is **necessary but not sufficient**: `a+kd` must
actually lie in the `Q o M` image. `pencil.py:as_QoM` decides that exactly (it
enumerates `u^2+uv+v^2 = A` and `= C` and matches the middle coefficient), so
every "hit" reported here is a verified polynomial identity, not a discriminant
coincidence.

## 2. THEOREM 1 (realizability) — the Pell mechanism is dead

> Let `(a,d)` be a real pencil of binary quadratic forms such that `a + k0 d` is
> positive definite for some `k0`, and let `S(k) = -disc(a+kd)/3`. Then
> `disc_k(S) = beta^2 - 4 alpha gamma >= 0`, with equality iff `d = lambda a`.

*Proof.* Shift so `k0 = 0`. Let `S1, S2` be the symmetric matrices of `a, d`; then
`S(k) = (4/3) det(S1 + k S2)`. Since `S1 > 0`, write `S1 = L L^T`, so
`det(S1+kS2) = det(S1) det(I + k L^{-1} S2 L^{-T})`. The matrix `L^{-1} S2 L^{-T}`
is real symmetric, hence has real eigenvalues `mu1, mu2`, so

    S(k) = (4/3) det(S1) (1 + k mu1)(1 + k mu2)

**factors over R**. A real quadratic that factors over R has non-negative
discriminant, zero iff `mu1 = mu2`, i.e. `S2 = mu S1`. `disc_k(S)` is invariant
under the shift. ∎

At any automatic index `a+kd = L1^2+L1L2+L2^2` with `det M_k != 0`, which is
positive definite — so **every** automatic family satisfies Theorem 1.

Consequences (`pell_obstruction.txt`, `sweep.txt`):

| regime | disc-level max (previous claim) | configs with >= 12 hits swept | **realizable max** |
|---|---|---|---|
| `alpha = 0` (linear `S`; the classical theory) | 13 | 186 | **13** |
| `alpha > 0` (Pell conic) | 15 | 1537 | **4** |
| `alpha < 0` (ellipse) | 17 | 2906 | **17** |

All 32 reported 15-hit `alpha > 0` configurations have
`beta^2 - 4 alpha gamma = -13728 < 0`. The 14-hit ones that *do* pass Theorem 1
(e.g. `alpha=21, beta=-630, gamma=3025`; `alpha=24, beta=-1368, gamma=18769`)
admit no integral pencil either (`pell.log`). The best `alpha > 0` pencil found
anywhere realizes 4 automatic indices: `a=(85,-40,220), d=(-3,0,-6)`,
`A = {8,18,24,26}`. For `alpha < 0` with `gamma > 0`, Theorem 1 is automatic,
which is why the elliptic branch survives intact.

## 3. F17 — the constructed 17-automatic pencil

Target `S(k) = -3 k^2 + 174 k + 25` (17 hits; `alpha = -3`, so `disc(d) = 9`).
`python3 pencil.py -3 174 25` enumerates **all** integral pencils with `a` reduced
positive definite — 24 of them, and in every one all 17 discriminant-level hits are
genuine `Q o M` identities. The cleanest:

    a(p,q) = 3 p^2 - 3 p q + 7 q^2         disc = -75   (positive definite)
    d(p,q) = 18 p^2 + 3 p q = 3 p (6p+q)   disc = +9    (indefinite, SQUARE disc)

    t_k = a + k d,  automatic at the 17 indices
    A = {0, 1, 2, 6, 8, 12, 17, 22, 25, 33, 36, 41, 46, 50, 52, 56, 57}

`family.py` prints the certificate for every index and asserts it coefficient by
coefficient; e.g.

    k= 0 : t_0  = (-2p+q)^2   + (-2p+q)(p-3q)    + (p-3q)^2        |det|= 5
    k=12 : t_12 = (-17p-q)^2  + (-17p-q)(7p-2q)  + (7p-2q)^2       |det|=41
    k=57 : t_57 = (-37p-3q)^2 + (-37p-3q)(17p+q) + (17p+q)^2       |det|=14

**Why `disc(d) = 9` is the whole point.** Every bad prime `p <= 53` must divide `d`
or else pay for its forced class (`parkability.txt` rechecks the parking argument
against this `A`: no residue class of 2, 5, 11, 17, 23, 29, 41, 47, 53 lies inside
`A` *and* has `p | y_k` throughout, because the `y_k` are
5,14,19,31,35,41,46,49,50 and the two indices with `41 | y_k` sit in classes that
leave `A`). In the pentagonal family `d = 24 m^2` is a **square** form, so
`M | d` requires `623645 | m` and costs `M^2`: `d_min = 24*623645^2 = 9.33e12`.
Here `disc(d)` is a perfect square, so `d` **splits into two integral linear
forms**, and divisibility is one linear congruence per prime:

    N = M1 * M2,   p = M1 P,   6p + q = M2 R    (P, R free integers)
    d = 3 M1 M2 P R = 3 N P R                ->  d_min = 3N
    T(P,R) = a + 57 d = 273 M1^2 P^2 + 84 M1 M2 P R + 7 M2^2 R^2
                        (positive definite, disc = -588 N^2)

so the family is **2-dimensional** and `d_min` is linear in `N`, not quadratic.
This is the one transferable lesson of the whole experiment.

Two realisations:

**Variant A**, `N = M = 2*5*11*17*23*29 = 1247290` (the minimum forced set).
`d_min = 3741870`; smallest member (`M1 = 493 = 17*29`, `M2 = 2530`, `P=R=1`)

    a = 2644447,  d = 3741870,  T = a+57d = 215931037       (9 digits)

against the pentagonal family's floor of `5.32e14`: **six orders of magnitude
lower**, and the covering experiment's claimed floor of `7.26e12` for coverage 17
is not binding, because 41/47/53 need not be in `d`.

**Variant B**, `N = M' = 2*5*11*17*23*29*41*47*53 = 127386974990`.
`d_min = 3M' = 382160924970`; smallest member (`M1 = 144478`, `M2 = 881705`)

    a = 57731773177,  d = 382160924970,  T = 21840904496467   (14 digits)

Variant B costs 1.3e4 in member density and buys 1.0e4 in per-member probability
(§5): the same expected yield at a given `T` to within 1.3x, but reached with
1.3e4 **fewer member tests**, so **B is the family to quote**.

### Verification (repo verifiers)

* `verify_auto.txt` — 8 members x 17 automatic indices = **136 terms, each passed
  by BOTH `src/verify.py` and `src/crosscheck.py`** (the latter constructs an
  explicit `(x,y)` with `x^2+xy+y^2 = t`, independent of the bad-prime criterion).
  0 failures.
* `python3 src/verify.py 2644447 3741870 3` -> `OVERALL: PASS`,
  `python3 src/crosscheck.py 2644447 3741870 3` -> `CROSS-CHECK: PASS`
  (indices 0,1,2 are automatic and contiguous, so `verify.py`'s AP interface
  applies directly).
* `rates.log` (314 500) + `variantB.log` (153 000) + `measure.py` (27 200) — in
  total **494 700 automatic terms over 29 100 random members spanning
  `T = 2.2e8 .. 1e18`, zero failures**; `tail.log` and `tailB.log` assert the same
  property on a further 280 000 members (4.8M terms) and never fire (verdicts from `build/apsearch --isl`,
  which agreed with `src/verify.py` on every spot check).

## 4. Member count — exact, not estimated (`count.py`, `variantB.log`)

A member is `(p,q)` up to sign with `N | d`, `d != 0`, largest term `<= T`. The
`d < 0` branch is a valid member too: read the AP backwards, automatic set
`57 - A`, also 17 indices. Because `N` is squarefree, `N | p(6p+q)` holds iff
`N' | (6p+q)` where `N'` is the product of the primes of `N` not dividing `p` —
one AP in `q` per `p` — so the count is **exact** in `O(sqrt T)`.

| T | variant A (exact) | A/T | variant B (exact) | B/T |
|---|---|---|---|---|
| 1e9  | 149 | 1.49e-7 | 0 | — |
| 1e11 | 45 869 | 4.59e-7 | 0 | — |
| 1e13 | 6 358 466 | 6.36e-7 | 0 | — |
| 1e14 | 67 402 801 | 6.74e-7 | 823 | 8.23e-12 |
| 1e15 | — | — | 18 216 | 1.82e-11 |
| 1e16 | — | — | 269 099 | 2.69e-11 |

Asymptotically `members(T) -> c T` with
`c = [2 pi / sqrt(588)] * [prod_{l|N} (2l-1)/l^2] * (wedge fractions)`:
`c_A = 7.1e-7`, `c_B = 5.38e-11`, both matched by the exact counts. For comparison
the pentagonal family has `9.8e-9 T` — **variant A has 72x more members and a
2.5e6 times lower floor.** `2^omega(N)/N` is the maximum density any binary
quadratic `d` can give (it has at most two linear factors), so this is optimal.

## 5. Pass rate, and the poisoning correction

Measured with `build/apsearch --isl` on random members (`rates.log`,
`poison.log`, `variantB.log`):

| T | 1e10 | 1e12 | 1e14 | 1e16 | 1e18 |
|---|---|---|---|---|---|
| variant A, average over the 41 non-automatic indices | 0.5565 | 0.5027 | 0.4647 | 0.4412 | 0.4130 |
| variant A, **clean** indices only | — | 0.5486 | — | 0.4760 | 0.4495 |
| variant A, **poisoned** indices (41/47/53 forced class) | — | 0.0119 | — | 0.0123 | 0.0104 |
| variant B (no poisoned index) | — | — | 0.5052 | 0.4736 | 0.4485 |

In variant A every member has `3.29` non-automatic indices (measured, 400 000
members) in the forced residue class of 41, 47 or 53, where the term needs
`p^2 | t` and the pass rate is `~rc/p`. The per-member probability is therefore

    P_A = rc^41 * E[ prod_{p in {41,47,53}} (1/p)^{m_p} ] = rc^41 * 9.97e-5

where `m_p` is the number of non-automatic indices in `p`'s forced class. **The
expectation must be taken over the distribution of `m`, not at its mean** — and
that distribution is strongly non-uniform, because `k0 = -a/d mod p` is a function
of `(p,q) mod p`, not a uniform residue (`poisondist.log`: `m = 0` occurs in
`2.3e-5` of members, versus `5.6e-3` under the uniform model). Using the mean
exponent (`rc^37.8 * 0.0115^3.3`) is **14x too pessimistic**; assuming uniform
forced classes is **74x too optimistic**. Independent confirmation: `tail.log`
(200 000 members, `T<=1e12`) has observed/binomial 0.8-1.26 up to 32/41 passes and
**0 observed versus 11.7 expected at >= 33/41** — members with `m >= 3` cannot
exceed 38 passes at all.

Variant B has 41, 47, 53 in `d` and no poisoned index, so `P_B = rho_B^41` is the
right shape. Its tail is **not** exactly binomial, though (`tailB.log`, 80 000
members at `T <= 1e15`, `rho = 0.48906`, 0 automatic failures):

| passes | 28/41 | 29 | 30 | 31 | 32 | 33 | 34+ |
|---|---|---|---|---|---|---|---|
| observed | 485 | 174 | 215 | 100 | 30 | 2 | 0 |
| binomial | 457 | 196 | 75 | 25.5 | 7.6 | 2.0 | 0.55 |
| ratio | 1.06 | 0.89 | **2.86** | **3.92** | **3.93** | 1.00 | — |

The 3-4x excess at 30-32 is heterogeneity, not correlation: members differ in
whether the bad primes `59..~200` happen to drop a forced index into the window,
and the sub-population with none of them passes at a higher rate. Heterogeneity
*raises* the upper tail (the same Jensen effect as in variant A, with the opposite
sign because here the favourable sub-population is common), so `rho_B^41` is if
anything **pessimistic** at the 58-term end. Offsetting that, the per-prime
argument for `p = 59..200` makes `rho_B^41` optimistic by an estimated factor ~3.
Net: treat `P_B` as carrying ~1 order of magnitude of uncertainty either way —
which does not touch a verdict that misses by 5 orders.

The same correction applies to the pentagonal family: its published
`0.4742^45 = 2.6e-15` per member is optimistic by ~1e4 for the poisoning and by a
further ~6e3 because 0.4742 was measured at `T = 5e14` and then applied at
`T = 1e20`.

## 6. The verdict table (`model.log`)

| T | A: members | A: P(member) | A: E[sols] | B: members | B: P(member) | B: E[sols] |
|---|---|---|---|---|---|---|
| 1e14 | 7.1e7  | 9.0e-17 | 6.4e-9 | 5.4e3  | 6.8e-13 | 3.7e-9 |
| 1e18 | 7.1e11 | 5.7e-19 | 4.1e-7 | 5.4e7  | 5.2e-15 | 2.8e-7 |
| 1e22 | 7.1e15 | 1.0e-20 | 7.1e-5 | 5.4e11 | 1.1e-16 | 5.6e-5 |
| 1e26 | 7.1e19 | 3.4e-22 | 2.4e-2 | 5.4e15 | 4.1e-18 | 2.2e-2 |
| 1e30 | 7.1e23 | 1.9e-23 | 1.4e1  | 5.4e19 | 2.5e-19 | 1.4e1  |

    variant A : break-even T* = 2.4e28,  members to test = 1.7e22
    variant B : break-even T* = 2.5e28,  members to test = 1.4e18     <-- THE NUMBER
    pentagonal 13-automatic, corrected : T* = 4.5e33, members = 4.5e25

A and B have **the same expected number of solutions below a given T** to within
1.3x — A trades 1.3e4 more members for 1.0e4 lower probability each, which is
nearly exactly a wash (this is the H1 "pinning is neutral" result reappearing).
They differ 1.3e4-fold in **work**: B reaches the same yield with 1.3e4 fewer
member tests, so B is the family to quote.

**Corrected coverage -> cost frontier** (this is the artefact asked for; the
previous table's `T` lower bounds were never realized and its 17-row floor of
7.26e12 assumed 41/47/53 in `d`, which is not required):

| family | automatic | realizable? | `d_min` | smallest member `T` | members(T) | break-even `T*` | members to test |
|---|---|---|---|---|---|---|---|
| `a=x^2, d=3m^2` | 8 | yes (classical) | — | — | ~T | — | hopeless |
| pentagonal `a=3x^2+m^2, d=24m^2` | 13 | yes (classical) | 9.33e12 | 5.32e14 | 9.8e-9 T | 4.5e33 | 4.5e25 |
| `alpha>0` Pell (claimed 15) | **4** | **NO — Theorem 1** | — | — | — | — | nonexistent |
| **F17/A** (this work) | **17** | **yes** | 3.74e6 | **2.16e8** | 7.1e-7 T | 2.4e28 | 1.7e22 |
| **F17/B** (this work) | **17** | **yes** | 3.82e11 | 2.18e13 | 5.38e-11 T | **2.5e28** | **1.4e18** |
| unconstrained, `M \| d` | 0 | — | 1.25e6 | 7.1e7 | `T^2/(114M)` | ~1e14 | — |

**Why it loses, in one line.** At equal term size and equal work per candidate,

    E_family / E_unconstrained = [c_B T rho^41] / [(T^2/(114 M)) rho^58]
                               = 2.4e-2 * rho(T)^-17 / T

The 17 automatic indices buy `rho^-17 ~ 3e7`, but the family is one parameter
poorer: `T` is a *quadratic* form in `(p,q)`, so the family has `~T` members while
unconstrained `(a,d)` has `~T^2` candidates. The ratio is already `8e-12` at
`T = 1e14` and `2e-23` at `T = 1e28`. Nothing is recoverable: the member density
is already maximal for a binary quadratic `d`, and 17 is the maximum coverage.

## 7. Pushing the 18+ impossibility (`scan.c`)

The old scan looped over `beta` with `|beta| <= 800`, which **cannot reach the
vertex-centred configurations for `n >~ 15`** (those have `beta ~ 57 n`), so its
box was not what it looked like. `scan.c` enumerates in the right coordinates.
Completing the square with `X_k = beta - 2 n k`, `alpha = -n`:

    X_k^2 + 4 n y_k^2 = Delta,     Delta := beta^2 + 4 n gamma = disc_k(S)

so a configuration **is** a pair `(n, Delta)` plus a choice of 58 consecutive terms
of the AP `X = beta - 2nk`; the hits are exactly the representations of `Delta` by
`X^2 + 4 n y^2` whose `X` falls in that window. Enumerating `(n, Delta)` directly
is complete and ~1e4 times faster per configuration than looping over `beta`.

    ./scan  300 400000000 18   ->  BEST=17   (scan_n300_D4e8.txt)
    ./scan 1200 400000000 18   ->  BEST=17   (scan_n1200_D4e8.txt)
    ./scan   60 1600000000 18  ->  still running at hand-off; see scan_n60_D16e8.txt
                                   (small n, 4x the Delta range)

`1 <= n <= 1200`, `Delta <= 4e8` contains every configuration all of whose hits
satisfy `X_k^2 + 4 n y_k^2 <= 4e8`; at the 17-hit optimum (`n = 12`) that allows
`y_k` up to 2887 versus the old limit of 800, and `|beta|` up to 2e4. **Max is 17**,
and both 17-hit configurations (`n=3, beta=174, gamma=25` and
`n=12, beta=672, gamma=784`) are reproduced.

**Why 17 is the ceiling, quantitatively.** A hit is a lattice point on the ellipse
`y^2 + n Z^2 = V` (`Z = k - k0`, `V = Delta/(4n)`) with `Z` confined to 58
consecutive integers. The total number of representations of `V` by `y^2 + n Z^2`
is `O(V^eps)` (bounded by a divisor function), while the window keeps only the
fraction `~ (114/2) / (sqrt(V/n)) * (2/pi)` of them, the representations being
equidistributed in angle. So

    E[#hits]  ~  r_{y^2+nZ^2}(V) * 36 sqrt(n/V),

which **decreases** in `V`: buying more representations costs more than it gains.
Eighteen hits would need `r(V) >~ 0.5 sqrt(V/n)`, impossible for large `V` because
`r(V) = O(V^eps)`; the scan then covers the finite remaining range. This is a
proof of finiteness *modulo* equidistribution (a theorem for fixed `n` and
`V -> infinity` along suitable sequences, but not uniformly in `n`), so an
unconditional uniform bound remains open.

**But it is now worth nothing.** §6 shows that even 25 automatic indices would not
change the verdict: the gain is `rho^-8 ~ 2e3` per member while the deficit against
unconstrained search at `T = 1e28` is `2e-23`. **Automatic terms are closed as a
direction for `n = 58`**, not because 18 is impossible but because the mechanism is
one dimension short no matter how many indices it frees.

## 8. What to carry forward

* **Delete the "Pell (`alpha>0`) gives 15, never used" lead** — Theorem 1 kills it,
  realizable max 4.
* **17 automatic indices are realizable and F17 is explicit and verified**, six
  orders of magnitude cheaper than the classical family and 3.2e7x better in final
  cost — and still 1e18 members from a laptop. Close the direction.
* **Fix the pentagonal family's published cost**: `0.4742^45` is wrong twice over.
  The per-member probability must treat the forced classes of 41, 47, 53 separately
  (pass rate 0.0115, not 0.47, and the expectation taken over the *distribution* of
  how many non-automatic indices they poison) — a ~1e4 correction that applies to
  any family leaving a bad prime `p < 58` out of `d` — and the rate must be
  evaluated at the `T` where it is used, not at `5e14`.
* **Transferable:** choose `d` with **square discriminant**. Then `d` splits into
  linear forms and forcing a squarefree `N` into `d` costs `N`, not `N^2`. That is
  worth 2.5e6 in term size here and would be worth the same in any future family
  that needs forced primes in `d`.

## Files

* `pencil.py` — enumerates ALL integral pencils realizing a given `S`, by reducing
  the problem to integer points on a finite `(A',B')` ellipse; `as_QoM` is the
  exact `Q o M` representability test. `python3 pencil.py <alpha> <beta> <gamma>`
* `family.py` — F17: the forms, the 17 `L1/L2` certificates, the parametrisation.
* `verify_auto.txt` — 136 automatic terms x both repo verifiers.
* `measure.py`, `rates.py`/`rates.log`, `variantB.log` — pass-rate measurements.
* `poison.log`, `poisondist.log`, `tail.log`, `tailB.py`/`tailB.log` — the
  41/47/53 poisoning correction (including the non-uniformity of the forced
  classes) and the independence/tail checks.
* `parkability.txt` — which bad primes must divide `d`, for this `A`.
* `count.py`/`count.log` — exact member counts.
* `model.py`/`model.log` — the cost table and the break-even solve.
* `sweep.py`/`sweep.txt`, `pell_obstruction.txt`, `pell.log` — Theorem 1 sweep.
* `scan.c`, `scan_n*.txt` — the extended 18+ scan in `(n, Delta)` coordinates.

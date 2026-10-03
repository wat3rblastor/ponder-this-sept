# LOSS — the best family constructible from automatic terms is F17 below; it needs **3.6e18 members tested** to yield one 58-term progression, at term size **T ~ 5e24**. That is 1e5 times past the laptop threshold (1e13) and 6e-16 of the unconstrained search's yield at equal term size and equal work per candidate.

Two results came out of closing this gap, and both are permanent:

1. **The alpha > 0 ("Pell conic") mechanism does not exist.** Theorem 1 below is a
   realizability obstruction that no discriminant-level search can see: a pencil
   with a positive definite member forces `disc_k(S) = beta^2 - 4*alpha*gamma >= 0`.
   **All 32 of the previously reported 15-hit alpha > 0 configurations violate it.**
   Sweeping every configuration with >= 12 discriminant-level hits and building the
   integral pencil explicitly, the maximum alpha > 0 coverage actually realizable is
   **4**, not 15 (`sweep.txt`).
2. **The alpha < 0 ("elliptic") 17-hit mechanism is real and is constructed here**
   (`F17`, below), the first automatic family in this project beating 13. It is
   **~3300x better than the classical 13-hit pentagonal family** and still loses by
   16 orders of magnitude. Its advantage is not the extra 4 automatic indices
   (worth 14x); it is that `disc(d) = 9` is a **perfect square**, so `d` splits into
   two integral linear forms and forcing the mandatory bad primes into `d` costs one
   congruence per prime instead of a square condition: `d_min = 3741870` instead of
   the pentagonal family's `24*623645^2 = 9.33e12`, and the family is 2-dimensional.

---

## 1. Correct framing (one normalisation fix)

`t` is Loeschian iff `t = x^2 + x y + y^2`. A term of a parametric family is
automatic iff that representation is a polynomial identity, i.e. iff the binary
quadratic form `a + k d` (in the family parameters `p, q`) satisfies

    a + k d = Q o M_k  for an integer matrix M_k,   Q := x^2 + x y + y^2,

equivalently `a+kd = L1^2 + L1 L2 + L2^2` for integral linear forms `L1, L2`.
Since `disc(Q) = -3`, the governing quadratic is

    S(k) := -disc(a + k d)/3 = (4(A+kA')(C+kC') - (B+kB')^2)/3,
    a = (A,B,C),  d = (A',B',C'),   alpha = -disc(d)/3, gamma = -disc(a)/3.

`experiments/2026-10-03-covering` used `P = -disc/12`, the norm form `x^2+3y^2`
of the index-2 suborder; that is `S/4` and silently assumes `det M_k` even. Its
*search* was over all integer `(alpha,beta,gamma)` so the hit counts carry over
verbatim; only the translation to `(a,d)` changes (`disc(d) = -3*alpha`,
`disc(a) = -3*gamma`). Everything below uses `S`.

`S(k)` a perfect square is **necessary but not sufficient**: `a+kd` must actually
lie in the `Q o M` image. `pencil.py:as_QoM` decides this exactly (it enumerates
`u^2+uv+v^2 = A` and `= C` and matches the middle coefficient), so every "hit"
reported here is a verified polynomial identity, not a discriminant coincidence.

## 2. THEOREM 1 (realizability). Why the Pell mechanism is dead.

> Let `(a,d)` be a real pencil of binary quadratic forms such that `a + k0 d` is
> positive definite for some `k0`, and let `S(k) = -disc(a+kd)/3`. Then
> `disc_k(S) = beta^2 - 4 alpha gamma >= 0`, with equality iff `d = lambda a`.

*Proof.* Shift so `k0 = 0`. Let `S1, S2` be the symmetric matrices of `a, d`;
`S(k) = (4/3) det(S1 + k S2)`. `S1 > 0`, so `S1 = L L^T` and
`det(S1+kS2) = det(S1) det(I + k L^{-1} S2 L^{-T})`. The matrix
`L^{-1} S2 L^{-T}` is real symmetric, hence has real eigenvalues `mu1, mu2`, so
`S(k) = (4/3) det(S1) (1+k mu1)(1+k mu2)` **factors over R**. A real quadratic
that factors over R has non-negative discriminant; it is zero iff `mu1 = mu2`,
i.e. `S2 = mu S1`. (The discriminant of `S` is invariant under the shift.) ∎

At any automatic index `a+kd = L1^2+L1L2+L2^2` with `det M_k != 0` is positive
definite, so **every** family of automatic terms satisfies Theorem 1.

Consequences, all checked in `pell_obstruction.txt` and `sweep.txt`:

| regime | disc-level max (prev. claim) | configs with >=12 hits | **realizable max** |
|---|---|---|---|
| `alpha = 0` (linear `S`; the classical theory) | 13 | 186 | **13** |
| `alpha > 0` (Pell conic) | 15 | 1537 | **4** |
| `alpha < 0` (ellipse) | 17 | 2906 | **17** |

All 32 reported 15-hit `alpha > 0` configurations have
`beta^2 - 4 alpha gamma = -13728 < 0`. The 14-hit ones that do pass Theorem 1
(e.g. `alpha=21, beta=-630, gamma=3025`) still admit no integral pencil. The best
`alpha > 0` pencil found anywhere is 4 automatic indices
(`a=(85,-40,220), d=(-3,0,-6)`, `A = {8,18,24,26}`). **The Pell mechanism is not a
lead; it never existed.**

For `alpha < 0` with `gamma > 0`, Theorem 1 is automatic, which is why the
elliptic branch survives intact.

## 3. F17 — the constructed family (this is the new object)

Target `S(k) = -3 k^2 + 174 k + 25` (17 hits; `alpha = -3` so `disc(d) = 9`).
`pencil.py -3 174 25` enumerates **all** integral pencils with `a` reduced
positive definite: 24 of them, and in every one all 17 discriminant-level hits
are genuine. The cleanest:

    a(p,q) = 3 p^2 - 3 p q + 7 q^2        disc = -75   (positive definite)
    d(p,q) = 18 p^2 + 3 p q = 3 p (6p+q)  disc = +9    (indefinite, SQUARE disc)

    t_k = a + k d,  automatic at the 17 indices
    A = {0, 1, 2, 6, 8, 12, 17, 22, 25, 33, 36, 41, 46, 50, 52, 56, 57}

`family.py` prints the certificate for every index; e.g.

    k= 0 : t_0  = (-2p+q)^2  + (-2p+q)(p-3q)   + (p-3q)^2
    k=12 : t_12 = (-17p-q)^2 + (-17p-q)(7p-2q) + (7p-2q)^2
    k=57 : t_57 = (-37p-3q)^2+ (-37p-3q)(17p+q)+ (17p+q)^2

each verified as a coefficient-by-coefficient polynomial identity (assertions in
`family.py:pencil_certificates`).

**Why `disc(d) = 9` is the whole point.** Every bad prime `p` with `2p <= 58`
must divide `d` (`2*5*11*17*23*29 = 1247290 = M`; the parking argument in
`2026-10-03-covering` §4 is unaffected by anything here — rechecked against this
`A`: no residue class of 2,5,11,17,23,29 lies inside `A`, and no `y_k` is
divisible by 41, 47 or 53 in a class contained in `A`). In the pentagonal family
`d = 24 m^2` is a *square* form, so `M | d` needs `623645 | m` and costs `M^2`.
Here `disc(d)` is a perfect square, so `d` **splits**, and `M | 3p(6p+q)` is
satisfied by one linear congruence per prime:

    M = M1 * M2,   p = M1 P,   6p + q = M2 R    (P, R free integers)
    d = 3 M1 M2 P R = 3 M P R           ->   d_min = 3M = 3741870
    T(P,R) = a + 57 d = 273 M1^2 P^2 + 84 M1 M2 P R + 7 M2^2 R^2

The smallest member (`M1 = 493 = 17*29`, `M2 = 2530`, `P = R = 1`):

    a = 2644447,  d = 3741870,  T = a + 57 d = 215931037   (9 digits)

against the pentagonal family's floor `T >= 5.32e14`. **Six orders of magnitude
lower, and the family is 2-dimensional in `(p,q)` instead of 1-dimensional.**

### Verification (repo verifiers only)

* `verify_auto.txt`: 8 members x 17 automatic indices = **136 terms, each passed
  by BOTH `src/verify.py` and `src/crosscheck.py`** (the latter constructs an
  explicit `(x,y)` with `x^2+xy+y^2 = t`). 0 failures.
* `python3 src/verify.py 2644447 3741870 3` -> `OVERALL: PASS`;
  `python3 src/crosscheck.py 2644447 3741870 3` -> `CROSS-CHECK: PASS`
  (indices 0,1,2 are automatic and contiguous).
* `rates.log`: 18 500 random members across `T ~ 1e10 .. 1e18`,
  **314 500 automatic terms, 0 failures** (verdicts from `build/apsearch --isl`).
* `measure.py` on the 1600 smallest members (`T` from 2.16e8): 27 200 automatic
  terms, 0 failures; longest overall run observed 18.

## 4. Cost (`count.py`, `rates.log`, `model.py`)

**Member count — exact lattice count, not an estimate.** A member is `(p,q)` up
to sign with `M | d`, `d != 0` and largest term `<= T`; the `d < 0` branch is a
valid member too (read the AP backwards, automatic set `57 - A`). Because `M` is
squarefree, `M | p(6p+q)` is one AP in `q` for each `p`, so the count is exact in
`O(sqrt T)`:

| T | members (exact) | members/T |
|---|---|---|
| 1e9  | 149 | 1.49e-7 |
| 1e11 | 45 869 | 4.59e-7 |
| 1e13 | 6 358 466 | 6.36e-7 |
| 1e14 | 67 402 801 | 6.74e-7 |

Asymptotically `members(T) -> c T` with
`c = [2 pi / sqrt(588)] * [prod_{l|M} (2l-1)/l^2] * (wedge fractions) = 7.1e-7`,
matched by the exact counts. Compare the pentagonal family: `9.8e-9 T`.

**Pass rate on the 41 non-automatic indices (measured, `rates.log`):**

| T | 1e10 | 1e12 | 1e14 | 1e16 | 1e18 |
|---|---|---|---|---|---|
| rho | 0.5565 | 0.5027 | 0.4647 | 0.4412 | 0.4130 |
| rho^41 | 3.67e-11 | 5.67e-13 | 2.27e-14 | 2.68e-15 | 1.80e-16 |

(Fit `rho(T) = 2.6405 (ln T)^-0.4981`, residuals < 0.8%. These are lower than the
0.4742 quoted for the pentagonal family because that figure was measured at
`T = 5e14` and because only `M` — not 41, 47, 53 — is forced into `d` here.)

**Expected 58-term solutions and members to test** (`model.log`):

| T | members | rho^41 | E[solutions] | members to test for one |
|---|---|---|---|---|
| 1e14 | 7.1e7  | 3.0e-14 | 2.2e-6 | 3.3e13 |
| 1e18 | 7.1e11 | 1.8e-16 | 1.3e-4 | 5.6e15 |
| 1e20 | 7.1e13 | 2.1e-17 | 1.5e-3 | 4.8e16 |
| 1e22 | 7.1e15 | 3.0e-18 | 2.1e-2 | 3.4e17 |
| **5.0e24** | **3.6e18** | **2.8e-19** | **1.0** | **3.6e18** |

**Break-even `T* = 5.0e24`, 3.6e18 members tested.** The threshold for a win was
1e13 members with a solution in reach: **missed by 5 orders of magnitude.**

**Corrected comparison table** (the useful artefact; the covering experiment's
figures for the pentagonal family used a pass rate measured at `T = 5e14` and
applied it at `T = 1e20`, which is optimistic by `(0.4742/0.3919)^45 = 5.9e3`):

| family | auto | `d_min` | `T` floor | members(T) | E[sols] below 1e20 | break-even T | members to test |
|---|---|---|---|---|---|---|---|
| pentagonal `a=3x^2+m^2, d=24m^2` (known) | 13 | 9.33e12 | 5.32e14 | 9.8e-9 T | 4.2e-7 (corrected) | ~1e28 | ~1e22 |
| **F17 (this work)** | **17** | **3.74e6** | **2.16e8** | **7.1e-7 T** | **1.5e-3** | **5.0e24** | **3.6e18** |
| alpha>0 Pell (claimed 15) | **4** | — | — | — | — | — | nonexistent |
| unconstrained, `M \| d` | 0 | 1.2e6 | 7.1e7 | `T^2/(114 M)` | 1.8e8 | ~3e13 | — |

**Why it loses, in one line.** At equal term size and equal work per candidate,

    E_family / E_unconstrained = [c T rho^41] / [(T^2/(114 M)) rho^58]
                               = 101 * rho(T)^-17 / T

The automatic indices buy `rho^-17 ~ 6e7`, but the family is one parameter
poorer: `T` is a *quadratic* form in `(p,q)` so the family has `~T` members,
while unconstrained `(a,d)` has `~T^2` candidates. The ratio is `1e-5` already at
`T = 1e12` and `5.5e-16` at `T = 1e25`. Nothing recoverable: the `2^6/M` density
of `M | p(6p+q)` is already the maximum a binary quadratic `d` can achieve (it
has at most 2 linear factors), and `17` is the maximum coverage.

## 5. Pushing the 18+ impossibility (`scan.c`)

The old scan looped over `beta`, and with `|beta| <= 800` it **could not reach the
vertex-centred configurations for `n >~ 15`** (those have `beta ~ 57 n`), so its
box was not what it looked like. `scan.c` enumerates in the right coordinates.
Completing the square with `X_k = beta - 2 n k`, `alpha = -n`:

    X_k^2 + 4 n y_k^2 = Delta,    Delta := beta^2 + 4 n gamma = disc_k(S)

so a configuration **is** a pair `(n, Delta)` plus a choice of 58 consecutive
terms of the AP `X = beta - 2nk`; the hits are the representations of `Delta` by
`X^2 + 4n y^2` whose `X` falls in that window. Enumerating `(n, Delta)` directly
is complete and 1e4 times faster per configuration.

    ./scan 300 400000000 18   ->   BEST=17     (no 18 anywhere in the box)
    ./scan 1200 400000000 18  ->   see scan_n1200_D4e8.txt

Box `1 <= n <= 300`, `Delta <= 4e8` contains every configuration all of whose
hits satisfy `X_k^2 + 4 n y_k^2 <= 4e8`; for the 17-hit optimum (`n=12`) that
allows `y_k` up to 2887 versus the old limit of 800, and it covers `|beta|` up to
2e4. **Max is 17.** Both 17-hit families (`n=3` and `n=12`) are reproduced.

**Why 17 is the ceiling, quantitatively.** A hit is a lattice point on the
ellipse `y^2 + n Z^2 = V` (`Z = k - k0`, `V = Delta/(4n)`) with `Z` confined to a
window of 58 consecutive integers. The total number of representations of `V` by
`y^2 + n Z^2` is `O(V^eps)` (it is bounded by a divisor function), while the
window keeps only the fraction `~ 114/(2 sqrt(V/n)) * (2/pi)` of them — the
representations are equidistributed in angle. So

    E[#hits] ~ r_{y^2+nZ^2}(V) * 36 sqrt(n/V),

which **decreases** in `V`: buying more representations costs more than it gains.
Getting 18 hits needs `r(V) >~ 0.5 sqrt(V/n)`, which fails for all large `V`
because `r(V) = O(V^eps)`. This is a proof of finiteness modulo
equidistribution (which is a theorem for `n` fixed and `V -> infinity` along
sequences with bounded class number, but not uniformly); the exhaustive scan then
covers the finite remaining range. **An unconditional uniform bound is still
open** — but it is now worth nothing, because §4 shows 18, or even 25, automatic
indices would not change the verdict: at 25 automatic indices the per-member gain
is `rho^-8 ~ 1e3` and break-even would still be at `T ~ 1e21` with `~1e15`
members, and the `rho^-17/T` argument still beats it.

## 6. What to carry forward

* **Delete the "Pell (alpha>0) gives 15, never used" lead.** Theorem 1 kills it;
  realizable max is 4.
* **The alpha<0 mechanism is real and 17 is realizable** — F17 is explicit,
  verified, and six orders of magnitude cheaper than the classical family. It is
  still 1e16 short. Automatic terms are closed as a direction for `n = 58`.
* The one transferable lesson: **square `disc(d)` makes `d` split, and forced
  prime divisibility then costs `M` instead of `M^2`.** If any future family
  needs forced primes in `d`, choose `d` with square discriminant.
* `src/verify.py` / `src/crosscheck.py` agreed on all 136 spot-checked automatic
  terms; `build/apsearch --isl` agreed with them and found 0 failures in 314 500
  automatic terms. No numeric claim here is unverified except the `rho`
  extrapolation beyond `T = 1e18` (explicitly a fit, shown with its residuals).

## Files

* `pencil.py` — enumerates ALL integral pencils realizing a given `S`, via the
  finite `(A',B')` ellipse; `as_QoM` is the exact representability test.
  `python3 pencil.py <alpha> <beta> <gamma>`
* `family.py` — F17: the forms, the 17 `L1/L2` certificates, the parametrisation.
* `verify_auto.txt` — 136 automatic terms x both repo verifiers.
* `measure.py`, `rates.py` / `rates.log` — pass-rate measurements.
* `count.py` / `count.log` — exact member counts.
* `model.py` / `model.log` — the cost table and the break-even solve.
* `sweep.py` / `sweep.txt` — realizable maximum per regime (the Theorem 1 sweep).
* `pell_obstruction.txt` — Theorem 1 applied to the 15-hit alpha>0 log.
* `scan.c` / `scan_n*.txt` — the extended 18+ scan in `(n, Delta)` coordinates.

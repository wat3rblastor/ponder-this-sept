# VERDICT: **CONFIRMED — the drift onset does NOT occur before `s = 0.707`.** The transport ratio `res(s~0.70)/res(s~0.44)`, measured at matched endpoints in four purpose-built families with 48 to 2.0e9 exact solutions per point, is **0.944 .. 1.045 (geometric mean 0.998)**: the `s = 0.44 -> 0.707` extrapolation that `groundtruth` flagged as its honest gap is worth **1.00 +- 0.05**, not the factor of several its `base = 330` series implied. **Residual at `s = 0.707` in the `base(58)` family = 0.958** (its own measured 0.96 at `s <= 0.44`, transported by 0.998 — an EXTRAPOLATION, not a measurement; band 0.88 .. 1.04). **Corrected n=58 even-odds cost: 1.64e5 core-hours, band 1.43e5 .. 1.85e5 — UNCHANGED** against the current 1.5e5-1.8e5, and inside the current 0.7e5-3.6e5 band, which this narrows.

**Scope, stated up front, because two very different claims are easy to confuse here.** I did **NOT** reach `s = 0.707` in the `base(58) = 3741870` family itself. That family's own exact reach is **unchanged at `s = 0.44` (n = 36, X = 2.5e9, `groundtruth`)**; my own `base(58)` run was at a deliberately *smaller* `X = 1e9` and reached only n = 34, `s = 0.41`. What I reached at `s = 0.707` and beyond (to `s = 0.955`) is in **smaller-`B` families, where large `s` is cheap** — chosen so that their structure matches `base(58)`'s at n=58 (several absent bad primes, the smallest of them multi-hit). Only a direct `base(58)` measurement at `s = 0.707` would close the literal gap; **that is still open, and section 2 shows it is not closable on this hardware at any `X`.** What is closed is the question the gap was *about* — whether the residual drifts between `s = 0.44` and `s = 0.707` — and the answer is no.

**The brief's route to n=42 does not exist.** The largest `X` I ran was **1e9** (measured peak RSS **0.61 GB**), not the requested ~1e10 / ~6.5 GB, because this machine's 8.2 GB swap was **already 6.6 GB full** from other agents before I started and never fell below that (load average 29-120 throughout). I began at `X = 2e9` (1.38 GB) and backed off to 1e9 when swap reached 7.4/8.2 GB. **No segment-size reduction was needed and this experiment caused no swapping.** The reduction cost nothing scientifically, because the requested computation could not have worked at any `X` that fits here: the expected number of 42-term APs in the `base(58)` family is **0.023 at X = 1e10** and reaches 1 only near **X = 1.6e11** (40 GB leanest, 110 GB with `exact.c` as written). `groundtruth`'s "X ~ 1e10 reaches n = 42 / ~6.5 GB / ~1 hour" was the condition that a legal `d` exists (`mmax >= 1`), not that a 42-term AP exists. Section 2 gives the arithmetic.

**Nothing of length >= 50 was found and nothing could have been.** The largest `n` with a nonzero exact count anywhere in this experiment is **n = 36**, at last terms below 1e9. No progression claim is made in this directory, and `src/verify.py` / `src/crosscheck.py` were not invoked on any new progression because there was none to invoke them on.

Laptop only (M2, 8 cores, 8 GB), 2026-10-03, ~2 h wall. Every run `nice -n 18/19`, `OMP_NUM_THREADS=3`. Per-core rates in section 6.

---

## 0. Conventions, stated before any number is used

`exact.c` reports three counts per `(n, X, B)`. `COUNT` and `MAXIMAL` differ by up to 4.4x here
and the repo has been corrupted by confusing them, so:

| column | exactly what it counts |
|---|---|
| **`COUNT`** | **PAIRS** `(a,d)`, `d = B*m` with `m >= 1`, `a >= 1`, `a = 1 (mod 3)`, `a` coprime to every bad prime dividing `B`, `a+(n-1)d <= X`, and `a+kd` Loeschian for all `0 <= k < n`. So: **APs of length AT LEAST `n`, counted once per starting pair** — an AP of length `L >= n` contributes `L-n+1`. **ALL APs, not maximal ones.** |
| `RUNSTART` | the same restricted to pairs with `a-d` not Loeschian (or `< 1`): **runs of length AT LEAST `n`, each once**. Not maximal (may extend right). |
| `MAXIMAL` | runs of length `>= n` extendable in neither direction inside `[1,X]`: **maximal runs of length AT LEAST `n`**. |

**Every model/exact ratio in this file is against `COUNT`** — length **at least**
`n`, **all** APs — because an expected-count model integrates over all `(a,d)`
with no left-end condition. **Nothing here is a count of APs of length exactly
`n`, and nothing here is a count of maximal APs.** `counts.txt` tabulates
`COUNT/RUNSTART` and `COUNT/MAXIMAL` per point so any number here can be
converted; `COUNT/RUNSTART` runs 1.00 .. 2.09 and `COUNT/MAXIMAL` 1.00 .. 4.37 across these families.

`s = n / (2 q_min)`, with `q_min` the smallest bad prime **not** dividing `d`.
`n > 2 q_min` is empty by the forced-divisibility theorem, so `s <= 1` always
and `s -> 1` is a hard wall. n=58 in the forced family has `q_min = 41`, hence
`s = 58/82 = 0.707`.

## 1. The counter is `groundtruth`'s, and it was re-verified before use

`exact.c` here is `experiments/2026-10-04-groundtruth/exact.c`, copied, with one
3-line change: `--dens` now honours `--base`, so class densities can be measured
for families other than the five hard-coded ones. **The counting path is
untouched.**

| check | result |
|---|---|
| the four published counts, `X = 2e6`, `n = 12..18` | **434143 / 95725 / 16975 / 2851** — exact |
| the four published counts, `X = 1.9e7`, `n = 22..28` | **1819 / 285 / 48 / 5** — exact |
| the PROVED `n = 35` minimum, `X = 311958331` | **exactly 1**, `a = 219830911, d = 2709630` = `records.json g1_min_last_term_35` |
| **NEW — the designed families brute-forced against `src/loeschian.is_loeschian`** (pure Python, exact factorisation; this validates the *family conventions* — which `d`, which `a` — not just the bitmap) | `B=748374, n=9, X=1e7`: **737 / 737** both. `B=162690, n=24, X=6e6`: **0 / 0** both. `B=340170, n=15, X=8e6`: **104 / 53** both. **0 disagreements** (`xcheck_newfam.txt`) |
| measured class density `rho(1e9)` for `B = 3741870` | **0.6157**, against `2026-10-04-rhoexp`'s independent `min6` ladder **0.6135** |
| the residual at fixed `s`, re-measured at `X = 2e9` | **0.907 / 0.907 / 0.905** against **0.910 / 0.910 / 0.903** at `X = 1e9` — agrees to 0.3% (section 4) |
| `cost58.py` with residual 1.0 | reproduces `groundtruth`'s forecast: `T* = 1.93e18`, 1.74e5 core-hours vs their 1.94e18, 1.77e5 (`cost58_R1.txt`) |

## 2. Why an `n = 42` count in the `base(58)` family is unreachable here

`exact.c` holds 5 bitmaps of `X` bits plus a prime sieve of `X/2` bits =
`0.6875 X` bytes; the leanest possible variant (`L` + one in-place working
buffer, a periodic admissible-`a` pattern) is `0.25 X` bytes. Expected `COUNT`
is the model validated below (`n42_reach.py`):

| X | `mmax` at n=42 | E[n=38] | E[n=40] | **E[n=42]** | `exact.c` GB | leanest GB |
|---|---|---|---|---|---|---|
| 2.5e9 (`groundtruth`'s best) | 16 | 0.254 | 0.029 | **0.003** | 1.7 | 0.6 |
| **1e10 (the brief's target)** | 65 | 1.66 | 0.199 | **0.023** | **6.9** | 2.5 |
| 2e10 | 130 | 4.24 | 0.523 | **0.063** | 13.8 | 5.0 |
| 4e10 | 260 | 10.6 | 1.31 | **0.161** | 27.5 | 10.0 |
| 1e11 | 651 | 34.6 | 4.20 | **0.511** | 68.8 | 25.0 |
| ~1.6e11 | ~1040 | ~55 | ~6.7 | **~0.8** | ~110 | ~40 |

Counts fall ~2x per term of `n` while supply grows only as `X^2`, so each extra
term of `n` costs 1.41x in `X` and 2x in memory: going from `groundtruth`'s
n = 36 to n = 42 costs 8x in `X`, i.e. **64x in memory**. **No `X` that fits in
8 GB gives an expected count anywhere near 1 at `n = 42` in this family** (at
the brief's 6.9 GB it is 0.023, a 2% chance of a single hit); even 32 GB with a
counter that does not exist yet reaches only `E ~ 0.65`. This is a correction to the follow-up `groundtruth`
named for itself, not a shortfall in execution.

## 3. What was done instead: an exact `s`-scan in purpose-built families

A family `d = B*m` is exactly countable whenever `base(n) | B`. Choosing `B` to
contain every forced bad prime **except one**, `q_min`, puts `s = n/(2 q_min)`
anywhere in `(0,1]` at a **moderate** `n`, because `s` is set by `q_min`, not by
`n`. Cost is `~7 X^2 / (64 B (n-1))` word ops — **inverse in `B`** — so these
families are cheap exactly where the forced family is unaffordable. Each has the
same qualitative structure as `base(58)` at n=58: several absent bad primes,
with the smallest one or two of them multi-hit (`q < n`).

All at **`X = 1e9`** (plus 3 points of the `q_min = 17` family re-run at
`X = 2e9`; see the X-independence subsection below). Exact counts: **`counts.txt`** (the standalone raw table).
Per-point model ratios: **`table.txt`**.

| `q_min` | `B` | factorisation | `n` run | `s` range | largest `n` with nonzero `COUNT` | `COUNT` range |
|---|---|---|---|---|---|---|
| 5 | 748374 | `2*3*11*17*23*29` | 4..9 | 0.400-0.900 | **9** | 2.70e7 .. 1.95e9 |
| 11 | 340170 | `2*3*5*17*23*29` | 10..21 | 0.455-0.955 | **21** | 1.74e4 .. 3.75e7 |
| 17 | 220110 | `2*3*5*11*23*29` | 16..33 | 0.471-0.971 | **33** | 4 .. 1.08e6 |
| 23 | 162690 | `2*3*5*11*17*29` | 24..45 | 0.522-0.978 | **36** | 1 .. 5991 |
| 29 | 129030 | `2*3*5*11*17*23` | 24..40 (step 2) | 0.414-0.690 | **34** | 4 .. 21455 |
| **41** | **3741870 = base(58)** | `2*3*5*11*17*23*29` | 24..42 | 0.293-0.512 | **34** | 1 .. 2503 |

### The residual, point by point (M3/exact, `COUNT` in parentheses)

`q_min = 5`, n = 4..9 — the series with overwhelming statistics:

| `s` | 0.400 | 0.500 | 0.600 | **0.700** | 0.800 | 0.900 |
|---|---|---|---|---|---|---|
| M3/exact | 0.977 | 0.962 | 0.972 | **0.983** | 0.995 | 1.007 |
| `COUNT` | 1.95e9 | 6.32e8 | 2.81e8 | 1.28e8 | 5.89e7 | 2.70e7 |

Poisson error here is below 0.02%, so this is not a statistical statement:
**the residual moves by 3.0% from `s = 0.400` to `s = 0.900`, monotonically and
upward.** The model very slightly *under*predicts at large `s`.

`q_min = 11`, n = 10..21 (`s` = 0.455 .. 0.955), counts 3.75e7 down to 1.74e4:
0.956, 0.934, 0.936, 0.942, 0.944, 0.953, 0.956, 0.963, 0.972, 0.980, 0.994, 1.010.

`q_min = 17`, n = 16..33 (`s` = 0.471 .. 0.971), counts 1.08e6 down to 4:
0.935, 0.908, 0.909, 0.910, 0.910, 0.910, 0.903, 0.887, **0.883** (`s=0.706`),
0.873, 0.886, 0.902, 0.938, 0.933, 1.001, 1.068, 1.202, 1.488 — the last four
points have 27, 11, 4 and 4 solutions.

`q_min = 23`, n = 24..36 (`s` = 0.522 .. 0.783), counts 5991 down to 1:
0.855, 0.834, 0.826, 0.811, 0.775, 0.772, 0.795, **0.894** (`s=0.674`, 48
solutions), 0.929, 0.958, 0.881, 0.984, 0.844 — the dip to 0.772 at `s = 0.630`
and the recovery to 0.894 at `s = 0.674` are a 15% wobble on counts falling from
231 to 48, i.e. Poisson, and they bracket rather than exceed the level.

### The `base = 330` "drift onset" is not reproduced at matched `q_min`

`groundtruth` put the drift onset between `s = 0.735` and `s = 0.794` on the
strength of its `base = 330` family (`q_min = 17`), where M3/exact runs 0.786,
0.742, 0.748, 0.746, 0.818, 0.903, 1.163, 1.612, 2.017, 4.738 for
`s = 0.647 .. 0.912` — on counts 19342, 6643, 3082, 1400, 555, 205, **59, 13,
5, 1**. My `q_min = 17` family is the matched control: same `q_min`, same absent
bad primes (17, 23, 29, 41, ...), the same two of them multi-hit at the top of
the range, and 2-12x the solutions at the same `s`:

| `s` | `base=330` M3/exact (`COUNT`) | **`B=220110` M3/exact (`COUNT`)** |
|---|---|---|
| 0.706 | 0.748 (3082) | **0.883 (4615)** |
| 0.735 | 0.746 (1400) | **0.873 (2389)** |
| 0.765 | 0.818 (555) | **0.886 (1189)** |
| 0.794 | 0.903 (205) | **0.902 (586)** |
| 0.824 | **1.163 (59)** | **0.938 (278)** |
| 0.853 | **1.612 (13)** | **0.933 (136)** |
| 0.882 | **2.017 (5)** | **1.001 (61)** |

**The upward excursion lives in the low-count tail and does not survive better
statistics at the same `s`.** The pooled `base = 330` deficit above `s = 0.82`
(78 observed against ~130 expected at a flat 0.80 residual) is a real ~4-sigma
fluctuation in that one series, but it is **not a function of `s`**: at
`s = 0.824` with 278 solutions the residual is 0.938, and at `s = 0.900` with
2.7e7 solutions it is 1.007.

## 4. The decisive number: the within-family transport ratio

`base(58)` can be measured only at `s <= 0.44`; n=58 sits at `s = 0.707`. The
quantity that licenses carrying its own measured residual across that gap is
`res(high s) / res(low s)` **within one family, at matched endpoints**:

| `q_min` | `B` | low anchor | high anchor | res(low) | res(high) | **ratio** | `COUNT`s |
|---|---|---|---|---|---|---|---|
| 5 | 748374 | n=4, `s=0.400` | n=7, `s=0.700` | 0.977 | 0.983 | **1.006** | 1.95e9, 1.28e8 |
| 11 | 340170 | n=10, `s=0.455` | n=16, `s=0.727` | 0.956 | 0.956 | **1.000** | 3.75e7, 5.33e5 |
| 17 | 220110 | n=16, `s=0.471` | n=24, `s=0.706` | 0.935 | 0.883 | **0.944** | 1.08e6, 4615 |
| 23 | 162690 | n=24, `s=0.522` | n=31, `s=0.674` | 0.855 | 0.894 | **1.045** | 5991, 48 |

> **Transport ratio 0.944 .. 1.045, geometric mean 0.998 over four families.**
> The `s = 0.44 -> 0.707` extrapolation is worth **1.00 +- 0.05**.

`base(58)`'s own count-weighted geometric-mean residual over its exact points
with `COUNT >= 30` is **0.942** (`X = 1e9`, 7 pts, `s = 0.29-0.37`), **0.996**
(`X = 1.5e9`, 10 pts, `s = 0.17-0.38`) and **0.935** (`X = 2.5e9`, 7 pts,
`s = 0.32-0.39`) — three independent `X`, mean **0.96**. Transporting:
**res(n=58, `s = 0.707`) = 0.958**, band 0.88 .. 1.04. **The model
`under`predicts by ~1.04x: there are slightly *more* 58s than it says.**

### Slopes, and why the nominal CIs are not the uncertainty

Per-family weighted LS of `ln(res)` on `s`, over each family's own `s` range
(points with `COUNT >= 30`):

| `q_min` | `B` | pts | `s` range | `d ln(res)/ds` | res(hi)/res(lo) |
|---|---|---|---|---|---|
| 5 | 748374 | 6 | 0.400-0.900 | `+0.012 +- 0.037` | 1.031 |
| 11 | 340170 | 12 | 0.455-0.955 | `-0.060 +- 0.051` | 1.056 |
| 17 | 220110 | 18 | 0.471-0.882 | `-0.212 +- 0.043` | 1.071 |
| 23 | 162690 | 8 | 0.522-0.674 | `-0.776 +- 0.174` | 1.046 |
| 29 | 129030 | 10 | 0.345-0.517 | `-0.616 +- 0.048` | 0.829 |
| **41** | **3741870** | 24 | 0.171-0.390 | `-0.169 +- 0.035` | **1.008** |
| 11 | 30 (legacy) | 4 | 0.545-0.818 | `-0.315 +- 0.039` | 0.912 |
| 17 | 330 (legacy) | 10 | 0.647-0.824 | `-0.156 +- 0.450` | 1.479 |

Pooled over all 92 points with `COUNT >= 30`:
`d ln(res)/ds = +0.008 +- 0.009` (all `s`), `-0.027 +- 0.011` (`s <= 0.75`),
`+0.118 +- 0.008` (`s in [0.6,1.0]`), giving `res(0.707) = 0.966 .. 0.984`.
**`chi2/dof` is 270-4000**, because with 1e9-solution points the Poisson error
is ~1e-5 while the real family-to-family scatter is ~10%. **The nominal 95% CIs
(0.95-1.00) are therefore far too tight to quote as the uncertainty.** The
honest uncertainty is the scatter: the per-family slopes run from `-0.78` to
`+0.01` per unit `s`, in **both** signs with no trend in `q_min`, and over the
0.267-wide gap from `s = 0.44` to `s = 0.707` the worst of them is worth
`e^{-0.78*0.267} = 0.81` to `e^{+0.01*0.267} = 1.00`. The four **direct**
transport measurements at matched endpoints (0.944-1.045) are a tighter and
better-posed estimate than any slope fit, and they are what the headline uses.
**Either way the answer to "does drift onset occur before `s = 0.707`" is no:
no family drifts upward across that gap, and the largest downward excursion is
19%.**

### Is `s` the right scaling variable? Partly — and that is handled

Count-weighted geometric-mean residual by `s` window and `q_min`:

| `s` window | `q_min=5` | `q_min=11` | `q_min=17` | `q_min=23` |
|---|---|---|---|---|
| [0.60, 0.70) | 0.972 | 0.943 | 0.889 | 0.783 |
| [0.70, 0.80) | 0.983 | 0.955 | 0.835 | — |
| [0.80, 0.90) | 0.995 | 0.972 | 0.967 | — |
| [0.90, 1.01) | 1.007 | 0.999 | — | — |

There **is** a `q_min` dependence of the absolute level (the residual falls with
`q_min` at fixed `s`), but **no `s` dependence within a fixed `q_min`**. That is
exactly the shape this argument needs: the absolute level for `q_min = 41` is
taken from `base(58)`'s *own* exact points, never extrapolated in `q_min`, and
only the flat `s` behaviour is transported.

### X-independence: the residual is a function of `s`, not of `X`

`run_xcheck.sh` re-measures the `q_min = 17` family at a second, independent
`X = 2e9`. It was **started without my intending it** (another process in this
workspace picked it up from section 8) and I **killed it at 3 of 32 points on
swap pressure** — swap free was falling, 740 -> 568 MB, and this machine is
shared. The three points it did finish carry 2.7-3.0x the solutions of the
`X = 1e9` run at the same `s`, and they settle the question:

| `s` | n | `X = 1e9`: M3/exact (`COUNT`) | **`X = 2e9`: M3/exact (`COUNT`)** |
|---|---|---|---|
| 0.588 | 20 | 0.910 (60512) | **0.907 (179555)** |
| 0.618 | 21 | 0.910 (31867) | **0.907 (93302)** |
| 0.647 | 22 | 0.903 (16848) | **0.905 (48356)** |

**Agreement to 0.3% at a 2x change in `X`.** The residual depends on `s` and on
the family, not on the range swept — which is what has to be true for any of
this to transport to `T ~ 1e18`. The remaining 29 points of that script are the
cheapest outstanding check in this directory (~30 min on an idle machine); they
would tighten `q_min = 17, 23` at `s = 0.7-0.95`, not change the headline.

Folding the `X = 2e9` points in moves no published number here by more than
0.001 (`table_X1e9only.txt` is the pre-fold table, kept for that comparison);
the transport ratio span moves from 0.945-1.046 to 0.944-1.045 and its
geometric mean is unchanged at 0.998.

## 5. HEADLINE: the n=58 cost is unchanged — and the drift could not have reached it anyway

Two questions the repo has mixed:

**(A) SUPPLY — how many 58-term APs exist below `T`.** Family `d = 3741870*m`,
`q_min = 41`, so **this** question sits at `s = 0.707` and is the one exposed to
`s`-drift. With `res = 0.958`:

| T | E[#58, all terms <= T], model | corrected |
|---|---|---|
| 1e16 | 2.06e-2 | 2.15e-2 |
| 1e17 | 3.44e-1 | 3.59e-1 |
| 1e18 | 6.64 | 6.93 |
| 1e19 | 1.42e2 | 1.48e2 |

Even-odds existence moves from `T = 1.75e17` to **`T = 1.69e17`** — a 3% shift.
(Divide by ~2 for distinct runs.) `assumptions` S4 and `groundtruth` 4A stand.

**(B) COST — core-hours for the production engine.** Family `d = K*D0`,
`D0 = 2*3*5*11*17*23*29*41*47*53`, cost `4.69e-32 T^2` core-hours. **The
production sieve `b2 = 2000` guarantees that no bad prime `<= 2000` divides ANY
of the 58 terms** (`src/c/apsearch.c`: tier B picks `a mod q` to push the hit
outside `[0,n-1]`; tier C bitmasks the rest up to `b2`; every bad prime `<= 58`
divides `D0` already, which I checked by factoring `D0`). So in the engine
family **every prime that could produce a hit is `> 2000 >> 58`. Its effective
`s` is `58/4006 = 0.014`, not 0.707, and it has no multi-hit or
class-exhaustion structure at all** — the engine deliberately forgoes exactly
the single-hit-with-`q^2` configurations whose statistics were in question.
**The `s`-drift, whatever its answer had been, cannot touch the engine cost.**
That is a result this experiment adds, and it is the deeper reason the DECIDES
answer is "unchanged".

With the low-`s` residual 0.96:

| T | E[#58] model | corrected | core-hours | **core-h per 58** |
|---|---|---|---|---|
| 1e17 | 1.58e-2 | 1.65e-2 | 4.69e2 | 2.85e4 |
| 1e18 | 2.95e-1 | 3.08e-1 | 4.69e4 | **1.52e5** |
| 1e19 | 6.20 | 6.46 | 4.69e6 | 7.26e5 |

> **Even-odds find: `T* = 1.87e18`, cost = 1.64e5 core-hours.**
> Band from the residual envelope 0.88 .. 1.04: **`T* = 1.75e18 .. 1.99e18`,
> cost 1.43e5 .. 1.85e5 core-hours.**
> Against `groundtruth`'s 1.5e5-1.77e5 and the plan's ~2.0e5: **UNCHANGED.**
> **The term-size band a production run should target is a last term near
> `T* = 1.9e18`** — the band already in the plan. **No reconfiguration is
> implied and none is recommended.**

1.43e5 .. 1.85e5 is **much narrower than the 0.7e5 .. 3.6e5 currently carried**,
because that band came from `groundtruth`'s inter-family slope envelope
`|b| <= 0.015` per term of `n` (a factor 0.4..2.4 at n=58); the `s`-scan
replaces that with a directly measured transport ratio of 0.998 +- 0.05 and
`base(58)`'s own residual measured at three independent `X`.

## 6. Rates and memory, for the next agent

Swap was **already 6.6 GB of 8.2 GB used** before I started (load average
29-120 from other agents) and never fell below it. I started at `X = 2e9`
(1.38 GB), **backed off to `X = 1e9`** when swap hit 7.4/8.2 GB, and measured
peak RSS **0.61 GB** thereafter. **This experiment caused no swapping.** No
segment-size reduction was needed, because I never attempted a large-`X` run:
**the right lever for this question is `B`, not `X`.**

| run | X | B | n | wall per n | RSS |
|---|---|---|---|---|---|
| bitmap build | 1e9 | — | — | ~20 s | 0.61 GB |
| `--dens` (single-threaded bit scan) | 1e9 | any | — | ~150 s each | 0.61 GB |
| `--base 748374` | 1e9 | 748374 | 4..9 | 15-75 s | 0.61 GB |
| `--base 340170` | 1e9 | 340170 | 10..21 | 25-33 s | 0.61 GB |
| `--base 220110` | 1e9 | 220110 | 16..33 | 8-16 s | 0.61 GB |
| `--base 162690` | 1e9 | 162690 | 24..45 | 11-14 s | 0.61 GB |
| `--base 3741870` | 1e9 | 3741870 | 24..42 | 2-20 s | 0.61 GB |
| `--base 220110` (X-check, killed) | 2e9 | 220110 | 20..22 | 105-118 s | 1.08 GB |

At `nice 18` with 3 OpenMP threads under load average 29-76 the counter
delivered **1.3e8 to 1.1e9 word-ops per core-second** — the low end whenever
`COUNT >= 1e8`, where the per-solution bit-enumeration loop (needed for
`RUNSTART` / `MAXIMAL`) dominates the bit-parallel AND chain.

**The cheap direction for `s` is a designed `B`, not a large `X`.** Reaching
`s = 0.9` with 2.7e7 exact solutions cost 15 s and 0.61 GB here; reaching
`s = 0.707` in the forced family would cost ~110 GB.

## 7. Scope — what this does NOT show

- **It does not extend the `base(58)` family's own exact reach.** Still n = 36
  at `X = 2.5e9`, `s = 0.44` (`groundtruth`). My `base(58)` run at `X = 1e9`
  reached n = 34 and is only a consistency check at a third `X` (residual 0.942,
  against 0.996 and 0.935 at the other two). **The literal gap — `base(58)`
  measured at `s = 0.707` — remains open and is not closable on this hardware
  (section 2).**
- The remaining assumption is that the `s`-flatness measured at
  `q_min = 5, 11, 17, 23` also holds at `q_min = 41`. The `q_min` trend in the
  *level* is measured and sidestepped (the level is taken from `base(58)`
  itself). The trend in the *slope* is not separately measurable here, because
  no family with `q_min >= 29` has affordable counts above `s = 0.7`: the
  `q_min = 29` family's counts die at n = 34 (`s = 0.59`) even at `X = 1e9`, and
  pushing it to `s = 0.707` needs `X ~ 1.6e10`. The six measured slopes show no
  trend in `q_min` and straddle zero, which is the evidence available.
- All ratios are for `COUNT`: pairs, length **at least** `n`, **all** APs.
  Nothing here is a maximal-AP or exactly-`n` count.
- Section 5B inherits `rho_s = 1.72409e-12` and `4.69e-32 T^2` core-hours from
  `2026-10-04-smallt` and `2026-10-03-laptopcost` unchanged. This experiment
  tested the probability side against exact truth; it did not re-measure the
  cost law or the window density.
- No statement about which regions are already swept; that is `smallt`'s
  coverage map, and it is what actually decides where to spend.
- `run_xcheck.sh` (a second `X = 2e9` series for `q_min = 17, 23`) ran only
  **3 of its 32 points** before I killed it on swap pressure; those 3 confirm
  `X`-independence to 0.3% (section 4). The other 29 are the cheapest remaining
  check, ~30 minutes on an idle machine, and would tighten `q_min = 17, 23` at
  `s = 0.7-0.95` without changing the headline.

## 8. Exact reproduction

```sh
cd experiments/2026-10-05-n42
gcc-15 -O3 -fopenmp -o exact exact.c -lm

# 1. re-verify the counter (the 8 published counts + the PROVED n=35 minimum)
OMP_NUM_THREADS=3 nice -n 18 ./exact 2000000    12 18 2        # 434143 95725 16975 2851
OMP_NUM_THREADS=3 nice -n 18 ./exact 19000000   22 28 2        # 1819 285 48 5
OMP_NUM_THREADS=3 nice -n 18 ./exact 311958331  35 35 1 --list # exactly 1
python3 xcheck_newfam.py > xcheck_newfam.txt   # brute force vs src/loeschian

# 2. the s-scan: class densities, then exact counts (6 families, X = 1e9, 0.61 GB)
X=1000000000 T=3 N=18 ./run_sscan.sh           # -> dens_*.txt, q5/q11/q17/q23/q29/q41.txt

# 3. the artefacts
python3 table_counts.py > counts.txt           # THE RAW EXACT COUNTS (standalone)
python3 analyse.py      > table.txt            # ratios, residual vs s, slopes, transport
python3 n42_reach.py                           # section 2: why n=42 is unreachable
python3 cost58.py 1.0   1.0  1.0  1.0  > cost58_R1.txt  # baseline; reproduces groundtruth
python3 cost58.py 0.958 0.96 0.88 1.04 > cost58.txt     # THE CORRECTED n=58 COST

# X-independence at a second X (~1.1 GB RSS; WATCH SWAP -- I killed this at
# 3 of 32 points when swap free fell below 600 MB).  q17b.txt holds those 3.
X=2000000000 T=3 N=18 ./run_xcheck.sh
```

`counts.txt` (raw exact counts) and `table.txt` (ratios and fits) are the
standalone artefacts; `cost58.txt` is the headline. Every number quoted above is
in one of them.

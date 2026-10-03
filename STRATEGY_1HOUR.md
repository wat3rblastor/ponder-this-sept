# Strategy: a verified n ≥ 58 within one hour of GPU wall-clock

Derived 2026-10-03. Numbers from the calibrated model in PROGRESS.md, re-derived
independently here, and cross-checked against the measured 0.4 hits/hour on the 8× RTX PRO 6000
box. Reproduce every table with the snippets in §6.

---

## 1. The bottom line

**Rent ~64 GPUs of RTX PRO 6000 Blackwell class for one hour, run the best remaining units
first, and stop on the first `n ≥ 58`.**

| GPUs for 1 h | E[hits] (model) | P(success) | P(success), 2× safety |
|---|---|---|---|
| 8 | 0.48 | 38% | 21% |
| 16 | 0.91 | 60% | 36% |
| 32 | 1.65 | 81% | 56% |
| **64** | **2.91** | **95%** | **77%** |

Cost at ~$1.5–2 per GPU-hour: **$95–130 for one hour.** The 2× safety column exists because the
942,000 units already completed carry `E[hits] = 1.32` by this model and produced no 58 — a 27%
outcome, so the model may be ~2× optimistic. Size for the safety column.

A6000-class cards are ~3–4× slower per GPU on this kernel, so the same confidence needs ~200+ of
them; prefer fewer, faster GPUs.

---

## 2. Why one hour is reachable at all: yield is wildly non-uniform

Expected hits per residue is `64 · ρ(T)^58 · corr`, and `ρ(T) = 2·0.638 / (√(ln T) · F)` falls
with term size. Across the plan the yield per residue falls about **16×**. Running best-first
instead of sweeping is worth **~2.5×** on its own:

- measured average over a long flat run: 0.05 hits per GPU-hour
- best-first (model): 0.124 hits per GPU-hour over the first 7 GPU-hours

**This lever is already implemented and does not need work.** `tools/plan_units.py` ranks units
by a score that is strictly better than the analytic one above: a Monte-Carlo per-term pass rate
`ρ(t) = 0.71992·(ln t/39.144)^-0.485` (1.2e7 samples at each of 9 sizes from 1e15 to 3e18)
applied to the actual terms rather than to `T_last`; an exact window-position offset for the
strided kernel; a correction for good primes `p ≡ 1 mod 3` dividing `K` (whose windows are
rescalings of already-ranked units); and an absolute calibration of windows per unit. It heaps
the units and skips every `(K, shift)` already present in `experiments/*/*.jsonl`.

Consequence for sizing: the measured 0.4 hits/hour on 8 GPUs was **already best-first**, so the
table in §1 is the honest marginal rate, not an optimistic one. Just run it.

---

## 3. The configuration is already optimal — proved, so don't spend the hour tuning

The engine's whole efficiency is *how selective stage 1 is at a given term size*, because
`hits/residue = 64 · P_raw / δ` and `P_raw` is fixed by nature. Term size is floored by
`T ≥ 64·MOD`, so pinning more primes costs size. Sweeping the pin set:

| pinned 59..P | MOD | T floor | 1/δ | hits/residue |
|---|---|---|---|---|
| 59 | 1.8e3 | 2.18e13 | 442 | 3.6e-15 |
| 89 | 9.3e8 | 2.18e13 | 23035 | 1.7e-13 |
| **101** | 9.4e10 | 2.18e13 | 54106 | **4.0e-13** |
| 107 | 1.0e13 | 6.4e14 | 118150 | 4.2e-14 |
| **113 (current)** | 1.13e15 | 7.3e16 | 242745 | 2.0e-15 |
| 131 | 1.5e17 | 9.5e18 | 435612 | 1.1e-16 |

Pinning 59..101 is 200× better per residue — **but only at the floor term size `T = 57·D0`,
which exists solely at `K = 1`, where there is no supply.** Once `d` pushes `T` above `64·MOD`
the two configurations coincide. Measured: a per-K adaptive `MOD` is worth **1.0×**. Keep
`MODCAP = 2e16`, pin 59..113.

Selectivity per unit of modulus is `ln(q/(q−58)) / ln q`: 59 scores **1.00**, 71 scores 0.40, 113
scores 0.15, 131 scores 0.12. That is why the set stops where it does.

### Ruled out — do not spend the hour on these

| idea | verdict |
|---|---|
| per-K adaptive MOD | **1.0×** (computed above) |
| `q²` "square option" at 59 (avoid ∪ 59²\|term) | modulus-inefficient: 1/δ falls 29.75 per 3481 of modulus vs 59 per 59. Running it as a second mode is **4.3× worse per residue** |
| square option at 41/47/53 (out of `d`) | reaches 5.35× more solutions but at **~95× lower** yield per residue at matched `T`; the `41²·47²·53²` stage-1 modulus alone forces `T ≥ 2e13` |
| dropping 59 from the pin set, letting tier C take it | 3.7 vs 0.99 hits per GPU-hour — **keep 59 pinned** |
| pentagonal automatic terms (`a = 3x²+m²`, `d = 24m²`) | measured today: 1400× better per candidate, 6·10⁸× fewer candidates. **Loses by ~10⁵** (`experiments/2026-10-03-pentagon`) |
| more pinned primes via 128-bit MOD | every row past 113 loses; 197 is 3·10⁷× worse |

---

## 4. Run book

```bash
make check                                  # python tests + CPU sieve
tools/build_here.sh                         # must print best=47 for K=205
MODCAP=2e16 python3 tools/plan_units.py --budget-res 2e16 --kmax 600000 --smax 3000 \
    --out experiments/remote_plan.txt       # value-ordered; finished units excluded
tools/multi_gpu.sh experiments/remote_plan.txt <tag>     # PPG=2, MPS on, --kernel 31 --t0 24
setsid nohup tools/autopromote.sh 60 >/dev/null 2>&1 </dev/null &
```

Across several machines, shard the one plan with `--slice i/N` so no unit is searched twice.

`MODCAP` is 2e15 by default and 2e16 in the GOAL.md run book; both select the same pin set
(59..113, `MOD = 1.133e15`), since adding 131 needs 1.48e17. Either is fine.

**Three things that will silently cost you the hour:**

1. **Stage-3 starvation.** Stage 3 runs in a `std::thread` that ignores `--threads` and needs
   **~25 cores per GPU**; `OMP_NUM_THREADS=8` measured **3.8× slower**. A 62-vCPU box with 8 GPUs
   is in exactly that regime. Either get ≥20 cores/GPU, or merge the parked Montgomery stage-3
   patch (2.2× less CPU per exact test — it was shelved as "not needed while GPU-bound", which
   stops being true on a core-starved box).
2. **A new `--tag` per machine**, or coverage logs collide.
3. **Never `pkill -f <pattern>` from a shell whose own command line contains the pattern.**

The engines stop themselves when any log shows a run ≥ 58. Then `tools/promote.py <a> <d> <n>`
(must pass both `src/verify.py` and `src/crosscheck.py`), confirm `records.json` / `ANSWER.md`,
append to PROGRESS.md, commit, push.

---

## 5. Why there is no cheaper route

Four results, each established independently, jointly say the hour has to be bought rather than
out-thought:

1. **No linear function of `k` is identically a norm** — `N(z₀+kw)` is quadratic in `k` with
   leading coefficient `N(w)`, zero only for `w = 0`. No parametrization makes every term a norm.
2. **Scaling cannot extend or repair.** With `π(n) = {q ≡ 2 mod 3 : v_q(n) odd}`, Loeschian means
   `π = ∅` and `π(cn) = π(c) △ π(n)`; scaling XORs every term identically, and in a primitive AP
   no two terms share a bad prime. Probability exactly 0.
3. **CRT can force at most ~2/58 of each term's digits** — `a` pinned mod `∏q²` is itself `~∏q²`.
   The terms are necessarily mostly unstructured.
4. **The automatic-terms identity is real but costs a dimension.** `3cd·t_{k0}` square makes
   `t_{k0±cj²}` free, but it forces `a = s²/(3cd) − k₀d`, dropping the candidate count from
   `~T²/D0` to `~T/√D0`. Break-even needs `T ≲ 1400·√D0 ≈ 4·10⁹`; the family forces
   `T ≥ 5.3·10¹⁴`. Measured, not just argued.

Exhaustive small cases (`n = 8..28`) show optimal APs look like random admissible ones.

---

## 6. Reproducing the numbers

The pin-set table, the cumulative-yield curve and the one-hour table are all produced by the
snippets recorded in this file's git commit message and in
`experiments/2026-10-03-pentagon/README.md`. Constants: `D0 = 382160924970`,
`MOD = 1.133e15` (pin 59..113), `R = 4.67e9` residues per unit, `9.3e10` residues/s per
RTX PRO 6000 Blackwell, `ρ(T) = 2·0.638/(√(ln T)·F)` with `F = ∏_{bad p ≤ 113} p/(p+1)`,
`corr = ∏_{bad q > 113} (1−58/(q+1))/e^{−58/(q+1)}`.

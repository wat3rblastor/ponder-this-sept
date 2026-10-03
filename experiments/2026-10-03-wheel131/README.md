**VERDICT: NOT WORTH DOING. The "free 1.46x from pinning 131" is not obtainable in the production CUDA engine. Upper bound on the gain is 4.9% of the plan's expected 58s (planner model), and that bound assumes three things that are all false. With the measured GPU cost per residue it is 0.38%, with the engine's real window offset it is 0.08%, and with both it is 0.00%. The change is also infeasible as written: pinned MOD' = 1.485e17 >= 2^56, and the largest candidate a = 2.8e19 > 2^64. 69% of production residues already pin 131. No patch.**

Scope: analysis plus a short probe of the unchanged engine. Nothing under src/ or tools/ was changed, and the production engines were not touched.

## 1. How the engine enumerates, and what "pin 131" would mean there

- `d = K*D0`. Stage 1 (`prepare` in `src/c/apsearch_cuda.cu`, around lines 924-950) builds a CRT wheel with modulus MOD. The wheel takes 3 (a=1), 2 (a=1) and 5 (a!=0; 4 residues). It then greedily adds the smallest bad primes q > 58 with q not dividing d (tier B), as long as `MOD*q <= modcap`. Each tier-B prime contributes `q-58` admissible residues.
- `thr` = the product of all component counts except the two largest, which run on-thread (the in1/in2 loops).
- Each residue R is one 64-bit word whose bit b is the candidate `a = R + (b + 64*shift)*MOD`. The kernel cost is per word.
- Tier C covers every other bad prime r <= b2 = 10000, with one AND per word in kill order. **131 is a tier-C row whenever it is not pinned.** So every stage-3 survivor already has 131 dividing none of its 58 terms. The per-candidate pass-probability gain from "131 avoided" is already realised. Pinning instead of filtering changes only how many a-values one word covers.
- Kernel 31 is strided: it walks unreduced `A < MOD*(1+c1+c2+nch)` and needs `MOD < 2^56` (guard at line 1041).

The pin set per MOD class, measured on the production plan `experiments/remote_plan_cost_f.txt` with MODCAP 2e16 (`plan_stats.out`):

| MOD | pins | share of plan residues | units |
|---|---|---|---|
| 1.1337e15 (main) | 59,71,83,89,101,107,113 (131 in tier C) | 30.0% | 1,657,874 |
| 2.5171e15 (59 divides K) | 71..113, **131** | 48.9% | 36,957 |
| 1.47e15-2.09e15 (one of 71..113 divides K) | 59.. minus that prime, **131** | about 6% | about 108k |
| 3.05e15-4.86e15 (59 and another prime divide K) | ..., **131, 137** | about 10% | about 2.8k |

By 131 status: **pinned 68.9%**, tier-C 28.8%, 131 divides K 2.3%. So MODCAP 2e16 already pins 131 for every class where it fits, which is 69% of the work. The remaining question concerns only the main class (28.8%).

## 2. Why the enumfloor "free 1.46x" does not carry over

1. **It does not fit.** enumfloor's wheel modulus `M = prod Q = 4.95e15` left out the factor 30 (2, 3 and 5). In the engine those are in the wheel, so the main class with 131 has `MOD' = 30*59*...*113*131 = 1.485e17`. That exceeds the term floor `X = 64*MOD = 7.26e16`, so by enumfloor's own rule (`prod Q <= X`), 131 is not free.
2. **The GPU lift is a 64-bit word.** Keeping the a-range at 64*MOD with MOD' = 131*MOD leaves 64/131 = 0.49 lifts per residue, i.e. less than one useful bit per word. Keeping 64 lifts per word means each unit covers a < 64*MOD' = 9.5e18. That is 131 old shifts' worth of range for 73 old units' worth of work (73 x 4.67e9 = 3.41e11 residues, which equals the existing 2.5171e15 class's residue count). The plan never takes more than 81 shifts of any K (median 8). No main-class K has a complete block of 131 shifts in the plan, so the pinned unit always covers a-ranges the planner rejected.
3. **The gain figure was the wrong model and is already banked.** `((q+1)/q)^58 = 1.554` is the independent-terms model. The exact enrichment of a 58-AP given that 131 is avoided is `1/(1-58/q+58/(q(q+1))) = 1.784`. That is equal to the cost `q/(q-58) = 1.795`, up to the 0.6% avoid-only loss. Filtering 131 in tier C already gives this enrichment per candidate at one table row per word. So pinning has nothing left to gain per candidate. Its only effect is 1.795x more a-coverage per word, at 131x the term range.
4. **64-bit arithmetic.** MOD' = 1.485e17 >= 2^56 = 7.2e16, so the engine refuses it (probe below). With the engine's default inner components (73, 55), the largest candidate is about `MOD'*(1+72+54+64) = 2.84e19 > 2^64`. Choosing the small components (13, 25) as the inner loops gives 1.60e19 + 57d, which is under 2^64 only for K < about 1.3e5, and it would cut the per-thread work 12x (c1*c2 = 325 vs 4015).

## 3. Numbers: expected 58s per GPU-second (planner model, `tools/plan_units.py` scoring)

`pin131_plan.py` samples 3000 main-class Ks from the production plan (it skips K where 131 divides K):

- the plan's own units for those K give E58/res = 9.75e-17;
- a 131-pinned unit (K, s'=0) gives E58/res = 7.06e-17, which is **0.72x**. This is optimistic: overflow is ignored, the window offset is wbar = 0.5, and GPU cost per residue is assumed equal. The per-K median is 0.75, and the pinned unit is better for 104 of the 3000 K.

`pin131_gain.py` is an upper bound over every main-class K in the plan (191,834 K). For each K, the pinned unit replaces that K's planned units if it beats them net of displaced work valued at the plan's marginal E58/res (CUT):

| assumptions | CUT (E58/res) | K that switch | upper-bound gain in plan E58 |
|---|---|---|---|
| no overflow, wbar 0.5, equal cost | 1.0e-16 | 6,522 | 1.20% |
| same | 8.5e-17 | 18,634 | 2.55% |
| same, CUT = the plan's last units (7.6e-17, measured from the plan tail) | 7.6e-17 | 39,930 | **4.93%** |
| + measured cost/res for res-3.4e11 units (0.76 of reference, PROGRESS.md solo benchmark) | 7.6e-17 | 3,316 | 0.38% |
| + the engine's real window offset for MOD' | 7.6e-17 | 780 | 0.08% |
| both | 7.6e-17 | 28 | 0.00% |

The plan's model E58 is 2.812. In the second check above, the planner itself with MODCAP 1.5e17 (which drops the overflowing main class) gives E58 = 2.958 at 2e16 residues, against **3.113 with MODCAP 2e16**. That is a 5% loss, so raising MODCAP is also not a cheap route to 131.

All the probability numbers above are the planner's model, not measured yields. The cost factor 0.76 is the measured solo-run figure from PROGRESS.md for res-3.4e11 units. A pinned main-class unit has exactly that shape (res 3.409e11, inner components 73 x 55, the same as the 2.5171e15 class).

## 4. Engine probe (unchanged production binary, GPU 7, about 1 s each; `probe.out`)

- `--modcap 2e16`, K=7: MOD=1133661268029390, res 4.67e9, tierC 608.
- `--modcap 2e16`, K=59: MOD=2517112306980510, res 3.409e11 (131 pinned).
- `--modcap 1.5e17`, K=7: `strided kernel needs b2 <= 16383, MOD < 2^56; skipping` (the pinned main class is refused).
- `--modcap 1.5e17`, K=59: unchanged (2.5171e15 x 137 > 1.5e17).

## 5. Reproduce

```
cd /workspace/ponder-this-sept
D=experiments/2026-10-03-wheel131
python3 $D/plan_stats.py                       # pin classes, shift structure, 131 status  -> plan_stats.out
python3 $D/pin131_model.py                     # per-K v(old) vs v(pinned), 64-bit headroom -> pin131_model.out
python3 $D/pin131_plan.py                      # E58/res, plan units vs pinned unit       -> pin131_plan.out
python3 $D/pin131_gain.py experiments/remote_plan_cost_f.txt 7.6e-17     # upper bound (~4 min)
PINCOST=1.316 python3 $D/pin131_gain.py experiments/remote_plan_cost_f.txt 7.6e-17
REALW=1 python3 $D/pin131_gain.py experiments/remote_plan_cost_f.txt 7.6e-17
for mc in 2e16 1.5e17; do MODCAP=$mc python3 tools/plan_units.py --budget-res 2e16 --kmax 600000 \
  --smax 3000 --done-glob /dev/null --old-plane0-kmax 0 --marks 1e15,3e15,1e16,2e16 \
  --out $D/plan_mc$mc.txt 2> $D/plan_mc$mc.log; done     # ~10 min each; plan .txt files were deleted (9.5 MB)
CUDA_VISIBLE_DEVICES=7 ./build/apsearch_cuda --nterms 58 --kmin 7 --kmax 7 --shifts 1 \
  --modcap 150000000000000000 --b2 10000 --threads 4 --kernel 31 --t0 24 --prep 8 --report 55 \
  --out $D/probe.out.txt                       # prints the 2^56 refusal
```

The probe outputs were renamed from `.jsonl` to `.out.txt` so the planner's `experiments/*/*.jsonl` done-glob does not pick them up.

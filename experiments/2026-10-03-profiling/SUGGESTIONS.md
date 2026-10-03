# Live profiling of the v2 GPU search (8x RTX PRO 6000, 2x EPYC 7742)

Raw samples in this directory: `rate_*.txt` (aggregate and per-engine res/s from the v2 jsonl
files, via `rate.sh N`), `smi_*.csv` / `dmon_3.txt` (GPU samples), `cgroup_*.txt` (CFS quota
counters), `thrcpu_*.txt` (per-thread CPU of one engine, via `thrcpu.py`), `bench/` (stage-3
microbenchmark, old vs new `loesch_core.h`).

## Cycle 1 — 2026-10-03 17:43–17:50 UTC (old binary, OpenMP threads not capped)

**Rate: 6.28e11 res/s** (92.7 s window, 2708 units). Per GPU 5.0e10 to 1.04e11.

Where the time went:
- The container has a CFS CPU quota of **245.76 cores** (`/sys/fs/cgroup/cpu.max` = 24575999/100000),
  not 256. Over 60 s it used 245.7 cores and was **throttled in 581 of 600 periods (97%)**:
  the whole container, including the threads that launch kernels, was frozen at the end of
  almost every 100 ms period.
- Load average ~1650. Each engine had **268 threads, 256 of them runnable**: the stage-3 worker
  is a std::thread, so `--threads 16` never reached it and every engine ran a 256-thread
  OpenMP team (16 x 256 = 4096 threads on 246 cores).
- GPUs: mean utilisation 76-94%, 7-23% of 2 s samples below 50%.
- Stage-3 cost (microbenchmark, thread CPU time, under load): old `is_loeschian` 39.2 us per
  call, new (Montgomery, multiply-by-inverse trial division) 7.5 us. About 5 calls per
  candidate (simulated early exit for `--report 55`).
- `ptrace` is not permitted in this container (gdb attach fails), so no stack sampling.

The orchestrator applied the fix (new `loesch_core.h` + `OMPT=8`) and restarted at ~17:50.

## Cycle 2 — 2026-10-03 17:51–17:58 UTC (new binary, OMP_NUM_THREADS=8)

**Rate: 7.63e11 res/s** (180.5 s window, 5558 units); 8.07e11 in the first 60 s window.
Per GPU (180 s): 7.8e10, 9.8e10, 9.5e10, 1.01e11, 1.07e11, 9.7e10, 1.03e11, 8.4e10. The spread
comes from unit-size variance within a 3-minute window, not from any one GPU being slow.

Measured before and after the restart:

| | before | after |
|---|---|---|
| aggregate res/s | 6.28e11 | 7.63e11 (+21.6%) |
| CPU cores used (cgroup) | 245.7 | 51 |
| CFS-throttled periods | 97% | 0% (0 of 627) |
| threads per engine | 268 | 20 |
| GPU SM busy (dmon, 1 s) | 76-94% | 100% in every sample on every GPU |
| GPU power | ~470-520 W | 549 W (cap 550 W) |

Where the time goes now: **the GPUs are the bottleneck, and they are capped by power.**
- `nvidia-smi dmon` gave sm = 100% in all 62 one-second samples on all 8 GPUs. The only clock
  event reason active is SW Power Cap (0x4), at 549 W.
- The power limit is **550 W. The default and the maximum are both 600 W.**
  `nvidia-smi -pl` returns "Insufficient Permissions" in this container. I tested it by
  re-setting the current value, which changes nothing.
- SM clocks run 2227-2381 MHz against a 3090 MHz maximum. The hotter cards run slower
  (GPU5: 86 C, 2227 MHz; GPU3: 72 C, 2381 MHz).
- Host pipeline (log fields since the restart, 10873 units): prep is 0.011 s per unit, up is
  0.024 s, and stage-3 wall time is 0.126 s on a worker thread. None of these block launches.
- `gpu=` interval: 36% of the summed `gpu=` belongs to small units that report < 3e9 res/s
  (gpu = 2-4 s for 4.67e9 residues). These are not idle gaps, because SM busy is 100%. They
  are small kernels queued behind the sibling MPS engine's large units (MOD 2.5e15-4.9e15,
  1-4 s each). The `gpu=` field cannot be used for per-unit efficiency.

### Ranked suggestions

1. **(Done, measured +21.6%)** Montgomery stage 3 and `OMPT=8`. Keep `OMPT` at or below about
   12 on this host. Total CPU demand must stay below the 245.76-core quota, or CFS throttling
   freezes the launch threads.
2. **Raise the GPU power limit to 600 W.** This needs the host (Vast) side; it cannot be done
   in the container. Expected gain: +3-6%. This is an estimate, not a measurement: +9% power
   at the cap, and throughput scales sub-linearly with power. To validate, compare `rate.sh 180`
   before and after; per-GPU clocks should rise. The only path is a different offer or host
   whose limit is already 600 W. Check with
   `nvidia-smi --query-gpu=power.limit,power.default_limit --format=csv`.
3. **Measure kernel time with events, not launch intervals.** This is instrumentation only, for
   the next tuning step. Record a start event before the launch, and log
   `cudaEventElapsedTime(start, done)` as `kern=` next to `gpu=`. Once energy per residue is the
   limit, the next gains come from kernel efficiency (instructions per residue) or from the
   unit mix. Neither can be measured from `gpu=`. Expected gain: 0% by itself.
4. Not worth doing: using the ~195 idle cores for CPU sieving. The CPU engine measured about
   6e6 res/s per core (PROGRESS.md), so 195 cores would add about 1e9 res/s, roughly 0.1%
   (estimate).

No new in-container suggestion with measured evidence of >= 3% exists in this state. The
machine is at 100% SM and at the power cap, and the CPU has 4.8x headroom. The target of
7.45e11 is exceeded (7.63e11 over 180 s).

## Cycle 3 — 2026-10-03 17:59–18:13 UTC (same 16 engines, OMPT=8; work per joule / per kernel-second)

Raw data: `cycle3/` (units_ts.txt = timestamped log lines, sys_ts.txt = cgroup + GPU every 5 s,
plan_classes.txt, plan_cost_{0,1}.txt), `ab_ksolo.txt`, `ab2_ksolo.txt`, `ab3_ksolo.txt`.
**GPU use by me:** about 30 s of KBENCH kernels in total, on GPUs 0-3, including one 15 s unit on
GPU 2. They ran with the patched copy `cycle3/patchwork/apsearch_cuda_kern`, writing their
output only into `cycle3/patchwork/`.

### Drift and system state
- CPU: 52.5 cores (cgroup, 278 s), **0 of 2784 periods throttled**. Load ~55.
- The GPUs heat-soaked: GPU0 went from 78 C / 2275 MHz (17:51) to 87 C / 2171 MHz, and GPU5
  runs at 90 C. Clocks are down ~3-5% at the same 550 W. This is hardware and cannot be fixed
  in the container. It explains part of any slow downward drift in res/s.
- New bests in the window were 48, 50, 49 and 42 (no 55+). These are per engine since its
  restart and do not affect the rate.

### 1. Unit mix: GPU cost per residue depends on the unit class (measured)
The engine was run with KBENCH, which synchronises after each kernel. Units were interleaved
with the reference class MOD=1.1337e15 on a loaded production GPU. The table gives kernel
res/s relative to the reference class measured in the same run:

| class (MOD) | plan residue share | thr | res/thr | rel. res per kernel-s |
|---|---|---|---|---|
| 1.1337e15 (reference) | 15.3% | 1.73M | 2695 | 1.00 |
| 2.5171e15 | 44.8% | 84.9M | 4015 | 0.85 |
| 2.0917e15 | 3.5% | 6.5M | 4014 | 0.83 |
| 1.6686e15 | 1.5% | 2.7M | 4016 | 0.80 |
| 1.4704e15 | 1.1% | 2.0M | 4015 | 0.85 |
| 1.3142e15 | — | 1.73M | 3577 | 0.81 |
| 1.3879e15 | 1.0% | 1.73M | 4015 | 0.72 |
| 1.7893e15 | 1.9% | 3.4M | 4016 | 0.65 |
| 3.0517e15 (one unit, 15 s) | 2.1% (big-MOD classes total ~23%) | 84.9M | 5767 | 0.68 |

Yield per residue in the planner is flat across classes (0.997-1.004 relative,
`cycle3/plan_classes.txt`). E58 per GPU-second is therefore proportional to the column above.
The planner budgets and ranks by residues, so it buys the expensive classes at the same price
as the cheap one.

Planner simulation (`cycle3/plan_cost.py` is a copy of tools/plan_units.py with a cost factor
per class; same done set, `MODCAP=2e16 --kmax 600000 --smax 3000`). GPU cost is in
reference-class residue-equivalents; 1e16 is about 2.9 h at the current rate.

| GPU budget | E58 residue-ranked (current) | E58 cost-ranked | gain | residues executed |
|---|---|---|---|---|
| 3e15 | 0.358 | 0.388 | +8.4% | 2.40e15 -> 2.74e15 |
| 6e15 | 0.686 | 0.724 | +5.5% | 4.8e15 -> 5.39e15 |
| 1e16 | 1.088 | 1.133 | +4.1% | 8.0e15 -> 8.75e15 |
| 1.5e16 | 1.552 | 1.605 | +3.4% | 1.2e16 -> 1.3e16 |
| 2e16 | 1.985 | 2.044 | +3.0% | 1.6e16 -> 1.72e16 |

On both metrics: **aggregate res/s also rises about 8-14%**, because the same GPU time executes
more of the cheap residues. E58 per hour rises 4-8% over the next 1-3 hours and about 3% over 6
hours. The new plan has ~1.9x more units; the host has the headroom (52 of 246 cores).

Caveats: one or two samples per class, measured under contention. The 0.68 factor for all
MOD >= 2.6e15 classes rests on a single 3.05e15 unit, and other classes are assumed 0.78.

PPG / pairing: SM busy is 100% in every 1 s sample. Small units queued behind the sibling's
large units under MPS are not lost time. Reordering units between the two engines, or changing
PPG, therefore has no measured upside (expected ~0). It cannot be tested without a restart.

### 2. `kern=` patch (ready, not applied): `kern_field.patch`
`git apply experiments/2026-10-03-profiling/kern_field.patch` (checked with `git apply --check`
against the current tree). It records a start event on the unit's stream just before the
launch. It logs `kern=` (start -> done, s) after `wait=` and `"kern_s"` in the jsonl. Nothing
else changes, and the resume and planner parsers ignore the new field.
- Compiled and run on K=206: `kern=0.0998s`, the same survivors and best as before.
- Under MPS the value includes time the kernel shares SMs with the other engine. Compare it
  per class or per GPU, not as solo time.

### Ranked suggestions
1. **Rank and budget plan units by GPU cost, not residues.** Estimated +3-8% E58/hour, plus
   ~+9% res/s, from the simulation above.
   - Command (writes a new plan file; it does not touch tools/):
     `COSTAWARE=1 MODCAP=2e16 python3 experiments/2026-10-03-profiling/cycle3/plan_cost.py --budget-res 2.2e16 --kmax 600000 --smax 3000 --out experiments/remote_plan_cost.txt`
   - Then relaunch with that plan under a new tag. Finished units are excluded through the done
     glob.
   - Better first: apply `kern_field.patch` at the same restart, collect ~20 min of `kern_s` by
     MOD, refit `cost()` in plan_cost.py, and re-plan.
   - Validation: per-class kern_s/res from production; aggregate res/s should rise.
   - Risk: the cost factors are thin (n=1-2 per class, under contention). If the big-MOD factor
     is wrong, the gain shrinks but should not turn negative, because yield per residue is flat
     across classes. Changing the plan changes the slice assignment; resume is per (K, shift),
     so that is safe.
2. Power limit 600 W (host side only), as in cycle 2.

## Cycle 4 — 2026-10-03 18:09–18:25 UTC (engines with kern_s, OMPT=8, plan v2)

Raw data: `cycle4/` (units_ts.txt, sys_ts.txt, fit.py, winfit.py, totfit.py, phases.py,
plan_cost.py, plan_cost_{0,1}.txt, solo_units.txt).

### Drift
591 s span: 51.0 cores, **0 of 5907 periods throttled**. GPUs are at 99-100% util and
544-550 W. SM clocks are flat at 2174-2322 MHz with no further drop since cycle 3, and
temperatures are 76-90 C. Four new n=55 lines appeared (K=41949, 26019, 38763, 141423); none
is >= 56.

### 1. Refitting cost per class from production kern_s: the data cannot do it
- **Per-unit kern_s under MPS cannot be attributed to a class.**
  - Reference class (MOD 1.1337e15, n=2095): median kern/res is 0.0155 ns, but pooled
    sum(kern)/sum(res) is 0.030-0.056 ns per GPU, 2-3.5x the median. Small kernels queue behind
    the sibling engine's 85M-thread grids, and that wait is counted in their kern_s.
  - The large classes come out at 0.2-0.6 of the reference (2.5171e15: n=81, pooled 0.17-0.32
    per GPU). That is the reverse of KBENCH, and it is the same artefact: block dispatch is
    first-come, so a large grid runs nearly alone while the small kernels wait.
  - `cycle4/cost_fit.json` and the fit.py output are kept only as evidence of this.
- **A conservation fit does not identify the costs.** The model is
  "per GPU and window, sum of res x cost = window x speed", with 30 s or 60 s windows, or whole
  spans of 899 s (cycle 3) and 490 s (cycle 4). It returned negative costs and negative speeds.
  - Per-GPU mixes differ too little (large-class share 24-41%), and large units are too lumpy
    (1 unit is up to ~60-100 s of GPU) for these spans.
  - Direct look: GPU1 had the highest large-class share in both spans. In cycle 3 it was the
    fastest GPU (9.65e10 vs mean 8.9e10); in cycle 4 it was below the mean (9.36e10 vs
    9.55e10). That is no evidence either way at the ~3% level.
- **Phase split** (small kernels only in flight vs a large kernel in flight, 0.05 s bins): it
  gives 4.0e10 res/s for small-only phases against 1.2e11 for large phases. It depends on the
  same kern intervals, so it inherits the queueing artefact. It is not a cost estimate either.

**Plain verdict:** the only cost numbers are still the cycle-3 KBENCH samples: 1-2 per class,
measured as a non-MPS context time-sliced against the production MPS server. Production data
neither confirms nor refutes them. **The +3-8% gain is unverified.**

### 2. Cost-ranked plan (generated; deploy only after the solo benchmark below)
`experiments/remote_plan_cost.txt`: 1,886,972 units with the KBENCH cost model
(ref 1.0, 2.5171e15 1/0.85, MOD >= 2.6e15 1/0.68, other 1/0.78), `--budget-res 2.6e16`.
- The format matches experiments/remote_plan.txt ("K shift\n", ASCII, trailing newline):
  0 format violations and 0 duplicates.
- 0 finished units at generation time. When I checked a few minutes later, 3380 of them had
  been finished by the running engines. Filter again immediately before launch (command below),
  because the engine's resume only reads its own tag/slice file.

Expected 58s at equal GPU time. GPU time is in reference-class residue-equivalents; at
~7.6e11 res/s and mean cost ~1.25, 1 h is about 3.42e15. This holds **only under the KBENCH
cost model**:

| horizon | residue-ranked (current) | cost-ranked | gain | residues |
|---|---|---|---|---|
| 1 h | 0.395 | 0.427 | +8.1% | 2.73e15 -> 3.13e15 |
| 3 h | 1.093 | 1.140 | +4.3% | 8.24e15 -> 9.04e15 |
| 6 h | 1.993 | 2.054 | +3.1% | 1.64e16 -> 1.76e16 |

If the true costs are flat (all ~1), the cost-ranked plan is **worse**. Measured from the same
two planner runs at equal residues: 0.410 vs 0.431 at 3e15 (-4.9%), 1.249 vs 1.298 at 1e16
(-3.8%), 2.286 vs 2.350 at 2e16 (-2.7%). That is roughly -5% / -4% / -3% at 1 h / 3 h / 6 h.
The plan buys reference-class units at lower yield per residue for no saving. Upside and
downside are about the same size, so the decision rests on the cost measurement.

### Ranked suggestions
1. **Settle the costs with a solo benchmark before switching plans** (about 1-2 min of one GPU).
   - At the next restart leave one GPU out, e.g. `GPUS="0 1 2 3 4 5 6" PPG=2` (14 engines), or
     pause its two engines. Then run, on that idle GPU:
     `OMP_NUM_THREADS=8 KBENCH=1 CUDA_VISIBLE_DEVICES=7 ./build/apsearch_cuda --nterms 58 --units experiments/2026-10-03-profiling/cycle4/solo_units.txt --slice 0 1 --modcap 20000000000000000 --b2 10000 --kernel 31 --t0 24 --prep 8 --report 55 --out /tmp/solo.jsonl 2>&1 | grep ksolo`
   - solo_units.txt holds 36 units from 9 classes, each interleaved with a reference-class unit.
   - If the classes come out within ~10% of the reference, keep the current plan. If they match
     KBENCH (0.65-0.85), deploy the cost plan: write the measured costs into a JSON
     `{"ref":1,"c2517":..,"big":..,"other":..}` and regenerate with
     `COSTFILE=that.json COSTAWARE=1 MODCAP=2e16 python3 experiments/2026-10-03-profiling/cycle4/plan_cost.py --budget-res 2.6e16 --kmax 600000 --smax 3000 --out experiments/remote_plan_cost.txt`.
2. Before any launch on a regenerated plan, drop the units finished since generation:
   `python3 -c "import json,glob;d={(j['K'],j['shift']) for f in glob.glob('experiments/*/*.jsonl') for l in open(f) if '\"covered\"' in l for j in [json.loads(l)]};L=[l for l in open('experiments/remote_plan_cost.txt') if tuple(map(int,l.split())) not in d];open('experiments/remote_plan_cost.txt','w').writelines(L);print(len(L))"`
3. Power limit 600 W (host side), unchanged from cycle 2.

## Cycle 5 — 2026-10-03 18:26–18:39 UTC (v3: cost-ranked plan, 16 engines, OMPT=8)

Raw data: `cycle5/` (units_ts.txt = v3 log lines with timestamps, sys_ts.txt, score.py,
epw.py, variants.sh, variants_sum.py, variant_units.txt).

### 1. Validating the switch: the gain is there, and larger than the planner predicted
Method: each logged unit gets its planner score E(K, shift) from tools/plan_units.py
(unit_shape, size_term, WCAL calibration). Units are summed over a wall-clock window taken from
the timestamped log tailers, skipping the first 60 s.

| run | window | units | res/s | E58 per wall-hour | E58 per 1e15 res |
|---|---|---|---|---|---|
| v2 (residue-ranked), 17:59-18:14 | 839 s | 25,494 | 7.09e11 | 0.400 | 0.157 |
| v2, 18:09-18:24 | 836 s | 25,004 | 7.08e11 | 0.392 | 0.154 |
| **v3 (cost-ranked), 18:27-18:37** | 614 s | 121,640 | **9.55e11** | **0.472** | 0.137 |

- **E58 per wall-hour: +18% to +21%. Aggregate res/s: +35%.** Yield per residue fell 11%, as
  the planner intended.
- The planner simulation predicted +8% E58/h at 1 h. Production does better, because the
  solo cost table under-predicts the gap.
  - The solo table (ref 1.0, 2.5171e15 at 1/0.76, MOD >= 2.6e15 at 1/0.62, others ~1/0.8)
    gives a v2 mean cost of ~1.25, i.e. 8.3e11 res/s at the solo reference rate of
    1.3e11 per GPU. v2 measured 7.08e11, 85% of that.
  - For v3 (46% reference, 54% 2.5171e15 residues) the same table gives ~8.9e11; v3 measured
    9.55e11, 107% of that.
  - So under MPS the large-MOD units ran worse than solo, and/or the small units ran better.
    That is consistent with two 85M-thread grids from two engines competing for L1/L2.
- Host side: 198 units/s (was ~30). CPU 60.4 cores (was 51), **0 of 6768 periods throttled**.

### Drift
676 s span: all GPUs at 100% util and 550 W. Clocks 2149-2292 MHz (about 25 MHz lower than
cycle 4), temperatures 78-90 C. This is stable heat soak; nothing changed.

### 2. Why the non-reference classes cost more per residue (from kernels_cmp.cuh kernel 31)
- **Work per thread.** Each thread walks c1 x c2 residues: the two largest pinned components,
  count q - 58 each, in groups of NCH = 7 adjacent words. Phase 1 runs T0 = 24 fully unrolled
  tier-C rows per group, with a per-word load guard from row 8. Surviving words are then
  compacted per warp, and the tail runs rows 24..607 one word per lane.
- **Reference class** (pins up to 113; inner pair 107, 113; 2695 residues per thread):
  tier-C starts at **131, 137, 149, ...**
- **2.5171e15 class = K divisible by 59.** 59 divides d, so it cannot be pinned (it is useless
  there), and 131 is pinned instead (MOD ratio 2.2203 = 131/59). Inner pair 113, 131;
  4015 residues per thread. Tier-C starts at **137**: the strongest row, 131 (kills ~44% of
  the bits), is gone. Each additional pin of this kind (MOD >= 2.6e15 classes: 131 and 137
  pinned) removes the next-strongest row.
- **Consequence.** After the 24 phase-1 rows, a larger fraction of words is still alive
  (about 1.8x more bits per removed strong row; model estimate). Fewer lanes go dead and skip
  their loads, so phase 1 does more loads. More words also reach the compacted tail, where each
  costs about one L2 load per row.
- Per residue the GPU therefore does more rows of work. This is the price of residues that
  pinning has already pre-filtered. It is not a bug, and cost-aware planning is the right
  response to it.
- **Padding (exact, minor).** c2 = 73 (q = 131) with NCH = 7 means 11 groups of 7 = 77 slots:
  94.8% useful, against 98.2% for the reference (c2 = 55, 56 slots). NCH = 5 would make it
  97.3%.
  - Residue-weighted over the v3 plan: NCH 5/6/7/8 give 0.986/0.927/0.964/0.945. Even a perfect
    per-unit NCH choice gains at most **+2.3% from padding alone**, below the 3% bar. It only
    matters if NCH = 5 is not slower for other reasons.

### 3. Kernel parameter sweep (needs a GPU without production engines; ready to run)
`experiments/2026-10-03-profiling/cycle5/variants.sh 7` runs kernel 31 solo with KBENCH on 7
units:
- reference K=11597, 513, 39519
- 2.5171e15 class K=135936, 17228 (3.4e11 residues each, ~3.2 s)
- 2.0917e15 class K=93152

It sweeps (nch, t0) over 7/24 (twice, as a drift check), 7/20, 7/28, 7/32, 7/36, 8/24,
8/28, 6/24, 6/28, 5/24. That is about 11 x ~9 s, roughly 2 minutes of GPU 7. It prints each
unit's time relative to 7/24 and checks that surv/bits/conf are identical across configs.
- If the 2.5171e15 units improve by >= 5% with a different (nch, t0) while the reference
  units do not get slower, the cheap change is a per-unit (nch, t0) choice by inner pair.
  That needs a small code change in prepare()/launch (pad = nch - 1 is set per unit in
  prepare). I will draft it as a patch once the numbers exist.
- If nothing beats 7/24 by >= 3% on that class, the kernel is done for this host.

### Ranked suggestions
1. Keep v3 (measured +18-21% E58/h, +35% res/s).
2. Run `cycle5/variants.sh 7` at the next convenient pause of GPU 7 (~2 min) and send me
   `cycle5/variants_out.txt`. Expected gain: unknown, at most a few percent; padding alone
   bounds the NCH part at +2.3%.
3. Power limit 600 W (host side), unchanged.

## Cycle 6 — 2026-10-03 18:43–18:56 UTC (v3 with --nch 8 --t0 28, restarted 18:42:12)

Raw data: `cycle6/` (units_ts.txt, sys_ts.txt, acct.py, large_units.txt, ab_mps.sh).

### 1. Did (8,28) regress in production? No sign of it; I cannot confirm the +2% either
| window | plan region (res share) | res/s | E58 per wall-hour | E58 per 1e15 res | ref-equiv GPU rate (c2517 cost 1.32 / 1.40) |
|---|---|---|---|---|---|
| v3 (7,24), 18:28-18:41, GPU7 excluded | ref 83%, c2517 17% | 9.43e11 | 0.467 | 0.138 | 1.227e11 / 1.242e11 |
| v3 (8,28), 18:44-18:54 | ref 36%, c2517 64% | 8.83e11 | 0.440 | 0.138 | 1.329e11 / 1.385e11 |

- Raw res/s and E58/h fell 6%. That is because the plan moved into a 59|K-heavy region: the
  2.5171e15 share went from 17% to 64% of residues, and those cost more per residue.
- Normalised by class cost, the GPU rate went **up** 8-11%. The sweep predicted about +2%, so
  the normalisation is the uncertain part. The conclusion holds unless the production cost of
  the 59|K class is below ~1.2, and the estimates below are 1.36-1.60.
- E58 per wall-hour is still +10-12% over v2 (0.39-0.40).

### 2. Host side is not a limiter
- Per-engine accounting (acct.py): sum(kern_s)/wall = **0.99 per engine** in both windows.
  Each engine has a kernel on the GPU 99% of the time, and the sibling engine covers the
  other 1%.
  - (7,24) window: mean kern 94-96 ms, up 2.4-3.0 ms (async, overlapped), wait ~92 ms.
  - (8,28) window: 70 units/s.
  - GPU SM busy is 99-100% in every sample.
- Expected gain from `experiments/2026-10-03-structure/botl/queued_launch.patch`: **<= 1%**
  (estimate; the per-engine gap is ~1% and is already filled by the sibling).
- It no longer applies cleanly: `patch -p0 --dry-run` in src/c gives 3 of 11 hunks FAILED.
  It conflicts with the kern_field change (both add a `start` event and a gpu_s/kern path).
  Not worth updating at < 3%. I did not prepare a new version.

### 3. Large units under MPS vs solo
- I solved two production regimes for two unknowns, per GPU-second: v3 (7,24) with ref plus
  2.5171e15, against v2 windows with ref, other, 2.5171e15 and MOD >= 2.6e15. I assumed the
  solo ratio for big vs 2.5171e15 (1.22, or 1.0 as a bound) and other = 1.0-1.25.
- Results (relative to ref):

  | | production (MPS) | solo |
  |---|---|---|
  | 2.5171e15 | 1.36-1.43 (1.54-1.60 if big = 2.5171e15) | 1.32 |
  | MOD >= 2.6e15 | 1.66-1.75 | 1.61 |
  | ref-equivalent GPU rate | 1.23-1.28e11 | ~1.3e11 |

- The **MPS penalty for large units is ~3-8% per large unit, and ~0% for ref**. Large units
  are ~60% of v3 GPU time, so dedicated single-engine GPUs for them would gain an estimated
  **~2-5% overall**.
  - This rests on a 2-equation solve across windows 30 min apart (clock drift ~1%). It is not
    measured evidence of >= 3%.
  - Splitting the plan by class across GPUs also costs flexibility. Single-engine GPUs lose
    the sibling that hides the ~1% launch gap, which is negligible for 3 s kernels.
- **Decisive A/B (~13 min of one paused GPU):** `experiments/2026-10-03-profiling/cycle6/ab_mps.sh 7`
  - 2 engines under MPS, then 1 engine, then 2 engines again, 240 s each, all on 200
    2.5171e15-class units from the END of the v3 plan.
  - Outputs go to /tmp only, so nothing is recorded as coverage.
  - Deploy a class split only if 1-engine res/s >= 1.05x the 2-engine runs.

### 4. Drift
683 s span: 49.6 cores, **0 of 6828 periods throttled**. GPUs at 99-100% and 549-550 W.
Clocks 2164-2309 MHz (about +15 MHz versus cycle 5), temperatures 78-90 C. Stable.

### Ranked suggestions
1. Keep v3 with (8,28). No regression is visible.
2. Optional: run `cycle6/ab_mps.sh 7`. A class split is worth doing only at >= 5% per large
   unit (expected overall 2-5%, estimate).
3. Do not spend effort on queued launches (<= 1%).
4. Power limit 600 W (host side), unchanged.

## Watchdog — from 2026-10-03 19:05 UTC (v3, --nch 8 --t0 28, 16 engines)

Each line covers the preceding window. Sources are `watch/units_ts.txt` (timestamped engine
log lines) and `watch/samples.txt` (5 s samples); the line is produced by `watch/watch.py`.
GPU busy/MHz/C are per-GPU means over the window.

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

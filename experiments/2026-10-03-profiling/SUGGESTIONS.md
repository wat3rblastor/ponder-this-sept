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

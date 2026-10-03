# 2026-10-03-cuda — G3 campaign on an NVIDIA GB10 (resumes the Metal campaign)

Engine: `build/apsearch_cuda` (`make cuda`), the CUDA twin of `src/c/apsearch_gpu.m`.
Same work units: `d = K * D0`, `D0 = 382160924970`, stage 1 pins `a` mod `{3,2,5}` and the
seven smallest bad primes above 58 (`59,71,83,89,101,107,113`, `MOD = 1133661268029390`),
stage 2 bitmasks every other bad prime `<= 10000`, stage 3 is the exact test in
`loesch_core.h` under OpenMP. A unit `(K, shift)` covers `a` in `[64*shift*MOD, 64*(shift+1)*MOD)`.

## Machine

- NVIDIA GB10 (Blackwell, cc 12.1, 48 SMs), CUDA 13.2, driver 595.99. Unified memory, 121 GB.
- 20 CPU cores. Stage 3 takes ~0.03 s per unit on them; the GPU kernel is the whole cost.
- Throughput: 3.45e9 residues/s = one 4.67e9-residue unit per 1.35 s (M2 Metal: 11-12 s).

## Validation (before launch)

- Identical hit sets to the CPU engine (`build/apsearch`) on four small units
  (`K=205,206`, shifts 0-1, `--modcap 1.1e13 --b2 2000 --report 8`): 2814 hits each, `diff` clean.
- Rediscovers the n=47 record at `K=205` (`a=2646171143023357`) and n=43 at `K=128` with
  `--b2 2000` (its 58-window tail is killed by tier-C primes above 2000, so at `--b2 10000`
  it is correctly not reported as a 58-window survivor).
- Magic-multiply modulo checked against `%` for every tier-C prime at startup.

## Partition and coverage

Shift plane 0, `K = 790..3365`, one process (`c1.jsonl`, `c1.log`), launched 2026-10-03 03:17.
`K <= 789` at shift 0 is covered by `experiments/2026-10-03-gpu` (`g2/g4/g5.jsonl`).
`--resume` reads `c1.jsonl` and skips recorded units, so a restart never re-searches.
`--report 36`: every run of length >= 36 is logged as a `hit` line, to calibrate the
per-term decay (hits at 36..47 per unit) and project the cost of 57.

## Why this ordering (search space ranking by term size)

Per residue the hit rate is `∝ ∏_B q/(q-58) · ρ(T)^57`, `ρ ∝ 1/sqrt(ln T)`, with
`T ≈ 64·MOD·(shift+1) + 57·K·D0`. So units are ranked by `T`: shift 0 ascending `K` first
(`T` from 7.3e16 at K≈200 to 1.5e17 at K=3365), then `{shift 0, K 3366..6730}` and
`{shift 1, K <= 3365}` which cost the same and have the same `T`, and so on. The 7-prime tier B
is right for every `K` below ~4e5: adding 131 makes `64·MOD = 9.5e18` and the `ρ^57`
penalty (×30) dwarfs the ×1.8 pinning gain; dropping 113 is only better for `K < ~300`, which
is already covered.

Model (PROGRESS.md 2026-10-02, calibrated to n=31..36 counts): expected 57-term runs are
~2.2e-4 per unit at K≈800 falling to ~1.1e-4 at K≈3365, i.e. ~0.5 over this plane, and
~0.35 per further plane/band of 3365 units (73 min each at 1.35 s/unit). The pessimistic
branch anchored on the single n=47 in 748 units is ~20x lower. The `--report 36` counts
decide which branch is real.

# Running the search on a rented multi-GPU machine

Requirements on the rented box: NVIDIA GPUs, `nvcc` (CUDA devel image; CUDA >= 12.8 for
RTX 50xx / RTX PRO Blackwell), `gcc` with OpenMP, about 4 CPU cores per GPU. No Python needed there.

## Driven from the orchestrator box (normal path)

The orchestrator's SSH key must be authorized on the rented box.

```
tools/deploy_remote.sh setup HOST PORT        # copy repo, build for that GPU, benchmark one unit
tools/deploy_remote.sh start HOST PORT 40     # replan, split 40:1 remote:local, launch all GPUs
tools/deploy_remote.sh pull  HOST PORT &      # fetch remote logs every 30 s
tools/autopromote.sh 60 &                     # verifies + records any new longest AP (already running)
```

`setup` prints the measured residues/s of one GPU; set the ratio in `start` to
(remote GPUs x that rate) / (this box's ~1.5e10). `stop` halts the remote engines; every
finished unit is a line in `experiments/remote/*.jsonl`, and rerunning the launcher resumes.

## Standalone on the rented box

```
git clone git@github.com:wat3rblastor/ponder-this-sept.git && cd ponder-this-sept
tools/build_here.sh                                          # must show best=47 for K=205
MODCAP=2e16 python3 tools/plan_units.py --budget-res 2e16 --kmax 600000 --smax 3000 --out experiments/remote_plan.txt
tools/multi_gpu.sh experiments/remote_plan.txt <new-tag>     # 2 engines per GPU under MPS
grep -h '\*\*\* n=' experiments/remote/*.log | sort -t= -k2 -n | tail -3
python3 src/verify.py <a> <d> <n> --maximal                  # verify a hit
```

A watcher stops all engines when any log shows a run of 58 or more. Nothing is a record until
`src/verify.py` and `src/crosscheck.py` both pass (`tools/promote.py <a> <d> <n>`).

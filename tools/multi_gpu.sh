#!/bin/bash
# Run a plan on the GPUs of THIS machine: PPG engines per GPU, disjoint slices.
#
#   [GPUS="0 1 2"] [PPG=3] [TH=8] tools/multi_gpu.sh <plan.txt> <tag>
#
# Engine j of N = (#GPUS * PPG) takes plan entries with index % N == j and runs
# on GPU GPUS[j / PPG]. Several engines per GPU keep the device busy while each
# one's host-side per-K setup runs (a single engine leaves the GPU idle then).
# Output goes to experiments/remote/<tag>_s<j>.jsonl/.log and is resumable
# (rerun the same command after an interruption). A watcher stops every engine
# as soon as any log reports a run of 58 or more.
# Build first:  tools/build_here.sh
cd "$(dirname "$0")/.." || exit 1
PLAN=$1; TAG=$2
GPUS=${GPUS:-$(nvidia-smi --query-gpu=index --format=csv,noheader | tr '\n' ' ')}
PPG=${PPG:-2}
# stage 3 (exact test, Montgomery since 2026-10-03) needs few cores; its worker thread ignores
# --threads, and an uncapped OpenMP team (one thread per core, per engine) oversubscribed the box
OMPT=${OMPT:-8}
set -- $GPUS; NG=$#; N=$(( NG * PPG ))
NCPU=$(nproc); TH=${TH:-$(( NCPU / N ))}; [ "$TH" -lt 2 ] && TH=2
E=experiments/remote; mkdir -p $E
# MPS=1 (default): engines on one GPU share it through the CUDA MPS server instead
# of time-slicing (measured ~+15% per GPU). The daemon is started once and reused.
if [ "${MPS:-1}" = 1 ]; then
  export CUDA_MPS_PIPE_DIRECTORY=$PWD/build/mps/pipe CUDA_MPS_LOG_DIRECTORY=$PWD/build/mps/log
  mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
  pgrep -x nvidia-cuda-mps >/dev/null || nvidia-cuda-mps-control -d
fi
j=0
for g in $GPUS; do
  for p in $(seq 1 $PPG); do
    env CUDA_VISIBLE_DEVICES=$g ${OMPT:+OMP_NUM_THREADS=$OMPT} OMP_WAIT_POLICY=passive setsid nohup ./build/apsearch_cuda --nterms 58 --units "$PLAN" \
        --slice $j $N --modcap ${MODCAP:-20000000000000000} --b2 10000 --report 44 \
        --threads $TH ${EXTRA:---kernel 31 --nch 8 --t0 28 --prep 8 --report 55} --out "$E/${TAG}_s${j}.jsonl" --resume \
        >> "$E/${TAG}_s${j}.log" 2>&1 < /dev/null &
    j=$(( j + 1 ))
  done
done
setsid nohup bash -c "while pgrep -x apsearch_cuda >/dev/null; do
    if grep -qE '\*\*\* n=(5[89]|[6-9][0-9])' $E/*.log 2>/dev/null; then
      for p in \$(pgrep -x apsearch_cuda); do kill -INT \$p; done; echo \"\$(date) FOUND >=58, engines stopped\" >> $E/watch.log; break; fi
    sleep 10; done" > /dev/null 2>&1 < /dev/null &
echo "launched $N engines on GPUs [$GPUS] ($PPG per GPU, $TH CPU threads each), logs in $E/${TAG}_s*.log"

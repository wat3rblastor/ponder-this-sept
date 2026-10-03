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
PPG=${PPG:-3}
set -- $GPUS; NG=$#; N=$(( NG * PPG ))
NCPU=$(nproc); TH=${TH:-$(( NCPU / N ))}; [ "$TH" -lt 2 ] && TH=2
E=experiments/remote; mkdir -p $E
j=0
for g in $GPUS; do
  for p in $(seq 1 $PPG); do
    CUDA_VISIBLE_DEVICES=$g OMP_NUM_THREADS=$TH OMP_WAIT_POLICY=passive setsid nohup ./build/apsearch_cuda --nterms 58 --units "$PLAN" \
        --slice $j $N --modcap ${MODCAP:-20000000000000000} --b2 10000 --report 44 \
        --threads $TH ${EXTRA:---kernel 31 --t0 24} --out "$E/${TAG}_s${j}.jsonl" --resume \
        >> "$E/${TAG}_s${j}.log" 2>&1 < /dev/null &
    j=$(( j + 1 ))
  done
done
setsid nohup bash -c "while pgrep -x apsearch_cuda >/dev/null; do
    if grep -qE '\*\*\* n=(5[89]|[6-9][0-9])' $E/*.log 2>/dev/null; then
      for p in \$(pgrep -x apsearch_cuda); do kill -INT \$p; done; echo \"\$(date) FOUND >=58, engines stopped\" >> $E/watch.log; break; fi
    sleep 10; done" > /dev/null 2>&1 < /dev/null &
echo "launched $N engines on GPUs [$GPUS] ($PPG per GPU, $TH CPU threads each), logs in $E/${TAG}_s*.log"

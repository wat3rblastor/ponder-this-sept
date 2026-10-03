#!/bin/bash
# Run a plan on every GPU of THIS machine: one engine per GPU, disjoint slices.
#
#   tools/multi_gpu.sh <plan.txt> <tag>
#
# GPU g takes plan entries with index % NGPU == g. Output goes to
# experiments/remote/<tag>_g<g>.jsonl/.log and is resumable (rerun the same
# command after an interruption). A watcher stops every engine as soon as any
# log reports a run of 58 or more.
# Build first:  tools/build_here.sh
cd "$(dirname "$0")/.." || exit 1
PLAN=$1; TAG=$2
NG=$(nvidia-smi --query-gpu=index --format=csv,noheader | wc -l)
NCPU=$(nproc); TH=$(( NCPU / NG )); [ "$TH" -lt 2 ] && TH=2
E=experiments/remote; mkdir -p $E
for g in $(seq 0 $((NG - 1))); do
  CUDA_VISIBLE_DEVICES=$g setsid nohup ./build/apsearch_cuda --nterms 58 --units "$PLAN" \
      --slice $g $NG --modcap ${MODCAP:-20000000000000000} --b2 10000 --report 44 \
      --threads $TH --out "$E/${TAG}_g${g}.jsonl" --resume \
      >> "$E/${TAG}_g${g}.log" 2>&1 < /dev/null &
done
setsid nohup bash -c "while pgrep -x apsearch_cuda >/dev/null; do
    if grep -qE '\*\*\* n=(5[89]|[6-9][0-9])' $E/*.log 2>/dev/null; then
      pkill -INT -x apsearch_cuda; echo \"\$(date) FOUND >=58, engines stopped\" >> $E/watch.log; break; fi
    sleep 10; done" > /dev/null 2>&1 < /dev/null &
echo "launched $NG engines ($TH CPU threads each), logs in $E/${TAG}_g*.log"

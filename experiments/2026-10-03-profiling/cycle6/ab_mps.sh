#!/bin/bash
# A/B on ONE paused GPU: 2 engines under MPS vs 1 engine, on 2.5171e15-class units only
# (taken from the END of the v3 plan, outputs to /tmp, not recorded as coverage).
# usage: experiments/2026-10-03-profiling/cycle6/ab_mps.sh [GPU] [SECONDS]
cd /workspace/ponder-this-sept || exit 1
G=${1:-7}; S=${2:-240}; D=experiments/2026-10-03-profiling/cycle6; U=$D/large_units.txt
export CUDA_MPS_PIPE_DIRECTORY=$PWD/build/mps/pipe CUDA_MPS_LOG_DIRECTORY=$PWD/build/mps/log
run() { # $1 = number of engines
  rm -f /tmp/abm_*.jsonl
  for j in $(seq 0 $(( $1 - 1 ))); do
    env CUDA_VISIBLE_DEVICES=$G OMP_NUM_THREADS=8 OMP_WAIT_POLICY=passive ./build/apsearch_cuda --nterms 58 --units $U \
      --slice $j $1 --modcap 20000000000000000 --b2 10000 --kernel 31 --nch 8 --t0 28 --prep 8 --report 55 \
      --out /tmp/abm_$j.jsonl > /tmp/abm_$j.log 2>&1 &
  done
  sleep 20; a=$(cat /tmp/abm_*.jsonl 2>/dev/null | python3 -c "import sys,json;print(sum(json.loads(l)['res'] for l in sys.stdin if '\"res\"' in l))")
  sleep $S;  b=$(cat /tmp/abm_*.jsonl 2>/dev/null | python3 -c "import sys,json;print(sum(json.loads(l)['res'] for l in sys.stdin if '\"res\"' in l))")
  for p in $(jobs -p); do kill -INT $p; done; wait
  python3 -c "print('engines=$1  res/s = %.4e' % (($b - $a) / $S))"
}
run 2; run 1; run 2

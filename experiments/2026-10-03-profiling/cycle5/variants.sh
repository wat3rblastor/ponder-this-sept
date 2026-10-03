#!/bin/bash
# Solo kernel-variant sweep. Run ONLY on a GPU with no production engines (GPU=7 default).
# usage: experiments/2026-10-03-profiling/cycle5/variants.sh [GPU]
cd /workspace/ponder-this-sept || exit 1
G=${1:-7}; D=experiments/2026-10-03-profiling/cycle5; OUT=$D/variants_out.txt; : > $OUT
for cfg in "7 24" "7 20" "7 28" "7 32" "7 36" "8 24" "8 28" "6 24" "6 28" "5 24" "7 24"; do
  set -- $cfg
  echo "=== nch=$1 t0=$2" >> $OUT
  OMP_NUM_THREADS=8 KBENCH=1 CUDA_VISIBLE_DEVICES=$G ./build/apsearch_cuda --nterms 58 \
    --units $D/variant_units.txt --slice 0 1 --modcap 20000000000000000 --b2 10000 \
    --kernel 31 --nch $1 --t0 $2 --prep 8 --report 55 --out /tmp/variant_$1_$2.jsonl 2>&1 \
    | grep -E "^ksolo|^K=" | sed -E 's/ prep=.*//' >> $OUT
  rm -f /tmp/variant_$1_$2.jsonl
done
python3 $D/variants_sum.py $OUT

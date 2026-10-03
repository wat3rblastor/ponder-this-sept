#!/bin/bash
# Run the plan on every GPU of this machine, one engine per GPU, disjoint slices.
#
#   tools/multi_gpu.sh <plan.txt> <tag> [offset] [total_slices]
#
# GPU g of this box takes plan entries with index % total_slices == offset + g.
# With several machines, give each a distinct offset range and the same
# total_slices (e.g. machine A: offset 0, machine B: offset 8, total 16).
# Output: experiments/2026-10-03-cuda/<tag>_g<slice>.jsonl/.log (resumable).
# Build first:  make cuda CUDA_ARCH=sm_120   (RTX 50xx; sm_89 for RTX 40xx)
cd "$(dirname "$0")/.." || exit 1
PLAN=$1; TAG=$2; OFF=${3:-0}
NG=$(nvidia-smi --query-gpu=index --format=csv,noheader | wc -l)
TOT=${4:-$NG}
NCPU=$(nproc); TH=$(( NCPU / NG )); [ "$TH" -lt 2 ] && TH=2
E=experiments/2026-10-03-cuda; mkdir -p $E
for g in $(seq 0 $((NG - 1))); do
  sl=$((OFF + g))
  CUDA_VISIBLE_DEVICES=$g setsid nohup ./build/apsearch_cuda --nterms 58 --units "$PLAN" \
      --slice $sl $TOT --modcap ${MODCAP:-20000000000000000} --b2 10000 --report 44 \
      --threads $TH --out "$E/${TAG}_g${sl}.jsonl" --resume \
      >> "$E/${TAG}_g${sl}.log" 2>&1 < /dev/null &
done
echo "launched $NG engines, slices $OFF..$((OFF + NG - 1)) of $TOT, $TH threads each"

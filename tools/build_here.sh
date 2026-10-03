#!/bin/bash
# Build the CUDA engine for whatever GPU this machine has, then run the
# selftest unit (K=205 must report best=47 at shift 0 with --kernel 4).
cd "$(dirname "$0")/.." || exit 1
cap=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader | sort | uniq -c | tr '\n' ' ') -> sm_$cap"
rm -f build/apsearch_cuda build/loesch_core_api.o
make cuda CUDA_ARCH=sm_$cap 2>&1 | grep -iE "error|nvcc" 
[ -x build/apsearch_cuda ] || { echo "BUILD FAILED"; exit 1; }
CUDA_VISIBLE_DEVICES=0 ./build/apsearch_cuda --kmin 205 --kmax 205 --kernel 4 --report 40 2>&1 | grep -E "^K=|\*\*\* n=47"
echo "--- production kernel, one standard unit (4.67e9 residues):"
CUDA_VISIBLE_DEVICES=0 ./build/apsearch_cuda --kmin 206 --kmax 206 --kernel 31 --t0 24 --prep 8 2>&1 | grep -E "^K=" | grep -oE "surv=[0-9]+|gpu=[0-9.]+s|\([0-9.e+]+ res/s\)" | tr '\n' ' '; echo

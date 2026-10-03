#!/bin/bash
# usage: run_small.sh K mode   -- engine (10 shifts, windows with 128 <= a0/MOD < 640) vs oracle
K=$1; M=$2; T=K${K}_m${M}
rm -f en_$T.jsonl en_$T.win
APS_DUMPWIN=en_$T.win ../../build/apsearch_cuda_sq --D0 129030 --kmin $K --kmax $K --mode $M --shifts 10 --mrange 128 640 \
   --nterms 24 --b2 2000 --modcap 10000000000 --report 14 --threads 16 --kernel 31 --t0 24 --prep 4 --out en_$T.jsonl > en_$T.log 2>&1
../../build/oracle_sq --D0 129030 --K $K --mode $M --nterms 24 --b2 2000 --modcap 1e10 --report 14 --mlo 128 --mhi 640 --win or_$T.win > or_$T.hits 2> or_$T.log
python3 compare.py $K $M

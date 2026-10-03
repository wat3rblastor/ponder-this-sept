#!/bin/bash
# sweep.sh <tag> <N> <dbase> <mmax> <nproc> <histmin>
# Exhaustive sweep of d = dbase*m, m in [1,mmax], all starts, terms <= N.
# Split into <nproc> disjoint m-ranges. Raw per-m JSONL in raw/<tag>.p*.jsonl
set -e
cd "$(dirname "$0")"
TAG=$1; N=$2; DB=$3; MMAX=$4; NP=$5; HM=${6:-20}
CH=$(( (MMAX + NP - 1) / NP ))
for ((i=0;i<NP;i++)); do
  LO=$(( i*CH + 1 )); HI=$(( (i+1)*CH )); (( HI > MMAX )) && HI=$MMAX
  (( LO > MMAX )) && continue
  ./mcurve --N "$N" --dbase "$DB" --mmin "$LO" --mmax "$HI" --hist-min "$HM" \
      --out "raw/${TAG}.p${i}.jsonl" > "raw/${TAG}.p${i}.log" 2>&1 &
done
wait
echo "DONE $TAG"

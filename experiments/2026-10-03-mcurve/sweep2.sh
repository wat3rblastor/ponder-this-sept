#!/bin/bash
# sweep2.sh <tag> <N> <dbase> <mmax> <nchunk> <histmin> <chunk-ids...>
set -e
cd "$(dirname "$0")"
TAG=$1; N=$2; DB=$3; MMAX=$4; NC=$5; HM=$6; shift 6
CH=$(( (MMAX + NC - 1) / NC ))
for i in "$@"; do
  LO=$(( i*CH + 1 )); HI=$(( (i+1)*CH )); (( HI > MMAX )) && HI=$MMAX
  (( LO > MMAX )) && continue
  ./mcurve --N "$N" --dbase "$DB" --mmin "$LO" --mmax "$HI" --hist-min "$HM" \
      --out "raw/${TAG}.p${i}.jsonl" > "raw/${TAG}.p${i}.log" 2>&1 &
done
wait
echo "DONE $TAG $*"

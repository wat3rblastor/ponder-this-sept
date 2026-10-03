#!/bin/bash
# 8 disjoint fresh K blocks (K >= 30001; the band 20001..27010 is already covered).
# Identical settings to experiments/2026-10-03-laptopcost so the two runs pool.
cd /Users/bowencheng/Projects/ponder-this-sept
D=experiments/2026-10-04-decay
i=0
for K in 30001 31001 32001 33001 34001 35001 36001 37001; do
  build/apsearch --nterms 58 --kmin $K --kmax $((K+9)) --shifts 3 --b2 2000 \
      --report 30 --out $D/w$i.jsonl > $D/w$i.out 2> $D/w$i.log &
  echo $! > $D/w$i.pid
  i=$((i+1))
done
wait

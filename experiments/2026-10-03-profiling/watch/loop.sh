#!/bin/bash
cd /workspace/ponder-this-sept/experiments/2026-10-03-profiling/watch
for i in 1 2 3; do
  for k in $(seq 1 60); do   # 10 min, check alerts every 10 s over last 90 s
    sleep 10
    a=$(python3 watch.py 90 | grep ALERT | grep -v "ENGINES min=-1")
    if [ -n "$a" ]; then echo "$(date -u +%H:%M:%S) $a" >> alerts.txt; fi
    # engine count instant check
    n=$(pgrep -x apsearch_cuda | wc -l); [ "$n" -lt 16 ] && echo "$(date -u +%H:%M:%S) ENGINES now $n" >> alerts.txt
    [ -s alerts.txt ] && exit 1
  done
  python3 watch.py 600 | grep -v ALERT | tee -a ../SUGGESTIONS.md >> lines.txt
done

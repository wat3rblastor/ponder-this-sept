#!/bin/bash
# every 5 s: epoch, engine count, cgroup counters, per-GPU util/clk/temp/power
end=$(( $(date +%s) + $1 ))
while [ $(date +%s) -lt $end ]; do
  echo "T $(date +%s) $(pgrep -x apsearch_cuda | wc -l) $(grep -E '^(usage_usec|nr_periods|nr_throttled) ' /sys/fs/cgroup/cpu.stat | awk '{print $2}' | tr '\n' ' ')"
  nvidia-smi --query-gpu=index,utilization.gpu,clocks.sm,temperature.gpu,power.draw --format=csv,noheader,nounits | sed 's/^/G /'
  sleep 5
done

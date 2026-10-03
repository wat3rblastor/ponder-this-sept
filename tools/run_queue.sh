#!/bin/bash
# Run the CUDA campaign through successive (K-band, shift) segments in order of
# term size T ~ 64*MOD*(shift+1) + 57*K*D0, never idling the GPU.
# Each segment has its own jsonl and is resumable; stops once records.json
# holds n >= TARGET (default 57). Usage: tools/run_queue.sh [target]
cd "$(dirname "$0")/.." || exit 1
TARGET=${1:-57}
E=experiments/2026-10-03-cuda
BAND=3365
recn() { python3 -c 'import json;print(json.load(open("records.json"))["g3_longest"]["n"])' 2>/dev/null || echo 0; }
# wait for any engine already running (the c1 segment) to finish
while pgrep -x apsearch_cuda >/dev/null; do sleep 5; done
# level L = shift + band index; all (shift, band) with the same L have the same T
for L in $(seq 1 60); do
  for s in $(seq 0 $L); do
    b=$((L - s))
    kmin=$((b * BAND + 1)); kmax=$(((b + 1) * BAND))
    tag="s${s}_k${kmin}-${kmax}"
    [ -f "$E/$tag.done" ] && continue
    [ "$(recn)" -ge "$TARGET" ] && { echo "$(date '+%F %T') target reached" >> "$E/queue.log"; exit 0; }
    echo "$(date '+%F %T') start $tag" >> "$E/queue.log"
    ./build/apsearch_cuda --nterms 58 --kmin $kmin --kmax $kmax --shift0 $s --shifts 1 \
        --modcap 2000000000000000 --b2 10000 --report 44 --out "$E/$tag.jsonl" --resume \
        >> "$E/$tag.log" 2>&1 &
    pid=$!
    while kill -0 $pid 2>/dev/null; do
      sleep 20
      [ "$(recn)" -ge "$TARGET" ] && kill -INT $pid
    done
    wait $pid
    if tail -1 "$E/$tag.log" | grep -q "INTERRUPTED"; then
      echo "$(date '+%F %T') interrupted $tag" >> "$E/queue.log"; exit 0
    fi
    echo "{\"segment\":\"$tag\",\"shift\":$s,\"kmin\":$kmin,\"kmax\":$kmax,\"exhausted\":true}" > "$E/$tag.done"
    echo "$(date '+%F %T') done $tag" >> "$E/queue.log"
  done
done

#!/bin/bash
# Run the value-ordered unit plan (tools/plan_units.py) until records.json
# holds n >= TARGET. Resumable: the engine skips units already in the jsonl.
# Usage: tools/run_plan.sh <plan.txt> <out.jsonl> [target]
cd "$(dirname "$0")/.." || exit 1
PLAN=$1; OUT=$2; TARGET=${3:-58}
recn() { python3 -c 'import json;print(json.load(open("records.json"))["g3_longest"]["n"])' 2>/dev/null || echo 0; }
./build/apsearch_cuda --nterms 58 --units "$PLAN" --modcap ${MODCAP:-20000000000000000} --b2 10000 \
    --report 44 --out "$OUT" --resume >> "${OUT%.jsonl}.log" 2>&1 &
pid=$!
while kill -0 $pid 2>/dev/null; do
  sleep 20
  [ "$(recn)" -ge "$TARGET" ] && kill -INT $pid
done

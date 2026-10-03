#!/bin/bash
# Write STATUS.md (a phone-readable progress page), commit and push it every N seconds.
# Usage: setsid nohup tools/status_loop.sh [seconds] > /dev/null 2>&1 < /dev/null &
cd "$(dirname "$0")/.." || exit 1
N=${1:-600}
tot() { cat experiments/remote/*.jsonl 2>/dev/null | python3 -c "
import sys,json
s=0
for l in sys.stdin:
    if '\"res\"' in l:
        try: s+=json.loads(l)['res']
        except Exception: pass
print(s)"; }
prev=$(tot); tp=$(date +%s)
while true; do
  sleep "${FIRST:-$N}"; FIRST=
  cur=$(tot); tn=$(date +%s)
  python3 - "$prev" "$cur" $(( tn - tp )) "$(pgrep -x apsearch_cuda | wc -l)" > STATUS.md <<'P'
import sys, json, datetime
prev, cur, dt, eng = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
r = json.load(open("records.json")); g = r["g3_longest"]
runs = r.get("g3_runs_55plus", [])
print("# Search status\n")
print(f"Updated {datetime.datetime.utcnow():%Y-%m-%d %H:%M} UTC (refreshes about every {dt // 60} min)\n")
print(f"- **Target n >= 58: {'MET' if g['n'] >= 58 else 'not yet'}**")
print(f"- Longest verified progression: **n = {g['n']}**, a = `{g['a']}`, d = `{g['d']}`")
print(f"- Engines running: {eng} of 16")
print(f"- Throughput: {(cur - prev) / dt:.2e} residues/s over the last {dt // 60} min")
print(f"- Residues covered on the GPUs, all runs: {cur:.2e}")
print(f"- Distinct progressions of 55+ terms: {len(runs)} (" + ", ".join(str(e['n']) for e in runs) + ")")
print("\nDetails: ANSWER.md (results), PROGRESS.md (log).")
P
  prev=$cur; tp=$tn
  git add STATUS.md 2>/dev/null && git commit -q -m "status: refresh STATUS.md" -- STATUS.md 2>/dev/null && git push -q origin HEAD 2>/dev/null
done

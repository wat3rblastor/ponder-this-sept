#!/bin/bash
# Watch the running campaigns and promote any new longest progression.
#
# tools/promote.py is the gate: it refuses anything that does not pass both
# src/verify.py and src/crosscheck.py, so this loop cannot record a bad record.
# Usage: tools/autopromote.sh [poll-seconds]
cd "$(dirname "$0")/.." || exit 1
POLL=${1:-60}
LOG=experiments/autopromote.log
echo "$(date '+%F %T') autopromote started (poll ${POLL}s)" >> "$LOG"

while true; do
  # current record
  cur=$(python3 -c 'import json;print(json.load(open("records.json"))["g3_longest"]["n"])' 2>/dev/null)
  [ -z "$cur" ] && cur=0

  # best candidate across every campaign log, as "n a d"
  best=$(grep -hoE '\*\*\* n=[0-9]+ a=[0-9]+ d=[0-9]+' \
           experiments/*/*.log 2>/dev/null \
         | sed -E 's/\*\*\* n=([0-9]+) a=([0-9]+) d=([0-9]+)/\1 \2 \3/' \
         | sort -k1,1nr | head -1)
  set -- $best
  n=$1; a=$2; d=$3

  if [ -n "$n" ] && [ "$n" -gt "$cur" ]; then
    echo "$(date '+%F %T') candidate n=$n a=$a d=$d (record was $cur)" >> "$LOG"
    if python3 tools/promote.py "$a" "$d" "$n" >> "$LOG" 2>&1; then
      echo "$(date '+%F %T') PROMOTED n=$n" >> "$LOG"
      git add -A >> "$LOG" 2>&1
      git commit -q -m "record: G3 n=$n at a=$a, d=$d

Found by the running campaign and promoted automatically by
tools/autopromote.sh, which records nothing that fails either src/verify.py or
the constructive cross-check in src/crosscheck.py." >> "$LOG" 2>&1
      echo "$(date '+%F %T') committed n=$n" >> "$LOG"
      git push -q origin HEAD >> "$LOG" 2>&1 && echo "$(date '+%F %T') pushed n=$n" >> "$LOG"
    else
      echo "$(date '+%F %T') promote REFUSED n=$n -- left alone" >> "$LOG"
    fi
  fi
  # list every distinct (primitive) progression of 55+ terms in ANSWER.md
  if python3 tools/list55.py >> "$LOG" 2>&1; then
    git add records.json ANSWER.md >> "$LOG" 2>&1
    git commit -q -m "answer: list every distinct progression of 55+ terms found so far" >> "$LOG" 2>&1
    git push -q origin HEAD >> "$LOG" 2>&1
    echo "$(date '+%F %T') 55+ list updated and pushed" >> "$LOG"
  fi
  # every ~10 polls, commit and push the coverage logs so another machine can resume
  tick=$(( ${tick:-0} + 1 ))
  if [ $(( tick % 10 )) -eq 0 ]; then
    git add experiments/remote/*.jsonl >> "$LOG" 2>&1
    if ! git diff --cached --quiet; then
      git commit -q -m "coverage: units searched so far (experiments/remote)" >> "$LOG" 2>&1
      git push -q origin HEAD >> "$LOG" 2>&1
    fi
  fi
  sleep "$POLL"
done

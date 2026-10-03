#!/bin/bash
# usage: rate.sh SECONDS -> aggregate & per-engine res/s from v2 jsonl
cd /workspace/ponder-this-sept/experiments/remote
snap(){ for i in $(seq 0 15); do python3 -c "
import json,sys
s=0;n=0
for l in open('v3_s$i.jsonl'):
  try: s+=json.loads(l)['res'];n+=1
  except: pass
print($i,s,n)"; done; }
snap > /tmp/claude-0/r0.$$; t0=$(date +%s.%N); sleep $1; snap > /tmp/claude-0/r1.$$; t1=$(date +%s.%N)
python3 - $t0 $t1 /tmp/claude-0/r0.$$ /tmp/claude-0/r1.$$ <<'P'
import sys
t0,t1=float(sys.argv[1]),float(sys.argv[2]);dt=t1-t0
a={int(l.split()[0]):(float(l.split()[1]),int(l.split()[2])) for l in open(sys.argv[3])}
b={int(l.split()[0]):(float(l.split()[1]),int(l.split()[2])) for l in open(sys.argv[4])}
tot=0
for i in range(16):
  r=(b[i][0]-a[i][0])/dt; tot+=r
  print(f"s{i:2d} gpu{i//2} {r:.3e} res/s units={b[i][1]-a[i][1]}")
for g in range(8): print(f"GPU{g} {sum((b[i][0]-a[i][0])/dt for i in (2*g,2*g+1)):.3e}")
print(f"AGG {tot:.4e} res/s over {dt:.1f}s units={sum(b[i][1]-a[i][1] for i in range(16))}")
P

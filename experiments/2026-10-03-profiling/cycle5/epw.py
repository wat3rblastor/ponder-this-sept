# expected-58 per wall hour from a tailer file: epw.py units_ts.txt [skip_s] [tail_trim_s]
import sys,re
sys.path.insert(0,'/workspace/ponder-this-sept/experiments/2026-10-03-profiling/cycle5'); from score import E
fn=sys.argv[1]; skip=float(sys.argv[2]) if len(sys.argv)>2 else 60; trim=float(sys.argv[3]) if len(sys.argv)>3 else 0
rows=[]
for l in open(fn):
    p=l.split(' ',2)
    if not p[2].startswith('K='): continue
    m=re.match(r'K=(\d+) sh=(\d+) .*?res=([\d.e+]+)',p[2]); rows.append((float(p[0]),int(m[1]),int(m[2]),float(m[3])))
t0=min(r[0] for r in rows)+skip; t1=max(r[0] for r in rows)-trim
sel=[r for r in rows if t0<=r[0]<=t1]
e=res=0
for t,K,s,r in sel:
    eu,ru=E(K,s); e+=eu; res+=r
dt=t1-t0
print(f'{fn}: window {dt:.0f}s units {len(sel)} res/s {res/dt:.4e}  E58/h {3600*e/dt:.4f}  E58 per 1e15 res {1e15*e/res:.4f}')

# watch.py WINDOW_S : summarise the last WINDOW_S seconds; print one line and ALERT lines
import sys,re,time,collections
sys.path.insert(0,'/workspace/ponder-this-sept/experiments/2026-10-03-profiling/cycle5'); from score import E
W=float(sys.argv[1]); now=time.time(); ta=now-W
REF=1133661268029390
res=e=0; n=0; mix=collections.Counter()
for l in open('units_ts.txt'):
    p=l.split(' ',2)
    if not p[2].startswith('K=') or float(p[0])<ta: continue
    m=re.match(r'K=(\d+) sh=(\d+) MOD=(\d+) .*?res=([\d.e+]+)',p[2])
    K,s,MOD,r=int(m[1]),int(m[2]),int(m[3]),float(m[4])
    eu,_=E(K,s); e+=eu; res+=r; n+=1
    mix['ref' if abs(MOD-REF)<1e6 else ('c59' if abs(MOD-2517112306980510)<1e6 else ('big' if MOD>=2.6e15 else 'oth'))]+=r
T=[];G=collections.defaultdict(list)
cur=None
for l in open('samples.txt'):
    f=l.split()
    if f[0]=='T':
        cur=int(f[1]); 
        if cur>=ta and len(f)>=6: T.append((cur,int(f[2]),int(f[3]),int(f[4]),int(f[5])))
    elif f[0]=='G' and cur and cur>=ta:
        g=[x.strip() for x in l[2:].split(',')]; G[int(g[0])].append((cur,float(g[1]),float(g[2]),float(g[3]),float(g[4])))
alerts=[]
eng=min(t[1] for t in T) if T else -1
if eng<16: alerts.append(f'ENGINES min={eng}')
cores=(T[-1][2]-T[0][2])/1e6/(T[-1][0]-T[0][0]) if len(T)>1 else 0
thr=T[-1][4]-T[0][4] if len(T)>1 else 0
if thr>0: alerts.append(f'THROTTLED {thr} periods')
gs=[]
for g in range(8):
    v=G[g]; run=0; worst=0
    for x in v:
        run = run+1 if x[1]<90 else 0; worst=max(worst,run)
    if worst*5>60: alerts.append(f'GPU{g} <90% busy for {worst*5}s')
    gs.append('%d:%.0f%%/%.0f/%.0fC'%(g,sum(x[1] for x in v)/len(v),sum(x[2] for x in v)/len(v),sum(x[4-1] for x in v)/len(v)))
R=sum(mix.values()) or 1
req=(mix['ref']+1.40*mix['c59']+1.70*mix['big']+1.25*mix['oth'])/W/8
if res>0 and W>=300 and req<0.9*1.33e11: alerts.append(f'REF-EQUIV RATE {req:.3e} < 90% of 1.33e11 per GPU')
line=(f"- {time.strftime('%H:%M',time.gmtime(now))} UTC [{W/60:.0f} min]: {res/W:.3e} res/s, E58/h {3600*e/W:.3f}, "
      f"E58/1e15res {1e15*e/max(res,1):.3f}, mix ref/59|K/big/other {100*mix['ref']/R:.0f}/{100*mix['c59']/R:.0f}/{100*mix['big']/R:.0f}/{100*mix['oth']/R:.0f}%, "
      f"ref-equiv/GPU {req:.3e}, units/s {n/W:.0f}, engines {eng}, CPU {cores:.1f} cores, throttled {thr}, GPU busy/MHz/C " + ' '.join(gs))
print(line)
for a in alerts: print('ALERT',a)

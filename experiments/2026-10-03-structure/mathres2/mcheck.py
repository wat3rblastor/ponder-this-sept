import sys,math,os,collections,statistics
sys.path.insert(0,'/workspace/mathres')
os.environ['MODCAP']='2e16'
from load import *
import importlib.util
spec=importlib.util.spec_from_file_location('np_','/workspace/ponder-this-sept/tools/plan_units.py'); N=importlib.util.module_from_spec(spec); spec.loader.exec_module(N)
N.MODCAP=int(2e16)
BADALL=[p for p in N.PR if p%3==2]
E_='/workspace/ponder-this-sept/experiments/'
U=load([E_+'remote/r*.jsonl'])
ms=[];ws=[];rel=[]
for u in U:
    if not u['hitrec']: continue
    K,s=u['K'],u['shift']; sh=N.unit_shape(K); MOD,res,ly,wbar,lg=sh
    if abs(res/u['res']-1)>0.01: continue
    d=K*D0
    for h in u['hitrec']:
        a=h['a']; n=h['n']; found=[]
        for j in range(-(58-44),n-44+1+14):
            a0=a+j*d   # window start
            if a0<=0: continue
            # window must contain >=44 of run [a, a+(n-1)d]: overlap
            lo=max(0,-j); hi=min(57,n-1-j)  # indices in window that are in run
            if hi-lo+1<44: continue
            ok=True
            for i in range(58):
                t=a0+i*d
                for q in BADALL:
                    if q>58 and t%q==0: ok=False;break
                if not ok: break
            if ok:
                m=a0//MOD-64*s
                if 0<=m<64+120: found.append(m)
        if found:
            ms.append(statistics.mean(found)); ws.append(wbar+32); rel.append(u['file'].split('/')[-1][:2])
print(len(ms),'mean m obs %.1f  model wbar+32 %.1f'%(statistics.mean(ms),statistics.mean(ws)))
for b in ['r1','r2','r3','r4','r5']:
    x=[m for m,r in zip(ms,rel) if r==b]; y=[w for w,r in zip(ws,rel) if r==b]
    if x: print(b,len(x),'%.1f %.1f'%(statistics.mean(x),statistics.mean(y)), 'min %.0f max %.0f'%(min(x),max(x)))

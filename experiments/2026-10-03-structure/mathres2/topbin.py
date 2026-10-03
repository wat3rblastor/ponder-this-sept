import sys,math,collections,os
sys.path.insert(0,'/workspace/mathres')
os.environ['MODCAP']='2e16'
from load import *
from model import *
E_='/workspace/ponder-this-sept/experiments/'
def prun(p,n,N=58):
    return p**n*(1+(N-n)*(1-p))
def fac(n):
    f=[];p=2
    while p*p<=n:
        while n%p==0: f.append(p); n//=p
        p+=1
    if n>1: f.append(n)
    return f
U=load([E_+'remote/r*.jsonl']); seen=set(); rows=[]
agg=collections.defaultdict(lambda:[0,0])
for u in U:
    K,s=u['K'],u['shift']
    if (K,s) in seen: continue
    seen.add((K,s)); P.MODCAP=int(2e16); sh=shape(K); MOD,r,ly,pinned=sh
    wbar=((pinned[-1]-58)+(pinned[-2]-58))/4; d=K*D0; e=0
    for i in range(4):
        a=(64*s+(i+.5)*16+wbar)*MOD
        p=math.exp(sum(math.log(rho(a+k*d)) for k in range(0,58,6))/10)
        e+=u['conf']/4*prun(p,44)*goodpen(K)
    o=sum(1 for h in u['hitrec'] if all(math.gcd(h['a'],h['d'])%p for p in GOOD))
    T=(64*s+62)*MOD+57*d
    bq=tuple(q for q in sorted(set(fac(K))) if q%3==2 and q>58)
    key=('nbad=%d'%len(bq), 'pins=%d..%d'%(pinned[0],pinned[-1]))
    agg[key][0]+=e; agg[key][1]+=o
    agg[('nbad=%d'%len(bq),)][0]+=e; agg[('nbad=%d'%len(bq),)][1]+=o
    for q in bq: agg[('q|K',q if q<200 else 'big')][0]+=e; agg[('q|K',q if q<200 else 'big')][1]+=o
for k in sorted(agg,key=str):
    e,o=agg[k]
    if e>4: print(k,'exp %.1f obs %d ratio %.2f z %+.1f'%(e,o,o/e,(o-e)/e**.5))

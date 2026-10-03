from load import *
from model import *
E_='/workspace/ponder-this-sept/experiments/'
U=load([E_+'remote/r1_g*.jsonl'])
def prun(p,n,N=58):
    # P(max run >= n in N iid bernoulli(p))
    # dp over current run length
    st=[0.0]*n; st[0]=1.0; hit=0.0
    for _ in range(N):
        ns=[0.0]*n
        for r,v in enumerate(st):
            if v==0: continue
            ns[0]+=v*(1-p)
            if r+1>=n: hit+=v*p
            else: ns[r+1]+=v*p
        st=ns
    return hit
import collections
ns=[44,46,48,50,52,55,58]
exp=collections.defaultdict(float); obs=collections.defaultdict(int)
expc=collections.defaultdict(lambda: collections.defaultdict(float)); obsc=collections.defaultdict(lambda: collections.defaultdict(int))
cache={}
def isprim(h):
    g=math.gcd(h['a'],h['d'])
    return all(g%p for p in GOOD)
tres=0
for u in U:
    K,s=u['K'],u['shift']
    sh=shape(K); MOD,res,ly,pinned=sh
    wbar=((pinned[-1]-58)+(pinned[-2]-58))/4
    d=K*D0
    W=u['conf']
    gp=goodpen(K)
    # class
    lnT=math.log((64*s+32+wbar)*MOD+57*d)
    cl='lnT<40' if lnT<40 else 'lnT 40-41' if lnT<41 else 'lnT 41-42' if lnT<42 else 'lnT>42'
    for i in range(4):
        b=(i+0.5)*16; a=(64*s+b+wbar)*MOD
        p=math.exp(sum(math.log(rho(a+k*d)) for k in range(0,58,6))/10)
        key=round(p,3)
        if key not in cache: cache[key]=[prun(key,n) for n in ns]
        for n,pr in zip(ns,cache[key]):
            exp[n]+=W/4*pr*gp; expc[cl][n]+=W/4*pr*gp
    for h in u['hitrec']:
        if not isprim(h): continue
        for n in ns:
            if h['n']>=n: obs[n]+=1; obsc[cl][n]+=1
    tres+=u['res']
print('units',len(U),'res %.3g'%tres)
print('n  expected(primitive, in-window)  observed(primitive)')
for n in ns: print(n,'%.3g'%exp[n],obs[n])
for cl in sorted(expc):
    print(cl,' '.join('%d:%.2f/%d'%(n,expc[cl][n],obsc[cl][n]) for n in ns))

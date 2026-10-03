from load import *
from model import *
import collections
E_='/workspace/ponder-this-sept/experiments/'
def prun(p,n,N=58):
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
def isprim(h):
    g=math.gcd(h['a'],h['d'])
    return all(g%p for p in GOOD)
def shape_for(u):
    for cap in (2e16,2e15,1.2e15,1.2e14):
        P.MODCAP=int(cap); sh=shape(u['K'])
        if sh and abs(sh[1]/u['res']-1)<0.01: return sh
    return None
cache={}
def run(files,ns,label):
    U=load(files); exp=collections.defaultdict(float); obs=collections.defaultdict(int); skipped=0; seen=set()
    for u in U:
        K,s=u['K'],u['shift']; sh=shape_for(u)
        if sh is None: skipped+=1; continue
        MOD,res,ly,pinned=sh
        wbar=((pinned[-1]-58)+(pinned[-2]-58))/4
        if 'c1' in u['file']: wbar=0
        d=K*D0; W=u['conf']; gp=goodpen(K)
        for i in range(4):
            b=(i+0.5)*16; a=(64*s+b+wbar)*MOD
            p=math.exp(sum(math.log(rho(a+k*d)) for k in range(0,58,6))/10)
            key=round(p,3)
            for n in ns:
                if (key,n) not in cache: cache[(key,n)]=prun(key,n)
                exp[n]+=W/4*cache[(key,n)]*gp
        for h in u['hitrec']:
            if not isprim(h): continue
            seen.add((h['a'],h['d']))
            for n in ns:
                if h['n']>=n: obs[n]+=1
    print(label,'units',len(U),'skipped',skipped)
    print('  '+'  '.join('%d: %.2f/%d'%(n,exp[n],obs[n]) for n in ns))
    return exp,obs
ns=[36,38,40,42,44,46,48,50,52,55,58]
a=run([E_+'2026-10-03-cuda/c1.jsonl'],ns,'c1 (report36)')
b=run([E_+'2026-10-03-cuda/c2.jsonl',E_+'2026-10-03-cuda/p*.jsonl',E_+'2026-10-03-cuda/s0*.jsonl'],ns,'c2,p1,p2,s0')
c=run([E_+'remote/r1_g*.jsonl'],ns,'remote')
print('combined n>=44:')
for n in ns:
    if n>=44: print(n,'%.2f'%(a[0][n]+b[0][n]+c[0][n]), a[1][n]+b[1][n]+c[1][n])

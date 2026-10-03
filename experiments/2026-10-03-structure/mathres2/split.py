import sys,math,collections,os
sys.path.insert(0,'/workspace/mathres')
os.environ['MODCAP']='2e16'
from load import *
from model import *
exec(open('cal.py').read().split('cache={}')[0].split("E_=")[0].split('from model import *')[1]) if 0 else None
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
cache={}
def feats(u,sh):
    K=u['K']; f=[]
    f.append('good|K' if goodpen(K)<1 else 'nogood')
    f.append('K%%2=%d'%(K%2)); f.append('K%%3=%d'%(min(K%3,1)))
    f.append('5|K' if K%5==0 else '5nK')
    f.append('shift'+('0' if u['shift']==0 else '1-2' if u['shift']<3 else '3-7' if u['shift']<8 else '8+'))
    f.append('K<3e3' if K<3000 else 'K<1e4' if K<10000 else 'K<3e4' if K<30000 else 'K>=3e4')
    f.append('pin0=%d'%sh[3][0]); f.append('pinlast=%d'%sh[3][-1])
    f.append('file='+u['file'].split('/')[-1][:2])
    return f
U=load([E_+'remote/r*.jsonl'])
exp=collections.defaultdict(float); obs=collections.defaultdict(int); seen=set()
for u in U:
    K,s=u['K'],u['shift']
    if (K,s) in seen: continue
    seen.add((K,s)); P.MODCAP=int(2e16); sh=shape(K)
    MOD,r,ly,pinned=sh
    wbar=((pinned[-1]-58)+(pinned[-2]-58))/4
    d=K*D0; W=u['conf']; gp=goodpen(K); e=0
    for i in range(4):
        b=(i+0.5)*16; a=(64*s+b+wbar)*MOD
        p=round(math.exp(sum(math.log(rho(a+k*d)) for k in range(0,58,6))/10),3)
        if p not in cache: cache[p]=prun(p,44)
        e+=W/4*cache[p]*gp
    o=0
    for h in u['hitrec']:
        g=math.gcd(h['a'],h['d'])
        if all(g%p for p in GOOD): o+=1
    for f in feats(u,sh): exp[f]+=e; obs[f]+=o
for f in sorted(exp): print('%-12s exp %.1f obs %d  ratio %.2f  z %.1f'%(f,exp[f],obs[f],obs[f]/exp[f],(obs[f]-exp[f])/exp[f]**.5))

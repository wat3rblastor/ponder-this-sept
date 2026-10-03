import sys,math,os,collections,glob,json
os.environ['MODCAP']='2e16'
sys.path.insert(0,'/workspace/mathres')
from model import *
E_='/workspace/ponder-this-sept/experiments/'
def tail(p,N=58):
    # P(longest run >= n) for n=0..58
    out=[1.0]
    for n in range(1,N+1):
        st=[0.0]*n; st[0]=1.0; hit=0.0
        for _ in range(N):
            ns=[0.0]*n
            for r,v in enumerate(st):
                if v==0: continue
                ns[0]+=v*(1-p)
                if r+1>=n: hit+=v*p
                else: ns[r+1]+=v*p
            st=ns
        out.append(hit)
    return out
cache={}
obs=[0.0]*64; exp=[0.0]*64; exp2=[0.0]*64; nu=0; W=0
units={}
for f in glob.glob(E_+'2026-10-03-cuda/*.jsonl'):
    reduced = 'c1' in f
    for line in open(f):
        try: j=json.loads(line)
        except: continue
        if j.get('hist'): units[(j['K'],j['shift'],reduced)]=j['runs']
for (K,s,reduced),runs in units.items():
    if K%59==0: continue          # windows isolated only when 59 is pinned
    sh=None
    for cap in (2e16,2e15,1.2e15,1.2e14):
        P.MODCAP=int(cap); sh=shape(K)
        if sh and len(sh[3])==7 and sh[3][0]==59: break
    if not sh or len(sh[3])!=7 or goodpen(K)<1: continue
    MOD,r,ly,pinned=sh; d=K*D0
    wbar=0 if reduced else ((pinned[-1]-58)+(pinned[-2]-58))/4
    n_w=sum(runs); nu+=1; W+=n_w
    for i in range(4):
        a=(64*s+(i+.5)*16+wbar)*MOD
        p=round(math.exp(sum(math.log(rho(a+k*d)) for k in range(0,58,6))/10),3)
        for pp,ex in ((p,exp),(round(p*0.995,4),exp2)):
            if pp not in cache: cache[pp]=tail(pp)
            for n in range(59): ex[n]+=n_w/4*cache[pp][n]
    c=0
    for n in range(63,-1,-1):
        c+=runs[n]; obs[n]+=c
print('units',nu,'windows',W)
for n in list(range(4,48,3)):
    print('n>=%2d obs %8d exp %10.1f ratio %.3f  z %+.1f | exp(p*0.995) ratio %.3f'%(n,obs[n],exp[n],obs[n]/exp[n],(obs[n]-exp[n])/exp[n]**.5,obs[n]/exp2[n]))

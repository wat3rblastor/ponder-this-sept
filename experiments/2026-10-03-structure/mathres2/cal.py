import sys,math,collections,os
sys.path.insert(0,'/workspace/mathres')
os.environ['MODCAP']='2e16'
from load import *
from model import *
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
def prim(h):
    g=math.gcd(h['a'],h['d']); 
    return all(g%p for p in GOOD), (h['a']//g, h['d']//g)
def shape_for(u):
    for cap in (2e16,2e15,1.2e15,1.2e14):
        P.MODCAP=int(cap); sh=shape(u['K'])
        if sh and abs(sh[1]/u['res']-1)<0.01: return sh
    return None
cache={}
ns=[44,46,48,50,52,54,55,56,57,58]
def run(files,label):
    U=load(files); exp=collections.defaultdict(float); obs=collections.defaultdict(int); uniq=collections.defaultdict(set); skipped=0; res=0;gs=0
    seenu=set()
    for u in U:
        K,s=u['K'],u['shift']
        if (K,s) in seenu: continue
        seenu.add((K,s))
        sh=shape_for(u)
        if sh is None: skipped+=1; continue
        MOD,r,ly,pinned=sh; res+=u['res']; gs+=u.get('gpu_s',0)
        wbar=((pinned[-1]-58)+(pinned[-2]-58))/4
        d=K*D0; W=u['conf']; gp=goodpen(K)
        for i in range(4):
            b=(i+0.5)*16; a=(64*s+b+wbar)*MOD
            p=math.exp(sum(math.log(rho(a+k*d)) for k in range(0,58,6))/10)
            key=round(p,3)
            for n in ns:
                if (key,n) not in cache: cache[(key,n)]=prun(key,n)
                exp[n]+=W/4*cache[(key,n)]*gp
        for h in u['hitrec']:
            ok,pr=prim(h)
            for n in ns:
                if h['n']>=n:
                    uniq[n].add(pr)
                    if ok: obs[n]+=1
    print(label,'units',len(seenu),'skipped',skipped,'res %.3g gpu_s %.0f'%(res,gs))
    for n in ns: print('  n>=%d: exp %.2f  obs_prim %d  uniq_all %d'%(n,exp[n],obs[n],len(uniq[n])))
    return exp,obs,uniq
run([E_+'remote/r*.jsonl'],'remote all')
for lab in ['r1','r2','r3','r4','r5']:
    run([E_+'remote/%s_*.jsonl'%lab],lab)
run([E_+'2026-10-03-cuda/c2.jsonl',E_+'2026-10-03-cuda/p*.jsonl',E_+'2026-10-03-cuda/s0*.jsonl'],'cuda c2,p,s0')

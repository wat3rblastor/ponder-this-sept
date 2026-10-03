import sys,math,collections,os,glob,json
os.environ['MODCAP']='2e16'
sys.path.insert(0,'/workspace/ponder-this-sept/tools')
import plan_units as N
exec(open('cal2.py').read().split("NS=[44")[0])
def fac(n):
    f=[];p=2
    while p*p<=n:
        while n%p==0: f.append(p); n//=p
        p+=1
    if n>1: f.append(n)
    return f
def FF(K,MOD,q):
    co=MOD//q; return (K*D0%q)*pow(co%q,-1,q)%q/q
def pins(K,MOD):
    return [q for q in N.BAD if MOD%q==0]
def F2(K,MOD): return FF(K,MOD,pins(K,MOD)[-1])
def F1(K,MOD): return FF(K,MOD,pins(K,MOD)[-2])
exp=collections.defaultdict(float); obs=collections.defaultdict(int); seen=set()
rows=[]
for f in sorted(glob.glob(E_+'remote/r[1-5]_*.jsonl')):
    b=f.split('/')[-1][:2]; pend=[]; lines=open(f).read().split('\n'); nun=sum(1 for l in lines if 'covered' in l); iu=0
    for line in lines:
        if not line: continue
        j=json.loads(line)
        if j.get('hit'): pend.append(j); continue
        if 'covered' not in j: continue
        K,s=j['K'],j['shift']; hs=[h for h in pend if h['K']==K]; pend=[]; iu+=1
        if (K,s) in seen: continue
        seen.add((K,s))
        sh=N.unit_shape(K); MOD,res,ly,wbar,lg=sh
        d=K*D0; W=j['conf']; e=0
        for bb in (8,24,40,56):
            a=(64*s+bb+wbar)*MOD
            e+=W/4*pwin(a,d,44)*goodpen(K,math.log(a+29*d))
        o=0
        for h in hs:
            g=math.gcd(h['a'],h['d'])
            if not any(g%p==0 for p in GOOD): o+=1
        fs=set(fac(K)); bq=[q for q in fs if q%3==2 and q>58]
        T=(64*s+32+wbar)*MOD+57*d
        ft=[b, b+('.q%d'%min(4,1+4*iu//(nun+1))), 'nbad=%d'%len(bq), 'unit '+('std 4.67e9' if res<5e9 else 'mid <1e11' if res<1e11 else 'heavy'),
            '59|K' if K%59==0 else '59 pinned', 'gpu_s/res '+('fast' if j['res']/max(j['gpu_s'],.01)>2.2e10 else 'slow'),
            'wbar %02d-%02d'%(int(wbar//8)*8,int(wbar//8)*8+8), 'heavy wbar %02d'%(int(wbar//8)*8) if res>1e11 else 'light wbar %02d'%(int(wbar//16)*16), 'Kmod64=%d'%(K%64>31), 'f2 %.1f'%(int(F2(K,MOD)*5)/5), 'f1 %.1f'%(int(F1(K,MOD)*5)/5), 'shift=%d'%min(s,4),
            'T '+('<5e17' if T<5e17 else '<1e18' if T<1e18 else '<2e18' if T<2e18 else '>=2e18'),
            'cpu_s '+('<.05' if j['cpu_s']<.05 else '>=.05'), 'conf/surv '+('=1' if j['conf']==j['surv'] else '<1')]
        for k in ft: exp[k]+=e; obs[k]+=o
        exp['ALL']+=e; obs['ALL']+=o
for k in sorted(exp):
    if exp[k]>3: print('%-22s exp %6.1f obs %4d ratio %.2f z %+.1f'%(k,exp[k],obs[k],obs[k]/exp[k],(obs[k]-exp[k])/exp[k]**.5))

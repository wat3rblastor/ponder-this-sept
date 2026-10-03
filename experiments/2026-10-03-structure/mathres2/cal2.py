import sys,math,collections,os,glob,json
os.environ['MODCAP']='2e16'
sys.path.insert(0,'/workspace/ponder-this-sept/tools')
import plan_units as N
N.MODCAP=int(2e16); D0=N.D0
GOOD=[p for p in N.sieve(2000) if p%3==1]
E_='/workspace/ponder-this-sept/experiments/'
def rho(T): return 0.71992*(math.log(T)/39.144)**-0.485
def pwin(a,d,n):
    # exact P(longest run>=n) with position-dependent p_k, n>29: sum over start j of [fail at j-1]*prod
    lp=[math.log(rho(a+k*d)) for k in range(58)]
    cs=[0]
    for x in lp: cs.append(cs[-1]+x)
    tot=0
    for j in range(0,58-n+1):
        pr=math.exp(cs[j+n]-cs[j])
        tot+=pr*(1 if j==0 else 1-math.exp(lp[j-1]))
    return tot
def goodpen(K,lnT):
    f=1.0
    for p in GOOD:
        if p>K: break
        if K%p==0: f*=(1-1/p)*math.exp(-(58/(p-1))*0.485*math.log(p)/lnT)
    return f
NS=[44,46,48,50,52,55,58]
exp=collections.defaultdict(float); obs=collections.defaultdict(int); seen=set()
files=sorted(glob.glob(E_+'remote/r*.jsonl'))
for f in files:
    b=f.split('/')[-1][:2]; pend=[]
    for line in open(f):
        j=json.loads(line)
        if j.get('hit'): pend.append(j); continue
        if 'covered' not in j: continue
        K,s=j['K'],j['shift']; hs=[h for h in pend if h['K']==K]; pend=[]
        if (K,s) in seen: continue
        seen.add((K,s))
        sh=N.unit_shape(K); MOD,res,ly,wbar,lg=sh
        assert abs(res/j['res']-1)<0.01
        d=K*D0; W=j['conf']
        for bb in (8,24,40,56):
            a=(64*s+bb+wbar)*MOD
            gp=goodpen(K,math.log(a+29*d))
            for n in NS:
                e=W/4*pwin(a,d,n)*gp
                exp[(b,n)]+=e; exp[('all',n)]+=e
        for h in hs:
            g=math.gcd(h['a'],h['d'])
            if any(g%p==0 for p in GOOD): continue
            for n in NS:
                if h['n']>=n: obs[(b,n)]+=1; obs[('all',n)]+=1
for b in ['r1','r2','r3','r4','r5','r7','r8','r9','all']:
    if (b,44) in exp: print(b,'  '.join('%d: %.1f/%d'%(n,exp[(b,n)],obs[(b,n)]) for n in NS))
e,o=exp[('all',44)],obs[('all',44)]
print('all n>=44 ratio %.3f z %.1f ; E58 model total %.3f'%(o/e,(o-e)/e**.5,exp[('all',58)]))

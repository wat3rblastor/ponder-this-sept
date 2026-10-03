import sys,math,os,collections,glob,json
os.environ['MODCAP']='2e16'
sys.path.insert(0,'/workspace/ponder-this-sept/tools')
import plan_units as N
D0=N.D0; GOOD=[p for p in N.sieve(100000) if p%3==1]
E_='/workspace/ponder-this-sept/experiments/'
hits={}; done=set()
for f in glob.glob(E_+'*/*.jsonl'):
    for line in open(f):
        try: j=json.loads(line)
        except: continue
        if j.get('hit') and j['n']>=44: hits[(j['a'],j['d'])]=j['n']
        elif 'covered' in j: done.add((j['K'],j['shift']))
prims=collections.defaultdict(list)
for (a,d),n in hits.items():
    g=math.gcd(a,d); m=1
    for p in GOOD:
        while g%p==0: g//=p; m*=p
    prims[(a//m,d//m)].append((m,n))
np_=[(k,v) for k,v in prims.items() if any(m>1 for m,_ in v)]
print('unique run reports',len(hits),'unique primitive APs',len(prims),'of which seen via a rescaled copy',len(np_))
both=sum(1 for k,v in np_ if any(m==1 for m,_ in v))
print('  rescaled AND primitive form also reported directly:',both,' only via image (new find through image):',len(np_)-both)
cnt=collections.Counter()
for (a,d),v in np_:
    K=d//D0; sh=N.unit_shape(K); MOD=sh[0]; mm=a//MOD
    direct=any(m==1 for m,_ in v)
    cov=((K,mm//64) in done) or ((K,max(mm//64-1,0)) in done) or K<=1949 and mm<64
    cnt[(direct,cov)]+=1
    if max(n for _,n in v)>=50: print('   n=%d prim a=%d K=%d mult m/MOD=%.1f  copies x%s direct=%s'%(max(n for _,n in v),a,K,a/MOD,sorted(m for m,_ in v),direct))
print('(direct_found, primitive unit in done-set):',dict(cnt))

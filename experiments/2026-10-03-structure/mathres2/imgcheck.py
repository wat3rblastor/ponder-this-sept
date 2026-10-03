import sys,math,os,collections,glob,json
D0=382160924970; NT=58
def sieve(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
PR=sieve(10000); BAD=[p for p in PR if p%3==2]; GOODP=[p for p in PR if p%3==1 and p<400]
def comps(K,cap=2e16):
    d=K*D0; C=[(3,1,0,1),(2,1,0,1),(5,1,1,4)]; MOD=30
    for q in BAD:
        if q<=NT or d%q==0: continue
        if MOD*q>cap: break
        C.append((q,d%q,d%q,q-NT)); MOD*=q
    return C,MOD
def wclass(K,a):
    """return (MOD, w) for the class of a in unit K under the strided walk, or None if a not admissible"""
    C,MOD=comps(K); d=K*D0
    R0=0; ss=[]
    for (m,s,t,c) in C:
        co=MOD//m; e=co*pow(co%m,-1,m)%MOD
        R0=(R0+s*e)%MOD; ss.append(t*e%MOD)
    idx=[]
    for (m,s,t,c) in C:
        if c==1: idx.append(0); continue
        if m==5: i=(a-1)%5
        else: i=(a*pow(d%m,-1,m)-1)%m
        if i>=c: return None
        idx.append(i)
    n=len(C); in2=n-1; in1=n-2
    Rpre=R0
    for j in range(n):
        if j in (in1,in2): continue
        Rpre=(Rpre+idx[j]*ss[j])%MOD
    A=Rpre+idx[in1]*ss[in1]+idx[in2]*ss[in2]
    assert A%MOD==a%MOD,(K,a)
    return MOD,(A-a%MOD)//MOD
def clean(a0,d):
    for q in BAD:
        if q>NT and d%q:
            r=a0%q; dq=d%q
            # zero index k = -a0/d mod q
            k=(-r*pow(dq,-1,q))%q
            if k<NT: return False
        elif a0%q==0: return False
    return True
E_='/workspace/ponder-this-sept/experiments/'
done={}; hits=collections.defaultdict(list)
for f in glob.glob(E_+'remote/r*.jsonl')+glob.glob(E_+'2026-10-03-cuda/[ps]*.jsonl')+glob.glob(E_+'2026-10-03-cuda/c2.jsonl'):
    b=f.split('/')[-1][:2]
    for line in open(f):
        try: j=json.loads(line)
        except: continue
        if j.get('hit'): hits[(j['a'],j['d'])].append((j['n'],b))
        elif 'covered' in j: done[(j['K'],j['shift'])]=b
print('units',len(done),'hit runs',len(hits))
exp=collections.Counter(); found=collections.Counter(); miss=[]
for (a,d),v in hits.items():
    n=max(x[0] for x in v); K=d//D0
    # all clean window starts overlapping the run by >=44
    for j in range(-14,n-44+15):
        a0=a+j*d
        if a0<=0: continue
        lo=max(0,-j); hi=min(57,n-1-j)
        if hi-lo+1<44: continue
        if not clean(a0,d): continue
        # same-K coverage and rescaled coverage
        for p in [1]+GOODP:
            K2=K*p; a2=a0*p
            if K2*D0*57+a2>=1.8e19: break
            r=wclass(K2,a2)
            if r is None: continue
            MOD,w=r; m=a2//MOD-w
            if m<0: continue
            s=m//64
            if (K2,s) in done:
                b=done[(K2,s)]; key=(b,'self' if p==1 else 'image')
                exp[key]+=1
                if (a*p,d*p) in hits: found[key]+=1
                else: miss.append((b,p,K2,s,a*p,d*p,n))
for k in sorted(exp): print(k,'expected reports',exp[k],'found',found[k])
print('missing examples:',miss[:12])

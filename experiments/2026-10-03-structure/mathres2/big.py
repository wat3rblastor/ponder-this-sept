import sys, math, time, numpy as np, flint
from multiprocessing import Pool
def sieve(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
PR=sieve(10000); BAD=[p for p in PR if p%3==2]
def rho(T): return 0.71992*(math.log(T)/39.144)**-0.485
def is_loe(t):
    if t<=1: return t>=0
    for p,e in flint.fmpz(t).factor():
        if int(p)%3==2 and e&1: return False
    return True
NMAX=30
G={}
def work(chunk):
    kind,d,par=G['fam']; out=[]
    for x in chunk:
        r=0
        for k in range(NMAX):
            if kind=='a': t=x+k*d
            elif kind=='K': t=par+k*x*d
            else: t=x*x+k*d
            if is_loe(t): r+=1
            else: break
        out.append(r)
    return out
def prod_bad(y): 
    r=3
    for p in BAD:
        if p<=y: r*=p
    return r
def run(name,kind,d,A,par=None,n0=10,free=()):
    t0=time.time()
    alive=np.ones(A,dtype=bool); alive[0]=False
    idx=None
    if kind=='a':
        alive[0::3]=False; alive[2::3]=False
        for p in BAD:
            if d%p==0: alive[0::p]=False
    elif kind=='s':
        m=par; alive[0::3]=False
        for p in BAD:
            if m%p==0: alive[0::p]=False
    cand=int(alive.sum()); lsieve=0.0
    for q in BAD:
        if d%q==0: continue
        kills=set()
        if kind=='a':
            for k in range(n0): kills.add((-k*d)%q)
        elif kind=='K':
            for k in range(1,n0): kills.add((-pow(k*d*par if False else k*d,-1,q)*par)%q)
        else:
            sq={}
            for r in range(q): sq.setdefault(r*r%q,[]).append(r)
            for k in range(1,n0):
                for r in sq.get((-k*d)%q,[]):
                    if r: kills.add(r)
        for r in kills: alive[r::q]=False
        lsieve+=math.log(1-len(kills)/q)
    surv=np.nonzero(alive)[0].tolist(); del alive
    G['fam']=(kind,d,par)
    ch=[surv[i:i+200] for i in range(0,len(surv),200)]
    with Pool(28) as P: res=[r for c in P.map(work,ch) for r in c]
    hist=[0]*(NMAX+1)
    for r in res: hist[r]+=1
    tail=lambda n: sum(hist[n:])
    # size of terms
    if kind=='a': Tk=lambda k: A/2+k*d+1
    elif kind=='K': Tk=lambda k: par+k*(A/2)*d
    else: Tk=lambda k: (A/2)**2+k*d
    # measured per-term rho among non-free tested terms with index<n0 (clean to 1e4)
    tested=passed=0; exp_pass=0.0
    for r in res:
        for k in range(min(r+1,n0)):
            if k in free: continue
            tested+=1; passed+=(k<r); exp_pass+=rho(Tk(max(k,1)))
    print('== %s  d=%.3e (%d digits)  T~%.2e  candidates=%d  clean windows(n0=%d)=%d  [%.0fs]'%(name,d,len(str(d)),Tk(n0),cand,n0,len(surv),time.time()-t0))
    print('   sieve survival: measured %.4e  model %.4e'%(len(surv)/cand,math.exp(lsieve)))
    print('   per-term rho (clean to 1e4): measured %.4f +- %.4f   model %.4f   ratio %.3f'%(passed/tested,(0.25/tested)**.5,exp_pass/tested,passed/exp_pass))
    scale=1e9/cand
    for n in (10,15,20,25,30):
        # model: generic independent terms
        ls=0.0
        for q in BAD:
            if d%q==0: continue
            if kind=='a': ls+=math.log(1-n/q)
            elif kind=='K': ls+=math.log(1-(n-1)/q)
            else:
                c=0
                for k in range(1,n):
                    if pow(-k*d%q,(q-1)//2,q)==1: c+=2
                ls+=math.log(max(1e-9,1-c/q))
        nf=sum(1 for k in range(n) if k in free)
        lr=sum(math.log(rho(Tk(k))) for k in range(n) if k not in free)
        mod=cand*math.exp(ls+lr)
        gen=cand*math.exp(sum(math.log(1-n/q) for q in BAD if d%q)+sum(math.log(rho(Tk(max(k,1)))) for k in range(n)))
        print('   runs>=%2d: observed %6d  model %9.2f  (per 1e9 cand: obs %.3g model %.3g; generic-AP model at this d,T: %.3g)'%(n,tail(n),mod,tail(n)*scale,mod*scale,gen*scale))
    sys.stdout.flush()
if __name__=='__main__':
    which=sys.argv[1]; A=int(float(sys.argv[2]))
    if which=='F1':
        for y in (53,59,71,113,200): run('F1 y=%d: d=3*prod(bad<=y), a<%.0e'%(y,A),'a',prod_bad(y),A)
    elif which=='F2':
        for y in (71,113): run('F2 y=%d: a=1, d=K*3*prod(bad<=y), K<%.0e'%(y,A),'K',prod_bad(y),A,par=1,free=(0,))
        run('F2 y=71: a=7, K<%.0e'%A,'K',prod_bad(71),A,par=7,free=(0,))
    elif which=='F3':
        M0=prod_bad(53)//3
        fr=(0,1,4,9,16,25)
        run('F3: a=s^2, d=3*M0^2, M0=prod(bad<=53), s<%.0e'%A,'s',3*M0*M0,A,par=M0,free=fr)
        M1=prod_bad(29)//3
        run('F3b: a=s^2, d=3*M1^2, M1=prod(bad<=29) (41,47,53 not in d), s<%.0e'%A,'s',3*M1*M1,A,par=M1,free=fr)
        run('F3c: a=s^2, d=3*(59*M0)^2, s<%.0e'%A,'s',3*(59*M0)**2,A,par=59*M0,free=fr)

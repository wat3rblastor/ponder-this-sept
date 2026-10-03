import sys, math, signal
from math import gcd, isqrt
from multiprocessing import Pool
P=sorted({j*(3*j-1)//2 for j in range(-60,61)})
def status(args):
    j,T,tmo=args
    import cypari2
    pari=cypari2.Pari(); pari.allocatemem(10**9,silent=True) if False else None
    f=pari.factor(T,10**6)
    bad=[];rest=1
    for p,e in zip(f[0],f[1]):
        p=int(p);e=int(e)
        if p<10**6 or pari.isprime(p):
            if p%3==2 and e&1: bad.append(p)
        else: rest*=p**e
    if bad: return (j,False,bad,len(str(T)))
    if rest==1: return (j,True,[],len(str(T)))
    def h(*a): raise TimeoutError
    signal.signal(signal.SIGALRM,h); signal.alarm(tmo)
    try:
        f=pari.factor(rest); signal.alarm(0)
        for p,e in zip(f[0],f[1]):
            p=int(p);e=int(e)
            if p%3==2 and e&1: bad.append(p)
        return (j,not bad,bad,len(str(T)))
    except BaseException as ex:
        signal.alarm(0)
        return (j,None,['unfactored c%d'%len(str(rest))],len(str(T)))
def run(X,m,e1,e2,jlo,jhi,tmo=120,procs=14):
    auto=set(b+p for b in (0,e1,e2) for p in P)
    jobs=[]; st={}
    for j in range(jlo,jhi+1):
        T=3*X*X+m*m*(24*j+1)
        if T<=0: continue
        if j in auto: st[j]=(True,['auto'])
        else: jobs.append((j,T,tmo))
    with Pool(procs) as pool:
        for j,ok,bad,dg in pool.imap_unordered(status,jobs): st[j]=(ok,bad)
    return st,auto
if __name__=="__main__":
    X,m,e1,e2=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
    jlo,jhi=int(sys.argv[5]),int(sys.argv[6]); tmo=int(sys.argv[7])
    print('digits of m',len(str(m)),'term digits',len(str(3*X*X+m*m*24*57)),'m mod 41,47,53:',[m%q==0 for q in (41,47,53)], 'm even',m%2==0)
    st,auto=run(X,m,e1,e2,jlo,jhi,tmo)
    s=''.join('A' if j in auto else ('L' if st[j][0] else ('?' if st[j][0] is None else '.')) for j in sorted(st))
    print('j from',min(st)); print(s)
    for j in sorted(st):
        if j not in auto: print(j,st[j][0],st[j][1][:4])

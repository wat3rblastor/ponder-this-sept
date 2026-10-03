import sys, itertools, math
from math import gcd, isqrt
import cypari2
pari=cypari2.Pari(); pari.allocatemem(2*10**9, silent=True) if hasattr(pari,'allocatemem') else None
Q=5*11*17*23*29
P=sorted({j*(3*j-1)//2 for j in range(-60,61)})
SMALLBAD=[p for p in range(2,20000) if p%3==2 and all(p%d for d in range(2,isqrt(p)+1))]
def loesch(n, full=True):
    """return (is_loeschian, list of odd-exp bad primes found)"""
    if n<=0: return (n==0, [])
    bad=[]
    for p in SMALLBAD:
        if n%p==0:
            e=0
            while n%p==0: n//=p; e+=1
            if e&1: bad.append(p)
    if bad: return (False,bad)
    if n==1: return (True,[])
    f=pari.factor(n)
    for p,e in zip(f[0],f[1]):
        p=int(p); e=int(e)
        if p%3==2 and e&1: bad.append(p)
    return (not bad, bad)
def double(X,m,Y,Z,M,N):
    X2=X**4-M*N*m**4; m2=2*X*m*Y*Z
    Y2=X**4+2*M*X*X*m*m+M*N*m**4; Z2=X**4+2*N*X*X*m*m+M*N*m**4
    g=gcd(X2,m2); X2//=g;m2//=g;Y2//=g;Z2//=g
    assert X2*X2+M*m2*m2==Y2*Y2 and X2*X2+N*m2*m2==Z2*Z2
    return abs(X2),abs(m2),abs(Y2),abs(Z2)
def analyse(X,m,e1,e2,jlo=-80,jhi=140,verbose=False):
    """progression T(j)=3X^2+m^2(24j+1); bases at j=0,e1,e2"""
    auto=set()
    for b in (0,e1,e2):
        for p in P: auto.add(b+p)
    st={}
    for j in range(jlo,jhi+1):
        T=3*X*X+m*m*(24*j+1)
        if T<=0: continue
        ok,bad=loesch(T)
        st[j]=(ok,bad)
        if j in auto: assert ok,(X,m,e1,e2,j)
    # longest run
    best=(0,None); run=0
    for j in range(jlo,jhi+2):
        if j in st and st[j][0]:
            run+=1
            if run>best[0]: best=(run,j-run+1)
        else: run=0
    # best window of 57 by loeschian count
    bw=(0,None,0)
    for s in range(jlo,jhi-56):
        if all((s+k) in st for k in range(57)):
            c=sum(st[s+k][0] for k in range(57)); a=sum((s+k) in auto for k in range(57))
            if c>bw[0]: bw=(c,s,a)
    return best,bw,st,auto
if __name__=="__main__":
    fn=sys.argv[1]; maxdig=int(sys.argv[2]) if len(sys.argv)>2 else 60
    needQ=int(sys.argv[3]) if len(sys.argv)>3 else Q
    seen=set(); out=[]
    for line in open(fn):
        v=list(map(int,line.split())); X,m,E=v[0],v[1],v[2:]
        if m==1 and X>60: continue
        bases=[0]+E
        if len(bases)>8: bases=bases[:8]
        for b0,b1,b2 in itertools.combinations(bases,3):
            if b2-b0>70: continue
            Xr=isqrt(X*X+8*b0*m*m); Y=isqrt(X*X+8*b1*m*m); Z=isqrt(X*X+8*b2*m*m)
            M=8*(b1-b0); N=8*(b2-b0)
            if (Xr*m*Y*Z)%needQ: continue
            X2,m2,Y2,Z2=double(Xr,m,Y,Z,M,N)
            if m2%needQ: continue
            key=(X2,m2,M,N)
            if key in seen: continue
            seen.add(key)
            dig=len(str(3*X2*X2+m2*m2*24*57))
            if dig>maxdig: continue
            out.append((dig,X2,m2,b1-b0,b2-b0,(X,m,b0,b1,b2)))
    out.sort()
    print(len(out),'candidate doubled points with Q|m', file=sys.stderr)
    for dig,X2,m2,e1,e2,src in out:
        best,bw,st,auto=analyse(X2,m2,e1,e2)
        print(f"dig={dig} e=({e1},{e2}) src={src} X={X2} m={m2} longest_run={best} best57(count,start,auto)={bw} m%41,47,53={[m2%q==0 for q in (41,47,53)]}",flush=True)

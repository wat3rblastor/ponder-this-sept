# Family: a = 3x^2 + m^2, d = 24 m^2, m = 623645*K, gcd(x,m)=1.
# i in A  <=>  x^2 + 8 i m^2 is a square ; automatic k = i + pentagonal.
import sys, math, time, itertools
from multiprocessing import Pool
sys.path.insert(0,'/workspace/ponder-this-sept/src')
from loeschian import is_loeschian, factorize
M0=5*11*17*23*29
BAD=[p for p in range(2,3000) if p%3==2 and all(p%q for q in range(2,int(p**.5)+1))]
def issq(n):
    if n<0: return False
    r=math.isqrt(n); return r*r==n
P=[p for p in range(0,600) if issq(24*p+1)]
def divisors(f):
    ds=[1]
    for p,e in f.items():
        ds=[d*p**j for d in ds for j in range(e+1)]
    return ds
def fastloe(n):
    for q in BAD:
        if n%q==0:
            e=0
            while n%q==0: n//=q; e+=1
            if e&1: return False
    return is_loeschian(n)
def runs(a,d):
    fl=[fastloe(a+k*d) for k in range(58)]
    best=cur=0
    for f in fl:
        cur=cur+1 if f else 0; best=max(best,cur)
    return best,sum(fl)
def work(K):
    m=M0*K; fm=factorize(m); out=[]
    xs={}
    XMAX=int(sys.argv[3])*m if len(sys.argv)>3 else 300*m
    for i in range(-200,58):
        if i==0: continue
        f=dict((p,2*e) for p,e in fm.items())
        for p,e in factorize(2*abs(i)).items(): f[p]=f.get(p,0)+e
        n=2*abs(i)*m*m; r=math.isqrt(n)
        for f1 in divisors(f):
            if f1>r: continue
            g1=n//f1
            x=g1-f1 if i>0 else g1+f1
            if x<=0 or x>XMAX or math.gcd(x,m)!=1: continue
            xs.setdefault(x,set()).add(i)
    for x,A in xs.items():
        A=A|{0}
        cov={i+p for i in A for p in P if 0<=i+p<=57}
        if len(cov)>=int(sys.argv[4]) if len(sys.argv)>4 else len(cov)>=23:
            a=3*x*x+m*m; d=24*m*m
            best,cnt=runs(a,d)
            out.append((best,cnt,len(cov),sorted(A),K,x,a,d))
    return out
if __name__=="__main__":
    k0,k1=int(sys.argv[1]),int(sys.argv[2]); t=time.time()
    res=[]
    with Pool(10) as pool:
        for o in pool.imap_unordered(work,range(k0,k1)): res+=o
    from collections import Counter
    print("candidates",len(res),"time",time.time()-t)
    print("coverage hist",sorted(Counter(r[2] for r in res).items()))
    print("run hist",sorted(Counter(r[0] for r in res).items()))
    print("count hist",sorted(Counter(r[1] for r in res).items()))
    for r in sorted(res,reverse=True)[:8]: print(r)
    for r in sorted(res,key=lambda r:-r[2])[:5]: print("maxcov",r)

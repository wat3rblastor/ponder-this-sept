import sys, itertools, math
from multiprocessing import Pool
from sympy import factorint, primerange
BP=[p for p in primerange(30,3000) if p%3==2]
FREE={0,1,2,4,6,9,12,16,20,25,30,36,42,49,56}
ORDER=[k for k in range(58) if k not in FREE]
def lo(n):
    for p in BP:
        if n%p==0:
            e=0
            while n%p==0: n//=p; e+=1
            if e&1: return False
    if n%3!=1 and n%3!=0: return False
    return all(e%2==0 for p,e in factorint(n).items() if p%3==2)
def run(a,d):
    k=0
    while lo(a+k*d): k+=1
    return k
P=[5,11,17,23,29]
def work(mask):
    U0=1;V0=1
    for i,p in enumerate(P):
        if mask>>i&1: U0*=p
        else: V0*=p
    best=(0,0,0); hist={}
    for up in range(1,60):
        u=U0*up
        for vp in range(1,60):
            v=V0*vp
            if math.gcd(u,v)!=1 or v%3==0 or (u+v)%2==0: continue
            A=3*u*u+v*v
            if math.gcd(A,2*5*11*17*23*29)!=1: continue
            a=A*A; d=48*u*u*v*v
            r=run(a,d); hist[r]=hist.get(r,0)+1
            if r>best[0]: best=(r,a,d)
    return best,hist
if __name__=="__main__":
    with Pool(10) as pool:
        res=pool.map(work,range(32))
    H={}
    for b,h in res:
        for k,v in h.items(): H[k]=H.get(k,0)+v
    print(max(r[0] for r in res)); print(max(res)[0])
    print(sorted(H.items()))

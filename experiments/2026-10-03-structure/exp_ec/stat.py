import sys, itertools
from math import isqrt
from collections import Counter
c=Counter(); tot=0; best=[]
for line in open(sys.argv[1]):
    v=list(map(int,line.split())); X,m,E=v[0],v[1],v[2:]
    if m==1: continue
    bases=[0]+E
    for b0,b1,b2 in itertools.combinations(bases,3):
        Xr=isqrt(X*X+8*b0*m*m); Y=isqrt(X*X+8*b1*m*m); Z=isqrt(X*X+8*b2*m*m)
        pr=Xr*m*Y*Z; tot+=1
        ds=tuple(q for q in (5,11,17,23,29) if pr%q==0)
        c[len(ds)]+=1
        for q in ds: c['q%d'%q]+=1
        if len(ds)>=4: best.append((X,m,b0,b1,b2,ds))
print(tot,sorted(c.items(),key=str)); print(best[:40])

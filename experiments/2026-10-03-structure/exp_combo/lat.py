# L+S lattice count: s_k^2 | a+k d for k=0..56 (pairwise coprime primes s_k), D | d.
from math import log10
def primes(lo,n):
    out=[];p=lo
    while len(out)<n:
        if all(p%q for q in range(2,int(p**.5)+1)): out.append(p)
        p+=1
    return out
def crt(r1,m1,r2,m2):
    return (r1+m1*((r2-r1)*pow(m1,-1,m2)%m2))%(m1*m2), m1*m2
def gauss(u,v):
    n=lambda w:w[0]*w[0]+w[1]*w[1]
    while True:
        if n(u)>n(v): u,v=v,u
        mu=(2*(u[0]*v[0]+u[1]*v[1])+n(u))//(2*n(u))
        if mu==0: return u,v
        v=(v[0]-mu*u[0],v[1]-mu*u[1])
D=382160924970
for lo in (59,1000,10**6):
    s=primes(lo,57)
    c,M=0,1
    for k,p in enumerate(s): c,M=crt(c,M,(-k*D)%(p*p),p*p)
    u,v=gauss((c,D),(M,0))   # vectors (a, d) with d=D*e, a=c*e mod M
    a,d=abs(u[0]),abs(u[1])
    T=max(a,56*d)
    print(f"s_k~{lo}: log10(index)={log10(M*D):.1f}; shortest (a,d): log10 a={log10(a):.1f}, log10 d={log10(d):.1f}; "
          f"known square part per term {2*log10(s[28]):.1f} of {log10(T):.1f} digits = {200*log10(s[28])/log10(T):.1f}%")

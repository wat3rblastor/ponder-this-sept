import sys, pickle, itertools, math
from math import gcd, isqrt
import cypari2
pari=cypari2.Pari(); pari.allocatemem(4*10**9)
Qs=(5,11,17,23,29); Q=math.prod(Qs)
# difference pairs from good triples
pairs={}
for n in (57,58):
    for w,a,b,c in pickle.load(open(f'triples{n}.pkl','rb')):
        key=(b-a,c-a)
        if key not in pairs or pairs[key][0]<w: pairs[key]=(w,n,a)
print(len(pairs),'difference pairs',file=sys.stderr)
lo,hi=int(sys.argv[1]),int(sys.argv[2])
for (e1,e2),(w,n,a) in sorted(pairs.items())[lo:hi]:
    M,N=8*e1,8*e2
    E=pari.ellinit([0,M+N,0,M*N,0])
    try:
        rk=pari.ellrank(E,1)
    except Exception as ex:
        print(e1,e2,'ellrank fail',ex); continue
    gens=list(rk[3]); tors=pari.elltors(E)
    r=len(gens)
    if r==0:
        print(f"e=({e1},{e2}) cov={w} rank0 [{rk[0]},{rk[1]}] tors={tors[1]}",flush=True); continue
    hs=[float(pari.ellheight(E,g)) for g in gens]
    # torsion points
    tgens=list(tors[2]); tords=[int(o) for o in tors[1]]
    tpts=[]
    for cs in itertools.product(*[range(o) for o in tords]):
        T=pari('[0]')
        for c,g in zip(cs,tgens): T=pari.elladd(E,T,pari.ellmul(E,g,c))
        tpts.append(T)
    B={1:10,2:5,3:3}.get(r,2)
    best=None; found=[]; seen=set()
    ordinfo={}
    for cs in itertools.product(range(-B,B+1),repeat=r):
        if all(c==0 for c in cs): continue
        first=[c for c in cs if c][0]
        if first<0: continue
        R0=pari('[0]')
        for c,g in zip(cs,gens): R0=pari.elladd(E,R0,pari.ellmul(E,g,c))
        for T in tpts:
            R=pari.elladd(E,R0,T)
            D=pari.ellmul(E,R,2)
            if len(D)<2: continue
            x=D[0]; num=int(pari.numerator(x)); den=int(pari.denominator(x))
            X=isqrt(num); m=isqrt(den)
            assert X*X==num and m*m==den
            if (X,m) in seen: continue
            seen.add((X,m))
            dv=tuple(q for q in Qs if m%q==0)
            if len(dv)==5: found.append((len(str(m)),cs,X,m))
    found.sort()
    print(f"e=({e1},{e2}) cov={w}(n={n},a={a}) rank={r} [{rk[0]},{rk[1]}] tors={tors[1]} hts={[round(h,2) for h in hs]} npts={len(seen)} withQ={len(found)} best={[(f[0],f[1]) for f in found[:3]]}",flush=True)
    for f in found[:3]:
        print("PT",e1,e2,f[2],f[3],flush=True)

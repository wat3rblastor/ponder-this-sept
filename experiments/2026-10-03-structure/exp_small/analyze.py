import sys,glob,random,math,collections
from math import gcd
def primes(n):
    s=bytearray([1])*n; s[0:2]=b'\0\0'
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(range(i*i,n,i)))
    return [i for i in range(n) if s[i]]
PR=primes(200000)
def fac(m):
    f={}
    for p in PR:
        if p*p>m: break
        while m%p==0: f[p]=f.get(p,0)+1; m//=p
    if m>1: f[m]=f.get(m,0)+1
    return f
def isL(m): return all(e%2==0 for p,e in fac(m).items() if p%3==2)
bad=[p for p in PR if p%3==2 and p<200]
good=[p for p in PR if p%3==1 and p<200]
TOP=int(sys.argv[1]) if len(sys.argv)>1 else 50
print("n | N | #runs>=n (primitive) | top used | min last | frac q|d for bad q in (n/2,n) [q:frac_div/frac_square] | bad q>n in d | P(3^2|d) P(7|d) P(13|d) | mean#primeterms/n (baseline) | mean omega (baseline) | longest L")
for n in range(8,60,2):
    try: lines=open(f'res_{n}.txt').read().split('\n')
    except: continue
    N=int(float(lines[0].split('=')[1]))
    sols=[]
    for l in lines[1:]:
        if not l: continue
        a,d,L=map(int,l.split())
        if gcd(a,d)!=1: continue
        sols.append((a+(n-1)*d,a,d,L))
    sols.sort()
    top=sols[:TOP]
    if not top: continue
    opt=[q for q in bad if n/2<q<n]
    res=[]
    for q in opt:
        div=sum(1 for _,a,d,L in top if d%q==0)
        res.append(f"{q}:{div/len(top):.2f}/{1-div/len(top):.2f}")
    # weight expectation for square option: candidates ratio (2q-n)/q vs 1 -> expected frac_div = 1/(1+(2q-n))... report
    extra=[]
    for q in [x for x in bad if n<=x<3*n][:5]:
        div=sum(1 for _,a,d,L in top if d%q==0)
        extra.append(f"{q}:{div/len(top):.2f}")
    f9=sum(1 for _,a,d,L in top if d%9==0)/len(top)
    f7=sum(1 for _,a,d,L in top if d%7==0)/len(top)
    f13=sum(1 for _,a,d,L in top if d%13==0)/len(top)
    # prime terms & omega
    npr=0;om=0;cnt=0
    bpr=0;bom=0;bc=0
    random.seed(n)
    for _,a,d,L in top:
        for k in range(n):
            t=a+k*d; f=fac(t); cnt+=1
            if len(f)==1 and list(f.values())[0]==1: npr+=1
            om+=sum(f.values())
        # baseline: random Loeschian numbers in the same residue class a mod d? use t' = a' + K d with a' random coprime-admissible, same size, conditional on Loeschian
        tries=0
        while tries<n:
            t=random.randrange(a,a+n*d)//6*6+1
            if gcd(t,d)!=1: continue
            f=fac(t)
            if any(e%2 for p,e in f.items() if p%3==2): continue
            tries+=1; bc+=1
            if len(f)==1 and list(f.values())[0]==1: bpr+=1
            bom+=sum(f.values())
    print(f"{n} | {N:.3g} | {len(sols)} | {len(top)} | {top[0][0]} (a={top[0][1]}, d={top[0][2]}={fac(top[0][2])}) | {' '.join(res)} | {' '.join(extra)} | {f9:.2f} {f7:.2f} {f13:.2f} | {npr/cnt:.3f} ({bpr/bc:.3f}) | {om/cnt:.2f} ({bom/bc:.2f}) | {max(s[3] for s in sols)}")

import random, subprocess, collections
from math import isqrt, comb
from family import A_FORM,D_FORM,AUTO,val
NON=[k for k in range(58) if k not in AUTO]
PR=[2,5,11,17,23,29,41,47,53]
def Mp(pp):
    r=1
    for l in PR:
        if pp%l: r*=l
    return r
rng=random.Random(31337); Tmax=10**15; mem=[]
while len(mem)<80000:
    pp=rng.randrange(1,isqrt(Tmax//21)+1); rem=Tmax-21*pp*pp
    if rem<=0: continue
    U=isqrt(rem//7)
    if U<=6*pp: continue
    m=Mp(pp); c=(U-6*pp)//m
    if c<1: continue
    u=6*pp+m*(rng.randrange(c)+1); qq=u-12*pp
    mem.append((val(A_FORM,pp,qq),val(D_FORM,pp,qq)))
hist=collections.Counter(); npass=0
for i in range(0,len(mem),20000):
    ch=mem[i:i+20000]
    lines=[A+k*D for A,D in ch for k in range(58)]
    r=subprocess.run(['../../build/apsearch','--isl'],input="\n".join(map(str,lines)),capture_output=True,text=True)
    v=[int(x.split()[1]) for x in r.stdout.strip().split("\n")]
    for j in range(len(ch)):
        w=v[j*58:(j+1)*58]
        assert all(w[k] for k in AUTO)
        c=sum(w[k] for k in NON); hist[c]+=1; npass+=c
N=len(mem); rho=npass/(N*41)
print(f"variant B, T<=1e15, {N} members, rho={rho:.5f}, 0 automatic failures")
print("tail of the non-automatic pass count: observed vs binomial(41,rho)")
for c in range(28,42):
    e=N*comb(41,c)*rho**c*(1-rho)**(41-c)
    print(f"  {c:2d}/41  obs {hist[c]:7d}  exp {e:10.2f}  ratio {hist[c]/e if e>0 else 0:.3f}")

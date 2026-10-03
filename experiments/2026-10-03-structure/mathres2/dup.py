import sys,math,os,collections
os.environ['MODCAP']='2e16'
sys.path.insert(0,'/workspace/ponder-this-sept/tools'); sys.path.insert(0,'/workspace/mathres')
import plan_units as N
from load import load
D0=N.D0; GOOD=[p for p in N.sieve(100000) if p%3==1]
def goodfac(K):
    r=[]
    for p in GOOD:
        if p>K: break
        if K%p==0: r.append(p)
    return r
LT=40.5
def pen(p): return math.exp(-(58/(p-1))*0.485*math.log(p)/LT)
def boost(p): return (1-math.log(p)/LT)**(-0.485*58)
# ---- plan analysis
plan=sys.argv[1]
marks=[1e14,3e14,1e15,3e15,1e16,2e16,1e99]; mi=0
res=Eprim=Edup=wnp=0; nu=ng=0; resg=0
shapes={}
for line in open(plan):
    K,s=map(int,line.split()[:2])
    if K not in shapes: shapes[K]=(N.unit_shape(K),goodfac(K))
    sh,gf=shapes[K]
    if sh is None: continue
    MOD,r,ly,wbar,lgood=sh
    st,T=N.size_term(MOD,wbar,K,s)
    base=N.WCAL*64*N.RHO_REF58*math.exp(ly+st)*r     # generic-window yield of unit, no good-prime factor
    prim=1.0; tot=1.0; fp=1.0
    for p in gf:
        prim*=(1-1/p)*pen(p); tot*=((1-1/p)*pen(p)+boost(p)/p); fp*=(1-1/p)
    res+=r; nu+=1; Eprim+=base*prim; Edup+=base*(tot-prim); wnp+=r*(1-fp)
    if gf: ng+=1; resg+=r
    while res>=marks[mi]:
        print('at %.0e res: units %d, with good p|K: %.1f%% of units, %.1f%% of residues; non-primitive windows %.2f%% of all; E58 prim %.3f, dup(raw) %.3f = %.1f%% of raw'%(marks[mi],nu,100*ng/nu,100*resg/res,100*wnp/res,Eprim,Edup,100*Edup/(Edup+Eprim))); mi+=1
print('END res %.3g: units %d, good p|K units %.1f%%, residues %.1f%%; non-prim windows %.2f%%; E58 prim %.3f dup %.3f (%.1f%% of raw)'%(res,nu,100*ng/nu,100*resg/res,100*wnp/res,Eprim,Edup,100*Edup/(Edup+Eprim)))

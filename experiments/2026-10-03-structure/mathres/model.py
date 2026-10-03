import math, os, sys
sys.path.insert(0,'/workspace/ponder-this-sept/tools')
os.environ.setdefault('MODCAP','2e16')
import plan_units as P
D0=P.D0
# MC: rho for t=1 mod 3, coprime to all bad primes<=1e4
PTS=[(34.539,0.76510),(35.637,0.75343),(36.841,0.74129),(37.940,0.73097),(39.144,0.71992),(40.243,0.71014),(41.447,0.70031),(42.545,0.69141),(43.749,None)]
PTS=[p for p in PTS if p[1]]
def rho(T):
    x=math.log(math.log(T))
    xs=[math.log(p[0]) for p in PTS]; ys=[math.log(p[1]) for p in PTS]
    i=0
    while i<len(xs)-2 and x>xs[i+1]: i+=1
    return math.exp(ys[i]+(ys[i+1]-ys[i])*(x-xs[i])/(xs[i+1]-xs[i]))
GOOD=[p for p in P.sieve(200) if p%3==1]
def shape(K):
    """MOD,res,ly,pinned"""
    d=K*D0; MOD,res,pinned,divs=30,4,[],[]; full=False
    for q in P.BAD:
        if d%q==0: divs.append(q); continue
        if not full and MOD*q<=P.MODCAP: MOD*=q;res*=(q-58);pinned.append(q)
        else: full=True
    if len(pinned)<4: return None
    ly=P.BASE_UNPINNED
    for q in pinned: ly-=P.LOGF[q]
    for q in divs: ly+=math.log((q-1)/q)-P.LOGF[q]
    return MOD,res,ly,pinned
WCAL=0.752
def goodpen(K,lnT=40.5):
    f=1.0
    for p in GOOD:
        if K%p==0:
            f*=(1-1/p)*math.exp(-(58/(p-1))*0.485*math.log(p)/lnT)
    return f
def unit_E(K,s,n=58,nb=8,sh=None,wbar=None,primitive=True):
    """expected # of n-term windows (all n of first n terms... full 58-window prob) for unit"""
    sh=sh or shape(K)
    MOD,res,ly,pinned=sh
    W=WCAL*64*res*math.exp(ly)
    if wbar is None: wbar=((pinned[-1]-58)+(pinned[-2]-58))/4
    d=K*D0; tot=0
    for i in range(nb):
        b=(i+0.5)*64/nb
        a=(64*s+b+wbar)*MOD
        lp=0
        for k in range(0,58,3): lp+=math.log(rho(a+(k+1)*d))
        tot+=math.exp(lp*58/20)
    E=W*tot/nb
    if primitive: E*=goodpen(K)
    return E,W,(64*s+32+wbar)*MOD+57*d

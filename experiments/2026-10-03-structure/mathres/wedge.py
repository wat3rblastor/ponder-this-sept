import os,math,sys
os.environ['MODCAP']='2e16'
import importlib.util
spec=importlib.util.spec_from_file_location('np_','/workspace/mathres/plan_units.py'); N=importlib.util.module_from_spec(spec); spec.loader.exec_module(N)
D0=N.D0
def shape_cap(K,cap):
    N.MODCAP=int(cap); return N.unit_shape(K)
def wparams(K,MOD,pinned2):
    d=K*D0; r=[]
    for q in pinned2:
        co=MOD//q; f=(d%q)*pow(co%q,-1,q)%q/q
        r.append((q-59)*f)
    return r   # ranges of the two uniform parts
def pinned_of(K,cap):
    d=K*D0; MOD=30; pinned=[]
    for q in N.BAD:
        if d%q==0: continue
        if MOD*q<=cap: MOD*=q; pinned.append(q)
        else: break
    return MOD,pinned
def frac_w_gt(x,r1,r2):
    # P(u1+u2+0.5 > x), u_j~U[0,r_j]; numeric
    n=24; c=0
    for i in range(n):
        for j in range(n):
            if (i+.5)/n*r1+(j+.5)/n*r2+0.5>x: c+=1
    return c/(n*n)
CONST=N.WCAL*64*N.RHO_REF58
out=[]
KMAX=int(sys.argv[1]) if len(sys.argv)>1 else 30000
for K in range(1950,KMAX):
    big=shape_cap(K,2e16)
    if big is None: continue
    MODb,resb,lyb,wbarb,lgb=big
    _,pb=pinned_of(K,2e16)
    r1,r2=wparams(K,MODb,pb[-2:])
    if r1+r2<4: continue
    for cap in (2e15,1.2e14,1.2e13):
        sm=shape_cap(K,cap)
        if sm is None: continue
        MODs,ress,lys,wbars,lgs=sm
        if MODs*20>MODb: continue
        s=0
        while True:
            alo=(64*s+wbars)*MODs; amid=(64*s+32+wbars)*MODs
            if alo>(r1+r2)*MODb or s>400: break
            fr=frac_w_gt(amid/MODb,r1,r2)
            st,T=N.size_term(MODs,wbars,K,s)
            y=CONST*math.exp(lys+lgs+st)*fr
            out.append((y,ress,K,s,cap,T,fr))
            s+=1
        break   # only the largest cap that gives a smaller MOD
out.sort(reverse=True)
cum=0;E=0
marks=[1e13,3e13,1e14,3e14,1e15,3e15]; mi=0
for y,res,K,s,cap,T,fr in out:
    cum+=res;E+=y*res
    if mi<len(marks) and cum>=marks[mi]:
        print('wedge units: res %.3g  E58 %.4f  marginal y/res %.3g (K=%d s=%d cap=%.2g T=%.3g frac_new=%.2f)'%(cum,E,y,K,s,cap,T,fr)); mi+=1
print('total',len(out),'res %.3g E %.4f'%(cum,E))
import pickle; pickle.dump(out,open('wedge.pkl','wb'))

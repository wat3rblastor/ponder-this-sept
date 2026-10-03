import re,sys,collections,numpy as np,itertools
REF=1133661268029390
def grp(m):
    if abs(m-REF)<1e6: return 0
    if abs(m-2517112306980510)<1e6: return 1
    if m>=2.6e15: return 2
    return 3
def load(fn,tmin=None):
    rows=[]
    for l in open(fn):
        p=l.split(' ',2)
        if not p[2].startswith('K='): continue
        rows.append((float(p[0]),int(p[1])//2,grp(float(re.search(r'MOD=(\d+)',p[2])[1])),float(re.search(r'res=([\d.e+]+)',p[2])[1])))
    return rows
spans=[]
for fn in sys.argv[1:]:
    rows=load(fn); t0=min(r[0] for r in rows); t1=max(r[0] for r in rows)
    M=np.zeros((8,4))
    for t,g,c,r in rows: M[g,c]+=r
    spans.append((fn,t1-t0,M))
    print(fn,'span %.0fs'%(t1-t0)); 
    for g in range(8): print('  GPU%d'%g,' '.join('%.3e'%x for x in M[g]),' total res/s %.3e'%(M[g].sum()/(t1-t0)))
# fit: per span s, per GPU g: M[g]·cost = T_s * v_s   (cost_ref=1; v_s per span)
ns=len(spans)
def fit(sel):
    X=[];y=[]
    for si,(fn,T,M) in enumerate(spans):
        for g in range(8):
            if (si,g) not in sel: continue
            row=np.zeros(3+ns); row[:3]=M[g,1:]; row[3+si]=-T; X.append(row); y.append(-M[g,0])
    s,*_=np.linalg.lstsq(np.array(X),np.array(y),rcond=None); return s
allsel={(si,g) for si in range(ns) for g in range(8)}
s=fit(allsel)
print('cost c2517 %.3f big %.3f other %.3f ; speed per span (ref-equiv res/s per GPU):'%tuple(s[:3]),' '.join('%.3e'%v for v in s[3:]))
lo=[fit(allsel-{(si,g)}) for si,g in allsel]
lo=np.array(lo)
print('leave-one-out range: c2517 %.3f-%.3f big %.3f-%.3f other %.3f-%.3f'%(lo[:,0].min(),lo[:,0].max(),lo[:,1].min(),lo[:,1].max(),lo[:,2].min(),lo[:,2].max()))

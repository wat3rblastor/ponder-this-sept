import re,collections,numpy as np,sys
rows=collections.defaultdict(list)
for l in open('units_ts.txt'):
    p=l.split(' ',2)
    if not p[2].startswith('K='): continue
    s=p[2]; g=lambda k: float(re.search(k+r'=([\d.e+]+)',s)[1])
    rows[int(p[1])].append(dict(t=float(p[0]),res=g('res'),thr=g('thr'),MOD=g('MOD'),kern=g('kern')))
t0=max(min(r['t'] for r in v) for v in rows.values())+30; t1=min(max(r['t'] for r in v) for v in rows.values())-30
dt=0.05; n=int((t1-t0)/dt)
out=[]
for g in range(8):
    large=np.zeros(n,bool)
    for e in (2*g,2*g+1):
        for r in rows[e]:
            if r['thr']>=8e7:
                a=int((r['t']-r['kern']-t0)/dt); b=int((r['t']-t0)/dt)+1
                large[max(a,0):max(min(b,n),0)]=True
    TL=large.sum()*dt; TS=(~large).sum()*dt
    resL=resS=0.0; refS=0.0
    for e in (2*g,2*g+1):
        for r in rows[e]:
            if not (t0<=r['t']<t1): continue
            # attribute each unit's residues uniformly over its kernel interval
            a=r['t']-r['kern']; b=r['t']; 
            ia=int((max(a,t0)-t0)/dt); ib=int((min(b,t1)-t0)/dt)+1
            frac=(min(b,t1)-max(a,t0))/max(r['kern'],1e-6)
            if ib<=ia: continue
            seg=large[ia:ib]; fL=seg.mean()
            resL+=r['res']*fL*min(1,frac/1.0 if frac<1 else 1); resS+=r['res']*(1-fL)
    out.append((g,TS,TL,resS/TS if TS else 0,resL/TL if TL else 0))
    print('GPU%d  small-only phase %5.0fs  %.3e res/s   large-in-flight phase %5.0fs  %.3e res/s'%out[-1])
a=np.array(out); 
print('pooled: small-only %.3e res/s (T=%.0f), large phase %.3e res/s (T=%.0f), ratio L/S %.3f'%((a[:,3]*a[:,1]).sum()/a[:,1].sum(),a[:,1].sum(),(a[:,4]*a[:,2]).sum()/a[:,2].sum(),a[:,2].sum(),((a[:,4]*a[:,2]).sum()/a[:,2].sum())/((a[:,3]*a[:,1]).sum()/a[:,1].sum())))
print('per-GPU ratio L/S:',' '.join('%.2f'%(x[4]/x[3]) for x in out if x[3]))

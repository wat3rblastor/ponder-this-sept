import re,sys,collections,numpy as np,json
W=float(sys.argv[1]) if len(sys.argv)>1 else 60
REF=1133661268029390
def grp(m):
    if abs(m-REF)<1e6: return 'ref'
    if abs(m-2517112306980510)<1e6: return 'c2517'
    if m>=2.6e15: return 'big'
    return 'other'
G=['ref','c2517','big','other']
rows=[]
for l in open('units_ts.txt'):
    p=l.split(' ',2)
    if not p[2].startswith('K='): continue
    t=float(p[0]); e=int(p[1]); s=p[2]
    rows.append((t,e//2,grp(float(re.search(r'MOD=(\d+)',s)[1])),float(re.search(r'res=([\d.e+]+)',s)[1]),float(re.search(r'kern=([\d.]+)',s)[1])))
t0=min(r[0] for r in rows); t1=max(r[0] for r in rows)
nw=int((t1-t0)//W)
A=collections.defaultdict(lambda: np.zeros(4))
for t,g,c,res,k in rows:
    w=int((t-t0)//W)
    if w>=nw: continue
    A[(g,w)][G.index(c)]+=res
keys=sorted(A)
# unknowns: cost_c for c in big,c2517,other (ref=1) and speed s_g (res-equiv per second) for 8 GPUs
X=[];y=[]
for (g,w) in keys:
    v=A[(g,w)]; row=np.zeros(3+8); row[0:3]=v[1:4]; row[3+g]=-W; X.append(row); y.append(-v[0])
X=np.array(X);y=np.array(y)
sol,*_=np.linalg.lstsq(X,y,rcond=None)
res=y-X@sol
# bootstrap over windows
rng=np.random.default_rng(1); bs=[]
for _ in range(400):
    i=rng.integers(0,len(y),len(y)); s2,*_=np.linalg.lstsq(X[i],y[i],rcond=None); bs.append(s2[:3])
bs=np.array(bs)
tot=collections.Counter()
for t,g,c,r,k in rows: tot[c]+=r
TR=sum(tot.values())
print(f'span {t1-t0:.0f}s, windows {len(y)} ({W:.0f}s x 8 GPUs), units {len(rows)}')
print('group  res-share  cost(rel ref)  boot 5-95%')
print(f'ref    {100*tot["ref"]/TR:5.1f}%     1.000')
for i,c in enumerate(G[1:]):
    print(f'{c:6s} {100*tot[c]/TR:5.1f}%     {sol[i]:.3f}        {np.percentile(bs[:,i],5):.3f}-{np.percentile(bs[:,i],95):.3f}')
print('GPU speed (ref-equiv res/s):',' '.join(f'{v:.3e}' for v in sol[3:]))
cw=sum(tot[c]*([1]+list(sol[:3]))[G.index(c)] for c in G)
print(f'mix-weighted cost per res {cw/TR:.3f}; rms resid/W·s {np.sqrt(np.mean(res**2))/np.mean(sol[3:])/W:.3f}')
json.dump({'ref':1.0,'c2517':float(sol[0]),'big':float(sol[1]),'other':float(sol[2])},open(f'costs_w{int(W)}.json','w'))

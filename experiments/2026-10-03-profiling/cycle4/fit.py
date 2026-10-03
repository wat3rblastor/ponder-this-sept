import json,collections,statistics as st,sys
REF=1133661268029390
D='/workspace/ponder-this-sept/experiments/remote'
rows=[]
for i in range(16):
    for l in open(f'{D}/v2_s{i}.jsonl'):
        if '"kern_s"' not in l: continue
        try: j=json.loads(l)
        except: continue
        if j.get('kern_s',0)<=0: continue
        j['g']=i//2; j['e']=i; rows.append(j)
# MOD is not in jsonl: derive from planner shape
import os,sys; os.environ.setdefault('MODCAP','2e16'); sys.path.insert(0,'/workspace/ponder-this-sept/tools'); import plan_units as P
sh={}
for r in rows:
    if r['K'] not in sh: sh[r['K']]=P.unit_shape(r['K'])
    r['MOD']=sh[r['K']][0]
    if abs(sh[r['K']][1]*1.0-r['res'])/r['res']>1e-3: r['MOD']=None
bad=sum(r['MOD'] is None for r in rows)
print('units with kern_s',len(rows),'shape mismatches',bad)
rows=[r for r in rows if r['MOD']]
# per-GPU reference: median kern/res and pooled sum kern / sum res
ref={}
for g in range(8):
    x=[r['kern_s']/r['res'] for r in rows if r['g']==g and r['MOD']==REF]
    rs=[r for r in rows if r['g']==g and r['MOD']==REF]
    ref[g]=(st.median(x), sum(r['kern_s'] for r in rs)/sum(r['res'] for r in rs), len(x))
print('per-GPU ref ns/res (median, pooled, n):',{g:(round(v[0]*1e9,4),round(v[1]*1e9,4),v[2]) for g,v in ref.items()})
cls=collections.defaultdict(list)
for r in rows: cls[r['MOD']].append(r)
TR=sum(r['res'] for r in rows)
out={}
print('%-11s %6s %6s %7s %7s %7s %7s  %s'%('MOD','n','res%','cost_md','cost_pl','q25','q75','per-GPU pooled'))
for m,rs in sorted(cls.items(),key=lambda x:-sum(r['res'] for r in x[1])):
    rel=[ (r['kern_s']/r['res'])/ref[r['g']][0] for r in rs]
    pg={}
    for g in range(8):
        q=[r for r in rs if r['g']==g]
        if q: pg[g]=(sum(r['kern_s'] for r in q)/sum(r['res'] for r in q))/ref[g][1]
    pooled=st.median(pg.values())
    qs=st.quantiles(rel,n=4) if len(rel)>=4 else [min(rel),0,max(rel)]
    out[m]=dict(n=len(rs),share=sum(r['res'] for r in rs)/TR,cost_median=st.median(rel),cost_pooled=pooled,q25=qs[0],q75=qs[2],pergpu=pg)
    if len(out)<=25: print('%.4e %6d %5.1f%% %7.3f %7.3f %7.3f %7.3f  %s'%(m,len(rs),100*out[m]['share'],out[m]['cost_median'],pooled,qs[0],qs[2],' '.join('%.2f'%v for v in pg.values())))
json.dump({str(k):v for k,v in out.items()},open('cost_fit.json','w'),indent=0)

import sys,math,collections,json,glob,os
os.environ.setdefault("MODCAP","2e16")
sys.path.insert(0,'tools'); import plan_units as P
done=set()
for f in glob.glob('experiments/*/*.jsonl'):
    for l in open(f):
        if '"covered"' in l:
            try: j=json.loads(l); done.add((j['K'],j['shift']))
            except: pass
units=[tuple(map(int,l.split())) for l in open('experiments/remote_plan.txt')]
rem=[u for u in units if u not in done]
print('plan units',len(units),'remaining',len(rem))
shapes={}
cls=collections.defaultdict(lambda:[0,0.0,0.0])  # n, res, E
order=[]
for K,s in rem:
    if K not in shapes: shapes[K]=P.unit_shape(K)
    MOD,res,ly,wbar,lgood=shapes[K]
    st,T=P.size_term(MOD,wbar,K,s)
    y=math.exp(ly+lgood+st)
    c=cls[MOD]; c[0]+=1; c[1]+=res; c[2]+=y*res
    order.append((MOD,res,y))
R=sum(c[1] for c in cls.values()); E=sum(c[2] for c in cls.values())
for m,c in sorted(cls.items(),key=lambda x:-x[1][1])[:14]:
    print('MOD %.4e n=%7d res share %5.1f%%  E share %5.1f%%  yield/res rel %.3f'%(m,c[0],100*c[1]/R,100*c[2]/E,(c[2]/c[1])/(E/R)))
import pickle; pickle.dump(order,open('experiments/2026-10-03-profiling/cycle3/plan_order.pkl','wb'))

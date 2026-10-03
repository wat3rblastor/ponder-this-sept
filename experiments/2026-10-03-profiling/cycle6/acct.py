import re,sys,collections
fn=sys.argv[1]; ta=float(sys.argv[2]); tb=float(sys.argv[3])
E=collections.defaultdict(list)
for l in open(fn):
    p=l.split(' ',2)
    if not p[2].startswith('K='): continue
    t=float(p[0])
    if not (ta<=t<=tb): continue
    g=lambda k: float(re.search(k+r'=([\d.e+]+)',p[2])[1])
    E[int(p[1])].append(dict(t=t,res=g('res'),kern=g('kern'),gpu=g('gpu'),up=g('up'),wait=g('wait'),prep=g('prep'),cpu=g('cpu'),thr=g('thr')))
span=tb-ta
print('eng units  sum(kern)/wall  mean kern  mean up  mean wait  mean gpu  small-unit share of units')
tot=collections.Counter()
for e in sorted(E):
    v=E[e]; n=len(v)
    sk=sum(x['kern'] for x in v)
    print(f'{e:3d} {n:5d}   {sk/span:6.3f}   {sk/n*1e3:7.2f}ms {sum(x["up"] for x in v)/n*1e3:6.2f}ms {sum(x["wait"] for x in v)/n*1e3:7.2f}ms {sum(x["gpu"] for x in v)/n*1e3:7.1f}ms  {sum(x["thr"]<8e7 for x in v)/n:.2f}')
    tot['n']+=n; tot['k']+=sk
print(f'all: units/s {tot["n"]/span:.1f}, mean sum(kern)/wall per engine {tot["k"]/span/len(E):.3f}')

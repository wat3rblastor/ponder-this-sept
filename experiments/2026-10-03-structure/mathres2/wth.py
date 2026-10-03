import sys,math,collections,os,glob,json
os.environ['MODCAP']='2e16'
sys.path.insert(0,'/workspace/ponder-this-sept/tools')
import plan_units as N
N.MODCAP=int(2e16)
E_='/workspace/ponder-this-sept/experiments/'
agg=collections.defaultdict(lambda:[0,0,0,0])
for f in glob.glob(E_+'remote/r*.jsonl'):
    for line in open(f):
        j=json.loads(line)
        if 'covered' not in j: continue
        sh=N.unit_shape(j['K'])
        if not sh or abs(sh[1]/j['res']-1)>0.01: continue
        MOD,res,ly,wbar,lg=sh
        th=64*res*math.exp(ly)*0.7371
        K=j['K']
        key='59|K' if K%59==0 else ('other bad q|K' if any(K%q==0 for q in N.BAD) else 'no bad')
        for k in (key,'ALL', 'res>1e11' if res>1e11 else 'res<1e11'):
            a=agg[k]; a[0]+=j['surv']; a[1]+=j['conf']; a[2]+=th; a[3]+=1
for k,(s,c,t,n) in agg.items(): print('%-14s units %6d surv/theory %.4f conf/theory %.4f'%(k,n,s/t,c/t))

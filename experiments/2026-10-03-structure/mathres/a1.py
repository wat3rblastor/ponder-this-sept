from load import *
import os
os.environ['MODCAP']='2e16'
import plan_units as P
U=load(['/workspace/ponder-this-sept/experiments/remote/r1_g*.jsonl'])
print(len(U),'units',sum(u['res'] for u in U),'res', sum(len(u['hits']) for u in U),'hits')
# survivors model
import collections
rows=[]
for u in U:
    sh=P.unit_shape(u['K'])
    MOD,res,ly=sh
    T=64*MOD*(u['shift']+0.5)+57*u['K']*D0
    u['MOD']=MOD;u['T']=T;u['ly']=ly
    u['pred_surv']=res*64*math.exp(ly)   # windows
    u['resm']=res
# check res match
bad=[u for u in U if abs(u['res']/u['resm']-1)>0.01]
print('res mismatch',len(bad), [(u['K'],u['res'],u['resm']) for u in bad[:5]])
# surv = words with surviving bits; conf ~ windows
tot_s=sum(u['conf'] for u in U); tot_p=sum(u['pred_surv'] for u in U)
print('conf total',tot_s,'pred',tot_p,'ratio',tot_s/tot_p)
# by class of bad divisors
def badf(K):
    return tuple(q for q in P.BAD if K%q==0)
g=collections.defaultdict(lambda:[0,0,0,0,0])
for u in U:
    b=badf(u['K']); key=(len(b), u['MOD'])
    key=b if len(b)<=1 else ('multi',len(b))
    g[key][0]+=u['conf'];g[key][1]+=u['pred_surv'];g[key][2]+=1;g[key][3]+=u['res'];g[key][4]+=len(u['hits'])
for k in sorted(g,key=lambda k:-g[k][3])[:25]:
    v=g[k]; print(k,'units',v[2],'res %.3g'%v[3],'conf/pred %.3f'%(v[0]/v[1]),'hits44',v[4],'hits/conf %.2e'%(v[4]/v[0]))

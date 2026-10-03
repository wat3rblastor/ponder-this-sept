from load import *
import os,collections,time
E='/workspace/ponder-this-sept/experiments/'
U=load([E+'remote/r1_g*.jsonl'])
g=collections.defaultdict(lambda:[0,0,0,0])
tot=0
for u in U:
    k=float('%.2g'%u['res'])
    g[k][0]+=1;g[k][1]+=u['res'];g[k][2]+=u['gpu_s'];g[k][3]+=u['cpu_s']
    tot+=u['gpu_s']
print('sum gpu_s',tot,'units',len(U))
for k in sorted(g):
    v=g[k]
    if v[1]>1e12: print('res/unit %.2g n=%d totres=%.3g gpu_s=%.0f res/gpu_s=%.3g'%(k,v[0],v[1],v[2],v[1]/v[2]))

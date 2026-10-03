from load import *
import os; os.environ['MODCAP']='2e16'
import importlib.util
spec=importlib.util.spec_from_file_location('np_','/workspace/mathres/plan_units.py'); N=importlib.util.module_from_spec(spec); spec.loader.exec_module(N)
E_='/workspace/ponder-this-sept/experiments/'
U=load([E_+'remote/r*.jsonl'])
import statistics
xs=[];ys=[];lo=[];hi=[]
for u in U:
    if not u['hitrec']: continue
    sh=N.unit_shape(u['K']); MOD,res,ly,wbar,lg=sh
    if abs(res/u['res']-1)>0.01: continue
    q2,q1=None,None
    for h in u['hitrec']:
        m=h['a']/MOD-64*u['shift']
        xs.append(m); ys.append(wbar+32)
        # max possible w
        print(u['K'],u['shift'],'m=%.1f pred wbar+32=%.1f'%(m,wbar+32)) if len(xs)<=12 else None
print(len(xs),'mean obs',statistics.mean(xs),'mean pred',statistics.mean(ys))
import numpy as np
print('corr',np.corrcoef(xs,ys)[0,1])

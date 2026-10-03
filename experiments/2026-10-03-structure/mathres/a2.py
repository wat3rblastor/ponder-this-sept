from load import *
import os
os.environ['MODCAP']='2e16'
import plan_units as P
E='/workspace/ponder-this-sept/experiments/'
U=load([E+'remote/r1_g*.jsonl',E+'2026-10-03-cuda/p*.jsonl',E+'2026-10-03-cuda/c2.jsonl'])
def fac(n):
    f={};p=2
    while p*p<=n:
        while n%p==0: f[p]=f.get(p,0)+1; n//=p
        p+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
seen=set()
for u in U:
    for h in u['hitrec']:
        if (h['a'],h['d']) in seen: continue
        seen.add((h['a'],h['d']))
        if h['n']>=47: print(h['n'],h['K'],fac(h['K']),'shift',u['shift'],'a/d=%.1f'%(h['a']/h['d']),'T=%.3g'%(h['a']+57*h['d']), u['file'][-12:], 'res=%.3g'%u['res'])

import sys,math,os,collections
sys.path.insert(0,'/workspace/mathres')
os.environ['MODCAP']='2e16'
from load import *
from model import *
E_='/workspace/ponder-this-sept/experiments/'
for lab in ['r1','r2','r3','r4','r5']:
    U=load([E_+'remote/%s_*.jsonl'%lab]); sw=sp=0; n=0; ov=0; kinds=collections.Counter(); sh0=collections.Counter()
    for u in U:
        P.MODCAP=int(2e16); sh=shape(u['K'])
        if not sh or abs(sh[1]/u['res']-1)>0.01: continue
        MOD,res,ly,pinned=sh
        sp+=WCAL*64*res*math.exp(ly); sw+=u['conf']; n+=1; ov+=u.get('overflow',False)
        kinds[len(pinned)]+=1; sh0[min(u['shift'],6)]+=1
    print(lab,n,'conf',sw,'pred %.0f ratio %.3f'%(sp,sw/sp),'overflow',ov,dict(kinds),dict(sh0), 'meanK',sum(u['K'] for u in U)/len(U))

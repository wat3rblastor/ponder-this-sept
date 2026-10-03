import os,sys,math,json,re
os.environ.setdefault('MODCAP','2e16'); sys.path.insert(0,'/workspace/ponder-this-sept/tools'); import plan_units as P
CONST=P.WCAL*64*P.RHO_REF58; _sh={}
def E(K,s):
    if K not in _sh: _sh[K]=P.unit_shape(K)
    MOD,res,ly,wbar,lgood=_sh[K]; st,T=P.size_term(MOD,wbar,K,s)
    return CONST*math.exp(ly+lgood+st)*res, res

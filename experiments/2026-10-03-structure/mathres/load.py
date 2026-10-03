import json,glob,math,sys
sys.path.insert(0,'/workspace/ponder-this-sept/tools')
D0=382160924970
def load(globs):
    units=[];hits=[]
    for g in globs:
        for f in sorted(glob.glob(g)):
            pend=[]
            for line in open(f):
                try: j=json.loads(line)
                except: continue
                if j.get('hit'): pend.append(j)
                elif 'covered' in j:
                    j['hits']=[h['n'] for h in pend if h['K']==j['K']]
                    j['hitrec']=[h for h in pend if h['K']==j['K']]
                    j['file']=f; pend=[]; units.append(j)
    return units

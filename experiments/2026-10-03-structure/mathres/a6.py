import json,numpy as np
rows=[json.loads(l) for l in open('samp.jsonl')]
X=[];y=[]
for g in ['g0','g1','g2','g3','g4','g5','g6','g7']:
    for i in range(2,len(rows)):
        dt=rows[i]['t']-rows[i-1]['t']; v=rows[i][g]
        X.append([v[0],v[2],v[3]/1e10,v[5]/1e10]); y.append(dt)
X=np.array(X);y=np.array(y)
# aggregate in blocks of 4 intervals to reduce boundary noise
c,res,_,_=np.linalg.lstsq(X,y,rcond=None)
print('n',len(y),'coef: s per small unit %.4f, s per mid unit %.4f, s per 1e10 mid res %.4f, s per 1e10 big res %.4f'%tuple(c))
print('totals small units',X[:,0].sum(),'mid',X[:,1].sum(),'midres',X[:,2].sum(),'bigres',X[:,3].sum(),'time',y.sum())
print('small: res/s = %.3g ; big res/s = %.3g'%(4.67e9/c[0],1e10/c[3]))

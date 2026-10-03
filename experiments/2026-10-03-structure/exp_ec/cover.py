# coverage of window [0,n) by union of (a + pentagonal) for base sets
import sys, itertools
P=sorted({j*(3*j-1)//2 for j in range(-40,41)})
def mask(a,n):
    v=0
    for p in P:
        k=a+p
        if 0<=k<n: v|=1<<k
    return v
for n in (57,58):
    R=200
    ms={a:mask(a,n) for a in range(-R,n)}
    ms={a:v for a,v in ms.items() if v}
    keys=sorted(ms)
    best1=max(bin(v).count('1') for v in ms.values())
    best2=[];best3=[]
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            v=ms[a]|ms[b]; best2.append((v.bit_count(),a,b))
    best2.sort(reverse=True)
    print(n,'best1',best1,'best2',best2[:8])
    res=[]
    for i,a in enumerate(keys):
        for j in range(i+1,len(keys)):
            b=keys[j]; v=ms[a]|ms[b]
            for c in keys[j+1:]:
                w=(v|ms[c]).bit_count()
                if w>=30: res.append((w,a,b,c))
    res.sort(reverse=True)
    print(n,'triples>=30:',len(res),'top',res[:25])
    from collections import Counter
    print(Counter(r[0] for r in res))
    import pickle; pickle.dump(res,open(f'triples{n}.pkl','wb'))
    # 4 bases: extend top triples
    b4=[]
    for (w,a,b,c) in res[:3000]:
        v=ms[a]|ms[b]|ms[c]
        for d in keys:
            if d in (a,b,c): continue
            b4.append(((v|ms[d]).bit_count(),)+tuple(sorted((a,b,c,d))))
    b4=sorted(set(b4),reverse=True)
    print(n,'4-base (greedy ext) top',b4[:15])

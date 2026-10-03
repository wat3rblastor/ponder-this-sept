import sys,random,math
from collections import Counter
from multiprocessing import Pool
from fam import runs,M0
def w(seed):
    r=random.Random(seed); out=[]
    for _ in range(150):
        K=r.randrange(1,2000); d=6*M0*41*47*53*K
        while True:
            a=r.randrange(10**17,10**18)
            if a%3==1 and math.gcd(a,d)==1: break
        out.append(runs(a,d))
    return out
if __name__=="__main__":
    sys.argv=sys.argv[:1]
    res=[]
    with Pool(10) as p:
        for o in p.map(w,range(10)): res+=o
    print("n",len(res),"mean count",sum(c for _,c in res)/len(res))
    print("run hist",sorted(Counter(b for b,_ in res).items()))
    print("count hist",sorted(Counter(c for _,c in res).items()))

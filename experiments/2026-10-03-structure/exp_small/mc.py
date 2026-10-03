import random,subprocess,sys,math
bad=[2,5,11,17,23,29,41,47,53,59,71,83,89,101,107,113,131,137,149]
random.seed(1)
def run(X,Q,ns):
    bq=[q for q in bad if q<=Q]
    nums=[]
    while len(nums)<ns:
        t=random.randrange(X,2*X)//3*3+1
        if all(t%q for q in bq): nums.append(t)
    out=subprocess.run(['factor']+[str(t) for t in nums],capture_output=True,text=True).stdout
    ok=0;pr=0
    for l in out.strip().split('\n'):
        fs=[int(x) for x in l.split(':')[1].split()]
        c={}
        for f in fs: c[f]=c.get(f,0)+1
        if all(e%2==0 for p,e in c.items() if p%3==2): ok+=1
    return ok/ns
print("X \\ Q   "+"  ".join(f"{Q:>6}" for Q in (23,53,59,101,149)))
for e,ns in ((9,6000),(12,6000),(18,6000),(24,3000),(30,400)):
    print(f"1e{e:<4}", "  ".join(f"{run(10**e,Q,ns):6.3f}" for Q in (23,53,59,101,149)),flush=True)

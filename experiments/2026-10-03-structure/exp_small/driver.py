import subprocess,sys,os,time
bad=[5,11,17,23,29,41,47,53]
def D0(n):
    d=6
    for q in bad:
        if 2*q<=n: d*=q
    return d
start={8:20000}
N=20000
os.environ['OMP_NUM_THREADS']='10'
for n in range(int(sys.argv[1]),int(sys.argv[2])+1,2):
    N=int(N)
    while True:
        t=time.time()
        out=subprocess.run(['./apfind','find',str(N),str(n),str(D0(n))],capture_output=True,text=True).stdout
        lines=[l for l in out.split('\n') if l]
        print(n,N,D0(n),len(lines),round(time.time()-t,1),flush=True)
        if len(lines)>=60 or time.time()-t>150: break
        N=N*1.5
    open(f'res_{n}.txt','w').write(f'# N={N}\n'+'\n'.join(lines)+'\n')

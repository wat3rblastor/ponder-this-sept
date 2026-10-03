import os,sys,time,re
pid=int(sys.argv[1]); T=float(sys.argv[2]); log=sys.argv[3]
HZ=os.sysconf('SC_CLK_TCK')
def snap():
    d={}
    for t in os.listdir(f'/proc/{pid}/task'):
        try:
            f=open(f'/proc/{pid}/task/{t}/stat').read().rsplit(')',1)[1].split()
            d[int(t)]=(int(f[11])+int(f[12]))/HZ
        except Exception: pass
    return d
def nlines(): return sum(1 for l in open(log) if l.startswith('K='))
a=snap(); la=nlines(); t0=time.time(); time.sleep(T); b=snap(); lb=nlines(); dt=time.time()-t0
lines=[l for l in open(log) if l.startswith('K=')][la:lb]
bits=sum(int(re.search(r'bits=(\d+)',l)[1]) for l in lines)
res=sum(float(re.search(r'res=([\d.e+]+)',l)[1]) for l in lines)
tids=sorted(b)
dl={t:b[t]-a.get(t,0) for t in tids}
main=dl.get(pid,0)
# OMP team = the large block of threads with similar tids created late; classify by rank
others=sorted(((v,t) for t,v in dl.items() if t!=pid),reverse=True)
tot=sum(dl.values())
print(f'pid {pid} dt={dt:.1f}s threads={len(tids)} total cores={tot/dt:.2f} main={main/dt:.2f} units={len(lines)} bits={bits} res={res:.3e}')
print(f'core-s per unit {tot/max(1,len(lines)):.2f}  core-us per candidate bit {1e6*tot/max(1,bits):.0f}')
print('top threads (cores):',[(t,round(v/dt,3)) for v,t in others[:12]])
low=[t for t in tids if t<min(tids)+600]

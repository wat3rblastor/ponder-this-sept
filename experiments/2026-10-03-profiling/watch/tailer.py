import os,time,sys
D='/workspace/ponder-this-sept/experiments/remote'; dur=float(sys.argv[1])
fs={i:open(f'{D}/v3_s{i}.log') for i in range(16)}
for f in fs.values(): f.seek(0,2)
buf={i:'' for i in range(16)}
out=open('units_ts.txt','w'); t0=time.time()
while time.time()-t0<dur:
    now=time.time()
    for i,f in fs.items():
        buf[i]+=f.read()
        *lines,buf[i]=buf[i].split('\n')
        for l in lines:
            if l.startswith('K=') or l.startswith('***'): out.write(f'{now:.2f} {i} {l}\n')
    out.flush(); time.sleep(1.0)

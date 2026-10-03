import time,glob,json
fs=sorted(glob.glob('/workspace/ponder-this-sept/experiments/remote/r1_g*.jsonl'))
pos={f:0 for f in fs}
out=open('/workspace/mathres/samp.jsonl','w')
def step():
    row={'t':time.time()}
    for f in fs:
        fh=open(f); fh.seek(pos[f]); data=fh.read(); 
        k=data.rfind('\n')+1; data=data[:k]; pos[f]+=k
        small=0;sres=0;big=0;bres=0;mid=0;mres=0
        for line in data.splitlines():
            if '"covered"' not in line: continue
            j=json.loads(line)
            if j['res']<5e9: small+=1;sres+=j['res']
            elif j['res']<1e11: mid+=1;mres+=j['res']
            else: big+=1;bres+=j['res']
        row[f[-8:-6]]=[small,sres,mid,mres,big,bres]
    out.write(json.dumps(row)+'\n'); out.flush()
step()
for i in range(40):
    time.sleep(15); step()

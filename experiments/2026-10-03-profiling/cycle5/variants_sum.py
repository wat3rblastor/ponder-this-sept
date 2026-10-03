import sys,re,collections
cfg=None; T=collections.defaultdict(lambda: collections.defaultdict(list)); S=collections.defaultdict(dict)
for l in open(sys.argv[1]):
    if l.startswith('==='): cfg=l[4:].strip(); continue
    m=re.match(r'ksolo K=(\d+) sh=(\d+) res=([\d.e+]+) ntc=\d+ ([\d.]+)s',l)
    if m: T[cfg][(m[1],m[2])].append(float(m[4])); continue
    m=re.match(r'K=(\d+) sh=(\d+) .*surv=(\d+) bits=(\d+) conf=(\d+)',l)
    if m: S[cfg][(m[1],m[2])]=(m[3],m[4],m[5])
base='nch=7 t0=24'
for cfg in T:
    same=all(S[cfg].get(k)==S[base].get(k) for k in S[base])
    print(cfg.ljust(14),' '.join(f'{k[0]}:{min(v)/min(T[base][k]):.3f}' for k,v in T[cfg].items()),' results identical' if same else ' RESULTS DIFFER')

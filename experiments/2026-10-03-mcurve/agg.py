#!/usr/bin/env python3
"""Aggregate mcurve raw JSONL into M(n) = min last term, per sweep."""
import json, sys, glob, collections

def agg(patterns):
    best = {}
    Nmax = None; mmax = 0; base = None
    for pat in patterns:
        for fn in sorted(glob.glob(pat)):
            for line in open(fn):
                line = line.strip()
                if not line: continue
                try: r = json.loads(line)
                except Exception: continue
                if 'm' not in r: continue
                Nmax = r['N']; mmax = max(mmax, r['m'])
                b = r['d'] // r['m']
                base = b if base is None else base
                for k, a in r.get('mina', {}).items():
                    n = int(k); a = int(a)
                    if a == 0: continue          # reject degenerate a=0
                    last = a + (n - 1) * r['d']
                    if last > Nmax: continue     # truncated-run artefact guard
                    if n not in best or last < best[n][0]:
                        best[n] = (last, a, r['d'])
    return best, Nmax, mmax, base

if __name__ == '__main__':
    best, N, mmax, base = agg(sys.argv[1:])
    print(f"# sweep: N={N} base={base} mmax={mmax}")
    for n in sorted(best):
        last, a, d = best[n]
        # completeness: need N >= last AND mmax*base*(n-1) >= last
        cov_m = mmax * base * (n - 1)
        proved = (N >= last) and (cov_m >= last)
        print(f"n={n:3d} M={last:<14d} a={a:<14d} d={d:<10d} "
              f"{'PROVED' if proved else 'BOUND '} (need m<= {last//(base*(n-1))+1}, have {mmax})")

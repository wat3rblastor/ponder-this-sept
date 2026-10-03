#!/usr/bin/env python3
"""For every configuration found by the discriminant-level search, decide
whether an INTEGRAL pencil (a,d) actually realises it, and report the
realizable maximum per sign of alpha.

Anchoring: pencils() needs a positive definite member at k=0, so the target is
re-anchored at the smallest hit k1: S'(j) = S(j+k1).  A pencil (a',d) for S'
gives the original pencil a = a' - k1*d, and true_hits() then counts the
genuinely automatic indices in the ORIGINAL window [0,57].
"""
from __future__ import annotations
import re, sys, math, collections
from math import isqrt
from pencil import pencils, true_hits, disc_hits, issq

def hits_of(alpha, beta, gamma):
    return [k for k in range(58)
            if (v := alpha*k*k+beta*k+gamma) >= 0 and issq(v)]

def realize(alpha, beta, gamma, hits):
    k1 = hits[0]
    g2 = alpha*k1*k1 + beta*k1 + gamma
    b2 = 2*alpha*k1 + beta
    best = 0; bestpen = None
    for (a2, d) in pencils(alpha, b2, g2):
        a = tuple(a2[i] - k1*d[i] for i in range(3))
        th = true_hits(a, d)
        if len(th) > best:
            best = len(th); bestpen = (a, d, sorted(th))
    return best, bestpen

if __name__ == "__main__":
    files = sys.argv[1:] or ['../2026-10-03-covering/exh_pos.txt',
                             '../2026-10-03-covering/exh_neg.txt',
                             '../2026-10-03-covering/exh_lin2.txt']
    cfg = {}
    for fn in files:
        try: fh = open(fn)
        except OSError: continue
        for L in fh:
            m = re.match(r'hits=(\d+) alpha=(-?\d+) beta=(-?\d+) gamma=(-?\d+)', L)
            if not m: continue
            h, al, be, ga = map(int, m.groups())
            if h < 12: continue
            cfg[(al, be, ga)] = h
    print(f"{len(cfg)} configurations with >=12 discriminant-level hits")
    stat = collections.defaultdict(lambda: [0, 0, 0])   # sign -> [n, maxdisc, maxreal]
    exemplar = {}
    for (al, be, ga), h in sorted(cfg.items(), key=lambda t: -t[1]):
        sg = 'alpha<0' if al < 0 else ('alpha>0' if al > 0 else 'alpha=0')
        D = be*be - 4*al*ga
        s = stat[sg]; s[0] += 1
        s[1] = max(s[1], h)
        if D <= 0:                      # Theorem 1: not realizable by any pencil
            continue
        hits = hits_of(al, be, ga)
        if len(hits) != h: h = len(hits)
        if h <= s[2]: continue
        r, pen = realize(al, be, ga, hits)
        if r > s[2]:
            s[2] = r; exemplar[sg] = (r, al, be, ga, pen)
        sys.stdout.flush()
    print()
    for sg in ('alpha=0', 'alpha>0', 'alpha<0'):
        if sg not in stat: continue
        n, md, mr = stat[sg]
        print(f"{sg}: {n:6d} configs, max disc-level hits {md}, "
              f"max REALIZED by an integral pencil {mr}")
        if sg in exemplar:
            r, al, be, ga, pen = exemplar[sg]
            a, d, th = pen
            print(f"   best: S={al}k^2{be:+d}k{ga:+d}  a={a} d={d}  |A|={r}  A={th}")

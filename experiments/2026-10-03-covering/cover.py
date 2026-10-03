#!/usr/bin/env python3
"""Step 1: the covering problem as literally posed.

Automatic indices from one condition (c, k0) are  {k0 + c*j^2} cap [0,N-1].
Compute exact minimum number of conditions to cover 20,30,40,50,58 indices
(and the max coverage with r conditions), by ILP-free exact branch and bound
on a set-cover / max-coverage instance.  Sets are tiny (<=8 elements) so this
is easy.
"""
import sys
from functools import lru_cache
from itertools import combinations

N = 58


def sets_all():
    out = {}
    for c in range(1, N + 1):
        for k0 in range(0, N):
            s = frozenset(k0 + c * j * j for j in range(0, N) if k0 + c * j * j < N)
            if len(s) >= 2:
                out.setdefault(s, (c, k0))
    # keep only maximal sets (a subset never helps)
    keys = list(out)
    keys.sort(key=lambda s: -len(s))
    maximal = []
    for s in keys:
        if not any(s < t for t in maximal):
            maximal.append(s)
    return [(s, out[s]) for s in maximal]


def max_coverage(sets, r):
    """exact max |union| over r sets, branch and bound"""
    masks = []
    for s, lbl in sets:
        m = 0
        for e in s:
            m |= 1 << e
        masks.append((m, bin(m).count('1'), lbl))
    masks.sort(key=lambda t: -t[1])
    best = [0, None]

    def rec(i, used, cov, chosen):
        pc = bin(cov).count('1')
        if pc > best[0]:
            best[0] = pc
            best[1] = list(chosen)
        if used == r:
            return
        # bound: remaining picks * largest set size
        if pc + (r - used) * masks[0][1] <= best[0]:
            return
        for j in range(i, len(masks)):
            m, sz, lbl = masks[j]
            if pc + (r - used) * sz <= best[0]:
                break
            rec(j + 1, used + 1, cov | m, chosen + [lbl])

    rec(0, 0, 0, [])
    return best


if __name__ == '__main__':
    S = sets_all()
    print(f"distinct maximal sets (c,k0): {len(S)}")
    sz = {}
    for s, lbl in S:
        sz[len(s)] = sz.get(len(s), 0) + 1
    print("set sizes:", dict(sorted(sz.items(), reverse=True)))
    big = [(sorted(s), lbl) for s, lbl in S if len(s) >= 7]
    for s, lbl in big:
        print("  size", len(s), "(c,k0)=", lbl, s)
    print()
    for r in range(1, 15):
        cov, chosen = max_coverage(S, r)
        print(f"r={r:2d} conditions -> max coverage {cov:2d} / {N}   e.g. {chosen}")
        if cov == N:
            break

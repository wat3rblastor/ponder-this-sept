#!/usr/bin/env python3
"""Measure F17: automatic-index failures (must be zero) and the pass rate on the
41 non-automatic indices.  Uses src/loeschian.is_loeschian (exact) only."""
from __future__ import annotations
import sys, math, random, time
from pathlib import Path
from multiprocessing import Pool
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from loeschian import is_loeschian, factorize  # noqa
from family import A_FORM, D_FORM, AUTO, M, val, member  # noqa

NONAUTO = [k for k in range(58) if k not in AUTO]
BAD = [p for p in range(2, 5000) if p % 3 == 2 and all(p % q for q in range(2, int(p**.5)+1))]


def fastloe(n):
    for q in BAD:
        if q*q > n: break
        if n % q == 0:
            e = 0
            while n % q == 0:
                n //= q; e += 1
            if e & 1: return False
    return is_loeschian(n)


def scan(args):
    M1, P, R = args
    a, d, p, q = member(M1, P, R)
    if a <= 0 or d <= 0: return None
    autofail = [k for k in AUTO if not fastloe(a + k*d)]
    na = [k for k in NONAUTO if fastloe(a + k*d)]
    fl = [True]*58
    for k in NONAUTO: fl[k] = k in na
    best = cur = 0
    for f in fl:
        cur = cur+1 if f else 0
        best = max(best, cur)
    return (a, d, a+57*d, autofail, len(na), best)


if __name__ == "__main__":
    M1 = int(sys.argv[1]) if len(sys.argv) > 1 else 493
    NP = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    NR = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    jobs = [(M1, P, R) for P in range(1, NP+1) for R in range(1, NR+1)]
    t0 = time.time()
    with Pool(8) as pool:
        res = [r for r in pool.map(scan, jobs, chunksize=4) if r]
    nterm = len(res)*41
    tot = sum(r[4] for r in res)
    af = sum(len(r[3]) for r in res)
    print(f"M1={M1}  members tested {len(res)}  time {time.time()-t0:.1f}s")
    print(f"automatic-index failures: {af} out of {len(res)*17} automatic terms")
    print(f"non-automatic pass rate: {tot}/{nterm} = {tot/nterm:.4f}")
    import collections
    print("best-run hist:", sorted(collections.Counter(r[5] for r in res).items()))
    print("non-auto count hist:", sorted(collections.Counter(r[4] for r in res).items()))
    res.sort(key=lambda r: -r[5])
    for r in res[:5]:
        print(f"  best: a={r[0]} d={r[1]} T={r[2]} ({len(str(r[2]))} digits) run={r[5]} nonauto={r[4]}")
    Ts = sorted(r[2] for r in res)
    print(f"  T range {Ts[0]} .. {Ts[-1]}")

#!/usr/bin/env python3
"""Step 2c: the best *realizable* configuration found, made explicit.

Pentagonal family, 13 automatic indices:
    d = 24 m^2,  a = 3 x^2 + m^2,  t_k = 3x^2 + m^2 (24k+1)
automatic at A = {k : 24k+1 is a square} = 13 indices.

badprimes.py shows d only needs the primes 2,5,11,17 (not 23,29,41,47,53,
whose forced hit classes can be parked on automatic indices).  So
m = 935 * w  with 935 = 5*11*17, instead of m = 623645 * w.

This script: picks the congruence classes, finds concrete members, verifies
the 13 automatic terms really are Loeschian, measures the pass rate of the 45
non-automatic terms, and counts the family members below a size bound.
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from loeschian import is_loeschian  # noqa: E402

N = 58
A = [k for k in range(N) if int((24 * k + 1) ** 0.5 + 0.5) ** 2 == 24 * k + 1]
NONAUTO = [k for k in range(N) if k not in A]
FORCED = [2, 5, 11, 17]          # must divide d
PARKED = [23, 29, 41, 47, 53]    # may miss d, hit class parked on A


def admissible_classes(p):
    return [r for r in range(p)
            if set(k for k in range(N) if k % p == r) <= set(A)
            and any(k % p == r for k in range(N))]


def setup(m):
    """choose, for each parked prime, a class r and the x-classes mod p"""
    d = 24 * m * m
    cond = {}
    for p in PARKED:
        assert d % p and m % p
        found = None
        for r in admissible_classes(p):
            # need 3x^2 + m^2 = -r*d  (mod p)
            rhs = (-r * d - m * m) % p
            inv3 = pow(3, -1, p)
            tgt = rhs * inv3 % p
            xs = [x for x in range(p) if x * x % p == tgt]
            if xs:
                found = (r, xs)
                break
        if found is None:
            return None
        cond[p] = found
    return d, cond


def members(m, cond, count, start=1):
    """x values meeting every congruence; yields (a, d)"""
    d = 24 * m * m
    x = start
    got = 0
    while got < count:
        x += 1
        if x % 2:
            continue
        if any(x % p == 0 for p in (5, 11, 17)):
            continue
        if any(x % p not in cond[p][1] for p in PARKED):
            continue
        yield 3 * x * x + m * m, d
        got += 1


if __name__ == '__main__':
    m = 935
    d, cond = setup(m)
    print(f"A ({len(A)} automatic indices) = {A}")
    print(f"m = {m} = 5*11*17, d = 24 m^2 = {d}")
    for p in PARKED:
        print(f"  prime {p}: hit class r={cond[p][0]} (indices "
              f"{[k for k in range(N) if k % p == cond[p][0]]}), x = {cond[p][1]} mod {p}")
    print(f"T_min = a+57d at x=0 would be {1369*m*m} ({len(str(1369*m*m))} digits)")
    print()

    # --- verify the automatic indices on real members, measure the rest
    nmem = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    autopass = 0
    autototal = 0
    hits = [0] * N
    tot = 0
    best = (0, None)
    for a, d_ in members(m, cond, nmem):
        tot += 1
        ok = 0
        for k in range(N):
            t = a + k * d_
            L = is_loeschian(t)
            if L:
                hits[k] += 1
            if k in A:
                autototal += 1
                autopass += L
            else:
                ok += L
        if ok > best[0]:
            best = (ok, (a, d_))
    print(f"members tested: {tot}  (smallest T = {3*2*2+m*m+57*d}...)")
    print(f"automatic indices: {autopass}/{autototal} Loeschian "
          f"({'ALL PASS' if autopass == autototal else 'FAILURE'})")
    nap = sum(hits[k] for k in NONAUTO) / (tot * len(NONAUTO))
    print(f"non-automatic pass rate (45 indices): {nap:.4f}")
    print(f"best member: {best[0]}/45 non-automatic terms Loeschian, (a,d)={best[1]}")
    print(f"per-member probability estimate: {nap}^45 = {nap**45:.3e}")

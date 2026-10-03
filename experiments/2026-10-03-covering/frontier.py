#!/usr/bin/env python3
"""Step 1+2: the frontier.  coverage  vs  minimum term size T.

For a quadratic P(k) = alpha k^2 + beta k + gamma the automatic set is
    A(P) = {k in [0,57] : P(k) = y_k^2 is a perfect square},
and any family realizing P satisfies
    (i)  T >= 2 * max_{k in A} y_k          (minimum of a pos. def. binary
         quadratic form of discriminant -12 y^2 is >= 2y)
    (ii) T >= 57 * d  and every bad prime p < 58 that is NOT parkable
         divides d, where p is parkable only if some class r mod p has all
         its indices in A *and* p | y_k for each of them (a bad prime divides
         X^2+3Y^2 only when it divides both X and Y, i.e. only when it
         divides the determinant y_k of the representing matrix).

This script enumerates quadratics directly (same loop as qsearch.c, in
Python, small boxes) or reads a qsearch/sym log, and prints the frontier.
"""
import re
import subprocess
import sys
from math import isqrt

N = 58
BADP = [2, 5, 11, 17, 23, 29, 41, 47, 53]


def hits(al, be, ga):
    out = {}
    for k in range(N):
        v = al * k * k + be * k + ga
        if v >= 0:
            y = isqrt(v)
            if y * y == v:
                out[k] = y
    return out


def parkable(p, H):
    """can bad prime p be left out of d, given automatic set H={k:y_k}?"""
    for r in range(p):
        cls = [k for k in range(N) if k % p == r]
        if not cls:
            continue
        if all(k in H and H[k] % p == 0 for k in cls):
            return r
    return None


def forced_product(H):
    prod = 1
    forced = []
    for p in BADP:
        if parkable(p, H) is None:
            forced.append(p)
            prod *= p
    return prod, forced


def report(cfgs):
    """cfgs: list of (al,be,ga); print frontier coverage -> min T bound"""
    best = {}
    for al, be, ga in cfgs:
        H = hits(al, be, ga)
        c = len(H)
        if c < 2:
            continue
        T1 = 2 * max(H.values())
        prod, forced = forced_product(H)
        T2 = 57 * prod
        T = max(T1, T2)
        if c not in best or T < best[c][0]:
            best[c] = (T, al, be, ga, T1, T2, forced, sorted(H))
    print(f"{'cov':>4} {'T_lower':>14} {'2*maxy':>9} {'57*prod':>10} "
          f"{'alpha':>7} {'beta':>8} {'gamma':>10}  forced primes")
    for c in sorted(best):
        T, al, be, ga, T1, T2, forced, A = best[c]
        print(f"{c:>4} {T:>14} {T1:>9} {T2:>10} {al:>7} {be:>8} {ga:>10}  {forced}")
    return best


if __name__ == '__main__':
    cfgs = []
    for fn in sys.argv[1:]:
        for line in open(fn):
            m = re.match(r'hits=(\d+) alpha=(-?\d+) beta=(-?\d+) gamma=(-?\d+)', line)
            if m:
                cfgs.append(tuple(int(x) for x in m.groups()[1:]))
            m = re.match(r'hits=(\d+) alpha=(-?\d+) h=(-?\d+) C=(-?\d+)', line)
            if m:
                _, al, h, C = (int(x) for x in m.groups())
                # P(k) = al*(2k-h)^2 + C
                cfgs.append((4 * al, -4 * al * h, al * h * h + C))
    print(f"configurations read: {len(cfgs)}")
    report(cfgs)

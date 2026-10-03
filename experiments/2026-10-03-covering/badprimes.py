#!/usr/bin/env python3
"""Step 2b: which bad primes must divide d, given an automatic set A?

A bad prime p (p = 2 mod 3) with p not dividing d hits the indices
{k in [0,57] : k = r (mod p)} for exactly one r, and at those indices the term
needs v_p even.  If every hit index is AUTOMATIC the term is Loeschian no
matter what divides it, so p may be left out of d at zero cost in size -- it
only costs one congruence a = -r*d (mod p).

For each bad prime p < 58 this prints the admissible r (hit class fully inside
A), hence whether p can be dropped from d.
"""
import sys

N = 58
def isprime(n):
    return n > 1 and all(n % i for i in range(2, int(n ** 0.5) + 1))


BAD = [p for p in range(2, N) if p % 3 == 2 and isprime(p)]


def squares_ap(b, g):
    out = []
    for k in range(N):
        v = b * k + g
        if v >= 0 and int(v ** 0.5 + 0.5) ** 2 == v:
            out.append(k)
    return out


def analyse(A, label):
    A = set(A)
    print(f"--- {label}: |A|={len(A)}  A={sorted(A)}")
    forced, droppable = [], []
    for p in BAD:
        ok = []
        for r in range(p):
            hits = [k for k in range(N) if k % p == r]
            if hits and set(hits) <= A:
                ok.append((r, hits))
        if ok:
            droppable.append(p)
            print(f"  p={p:2d}: DROPPABLE, classes {[(r, h) for r, h in ok]}")
        else:
            forced.append(p)
            mn = min(len([k for k in range(N) if k % p == r]) for r in range(p))
            print(f"  p={p:2d}: forced into d (min #hits in a class = {mn})")
    print(f"  forced primes: {forced}")
    print(f"  droppable:     {droppable}")
    return forced, droppable


if __name__ == '__main__':
    pent = squares_ap(24, 1)          # pentagonal family, 13 automatic
    sq = [j * j for j in range(8)]     # a=x^2, d=3m^2, 8 automatic
    analyse(sq, "c=1 family (a=x^2, d=3m^2), A={j^2}")
    print()
    analyse(pent, "pentagonal family (a=3x^2+m^2, d=24m^2), A={k:24k+1 square}")

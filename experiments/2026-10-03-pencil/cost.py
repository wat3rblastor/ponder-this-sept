#!/usr/bin/env python3
"""Cost model for F17, all constants measured.

Member count.  A member is a pair (p,q) with
    d(p,q) = 3 p (6p+q) > 0,  M = 1247290 | d,  T = a+57d <= Tmax.
T(p,q) = a + 57 d = 1029 p^2 + 168 p q + 7 q^2, positive definite, disc = -588.
(p,q) and (-p,-q) give the same (a,d), so members are counted up to sign.

The density of (p,q) mod M with M | p(6p+q) is prod_{l|M} (2l-1)/l^2.
"""
from __future__ import annotations
import sys, math, random
from math import isqrt, gcd
from pathlib import Path
from multiprocessing import Pool
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from loeschian import is_loeschian  # noqa
from family import A_FORM, D_FORM, AUTO, M, val  # noqa

NONAUTO = [k for k in range(58) if k not in AUTO]
PRIMES_M = [2, 5, 11, 17, 23, 29]
BAD = [p for p in range(2, 10000) if p % 3 == 2 and all(p % q for q in range(2, int(p**.5)+1))]


def T_of(p, q):
    return 1029*p*p + 168*p*q + 7*q*q


def exact_member_count(Tmax):
    """brute-force exact count of members with T <= Tmax (small Tmax only)."""
    n = 0
    # T >= (1029*7-84^2/...)  bound p: min over q of T at fixed p is
    # T = 7(q + 12p)^2 + (1029-1008)p^2 = 7(q+12p)^2 + 21 p^2
    pmax = isqrt(Tmax//21)
    for p in range(1, pmax+1):
        rem = Tmax - 21*p*p
        if rem < 0: continue
        s = isqrt(rem//7)
        for u in range(-s, s+1):          # u = q + 12 p
            q = u - 12*p
            d = 3*p*(6*p+q)
            if d <= 0 or d % M: continue
            if 7*u*u + 21*p*p > Tmax: continue
            n += 1
    # p<0 is the (-p,-q) image of p>0, already excluded by counting up to sign
    return n


def asymptotic_density():
    dens = 1.0
    for l in PRIMES_M:
        dens *= (2*l - 1) / l**2
    area = 2*math.pi/math.sqrt(588)          # area of {T<=1} = 2 pi / sqrt|disc|
    return area * dens / 2                   # /2 for the +-(p,q) identification


def fastloe(n):
    for q in BAD:
        if q*q > n: break
        if n % q == 0:
            e = 0
            while n % q == 0:
                n //= q; e += 1
            if e & 1: return False
    return is_loeschian(n)


def sample_member(rng, Tmax, M1):
    """random member with T ~ Tmax, via the M1 | p, (M/M1) | 6p+q sublattice"""
    M2 = M // M1
    for _ in range(400):
        # T = 273 M1^2 P^2 + 84 M1 M2 P R + 7 M2^2 R^2
        Pb = isqrt(Tmax // (273*M1*M1)) or 1
        Rb = isqrt(Tmax // (7*M2*M2)) or 1
        P = rng.randrange(1, Pb+1); R = rng.randrange(1, Rb+1)
        p = M1*P; q = M2*R - 6*M1*P
        T = T_of(p, q)
        if T <= Tmax and T >= Tmax//10:
            return p, q
    return None


def work(args):
    seed, Tmax, M1, nmem = args
    rng = random.Random(seed)
    npass = ntot = autofail = 0
    runs = []
    for _ in range(nmem):
        m = sample_member(rng, Tmax, M1)
        if not m: continue
        p, q = m
        a = val(A_FORM, p, q); d = val(D_FORM, p, q)
        assert d % M == 0 and d > 0 and a > 0
        for k in AUTO:
            if not fastloe(a+k*d): autofail += 1
        fl = [True]*58
        for k in NONAUTO:
            ok = fastloe(a+k*d)
            fl[k] = ok
            npass += ok; ntot += 1
        best = cur = 0
        for f in fl:
            cur = cur+1 if f else 0
            best = max(best, cur)
        runs.append(best)
    return npass, ntot, autofail, runs


if __name__ == "__main__":
    print("asymptotic members-below-T constant c =", f"{asymptotic_density():.4e}",
          " (members(T) ~ c*T)")
    for e in (9, 10, 11, 12):
        Tm = 10**e
        n = exact_member_count(Tm)
        print(f"  exact members(T<=1e{e}) = {n:12d}   c_eff = {n/Tm:.3e}")
    print()
    mode = sys.argv[1] if len(sys.argv) > 1 else "rate"
    if mode == "rate":
        for e, M1, nm in [(12, 493, 120), (15, 493, 120), (18, 493, 90),
                          (21, 493, 60), (24, 493, 40)]:
            jobs = [(1000*e+i, 10**e, M1, nm) for i in range(8)]
            with Pool(8) as pool:
                res = pool.map(jobs.__getitem__ if False else work, jobs)
            npass = sum(r[0] for r in res); ntot = sum(r[1] for r in res)
            af = sum(r[2] for r in res); runs = sum((r[3] for r in res), [])
            rho = npass/ntot
            import math as _m
            print(f"T~1e{e}: members {len(runs):5d}  autofail {af}  "
                  f"rho_nonauto = {rho:.4f}  rho^41 = {rho**41:.3e}  "
                  f"maxrun {max(runs)}  meancount {npass/len(runs):.2f}")
            sys.stdout.flush()

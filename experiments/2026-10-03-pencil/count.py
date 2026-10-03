#!/usr/bin/env python3
"""EXACT member count for F17 below a term-size bound T.

A member is (p,q) with d = 3p(6p+q) > 0 and M = 1247290 | d and
T(p,q) = a + 57 d = 1029 p^2 + 168 p q + 7 q^2 <= Tmax.
Substituting u = q + 12 p gives T = 7 u^2 + 21 p^2 and 6p+q = u - 6p, so

    p >= 1,  u > 6p,  7u^2 + 21p^2 <= Tmax,  M | p (u - 6p).

Since M is squarefree, M | p(u-6p) iff M' | (u-6p) where M' is the product of
the primes of M that do not divide p.  So for each p the admissible u form a
single arithmetic progression and the count is exact in O(sqrt Tmax).
(p,q) and (-p,-q) give the same (a,d); the constraint 6p+q>0 with p>0 picks
exactly one of the pair.
"""
import sys, math
from math import isqrt, gcd
PRIMES_M = [2, 5, 11, 17, 23, 29]
M = 1
for l in PRIMES_M: M *= l

def Mprime(p):
    Mp = 1
    for l in PRIMES_M:
        if p % l: Mp *= l
    return Mp

def members_pos(Tmax):
    """d > 0 branch: largest term is a+57d = 7u^2+21p^2, u = q+12p > 6p."""
    tot = 0
    for p in range(1, isqrt(Tmax // 21) + 1):
        rem = Tmax - 21*p*p
        if rem < 0: break
        Umax = isqrt(rem // 7)
        if Umax <= 6*p: continue
        tot += (Umax - 6*p) // Mprime(p)
    return tot

def members_neg(Tmax):
    """d < 0 branch (AP read backwards; automatic set is 57-A, also 17 indices):
    w = -q > 6p, largest term is a = 3p^2+3pw+7w^2."""
    tot = 0
    for p in range(1, isqrt(Tmax // 273) + 1):
        # 7w^2+3pw+3p^2 <= Tmax
        disc = 9*p*p - 28*(3*p*p - Tmax)
        Wmax = (-3*p + isqrt(disc)) // 14
        while 7*Wmax*Wmax + 3*p*Wmax + 3*p*p > Tmax: Wmax -= 1
        while 7*(Wmax+1)**2 + 3*p*(Wmax+1) + 3*p*p <= Tmax: Wmax += 1
        if Wmax <= 6*p: continue
        tot += (Wmax - 6*p) // Mprime(p)
    return tot

def members(Tmax):
    return members_pos(Tmax) + members_neg(Tmax)

def asym():
    dens = 1.0
    for l in PRIMES_M: dens *= (2*l-1)/l**2
    return 2*math.pi/math.sqrt(588) * dens / 2

if __name__ == "__main__":
    print(f"asymptotic constant c = {asym():.5e}   (members(T) ~ c T)")
    for e in range(8, 23):
        T = 10**e
        if e <= 16:
            n = members(T)
            print(f"  members(T<=1e{e:2d}) = {n:>18,d}    n/T = {n/T:.4e}")
        else:
            print(f"  members(T<=1e{e:2d}) ~ {asym()*T:>18.4e}  (asymptotic)")

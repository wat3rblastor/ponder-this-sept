#!/usr/bin/env python3
"""The constructed 17-automatic family F17, and its verification.

PENCIL (constructed by pencil.py, target S(k) = -3k^2 + 174k + 25):

    a(p,q) = 3 p^2 - 3 p q + 7 q^2          (positive definite, disc -75)
    d(p,q) = 18 p^2 + 3 p q = 3 p (6 p + q) (INDEFINITE, disc +9 = a square!)

    t_k = a + k d,  automatic at the 17 indices
    A = {0,1,2,6,8,12,17,22,25,33,36,41,46,50,52,56,57}

The point of disc(d) = 9 being a perfect square: d SPLITS into linear forms, so
forcing the mandatory bad primes 2,5,11,17,23,29 into d costs only ONE linear
congruence per prime (index M in the (p,q) lattice) instead of the square
condition 623645 | m of the classical pentagonal family d = 24 m^2.  Hence
d_min = 3M = 3741870 instead of 24*623645^2 = 9.33e12, and the family is
2-dimensional.

PARAMETRISATION of the sub-family with M | d:
    M = 2*5*11*17*23*29 = 1247290 = M1 * M2,
    p = M1 * P,  6p + q = M2 * R   (so q = M2 R - 6 M1 P),  P,R free integers.
    d = 3 M1 M2 P R = 3 M P R
    T(P,R) = a + 57 d = 273 M1^2 P^2 + 84 M1 M2 P R + 7 M2^2 R^2
             (positive definite, disc = -588 M^2)
"""
from __future__ import annotations
import sys, math
from math import isqrt, gcd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from loeschian import is_loeschian, factorize   # noqa

A_FORM = (3, -3, 7)
D_FORM = (18, 3, 0)
AUTO = [0,1,2,6,8,12,17,22,25,33,36,41,46,50,52,56,57]
M = 2*5*11*17*23*29      # 1247290


def val(F, p, q):
    return F[0]*p*p + F[1]*p*q + F[2]*q*q


def reps_Q(n):
    out = []
    if n == 0: return [(0,0)]
    for v in range(-isqrt(4*n//3)-1, isqrt(4*n//3)+2):
        r = 4*n - 3*v*v
        if r < 0: continue
        s = isqrt(r)
        if s*s != r: continue
        for sg in (s,-s):
            if (sg - v) % 2: continue
            u = (sg - v)//2
            if u*u+u*v+v*v == n: out.append((u,v))
    return sorted(set(out))


def as_QoM(A,B,C):
    if A < 0 or C < 0: return None
    for (u1,v1) in reps_Q(A):
        for (u2,v2) in reps_Q(C):
            if 2*u1*u2 + u1*v2 + u2*v1 + 2*v1*v2 == B:
                return (u1,u2,v1,v2)
    return None


def pencil_certificates():
    """L1,L2 for every automatic index, checked as polynomial identities."""
    cert = {}
    for k in range(58):
        F = tuple(A_FORM[i] + k*D_FORM[i] for i in range(3))
        m = as_QoM(*F)
        if m is None: continue
        u1,u2,v1,v2 = m
        # identity check: (u1 p+u2 q)^2 + (..)(..) + (v1 p+v2 q)^2 == F(p,q)
        assert u1*u1+u1*v1+v1*v1 == F[0]
        assert u2*u2+u2*v2+v2*v2 == F[2]
        assert 2*u1*u2+u1*v2+u2*v1+2*v1*v2 == F[1]
        cert[k] = m
    return cert


def member(M1, P, R):
    M2 = M // M1
    p = M1*P
    q = M2*R - 6*M1*P
    a = val(A_FORM, p, q)
    d = val(D_FORM, p, q)
    return a, d, p, q


def enumerate_members(M1, Tmax, limit=None):
    """all (P,R) with P,R>=1 (so d>0) and T = a+57d <= Tmax, deduplicated."""
    M2 = M // M1
    out = []
    P = 1
    while 273*M1*M1*P*P + 84*M1*M2*P + 7*M2*M2 <= Tmax:
        R = 1
        while True:
            T = 273*M1*M1*P*P + 84*M1*M2*P*R + 7*M2*M2*R*R
            if T > Tmax: break
            out.append((P,R,T))
            R += 1
        P += 1
        if limit and len(out) > limit: break
    return out


if __name__ == "__main__":
    cert = pencil_certificates()
    print("automatic indices (polynomial identity verified):", sorted(cert))
    assert sorted(cert) == AUTO, sorted(cert)
    for k in AUTO:
        u1,u2,v1,v2 = cert[k]
        print(f"  k={k:3d}  t_k = L1^2+L1L2+L2^2  L1={u1}p{u2:+d}q  L2={v1}p{v2:+d}q"
              f"   det={u1*v2-u2*v1}")

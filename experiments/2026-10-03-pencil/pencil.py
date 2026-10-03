#!/usr/bin/env python3
"""Construct explicit INTEGRAL pencils of binary quadratic forms (a,d) whose
automatic set is prescribed.

Conventions.  A form F = (A,B,C) means A p^2 + B p q + C q^2.  Q = (1,1,1) is
x^2+xy+y^2, so F is automatically Loeschian for every (p,q) iff

    F = Q o M  for an integer 2x2 matrix M = [[u1,u2],[v1,v2]],
    i.e.  F(p,q) = L1^2 + L1 L2 + L2^2  with L1 = u1 p + u2 q, L2 = v1 p + v2 q.

Then disc(F) = -3 det(M)^2, so with

    S(k) := -disc(a + k d)/3 = (4 (A+kA')(C+kC') - (B+kB')^2) / 3

a NECESSARY condition for k automatic is that S(k) be a perfect square.  That is
the discriminant-level screen used by experiments/2026-10-03-covering.  It is
NOT sufficient: F must actually lie in the Q o M orbit-image, which is a real
extra condition (checked exactly here by `as_QoM`).

NOTE on normalisation: the covering experiment wrote P(k) = -disc/12, which is
S(k)/4 and silently restricts to det(M) even.  Its (alpha,beta,gamma) loops run
over all integers so the *search* was complete, but the translation to (a,d)
differs: here disc(d) = -3*alpha, disc(a) = -3*gamma.
"""
from __future__ import annotations
import sys, math
from math import isqrt, gcd

def issq(n):
    if n < 0: return False
    r = isqrt(n); return r*r == n

def reps_Q(n):
    """all (u,v) with u^2+uv+v^2 == n, n>=0"""
    out = []
    if n == 0: return [(0,0)]
    # u^2+uv+v^2 = n  =>  (2u+v)^2 + 3v^2 = 4n
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
    """Return M=[[u1,u2],[v1,v2]] with Q o M == (A,B,C), or None."""
    if A < 0 or C < 0: return None
    for (u1,v1) in reps_Q(A):
        for (u2,v2) in reps_Q(C):
            if 2*u1*u2 + u1*v2 + u2*v1 + 2*v1*v2 == B:
                return (u1,u2,v1,v2)
    return None

def Sval(a,d,k):
    A,B,C = a; A2,B2,C2 = d
    num = 4*(A+k*A2)*(C+k*C2) - (B+k*B2)**2
    assert num % 3 == 0 or True
    return num/3 if num % 3 else num//3

def disc_hits(a,d,n=58):
    """indices where S(k) is a non-negative perfect square (necessary cond.)"""
    out = []
    for k in range(n):
        A,B,C = a; A2,B2,C2 = d
        num = 4*(A+k*A2)*(C+k*C2) - (B+k*B2)**2
        if num < 0 or num % 3: continue
        if issq(num//3): out.append(k)
    return out

def true_hits(a,d,n=58):
    """indices where a+kd is genuinely Q o M (sufficient for automatic)."""
    out = {}
    A,B,C = a; A2,B2,C2 = d
    for k in range(n):
        F = (A+k*A2, B+k*B2, C+k*C2)
        num = 4*F[0]*F[2]-F[1]**2
        if num < 0: continue
        m = as_QoM(*F)
        if m is not None: out[k] = m
    return out

def reduced_forms(negdisc3):
    """all reduced positive definite (A,B,C) with 4AC-B^2 == 3*gamma"""
    tgt = 3*negdisc3
    out = []
    Amax = isqrt(tgt//3)+1
    for A in range(1, Amax+1):
        for B in range(-A, A+1):
            num = B*B + tgt
            if num % (4*A): continue
            C = num//(4*A)
            if C < A: continue
            out.append((A,B,C))
    return out

def pencils(alpha, beta, gamma, verbose=False):
    """all integral pencils (a,d) with a reduced pos.def., realising
       S(k) = alpha k^2 + beta k + gamma  (requires gamma a square, 0 a hit)."""
    res = []
    for (A,B,C) in reduced_forms(gamma):
        # A B'^2 - 2 B A' B' + (4C A'^2 - 3 beta A' + 3 alpha A) = 0
        # disc_B' = 12 [ -gamma A'^2 + A beta A' - alpha A^2 ]
        # solve the A'-range
        # gamma A'^2 - A beta A' + alpha A^2 <= 0
        if gamma == 0: continue
        disc = (A*beta)**2 - 4*gamma*alpha*A*A
        if disc < 0: continue
        r = math.sqrt(disc)
        lo = (A*beta - r)/(2*gamma); hi = (A*beta + r)/(2*gamma)
        if lo > hi: lo, hi = hi, lo
        for Ap in range(math.floor(lo)-1, math.ceil(hi)+2):
            dB = 12*(-gamma*Ap*Ap + A*beta*Ap - alpha*A*A)
            if dB < 0: continue
            s = isqrt(dB)
            if s*s != dB: continue
            for sg in (s,-s):
                num = 2*B*Ap + sg
                if num % (2*A): continue
                Bp = num//(2*A)
                if 4*Ap*Bp == 0 and False: pass
                # C' from the linear equation 4A C' = 3 beta + 2 B B' - 4 A' C
                num2 = 3*beta + 2*B*Bp - 4*Ap*C
                if num2 % (4*A): continue
                Cp = num2//(4*A)
                a = (A,B,C); d = (Ap,Bp,Cp)
                # verify S exactly
                ok = True
                for k in (0,1,2):
                    A2,B2,C2 = d
                    num3 = 4*(A+k*A2)*(C+k*C2)-(B+k*B2)**2
                    if num3 != 3*(alpha*k*k+beta*k+gamma): ok = False
                if not ok: continue
                res.append((a,d))
    return res

if __name__ == "__main__":
    alpha, beta, gamma = (int(x) for x in sys.argv[1:4])
    print(f"target S(k) = {alpha}k^2 + {beta}k + {gamma}")
    dh = None
    ps = pencils(alpha,beta,gamma)
    print(f"{len(ps)} integral pencils with a reduced positive definite")
    best = []
    for a,d in ps:
        if dh is None:
            dh = disc_hits(a,d); print("disc-level hits:", dh, len(dh))
        th = true_hits(a,d)
        best.append((len(th), a, d, sorted(th)))
    best.sort(key=lambda t: -t[0])
    seen = set()
    for n,a,d,th in best[:40]:
        print(f"true hits={n:3d}  a={a} d={d}  disc(d)={d[1]**2-4*d[0]*d[2]}  A={th}")

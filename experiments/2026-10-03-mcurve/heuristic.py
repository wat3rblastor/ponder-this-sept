#!/usr/bin/env python3
"""Heuristic first-existence threshold T(n) for an n-term Loeschian AP.

Model (re-derived from scratch; see README).  For a set D of bad primes
dividing d (d = 3*prod(D)*m), the expected number of n-term APs with all
terms <= T is

  E_D(T) = T^2 * prod_{p in D}(1-1/p) / (18 (n-1) prod D)
           * rho_D(T)^n * prod_{q bad, q not in D} corr_q(n)

  rho_D(T) = rho1(T) / prod_{p in D} p/(p+1)        (local-factor removal)
  rho1(T)  = density of Loeschian among integers == 1 mod 3
           = 2*0.638909/sqrt(ln T) * (1 + 0.36/ln T)   [fitted to sieve data]

  corr_q(n): q bad, q not | d.  P(all n terms have even v_q):
     q >= n        : 1 - n/(q+1)
     n/2 <= q < n  : (1 - (n-q)/q)/(q+1)      (>=1 term always hit)
   divided by the (q/(q+1))^n already inside rho1^n.

E(T) = sum over subsets D (forced primes always in) ; solve E(T) = 1.
"""
import math
from itertools import combinations

K2 = 2 * 0.638909

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

BAD = [p for p in primes_upto(300_000) if p % 3 == 2]

def rho1(T, refined=True):
    L = math.log(T)
    r = K2 / math.sqrt(L)
    return r * (1 + 0.36 / L) if refined else r

def forced(n):
    return [p for p in BAD if 2 * p <= n]

def base_of(n):
    b = 3
    for p in forced(n): b *= p
    return b

def corr_q(q, n, refined=True):
    """P(all n terms even v_q) / (q/(q+1))^n  for bad q not dividing d."""
    base = (q / (q + 1.0)) ** n
    if q >= n:
        p = 1 - n / (q + 1.0)
    elif 2 * q >= n:
        p = (1 - (n - q) / q) / (q + 1.0)
    else:
        return None                      # q is forced, cannot be absent
    if not refined:                      # the "stated" formula, verbatim
        p = 1 - n / (q + 1.0)
    return p / base

_CACHE = {}

def coefs(n, refined=True, triangular=True, coprime_a=True):
    """Per-(n,variant) list of (A_D, |D|) with logE_D(T) = A_D + 2lnT + n*ln rho1."""
    key = (n, refined, triangular, coprime_a)
    if key in _CACHE: return _CACHE[key]
    f = forced(n)
    # only q < ~1.4n can ever pay for itself: including q costs a factor q in
    # the number of available d, and gains at most 1/corr_q(n).
    CUT = 1.4 * n
    opts = [q for q in BAD if q not in f and q <= CUT][:8]
    big = [q for q in BAD if q not in f and q not in opts]
    tail = 0.0
    for q in big:
        c = corr_q(q, n, refined)
        tail += math.log(c) if (c and c > 0) else -1e9
    out = []
    for r in range(len(opts) + 1):
        for extra in combinations(opts, r):
            D = f + list(extra)
            A = -math.log(n - 1.0) - math.log(3.0) - math.log(3.0) + tail
            prodD = 1
            for p_ in D: prodD *= p_
            A -= math.log(float(prodD))
            if triangular: A -= math.log(2.0)
            if coprime_a:
                for p_ in D: A += math.log(1 - 1.0 / p_)
            F = 1.0
            for p_ in D: F *= p_ / (p_ + 1.0)
            A += -n * math.log(F)
            bad_cfg = False
            for q in opts:
                if q in extra: continue
                c = corr_q(q, n, refined)
                if c is None or c <= 0: bad_cfg = True; break
                A += math.log(c)
            if bad_cfg: continue
            out.append(A)
    _CACHE[key] = out
    return out

def logE(T, n, refined=True, **kw):
    cs = coefs(n, refined=refined, **kw)
    if not cs: return None
    base = 2 * math.log(T) + n * math.log(rho1(T, refined))
    m = max(cs)
    tot = m + math.log(sum(math.exp(c - m) for c in cs))
    return tot + base

def solve(n, **kw):
    lo, hi = 1e2, 1e60
    for _ in range(120):
        mid = math.sqrt(lo * hi)
        v = logE(mid, n, **kw)
        if v is None or v < 0: lo = mid
        else: hi = mid
    return math.sqrt(lo * hi)

if __name__ == '__main__':
    print(f"{'n':>3} {'base':>9} {'T_refined':>12} {'T_stated':>12}")
    for n in sorted(set(list(range(10, 60, 2)) + [27, 28, 35, 58])):
        tr = solve(n)
        ts = solve(n, refined=False, triangular=False, coprime_a=False)
        print(f"{n:>3} {base_of(n):>9} {tr:>12.3e} {ts:>12.3e}")

#!/usr/bin/env python3
"""Corrected unit scorer for the plan-ranking audit (2026-10-03).

Computes, for a unit (K, shift):
  v_old : exactly what tools/plan_units.py computes (imported, not re-derived)
  v_new : corrected score, see CORRECTIONS below.

Both are "log expected number of 58-term Loeschian APs per stage-1 residue,
up to a unit-independent constant".

CORRECTIONS relative to tools/plan_units.py
  C1 window offset.  plan_units adds `wbar`, the mean window offset w of the
     CUDA kernel's *unreduced* class representative A = Rred + w*MOD
     (src/c/apsearch_cuda.cu:419-422, w < c1+c2 so w can reach ~50).  The CPU
     engine src/c/apsearch.c keeps R reduced mod MOD (a0 = R + (64*shift+b)*MOD,
     R < MOD), so for it w == 0 and the only offset is E[R]/MOD = 0.5.
     v_new uses wbar = 0.5 (option --wbar old keeps the CUDA value).
  C2 term indexing.  plan_units evaluates rho at a + k*d + d, i.e. the terms
     a+d .. a+58d; the unit's terms are a .. a+57d.
  C3 term averaging.  plan_units averages ln ln t over the 6-point grid
     k in {0,11.4,22.8,34.2,45.6,57} and b over 8 points; v_new averages over
     all 58 k and all 64 b exactly.
  C4 divisor bookkeeping.  plan_units scans only its hard-coded bad-prime list
     (q <= 10000) for q | K and good primes p < 1000 for p | K; v_new factors K.
  C5 good-prime size scale.  plan_units hard-codes ln T = 40.5 in the good-prime
     rho correction; v_new uses the unit's own mean ln t.
Audit notes on factors that are NOT corrected because they are
unit-independent (they shift E58 but cannot change the ranking):
  * the missing prod_{p | D0} (1-1/p) ("p does not divide a" for the tier-A
    primes 2,5,11,...,53) -- plan_units has no factor for them at all, while
    rho^58 carries the generic (p/(p+1))^58; the absolute E58 is therefore
    wrong by ~1e14, constant across units.
  * the O(1/q^2) "q | t but q^2 | t" rescue for unfiltered bad primes,
    prod_q (1 + 58/(q(q-58))) ~ 2.2 at q=59... -- same set of primes for every
    unit of the same shape.
  * the truncation of the bad-prime product at q <= 10000.
"""
import math, sys, os, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "tools"))
import plan_units as PU

D0 = PU.D0
NT = 58
EXPO = PU.EXPO


def factor(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def shape_new(K, modcap):
    """Replicates the engine's stage-1 pinning (apsearch.c:190-198) exactly."""
    d = K * D0
    MOD, res, pinned = 30, 4, []
    for q in PU.PR:
        if q <= NT or q % 3 != 2 or d % q == 0:
            continue
        if MOD > 4e18 / q or MOD * q > modcap:
            break
        MOD *= q
        res *= q - NT
        pinned.append(q)
    if len(pinned) < 4:
        return None
    return MOD, res, pinned


def v_new(K, shift, modcap, wbar=0.5, use_dedup=True):
    sh = shape_new(K, modcap)
    if sh is None:
        return None
    MOD, res, pinned = sh
    d = K * D0
    if 64 * MOD * (shift + 1) + 105 * MOD + 57 * d >= 1.8e19:
        return None
    fk = factor(K)
    # --- bad primes: corrections relative to the generic AP factor (q-58)/q
    ly = 0.0
    for q in pinned:                      # pinned: free, factor 1
        ly -= math.log((q - NT) / q)
    for q in fk:                          # q | d: need q does not divide a
        if q > NT and q % 3 == 2:
            ly += math.log((q - 1) / q) - math.log((q - NT) / q)
    # --- size term: exact mean over the 64 window multiples of prod_k rho(t_k)
    acc, lnt_acc = 0.0, 0.0
    for b in range(64):
        a = (64 * shift + b + wbar) * MOD
        m = 0.0
        for k in range(NT):
            m += math.log(math.log(a + k * d))
        acc += math.exp(-EXPO * (m - NT * math.log(39.144)))
        lnt_acc += m / NT
    st = math.log(acc / 64)
    lnT = lnt_acc / 64                    # mean ln(ln t) ... use exp for scale
    lnTbar = math.exp(lnT)
    # --- good primes p | d: no term is divisible by p, which removes the
    #     free good-prime factors and lowers rho
    lgood = 0.0
    for p in fk:
        if p % 3 == 1:
            lgood += -(NT / (p - 1)) * EXPO * math.log(p) / lnTbar
            if use_dedup:
                lgood += math.log(1 - 1 / p)
    return ly + lgood + st, res


def v_old(K, shift):
    r = PU.E_unit_direct(K, shift) if hasattr(PU, "E_unit_direct") else None
    sh = PU.unit_shape(K)
    if sh is None:
        return None
    MOD, res, ly, wbar, lgood = sh
    if 64 * MOD * (shift + 1) + 105 * MOD + 57 * K * D0 >= 1.8e19:
        return None
    st, T = PU.size_term(MOD, wbar, K, shift)
    return ly + lgood + st, res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("units", help="file with 'K shift' lines, or K:shift list")
    ap.add_argument("--modcap", type=float, default=2e16)
    args = ap.parse_args()
    out = []
    for line in open(args.units):
        if not line.strip():
            continue
        K, s = map(int, line.split()[:2])
        o = v_old(K, s)
        n = v_new(K, s, args.modcap)
        n_cuda = v_new(K, s, args.modcap, wbar=PU.unit_shape(K)[3])
        rec = dict(K=K, shift=s,
                   v_old=o[0] if o else None, res_old=o[1] if o else None,
                   v_new=n[0] if n else None, res_new=n[1] if n else None,
                   v_new_cudaw=n_cuda[0] if n_cuda else None,
                   wbar=PU.unit_shape(K)[3])
        out.append(rec)
        print(json.dumps(rec))

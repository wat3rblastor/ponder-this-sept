#!/usr/bin/env python3
"""Re-derive the pin-set trade-off table from scratch (no inherited numbers).

Definitions used (all re-derived here, nothing trusted):
  bad prime      p = 2 (mod 3)
  n Loeschian    iff every bad prime divides n to an even power
  AP             a, a+d, ..., a+(L-1)d   with L = 58 terms
  forced         3 | d, a = 1 (mod 3); every bad p with 2p <= L must divide d
  D0             working step base
  admissible     for a bad prime q not dividing d: q divides no term iff
                 a = j*d (mod q) with 1 <= j <= q-L.  (q-L residues out of q.)
"""
from math import prod, log, sqrt

L = 58

def is_bad(p):
    return p % 3 == 2

def primes(n):
    s = [True]*(n+1); s[0]=s[1]=False
    for i in range(2,int(n**.5)+1):
        if s[i]:
            for j in range(i*i,n+1,i): s[j]=False
    return [i for i,b in enumerate(s) if b]

P = primes(5000)
BAD = [p for p in P if is_bad(p)]
print("bad primes <= 140:", [p for p in BAD if p <= 140])

# --- forced divisors of d -------------------------------------------------
forced = [p for p in BAD if 2*p <= L]
print("bad primes forced to divide d (2p <= 58):", forced)
# The working base also absorbs 41,47,53 (each would otherwise hit 2 terms
# only if 2p<=58; for 58<2p<=116 a bad prime q not dividing d hits exactly one
# term unless excluded -- including them in d is a choice, not forced).
D0 = 3 * prod([2,5,11,17,23,29,41,47,53])
print("D0 =", D0, " (3 * 2*5*11*17*23*29*41*47*53)")
assert D0 == 382160924970, D0
print("term floor at K=1 (57*D0) =", 57*D0)

# --- the F factor: bad primes that can never divide any term --------------
# A bad prime p | d with p coprime to a divides no term at all.
# A bad prime q with q-L <= 0, i.e. q < L, that does not divide d WOULD have to
# hit a term; so all bad q < 58 must divide d (that is the 2p<=58 rule extended:
# in fact q <= 57 => q-L <= 1).  We take the set of bad primes excluded from all
# terms = divisors of d (coprime to a) + pinned primes q in Q.
def F_of(excluded):
    return prod([p/(p+1) for p in excluded])

K_LR = 0.638909405  # Landau-Ramanujan-type constant for x^2+xy+y^2 (inherited)

def rho(T, excluded):
    """INHERITED model (not verified here): per-term pass probability."""
    return 2*K_LR/(sqrt(log(T))*F_of(excluded))

# --- the pin-set table ----------------------------------------------------
cand = [q for q in BAD if q > L and q not in forced]   # 59,71,83,89,101,...
print("\npinnable primes:", cand[:10])

MODSMALL = 3*2*5   # the a mod 3 / a odd / 5 ! a part of the engine wheel
rows = []
for k in range(1, 8):
    Q = cand[:k]
    M = MODSMALL*prod(Q)
    delta = prod([(q-L)/q for q in Q])
    # the engine's own floor: it lifts a residue mod M by b*M, b=0..63
    engine_floor = max(57*D0, 64*M)
    nat_floor = 57*D0
    excl = [2,5,11,17,23,29,41,47,53] + Q
    rows.append(dict(Q=Q, M=M, delta=delta, inv=1/delta,
                     engine_floor=engine_floor, nat_floor=nat_floor,
                     hit_engine=rho(engine_floor, excl)**L,
                     hit_nat=rho(nat_floor, excl)**L))

print("\n%-5s %-11s %-11s %-10s %-11s %-11s %-11s" %
      ("maxq","M","1/delta","64*M","engine floor","P58@engine","P58@2.2e13"))
for r in rows:
    print("%-5d %-11.3g %-11.0f %-10.3g %-11.3g %-11.4g %-11.4g" %
          (r["Q"][-1], r["M"], r["inv"], 64*r["M"], r["engine_floor"],
           r["hit_engine"], r["hit_nat"]))

print("\n-- ratio of P58 at the natural floor (2.18e13) to P58 at the engine floor --")
base = rows[-1]["hit_engine"]           # what we actually run today (pin to 113)
for r in rows:
    print("  pin 59..%-4d  P58@nat/P58@today = %9.1f   (1/delta = %.0f)"
          % (r["Q"][-1], r["hit_nat"]/base, r["inv"]))

print("\n-- expected hits per candidate emitted, per unit enumeration cost --")
print("   today: 1 candidate per ~1 loop iteration, P58 = %.4g" % base)
for r in rows:
    print("  pin 59..%-4d at term size 2.18e13: P58 = %.4g" % (r["Q"][-1], r["hit_nat"]))

# --- supply: how many admissible a are there below X at all? --------------
X = 57*D0
print("\n-- supply of admissible a below X = 57*D0 = %.3g --" % X)
# density from the small part: a = 1 mod 3 (1/3), a odd (1/2), 5 ! a (4/5)
SMALL_DENS = (1/3)*(1/2)*(4/5)
print("   small-part density (a=1 mod 3, odd, 5!a) = %.5f  => 1/%.1f"
      % (SMALL_DENS, 1/SMALL_DENS))
for r in rows:
    n = r["delta"]*SMALL_DENS*X
    print("  pin 59..%-4d : 1/delta_total = %8.0f, admissible a < X = %.4g,"
          " expected 58s = %.3g"
          % (r["Q"][-1], 1/(r["delta"]*SMALL_DENS), n, n*r["hit_nat"]))

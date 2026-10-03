#!/usr/bin/env python3
"""Where the 5.33x goes, and what term size each mode can actually reach.

n = 58.  Forced bad primes (2q <= n): 2,5,11,17,23,29  -> always in d.
Bracket primes 41,47,53: either IN d, or SQUARE (q^2 divides exactly one term).

For a fixed d, the density of a that survive every bad prime up to B is
    rho(d) = prod_{r | d, r bad} (1 - 1/r)                 # r must not divide a
           * prod_{r bad, r > n, r !| d} (1 - n/r)         # r must miss window
           * prod_{q square} (2q - n)/q^2                  # q^2 hits one term
and the number of d of a given mode below D is D / dmin(mode).

Counting (a, d) pairs with last term <= X:
    N(X) = (X / ((n-1) * dmin)) * X * rho
This is the quantity the "5.33x" refers to (at fixed X), and 1/rho is the
smallest a at which the mode has ANY admissible candidate at all -- the number
that the "term size drops by 1e5" claim leaves out.
"""
n = 58
B = 2000


def primes(b):
    return [x for x in range(2, b) if all(x % y for y in range(2, int(x**.5) + 1))]


BAD = [p for p in primes(B) if p % 3 == 2]
FORCED = [2, 5, 11, 17, 23, 29]
BRACKET = [41, 47, 53]
BASE = 3 * 2 * 5 * 11 * 17 * 23 * 29


def mode(sq):
    d = BASE
    for q in BRACKET:
        if q not in sq:
            d *= q
    rho = 1.0
    for r in BAD:
        if d % r == 0:
            rho *= 1 - 1.0 / r          # r | d: need r !| a
        elif r in sq:
            rho *= (2.0 * r - n) / (r * r)   # square option
        elif r > n:
            rho *= 1 - float(n) / r     # must miss the window
        else:
            pass                        # r < n and r !| d: impossible, excluded
    return d, rho


def N(X, d, rho):
    return (X / ((n - 1.0) * d)) * X * rho


base_d, base_rho = mode([])
print(f"{'mode':<16} {'dmin':>14} {'rho':>12} {'1/rho':>10} "
      f"{'(n-1)*dmin':>12} {'floor':>10} {'N(X)/N_alld(X)':>15}")
rows = [[], [41], [47], [53], [41, 47], [41, 53], [47, 53], [41, 47, 53]]
tot = 0.0
for sq in rows:
    d, rho = mode(sq)
    ratio = N(1e16, d, rho) / N(1e16, base_d, base_rho)
    tot += ratio
    floor = max((n - 1) * d, 1 / rho)
    print(f"{('sq=' + ','.join(map(str, sq)) if sq else 'all-in-d'):<16} "
          f"{d:>14} {rho:>12.4g} {1/rho:>10.4g} {(n-1)*d:>12} "
          f"{floor:>10.4g} {ratio:>15.4f}")
print(f"\nsum over all 8 modes / all-in-d = {tot:.4f}"
      "   <-- this is the '5.33x'")
print("per-prime factors 1 + (2q-n)/q:",
      ", ".join(f"q={q}: {1 + (2.0*q-n)/q:.4f}" for q in BRACKET),
      f"-> product {(1+24/41.)*(1+36/47.)*(1+48/53.):.4f}")
print("pure square-family share of the solution space: "
      f"{(24/41.)*(36/47.)*(48/53.):.4f} x all-in-d")

print("\nsmallest last term with N(X) = 1 (the honest 'reachable term size'):")
for sq in rows:
    d, rho = mode(sq)
    X = ((n - 1.0) * d / rho) ** 0.5
    print(f"  {('sq=' + ','.join(map(str, sq)) if sq else 'all-in-d'):<16} "
          f"X = {X:.4g}")

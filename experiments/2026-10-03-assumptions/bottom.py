"""How many n-term Loeschian APs exist with all terms <= X?

The forced structure (GOAL.md) is: 3 | d, and every bad prime q with 2q <= n
divides d.  That is ALL that is forced.  base(n) = 3 * prod(bad q <= n/2).
For n = 58 base = 3*2*5*11*17*23*29 = 3741870, NOT
D0 = 3*2*5*11*17*23*29*41*47*53 = 382160924970 (41,47,53 are optional: each may
instead take q^2 on its single hit term).

Model.  Family F(n,X) = {(a,d) : d = base*m, a = 1 mod 3, a != 0 mod q for the
forced bad q, a+(n-1)d <= X}.  |F| is counted exactly by the formula below.
The per-term Loeschian rate rho is measured by Monte Carlo on exactly that
sampling distribution, and E(X) = |F| * rho^n (terms treated as independent).

Validation: n = 35 has a PROVED unique minimum (records.json, last term
311958331, exactly one AP at or below it).  If the model says E_35(3.1e8) ~ 1
it is calibrated, and its n = 58 numbers can be believed to an order of
magnitude.
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))
from loeschian import is_loeschian, is_prime


def bad_primes(lo, hi):
    return [q for q in range(lo, hi + 1) if is_prime(q) and q % 3 == 2]


def base_of(n):
    b = 3
    for q in bad_primes(2, n // 2):
        b *= q
    return b


def count_pairs(n, X, base, forced):
    """|F(n,X)| exactly (sum over d of the count of admissible a)."""
    dens = 1.0 / 3.0
    for q in forced:
        dens *= (q - 1) / q
    mmax = (X - 1) // ((n - 1) * base)
    tot = 0.0
    for m in range(1, mmax + 1):
        tot += (X - (n - 1) * base * m) * dens
    return tot, mmax


def rho_mc(n, X, base, forced, trials, seed=1):
    """Per-term Loeschian rate, sampling (a,d) uniformly from F(n,X)."""
    rng = random.Random(seed)
    mmax = (X - 1) // ((n - 1) * base)
    good = tot = 0
    runs = []
    while tot < trials:
        m = rng.randint(1, mmax)
        d = base * m
        hi = X - (n - 1) * d
        a = rng.randint(1, hi)
        a -= (a - 1) % 3                      # a = 1 mod 3
        if a < 1:
            continue
        if any(a % q == 0 for q in forced):
            continue
        hits = [is_loeschian(a + k * d) for k in range(n)]
        good += sum(hits)
        tot += n
        # longest run, for a sanity check on clustering
        r = c = 0
        for h in hits:
            c = c + 1 if h else 0
            r = max(r, c)
        runs.append(r)
    return good / tot, runs


def report(n, X, base=None, label=''):
    forced = bad_primes(2, n // 2)
    if base is None:
        base = base_of(n)
    pairs, mmax = count_pairs(n, X, base, forced)
    trials = 20000 if n > 40 else 20000
    rho, runs = rho_mc(n, X, base, forced, trials)
    E = pairs * rho ** n
    print(f"n={n:3d} X=1e{len(str(X))-1:<3d} base={base:<14d} {label}")
    print(f"      #d={mmax:<8d} |F|={pairs:.3e}  rho={rho:.4f}  "
          f"rho^n={rho**n:.3e}   E(X)={E:.3e}   E*0.82={E*0.82:.3e}")
    return E


if __name__ == '__main__':
    print("--- CALIBRATION: n=35, where the minimum is proved unique ---")
    for X in (1e8, 3.1e8, 1e9, 1e10):
        report(35, int(X))
    print("\n--- n=58, truly forced base 3741870 (41,47,53 by square option) ---")
    for X in (1e10, 1e11, 1e12, 1e13, 1e14, 1e16):
        report(58, int(X))
    print("\n--- n=58, conventional D0 (41,47,53 all in d) ---")
    D0 = 382160924970
    for X in (1e13, 1e14, 1e16, 1e18, 1e20):
        report(58, int(X), base=D0, label='(D0 family)')

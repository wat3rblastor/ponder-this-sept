"""Ground truth for the counting model, via a fast exact Loeschian bitmap.

Builds a Loeschian indicator for 0..X with X < 4472^2 using
  v_q(t) mod 2 = XOR over k>=1 of [q^k | t]
for bad primes q <= sqrt(X), plus the single large cofactor.

Then, for the family  d = base(n)*m, a = 1 mod 3, a != 0 mod (bad q <= n/2):
  * EXACT count of n-term APs with last term <= X  (ground truth),
  * the independence model rho^n,
  * the per-prime local-density model of count.py.

Cross-checked: the bitmap is compared against src/loeschian.is_loeschian on
200k random points before use.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))
from loeschian import is_loeschian
from count import bad_primes, base_of, local_product, rho_tail

X = 19_000_000           # < 4472^2 = 19998784, so the cofactor is 1 or one prime


def build(X):
    n = X + 1
    S = int(X ** 0.5)
    sieve = np.ones(S + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(S ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    primes = np.flatnonzero(sieve)
    res = np.arange(n, dtype=np.int64)
    par = np.zeros(n, dtype=np.uint8)
    for p in primes:
        p = int(p)
        pk = p
        while pk <= X:
            res[pk::pk] //= p
            if p % 3 == 2:
                par[pk::pk] ^= 1
            pk *= p
    # leftover cofactor: 1 or a single prime > S, valuation 1
    big = res > 1
    par[big] ^= ((res[big] % 3) == 2).astype(np.uint8)
    loe = par == 0
    loe[0] = True
    return loe


def verify(loe, trials=200000, seed=11):
    rng = np.random.default_rng(seed)
    pts = rng.integers(1, X, trials)
    bad = [int(t) for t in pts if bool(loe[t]) != is_loeschian(int(t))]
    print(f"bitmap cross-check vs src/loeschian.is_loeschian on {trials} points: "
          f"{len(bad)} mismatches {bad[:5]}")
    assert not bad


def admissible_density(n):
    dens = 1 / 3
    for q in bad_primes(2, n // 2):
        dens *= (q - 1) / q
    return dens


def exact_count(loe, n, X):
    """Exact number of (a,d) in the family with last term <= X, plus per-term rate."""
    base = base_of(n)
    forced = bad_primes(2, n // 2)
    mmax = (X - 1) // ((n - 1) * base)
    total_pairs = 0
    total_terms = 0
    total_ok = 0
    found = []
    for m in range(1, mmax + 1):
        d = base * m
        hi = X - (n - 1) * d
        a = np.arange(1, hi + 1, 3, dtype=np.int64)   # a = 1 mod 3
        for q in forced:
            a = a[a % q != 0]
        total_pairs += a.size
        mask = loe[a]
        total_terms += a.size
        total_ok += int(mask.sum())
        for k in range(1, n):
            mask &= loe[a + k * d]
            if not mask.any():
                break
        else:
            for aa in a[mask]:
                found.append((int(aa), int(d)))
    return total_pairs, total_ok / max(total_terms, 1), found, mmax, base


if __name__ == '__main__':
    print(f"X = {X}; building bitmap ...")
    loe = build(X)
    verify(loe)
    print(f"Loeschian density below X: {loe.sum()/X:.4f}\n")
    rt = rho_tail(X, 1000, 3000)
    print(f"rho_tail(X, Q=1000) = {rt:.4f}\n")
    for n in (10, 12, 14, 16, 18, 20, 22):
        pairs, rho, found, mmax, base = exact_count(loe, n, X)
        indep = pairs * rho ** n
        loc = local_product(n, 1000) / admissible_density(n)
        mod2 = pairs * loc * rt ** n
        print(f"n={n:3d} base={base:<7d} #d={mmax:<6d} pairs={pairs:.4e} "
              f"rho1={rho:.4f}")
        print(f"      EXACT count = {len(found):<6d}   model I (rho^n) = {indep:10.4g}"
              f"  ratio={len(found)/indep if indep else 0:8.3g}")
        print(f"                                model II (local)= {mod2:10.4g}"
              f"  ratio={len(found)/mod2 if mod2 else 0:8.3g}")
        if found[:3]:
            print(f"      examples {found[:3]}")
        sys.stdout.flush()

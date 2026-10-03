"""Exact-local-density count of n-term Loeschian APs with all terms <= X.

E(X) = |F(n,X)| * prod_{bad q <= Q} P_q(n) * rho_Q(X)^n

|F(n,X)| = #{(a,d): d = base(n)*m, a = 1 mod 3, a + (n-1)d <= X}
           base(n) = 3 * prod(bad q <= n/2)   <-- all that is forced

P_q(n) = P over (a, d mod q^inf), q | d allowed, of "every term has even v_q".
  q <= n/2 (forced q | d):          P_q = (q-1)/q
  n/2 < q <= n:  s = 2q-n single-hit classes; classes with two hits are
      impossible when q !| d (their difference is q*d).
      P_q = (1/q)(q-1)/q + ((q-1)/q)(s/q)(1/(q+1))
  q > n:
      P_q = (1/q)(q-1)/q + ((q-1)/q)[(1 - n/q) + (n/q)/(q+1)]
  using P(v_q even and > 0 | q | t) = sum_{j>=1} q^-(2j-1)(1-1/q) = 1/(q+1).

rho_Q(X) = P(every bad prime q > Q has even valuation in t), t ~ X: Monte Carlo
by exact factorisation.  Terms are treated as independent ONLY for q > Q, where
the error term is n^2/(2 q^2) summed over bad q > Q (< 0.12 for Q = 1000).

Calibration target: n = 35 has a proved unique minimum with last term
311958331 (records.json), so E_35(3.1e8) must come out near 1.
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))
from loeschian import factorize, is_prime


def bad_primes(lo, hi):
    return [q for q in range(max(lo, 2), hi + 1) if is_prime(q) and q % 3 == 2]


def base_of(n):
    b = 3
    for q in bad_primes(2, n // 2):
        b *= q
    return b


def P_q(q, n):
    if 2 * q <= n:
        return (q - 1) / q                     # q | d forced
    if q <= n:
        s = 2 * q - n                          # single-hit classes
        return (1 / q) * (q - 1) / q + ((q - 1) / q) * (s / q) * (1 / (q + 1))
    return (1 / q) * (q - 1) / q + ((q - 1) / q) * ((1 - n / q) + (n / q) / (q + 1))


def local_product(n, Q):
    p = 1.0
    for q in bad_primes(2, Q):
        p *= P_q(q, n)
    return p


_RHO_CACHE = {}


def rho_tail(X, Q, trials=4000, seed=7):
    """P(all bad primes > Q have even valuation), t uniform near X."""
    key = (X, Q, trials, seed)
    if key in _RHO_CACHE:
        return _RHO_CACHE[key]
    rng = random.Random(seed)
    ok = 0
    for _ in range(trials):
        t = rng.randrange(X // 2, X)
        f = factorize(t)
        if all(not (p % 3 == 2 and p > Q and e % 2) for p, e in f.items()):
            ok += 1
    r = ok / trials
    _RHO_CACHE[key] = r
    return r


def pairs(n, X, base):
    mmax = (X - 1) // ((n - 1) * base)
    M, s = float(mmax), float((n - 1) * base)
    tot = (M * X - s * M * (M + 1) / 2) / 3.0
    return tot, mmax


def E(n, X, base=None, Q=1000, label='', trials=4000):
    if base is None:
        base = base_of(n)
    F, mmax = pairs(n, X, base)
    L = local_product(n, Q)
    # a forced base larger than base_of(n) means those primes are in d: their
    # local factor becomes (q-1)/q instead of the mixed value.
    extra = 1.0
    for q in bad_primes(2, n):
        if base % q == 0 and 2 * q > n:
            extra *= ((q - 1) / q) / P_q(q, n)
    r = rho_tail(X, Q, trials)
    val = F * L * extra * r ** n
    print(f"  n={n:3d} X=1e{len(str(X))-1:<3d} base={base:<13d} #d={mmax:<9d} "
          f"|F|={F:.3e} loc={L*extra:.3e} rho={r:.4f} rho^n={r**n:.3e}  E={val:.3e}  {label}")
    return val


if __name__ == '__main__':
    print("--- CALIBRATION n=35 (exhaustively proved: exactly ONE AP with last term <= 3.1e8) ---")
    for X in (1e8, 311958331, 1e9, 1e10):
        E(35, int(X))
    print("--- CALIBRATION n=36,43,47,55: records.json term sizes, E should be >> 1 ---")
    E(36, 1383069883)
    E(43, 7468034210927227)
    E(47, 6249948665490457)
    E(55, 22645664238935173)
    print("\n--- n=58, truly forced base 3741870 (41,47,53 by square option) ---")
    for X in (1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e18, 1e20):
        E(58, int(X))
    print("\n--- n=58, conventional D0 = 3*2*5*11*17*23*29*41*47*53 ---")
    D0 = 382160924970
    for X in (1e14, 1e15, 1e16, 1e18, 1e20):
        E(58, int(X), base=D0, label='(D0)')
    print("\n--- n=57 for reference, forced base ---")
    for X in (1e11, 1e12, 1e13):
        E(57, int(X))

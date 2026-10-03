#!/usr/bin/env python3
"""Residue / hit-index / d-structure forensics on the record APs."""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from loeschian import bad_primes, factorize  # noqa: E402

RECORDS = [
    ("55a", 11687581876345393, 202927451159070, 55),
    ("55b", 296246969176050787, 11950172123811900, 55),
    ("47", 2646171143023357, 78342989618850, 47),
    ("43", 5413537078288507, 48916598396160, 43),
    ("35opt", 219830911, 2709630, 35),
]

QMAX = 1000


def fmt(f):
    return " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))


def main():
    allhits = []
    for name, a, d, n in RECORDS:
        fd = factorize(d)
        badd = sorted(p for p in fd if p % 3 == 2)
        goodd = sorted(p for p in fd if p % 3 != 2)
        forced = [p for p in bad_primes(n) if 2 * p <= n]
        bp_lt_n = [p for p in bad_primes(n + 1)]
        nxt = next(p for p in bad_primes(10000) if p > max(badd))
        print(f"=== {name}  n={n}")
        print(f"  d = {fmt(fd)}")
        print(f"  bad primes | d        : {badd}")
        print(f"  forced (2p<=n)        : {forced}")
        print(f"  extra bad primes in d : {[p for p in badd if p not in forced]}")
        print(f"  bad primes <= n not | d: {[p for p in bp_lt_n if p not in badd]}")
        print(f"  good primes | d       : {goodd}")
        print(f"  v: {[(p, fd[p]) for p in badd]}")
        print(f"  largest bad p | d = {max(badd)}; next bad prime = {nxt} "
              f"({'>' if nxt > n - 1 else '<='} n-1={n-1})")
        # hit indices for bad primes not dividing d
        exp = 0.0
        obs = 0
        rows = []
        for q in bad_primes(QMAX):
            if d % q == 0:
                continue
            k0 = (-a * pow(d, -1, q)) % q
            hits = [k for k in range(n) if (k - k0) % q == 0]
            exp += min(1.0, n / q) if q >= n else n / q
            if hits:
                obs += 1
                vs = [(k, max(e for p, e in factorize(a + k * d).items() if p == q)
                       if (a + k * d) % q == 0 else 0) for k in hits]
                rows.append((q, k0, vs))
        print(f"  bad q<{QMAX}, q!|d: #with a hit in window = {obs}, "
              f"expected (uniform k0) = {exp:.2f}")
        for q, k0, vs in rows:
            print(f"     q={q} k0={k0} v_q(term)={vs}")
        allhits.append((name, n, obs, exp, rows))
    # aggregate binomial test
    tot_obs = sum(x[2] for x in allhits)
    tot_exp = sum(x[3] for x in allhits)
    sd = math.sqrt(sum(
        sum(min(1.0, x[1] / q) * (1 - min(1.0, x[1] / q))
            for q in bad_primes(QMAX)) for x in allhits))
    print(f"\nAGGREGATE hit-prime count: obs={tot_obs} exp={tot_exp:.2f} "
          f"(ratio {tot_obs / tot_exp:.2f}, crude sd~{sd:.2f}, "
          f"z={(tot_obs - tot_exp) / sd:.2f})")

    # relations between records
    print("\n=== relations")
    import itertools
    for (n1, a1, d1, _), (n2, a2, d2, _) in itertools.combinations(RECORDS, 2):
        print(f"  {n1} vs {n2}: gcd(d)={math.gcd(d1, d2)} gcd(a)={math.gcd(a1, a2)} "
              f"gcd(a1,d2)={math.gcd(a1, d2)} gcd(a2,d1)={math.gcd(a2, d1)}")
    for name, a, d, n in RECORDS:
        print(f"  {name}: a/d = {a / d:.3f}  (centre index {a / d:.1f}); "
              f"a mod d = {a % d}; gcd(a,d)={math.gcd(a, d)}")
        fa = factorize(a)
        print(f"      a = {fmt(fa)}")


if __name__ == "__main__":
    main()

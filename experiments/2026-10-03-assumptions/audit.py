"""Structure audit of every known record progression.

For each (a, d, n): which bad primes q with n/2 < q <= n are NOT in d (so the
progression must be paying q^2 on a single hit term -- the "square option"),
and which bad primes q > n actually hit a term with q^2.

Run: python3 experiments/2026-10-03-assumptions/audit.py
"""
import json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))
from loeschian import is_loeschian, factorize, is_bad_prime


def bad_primes_upto(x):
    out = []
    n = 2
    while n <= x:
        from loeschian import is_prime
        if is_prime(n) and n % 3 == 2:
            out.append(n)
        n += 1
    return out


def audit(a, d, n, label):
    assert all(is_loeschian(a + k * d) for k in range(n)), label
    print(f"\n=== {label}: n={n} a={a} d={d}")
    print(f"    terms {a} .. {a+(n-1)*d}  (~1e{len(str(a+(n-1)*d))-1})  d={factorize(d)}")
    fd = factorize(d)
    # which bad primes q <= n are absent from d
    absent = [q for q in bad_primes_upto(n) if q not in fd]
    print(f"    bad q <= n absent from d: {absent}")
    for q in absent:
        hits = [k for k in range(n) if (a + k * d) % q == 0]
        info = []
        for k in hits:
            v = 0
            t = a + k * d
            while t % q == 0:
                t //= q
                v += 1
            info.append((k, v))
        print(f"      q={q}: hit indices/valuations {info}  (needs even valuation)")
    # bad primes > n that hit a term to an even power (>=2)
    sq = {}
    for k in range(n):
        t = a + k * d
        for p, e in factorize(t).items():
            if p % 3 == 2 and p > n and e >= 2:
                sq.setdefault(p, []).append((k, e))
    print(f"    bad q > n appearing squared in some term: "
          f"{ {p: v for p, v in sorted(sq.items())} }")
    return absent, sq


def main():
    recs = json.load(open(os.path.join(os.path.dirname(__file__), '..', '..', 'records.json')))
    items = []
    g = recs['g3_longest']
    items.append((g['a'], g['d'], g['n'], 'g3_longest (best known here)'))
    g = recs['g1_min_last_term_35']
    items.append((g['a'], g['d'], g['n'], 'n=35 PROVED MINIMAL'))
    for h in recs.get('history', []):
        if 'a' in h:
            items.append((h['a'], h['d'], h['n'], f"history n={h['n']}"))
    stats = []
    for a, d, n, lab in items:
        stats.append((lab, n) + audit(a, d, n, lab))
    print("\n--- summary: square option used? ---")
    for lab, n, absent, sq in stats:
        print(f"  {lab:34s} n={n:3d} square-option primes in (n/2,n]: {absent}"
              f"   bad squares > n: {sorted(sq)}")


if __name__ == '__main__':
    main()

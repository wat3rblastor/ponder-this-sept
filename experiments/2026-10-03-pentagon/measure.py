#!/usr/bin/env python3
"""Matched measurement: pentagonal automatic-terms family vs. generic baseline.

Same d, same term size, same congruence conditions; only the shape of a differs.
Usage: python3 experiments/2026-10-03-pentagon/measure.py [N]
Requires build/apsearch (make apsearch) for the exact Loeschian test.
"""
import subprocess, sys, random
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
m = 5 * 11 * 17 * 23 * 29          # 623645: puts 5,11,17,23,29 into d = 24m^2
assert m % 2 == 1 and m % 3 != 0
d = 24 * m * m
AUTO = {k for k in range(58) if isqrt(24 * k + 1) ** 2 == 24 * k + 1}

xs, x = [], 2
while len(xs) < N:                  # x even -> a odd; x coprime to the forced primes
    if all(x % p for p in (5, 11, 17, 23, 29)):
        xs.append(x)
    x += 2
rows = [3 * x * x + m * m for x in xs]

rng = random.Random(7)
lo, hi = rows[0], rows[-1]
base = []
while len(base) < N:                # generic a, same magnitude and congruences
    a = rng.randrange(lo, hi + 1)
    if a % 2 and a % 3 == 1 and all(a % p for p in (5, 11, 17, 23, 29)):
        base.append(a)

terms = [a + k * d for a in rows + base for k in range(58)]
p = subprocess.run([str(ROOT / "build/apsearch"), "--isl"], input="\n".join(map(str, terms)),
                   capture_output=True, text=True)
v = [int(l.split()[1]) for l in p.stdout.splitlines()]


def stats(label, chunk, skip=frozenset()):
    ok = tot = na_ok = na_tot = 0
    longest = [0] * 60
    for i in range(N):
        r = chunk[i * 58:(i + 1) * 58]
        best = cur = 0
        for k, b in enumerate(r):
            cur = cur + 1 if b else 0
            best = max(best, cur)
            if k not in skip:
                na_tot += 1; na_ok += b
        longest[best] += 1
        ok += sum(r); tot += 58
    ge = [sum(longest[j:]) for j in range(60)]
    print(f"{label}: per-term {ok/tot:.4f}  non-automatic {na_ok/na_tot:.4f}  " +
          "  ".join(f"P(>={L})={ge[L]/N:.4f}" for L in (10, 12, 14, 16)))
    return na_ok / na_tot


pe = stats("PENTAGON", v[:N * 58], AUTO)
ba = stats("BASELINE", v[N * 58:], frozenset())
print(f"\nautomatic positions: {sorted(AUTO)} ({len(AUTO)} of 58)")
print(f"extrapolated P(58-run) ratio = {pe**45 / ba**58:.3g}x per candidate")

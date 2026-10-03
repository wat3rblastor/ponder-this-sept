#!/usr/bin/env python3
"""Analytic power analysis for the 'automatic family' feature.

3*c*d*t is a perfect square  <=>  t = S_c * u^2  with S_c = squarefree kernel
of 3*c*d.  So the automatic-family values form a sparse set of density
1/(2*sqrt(X*S_c)) near X.  This computes the exact expected number of
automatic-family terms in each record's window under the null, i.e. the power
the Monte-Carlo test could ever have had, and the exact observed number.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from loeschian import factorize  # noqa: E402

RECORDS = [
    ("55a", 11687581876345393, 202927451159070, 55),
    ("55b", 296246969176050787, 11950172123811900, 55),
    ("47", 2646171143023357, 78342989618850, 47),
    ("43", 5413537078288507, 48916598396160, 43),
    ("35opt", 219830911, 2709630, 35),
]
CMAX = 60


def sqfree(n):
    out = 1
    for p, e in factorize(n).items():
        if e % 2:
            out *= p
    return out


print(f"{'rec':<6} {'E[auto terms] (null)':>22} {'observed':>9} "
      f"{'min S_c':>14} {'last term':>22}")
for name, a, d, n in RECORDS:
    last = a + (n - 1) * d
    exp = 0.0
    obs = 0
    minS = None
    for c in range(1, CMAX + 1):
        S = sqfree(3 * c * d)
        minS = S if minS is None else min(minS, S)
        # expected count of t in the AP with t = S*u^2:
        # density of {S u^2} near X is 1/(2 sqrt(X S)); n terms
        exp += sum(1.0 / (2 * math.sqrt((a + k * d) * S)) for k in range(n))
        # exact check
        for k in range(n):
            t = a + k * d
            if t % S == 0:
                u = math.isqrt(t // S)
                if u * u * S == t:
                    obs += 1
    print(f"{name:<6} {exp:>22.3e} {obs:>9} {minS:>14} {last:>22}")

print("\nWith E[auto terms] this small the Monte-Carlo test has no power at "
      "all; the analytic null is the test, and the records are consistent "
      "with it (0 observed).")
print("For an automatic-family term to exist at all the window would have to "
      "reach down to ~S_c, i.e. terms of order sqfree(3*c*d) * u^2 with small "
      "u -- that means terms far smaller than d^1, which no long AP has.")

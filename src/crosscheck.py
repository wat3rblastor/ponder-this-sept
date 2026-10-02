#!/usr/bin/env python3
"""Independent, CONSTRUCTIVE cross-check of a Loeschian AP (GOAL.md §8 step 2).

`src/verify.py` decides "is this Loeschian?" from the bad-prime criterion: every
prime p = 2 (mod 3) to an even power. That criterion is the thing our whole
search is built on, so verifying a record with it alone is close to circular.

This script does not use the criterion at all. For each term t it *constructs*
integers x, y with

    x^2 + x*y + y^2 == t

and checks that identity by direct integer arithmetic. A representation is
positive proof that t is Loeschian, independent of any theory about bad primes:

  1. factor t (re-multiplied to confirm the factorization),
  2. represent each prime power as the norm of an Eisenstein integer:
       - p = 3                      -> 3 = N(2 + w)
       - p = 1 (mod 3)              -> Cornacchia: p = u^2 + 3v^2, then
                                       p = N((u+v) + 2v*w)
       - p = 2 (mod 3), even exp 2f -> p^2 = N(p + 0*w), used f times
  3. multiply the Eisenstein integers together (the norm is multiplicative),
  4. read off (x, y) and check x^2 + x*y + y^2 == t exactly.

Z[w] with w = (-1 + sqrt(-3))/2 has N(a + b*w) = a^2 - a*b + b^2, so the
reported pair is (x, y) = (a, -b).

If a prime p = 2 (mod 3) ever appears to an ODD power, no representation exists
and this script reports FAIL for that term -- which is exactly the independent
check we want.

Usage:  crosscheck.py A D N
Exit status 0 iff every term got a verified representation.
"""

from __future__ import annotations

import sys
from math import isqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from loeschian import factorize, is_prime  # noqa: E402


def sqrt_mod(a: int, p: int) -> int | None:
    """Tonelli-Shanks: a square root of a mod odd prime p, or None."""
    a %= p
    if a == 0:
        return 0
    if p == 2:
        return a
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t = t * c % p
        r = r * b % p
    return r


def cornacchia_3(p: int) -> tuple[int, int] | None:
    """Solve u^2 + 3*v^2 == p for a prime p = 1 (mod 3). Verified before return."""
    if p == 3:
        return (0, 1)
    r = sqrt_mod(-3 % p, p)
    if r is None:
        return None
    lim = isqrt(p)
    for r0 in (r, p - r):
        a, b = p, r0
        while b > lim:
            a, b = b, a % b
        rem = p - b * b
        if rem >= 0 and rem % 3 == 0:
            v2 = rem // 3
            v = isqrt(v2)
            if v * v == v2 and b * b + 3 * v * v == p:
                return (b, v)
    # last resort for tiny p
    v = 0
    while 3 * v * v <= p:
        u2 = p - 3 * v * v
        u = isqrt(u2)
        if u * u == u2:
            return (u, v)
        v += 1
    return None


def eis_mul(z: tuple[int, int], w: tuple[int, int]) -> tuple[int, int]:
    """(a + b*w) * (c + d*w) with w^2 = -1 - w."""
    a, b = z
    c, d = w
    return (a * c - b * d, a * d + b * c - b * d)


def eis_norm(z: tuple[int, int]) -> int:
    a, b = z
    return a * a - a * b + b * b


def represent(t: int) -> tuple[int, int] | None:
    """Return (x, y) with x^2 + x*y + y^2 == t, or None if t is not a norm."""
    if t == 0:
        return (0, 0)
    if t == 1:
        return (1, 0)
    f = factorize(t)
    prod = 1
    for p, e in f.items():
        if not is_prime(p):
            return None
        prod *= p**e
    if prod != t:
        raise AssertionError(f"factorization of {t} does not re-multiply")

    z = (1, 0)
    for p, e in sorted(f.items()):
        if p == 3:
            unit = (2, 1)                 # N(2 + w) = 4 - 2 + 1 = 3
            reps = e
        elif p % 3 == 1:
            uv = cornacchia_3(p)
            if uv is None:
                return None
            u, v = uv                     # p = u^2 + 3 v^2
            unit = (u + v, 2 * v)         # N = (u+v)^2 - (u+v)2v + 4v^2 = u^2+3v^2
            if eis_norm(unit) != p:
                return None
            reps = e
        else:                             # p = 2 (mod 3): needs an even exponent
            if e % 2 == 1:
                return None
            unit = (p, 0)                 # N(p) = p^2
            reps = e // 2
        for _ in range(reps):
            z = eis_mul(z, unit)

    if eis_norm(z) != t:
        return None
    a, b = z
    return (a, -b)                        # a^2 - a*b + b^2 = x^2 + x*y + y^2


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__)
        return 2
    a, d, n = (int(v) for v in sys.argv[1:4])
    print(f"Independent constructive cross-check: a={a} d={d} n={n}")
    print("Each term is certified by an explicit (x, y) with x^2+x*y+y^2 = t.\n")
    ok = True
    for k in range(n):
        t = a + k * d
        rep = represent(t)
        if rep is None:
            print(f"  k={k:3d}  t={t}   FAIL: no representation exists")
            ok = False
            continue
        x, y = rep
        lhs = x * x + x * y + y * y
        if lhs != t:
            print(f"  k={k:3d}  t={t}   FAIL: {x}^2+{x}*{y}+{y}^2 = {lhs} != {t}")
            ok = False
            continue
        print(f"  k={k:3d}  t={t} = ({x})^2 + ({x})*({y}) + ({y})^2   OK")
    print()
    print(f"CROSS-CHECK: {'PASS' if ok else 'FAIL'}  "
          f"({n} terms, each with a verified x^2+x*y+y^2 representation)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

"""Loeschian numbers: n = x^2 + x*y + y^2 for some x, y in Z.

Core fact (GOAL.md section 2): n >= 1 is Loeschian iff every prime p = 2 (mod 3)
divides n to an even power. Primes 3 and p = 1 (mod 3) are unconstrained.
n = 0 is Loeschian (x = y = 0).

Everything here uses Python ints. No fixed-width arithmetic anywhere, so the
int64 overflow trap in GOAL.md section 2b cannot bite.
"""

from __future__ import annotations

import random

# --- primality / factorization (stdlib only, exact for arbitrary size) -------

_SMALL_PRIMES = None


def small_primes(limit: int = 100_000) -> list[int]:
    """Primes < limit, cached for the default limit."""
    global _SMALL_PRIMES
    if limit == 100_000 and _SMALL_PRIMES is not None:
        return _SMALL_PRIMES
    sieve = bytearray([1]) * limit
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(range(i * i, limit, i)))
    out = [i for i in range(limit) if sieve[i]]
    if limit == 100_000:
        _SMALL_PRIMES = out
    return out


def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin for n < 3.3e24, strong-probable beyond."""
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _pollard_brent(n: int, rng: random.Random) -> int:
    """Return a nontrivial factor of composite, odd, non-prime n."""
    if n % 2 == 0:
        return 2
    while True:
        y = rng.randrange(1, n)
        c = rng.randrange(1, n)
        m = 128
        g = r = q = 1
        x = ys = y
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = _gcd(q, n)
                k += m
            r *= 2
        if g == n:
            g = 1
            y = ys
            while g == 1:
                y = (y * y + c) % n
                g = _gcd(abs(x - y), n)
        if g != n:
            return g


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


def factorize(n: int) -> dict[int, int]:
    """Full prime factorization of n >= 1 as {prime: exponent}."""
    if n < 1:
        raise ValueError(f"factorize needs n >= 1, got {n}")
    out: dict[int, int] = {}
    if n == 1:
        return out
    for p in small_primes():
        if p * p > n:
            break
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    if n == 1:
        return out
    rng = random.Random(0xC0FFEE ^ (n & 0xFFFFFFFF))
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            out[m] = out.get(m, 0) + 1
            continue
        r = _isqrt_exact(m)
        if r is not None:
            # m is a perfect square: factor the root, double the exponents.
            for p, e in factorize(r).items():
                out[p] = out.get(p, 0) + 2 * e
            continue
        f = _pollard_brent(m, rng)
        stack.append(f)
        stack.append(m // f)
    return out


def _isqrt_exact(n: int) -> int | None:
    import math

    r = math.isqrt(n)
    return r if r * r == n else None


# --- the Loeschian test ------------------------------------------------------


def is_bad_prime(p: int) -> bool:
    """A 'bad' prime is p = 2 (mod 3): it must occur to an even power."""
    return p % 3 == 2


def is_loeschian(n: int) -> bool:
    """True iff n is of the form x^2 + x*y + y^2 with x, y in Z.

    Exact for all n >= 0, any size.
    """
    if n < 0:
        return False
    if n == 0:
        return True
    for p, e in factorize(n).items():
        if p % 3 == 2 and e % 2 == 1:
            return False
    return True


def is_loeschian_bruteforce(n: int) -> bool:
    """Independent cross-check: search x, y >= 0 with x^2 + x*y + y^2 == n.

    Non-negative x, y suffice: the form's value set over Z equals its value set
    over the non-negative cone (A003136). Only for small n -- it is O(sqrt(n)).
    """
    if n < 0:
        return False
    x = 0
    while x * x <= n:
        # solve x^2 + x*y + y^2 = n for y >= 0
        disc = x * x - 4 * (x * x - n)
        if disc >= 0:
            import math

            r = math.isqrt(disc)
            for rr in (r, r + 1):
                if rr * rr == disc and (-x + rr) % 2 == 0 and (-x + rr) >= 0:
                    return True
        x += 1
    return False


def loeschian_sieve(limit: int) -> bytearray:
    """Bitmap-free sieve: out[n] == 1 iff n is Loeschian, for 0 <= n <= limit.

    Built by direct enumeration of x^2 + x*y + y^2 <= limit with x, y >= 0,
    which is O(limit) pairs. Pure Python -- use the C searcher for big limits.
    """
    if limit < 0:
        raise ValueError("limit must be >= 0")
    out = bytearray(limit + 1)
    x = 0
    while x * x <= limit:
        y = 0
        while True:
            v = x * x + x * y + y * y
            if v > limit:
                break
            out[v] = 1
            y += 1
        x += 1
    return out


# --- arithmetic progressions -------------------------------------------------

#: Bad primes p = 2 (mod 3) in increasing order, as far as asked for.
def bad_primes(limit: int) -> list[int]:
    return [p for p in small_primes() if p < limit and p % 3 == 2]


def forced_bad_primes(n: int) -> list[int]:
    """Bad primes that MUST divide d for an n-term AP to exist (GOAL.md 2b).

    If a bad prime p does not divide d, the terms divisible by p form a
    sub-progression with step p*d. If two or more terms are divisible by p,
    both need v_p >= 2, i.e. p^2 | t and p^2 | t', so p^2 | (t' - t) = p*j*d
    with j < p, forcing p | d -- contradiction. Two or more hits are guaranteed
    exactly when the window of n consecutive indices covers some residue class
    mod p twice, i.e. when n >= 2*p (n > p + (p-1), i.e. n >= 2p).
    """
    return [p for p in bad_primes(n) if 2 * p <= n]


def ap_run_length(a: int, d: int, cap: int = 10_000) -> int:
    """Number of consecutive Loeschian terms starting at a with step d."""
    if d < 1:
        raise ValueError("d must be >= 1")
    k = 0
    while k < cap and is_loeschian(a + k * d):
        k += 1
    return k

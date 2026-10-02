"""Unit tests for the Loeschian core. Encodes the GOAL.md section 2b invariants.

Run with:  python3 -m unittest discover -s tests -v    (or: make check)
"""

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from loeschian import (  # noqa: E402
    ap_run_length,
    bad_primes,
    factorize,
    forced_bad_primes,
    is_loeschian,
    is_loeschian_bruteforce,
    is_prime,
    loeschian_sieve,
)

# A003136, from OEIS / GOAL.md section 1.
KNOWN_PREFIX = [0, 1, 3, 4, 7, 9, 12, 13, 16, 19, 21, 25, 27, 28, 31, 36, 37,
                39, 43, 48, 49, 52, 57, 61, 63, 64, 67, 73, 75, 76, 79, 81,
                84, 91, 93, 97, 100]


class TestBadPrimes(unittest.TestCase):
    def test_modulus_is_3_not_4(self):
        # The classic muscle-memory bug is p = 3 (mod 4). Pin the real list.
        self.assertEqual(
            bad_primes(110),
            [2, 5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83, 89, 101, 107],
        )

    def test_bad_primes_are_not_loeschian_but_squares_are(self):
        for p in bad_primes(200):
            self.assertFalse(is_loeschian(p), f"{p} must not be Loeschian")
            self.assertTrue(is_loeschian(p * p), f"{p}^2 must be Loeschian")
            self.assertFalse(is_loeschian(p**3))
            self.assertTrue(is_loeschian(p**4))

    def test_good_primes_are_loeschian(self):
        for p in (3, 7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103, 109):
            self.assertTrue(is_loeschian(p), f"{p} should be Loeschian")

    def test_forced_bad_primes_for_target_lengths(self):
        # GOAL.md: for n = 58 these must divide d.
        self.assertEqual(forced_bad_primes(58), [2, 5, 11, 17, 23, 29])
        self.assertEqual(forced_bad_primes(4), [2])
        self.assertEqual(forced_bad_primes(35), [2, 5, 11, 17])
        self.assertEqual(forced_bad_primes(120), [2, 5, 11, 17, 23, 29, 41, 47, 53, 59])


class TestEdgeValues(unittest.TestCase):
    def test_zero_and_one_are_loeschian(self):
        self.assertTrue(is_loeschian(0))
        self.assertTrue(is_loeschian(1))

    def test_negatives_are_not(self):
        self.assertFalse(is_loeschian(-1))
        self.assertFalse(is_loeschian(-3))


class TestAgainstKnownPrefix(unittest.TestCase):
    def test_prefix_matches(self):
        got = [n for n in range(101) if is_loeschian(n)]
        self.assertEqual(got, KNOWN_PREFIX)

    def test_sieve_matches_prefix(self):
        s = loeschian_sieve(100)
        self.assertEqual([n for n in range(101) if s[n]], KNOWN_PREFIX)


class TestCrossImplementation(unittest.TestCase):
    """Three independent implementations must agree: factor test, brute-force
    x^2+xy+y^2 search, and the enumeration sieve."""

    def test_all_three_agree_to_20000(self):
        limit = 20_000
        sieve = loeschian_sieve(limit)
        for n in range(limit + 1):
            a = bool(sieve[n])
            b = is_loeschian(n)
            self.assertEqual(a, b, f"sieve vs factor-test disagree at n={n}")

    def test_bruteforce_agrees_on_sample(self):
        rng = random.Random(12345)
        for _ in range(400):
            n = rng.randrange(0, 50_000)
            self.assertEqual(
                is_loeschian(n), is_loeschian_bruteforce(n),
                f"factor-test vs brute-force disagree at n={n}")

    def test_bruteforce_agrees_on_structured_sample(self):
        # numbers built to stress the even-power rule
        for n in (4, 25, 100, 121, 4 * 7, 2 * 7, 2 * 2 * 5 * 5, 5 * 5 * 11,
                  5 * 11, 2**6 * 3, 2**5 * 3, 11**2 * 17**2, 11**2 * 17):
            self.assertEqual(is_loeschian(n), is_loeschian_bruteforce(n), f"n={n}")


class TestFactorize(unittest.TestCase):
    def test_small(self):
        self.assertEqual(factorize(1), {})
        self.assertEqual(factorize(2), {2: 1})
        self.assertEqual(factorize(360), {2: 3, 3: 2, 5: 1})

    def test_big_semiprime_and_square(self):
        p, q = 1_000_000_007, 1_000_000_009
        self.assertTrue(is_prime(p) and is_prime(q))
        self.assertEqual(factorize(p * q), {p: 1, q: 1})
        self.assertEqual(factorize(p * p), {p: 2})
        # products always re-multiply to the input
        rng = random.Random(7)
        for _ in range(25):
            n = rng.randrange(10**14, 10**16)
            prod = 1
            for pp, e in factorize(n).items():
                self.assertTrue(is_prime(pp))
                prod *= pp**e
            self.assertEqual(prod, n)

    def test_beyond_2_63(self):
        # GOAL.md 2b: must be exact past the int64 cliff.
        n = (2**64 + 1) * 3
        prod = 1
        for p, e in factorize(n).items():
            prod *= p**e
        self.assertEqual(prod, n)
        self.assertTrue(is_loeschian(3 * 7 * 13 * 2**64))  # 2^64 is an even power


class TestProgressions(unittest.TestCase):
    def test_puzzle_example(self):
        # 1, 7, 13, 19: n=4, a=1, d=6
        for t in (1, 7, 13, 19):
            self.assertTrue(is_loeschian(t))
        self.assertGreaterEqual(ap_run_length(1, 6), 4)
        # The puzzle's 4-term example is not maximal: it runs
        # 1,7,13,19,25,31,37,43,49 and stops at 55 = 5*11 (two bad primes to
        # the first power). Pin the true length so an off-by-one shows up.
        self.assertEqual(ap_run_length(1, 6), 9)
        self.assertFalse(is_loeschian(1 + 9 * 6))  # 55 = 5 * 11

    def test_run_length_counts_terms_not_steps(self):
        # n terms means last index n-1 (GOAL.md 2b off-by-one guard)
        self.assertEqual(ap_run_length(0, 1), 2)  # 0, 1 Loeschian; 2 is not
        self.assertEqual(ap_run_length(3, 1), 2)  # 3, 4 Loeschian; 5 is not

    def test_d_must_be_positive(self):
        with self.assertRaises(ValueError):
            ap_run_length(1, 0)

    def test_square_scaling_preserves_aps(self):
        # (a, d) -> (m^2 a, m^2 d) keeps every term Loeschian (GOAL.md 2-3)
        a, d, n = 1, 6, 5
        for m in (2, 3, 5, 7, 11):
            for k in range(n):
                self.assertTrue(is_loeschian(m * m * (a + k * d)))


class TestSieveBounds(unittest.TestCase):
    def test_sieve_is_dense_enough_to_be_believable(self):
        s = loeschian_sieve(100_000)
        count = sum(s)
        # density ~ K/sqrt(log N), K ~ 0.638 -> ~0.19 at 1e5; keep a loose band
        self.assertTrue(0.15 < count / 100_000 < 0.30, count / 100_000)



class TestConstructiveCrossCheck(unittest.TestCase):
    """src/crosscheck.py must agree with the bad-prime criterion, and every
    representation it returns must satisfy x^2+x*y+y^2 = t exactly. This is the
    independent check GOAL.md §8 requires before a record is claimed."""

    def test_agrees_with_criterion_and_arithmetic(self):
        import crosscheck
        for t in range(0, 3000):
            rep = crosscheck.represent(t)
            self.assertEqual(rep is not None, is_loeschian(t),
                             f"representability disagrees with criterion at t={t}")
            if rep is not None:
                x, y = rep
                self.assertEqual(x * x + x * y + y * y, t, f"bad rep for t={t}")

    def test_large_values(self):
        import crosscheck
        import random
        rng = random.Random(4242)
        checked = 0
        for _ in range(200):
            t = rng.randrange(10**14, 10**15)
            rep = crosscheck.represent(t)
            self.assertEqual(rep is not None, is_loeschian(t), f"t={t}")
            if rep is not None:
                x, y = rep
                self.assertEqual(x * x + x * y + y * y, t, f"t={t}")
                checked += 1
        self.assertGreater(checked, 10, "too few large Loeschian samples")

if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Authoritative verifier for Loeschian arithmetic progressions.

Nothing becomes a record until this script says PASS (GOAL.md section 3).
It factors every term from scratch with Python ints -- it never consults a
sieve, a searcher's verdict, or any fixed-width arithmetic.

Usage:
    verify.py A D N            # check the N-term AP a=A, d=D
    verify.py A D N --quiet    # only the summary lines
    verify.py A D N --maximal  # also report whether the run extends

Exit status 0 on PASS, 1 on FAIL.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from loeschian import factorize, is_loeschian  # noqa: E402


def fmt_factors(f: dict[int, int]) -> str:
    return " * ".join(
        f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items())
    ) or "1"


def verify(a: int, d: int, n: int, quiet: bool = False, maximal: bool = False) -> bool:
    ok = True
    if d < 1:
        print(f"FAIL: d must be >= 1, got d={d}  (a d=0 'progression' is degenerate)")
        return False
    if a < 0:
        print(f"FAIL: a must be >= 0, got a={a}")
        return False
    if n < 1:
        print(f"FAIL: n must be >= 1, got n={n}")
        return False

    last = a + (n - 1) * d
    print(f"AP: a = {a}")
    print(f"    d = {d}")
    print(f"    n = {n} terms  (indices k = 0 .. {n - 1})")
    print(f"    last term a + {n - 1}*d = {last}")
    print(f"    last term has {len(str(last))} decimal digits")
    print()

    for k in range(n):
        t = a + k * d
        if t == 0:
            if not quiet:
                print(f"  k={k:3d}  t=0  = 0^2+0*0+0^2            PASS (0 is Loeschian)")
            continue
        f = factorize(t)
        # independent re-multiplication guards against a factorization bug
        prod = 1
        for p, e in f.items():
            prod *= p**e
        if prod != t:
            print(f"  k={k:3d}  t={t}  FACTORIZATION BUG: product {prod} != {t}")
            return False
        bad_odd = sorted(p for p, e in f.items() if p % 3 == 2 and e % 2 == 1)
        good = not bad_odd
        ok &= good
        if not quiet or not good:
            verdict = "PASS" if good else f"FAIL (odd power of bad prime(s) {bad_odd})"
            print(f"  k={k:3d}  t={t}  = {fmt_factors(f)}   {verdict}")

    print()
    if maximal and ok:
        before = a - d
        ext_lo = before >= 0 and is_loeschian(before)
        ext_hi = is_loeschian(a + n * d)
        print(f"  maximality: a-d = {before} "
              f"{'IS' if ext_lo else 'is NOT'} Loeschian"
              f"{' (negative, n/a)' if before < 0 else ''}")
        print(f"  maximality: a+{n}*d = {a + n * d} "
              f"{'IS' if ext_hi else 'is NOT'} Loeschian")
        if ext_lo or ext_hi:
            print("  NOTE: this run EXTENDS -- reported n understates the record.")
        else:
            print("  run is maximal in both directions")
        print()

    print(f"OVERALL: {'PASS' if ok else 'FAIL'}  (n={n}, a={a}, d={d})")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", type=int)
    ap.add_argument("d", type=int)
    ap.add_argument("n", type=int)
    ap.add_argument("--quiet", action="store_true", help="suppress per-term lines")
    ap.add_argument("--maximal", action="store_true",
                    help="also test whether the run extends either way")
    args = ap.parse_args()
    return 0 if verify(args.a, args.d, args.n, args.quiet, args.maximal) else 1


if __name__ == "__main__":
    raise SystemExit(main())

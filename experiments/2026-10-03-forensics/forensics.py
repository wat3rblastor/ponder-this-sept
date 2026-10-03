#!/usr/bin/env python3
"""Forensics on the long Loeschian APs we possess, vs a matched control.

Control (per record): random starts a' with a' = a (mod M), where
M = prod_{bad p | d} p^(v_p(d)+2).  Congruence mod M fixes exactly the local
conditions at the bad primes dividing d (which is what makes a start
admissible at all), so the control terms t'_k = a'+k*d satisfy the same local
constraints as the real terms, are the same size, and lie in the same
arithmetic class mod d.  We keep only those t'_k that are Loeschian: the null
hypothesis is "the terms of a long AP look like independent random Loeschian
numbers of that size in that class".
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from loeschian import factorize, is_loeschian  # noqa: E402

RECORDS = [
    ("55a", 11687581876345393, 202927451159070, 55),
    ("55b", 296246969176050787, 11950172123811900, 55),
    ("47", 2646171143023357, 78342989618850, 47),
    ("43", 5413537078288507, 48916598396160, 43),
    ("35opt", 219830911, 2709630, 35),
]

SMALL_C = list(range(1, 61))


def issq(n: int) -> bool:
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n


def bad_part_modulus(d: int) -> int:
    M = 1
    for p, e in factorize(d).items():
        if p % 3 == 2:
            M *= p ** (e + 2)
    return M


def term_features(t: int, d: int) -> dict:
    f = factorize(t)
    lp = max(f)
    # square part contributed by GOOD primes (bad primes are forced even)
    good_sq = 1
    good_distinct = 0
    bad_distinct = 0
    for p, e in f.items():
        if p % 3 == 2:
            bad_distinct += 1
        else:
            good_distinct += 1
            if e >= 2:
                good_sq *= p ** (e - e % 2)
    rad = 1
    for p in f:
        rad *= p
    feats = {
        "lp": lp,
        "lp_ratio": math.log(lp) / math.log(t),
        "nfac": sum(f.values()),
        "ndistinct": len(f),
        "good_sq": good_sq,
        "good_sq_big": int(good_sq > 1),
        "rad_ratio": math.log(rad) / math.log(t),
        "smooth1e4": int(lp <= 10**4),
        "smooth1e6": int(lp <= 10**6),
        "bad_distinct": bad_distinct,
    }
    # t = c * perfect square, for small c
    cs = [c for c in SMALL_C if t % c == 0 and issq(t // c)]
    feats["csquare"] = cs[0] if cs else 0
    feats["is_csquare"] = int(bool(cs))
    # automatic family: 3*c*d*t is a perfect square
    auto = [c for c in SMALL_C if issq(3 * c * d * t)]
    feats["auto"] = auto[0] if auto else 0
    feats["is_auto"] = int(bool(auto))
    # weaker variant: c*d*t a perfect square
    feats["is_auto2"] = int(any(issq(c * d * t) for c in SMALL_C))
    return feats


def collect_real(name, a, d, n):
    out = []
    for k in range(n):
        t = a + k * d
        assert is_loeschian(t), (name, k)
        fe = term_features(t, d)
        fe["k"] = k
        fe["rec"] = name
        out.append(fe)
    return out


def collect_control(name, a, d, n, target, rng):
    M = bad_part_modulus(d)
    out = []
    tries = 0
    while len(out) < target and tries < 200000:
        tries += 1
        # a' = a + M*j, same magnitude as a
        j = rng.randrange(-(a // (2 * M)), a // (2 * M) + 1)
        ap = a + M * j
        if ap <= 0:
            continue
        k = rng.randrange(n)
        t = ap + k * d
        if not is_loeschian(t):
            continue
        fe = term_features(t, d)
        fe["k"] = k
        fe["rec"] = name
        out.append(fe)
    return out


def main():
    target = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    rng = random.Random(12345)
    real, ctrl = [], []
    for name, a, d, n in RECORDS:
        r = collect_real(name, a, d, n)
        c = collect_control(name, a, d, n, target, rng)
        real += r
        ctrl += c
        print(f"{name}: real {len(r)} ctrl {len(c)}", flush=True)
    Path("real_terms.json").write_text(json.dumps(real))
    Path("ctrl_terms.json").write_text(json.dumps(ctrl))


if __name__ == "__main__":
    main()

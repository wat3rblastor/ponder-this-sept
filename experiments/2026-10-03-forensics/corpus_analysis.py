#!/usr/bin/env python3
"""Analysis of the exhaustive corpus from build/loeschsearch.

Part 1: which d are good?  Regress log(#starts with run >= k) on log d and
        indicators for each bad prime dividing d.  Exhaustive counts, so the
        only noise is the arithmetic itself.
Part 2: per-term forensics on the (exhaustively optimal) longest run for each
        d, against the same matched control as the records.
Part 3: hit-index statistics against the CORRECT conditional null.
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from loeschian import bad_primes, factorize, is_loeschian  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
from forensics import bad_part_modulus, term_features  # noqa: E402

BADS = [2, 5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83, 89]


def load(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        o = json.loads(line)
        if "g1" in o:
            continue
        rows.append(o)
    return rows


def lstsq(X, y):
    """Ordinary least squares via normal equations (small p)."""
    p = len(X[0])
    A = [[sum(X[i][a] * X[i][b] for i in range(len(X))) for b in range(p)]
         for a in range(p)]
    b = [sum(X[i][a] * y[i] for i in range(len(X))) for a in range(p)]
    # gaussian elimination with partial pivot
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(p):
        piv = max(range(c, p), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        if abs(M[c][c]) < 1e-12:
            continue
        for r in range(p):
            if r == c:
                continue
            f = M[r][c] / M[c][c]
            for cc in range(c, p + 1):
                M[r][cc] -= f * M[c][cc]
    beta = [M[i][p] / M[i][i] if abs(M[i][i]) > 1e-12 else 0.0 for i in range(p)]
    resid = [y[i] - sum(X[i][j] * beta[j] for j in range(p)) for i in range(len(X))]
    rss = sum(r * r for r in resid)
    dof = max(1, len(X) - p)
    s2 = rss / dof
    # standard errors from (X'X)^-1 diag -- recompute inverse diag crudely
    # invert A
    n = p
    aug = [A[i][:] + [1.0 if j == i else 0.0 for j in range(n)] for i in range(n)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(aug[r][c]))
        aug[c], aug[piv] = aug[piv], aug[c]
        if abs(aug[c][c]) < 1e-12:
            continue
        pv = aug[c][c]
        aug[c] = [v / pv for v in aug[c]]
        for r in range(n):
            if r == c:
                continue
            f = aug[r][c]
            if f:
                aug[r] = [aug[r][k] - f * aug[c][k] for k in range(2 * n)]
    se = [math.sqrt(max(0.0, s2 * aug[i][n + i])) for i in range(n)]
    return beta, se, math.sqrt(s2)


def part1(rows, label, kset):
    print(f"\n===== PART 1 ({label}): which d give many long runs?  "
          f"{len(rows)} exhaustively-scanned steps d")
    for k in kset:
        data = [(r["d"], r["ge"].get(str(k), 0)) for r in rows
                if str(k) in r["ge"] or True]
        data = [(d, c) for d, c in data if c > 0]
        if len(data) < 20:
            print(f"  k={k}: only {len(data)} steps with any run >= k -- skip")
            continue
        used = [p for p in BADS
                if 3 <= sum(1 for d, _ in data if d % p == 0) <= len(data) - 3]
        X, y = [], []
        for d, c in data:
            row = [1.0, math.log(d)]
            for p in used:
                row.append(1.0 if d % p == 0 else 0.0)
            X.append(row)
            y.append(math.log(c))
        beta, se, s = lstsq(X, y)
        names = ["const", "log d"] + [f"p={p}" for p in used]
        print(f"  --- k={k}: n={len(data)} steps, resid sd={s:.2f} (log scale)")
        for nm, b, e in zip(names, beta, se):
            star = "  ***" if abs(b) > 2.5 * e and nm not in ("const",) else ""
            mult = f"  x{math.exp(b):.2f}" if nm.startswith("p=") else ""
            print(f"      {nm:<8} {b:+8.3f} +- {e:.3f}{mult}{star}")
        # raw group means for the headline primes
        for p in used:
            w = [c for d, c in data if d % p == 0]
            wo = [c for d, c in data if d % p != 0]
            if w and wo:
                gm = math.exp(sum(math.log(c) for c in w) / len(w))
                gmo = math.exp(sum(math.log(c) for c in wo) / len(wo))
                print(f"      raw geo-mean count p={p}: with {gm:9.1f} (n={len(w)})"
                      f"  without {gmo:9.1f} (n={len(wo)})  ratio {gm/gmo:7.2f}x")


def part2(rows, nmin, target_ctrl, rng):
    sel = [r for r in rows if r["best_n"] >= nmin and r["best_a"] > 0]
    print(f"\n===== PART 2: per-term forensics on {len(sel)} exhaustively-optimal "
          f"runs with n >= {nmin}")
    real, ctrl = [], []
    for r in sel:
        a, d, n = r["best_a"], r["d"], r["best_n"]
        for k in range(n):
            t = a + k * d
            assert is_loeschian(t), (a, d, k)
            fe = term_features(t, d)
            fe["k"] = k
            fe["n"] = n
            real.append(fe)
        M = bad_part_modulus(d)
        got = 0
        tries = 0
        while got < target_ctrl and tries < 100000:
            tries += 1
            j = rng.randrange(-(a // (2 * M)) if a > 2 * M else -1, a // (2 * M) + 1)
            ap = a + M * j
            if ap <= 0:
                continue
            k = rng.randrange(n)
            t = ap + k * d
            if not is_loeschian(t):
                continue
            fe = term_features(t, d)
            fe["k"] = k
            fe["n"] = n
            ctrl.append(fe)
            got += 1
    Path("corpus_real_terms.json").write_text(json.dumps(real))
    Path("corpus_ctrl_terms.json").write_text(json.dumps(ctrl))
    print(f"  real terms {len(real)}  control terms {len(ctrl)}")
    return real, ctrl


def part3(rows, nmin):
    sel = [r for r in rows if r["best_n"] >= nmin and r["best_a"] > 0]
    print(f"\n===== PART 3: hit indices for bad q not dividing d, "
          f"{len(sel)} optimal runs, q < 1000")
    obs = 0
    exp_uncond = 0.0
    exp_cond = 0.0
    kpos = []
    for r in sel:
        a, d, n = r["best_a"], r["d"], r["best_n"]
        for q in bad_primes(1000):
            if d % q == 0:
                continue
            k0 = (-a * pow(d, -1, q)) % q
            hits = [k for k in range(n) if (k - k0) % q == 0]
            nh = len(hits)
            # unconditional expected number of primes with >=1 hit
            exp_uncond += min(1.0, n / q)
            # conditional on the AP being all-Loeschian.  Over k0 uniform mod q:
            #   hits==0 -> weight 1; hits==1 -> weight 1/q (needs q^2 | t);
            #   hits>=2 -> weight 0 (would force q | d).
            A = sum(1 for z in range(q) if not any((k - z) % q == 0
                                                   for k in range(n)))
            B = sum(1 for z in range(q)
                    if sum(1 for k in range(n) if (k - z) % q == 0) == 1)
            exp_cond += (B / q) / (A + B / q) if (A + B / q) > 0 else 1.0
            if nh:
                obs += 1
                kpos.append((q, k0, n, hits))
    print(f"  primes with >=1 hit: obs={obs}  exp_uncond={exp_uncond:.1f}  "
          f"exp_cond(all-Loeschian)={exp_cond:.2f}")
    if exp_cond:
        print(f"  obs/exp_cond = {obs/exp_cond:.2f}x")
    for q, k0, n, hits in kpos[:40]:
        print(f"    q={q} n={n} k0={k0} hits={hits} "
              f"({'END' if min(hits) < 3 or max(hits) > n - 4 else 'mid'})")
    # edge clustering of hit index
    if kpos:
        rel = [h / (n - 1) for q, k0, n, hits in kpos for h in hits]
        print(f"  relative hit position mean={sum(rel)/len(rel):.3f} "
              f"(uniform null 0.5), n={len(rel)}")


if __name__ == "__main__":
    rows = load(sys.argv[1])
    kset = [int(x) for x in sys.argv[2].split(",")]
    nmin = int(sys.argv[3])
    part1(rows, Path(sys.argv[1]).name, kset)
    rng = random.Random(777)
    part2(rows, nmin, 60, rng)
    part3(rows, nmin)

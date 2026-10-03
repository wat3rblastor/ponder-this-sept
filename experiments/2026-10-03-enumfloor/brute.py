#!/usr/bin/env python3
"""Independent brute-force check of enum.c: scan every a < X and compare sets."""
import subprocess, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
L = 58

def brute(d, X, Q):
    """all a in [0,X) with a = j*d (mod q), 1<=j<=q-58, for every q in Q"""
    ok = []
    inv = {q: pow(d % q, -1, q) for q in Q}
    for a in range(X):
        good = True
        for q in Q:
            j = (a % q) * inv[q] % q
            if not (1 <= j <= q - L): good = False; break
        if good: ok.append(a)
    return ok

def run(d, X, Q, split=None):
    cmd = [os.path.join(HERE, "enum"), "--d", str(d), "--X", str(X),
           "--Q", ",".join(map(str, Q)), "--emit", "--verify"]
    if split is not None: cmd += ["--split", str(split)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    got = sorted(int(l.split()[1]) for l in r.stdout.splitlines() if l.startswith("A "))
    stat = [l for l in r.stdout.splitlines() if l.startswith("out=")][0]
    return got, stat

D0 = 382160924970
cases = [
    (D0, 2_000_000, [59, 71, 83]),
    (D0, 2_000_000, [59, 71, 83, 89]),
    (D0, 20_000_000, [59, 71, 83, 89, 101]),
    (D0, 20_000_000, [59, 71, 83, 89, 101, 107, 113]),
    (D0 * 7, 20_000_000, [59, 71, 83, 89, 101, 107, 113]),
    (D0 * 1234567, 5_000_000, [59, 71, 83, 89, 101, 107, 113, 131, 137]),
]
fail = 0
for d, X, Q in cases:
    exp = brute(d, X, Q)
    for split in (None, 0, 1, 2, 3):
        if split is not None and split > len(Q): continue
        got, stat = run(d, X, Q, split)
        same = got == exp
        fail += not same
        print("d=%-16d X=%-9d Q=%-28s split=%-4s brute=%-7d enum=%-7d %s"
              % (d, X, ",".join(map(str, Q)), split, len(exp), len(got),
                 "MATCH" if same else "*** MISMATCH ***"))
    print("    " + stat)
print("FAILURES:", fail)
sys.exit(1 if fail else 0)

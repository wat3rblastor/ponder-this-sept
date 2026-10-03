#!/usr/bin/env python3
"""Measured pass rate of F17 on the 41 non-automatic indices, by term size,
and a zero-failure check on the 17 automatic indices.

Verdicts come from build/apsearch --isl (exact, 64-bit), so T <= ~9e18 only.
A sample of the automatic terms is independently re-verified by src/verify.py
and src/crosscheck.py in verify_auto.txt.
"""
import random, subprocess, sys
from math import isqrt
sys.path.insert(0, '.')
from family import A_FORM, D_FORM, AUTO, M, val

NON = [k for k in range(58) if k not in AUTO]
PR = [2, 5, 11, 17, 23, 29]
ROOT = '../..'


def Mprime(p):
    r = 1
    for l in PR:
        if p % l: r *= l
    return r


def sample(Tmax, rng, nm):
    out = []
    while len(out) < nm:
        p = rng.randrange(1, isqrt(Tmax//21)+1)
        rem = Tmax - 21*p*p
        if rem <= 0: continue
        U = isqrt(rem//7)
        if U <= 6*p: continue
        m = Mprime(p); c = (U - 6*p)//m
        if c < 1: continue
        u = 6*p + m*(rng.randrange(c)+1); q = u - 12*p
        a = val(A_FORM, p, q); d = val(D_FORM, p, q)
        assert d > 0 and d % M == 0
        T = a + 57*d
        if T > Tmax or T < Tmax//10: continue
        out.append((a, d))
    return out


if __name__ == "__main__":
    for e in [10, 12, 14, 16, 18]:
        rng = random.Random(777+e)
        nm = 4000 if e < 18 else 2500
        mem = sample(10**e, rng, nm)
        lines = [a+k*d for a, d in mem for k in range(58)]
        r = subprocess.run([ROOT+'/build/apsearch', '--isl'],
                           input="\n".join(map(str, lines)),
                           capture_output=True, text=True)
        v = [int(x.split()[1]) for x in r.stdout.strip().split("\n")]
        assert len(v) == len(lines)
        af = np_ = nt = 0; runs = []
        for i in range(len(mem)):
            w = v[i*58:(i+1)*58]
            af += sum(1 for k in AUTO if not w[k])
            for k in NON:
                np_ += w[k]; nt += 1
            best = cur = 0
            for x in w:
                cur = cur+1 if x else 0; best = max(best, cur)
            runs.append(best)
        rho = np_/nt
        print(f"T~1e{e}: members {len(mem)}  automatic failures {af}/{len(mem)*17}"
              f"  rho_nonauto={rho:.4f}  rho^41={rho**41:.3e}"
              f"  longest run {max(runs)}  mean {np_/len(mem):.2f}/41")
        sys.stdout.flush()

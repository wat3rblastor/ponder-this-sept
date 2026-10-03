"""Independent cross-check of the new designed families: brute-force the SAME
count with src/loeschian.is_loeschian (pure Python, exact factorisation) and
compare to exact.c.  Validates the family conventions (which d, which a are
admissible), not just the Loeschian bitmap."""
import os, subprocess, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
from loeschian import is_loeschian

HERE = os.path.dirname(os.path.abspath(__file__))
BAD = [2, 5, 11, 17, 23, 29]

CASES = [(748374, 9, 10_000_000), (162690, 24, 6_000_000), (340170, 15, 8_000_000)]
for B, n, X in CASES:
    qs = [q for q in BAD if B % q == 0]
    tot = rs = 0
    m = 1
    while (n - 1) * B * m < X:
        d = B * m
        for a in range(1, X - (n - 1) * d + 1):
            if a % 3 != 1:
                continue
            if any(a % q == 0 for q in qs):
                continue
            if all(is_loeschian(a + k * d) for k in range(n)):
                tot += 1
                if a - d < 1 or not is_loeschian(a - d):
                    rs += 1
        m += 1
    out = subprocess.run([os.path.join(HERE, 'exact'), str(X), str(n), str(n),
                          '1', '--base', str(B)], capture_output=True, text=True,
                         env=dict(os.environ, OMP_NUM_THREADS='2'))
    line = [l for l in out.stdout.splitlines() if l.startswith('n=')][0]
    f = dict(kv.split('=', 1) for kv in line.split() if '=' in kv)
    ok = (int(f['COUNT']) == tot and int(f['RUNSTART']) == rs)
    print(f"base={B} n={n} X={X}: python COUNT={tot} RUNSTART={rs} | "
          f"exact.c COUNT={f['COUNT']} RUNSTART={f['RUNSTART']}  "
          f"{'AGREE' if ok else '*** DISAGREE ***'}")

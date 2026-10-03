"""The raw exact counts of this experiment, as a standalone table.

CONVENTIONS (stated explicitly; this repo has been burned by them before):
  COUNT    = number of PAIRS (a,d) with d = B*m (m >= 1), a >= 1, a = 1 mod 3,
             a coprime to every bad prime dividing B, a+(n-1)d <= X, and
             a+kd Loeschian for all 0 <= k < n.
             So COUNT counts APs of length AT LEAST n, once per starting pair:
             an AP of length L >= n contributes L-n+1.  ALL APs, not maximal.
  RUNSTART = the same, restricted to pairs with a-d not Loeschian (or < 1):
             runs of length >= n, each counted once.
  MAXIMAL  = runs of length >= n not extendable in either direction in [1,X].
These are exact integers, not estimates.  Produced by exact.c (the counter of
experiments/2026-10-04-groundtruth, reused unmodified except that --dens now
accepts --base).
"""
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))
QMIN = {748374: 5, 340170: 11, 220110: 17, 162690: 23, 129030: 29, 3741870: 41}
FACT = {748374: '2*3*11*17*23*29', 340170: '2*3*5*17*23*29',
        220110: '2*3*5*11*23*29', 162690: '2*3*5*11*17*29',
        129030: '2*3*5*11*17*23', 3741870: '2*3*5*11*17*23*29 = base(58)'}

print(__doc__)
rows = []
for f in sorted(glob.glob(os.path.join(HERE, 'q*.txt'))):
    X = None
    for line in open(f):
        if line.startswith('# X='):
            X = int(line.split('=')[1])
        if line.startswith('n='):
            d = dict(kv.split('=', 1) for kv in line.split() if '=' in kv)
            if 'COUNT' in d:
                rows.append((int(d['base']), X, int(d['n']), int(d['nd']),
                             int(d['COUNT']), int(d['RUNSTART']),
                             int(d['MAXIMAL']), float(d['t'][:-1])))
rows.sort(key=lambda r: (QMIN.get(r[0], 0), r[1], r[2]))
cur = None
for base, X, n, nd, c, rs, mx, t in rows:
    if (base, X) != cur:
        cur = (base, X)
        q = QMIN.get(base)
        print(f"\nFAMILY d = {base}*m   ({FACT.get(base,'')}),  q_min = {q} "
              f"(smallest bad prime not dividing d),  s = n/{2*q},  X = {X:.3e}")
        print(f"  non-empty only for n <= {2*q-1} (every bad prime p with "
              f"2p <= n must divide d)")
        print(f"{'n':>4} {'s':>6} {'#d':>6} {'COUNT':>13} {'RUNSTART':>13} "
              f"{'MAXIMAL':>13} {'C/RS':>6} {'C/MX':>6} {'sec':>7}")
    q = QMIN.get(base)
    print(f"{n:4d} {n/(2.0*q):6.3f} {nd:6d} {c:13d} {rs:13d} {mx:13d} "
          f"{(c/rs if rs else 0):6.2f} {(c/mx if mx else 0):6.2f} {t:7.1f}")

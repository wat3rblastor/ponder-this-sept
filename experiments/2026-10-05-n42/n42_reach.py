"""Why an exact NONZERO count at n=42 in the base(58) family is out of reach on
this laptop.  exact.c's memory is 5 bitmaps of X bits + a prime sieve of X/2
bits = 0.6875*X bytes; the leanest possible variant (L + one working buffer,
periodic admissible-a pattern, no second scratch) is 0.25*X bytes.  Expected
COUNT is the validated model M2/M3 at n=42."""
import math, os, sys, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '2026-10-04-groundtruth'))
import model as M
HERE = os.path.dirname(os.path.abspath(__file__))
for f in sorted(glob.glob(os.path.join(HERE, 'dens_*.txt'))):
    for b, pts in M.load_dens(f).items():
        if pts:
            lo = min(t for t, r, c in pts)
            M.DENS[b] = sorted([p for p in M.DENS.get(b, []) if p[0] < lo] + pts)
B = 3741870
qs = [q for q in M.BAD if B % q == 0]
rho = M.Rho(B, extend='min6')
print(f"{'X':>10} {'mmax':>6} {'E[n=38]':>10} {'E[n=40]':>10} {'E[n=42]':>10} "
      f"{'exact.c GB':>11} {'lean GB':>9}")
for X in (2.5e9, 1e10, 2e10, 4e10, 1e11, 4e11):
    r = []
    for n in (38, 40, 42):
        C = M.corr_family(n, int(X), B, qs)
        r.append(M.term_integral(n, int(X), B, qs, rho) * C)
    print(f"{X:10.1e} {int((X-1)//(41*B)):6d} {r[0]:10.3f} {r[1]:10.3f} "
          f"{r[2]:10.3f} {0.6875*X/1e9:11.2f} {0.25*X/1e9:9.2f}")

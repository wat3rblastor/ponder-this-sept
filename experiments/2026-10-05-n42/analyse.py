"""Residual (model/exact) as a function of s = n/(2 q_min), pooled over the
new s-scan families and the 2026-10-04-groundtruth points.

Model = M2 and M3 of experiments/2026-10-04-groundtruth/model.py, i.e. the
exact L_q product form (corr_family / M3's exact d-sum) with the per-term rate
taken at each term's actual size and integrated over the (a,d) supply.  NOT
rho^n * 0.82^n.

Convention, throughout: COUNT = pairs (a,d) with a..a+(n-1)d all Loeschian and
a+(n-1)d <= X, i.e. APs of length AT LEAST n counted by starting pair (an AP of
length L >= n contributes L-n+1).  All APs, not maximal ones.  RUNSTART and
MAXIMAL are carried along for reference only.
"""
import math, os, sys, glob

HERE = os.path.dirname(os.path.abspath(__file__))
GT = os.path.join(os.path.dirname(HERE), '2026-10-04-groundtruth')
sys.path.insert(0, GT)
import model

# ---- splice our own measured class densities in (keyed by base) ----
for f in sorted(glob.glob(os.path.join(HERE, 'dens_*.txt'))):
    for base, pts in model.load_dens(f).items():
        old = model.DENS.get(base, [])
        # prefer the new, wider-range measurement; keep old points below its floor
        if pts:
            lo = min(t for t, r, c in pts)
            model.DENS[base] = sorted([p for p in old if p[0] < lo] + pts)

BADSET = set(model.BAD)


def qmin_of(base):
    for q in model.BAD:
        if base % q:
            return q


def parse(path):
    X = None
    rows = []
    for line in open(path):
        if line.startswith('# X='):
            X = int(line.split('=')[1])
        if line.startswith('n='):
            f = dict(kv.split('=', 1) for kv in line.split() if '=' in kv)
            if 'COUNT' not in f:
                continue
            rows.append(dict(n=int(f['n']), X=X, base=int(f['base']),
                             pairs=int(f['pairs']), rho1=float(f['rho1']),
                             count=int(f['COUNT']), rs=int(f['RUNSTART']),
                             mx=int(f['MAXIMAL']), src=os.path.basename(path)))
    return rows


def wls(xs, ys, ws):
    """weighted least squares y = a + b x; returns b, sigma(b), a, chi2/dof"""
    S = sum(ws); Sx = sum(w * x for w, x in zip(ws, xs))
    Sy = sum(w * y for w, y in zip(ws, ys))
    Sxx = sum(w * x * x for w, x in zip(ws, xs))
    Sxy = sum(w * x * y for w, x, y in zip(ws, xs, ys))
    D = S * Sxx - Sx * Sx
    if D == 0:
        return 0.0, float('inf'), 0.0, 0.0
    b = (S * Sxy - Sx * Sy) / D
    a = (Sy - b * Sx) / S
    dof = max(len(xs) - 2, 1)
    chi2 = sum(w * (y - a - b * x) ** 2 for w, x, y in zip(ws, xs, ys))
    s2 = chi2 / dof
    return b, math.sqrt(s2 * S / D), a, s2


def main():
    files = sorted(glob.glob(os.path.join(HERE, 'q*.txt')))
    files += [os.path.join(GT, f) for f in
              ('fam3741870_15e8.txt', 'fam3741870b.txt', 'fam3741870c.txt',
               'fam129030_1e9.txt', 'ext_1e8.txt', 'ext_1e9.txt',
               'repro_2e6.txt', 'repro_19e6.txt')]
    rows = []
    for f in files:
        if os.path.exists(f):
            rows += parse(f)

    out = []
    print('#' + '=' * 134)
    print('# EXACT COUNTS AND MODEL RATIOS.  COUNT = APs of length >= n counted by '
          'starting pair (all APs, not maximal).')
    print('#' + '=' * 134)
    hdr = (f"{'qmin':>4} {'base':>8} {'X':>10} {'n':>3} {'s':>6} {'#d':>6} "
           f"{'COUNT':>10} {'RUNSTART':>9} {'MAXIMAL':>8} {'Corr':>10} "
           f"{'M2':>11} {'M2/ex':>7} {'M2s/ex':>7} {'M3':>11} {'M3/ex':>7}  src")
    print(hdr)
    for r in rows:
        if r['count'] == 0:
            continue
        base, n, X = r['base'], r['n'], r['X']
        if base not in model.DENS:
            continue
        q = qmin_of(base)
        s = n / (2.0 * q)
        qs = [p for p in model.BAD if base % p == 0]
        rho = model.Rho(base)
        m = model.models(n, X, base=base, rho=rho, rho1_meas=r['rho1'])
        M3 = model.m3(n, X, base, qs, rho)
        # M2 uses the asymptotic supply delta*X^2/(2(n-1)base), which is up to
        # 9% high when mmax is small.  M2s divides that discretisation error out
        # by using the EXACT pair count, isolating the probability factor.
        M2s = m['M2'] * r['pairs'] / m['F']
        r.update(s=s, qmin=q, M2=m['M2'], M2s=M2s, M3=M3, Corr=m['Corr'], F=m['F'])
        out.append(r)
        print(f"{q:4d} {base:8d} {X:10.3e} {n:3d} {s:6.3f} "
              f"{int((X-1)//((n-1)*base)):6d} {r['count']:10d} {r['rs']:9d} "
              f"{r['mx']:8d} {m['Corr']:10.3e} {m['M2']:11.4e} "
              f"{m['M2']/r['count']:7.3f} {M2s/r['count']:7.3f} "
              f"{M3:11.4e} {M3/r['count']:7.3f}  {r['src']}")

    # ---------------- supply check ----------------
    worst = max(out, key=lambda r: abs(math.log(r['F'] / r['pairs'])))
    print(f"\nsupply law check: worst |ln(model pairs/exact pairs)| = "
          f"{abs(math.log(worst['F']/worst['pairs'])):.5f} "
          f"(base {worst['base']}, n={worst['n']})")

    # ---------------- residual vs s ----------------
    print('\n' + '=' * 90)
    print('RESIDUAL vs s = n/(2 q_min).   binned, weight = exact COUNT '
          '(Poisson), points with COUNT >= 30 only')
    print('=' * 90)
    good = [r for r in out if r['count'] >= 30]
    print(f"{'s bin':>12} {'pts':>4} {'sum COUNT':>11} {'<M2s/ex>':>9} "
          f"{'<M3/ex>':>8} {'spread M3':>18}")
    edges = [0, .2, .3, .4, .5, .6, .65, .7, .75, .8, .85, .9, 1.01]
    for lo, hi in zip(edges, edges[1:]):
        bs = [r for r in good if lo <= r['s'] < hi]
        if not bs:
            continue
        W = sum(r['count'] for r in bs)
        g2 = math.exp(sum(r['count'] * math.log(r['M2s'] / r['count']) for r in bs) / W)
        g3 = math.exp(sum(r['count'] * math.log(r['M3'] / r['count']) for r in bs) / W)
        rr = sorted(r['M3'] / r['count'] for r in bs)
        print(f"  [{lo:.2f},{hi:.2f}) {len(bs):4d} {W:11d} {g2:9.3f} {g3:8.3f} "
              f"   {rr[0]:.3f} .. {rr[-1]:.3f}")

    # per q_min, at s nearest 0.707
    print('\nresidual at the s nearest 0.707, per family (q_min):')
    print(f"{'qmin':>5} {'n':>4} {'s':>6} {'COUNT':>9} {'M2s/ex':>7} {'M3/ex':>7}")
    for q in sorted({r['qmin'] for r in good}):
        bs = [r for r in good if r['qmin'] == q]
        b = min(bs, key=lambda r: abs(r['s'] - 0.707))
        print(f"{q:5d} {b['n']:4d} {b['s']:6.3f} {b['count']:9d} "
              f"{b['M2s']/b['count']:7.3f} {b['M3']/b['count']:7.3f}")

    # ---------------- slope fits ----------------
    print('\n' + '=' * 90)
    print('SLOPE OF ln(residual) IN s   (weighted LS, weight = exact COUNT; '
          '95% CI = +-1.96 sigma)')
    print('=' * 90)
    for name, key in (('M2s', 'M2s'), ('M3', 'M3')):
        for lab, sel in (('all s, COUNT>=30', lambda r: True),
                         ('s <= 0.75', lambda r: r['s'] <= 0.75),
                         ('s in [0.4,1.0]', lambda r: r['s'] >= 0.4),
                         ('s in [0.6,1.0]', lambda r: r['s'] >= 0.6)):
            bs = [r for r in good if sel(r)]
            if len(bs) < 3:
                continue
            xs = [r['s'] for r in bs]
            ys = [math.log(r[key] / r['count']) for r in bs]
            ws = [float(r['count']) for r in bs]
            b, sb, a, c2 = wls(xs, ys, ws)
            print(f"  {name} {lab:18s} pts={len(bs):3d}  d ln(res)/ds = "
                  f"{b:+.4f} +- {sb:.4f}   res(s=0.707) = "
                  f"{math.exp(a + b * 0.707):.3f}   [95% CI "
                  f"{math.exp(a + (b - 1.96 * sb) * 0.707):.3f}.."
                  f"{math.exp(a + (b + 1.96 * sb) * 0.707):.3f}]  chi2/dof={c2:.1f}")

    # ---------------- q_min dependence at matched s ----------------
    print('\nq_min dependence at matched s (is s the right scaling variable?)')
    print(f"{'s window':>14} " + ' '.join(f"{q:>8}" for q in
                                         sorted({r['qmin'] for r in good})))
    for lo, hi in ((.6, .7), (.7, .8), (.8, .9), (.9, 1.01)):
        line = f"  [{lo:.2f},{hi:.2f})   "
        for q in sorted({r['qmin'] for r in good}):
            bs = [r for r in good if r['qmin'] == q and lo <= r['s'] < hi]
            if not bs:
                line += f"{'-':>8} "
                continue
            W = sum(r['count'] for r in bs)
            g = math.exp(sum(r['count'] * math.log(r['M3'] / r['count'])
                             for r in bs) / W)
            line += f"{g:8.3f} "
        print(line)

    # ---------------- the logically decisive ratio ----------------
    # base(58) can only be measured at s <= 0.44, but n=58 sits at s = 0.707.
    # The quantity that licenses carrying its measured residual across that gap
    # is  res(s=0.707) / res(s=0.44)  WITHIN a family.  Measure it in every
    # family that spans both.
    print('\n' + '=' * 90)
    print('WITHIN-FAMILY TRANSPORT RATIO  res(s=0.707)/res(s=0.44)')
    print('  (= the factor by which carrying base(58)\'s own measured residual')
    print('   from its last exact point, s=0.44, out to n=58 at s=0.707 is wrong)')
    print('=' * 90)
    print(f"{'qmin':>5} {'base':>8} {'n@0.44':>7} {'n@0.707':>8} {'res@0.44':>9} "
          f"{'res@0.707':>10} {'ratio':>7} {'COUNTs':>18}")
    rats = []
    for q in sorted({r['qmin'] for r in good}):
        bs = [r for r in good if r['qmin'] == q]
        lo = min(bs, key=lambda r: abs(r['s'] - 0.44))
        hi = min(bs, key=lambda r: abs(r['s'] - 0.707))
        if abs(lo['s'] - 0.44) > 0.09 or abs(hi['s'] - 0.707) > 0.09 or lo is hi:
            continue
        rl = lo['M3'] / lo['count']; rh = hi['M3'] / hi['count']
        rats.append(rh / rl)
        print(f"{q:5d} {lo['base']:8d} {lo['n']:7d} {hi['n']:8d} {rl:9.3f} "
              f"{rh:10.3f} {rh/rl:7.3f} {lo['count']:9d},{hi['count']:8d}")
    print('\nPER-FAMILY SLOPE over the family\'s OWN s range (COUNT >= 30):')
    print(f"{'qmin':>5} {'base':>8} {'pts':>4} {'s range':>14} {'d ln(res)/ds':>13} "
          f"{'sigma':>7} {'res lo':>7} {'res hi':>7} {'hi/lo':>6}")
    for q, bb in sorted({(r['qmin'], r['base']) for r in good}):
        bs = [r for r in good if r['qmin'] == q and r['base'] == bb]
        if len(bs) < 3:
            continue
        bs.sort(key=lambda r: r['s'])
        xs = [r['s'] for r in bs]
        ys = [math.log(r['M3'] / r['count']) for r in bs]
        ws = [float(r['count']) for r in bs]
        b, sb, a, c2 = wls(xs, ys, ws)
        lo = bs[0]['M3'] / bs[0]['count']; hi = bs[-1]['M3'] / bs[-1]['count']
        print(f"{q:5d} {bb:8d} {len(bs):4d} {xs[0]:6.3f}-{xs[-1]:6.3f} "
              f"{b:+13.4f} {sb:7.4f} {lo:7.3f} {hi:7.3f} {hi/lo:6.3f}")
    if rats:
        print(f"  -> transport ratio spans {min(rats):.3f} .. {max(rats):.3f} "
              f"over {len(rats)} families (geometric mean "
              f"{math.exp(sum(math.log(v) for v in rats)/len(rats)):.3f})")
    return out





def famsum():
    """count-weighted geometric-mean residual per family, for the README."""
    out = main()
    print('\n' + '=' * 90)
    print('COUNT-WEIGHTED GEOMETRIC-MEAN RESIDUAL PER FAMILY (COUNT >= 30)')
    print('=' * 90)
    print(f"{'qmin':>5} {'base':>8} {'X':>10} {'pts':>4} {'s range':>14} "
          f"{'sum COUNT':>12} {'<M2s/ex>':>9} {'<M3/ex>':>8}")
    key = {}
    for r in out:
        if r['count'] >= 30:
            key.setdefault((r['qmin'], r['base'], r['X']), []).append(r)
    for k in sorted(key):
        bs = key[k]
        W = sum(r['count'] for r in bs)
        g2 = math.exp(sum(r['count'] * math.log(r['M2s'] / r['count']) for r in bs) / W)
        g3 = math.exp(sum(r['count'] * math.log(r['M3'] / r['count']) for r in bs) / W)
        print(f"{k[0]:5d} {k[1]:8d} {k[2]:10.2e} {len(bs):4d} "
              f"{min(r['s'] for r in bs):6.3f}-{max(r['s'] for r in bs):6.3f} "
              f"{W:12d} {g2:9.3f} {g3:8.3f}")


if __name__ == '__main__':
    famsum()

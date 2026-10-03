#!/usr/bin/env python3
"""Build the M(n) table from all sweeps, compare with the heuristic, fit, extrapolate."""
import json, glob, math, os, collections
import heuristic as H

RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw')

def load():
    """-> {tag: {'base':b,'N':N,'mmax':m,'best':{n:(last,a,d)}}}"""
    sweeps = {}
    for fn in sorted(glob.glob(os.path.join(RAW, '*.jsonl'))):
        tag = os.path.basename(fn).split('.')[0]
        s = sweeps.setdefault(tag, {'base': None, 'N': 0, 'mmax': 0,
                                    'best': {}, 'ms': set()})
        for line in open(fn):
            line = line.strip()
            if not line: continue
            try: r = json.loads(line)
            except Exception: continue
            if 'm' not in r: continue
            s['N'] = r['N']; s['ms'].add(r['m'])
            s['base'] = r['d'] // r['m']
            for k, a in r.get('mina', {}).items():
                n = int(k); a = int(a)
                if a == 0: continue
                last = a + (n - 1) * r['d']
                if last > r['N']: continue
                if n not in s['best'] or last < s['best'][n][0]:
                    s['best'][n] = (last, a, r['d'])
    # mmax must be the CONTIGUOUS prefix 1..mmax actually scanned: parallel
    # chunks finish out of order, so max(m) would overstate the coverage.
    for s in sweeps.values():
        k = 0
        while (k + 1) in s['ms']: k += 1
        s['mmax'] = k
        # drop results whose d lies outside the contiguous prefix
        s['best'] = {n: v for n, v in s['best'].items()
                     if v[2] // s['base'] <= k}
    return sweeps

def table(sweeps):
    rows = {}
    for n in range(14, 60):   # below 14 the sweeps' --hist-min suppresses output
        bf = H.base_of(n)
        proved = None; bound = None
        for tag, s in sweeps.items():
            if n not in s['best']: continue
            last, a, d = s['best'][n]
            if bound is None or last < bound[0]: bound = (last, a, d, tag)
            complete = (bf % s['base'] == 0                   # sweep base divides forced base
                        and s['N'] >= last
                        and s['mmax'] * s['base'] * (n - 1) >= last)
            if complete and (proved is None or last < proved[0]):
                proved = (last, a, d, tag, s['N'], s['mmax'], s['base'])
        # a lower bound from any complete sweep that found nothing for this n
        lower = 0
        for tag, s in sweeps.items():
            if bf % s['base']: continue
            cov = min(s['N'], s['mmax'] * s['base'] * (n - 1))
            if n not in s['best'] or s['best'][n][0] > cov:
                lower = max(lower, cov)
        if proved or bound or lower:
            rows[n] = dict(base=bf, proved=proved, bound=bound, lower=lower)
    return rows

if __name__ == '__main__':
    sw = load()
    print("== sweeps ==")
    for t, s in sorted(sw.items()):
        print(f"  {t:9s} base={s['base']:<8d} N={s['N']:<12d} mmax={s['mmax']}")
    rows = table(sw)
    print("\n== M(n) ==")
    print(f"{'n':>3} {'base':>8} {'status':>6} {'M(n)':>14} {'a':>14} {'d':>10} "
          f"{'heur_ref':>11} {'ratio':>7} {'heur_stated':>12} {'r_st':>7}")
    data = []
    for n in sorted(rows):
        r = rows[n]
        hr = H.solve(n); hs = H.solve(n, refined=False, triangular=False, coprime_a=False)
        if r['proved']:
            last, a, d = r['proved'][:3]; st = 'PROVED'
            data.append((n, last, hr, hs))
        elif r['bound']:
            last, a, d = r['bound'][:3]; st = 'BOUND'
        elif r['lower']:
            print(f"{n:>3} {r['base']:>8} {'LOWER':>6} {'>'+str(r['lower']):>14} "
                  f"{'-':>14} {'-':>10} {H.solve(n):>11.3e} {'-':>7} "
                  f"{H.solve(n, refined=False, triangular=False, coprime_a=False):>12.3e}")
            continue
        else:
            continue
        print(f"{n:>3} {r['base']:>8} {st:>6} {last:>14d} {a:>14d} {d:>10d} "
              f"{hr:>11.3e} {last/hr:>7.2f} {hs:>12.3e} {last/hs:>7.2f}"
              + (f"   [>{r['lower']:.3g} excluded]" if r['lower'] and not r['proved'] else ""))
    json.dump({'rows': {str(n): rows[n] for n in rows}}, open('table.json', 'w'),
              indent=1, default=str)

    # ---- fits -------------------------------------------------------------
    print("\n== fits (proved points only) ==")
    import statistics
    def lsq(xs, ys):
        mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
        sxx = sum((x-mx)**2 for x in xs); sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
        b = sxy/sxx; a = my - b*mx
        res = [y - (a + b*x) for x, y in zip(xs, ys)]
        s = math.sqrt(sum(r*r for r in res)/max(1, len(xs)-2))
        seb = s/math.sqrt(sxx)
        return a, b, s, seb, mx, sxx, len(xs)

    for lo in (22, 26, 28):
        pts = [(n, math.log10(m)) for n, m, _, _ in data if n >= lo]
        if len(pts) < 3: continue
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        a, b, s, seb, mx, sxx, k = lsq(xs, ys)
        pred = a + b*58
        sep = s*math.sqrt(1 + 1/k + (58-mx)**2/sxx)   # prediction interval
        print(f"  log10 M ~ a+b*n, n>={lo}: b={b:.4f} dex/term "
              f"(x{10**b:.3f}/term), resid sd={s:.3f} dex, n_pts={k}")
        print(f"      -> log10 M(58) = {pred:.3f} +- {sep:.3f} (1 sd pred. int.)  "
              f"M(58) = {10**pred:.3e}  [{10**(pred-2*sep):.2e}, {10**(pred+2*sep):.2e}]")

    # ratio-to-heuristic fit (structure-aware: heuristic carries the base jumps)
    pts = [(n, math.log10(m/h)) for n, m, h, _ in data if n >= 22]
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    a, b, s, seb, mx, sxx, k = lsq(xs, ys)
    h58 = H.solve(58)
    pred = a + b*58; sep = s*math.sqrt(1 + 1/k + (58-mx)**2/sxx)
    print(f"\n  log10(M/heur_refined) ~ a+b*n, n>=22: b={b:+.4f} dex/term, "
          f"resid sd={s:.3f}, n_pts={k}")
    print(f"      mean ratio = {10**statistics.fmean(ys):.2f}x")
    print(f"      -> log10 M(58) = log10(heur={h58:.3e}) + {pred:.3f} +- {sep:.3f}")
    print(f"      M(58) = {h58*10**pred:.3e}  "
          f"[{h58*10**(pred-2*sep):.2e}, {h58*10**(pred+2*sep):.2e}]")
    # constant-ratio variant
    cm = statistics.fmean(ys); cs = statistics.stdev(ys)
    print(f"      constant-ratio variant: M(58) = {h58*10**cm:.3e} "
          f"[{h58*10**(cm-2*cs):.2e}, {h58*10**(cm+2*cs):.2e}]")

    # ---- why the structure-blind fit fails ---------------------------------
    print("\n== heuristic growth rate (dex per extra term) ==")
    prev = None
    for n in range(22, 60):
        t = H.solve(n)
        if prev: print(f"  n={n:2d} base={H.base_of(n):>8} heur={t:.3e} "
                       f"d(log10)/dn = {math.log10(t/prev):+.4f}")
        prev = t

    print("\n== three extrapolation routes for M(58) ==")
    h35, h58 = H.solve(35), H.solve(58)
    m35 = dict((n, m) for n, m, _, _ in data)[35]
    print(f"  (A) structure-blind log-linear fit on n>=26      : 6.0e14   "
          f"[REJECTED: assumes constant 0.27 dex/term, but the measured+modelled\n"
          f"       slope rises to 0.40-0.43 dex/term by n=55-58]")
    cm = statistics.fmean([math.log10(m/H.solve(n)) for n, m, _, _ in data if n >= 22])
    cs = statistics.stdev([math.log10(m/H.solve(n)) for n, m, _, _ in data if n >= 22])
    print(f"  (B) this model x constant measured bias ({10**cm:.2f}x)    : "
          f"{h58*10**cm:.2e}  [{h58*10**(cm-2*cs):.1e}, {h58*10**(cm+2*cs):.1e}]")
    print(f"  (C) this model x bias trend-extrapolated to n=58 : "
          f"{h58*10**(cm + (-0.0164)*(58-27.5)):.2e}")
    print(f"  (D) project heuristic (1e18 at 58) x measured/pred at n=35 "
          f"({m35/1.3e9:.3f}) : {1e18*m35/1.3e9:.2e}")
    print(f"\n  measured within-plateau growth n=26..35: "
          f"{math.log10(m35/dict((n,m) for n,m,_,_ in data)[26])/9:.4f} dex/term "
          f"(x{10**(math.log10(m35/dict((n,m) for n,m,_,_ in data)[26])/9):.3f})")

#!/usr/bin/env python3
"""Plan-wide comparison of the old (tools/plan_units.py) and corrected scores.
Random sample of units over K <= 60000, shift <= 32.  Stdlib only.
Run: MODCAP=2e16 python3 experiments/2026-10-03-planaudit/plan_compare.py
"""
import sys, os, math, random
sys.path.insert(0, 'tools'); sys.path.insert(0, os.path.dirname(__file__))
import plan_units as PU, score2 as S2

def spearman(x, y):
    def rk(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
        while i < len(o):
            j = i
            while j+1 < len(o) and v[o[j+1]] == v[o[i]]: j += 1
            for t in range(i, j+1): r[o[t]] = (i+j)/2 + 1
            i = j+1
        return r
    a, b = rk(x), rk(y); n = len(x); ma, mb = sum(a)/n, sum(b)/n
    return sum((p-ma)*(q-mb) for p, q in zip(a, b))/math.sqrt(
        sum((p-ma)**2 for p in a)*sum((q-mb)**2 for q in b))

random.seed(11)
U = set()
while len(U) < 3000:
    U.add((random.randint(1, 60000), random.choice([0, 0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32])))
rows = []
for K, s in sorted(U):
    o = S2.v_old(K, s); n = S2.v_new(K, s, 2e16)
    if o is None or n is None: continue
    assert o[1] == n[1]
    rows.append(dict(K=K, shift=s, v_old=o[0], v_new=n[0], res=float(o[1]),
                     wbar=PU.unit_shape(K)[3]))
vo = [r["v_old"] for r in rows]; vn = [r["v_new"] for r in rows]
res = [r["res"] for r in rows]; tot = sum(res)
d = [a-b for a, b in zip(vo, vn)]; md = sum(d)/len(d)
sd = math.sqrt(sum((x-md)**2 for x in d)/len(d))
print("units %d   spearman(v_old, v_new) = %.6f" % (len(rows), spearman(vo, vn)))
print("v_old spread %.2f log (= %.1fx)" % (max(vo)-min(vo), math.exp(max(vo)-min(vo))))
print("v_old - v_new: constant %.3f, sd %.4f, worst unit %.3f log (= %.2fx)"
      % (md, sd, max(abs(x-md) for x in d), math.exp(max(abs(x-md) for x in d))))
for sel in (0, None):
    idx = [i for i, r in enumerate(rows) if sel is None or r["shift"] == sel]
    dd = [d[i]-md for i in idx]
    print("  shift==%s: n=%d sd %.4f worst %.3f" % (sel, len(idx),
          math.sqrt(sum(x*x for x in dd)/len(dd)), max(abs(x) for x in dd)))

def top(v, frac):
    """fractional knapsack: exactly frac*tot residues, last unit split"""
    idx = sorted(range(len(rows)), key=lambda i: -v[i]); b = frac*tot; t = 0.0; sel = []
    for i in idx:
        if t >= b: break
        w = min(res[i], b - t); sel.append((i, w)); t += w
    return sel

for frac in (0.05, 0.10, 0.25):
    to, tn = top(vo, frac), top(vn, frac)
    so, sn = set(i for i, w in to), set(i for i, w in tn)
    yo = sum(math.exp(vn[i])*w for i, w in to)
    yn = sum(math.exp(vn[i])*w for i, w in tn)
    print("top %4.0f%% of residues: |old|=%d |new|=%d  overlap %.1f%%  "
          "model yield ratio new/old = %.4f"
          % (100*frac, len(so), len(sn), 100*len(so & sn)/max(len(so), len(sn)), yn/yo))
print("\nhead of each ranking (K, shift):")
ro = sorted(range(len(rows)), key=lambda i: -vo[i])[:8]
rn = sorted(range(len(rows)), key=lambda i: -vn[i])[:8]
print("  old:", [(rows[i]["K"], rows[i]["shift"]) for i in ro])
print("  new:", [(rows[i]["K"], rows[i]["shift"]) for i in rn])

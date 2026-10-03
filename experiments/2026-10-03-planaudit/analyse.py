#!/usr/bin/env python3
"""Correlate the old (tools/plan_units.py) and corrected (score2.py) unit scores
with the measured per-unit yield from probe.c.  Stdlib only.

measured log yield per stage-1 residue:
    M = log(win_ok / res) + sum_k log p_k      (probe.c "sumlogp" + the first term)
sigma(M) is propagated from the binomial counts.

Run:  MODCAP=2e16 python3 experiments/2026-10-03-planaudit/analyse.py
"""
import json, math, os, sys
sys.path.insert(0, 'tools'); sys.path.insert(0, os.path.dirname(__file__))
import plan_units as PU, score2 as S2

HERE = os.path.dirname(__file__)
meas = [json.loads(l) for l in open(os.path.join(HERE, "measured.jsonl")) if l.strip().startswith("{")]


def spearman(x, y):
    def rk(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]: j += 1
            for t in range(i, j + 1): r[o[t]] = (i + j) / 2 + 1
            i = j + 1
        return r
    a, b = rk(x), rk(y); n = len(x)
    ma, mb = sum(a) / n, sum(b) / n
    return sum((p - ma) * (q - mb) for p, q in zip(a, b)) / math.sqrt(
        sum((p - ma) ** 2 for p in a) * sum((q - mb) ** 2 for q in b))


def pearson(x, y):
    n = len(x); mx, my = sum(x) / n, sum(y) / n
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(
        sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))


rows = []
for m in meas:
    K, s = m["K"], m["shift"]
    if not m["win_ok"] or m["sumlogp"] != m["sumlogp"]:
        continue
    o = S2.v_old(K, s); n = S2.v_new(K, s, 2e16); nc = S2.v_new(K, s, 2e16, wbar=PU.unit_shape(K)[3])
    if o is None or n is None:
        continue
    M = math.log(m["win_ok"] / m["res"]) + m["sumlogp"]
    p = m["nloe"] / m["nterm"]
    nwin = m["win_ok"]
    sig = math.sqrt(58 * (1 - p) / (p * nwin) + 1.0 / nwin)
    rows.append(dict(K=K, shift=s, M=M, sig=sig, v_old=o[0], v_new=n[0],
                     v_new_cudaw=nc[0], res=float(o[1]), nwin=nwin,
                     nterm=m["nterm"], p=p, MOD=m["MOD"]))

print("units measured: %d   median sigma(M) = %.3f   (score spread %.2f log)" % (
    len(rows), sorted(r["sig"] for r in rows)[len(rows) // 2],
    max(r["v_old"] for r in rows) - min(r["v_old"] for r in rows)))
M = [r["M"] for r in rows]
for nm in ("v_old", "v_new", "v_new_cudaw"):
    v = [r[nm] for r in rows]
    print("  spearman(M, %-12s) = %+.4f   pearson = %+.4f  sd(M - %s) = %.3f" % (
        nm, spearman(M, v), pearson(M, v), nm,
        (lambda d: math.sqrt(sum((x - sum(d) / len(d)) ** 2 for x in d) / len(d)))(
            [a - b for a, b in zip(M, v)])))

# slope of M on v (1.0 == the model's size/prime factors have the right scale)
def ols(x, y):
    n = len(x); mx, my = sum(x) / n, sum(y) / n
    b = sum((a - mx) * (c - my) for a, c in zip(x, y)) / sum((a - mx) ** 2 for a in x)
    r = [c - (my + b * (a - mx)) for a, c in zip(x, y)]
    se = math.sqrt(sum(z * z for z in r) / (n - 2) / sum((a - mx) ** 2 for a in x))
    return b, se

for nm in ("v_old", "v_new"):
    b, se = ols([r[nm] for r in rows], M)
    print("  OLS slope M on %s: %.3f +- %.3f" % (nm, b, se))

# does the residual still depend on something the score ignores?
res_v = [r["M"] - r["v_new"] for r in rows]
for nm, f in (("shift", lambda r: r["shift"]), ("ln MOD", lambda r: math.log(r["MOD"])),
              ("ln res", lambda r: math.log(r["res"])), ("wbar", lambda r: PU.unit_shape(r["K"])[3]),
              ("ln K", lambda r: math.log(r["K"]))):
    print("  pearson(residual M-v_new, %-7s) = %+.3f" % (nm, pearson(res_v, [f(r) for r in rows])))

# ---- the DECIDES statistic: measured yield of the top 10% of residues by each score
tot = sum(r["res"] for r in rows)
def top(key, frac=0.10):
    """fractional knapsack: exactly frac*tot residues, last unit split"""
    idx = sorted(rows, key=lambda r: -r[key]); b = frac * tot; t = 0.0; sel = []
    for r in idx:
        if t >= b: break
        w = min(r["res"], b - t); sel.append((r, w)); t += w
    return sel, t
print("\nmeasured yield of the best 10%% of the sampled residue budget:")
out = {}
for nm in ("v_old", "v_new", "M"):
    sel, spent = top(nm)
    y = sum(math.exp(r["M"]) * w for r, w in sel)
    out[nm] = y
    print("  ranked by %-6s: units=%2d residues=%.3g measured_yield=%.4g" % (nm, len(sel), spent, y))
print("  RATIO new/old = %.3f   (oracle M/old = %.3f)" % (out["v_new"] / out["v_old"], out["M"] / out["v_old"]))
for frac in (0.2, 0.3, 0.5):
    a = sum(math.exp(r["M"]) * w for r, w in top("v_new", frac)[0])
    b = sum(math.exp(r["M"]) * w for r, w in top("v_old", frac)[0])
    print("  at frac %.1f: ratio new/old = %.3f" % (frac, a / b))
json.dump(rows, open(os.path.join(HERE, "joined.json"), "w"), indent=0)

#!/usr/bin/env python3
"""Real-vs-control comparison of per-term features, with p-values."""
from __future__ import annotations

import json
import math
from collections import defaultdict

import sys
RF = sys.argv[1] if len(sys.argv)>1 else "real_terms.json"
CF = sys.argv[2] if len(sys.argv)>2 else "ctrl_terms.json"
real = json.loads(open(RF).read())
ctrl = json.loads(open(CF).read())
for x in real+ctrl:
    x.setdefault("rec", "corpus")

BIN = ["is_auto", "is_auto2", "is_csquare", "good_sq_big", "smooth1e4", "smooth1e6"]
NUM = ["lp_ratio", "rad_ratio", "ndistinct", "nfac", "bad_distinct"]


def binom_tail(k, n, p):
    """P[X >= k] (or <= k if k below mean) for Binomial(n,p); two-sided-ish."""
    if n == 0:
        return 1.0
    lo = k <= n * p
    tot = 0.0
    rng = range(0, k + 1) if lo else range(k, n + 1)
    for i in rng:
        lg = (math.lgamma(n+1)-math.lgamma(i+1)-math.lgamma(n-i+1)
               + i*math.log(p) + (n-i)*math.log1p(-p))
        tot += math.exp(lg)
    return min(1.0, 2 * tot)


print(f"{'feature':<14} {'real':>14} {'ctrl':>14} {'ratio':>7} {'p(2s)':>9}")
print("-" * 64)
for f in BIN:
    rk = sum(x[f] for x in real)
    rn = len(real)
    ck = sum(x[f] for x in ctrl)
    cn = len(ctrl)
    pc = ck / cn
    ratio = (rk / rn) / pc if pc > 0 else float("inf")
    p = binom_tail(rk, rn, pc) if pc > 0 else (1.0 if rk == 0 else 0.0)
    print(f"{f:<14} {rk}/{rn} = {rk/rn:.4f} {ck}/{cn} = {pc:.4f} "
          f"{ratio:>6.2f}x {p:>9.2g}")

print()
print(f"{'feature':<14} {'real mean':>10} {'ctrl mean':>10} {'diff/sd':>9} {'z':>8}")
print("-" * 56)
for f in NUM:
    rv = [x[f] for x in real]
    cv = [x[f] for x in ctrl]
    rm = sum(rv) / len(rv)
    cm = sum(cv) / len(cv)
    cs = math.sqrt(sum((v - cm) ** 2 for v in cv) / (len(cv) - 1))
    z = (rm - cm) / (cs / math.sqrt(len(rv)))
    print(f"{f:<14} {rm:>10.4f} {cm:>10.4f} {(rm-cm)/cs:>9.3f} {z:>8.2f}")

# per-record breakdown of the headline binaries
print("\nper-record:")
for f in ("is_auto", "is_csquare", "good_sq_big"):
    for rec in sorted({x["rec"] for x in real}):
        rr = [x for x in real if x["rec"] == rec]
        cc = [x for x in ctrl if x["rec"] == rec]
        print(f"  {f:<12} {rec:<6} real {sum(x[f] for x in rr)}/{len(rr)}  "
              f"ctrl {sum(x[f] for x in cc)}/{len(cc)} = "
              f"{sum(x[f] for x in cc)/len(cc):.4f}")

# which c values show up
cnt = defaultdict(int)
for x in real:
    if x["csquare"]:
        cnt[x["csquare"]] += 1
print("\nreal c with t = c*square:", dict(sorted(cnt.items())))
cnt = defaultdict(int)
for x in ctrl:
    if x["csquare"]:
        cnt[x["csquare"]] += 1
print("ctrl c with t = c*square:", dict(sorted(cnt.items())))
cnt = defaultdict(int)
for x in real:
    if x["auto"]:
        cnt[x["auto"]] += 1
print("real c with 3*c*d*t square:", dict(sorted(cnt.items())))
cnt = defaultdict(int)
for x in ctrl:
    if x["auto"]:
        cnt[x["auto"]] += 1
print("ctrl c with 3*c*d*t square:", dict(sorted(cnt.items())))

# index position of the special terms
print("\nindices k of real terms with good_sq_big / is_csquare:")
for rec in sorted({x["rec"] for x in real}):
    rr = [x for x in real if x["rec"] == rec]
    print(f"  {rec}: good_sq k={[x['k'] for x in rr if x['good_sq_big']]} "
          f"csq k={[x['k'] for x in rr if x['is_csquare']]}")

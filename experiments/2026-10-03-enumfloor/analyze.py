#!/usr/bin/env python3
"""Combine the measured enumeration cost (bench.out) with the 58-term pass model
to get expected 58-term hits per CPU-second.  The cost numbers are MEASURED; the
P58 model is INHERITED from the briefing and is not verified here."""
import re, os
from math import prod, log, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
L = 58
K_LR = 0.638909405
DIVD = [2,5,11,17,23,29,41,47,53]      # bad primes dividing D0 (coprime to a)

def P58(T, Q):
    F = prod([p/(p+1) for p in DIVD + list(Q)])
    return (2*K_LR/(sqrt(log(T))*F))**L

# --- measured ns/output, pulled out of bench.out ---------------------------
rows = []
txt = open(os.path.join(HERE, "bench.out")).read().split("### ")
for blk in txt[1:]:
    cmd = blk.splitlines()[0]
    m = re.search(r"out=(\d+) cand=(\d+) leaves=(\d+) sec=([\d.]+) work_per_out=([\d.]+) ns_per_out=([\d.]+)", blk)
    if not m: continue
    if "--verify" in cmd: continue
    Q = [int(x) for x in re.search(r"--Q ([\d,]+)", cmd).group(1).split(",")]
    X = int(re.search(r"--X (\d+)", cmd).group(1))
    rows.append(dict(Q=Q, X=X, cmd=cmd, out=int(m.group(1)),
                     work=float(m.group(5)), ns=float(m.group(6)),
                     split="--split" in cmd))

D0 = 382160924970
XMIN = 57*D0            # terms are >= 57*d and d >= D0, so X < XMIN is infeasible
print("measured cost, and expected 58-term hits per CPU-second (1 core)")
print("(rows with X < 57*D0 = %.3g are INFEASIBLE, marked *, kept only to show"
      " the cost curve)" % XMIN)
print("%-24s %-10s %-9s %-8s %-8s %-11s %-11s" %
      ("pins", "X (term)", "1/delta", "work/out", "ns/out", "P58", "hits/core-s"))
best = None
for r in rows:
    if r["split"]: continue
    delta = prod([(q-L)/q for q in r["Q"]])
    p = P58(r["X"], r["Q"])
    hps = p / (r["ns"]*1e-9)
    tag = "59..%d" % r["Q"][-1] + ("*" if r["X"] < XMIN else "")
    print("%-24s %-10.3g %-9.0f %-8.2f %-8.2f %-11.4g %-11.4g" %
          (tag, r["X"], 1/delta, r["work"], r["ns"], p, hps))
    if r["X"] >= XMIN and (best is None or hps > best[0]):
        best = (hps, tag, r["X"], r["ns"], p)

# the current engine: same pins, term floor 64*MOD, 1 candidate per iteration
eng = [r for r in rows if r["X"] == 72554321153880960][0]
delta_e = prod([(q-L)/q for q in eng["Q"]])
p_e = P58(eng["X"], eng["Q"])
hps_e = p_e/(eng["ns"]*1e-9)
print("\nbaseline (today's engine structure: pins 59..113, term floor 64*MOD = %.3g)" % eng["X"])
print("  ns/out = %.2f  P58 = %.4g  hits/core-s = %.4g" % (eng["ns"], p_e, hps_e))
print("\nbest measured configuration: pins %s at X = %.3g" % (best[1], best[2]))
print("  ns/out = %.2f  P58 = %.4g  hits/core-s = %.4g" % (best[3], best[4], best[0]))
print("  SPEEDUP over today's structure = %.1f x" % (best[0]/hps_e))

# ---------------------------------------------------------------------------
# What the win is actually worth: sweep K (d = K*D0, term size T = 57*K*D0).
# At every K the wheel modulus can stay at M1 = 59*71*83*89*101*107 = 3.344e11
# (<= T for all K >= 1), so the cost stays ~1.6 ns per candidate.  The engine's
# structure instead floors the term size at 64*MOD = 7.26e16 regardless of K.
# Approximation (flagged): K's own prime factors are assumed to add no new
# constraints, and a-supply at K is delta_total * 57*K*D0.
SMALL = (1/3)*(1/2)*(4/5)
Q107 = [59,71,83,89,101,107]
delta107 = prod([(q-L)/q for q in Q107]) * SMALL
NS_NEW, NS_OLD = 1.64, 1.43          # measured
Q113 = [59,71,83,89,101,107,113]
delta113 = prod([(q-L)/q for q in Q113]) * SMALL
ENG_X = 72554321153880960

print("\n-- value of the win: exhaustive sweep of a < 57*K*D0 at ~1.64 ns/cand --")
print("%-9s %-10s %-11s %-11s %-11s %-11s" %
      ("Kmax","term size","cands<=Kmax","core-hours","E[58s] new","E[58s] old"))
cum_n = cum_hit = cum_old = 0.0
for K in range(1, 400001):
    T = 57*K*D0
    n = delta107*57*D0            # increment of supply when K -> K+1
    cum_n += n
    cum_hit += n*P58(T, Q107)
    cum_old += n*P58(max(T, ENG_X), Q113)*(delta113/delta107)
    if K in (1, 10, 100, 1000, 3331, 10000, 100000, 400000):
        print("%-9d %-10.3g %-11.4g %-11.4g %-11.4g %-11.4g" %
              (K, T, cum_n, cum_n*NS_NEW*1e-9/3600, cum_hit, cum_old))
print("\n(last column = the same a-supply priced at the engine's term floor and"
      " pin set, i.e. what today's structure would get from the same work)")

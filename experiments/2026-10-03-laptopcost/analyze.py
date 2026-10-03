#!/usr/bin/env python3
"""Measured extrapolation of the apsearch 58-term cost on this laptop.

Inputs:
  w*.jsonl / w*.log  -- 8-core run, nterms=58, report=30 (throughput, window
                        count, and the empirical tail histogram of run lengths)
  q_scan.txt         -- mc_win: per-term pass probability inside an admissible
                        window, measured on the real AP terms

Model used for the extrapolation (validated against the tail histogram):
  inside an admissible window every term is free of every bad prime <= b2 and
  is 1 mod 3; each bad prime q > b2 exceeds 58 so it divides at most one of the
  58 indices, hence the 58 Loeschian events are independent Bernoulli(q).
  Terms immediately outside the window are only "1 mod 3 + coprime to the bad
  primes dividing d", pass probability p_out, and the engine extends outward.
"""
import glob, json, re, sys, math

HERE = __import__("os").path.dirname(__import__("os").path.abspath(__file__))

# ---------------------------------------------------------------- engine run
res = surv = 0
secs = 0.0
hist = {}
hits = set()
for f in sorted(glob.glob(HERE + "/w*.jsonl")):
    last = None          # last progress/unit record per (K,shift)
    per = {}
    for line in open(f):
        r = json.loads(line)
        if r.get("hit"):
            hits.add((r["a"], r["d"], r["n"]))
            continue
        per[(r["K"], r["shift"])] = r
    for r in per.values():
        res += r["res"]; surv += r["surv"]; secs += r["secs"]

for a, d, n in hits:
    for L in range(1, n + 1):
        hist[L] = hist.get(L, 0) + 1

wall = max(float(re.search(r"(\d+\.?\d*)s", l).group(1))
           for f in glob.glob(HERE + "/w*.log")
           for l in open(f) if l.startswith("  ..") or l.startswith("K="))

# ---------------------------------------------------------------- mc_win q
qs, ws, pouts = [], [], []
for l in open(HERE + "/q_scan.txt"):
    m = dict(kv.split("=") for kv in l.split())
    fr, lo = int(m["free"]), int(m["loe"])
    qs.append(lo / fr); ws.append(fr)
    pouts.append(lo / int(m["terms"]))      # unprotected-term pass rate
q = sum(a * b for a, b in zip(qs, ws)) / sum(ws)
sq = math.sqrt(q * (1 - q) / sum(ws))
p_out = sum(pouts) / len(pouts)

# ------------------------------------------------- exact run-length DP model
def p_run(L, q, p_out, nwin=58):
    """P(the engine reports a run >= L from one admissible window).

    Bi-infinite AP: indices 0..nwin-1 pass with prob q, all others with p_out.
    Success = some run of >= L consecutive passes whose span meets the window.
    DP over the index line carrying the current run length (capped at L)."""
    st = [0.0] * (L + 1); st[0] = 1.0; done = 0.0
    for i in range(-L - 1, nwin + L + 1):
        p = q if 0 <= i < nwin else p_out
        nxt = [0.0] * (L + 1)
        for r, m in enumerate(st):
            if m == 0.0: continue
            nxt[0] += m * (1 - p)                      # fail -> run resets
            r2 = r + 1
            if r2 >= L:
                # run span [i-L+1, i]; counts iff it meets [0, nwin-1]
                if i >= 0 and i - L + 1 <= nwin - 1: done += m * p
                else: nxt[L] += m * p
            else:
                nxt[r2] += m * p
        st = nxt
    return done

# ------------------------------------------------------------------- report
print("=== measured engine throughput (8 cores, nterms=58, report=30) ===")
print(f"residues processed      {res:.4g}")
print(f"core-seconds            {secs:.1f}   wall {wall:.0f}s")
print(f"residues/s/core         {res/secs:.4g}")
print(f"residues/s (8 cores)    {8*res/secs:.4g}")
print(f"admissible windows      {surv}   = {surv/res:.4g} per residue")
print(f"windows/s (8 cores)     {8*surv/secs:.4g}")
print()
print("=== measured per-term decay inside an admissible window ===")
print(f"q      = {q:.5f} +- {sq:.5f}   (mc_win, {sum(ws)} real AP terms, 4 K values)")
print(f"p_out  = {p_out:.5f}           (term just outside the sieved window)")
print()
print("=== empirical tail histogram vs the Bernoulli(q) model ===")
print(f"{'L':>4} {'obs N(>=L)':>11} {'pred N(>=L)':>12} {'obs/pred':>9} "
      f"{'r=pred(L+1)/pred(L)':>21}")
rows = []
for L in range(30, 37):
    pL = p_run(L, q, p_out)
    pred = surv * pL
    obs = hist.get(L, 0)
    r = p_run(L + 1, q, p_out) / pL
    rows.append((L, obs, pred, r))
    print(f"{L:>4} {obs:>11} {pred:>12.2f} "
          f"{(obs/pred if pred else float('nan')):>9.2f} {r:>21.4f}")
print()
print("decay ratio r(L) = P(>=L+1)/P(>=L) across the extrapolation range:")
for L in list(range(30, 58, 4)) + [56, 57]:
    print(f"   L={L:>3}  r={p_run(L+1,q,p_out)/p_run(L,q,p_out):.4f}")
print()

# ------------------------------------------------------------- extrapolation
p58 = p_run(58, q, p_out)
p58_pure = q ** 58
print("=== extrapolation to L = 58 ===")
print(f"P(all 58 in-window Loeschian)  = q^58            = {p58_pure:.4g}")
print(f"P(engine reports run >= 58)    = DP with extension = {p58:.4g}"
      f"   ({p58/p58_pure:.2f}x q^58)")
wps8 = 8 * surv / secs
rate8 = wps8 * p58
print(f"58-term solutions per second (8 laptop cores) = {rate8:.4g}")
print(f"58-term solutions per laptop-core-hour        = {rate8*3600/8:.4g}")
print(f"EXPECTED WALL-CLOCK, this laptop              = {1/rate8/3600:.4g} hours"
      f"  = {1/rate8/3600/24/365.25:.3g} years")
print()
# 1-sigma band from q alone
for sgn, lab in ((+1, "q+1sigma"), (-1, "q-1sigma")):
    qq = q + sgn * sq
    h = 1 / (wps8 * p_run(58, qq, p_out)) / 3600
    print(f"  {lab}: q={qq:.5f} -> {h:.4g} hours")
# band from the engine tail alone: fit q so pred N(>=30) == obs
obs30 = hist.get(30, 0)
lo = 0.55; hi = 0.75
for _ in range(60):
    mid = (lo + hi) / 2
    if surv * p_run(30, mid, p_out) < obs30: lo = mid
    else: hi = mid
qfit = (lo + hi) / 2
sig = math.sqrt(obs30) / obs30 if obs30 else 1.0
print(f"\n  q fitted from the engine's own N(>=30)={obs30}: q_fit={qfit:.4f}"
      f"  (Poisson 1sigma on the count = {100*sig:.0f}% -> "
      f"dq/q = {100*sig/30:.2f}%)")
for sgn in (+1, -1):
    qq = qfit * (1 + sgn * sig / 30)
    print(f"    q_fit{'+' if sgn>0 else '-'}1sigma = {qq:.4f} -> "
          f"{1/(wps8*p_run(58,qq,p_out))/3600:.4g} hours")

# --------------------------------------------------------- GPU cross-check
GPU_RES_S = 7.45e11
gpu_rate = GPU_RES_S * (surv / res) * p58
print()
print("=== cross-check against the 8x RTX PRO 6000 datapoint ===")
print(f"GPU residues/s (given)            {GPU_RES_S:.4g}")
print(f"laptop windows per residue        {surv/res:.4g}  (measured)")
print(f"=> GPU windows/s                  {GPU_RES_S*surv/res:.4g}")
print(f"=> predicted 58s per GPU-hour     {gpu_rate*3600:.3f}")
print(f"   claimed  58s per GPU-hour      0.4")
print(f"   ratio predicted/claimed        {gpu_rate*3600/0.4:.2f}")
print(f"laptop/GPU throughput ratio       {8*res/secs/GPU_RES_S:.4g}")
print(f"=> laptop hours implied by the 0.4/h figure: "
      f"{1/(0.4*(8*res/secs)/GPU_RES_S):.4g}")

# ------------------------------------------------- planner rho comparison
print()
print("=== the planner's per-term rho vs measurement ===")
lnt = 40.4   # mean ln(term) over the window actually searched
rho = 0.71992 * (lnt / 39.144) ** -0.485
print(f"tools/plan_units.py (mathres) rho at ln t={lnt}: {rho:.5f}")
print(f"measured q                                      : {q:.5f}")
print(f"ratio (meas/planner) = {q/rho:.4f}  -> ^58 = {(q/rho)**58:.3g}")
print(f"=> if rho were right, N(>=30) would have been "
      f"{surv*p_run(30, rho, p_out):.0f}; observed {hist.get(30,0)}.")

# ------------------------------------------------- tail-only sensitivity
print()
print("=== how much the engine tail alone constrains the answer ===")
n30, n33 = hist.get(30, 0), hist.get(33, 0)
r_obs = (n33 / n30) ** (1 / 3)
print(f"observed cumulative decay 30->33: r_obs = ({n33}/{n30})^(1/3) = {r_obs:.3f}")
print(f"model (r = q):                            {(p_run(33,q,p_out)/p_run(30,q,p_out))**(1/3):.3f} per term")
for r in (r_obs, 0.6247, 0.60, r_obs*0.93):
    p = (n33 / surv) * r ** 25
    print(f"  extrapolate P(>=33) by r={r:.3f}^25 -> P(>=58)={p:.3g} -> "
          f"{1/(wps8*p)/3600:.4g} laptop-hours")

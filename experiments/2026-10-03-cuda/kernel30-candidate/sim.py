"""CPU-only statistical model of sieve_strided: per-group death row, warp max, etc.
Reproduces the host table construction of apsearch_cuda.cu for one (K, shift)."""
import numpy as np, sys
from math import prod
K = int(sys.argv[1]) if len(sys.argv) > 1 else 894
SHIFT = int(sys.argv[2]) if len(sys.argv) > 2 else 1
NCH = 7; U = 4; NT = 58; D0 = 382160924970; MODCAP = 20000000000000000; B2 = 10000
WARPS = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
rng = np.random.default_rng(1)
def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
P = primes(B2)
d = K * D0
C = [(3, 1, 0, 1), (2, 1, 0, 1), (5, 1, 1, 4)]
MOD = 30
for q in P:
    if q <= NT or q % 3 != 2 or d % q == 0: continue
    if MOD * q > MODCAP: break
    C.append((q, d % q, d % q, q - NT)); MOD *= q
R0 = 0; sst = []
for (m, s, t, c) in C:
    co = MOD // m; e = co * pow(co % m, -1, m) % MOD
    R0 = (R0 + s * e) % MOD; sst.append(t * e % MOD)
idx = sorted(range(len(C)), key=lambda i: C[i][3])
in2, in1 = idx[-1], idx[-2]
c1, c2 = C[in1][3], C[in2][3]; s1, s2 = sst[in1], sst[in2]
pinned = {c[0] for c in C}
tc = []
for r in P:
    if r in pinned: continue
    divD = D0 % r == 0
    if not divD and (r % 3 != 2 or d % r == 0 or r <= NT): continue
    if s2 % r == 0: continue
    tc.append(r)
tc.sort(key=lambda r: -((1.0 / r) if D0 % r == 0 else NT / r))
print(f"K={K} MOD={MOD:.4g} c1={c1} c2={c2} ntc={len(tc)} first={tc[:12]}")
boff = SHIFT * 64
tabs = []
for r in tc:
    ok = np.ones(r, bool)
    if D0 % r == 0: ok[0] = False
    else:
        dm = d % r
        for k in range(NT): ok[(-k * dm) % r] = False
    modr = MOD % r; b0 = (boff % r) * modr % r
    x = np.arange(r)
    W = np.zeros(r, np.uint64)
    for b in range(64):
        W |= ok[(x + b0 + b * modr) % r].astype(np.uint64) << np.uint64(b)
    ys = (np.arange(r + NCH - 1) * (s2 % r)) % r
    tabs.append((r, W[ys], pow(s2 % r, -1, r)))
# sample warps: 32 lanes with independent random prefixes, common (i1, group)
others = [i for i in range(len(C)) if i not in (in1, in2) and C[i][3] > 1]
NG = WARPS * 32
R = np.full(NG, R0, dtype=object)
for i in others:
    R = (R + rng.integers(0, C[i][3], NG).astype(object) * sst[i]) % MOD
i1 = np.repeat(rng.integers(0, c1, WARPS), 32).astype(object)
g = np.repeat(rng.integers(0, (c2 + NCH - 1) // NCH, WARPS), 32).astype(object)
A = R + i1 * s1 + g * NCH * s2
nvalid = np.minimum(NCH, c2 - np.array(g, dtype=np.int64) * NCH)
death = np.full(NG, len(tabs), np.int64)
alive = np.ones(NG, bool)
acc = np.zeros((NG, NCH), np.uint64)
for c in range(NCH): acc[:, c] = np.where(c < nvalid, ~np.uint64(0), np.uint64(0))
Amod = {}
for t, (r, Wr, inv) in enumerate(tabs):
    am = np.array([int(a) % r for a in A], np.int64)  # (A mod r)
    y = am * inv % r
    for c in range(NCH): acc[:, c] &= Wr[y + c]
    dead = (acc == 0).all(1) & alive
    death[dead] = t; alive &= ~dead
    if not alive.any(): break
    if t > 200 and alive.sum() < 3: pass
print("groups surviving all rows:", alive.sum())
rows = death + 1                               # rows actually needed
iters = (rows + U - 1) // U * U               # rounded to U
wm = iters.reshape(WARPS, 32).max(1)
print(f"mean rows/group={rows.mean():.2f} (U-rounded {iters.mean():.2f}), median={np.median(rows)}")
print(f"warp max rows: mean={wm.mean():.2f}  p10={np.percentile(wm,10)} p90={np.percentile(wm,90)}")
print(f"divergence waste factor (warp max / lane mean) = {wm.mean()/iters.mean():.3f}")
for q in [4, 8, 12, 16, 20, 24, 32, 40, 48, 64]:
    print(f"  frac groups alive after {q:3d} rows: {(rows > q).mean():.4f}")
# two-phase / compaction models (cost in warp-rows per 32 groups)
cur = wm.mean()
print(f"current cost per 32 groups: {cur:.2f} warp-rows")
for T0 in [8, 12, 16, 20, 24, 28, 32]:
    # ideal stage compaction at T0: all groups run T0 (warp exits early if all die);
    # survivors are re-packed 32-dense and continue from T0 (warp max over survivors)
    ph1 = np.minimum(wm, T0).mean()
    surv = rows > T0
    rest = iters[surv] - T0
    rng.shuffle(rest)
    n = len(rest) // 32 * 32
    ph2 = rest[:n].reshape(-1, 32).max(1).sum() / WARPS if n else 0
    # restart variant: survivors recompute from row 0
    print(f"  compaction at T0={T0:2d}: phase1={ph1:.2f} phase2={ph2:.2f} total={ph1+ph2:.2f}  gain x{cur/(ph1+ph2):.2f}")
# multi-stage compaction every S rows
for S in [8, 12, 16]:
    tot = 0.0; live = iters.copy(); base = 0
    while len(live):
        rng.shuffle(live)
        m = (len(live) + 31) // 32
        pad = np.zeros(m * 32, np.int64); pad[:len(live)] = live
        tot += np.minimum(pad.reshape(m, 32).max(1) - base, S).clip(0).sum()
        base += S; live = live[live > base]
    print(f"  re-compaction every {S} rows: total={tot/WARPS:.2f} gain x{cur/(tot/WARPS):.2f}")
np.save(f"rows_K{K}.npy", rows)
# two compaction points (T0, T1): phase1 T0 rows on all, phase2 T0->T1 on 32-dense survivors, phase3 tail
for T0, T1 in [(20, 36), (24, 40), (16, 28), (20, 32), (24, 36)]:
    ph1 = T0
    s1_ = iters[iters > T0]; rng.shuffle(s1_)
    ph2 = len(s1_) / 32 * (T1 - T0) / WARPS
    s2_ = iters[iters > T1] - T1; rng.shuffle(s2_)
    n = len(s2_) // 32 * 32
    ph3 = s2_[:n].reshape(-1, 32).max(1).sum() / WARPS
    tot = ph1 + ph2 + ph3
    print(f"  two-point compaction T0={T0} T1={T1}: {ph1:.1f}+{ph2:.1f}+{ph3:.1f}={tot:.2f} gain x{cur/tot:.2f}")
print("words alive at row T:", {T: float((rows > T).mean()) for T in [24, 28, 32, 40]})

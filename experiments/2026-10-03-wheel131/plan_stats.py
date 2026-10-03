"""Pin-set classes and shift structure of the production plan (remote_plan_cost_f.txt)."""
import sys, collections
sys.path.insert(0, "tools")
import os; os.environ["MODCAP"] = "2e16"
import importlib.util
spec = importlib.util.spec_from_file_location("pu", "tools/plan_units.py"); pu = importlib.util.module_from_spec(spec); spec.loader.exec_module(pu)
plan = sys.argv[1] if len(sys.argv) > 1 else "experiments/remote_plan_cost_f.txt"
D0, NT = pu.D0, 58
def pins(K):
    d = K * D0; MOD = 30; P = []
    for q in pu.BAD:
        if d % q == 0: continue
        if MOD * q <= pu.MODCAP: MOD *= q; P.append(q)
        else: break
    return MOD, P
byK = collections.defaultdict(list); order = []
for line in open(plan):
    K, s = map(int, line.split()[:2]); byK[K].append(s); order.append((K, s))
cls = collections.Counter(); clsres = collections.Counter(); clspins = {}
for K, ss in byK.items():
    MOD, P = pins(K); res = 4
    for q in P: res *= q - NT
    cls[MOD] += len(ss); clsres[MOD] += res * len(ss); clspins[MOD] = P
tot = sum(clsres.values())
print(f"plan {plan}: units={len(order)} Ks={len(byK)} residues={tot:.4g}")
for MOD, n in sorted(cls.items(), key=lambda x: -clsres[x[0]])[:12]:
    print(f"  MOD={MOD:.5g} units={n} res_share={clsres[MOD]/tot:.3f} pins={clspins[MOD]} 131pinned={131 in clspins[MOD]}")
# shift structure
ns = sorted(len(v) for v in byK.values())
import statistics
print("shifts per K: median", statistics.median(ns), "mean %.1f" % statistics.mean(ns), "max", ns[-1])
for thr in (1, 10, 50, 131, 262):
    print(f"  Ks with >= {thr} shifts: {sum(1 for x in ns if x >= thr)}  units in them: {sum(x for x in ns if x >= thr)}")
mx = collections.Counter(max(v) for v in byK.values())
contig = sum(1 for v in byK.values() if sorted(v) == list(range(min(v), max(v)+1)))
print("Ks whose plan shifts are contiguous:", contig, "of", len(byK))
# fraction of units that would sit in COMPLETE blocks of 131 shifts [131j,131j+131) (for the 1.1337e15 class only)
full = 0; ref = 0
for K, ss in byK.items():
    MOD, P = pins(K)
    if 131 in P: continue
    ref += len(ss); st = set(ss)
    blocks = collections.Counter(s // 131 for s in ss)
    full += sum(c for b, c in blocks.items() if c == 131)
print(f"units in classes without 131 pinned: {ref}; of those in complete 131-shift blocks: {full} ({full/max(ref,1):.3f})")
# share of plan residues by 131 status
sh = collections.Counter()
for K, ss in byK.items():
    MOD, P = pins(K); res = 4
    for q in P: res *= q - NT
    st = "131|K" if K % 131 == 0 else ("pinned" if 131 in P else "tierC")
    sh[st] += res * len(ss)
print("residue share by 131 status:", {k: round(v / tot, 4) for k, v in sh.items()})

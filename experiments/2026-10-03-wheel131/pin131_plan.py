"""For main-class K (MOD = 1.1337e15, 131 not pinned) in the production plan, compare in the
planner's own model: E58 per residue of the plan's units for that K vs. a 131-pinned unit
(K, s'=0) [MOD' = 131*MOD], which covers the a-range of old shifts 0..130.
64-bit overflow is IGNORED here (upper bound on the pinned option), and the pinned unit is
given wbar = 0.5 (no window offset) as a further optimistic bound. Equal GPU cost/residue assumed."""
import os, sys, math, random, collections, importlib.util
os.environ["MODCAP"] = "2e16"
spec = importlib.util.spec_from_file_location("pu", "tools/plan_units.py"); pu = importlib.util.module_from_spec(spec); spec.loader.exec_module(pu)
D0 = pu.D0; CONST = pu.WCAL * 64 * pu.RHO_REF58
plan = sys.argv[1] if len(sys.argv) > 1 else "experiments/remote_plan_cost_f.txt"
byK = collections.defaultdict(list)
for line in open(plan):
    K, s = map(int, line.split()[:2]); byK[K].append(s)
random.seed(1)
Ks = list(byK); random.shuffle(Ks)
n = 0; Eo = Wo = En = Wn = 0.0; better = 0; ratios = []
for K in Ks:
    pu.MODCAP = 2e16; sh = pu.unit_shape(K)
    if sh is None or sh[0] != 1133661268029390: continue
    MOD, res, ly, wbar, lg = sh
    eo = sum(CONST * math.exp(ly + lg + pu.size_term(MOD, wbar, K, s)[0]) * res for s in byK[K])
    wo = res * len(byK[K])
    pu.MODCAP = 1.5e17; sh2 = pu.unit_shape(K)
    M2, r2, ly2, wb2, lg2 = sh2
    if M2 != MOD * 131: continue   # 131 | K: 131 already free
    en = CONST * math.exp(ly2 + lg2 + pu.size_term(M2, 0.5, K, 0)[0]) * r2
    Eo += eo; Wo += wo; En += en; Wn += r2
    ratios.append((en / r2) / (eo / wo)); better += (en / r2) > (eo / wo)
    n += 1
    if n >= 3000: break
ratios.sort()
print(f"sampled main-class K: {n}")
print(f"old plan units for these K: E58/res = {Eo/Wo:.4g}   (E58 {Eo:.4g} over {Wo:.4g} res)")
print(f"131-pinned unit s'=0 (optimistic): E58/res = {En/Wn:.4g}   (E58 {En:.4g} over {Wn:.4g} res)")
print(f"ratio new/old per residue: aggregate {(En/Wn)/(Eo/Wo):.3f}; per-K median {ratios[len(ratios)//2]:.3f}, "
      f"p90 {ratios[int(.9*len(ratios))]:.3f}, max {ratios[-1]:.3f}; K with new better: {better}")

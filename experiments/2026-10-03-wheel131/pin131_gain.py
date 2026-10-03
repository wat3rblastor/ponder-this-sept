"""Upper bound, in the planner's model, on the E58 gain from offering a 131-pinned unit
(K, s'=0; MOD' = 131*MOD) for EVERY main-class K of the production plan, at the plan's marginal
value per residue CUT (displaced work). Optimistic: ignores 64-bit overflow, wbar = 0.5,
equal GPU cost per residue.  gain_K = max(0, [E_new - CUT*W_new] - [E_old - CUT*W_old])."""
import os, sys, math, collections, importlib.util
os.environ["MODCAP"] = "2e16"
spec = importlib.util.spec_from_file_location("pu", "tools/plan_units.py"); pu = importlib.util.module_from_spec(spec); spec.loader.exec_module(pu)
D0 = pu.D0; CONST = pu.WCAL * 64 * pu.RHO_REF58
plan = sys.argv[1]; CUT = float(sys.argv[2])
COST = float(os.environ.get('PINCOST', '1'))   # GPU cost/residue of pinned unit vs reference
REALW = os.environ.get('REALW') == '1'           # use the engine's real window offset for MOD'
byK = collections.defaultdict(list)
for line in open(plan):
    K, s = map(int, line.split()[:2]); byK[K].append(s)
Etot = gain = Wextra = 0.0; nK = nwin = 0
for K, ss in byK.items():
    pu.MODCAP = 2e16; sh = pu.unit_shape(K)
    if sh is None: continue
    MOD, res, ly, wbar, lg = sh
    eo = sum(CONST * math.exp(ly + lg + pu.size_term(MOD, wbar, K, s)[0]) * res for s in ss)
    Etot += eo
    if MOD != 1133661268029390: continue
    pu.MODCAP = 1.5e17; M2, r2, ly2, wb2, lg2 = pu.unit_shape(K)
    if M2 != MOD * 131: continue
    nK += 1
    en = CONST * math.exp(ly2 + lg2 + pu.size_term(M2, wb2 if REALW else 0.5, K, 0)[0]) * r2 / COST
    g = (en - CUT * r2) - (eo - CUT * res * len(ss))
    if g > 0: gain += g; nwin += 1; Wextra += r2 - res * len(ss)
print(f"plan {plan}: model E58 of plan = {Etot:.4f}; main-class K (131 pinnable) = {nK}")
print(f"PINCOST={COST} REALW={REALW} cut {CUT:.3g}/res: K where pinned unit beats old units + displaced work: {nwin}; "
      f"upper-bound E58 gain = {gain:.4f} ({100*gain/Etot:.2f}% of plan); extra residues {Wextra:.3g}")

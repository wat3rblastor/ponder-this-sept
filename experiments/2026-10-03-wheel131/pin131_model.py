"""Planner-model value of pinning 131 for the K classes that do not pin it at MODCAP 2e16.

For each K, compare (in the production planner's own yield model, tools/plan_units.py):
  old : units (K, s), MOD = 1.1337e15 class, 4.67e9 residues each
  new : units (K, s'), MOD' = 131*MOD, 73x the residues, a-range 131x
Reports E58 per residue (log v) for both, and the 64-bit headroom of the new units.
"""
import os, sys, math, importlib.util
os.environ["MODCAP"] = "2e16"
spec = importlib.util.spec_from_file_location("pu", "tools/plan_units.py"); pu = importlib.util.module_from_spec(spec); spec.loader.exec_module(pu)
D0 = pu.D0
def E(K, s, modcap):
    pu.MODCAP = modcap
    sh = pu.unit_shape(K)
    MOD, res, ly, wbar, lg = sh
    st, T = pu.size_term(MOD, wbar, K, s)
    over = 64 * MOD * (s + 1) + 105 * MOD + 57 * K * D0 >= 1.8e19
    return MOD, res, ly + lg + st, T, over, wbar
CONST = pu.WCAL * 64 * pu.RHO_REF58
print("K      old: MOD res  v(s=0) v(s=7)  | new: MOD res v(s'=0)  T_last(s'=0)  64bit-overflow  amax=MOD'(1+c1+c2+64)")
for K in (7, 101, 1000, 3331, 10007, 50021, 100003, 200003, 400009):
    o0 = E(K, 0, 2e16); o7 = E(K, 7, 2e16); n0 = E(K, 0, 1.5e17)
    if 131 not in [q for q in pu.BAD if n0[0] % q == 0]: print(K, "131 not pinned in new"); continue
    amax = n0[0] * (1 + 73 + 55 + 64) + 57 * K * D0
    print(f"{K:7d} {o0[0]:.4g} {o0[1]:.3g} {o0[2]:+.3f} {o7[2]:+.3f} | {n0[0]:.4g} {n0[1]:.3g} {n0[2]:+.3f} {n0[3]:.3g} {n0[4]} {amax:.3g} {'>2^64' if amax >= 2**64 else ''}")

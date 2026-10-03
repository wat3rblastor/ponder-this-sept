#!/usr/bin/env python3
"""Cost model for the two realisations of the 17-automatic pencil F17.
All inputs are measured; nothing here is assumed.

  variant A  (M = 2*5*11*17*23*29 | d only)
     members(T) ~ 7.1e-7 T            (count.py, exact to 1e14)
     per-member success is NOT rho^41: 41, 47 and 53 are not in d, so every
     member has ~3.2 non-automatic indices in their forced residue classes
     where the term must be divisible by p^2.  Measured separately
     (poison.log): clean-index rate rc, poisoned-index rate rp ~ 0.011,
     P = rc^37.8 * rp^3.2, which is ~5e3 times smaller than rho^41.

  variant B  (M' = 2*5*11*17*23*29*41*47*53 | d)   <-- the better family
     members(T) ~ 5.38e-11 T          (variantB.log, exact to 1e16)
     no poisoned index, so P = rho_B^41 with rho_B measured directly.
"""
import math

D0 = 2*5*11*17*23*29

def fit(pairs):
    xs = [math.log(math.log(T)) for T, _ in pairs]
    ys = [math.log(r) for _, r in pairs]
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    s = -sum((x-mx)*(y-my) for x, y in zip(xs, ys))/sum((x-mx)**2 for x in xs)
    lnC = my + s*mx
    return (lambda T: math.exp(lnC)*math.log(T)**(-s)), math.exp(lnC), s

# ---- variant A -------------------------------------------------------------
C_A = 7.1e-7
RHO_A = [(1e10, 0.5565), (1e12, 0.5027), (1e14, 0.4647),
         (1e16, 0.4412), (1e18, 0.4130)]            # rates.log
RC_A = [(1e12, 0.5486), (1e16, 0.4760), (1e18, 0.4495)]   # poison.log (clean)
RP_A = 0.0115                                        # poison.log (poisoned)
MPOIS = 3.2
rho_A, _, _ = fit(RHO_A)
rc_A, _, _ = fit(RC_A)
def P_A(T): return rc_A(T)**(41-MPOIS) * RP_A**MPOIS

# ---- variant B -------------------------------------------------------------
C_B = 5.38e-11
RHO_B = [(1e14, 0.5052), (1e16, 0.4736), (1e18, 0.4485)]  # variantB.log
rho_B, CB, sB = fit(RHO_B)
def P_B(T): return rho_B(T)**41

def breakeven(c, P):
    lo, hi = 1e13, 1e45
    for _ in range(300):
        mid = math.sqrt(lo*hi)
        if c*mid*P(mid) < 1: lo = mid
        else: hi = mid
    return math.sqrt(lo*hi)

if __name__ == "__main__":
    print(f"variant B fit: rho_B(T) = {CB:.4f} (ln T)^(-{sB:.4f})")
    for T, r in RHO_B:
        print(f"   T=1e{round(math.log10(T))}: measured {r:.4f}  fit {rho_B(T):.4f}")
    print()
    print("  T     | A: members  P(member)   E[sols]  | B: members  P(member)   E[sols]")
    for e in range(14, 33, 2):
        T = 10.0**e
        ma, pa = C_A*T, P_A(T)
        mb, pb = C_B*T, P_B(T)
        print(f" 1e{e:<4d} | {ma:9.2e} {pa:10.2e} {ma*pa:9.2e}  | "
              f"{mb:9.2e} {pb:10.2e} {mb*pb:9.2e}")
    print()
    for nm, c, P in (("A", C_A, P_A), ("B", C_B, P_B)):
        T = breakeven(c, P)
        print(f"variant {nm}: break-even T* = {T:.2e}   members to test = {1/P(T):.2e}")
    print()
    print("pentagonal 13-automatic family under the same model "
          "(members 9.8e-9 T, 45 non-automatic, 41/47/53 not in d either):")
    Pp = lambda T: rc_A(T)**(45-MPOIS) * RP_A**MPOIS
    T = breakeven(9.8e-9, Pp)
    print(f"   E(1e20) = {9.8e-9*1e20*Pp(1e20):.2e}   break-even T* = {T:.2e}"
          f"   members to test = {1/Pp(T):.2e}")
    print()
    print("same-work ratio against the unconstrained search at equal T")
    print("  (unconstrained: ~T^2/(114*D0) candidates, 58 terms each)")
    for e in (14, 18, 22, 28):
        T = 10.0**e; r = rho_B(T)
        print(f"   T=1e{e:<3d}  variantB/unconstrained = "
              f"{(C_B*T*r**41)/((T*T/(114*D0))*r**58):.2e}")

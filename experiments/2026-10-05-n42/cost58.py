"""n=58 cost under the residual measured in this experiment.

Two distinct questions, which this repo has mixed before:

(A) SUPPLY -- how many 58-term APs EXIST below T.  Family d = 3741870*m.
    Smallest bad prime not dividing d is q_min = 41, so this question sits at
    s = 58/82 = 0.707 and IS exposed to any s-drift.

(B) COST -- core-hours for the production engine to find one.  Family
    d = K*D0, D0 = 2*3*5*11*17*23*29*41*47*53, with the production sieve
    b2 = 2000 which guarantees NO bad prime <= 2000 divides ANY of the 58
    terms.  Every bad prime that could produce a hit is therefore > 2000 >> 58,
    so the engine family has NO multi-hit / class-exhaustion structure at all:
    its effective s is 58/(2*2003) = 0.014, not 0.707.  The s-drift cannot
    reach it.  (This is the point this experiment adds to forecast.py.)

Residuals are read from the command line: R_A (model/exact at s ~ 0.707) and
R_B (model/exact at low s), so cost = COST_C * T*^2 with E_model(T*) = R*ln2.
"""
import math
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), '2026-10-04-groundtruth'))
import model as M
import forecast as F


_CACHE = {}


def solve(fn, target, tag):
    """E(T) is a smooth near-power-law in T, so sample ln E on a coarse log
    grid of T (cached) and invert by local log-log interpolation.  A 70-step
    bisection on fn would cost 70 full term integrals, which is minutes."""
    if tag not in _CACHE:
        grid = [10 ** e for e in (16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20)]
        _CACHE[tag] = [(math.log(T), math.log(max(fn(T), 1e-300))) for T in grid]
    pts = _CACHE[tag]
    y = math.log(target)
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if y0 <= y <= y1:
            return math.exp(x0 + (x1 - x0) * (y - y0) / (y1 - y0))
    return float('nan')


def main():
    RA = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    RB = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    band = [float(v) for v in sys.argv[3:5]] if len(sys.argv) > 4 else [RB, RB]
    print(f"residual applied: R_A(s=0.707) = {RA:.3f}, R_B(low s) = {RB:.3f}, "
          f"band on R_B = {band[0]:.3f}..{band[1]:.3f}")
    print("exact truth = model / R, so E_true(T) = E_model(T)/R.\n")

    print("(A) SUPPLY: expected # of 58-term APs (pairs) with all terms <= T, "
          "forced family d=3741870*m")
    print(f"{'T':>8} {'E model':>11} {'E corrected':>12}")
    for T in (1e16, 1e17, 1e18, 1e19):
        E = F.famA(T)
        print(f"{T:8.0e} {E:11.4e} {E/RA:12.4e}")
    TA = solve(F.famA, math.log(2) * RA, 'A')
    print(f"    even-odds existence (E_true >= ln2): T = {TA:.3e}  "
          f"(model alone: {solve(F.famA, math.log(2), chr(65)):.3e})")

    print("\n(B) COST: engine family, cost = 4.69e-32 T^2 core-hours")
    print(f"{'T':>8} {'E model':>11} {'E corrected':>12} {'core-hours':>12} "
          f"{'core-h per 58':>14}")
    for T in (1e17, 1e18, 1e19):
        E0 = F.famB(T)[0]
        E = E0 / RB
        c = F.COST_C * T * T
        print(f"{T:8.0e} {E0:11.4e} {E:12.4e} {c:12.4e} {c/E:14.4e}")
    out = []
    for R in (RB, band[0], band[1]):
        T = solve(lambda T: F.famB(T)[0], math.log(2) * R, 'B')
        out.append((R, T, F.COST_C * T * T))
    print(f"    even-odds find: T* = {out[0][1]:.3e}, "
          f"cost = {out[0][2]:.3e} core-hours")
    print(f"    band from R_B = {band[0]:.3f}..{band[1]:.3f}: "
          f"T* = {out[1][1]:.3e}..{out[2][1]:.3e}, "
          f"cost = {out[1][2]:.3e}..{out[2][2]:.3e} core-hours")
    print(f"    (groundtruth's figure: 1.5e5-1.77e5; plan's: ~2.0e5)")
    print(f"\n    term-size band a production run should target: "
          f"last term near T* = {out[0][1]:.2e}")


if __name__ == '__main__':
    main()

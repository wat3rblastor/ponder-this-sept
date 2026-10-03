#!/usr/bin/env python3
"""Patched copy of tools/plan_units.py (mathres, 2026-10-03). Same CLI and output format.

Changes to the score (all measured, see /workspace/mathres notes in the report):
  1. Size term: per-term pass probability rho(t) = 0.71992*(ln t/39.144)^-0.485 (Monte Carlo,
     1.2e7 samples at each of 9 sizes 1e15..3e18), applied to the actual terms a+k*d averaged
     over the unit's 64 window multiples, instead of (ln T_last)^-29.
  2. Window position: the strided kernel covers multiples [64s+w, 64s+64+w) with
     w = floor((R + i1*s1 + i2*s2)/MOD); mean w = ((c1-1)f1+(c2-1)f2)/2 with f_j = sstep_j/MOD
     computed exactly per K (the old score assumed w = 0).
  3. Good primes p = 1 (mod 3) dividing K: 1/p of the windows (p | a) are rescalings of unit
     K/p at lower a (already covered / ranked earlier), and the remaining windows have no term
     divisible by p, which lowers rho. Factor (1-1/p)*exp(-(58/(p-1))*0.485*ln p/ln T).
Absolute calibration: windows per unit = 0.752*64*res*exp(ly) (measured, constant across 12k units).
"""
import argparse, glob, json, math, sys, heapq, os

D0 = 382160924970
NT = 58
MODCAP = int(float(os.environ.get("MODCAP", "2e15")))
EXPO = 0.485
WCAL = 0.752

def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

PR = sieve(10000)
BAD = [p for p in PR if p % 3 == 2 and p > NT]
GOOD = [p for p in PR if p % 3 == 1 and p < 1000]
LOGF = {q: math.log((q - NT) / q) for q in BAD}
BASE_UNPINNED = sum(LOGF.values())

def unit_shape(K):
    """(MOD, res, ly, wbar, lgood) ; ly excludes the size term."""
    d = K * D0
    MOD, res, pinned, divs = 30, 4, [], []
    full = False
    for q in BAD:
        if d % q == 0:
            divs.append(q); continue
        if not full and MOD * q <= MODCAP:
            MOD *= q; res *= (q - NT); pinned.append(q)
        else:
            full = True
            if q > 600 and K < q: break
    if len(pinned) < 4: return None
    ly = BASE_UNPINNED
    for q in pinned: ly -= LOGF[q]
    for q in BAD:
        if q > K: break
        if K % q == 0: ly += math.log((q - 1) / q) - LOGF[q]
    # mean window offset of the strided walk: two largest-count components
    q2, q1 = pinned[-1], pinned[-2]
    WM = os.environ.get("WMODE", "0")
    wbar = 0.5
    for q in (q1, q2):
        co = MOD // q
        f = (d % q) * pow(co % q, -1, q) % q / q
        if WM in ("1","2"): f = min(f, 1 - f)
        if WM == "2" and q == q1: f = 0
        if WM == "3": f = 0
        wbar += (q - NT - 1) * f / 2
    lgood = 0.0
    for p in GOOD:
        if p > K: break
        if K % p == 0:
            lgood += math.log(1 - 1 / p) - (58 / (p - 1)) * EXPO * math.log(p) / 40.5
    return MOD, res, ly, wbar, lgood

BS = [4, 12, 20, 28, 36, 44, 52, 60]
KS = [0, 11.4, 22.8, 34.2, 45.6, 57]
def size_term(MOD, wbar, K, s):
    """ln of mean over window multiples of prod_k rho(t_k)/rho_ref^58, and mean last term."""
    d = K * D0
    acc = 0.0
    for b in BS:
        a = (64 * s + b + wbar) * MOD
        m = 0.0
        for k in KS: m += math.log(math.log(a + k * d + d))
        acc += math.exp(-EXPO * NT * (m / len(KS) - math.log(39.144)))
    return math.log(acc / len(BS)), (64 * s + 32 + wbar) * MOD + 57 * d

RHO_REF58 = 0.71992 ** NT

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget-res", type=float, default=7e14)
    ap.add_argument("--kmax", type=int, default=200000)
    ap.add_argument("--smax", type=int, default=64)
    ap.add_argument("--out", default=None)
    ap.add_argument("--done-glob", default="experiments/*/*.jsonl")
    ap.add_argument("--old-plane0-kmax", type=int, default=1949)
    ap.add_argument("--k-keep-mod", type=int, default=0)
    ap.add_argument("--k-drop-mod", type=int, default=0)
    ap.add_argument("--score-plan", default=None, help="evaluate an existing 'K shift' plan file under this model")
    ap.add_argument("--marks", default="1e15,3e15,6e15,1e16,1.5e16,2e16")
    args = ap.parse_args()

    done = set()
    for f in glob.glob(args.done_glob):
        for line in open(f):
            if '"covered"' in line:
                try: j = json.loads(line)
                except Exception: continue
                done.add((j["K"], j["shift"]))
    marks = [float(x) for x in args.marks.split(",")]

    shapes = {}
    def E_unit(K, s):
        """(E58 of the unit, residues) or None"""
        if K not in shapes: shapes[K] = unit_shape(K)
        sh = shapes[K]
        if sh is None: return None
        MOD, res, ly, wbar, lgood = sh
        if 64 * MOD * (s + 1) + 105 * MOD + 57 * K * D0 >= 1.8e19: return None
        st, T = size_term(MOD, wbar, K, s)
        v = ly + lgood + st                     # log yield per residue (up to constants)
        return v, res

    CONST = WCAL * 64 * RHO_REF58               # E58 per residue = CONST * exp(v)

    if args.score_plan:
        spent = E = 0.0; mi = 0; n = 0
        for line in open(args.score_plan):
            K, s = map(int, line.split()[:2])
            if (K, s) in done: continue
            r = E_unit(K, s)
            if r is None: continue
            v, res = r
            spent += res; E += CONST * math.exp(v) * res; n += 1
            while mi < len(marks) and spent >= marks[mi]:
                print(f"  at {marks[mi]:.3g} res: E58={E:.3f} (units {n}, marginal/res {CONST*math.exp(v):.3g})", file=sys.stderr); mi += 1
        print(f"plan {args.score_plan}: units={n} residues={spent:.3g} E58={E:.3f}", file=sys.stderr)
        return

    heap = []
    for K in range(1, args.kmax + 1):
        r = E_unit(K, 0)
        if r is None: continue
        heapq.heappush(heap, (-r[0], K, 0))
    plan, spent, E, mi = [], 0.0, 0.0, 0
    while heap and spent < args.budget_res:
        nv, K, s = heapq.heappop(heap)
        if s + 1 < args.smax:
            r2 = E_unit(K, s + 1)
            if r2 is not None: heapq.heappush(heap, (-r2[0], K, s + 1))
        if (K, s) in done or (s == 0 and K <= args.old_plane0_kmax): continue
        if args.k_keep_mod and K % args.k_keep_mod: continue
        if args.k_drop_mod and K % args.k_drop_mod == 0: continue
        res = shapes[K][1]
        plan.append((K, s))
        spent += res
        E += CONST * math.exp(-nv) * res
        while mi < len(marks) and spent >= marks[mi]:
            print(f"  at {marks[mi]:.3g} res: E58={E:.3f} (units {len(plan)}, marginal/res {CONST*math.exp(-nv):.3g})", file=sys.stderr); mi += 1
    print(f"units={len(plan)} residues={spent:.3g} E58={E:.3f} maxK={max(p[0] for p in plan)} "
          f"maxshift={max(p[1] for p in plan)}", file=sys.stderr)
    if args.out:
        with open(args.out, "w") as f:
            for K, s in plan: f.write(f"{K} {s}\n")

if __name__ == "__main__":
    main()

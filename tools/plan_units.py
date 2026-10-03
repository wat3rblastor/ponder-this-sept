#!/usr/bin/env python3
"""Rank (K, shift) work units by expected 58-term hits per residue.

For d = K*D0 the hit probability of an enumerated residue factors as
    prod over bad q > 58:   q | K      -> (q-1)/q        (only a != 0 mod q needed)
                            pinned     -> 1              (enumerated, cost-free)
                            otherwise  -> (q-58)/q       (tier-C filter, avoid-only)
times a large-prime density term ~ (ln T)^(-nterms/2), T = last term.
Dividing K by a bad prime q therefore multiplies the per-residue yield by
(q-1)/(q-58) relative to a generic K (and frees a pin slot if q was pinned),
at the price of whatever it does to T. Units are sorted by yield per residue.

Usage: plan_units.py --budget-res 7e14 [--kmax 200000] [--smax 60] --out plan.txt
Prints the model's summed expectation, calibrated so a generic K~1000, shift-0
unit yields E58 = 1e-5 (measured, PROGRESS.md 2026-10-03).
"""
import argparse, glob, json, math, sys, heapq

D0 = 382160924970
NT = 58
MODCAP = int(float(__import__("os").environ.get("MODCAP", "2e15")))

def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

BAD = [p for p in sieve(10000) if p % 3 == 2 and p > NT]
LOGF = {q: math.log((q - NT) / q) for q in BAD}            # unpinned, not dividing K
BASE_UNPINNED = sum(LOGF.values())

def unit_shape(K):
    """(MOD, residues, log yield factor excluding the T term) for this K."""
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
    if len(pinned) < 4: return None
    ly = BASE_UNPINNED
    for q in pinned: ly -= LOGF[q]
    for q in divs: ly += math.log((q - 1) / q) - LOGF[q]
    # large bad prime factors of K beyond 10000 change nothing measurable
    return MOD, res, ly

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget-res", type=float, default=7e14)
    ap.add_argument("--kmax", type=int, default=200000)
    ap.add_argument("--smax", type=int, default=64)
    ap.add_argument("--out", default=None)
    ap.add_argument("--done-glob", default="experiments/2026-10-03-cuda/*.jsonl")
    ap.add_argument("--old-plane0-kmax", type=int, default=1949,
                    help="shift 0 is covered for K <= this by the reduced-kernel campaigns")
    args = ap.parse_args()

    done = set()
    for f in glob.glob(args.done_glob):
        for line in open(f):
            if '"covered"' in line:
                j = json.loads(line); done.add((j["K"], j["shift"]))

    ref = unit_shape(1009)                       # generic K: nothing bad divides it
    refT = 64 * ref[0] * 0.5 + 57 * 1009 * D0
    def val(shape, K, s):                        # log yield per residue, relative
        MOD, res, ly = shape
        T = 64 * MOD * (s + 0.5) + 57 * K * D0
        if 64 * MOD * (s + 1) + 105 * MOD + 57 * K * D0 >= 1.8e19: return None, None
        return ly - (NT / 2) * math.log(math.log(T)), T
    refv, _ = val(ref, 1009, 0)
    E_REF = 1e-5                                 # E58 of the reference unit
    ref_per_res = E_REF / ref[1]

    # lazy best-first over shifts: value is decreasing in s for fixed K
    heap = []
    shapes = {}
    for K in range(1, args.kmax + 1):
        sh = unit_shape(K)
        if sh is None: continue
        shapes[K] = sh
        v, _ = val(sh, K, 0)
        if v is not None: heapq.heappush(heap, (-v, K, 0))
    plan, spent, E, Elevel = [], 0.0, 0.0, 0.0
    while heap and spent < args.budget_res:
        nv, K, s = heapq.heappop(heap)
        if s + 1 < args.smax:
            v2, _ = val(shapes[K], K, s + 1)
            if v2 is not None: heapq.heappush(heap, (-v2, K, s + 1))
        if (K, s) in done or (s == 0 and K <= args.old_plane0_kmax): continue
        res = shapes[K][1]
        plan.append((K, s, -nv - refv))
        spent += res
        E += ref_per_res * math.exp(-nv - refv) * res
    print(f"units={len(plan)} residues={spent:.3g} E58={E:.3f} "
          f"marginal rel. yield={math.exp(plan[-1][2]):.3f} maxK={max(p[0] for p in plan)} "
          f"maxshift={max(p[1] for p in plan)}", file=sys.stderr)
    if args.out:
        with open(args.out, "w") as f:
            for K, s, v in plan: f.write(f"{K} {s}\n")

if __name__ == "__main__":
    main()

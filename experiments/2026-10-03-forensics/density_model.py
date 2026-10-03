#!/usr/bin/env python3
"""Is the exhaustive count of long runs per step d fully explained by the
independent-local-densities model?

For a start a chosen at random and a window of n terms t_k = a+k*d, the
probability that no bad prime q spoils the window factors (heuristically) over
q.  For q not dividing d, let A = #{k0 mod q : no index in [0,n) hits q} and
B = #{k0 : exactly one index hits}.  Two or more hits is fatal (it would force
q | d), one hit survives with probability 1/q (needs q^2 | t).  So

    P_q(q !| d) = (A + B/q)/q,        P_q(q | d) ~ 1 - 1/q.

predicted(d, n)  =  (#starts in range) * prod_{bad q <= QCUT} P_q(d, n)

The q > QCUT tail is almost the same for every d in a corpus (it depends on d
only through the term size), so observed/predicted should be a CONSTANT across
d if the model is complete.  Scatter in that ratio is the only room left for
hidden structure in the choice of d.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from loeschian import bad_primes  # noqa: E402

QCUT = 300


def pq(q, n, divides):
    if divides:
        return 1.0 - 1.0 / q
    A = 0
    B = 0
    for k0 in range(q):
        h = (n - 1 - k0) // q + 1 if k0 < n else 0
        if h == 0:
            A += 1
        elif h == 1:
            B += 1
    return (A + B / q) / q


def predicted(d, n, N):
    R = N - (n - 1) * d
    if R <= 0:
        return 0.0
    logp = math.log(R)
    for q in bad_primes(QCUT):
        p = pq(q, n, d % q == 0)
        if p <= 0:
            return 0.0
        logp += math.log(p)
    return math.exp(logp)


def main(paths, ks):
    for path in paths:
        rows = [json.loads(ln) for ln in Path(path).read_text().splitlines()
                if ln.strip() and '"g1"' not in ln]
        print(f"\n===== {Path(path).name}: {len(rows)} exhaustive steps")
        for k in ks:
            data = []
            for r in rows:
                c = r["ge"].get(str(k), 0)
                pr = predicted(r["d"], k, r["N"])
                if pr <= 0:
                    continue
                data.append((r["d"], r["m"], c, pr))
            if len(data) < 10:
                print(f"  k={k}: {len(data)} usable steps -- skip")
                continue
            data_nz = [z for z in data if z[2] > 0]
            lr = [math.log(c / pr) for _, _, c, pr in data_nz]
            mu = sum(lr) / len(lr)
            sd = math.sqrt(sum((x - mu) ** 2 for x in lr) / (len(lr) - 1))
            print(f"  k={k}: {len(data)} steps, log(obs/pred) mean={mu:+.3f} "
                  f"sd={sd:.3f}  => obs/pred spread = "
                  f"x{math.exp(sd):.2f} (1 sd)")
            rng = sorted(data_nz, key=lambda z: z[2] / z[3])
            print(f"      worst  d={rng[0][0]} m={rng[0][1]} obs={rng[0][2]} "
                  f"pred={rng[0][3]:.1f} ratio={rng[0][2]/rng[0][3]:.2f}")
            print(f"      best   d={rng[-1][0]} m={rng[-1][1]} obs={rng[-1][2]} "
                  f"pred={rng[-1][3]:.1f} ratio={rng[-1][2]/rng[-1][3]:.2f}")
            # does the residual correlate with an extra bad prime in d?
            # censoring-free Poisson group comparison: sum(obs)/sum(pred)
            for q in (23, 29, 41, 47, 53, 59, 71):
                w = [z for z in data if z[0] % q == 0]
                wo = [z for z in data if z[0] % q != 0]
                if len(w) >= 3 and len(wo) >= 3:
                    rw = sum(z[2] for z in w) / sum(z[3] for z in w)
                    ro = sum(z[2] for z in wo) / sum(z[3] for z in wo)
                    nw = sum(z[2] for z in w)
                    no = sum(z[2] for z in wo)
                    rse = math.sqrt(1 / max(1, nw) + 1 / max(1, no))
                    if rw <= 0 or ro <= 0:
                        continue
                    print(f"      resid q={q}: obs/pred with={rw:.3e} "
                          f"({len(w)} steps, {nw} runs) without={ro:.3e} "
                          f"({len(wo)} steps, {no} runs)  model-corrected "
                          f"excess x{rw/ro:.2f} (logse {rse:.2f})"
                          f"{'  ***' if abs(math.log(rw/ro)) > 2.5 * rse else ''}")
            # raw prime effect on the PREDICTED count, for reference
            for q in (23, 29, 41, 47, 53, 59, 71):
                print(f"      model P_q(k={k}, q={q}): not|d {pq(q, k, False):.4f} "
                      f"| d {pq(q, k, True):.4f}  gain "
                      f"x{pq(q, k, True)/pq(q, k, False):.1f}")


if __name__ == "__main__":
    main(sys.argv[1:-1], [int(x) for x in sys.argv[-1].split(",")])

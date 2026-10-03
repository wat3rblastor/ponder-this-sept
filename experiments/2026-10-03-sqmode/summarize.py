#!/usr/bin/env python3
"""Per-unit-of-work yield from the .jsonl unit records of a matched run."""
import json
import sys

for path in sys.argv[1:]:
    tot = dict(secs=0.0, res=0, surv=0, conf=0, ge20=0, ge25=0, ge30=0,
               ge35=0, ge40=0)
    units = 0
    minlast = 0
    best = 0
    ind = sq = ""
    for line in open(path):
        r = json.loads(line)
        if r.get("progress") or r.get("hit"):
            continue
        units += 1
        ind, sq = r["ind"], r["sq"]
        for k in tot:
            tot[k] += r[k]
        if r["minlast"] and (not minlast or r["minlast"] < minlast):
            minlast = r["minlast"]
        best = max(best, r["best"])
    s = tot["secs"] or 1
    print(f"{path}")
    print(f"  mode ind=[{ind}] sq=[{sq}]  units={units}  secs={s:.1f}")
    print(f"  residues={tot['res']:.4g} ({tot['res']/s:.4g}/s)  "
          f"survivors={tot['surv']} ({tot['surv']/s:.4g}/s)")
    for k in ("ge20", "ge25", "ge30", "ge35", "ge40"):
        print(f"  {k}: {tot[k]:>8}   per 1000s: {1000*tot[k]/s:9.3f}   "
              f"per 1e6 survivors: {1e6*tot[k]/max(1,tot['surv']):9.2f}")
    print(f"  best run = {best}   smallest last term (runs >= 30) = {minlast}")

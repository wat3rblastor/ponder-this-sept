#!/usr/bin/env python3
"""Collect every distinct progression of 55+ terms found by any campaign.

Scans experiments/*/*.log and *.jsonl for hits with n >= 55, reduces each to its
primitive form (a/g, d/g with g = gcd(a, d); a rescaled copy m*(a, d) of a known
progression is the same progression and is listed once), verifies each with
src/verify.py and src/crosscheck.py, stores the list in records.json under
"g3_runs_55plus" and regenerates ANSWER.md. Exit status 0 if the list changed,
3 if it did not.
"""
import glob, json, re, subprocess, sys
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import promote  # noqa: E402

NMIN = 55


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    return p.returncode, p.stdout + p.stderr


def main():
    cands = {}
    pat = re.compile(r"n=(\d+) a=(\d+) d=(\d+)")
    for f in glob.glob(str(ROOT / "experiments/*/*.log")):
        for line in open(f, errors="replace"):
            if "*** n=" in line:
                m = pat.search(line)
                if m and int(m.group(1)) >= NMIN:
                    cands[(int(m.group(2)), int(m.group(3)))] = int(m.group(1))
    for f in glob.glob(str(ROOT / "experiments/*/*.jsonl")):
        for line in open(f, errors="replace"):
            if '"hit"' in line:
                try: j = json.loads(line)
                except ValueError: continue
                if j.get("n", 0) >= NMIN:
                    cands[(j["a"], j["d"])] = j["n"]
    recs = json.loads((ROOT / "records.json").read_text())
    g3 = recs["g3_longest"]
    cands[(g3["a"], g3["d"])] = g3["n"]

    prim = {}
    for (a, d), n in cands.items():
        g = gcd(a, d)
        key = (a // g, d // g)
        prim[key] = max(prim.get(key, 0), n)

    old = {(e["a"], e["d"]): e for e in recs.get("g3_runs_55plus", [])}
    out = []
    for (a, d), n in sorted(prim.items(), key=lambda kv: (-kv[1], kv[0][0] + (kv[1] - 1) * kv[0][1])):
        e = old.get((a, d))
        if e and e["n"] >= n:
            out.append(e); continue
        rc1, o1 = run(["python3", "src/verify.py", str(a), str(d), str(n), "--quiet", "--maximal"])
        rc2, o2 = run(["python3", "src/crosscheck.py", str(a), str(d), str(n)])
        if rc1 or "OVERALL: PASS" not in o1 or rc2 or "CROSS-CHECK: PASS" not in o2:
            print(f"SKIP (failed verification): a={a} d={d} n={n}"); continue
        out.append({"a": a, "d": d, "n": n, "last": a + (n - 1) * d,
                    "d_factored": promote.fmt_factored(d), "verified": True,
                    "maximal_both_ends": "maximal in both directions" in o1})
        print(f"listed n={n} a={a} d={d}")
    # rescaled copies exactly as the search found them (g*a, g*d), listed separately
    oldc = {(e["a"], e["d"]): e for e in recs.get("g3_copies_55plus", [])}
    listed = {(e["a"], e["d"]) for e in out}
    copies = []
    for (a, d), n in sorted(cands.items(), key=lambda kv: (-kv[1], kv[0][0] + (kv[1] - 1) * kv[0][1])):
        g = gcd(a, d)
        if g == 1 or (a // g, d // g) not in listed:
            continue
        e = oldc.get((a, d))
        if e and e["n"] >= n:
            copies.append(e); continue
        rc1, o1 = run(["python3", "src/verify.py", str(a), str(d), str(n), "--quiet"])
        rc2, o2 = run(["python3", "src/crosscheck.py", str(a), str(d), str(n)])
        if rc1 or "OVERALL: PASS" not in o1 or rc2 or "CROSS-CHECK: PASS" not in o2:
            print(f"SKIP copy (failed verification): a={a} d={d} n={n}"); continue
        copies.append({"a": a, "d": d, "n": n, "last": a + (n - 1) * d, "multiplier": g,
                       "primitive_a": a // g, "primitive_d": d // g, "verified": True})
        print(f"listed copy n={n} a={a} d={d} (x{g})")
    if out == recs.get("g3_runs_55plus") and copies == recs.get("g3_copies_55plus", []):
        return 3
    recs["g3_runs_55plus"] = out
    recs["g3_copies_55plus"] = copies
    (ROOT / "records.json").write_text(json.dumps(recs, indent=2) + "\n")
    promote.write_answer(recs)
    print(f"{len(out)} distinct progressions of {NMIN}+ terms; ANSWER.md regenerated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

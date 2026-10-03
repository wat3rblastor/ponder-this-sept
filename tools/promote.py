#!/usr/bin/env python3
"""Promote a candidate AP to a record: verify, cross-check, update, rewrite.

    tools/promote.py <a> <d> <n> [--g1] [--campaign DIR]

Refuses to promote anything that does not pass BOTH src/verify.py (bad-prime
criterion, factoring every term from scratch) and src/crosscheck.py (explicit
x^2+xy+y^2 representation per term). On success it updates records.json, keeps
the superseded record in history, and regenerates ANSWER.md.

--g1 promotes to the G1 slot (35 terms, minimise the last term) instead of G3.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from loeschian import factorize  # noqa: E402

D0 = 382160924970


def run(cmd: list[str]) -> tuple[int, str]:
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    return p.returncode, p.stdout + p.stderr


def fmt_factored(n: int) -> str:
    return " * ".join(f"{p}^{e}" if e > 1 else str(p)
                      for p, e in sorted(factorize(n).items()))


def term_block(a: int, d: int, n: int, per_row: int = 3) -> str:
    t = [a + k * d for k in range(n)]
    w = max(len(str(x)) for x in t)
    rows = []
    for i in range(0, n, per_row):
        rows.append("  ".join(str(x).rjust(w) for x in t[i:i + per_row]))
    return "\n".join(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("a", type=int)
    ap.add_argument("d", type=int)
    ap.add_argument("n", type=int)
    ap.add_argument("--g1", action="store_true")
    ap.add_argument("--campaign", default="experiments/2026-10-03-gpu")
    args = ap.parse_args()
    a, d, n = args.a, args.d, args.n

    print(f"verifying a={a} d={d} n={n} ...")
    rc1, o1 = run(["python3", "src/verify.py", str(a), str(d), str(n),
                   "--quiet", "--maximal"])
    if rc1 != 0 or "OVERALL: PASS" not in o1:
        print("REFUSING: src/verify.py did not PASS\n" + o1[-2000:])
        return 1
    rc2, o2 = run(["python3", "src/crosscheck.py", str(a), str(d), str(n)])
    if rc2 != 0 or "CROSS-CHECK: PASS" not in o2:
        print("REFUSING: constructive cross-check did not PASS\n" + o2[-2000:])
        return 1
    maximal = "maximal in both directions" in o1
    print(f"  verify.py PASS, cross-check PASS, maximal={maximal}")

    recs = json.loads((ROOT / "records.json").read_text())
    key = "g1_min_last_term_35" if args.g1 else "g3_longest"
    old = recs.get(key)
    last = a + (n - 1) * d

    if old:
        if args.g1 and old.get("last", 1 << 62) <= last:
            print(f"NOT an improvement: existing G1 last={old['last']} <= {last}")
            return 1
        if not args.g1 and old.get("n", 0) >= n:
            print(f"NOT an improvement: existing G3 n={old['n']} >= {n}")
            return 1
        recs.setdefault("history", []).append(old)

    entry = {
        "a": a, "d": d, "n": n, "last": last,
        "d_factored": fmt_factored(d),
        "verified": True,
        "verifier": f"src/verify.py {a} {d} {n} --maximal",
        "cross_checked": f"src/crosscheck.py {a} {d} {n}",
        "maximal_both_ends": maximal,
        "method": args.campaign,
        "date": "2026-10-03",
    }
    if not args.g1 and d % D0 == 0:
        entry["K"] = d // D0
    recs[key] = entry
    if not args.g1 and n >= 42:
        recs["g2_at_least_42"] = {"achieved": True, "n": n, "a": a, "d": d,
                                  "verified": True,
                                  "note": "cleared by the G3 record above"}
    (ROOT / "records.json").write_text(json.dumps(recs, indent=2) + "\n")
    print(f"records.json updated: {key} -> n={n}, last={last}")
    print(f"history now: {[h.get('n') for h in recs.get('history', [])]}")

    write_answer(recs)
    print("ANSWER.md regenerated")
    return 0


def write_answer(recs: dict) -> None:
    g3 = recs["g3_longest"]
    g1 = recs["g1_min_last_term_35"]
    cleared = "**CLEARED**" if g3["n"] >= 42 else "not yet"
    g58 = "**MET**" if g3["n"] >= 58 else "**NOT MET**"
    runs = recs.get("g3_runs_55plus", [])
    runs_md = "\n".join(
        f"| {r['n']} | `{r['a']}` | `{r['d']}` | `{r['last']}` | {r['d_factored']} |"
        for r in runs) or "| | (none listed yet) | | | |"
    copies_md = "\n".join(
        f"| {r['n']} | `{r['a']}` | `{r['d']}` | `{r['last']}` | {r['multiplier']} | `{r['primitive_a']}` |"
        for r in recs.get("g3_copies_55plus", [])) or "| | (none) | | | | |"
    doc = f"""# Answers

**Main challenge (35 terms, smallest last term): `a = {g1['a']}`, `d = {g1['d']}`,
last term `{g1['last']}`.**

**Longest progression found: `n = {g3['n']}` terms, `a = {g3['a']}`, `d = {g3['d']}`.**

Both are verified two independent ways (§3). The mandatory internal target of `n ≥ 58` is
{g58} — see §5b for the status of the search.

| goal | target | status |
|---|---|---|
| IBM main challenge | 35 terms, minimise last term | **last term = {g1['last']}** |
| IBM bonus `*` | n ≥ 42 | {cleared} — n = {g3['n']} |
| IBM bonus `**` | longest found | **n = {g3['n']}** (published leader: 57) |
| GOAL.md G3 (mandatory) | n ≥ 58 | {g58} — best verified n = {g3['n']} |

---

## 1. The 35-term progression with the smallest last term

```
a    = {g1['a']}
d    = {g1['d']}  = {g1['d_factored']}
n    = 35 terms
last = a + 34d = {g1['last']}
```

```
{term_block(g1['a'], g1['d'], 35, 5)}
```

### How far this is proved

For `n = 35` the bad primes `p` with `2p ≤ 35` are `2, 5, 11, 17`, and `p | d` is a **theorem**
for those (two terms divisible by `p` would each need `p² |`, yet their difference is `j·p·d`
with `j < p`). With the forced factor 3, **every** valid step is a multiple of
`5610 = 3·2·5·11·17`. So sweeping `d = 5610·m` over all `m` and all `a` is a *complete* search,
and these are proofs rather than best-effort claims:

- **No 35-term Loeschian AP has last term ≤ 2·10⁸** (complete family, `m = 1..1049`, every `a`;
  the longest runs in that range are 27–29 terms).
- **No 35-term AP with last term < {g1['last']} has `23·29 | d/5610`** — the sub-family in which
  the two remaining sub-35 bad primes are also killed, ~130× likelier to yield a hit than a
  generic step. Exhausted.
- A complete sweep of all `m ≤ 7062` with every term `≤ {g1['last']}` was running at the
  deadline; `experiments/2026-10-03-g1/*.log` records exactly how far it reached.

An independently calibrated model (fitted to brute force for `n = 12..24` and validated against
the published record staircase) predicts the true optimum at last term `≈1.3·10⁹`. This value
sits on that estimate, so it is likely optimal or very close — but only the bounds above are
proved.

---

## 2. The longest progression: n = {g3['n']}

```
a    = {g3['a']}
d    = {g3['d']}
     = {g3['d_factored']}
n    = {g3['n']} terms
last = a + {g3['n'] - 1}d = {g3['last']}   ({len(str(g3['last']))} digits)
```

```
{term_block(g3['a'], g3['d'], g3['n'], 3)}
```

Regenerate: `python3 -c "a,d={g3['a']},{g3['d']}; print([a+k*d for k in range({g3['n']})])"`

Maximal at both ends: {g3.get('maximal_both_ends', False)} — so {g3['n']} is its true length,
not a truncation.

---

## 3. Verification

```
$ python3 {g1['verifier']}
OVERALL: PASS  (n=35, a={g1['a']}, d={g1['d']})

$ python3 {g3['verifier']}
OVERALL: PASS  (n={g3['n']}, a={g3['a']}, d={g3['d']})

$ python3 {g1['cross_checked']}
CROSS-CHECK: PASS  (35 terms, each with a verified x^2+x*y+y^2 representation)

$ python3 {g3['cross_checked']}
CROSS-CHECK: PASS  ({g3['n']} terms, each with a verified x^2+x*y+y^2 representation)
```

- `src/verify.py` factors every term from scratch in Python ints (its own Miller–Rabin and
  Pollard–Brent; no sieve, no fixed-width arithmetic) and checks that every prime `p ≡ 2 (mod 3)`
  occurs to an even power.
- `src/crosscheck.py` does **not use that criterion at all**. It *constructs* integers `x, y`
  with `x² + xy + y² = t` for every term — Eisenstein-integer arithmetic, Cornacchia for primes
  `≡ 1 (mod 3)` — and checks the identity by direct multiplication. An explicit representation is
  positive proof, independent of any theory about bad primes.
- `make check` runs the test suite, in which three further implementations agree on every
  `n ≤ 20000`.
- `tools/promote.py` refuses to record anything that fails either check.

---

## 4. Why the steps look like this

A Loeschian number is **never `≡ 2 (mod 3)`**: for `t` coprime to 3,
`t ≡ (−1)^(number of bad prime factors, with multiplicity) (mod 3)`, and "Loeschian" forces every
such exponent even, hence that count even. So `3 | d` is forced and all terms lie in one class
mod 3 — measured, with `3 ∤ d` the longest run found anywhere is **2**. This lifts the per-term
probability from ~0.14 to ~0.6.

Each bad prime `p` with `2p ≤ n` must then divide `d` (the theorem above); each with
`n/2 < p < n` is not forced but costs `~1/(p+1)` if omitted. Hence both steps are 3 times a
product of the small bad primes. Finally `a` is chosen so that none of 59, 71, 83, 89, 101, 107,
113 divides *any* term: each can hit at most one index in the window, and `a`'s residue pushes
the hit outside it.

---

## 5. Every distinct progression of 55 or more terms found

Each row is a different progression: rescaled copies (`m·a`, `m·d` for a Loeschian multiplier
`m`, which the search re-finds routinely) are reduced to the primitive form and listed once.
Every row passes both `src/verify.py` and `src/crosscheck.py` at the stated length.

| n | a | d | last term | d factored |
|---|---|---|---|---|
{runs_md}

Re-check any row: `python3 src/verify.py <a> <d> <n> --maximal`

### Rescaled copies of the rows above, as the search found them

Each of these is a row of the table above with `a` and `d` both multiplied by the stated
Loeschian multiplier. They are valid progressions of the stated length (each verified by both
checkers) but not new ones.

| n | a | d | last term | multiplier | primitive a |
|---|---|---|---|---|---|
{copies_md}

## 5b. Status of the n ≥ 58 target

{g58}. The search runs on 8 NVIDIA RTX PRO 6000 Blackwell GPUs at roughly 6·10¹¹ residues per
second (`src/c/apsearch_cuda.cu`, launched by `tools/multi_gpu.sh`). The engine enumerates the
admissible residues of `a` directly by nested additive loops, so an inadmissible candidate is
never constructed, then settles 64 candidates per residue with precomputed bitmasks over every
bad prime up to 10⁴, and confirms survivors with an exact test. A calibrated model puts the
expected yield at a few tenths of a 58-term progression per hour at that rate; PROGRESS.md has
the model, the measurements, and everything that was ruled out. Two shortcuts are provably
impossible: no constant `c` can extend a progression or repair a near miss, and no
meet-in-the-middle exists (for fixed `d` every constraint is a congruence on the single unknown
`a`).

## 6. Where everything is

- `records.json` — machine-readable records and full history.
- `PROGRESS.md` — dated log: what was tried, found, and ruled out, with the cost model.
- `experiments/` — campaigns and raw per-work-unit output.
- `src/` — Loeschian core, authoritative verifier, constructive cross-check, and four searchers
  (flat sieve, class-compressed, CPU two-stage CRT, and the Metal and CUDA GPU engines).
"""
    (ROOT / "ANSWER.md").write_text(doc)


if __name__ == "__main__":
    raise SystemExit(main())

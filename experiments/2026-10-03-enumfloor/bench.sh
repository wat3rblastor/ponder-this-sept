#!/bin/bash
# Cost-per-output sweep.  Column meanings are in README.md.
cd "$(dirname "$0")"
D0=382160924970
X0=21783172723290      # 57*D0, the natural term floor at K=1
E=./enum
run() { echo "### $*"; $E "$@" 2>&1; echo; }

echo "=== A. headline: pin 59..113, X = 57*D0 = $X0, auto split ==="
run --Q 59,71,83,89,101,107,113 --X $X0
echo "=== A'. same, 113 pulled into the wheel (M1 > X, i.e. pure CRT walk w/ prune) ==="
run --Q 59,71,83,89,101,107,113 --X $X0 --split 7
echo "=== A''. same, verified sample (brute-force inverse verifier, 2e6 outputs) ==="
run --Q 59,71,83,89,101,107,113 --X $X0 --verify --maxout 2000000

echo "=== B. X sweep, pin 59..113, auto split ==="
for X in 1000000000000 10000000000000 21783172723290 100000000000000 1000000000000000 10000000000000000 100000000000000000; do
  run --Q 59,71,83,89,101,107,113 --X $X --maxout 30000000
done

echo "=== C. pin-set sweep at X = 57*D0 (each extra pinned prime costs ~q/(q-58)) ==="
run --Q 59,71,83,89,101 --X $X0 --maxout 30000000
run --Q 59,71,83,89,101,107 --X $X0 --maxout 30000000
run --Q 59,71,83,89,101,107,113 --X $X0 --maxout 30000000
run --Q 59,71,83,89,101,107,113,131 --X $X0 --maxout 30000000
run --Q 59,71,83,89,101,107,113,131,137 --X $X0 --maxout 30000000
run --Q 59,71,83,89,101,107,113,131,137,149 --X $X0 --maxout 30000000
run --Q 59,71,83,89,101,107,113,131,137,149,167 --X $X0 --maxout 30000000
run --Q 59,71,83,89,101,107,113,131,137,149,167,173 --X $X0 --maxout 30000000

echo "=== D. naive full-residue-enumeration baseline (split = all primes, no mid-split) ==="
echo "    cost per output should be ~ M/X ; shown where it is cheap enough to run"
run --Q 59,71,83,89,101,107,113 --X 100000000000 --split 7 --maxout 3000000
run --Q 59,71,83,89,101,107,113 --X 10000000000 --split 7 --maxout 3000000
run --Q 59,71,83,89,101,107,113 --X 1000000000 --split 7 --maxout 3000000

echo "=== E. engine-equivalent baseline: same pins, term floor 64*MOD ==="
echo "    (MOD = 30*59*71*83*89*101*107*113 = 1133661268029390, X = 64*MOD)"
run --Q 59,71,83,89,101,107,113 --X 72554321153880960 --maxout 30000000

echo "=== F. at the engine's own term floor X = 64*MOD = 7.2554e16: does pinning MORE pay? ==="
for Q in 59,71,83,89,101,107,113 59,71,83,89,101,107,113,131 \
         59,71,83,89,101,107,113,131,137 59,71,83,89,101,107,113,131,137,149,167,173; do
  run --Q $Q --X 72554321153880960 --maxout 30000000
done

echo "=== G. large verified samples (independent brute-force-inverse verifier) ==="
run --Q 59,71,83,89,101,107 --X $X0 --verify --maxout 20000000
run --Q 59,71,83,89,101,107,113 --X $X0 --verify --maxout 20000000
run --Q 59,71,83,89,101,107,113,131,137,149,167,173 --X $X0 --verify --maxout 5000000

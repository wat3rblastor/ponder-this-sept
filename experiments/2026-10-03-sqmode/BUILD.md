# Build / run

```sh
cc -O3 -march=native -std=gnu11 -Wall -Wextra -I src/c \
   -o build/apsearch_sq experiments/2026-10-03-sqmode/apsearch_sq.c -lm
```

`src/c/apsearch.c` is untouched; `apsearch_sq.c` is a copy plus the square
option.

## New flags

| flag | meaning |
|---|---|
| `--ind q,q,..`   | bracket primes put IN d: `d = K * D0 * prod(ind)` |
| `--square q,q,..`| bracket primes put in SQUARE mode: q must not divide d, and stage 1 gains a component (modulus `q^2`, start `-(max(0,n-q))*d`, step `-d`, count `min(n-1,q-1)-max(0,n-q)+1`, i.e. `2q-n` when `q<n`) |
| `--mintierb N`   | minimum pinned tier-B primes before a K is skipped (square modes eat the modulus budget, so this has to be loosened from the hard-coded 4) |
| `--histmin N`    | smallest run length whose last term is tracked (`minlast`) |
| `--tlimit S`     | stop after S seconds (matched wall-clock measurement) |
| `--atest A`      | do not search; report whether this exact `a` passes stage 1 / stage 2 and how long its Loeschian run is. Cheap form of a record regression: a full unit at the record's modulus is ~5e9 residues. |

The mode is part of the unit identity: every `--out` record carries
`"ind"` and `"sq"`, and the `*** n=` lines print them too.

Square-mode primes are excluded from tier B **and from tier C** (tier C's mask
is "q misses the whole window", the exact complement of the square family;
applying it would reject every candidate the mode generates). Both exclusions
are explicit, not incidental.

Two self-checks run per unit: the original CRT check (now covering the `q^2`
components), plus a new semantic check that every residue of a `q^2` component
really makes `q^2` divide exactly one in-window term and `q` divide no other.

## The runs in raw/

```sh
# decisive test: the PROVED-optimal 35-term AP, which lives in the square family
./build/apsearch_sq --nterms 35 --D0 5610 --ind 23 --square 29 \
    --kmin 21 --kmax 21 --shifts 1 --report 33 --histmin 33 \
    --modcap 200000000000

# record regressions (admissibility probe)
./build/apsearch_sq --nterms 47 --D0 3741870 --ind 41,47,53 --kmin 205 --kmax 205 \
    --atest 2646171143023357
./build/apsearch_sq --nterms 35 --D0 5610 --ind 23 --square 29 --kmin 21 --kmax 21 \
    --modcap 200000000000 --atest 219830911

# matched wall clock at n = 58 (900 s each, same modcap, same shifts, same b2)
./build/apsearch_sq --nterms 58 --D0 3741870 --ind 41,47,53 ...   # raw/n58_alld.*
./build/apsearch_sq --nterms 58 --D0 3741870 --ind 41,47 --square 53 ...  # raw/n58_sq53.*
./build/apsearch_sq --nterms 58 --D0 3741870 --square 41,47,53 ...       # raw/n58_sq.*
```

Equivalence with production (`build/apsearch`), identical args, n=58 K=7:
both report `MOD=93760753290 R=1732900 surv=19 conf=19 best=18
a=22225747809637 d=2675126474790` — the all-in-d path is unchanged.

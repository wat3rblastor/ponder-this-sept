/* enum.c -- enumerate CRT-admissible starting values a < X at O(1) cost each.
 *
 * Problem.  Fix d and a set Q of bad primes q > 58 with q ! d.  Call a
 * admissible at q iff  a = j*d (mod q)  for some 1 <= j <= q-58  (so q divides
 * no term of a, a+d, ..., a+57d).  Enumerate every a < X admissible at every
 * q in Q.  Hard case: X << prod(Q).
 *
 * Algorithm (wheel + filter, = meet-in-the-middle with the split chosen so one
 * side fits under X).  Split Q = Q1 u Q2 with M1 = prod(Q1) <= X.
 *   - Q1: enumerate the delta1*M1 = prod_{Q1}(q-58) admissible residues mod M1
 *     by a nested additive CRT loop, O(1) amortised per residue (no mod, no
 *     multiply: one add and one conditional subtract per level).
 *   - lift each residue r by t*M1 while r + t*M1 < X  (floor(X/M1) lifts).
 *   - Q2: reject with a table lookup per prime, a mod q maintained by addition.
 * Cost per emitted a = prod_{q in Q2} q/(q-58) * (|Q2| cheap ops), i.e. O(1)
 * when Q2 is small, and never worse than q_max/(q_max-58) ~ 2 per leftover
 * prime.
 *
 * Build: cc -O3 -march=native -o enum enum.c
 * Usage: ./enum --d <d> --X <X> --Q 59,71,83,89,101,107,113 [--split k]
 *               [--emit] [--verify] [--maxout N]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef unsigned long long u64;
typedef __uint128_t u128;

#define MAXP 16

static int nq, Q[MAXP], split;
static u64 D, X, M1;
/* wheel side */
static int w_n, w_q[MAXP], w_lo[MAXP], w_hi[MAXP];
static u64 w_step[MAXP];          /* (e_q * d) mod M1, the additive increment */
static u64 w_base;                /* value with all j_q = lo  */
/* filter side */
static int f_n, f_q[MAXP];
static unsigned char f_ok[MAXP][512];   /* f_ok[i][a % q] = admissible? */

static u64 cand = 0, out = 0, leaves = 0;     /* candidates tested (= inner iterations) */
static int do_emit = 0, do_verify = 0;
static long long maxout = -1;
static u64 verified = 0, badver = 0;
static u64 *emitbuf = NULL; static u64 emitn = 0, emitcap = 0;

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }

static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1 % m; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
/* independent, dumb inverse: brute-force search.  Used only by the verifier. */
static u64 inv_brute(u64 a, u64 p) {
    a %= p;
    for (u64 x = 1; x < p; x++) if ((a * x) % p == 1) return x;
    fprintf(stderr, "no inverse of %llu mod %llu\n", a, p); exit(1);
}

/* verifier: completely independent of the tables used by the search. */
static void verify(u64 a) {
    for (int i = 0; i < nq; i++) {
        u64 q = Q[i];
        u64 j = ((a % q) * inv_brute(D % q, q)) % q;   /* a = j*d (mod q) */
        if (!(j >= 1 && j <= q - 58)) { badver++; return; }
    }
    verified++;
}

static inline void leaf(u64 r) {
    leaves++;
    /* lift r by t*M1 while < X, filter by Q2 */
    u64 a = r;
    u64 m2[MAXP];
    int rq[MAXP];
    for (int i = 0; i < f_n; i++) { rq[i] = (int)(a % f_q[i]); m2[i] = M1 % f_q[i]; }
    while (a < X) {
        cand++;
        int ok = 1;
        for (int i = 0; i < f_n; i++) if (!f_ok[i][rq[i]]) { ok = 0; break; }
        if (ok) {
            out++;
            if (do_verify) verify(a);
            if (do_emit) {
                if (emitn == emitcap) { emitcap = emitcap ? emitcap*2 : 1024;
                    emitbuf = realloc(emitbuf, emitcap*sizeof(u64)); }
                emitbuf[emitn++] = a;
            }
            if (maxout >= 0 && (long long)out >= maxout) return;
        }
        a += M1;
        for (int i = 0; i < f_n; i++) {
            rq[i] += (int)m2[i];
            if (rq[i] >= f_q[i]) rq[i] -= f_q[i];
        }
    }
}

/* nested additive CRT walk over the wheel side */
static void walk(int lvl, u64 v) {
    if (lvl == w_n) { leaf(v); return; }
    u64 s = w_step[lvl];
    for (int j = w_lo[lvl]; j <= w_hi[lvl]; j++) {
        walk(lvl + 1, v);
        v += s; if (v >= M1) v -= M1;
        if (maxout >= 0 && (long long)out >= maxout) return;
    }
}

int main(int argc, char **argv) {
    D = 382160924970ULL; X = 21783172723290ULL; split = -1;
    const char *qs = "59,71,83,89,101,107,113";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--d")) D = strtoull(argv[++i], 0, 10);
        else if (!strcmp(argv[i], "--X")) X = strtoull(argv[++i], 0, 10);
        else if (!strcmp(argv[i], "--Q")) qs = argv[++i];
        else if (!strcmp(argv[i], "--split")) split = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--emit")) do_emit = 1;
        else if (!strcmp(argv[i], "--verify")) do_verify = 1;
        else if (!strcmp(argv[i], "--maxout")) maxout = atoll(argv[++i]);
        else { fprintf(stderr, "bad arg %s\n", argv[i]); return 1; }
    }
    { char buf[256]; strncpy(buf, qs, 255); buf[255]=0;
      for (char *t = strtok(buf, ","); t; t = strtok(0, ",")) Q[nq++] = atoi(t); }

    /* choose the split: largest prefix with prod <= X */
    if (split < 0) {
        u128 p = 1; split = 0;
        for (int i = 0; i < nq; i++) { if ((u128)p * Q[i] > (u128)X) break; p *= Q[i]; split = i+1; }
    }
    M1 = 1; for (int i = 0; i < split; i++) M1 *= Q[i];
    w_n = split; f_n = nq - split;
    for (int i = 0; i < split; i++) {
        w_q[i] = Q[i]; w_lo[i] = 1; w_hi[i] = Q[i] - 58;
        if (w_hi[i] < w_lo[i]) { fprintf(stderr, "prime %d has no admissible residue\n", Q[i]); return 1; }
    }
    for (int i = 0; i < f_n; i++) {
        int q = f_q[i] = Q[split + i];
        if (q >= 512) { fprintf(stderr, "prime %d exceeds table\n", q); return 1; }
        u64 id = powmod(D % q, q - 2, q);          /* Fermat inverse */
        for (int x = 0; x < q; x++) {
            u64 j = ((u64)x * id) % q;
            f_ok[i][x] = (j >= 1 && j <= (u64)(q - 58));
        }
    }
    /* additive increments: e_q*d mod M1 where e_q = 1 mod q, 0 mod M1/q */
    w_base = 0;
    for (int i = 0; i < w_n; i++) {
        u64 q = w_q[i], c = M1 / q;
        u64 e = mulmod(c, powmod(c % q, q - 2, q), M1);   /* idempotent */
        w_step[i] = mulmod(e, D % M1, M1);
        w_base = (w_base + mulmod(w_step[i], (u64)w_lo[i], M1)) % M1;
    }
    double wheelres = 1; for (int i = 0; i < w_n; i++) wheelres *= (w_hi[i]-w_lo[i]+1);
    double delta = 1; for (int i = 0; i < nq; i++) delta *= (double)(Q[i]-58)/Q[i];

    fprintf(stderr, "d=%llu X=%llu  Q1(wheel)=%d primes M1=%llu  Q2(filter)=%d primes\n",
            D, X, w_n, M1, f_n);
    fprintf(stderr, "wheel residues=%.0f  lifts=%llu  1/delta=%.0f  expected out=%.4g\n",
            wheelres, X / M1 + 1, 1/delta, delta*X);

    struct timespec t0, t1; clock_gettime(CLOCK_MONOTONIC, &t0);
    walk(0, w_base);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sec = (t1.tv_sec-t0.tv_sec) + 1e-9*(t1.tv_nsec-t0.tv_nsec);

    double work = (double)leaves + (double)cand;
    printf("out=%llu cand=%llu leaves=%llu sec=%.3f work_per_out=%.4f"
           " ns_per_out=%.2f expected_out=%.6g verified=%llu badver=%llu\n",
           out, cand, leaves, sec, out ? work/out : 0.0,
           out ? 1e9*sec/out : 0.0, delta*X, verified, badver);
    if (do_emit) { for (u64 i = 0; i < emitn; i++) printf("A %llu\n", emitbuf[i]); }
    return badver ? 2 : 0;
}

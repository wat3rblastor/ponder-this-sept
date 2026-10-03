/* exact.c -- EXACT counts of n-term Loeschian arithmetic progressions.
 *
 * Independent bit-parallel re-implementation of
 * experiments/2026-10-03-assumptions/joint.c, used to (a) reproduce its eight
 * published exact counts, (b) extend the exact table to larger n and larger X.
 *
 * CONVENTION (stated because the repo has been burned by it):
 *   COUNT  = #PAIRS (a,d), d = base(n)*m with m >= 1, a >= 1,
 *            a = 1 (mod 3), a != 0 (mod q) for every bad prime q | base(n),
 *            a + (n-1)d <= X, and  a + k d  Loeschian for all 0 <= k < n.
 *            An AP of length L >= n contributes L-n+1 to COUNT, so COUNT is
 *            "APs of length >= n counted by starting pair".  This is exactly
 *            what joint.c reports.  It is NOT a count of maximal APs.
 *   RUNSTART = same, restricted to pairs whose run cannot be extended left
 *            (a-d not Loeschian or < 1): i.e. #runs of length >= n.
 *   MAXIMAL  = runs of length >= n that also cannot be extended right inside
 *            [1,X]: i.e. #maximal runs of length >= n.
 *
 * base(n) = 3 * prod{bad q : 2q <= n}: the fully forced divisors of d.
 *
 * Loeschian bitmap: m >= 1 is Loeschian (m = x^2+xy+y^2) iff every prime
 * p = 2 (mod 3) divides m to an even power.  Built by parity-toggling over
 * bad p <= sqrt(X), plus a single pass over bad p > sqrt(X) (which can occur
 * only to the first power, since (S+1)^2 > X).
 *
 * build: gcc-15 -O3 -fopenmp -o exact exact.c -lm
 * usage: ./exact X nlo nhi [nstep] [--check] [--list]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#else
static double omp_get_wtime(void) { return 0; }
static int omp_get_max_threads(void) { return 1; }
#endif

typedef uint64_t u64;

static long X, W;
static u64 *L, *AD, *B, *C, *D;

#define GET(p,i) (((p)[(i)>>6] >> ((i)&63)) & 1ULL)
#define SET(p,i) ((p)[(i)>>6] |= 1ULL << ((i)&63))
#define CLR(p,i) ((p)[(i)>>6] &= ~(1ULL << ((i)&63)))
#define TOG(p,i) ((p)[(i)>>6] ^= 1ULL << ((i)&63))

static int isbad(long p) { return p % 3 == 2; }

static u64 *alloc_bits(long nw) {
    u64 *p = calloc(nw + 40, sizeof(u64));
    if (!p) { fprintf(stderr, "OOM\n"); exit(1); }
    return p;
}

/* ---- odd-only bit-packed prime sieve up to N ---- */
static u64 *psv;
static void prime_sieve(long N) {
    long half = (N >> 1) + 1, nw = (half >> 6) + 2;
    psv = alloc_bits(nw);
    for (long w = 0; w < nw; w++) psv[w] = ~0ULL;   /* bit i -> 2i+1 prime */
    CLR(psv, 0);                                     /* 1 */
    for (long i = 1; (2 * i + 1) * (2 * i + 1) <= N; i++) {
        if (!GET(psv, i)) continue;
        long p = 2 * i + 1;
        for (long j = p * p; j <= N; j += 2 * p) CLR(psv, j >> 1);
    }
}
static int isprime_s(long p) {
    if (p == 2) return 1;
    if (p < 2 || !(p & 1)) return 0;
    return (int)GET(psv, p >> 1);
}

static void build_loeschian(void) {
    long S = 1; while ((S + 1) * (S + 1) <= X) S++;
    u64 *sc = alloc_bits(W);
    u64 *od = alloc_bits(W);
    for (long p = 2; p <= S; p++) {
        if (!isprime_s(p) || !isbad(p)) continue;
        for (long pk = p;; pk *= p) {
            for (long m = pk; m <= X; m += pk) TOG(sc, m);
            if (pk > X / p) break;
        }
        for (long m = p; m <= X; m += p)
            if (GET(sc, m)) { SET(od, m); CLR(sc, m); }
    }
    for (long q = S + 1; q <= X; q++) {
        if (!isprime_s(q) || !isbad(q)) continue;
        for (long m = q; m <= X; m += q) SET(od, m);
    }
    free(sc);
    L = od;
    for (long w = 0; w < W; w++) L[w] = ~L[w];
    for (long i = X + 1; i < (W + 40) * 64; i++) CLR(L, i);
    SET(L, 0);
}

static long badlist[16]; static int nbad;
static long FORCEBASE = 0;      /* --base B: use the sub-family d = B*m instead
                                 * of the minimal forced base(n).  Legitimate:
                                 * any B with base(n) | B defines an exactly
                                 * countable sub-family, and B = 3741870 is the
                                 * actual n=58 family. */
static long base_of(int n) {
    long b = 3; nbad = 0;
    for (long q = 2; 2 * q <= n; q++)
        if (isprime_s(q) && isbad(q)) { b *= q; badlist[nbad++] = q; }
    if (FORCEBASE) {
        if (FORCEBASE % b) { fprintf(stderr, "base %ld not a multiple of base(%d)=%ld\n", FORCEBASE, n, b); exit(1); }
        b = FORCEBASE; nbad = 0;
        for (long q = 2; q <= 200; q++)
            if (isprime_s(q) && isbad(q) && b % q == 0) badlist[nbad++] = q;
    }
    return b;
}

static void build_adm(void) {
    memset(AD, 0, (size_t)(W + 40) * sizeof(u64));
    for (long a = 1; a <= X; a += 3) SET(AD, a);
    for (int i = 0; i < nbad; i++)
        for (long a = badlist[i]; a <= X; a += badlist[i]) CLR(AD, a);
    CLR(AD, 0);
}

/* exact #admissible a in [1,hi] by inclusion-exclusion */
static long long adm_count(long hi) {
    if (hi < 1) return 0;
    long long tot = 0;
    for (int S = 0; S < (1 << nbad); S++) {
        long Q = 1; int bits = 0;
        for (int i = 0; i < nbad; i++) if (S >> i & 1) { Q *= badlist[i]; bits++; }
        long M = 3 * Q, r = -1;
        for (long x = 0; x < M; x++) if (x % 3 == 1 && x % Q == 0) { r = x; break; }
        if (r == 0) r = M;
        long long c = (r > hi) ? 0 : (hi - r) / M + 1;
        tot += (bits & 1) ? -c : c;
    }
    return tot;
}

/* dst = src1 AND (src2 >> off bits), words [0,nw) */
static void shand(u64 *dst, const u64 *src1, const u64 *src2, long off, long nw) {
    long ws = off >> 6, bs = off & 63;
    if (bs == 0) {
        #pragma omp parallel for schedule(static)
        for (long w = 0; w < nw; w++) dst[w] = src1[w] & src2[w + ws];
    } else {
        #pragma omp parallel for schedule(static)
        for (long w = 0; w < nw; w++)
            dst[w] = src1[w] & ((src2[w + ws] >> bs) | (src2[w + ws + 1] << (64 - bs)));
    }
}

/* popcount of (src AND AD) over bits [1,hi] */
static long long pc_masked(const u64 *src, long hi) {
    if (hi < 1) return 0;
    long nw = (hi >> 6) + 1, b = hi & 63;
    long long t = 0;
    #pragma omp parallel for schedule(static) reduction(+:t)
    for (long w = 0; w < nw; w++) {
        u64 v = src[w] & AD[w];
        if (w == nw - 1 && b < 63) v &= (1ULL << (b + 1)) - 1;
        t += (long long)__builtin_popcountll(v);
    }
    return t;
}
/* popcount of (AD >> off) AND AD over [1,hi]  -- not needed; kept out */

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: exact X nlo nhi [step] [--check] [--list]\n"); return 1; }
    X = atol(argv[1]);
    int nlo = argc > 2 ? atoi(argv[2]) : 12, nhi = argc > 3 ? atoi(argv[3]) : 18;
    int step = (argc > 4 && argv[4][0] != '-') ? atoi(argv[4]) : 2;
    int dolist = 0, docheck = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--list")) dolist = 1;
        if (!strcmp(argv[i], "--check")) docheck = 1;
        if (!strcmp(argv[i], "--base")) FORCEBASE = atol(argv[i + 1]);
    }
    W = (X >> 6) + 2;
    double t0 = omp_get_wtime();
    prime_sieve(X);
    build_loeschian();
    AD = alloc_bits(W); B = alloc_bits(W); C = alloc_bits(W); D = alloc_bits(W);
    fprintf(stderr, "# build %.1fs, threads %d\n", omp_get_wtime() - t0, omp_get_max_threads());

    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "--spot")) {
        long v; while (scanf("%ld", &v) == 1) printf("%ld %d\n", v, (int)GET(L, v));
        return 0;
    }
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "--dens")) {
        /* EXACT per-term Loeschian density in the conditioning class of each
         * forced family, in 64 log-spaced bins of t.  Every term of every AP in
         * family base lies in exactly this class, so this IS the model's rho. */
        int ns[5] = {12, 22, 34, 46, 58};
        int nn = 5;
        if (FORCEBASE) { ns[0] = 4; nn = 1; }   /* base(4)=6 divides every base used here */
        for (int j = 0; j < nn; j++) {
            long base = base_of(ns[j]);
            build_adm();
            printf("DENS base=%ld n=%d:", base, ns[j]);
            for (int b = 0; b < 64; b++) {
                double f0 = pow((double)X, b / 64.0), f1 = pow((double)X, (b + 1) / 64.0);
                long lo = (long)f0, hi = (long)f1; if (hi > X) hi = X;
                if (hi <= lo) continue;
                long long nc = 0, nl = 0;
                for (long i = lo; i < hi; i++) if (GET(AD, i)) { nc++; nl += GET(L, i); }
                if (nc) printf(" %ld,%ld,%lld,%lld", lo, hi, nc, nl);
            }
            printf("\n"); fflush(stdout);
        }
        return 0;
    }
    printf("# X=%ld\n# first 60 Loeschian:", X);
    { int c = 0; for (long i = 0; i <= X && c < 60; i++) if (GET(L, i)) { printf(" %ld", i); c++; } }
    printf("\n");
    long long cnt = 0;
    for (long w = 0; w < W; w++) cnt += __builtin_popcountll(L[w]);
    cnt -= 1;
    printf("# Loeschian count in [1,X] = %lld  density %.6f  K=count*sqrt(lnX)/X = %.6f\n",
           cnt, (double)cnt / X, (double)cnt * sqrt(log((double)X)) / X);

    for (int n = nlo; n <= nhi; n += step) {
        long base = base_of(n);
        long mmax = (X - 1) / ((long)(n - 1) * base);
        if (mmax < 1) { printf("n=%d base=%ld NO_D\n", n, base); continue; }
        build_adm();
        double t1 = omp_get_wtime();

        if (docheck) {
            /* EXACT class density rho_class(t) in log-spaced bins: the terms of
             * every AP in this family lie in exactly this class. */
            printf("# CLASSDENS n=%d base=%ld bins(lo,hi,nclass,nloe):", n, base);
            for (int b = 0; b < 48; b++) {
                double f0 = pow((double)X, b / 48.0), f1 = pow((double)X, (b + 1) / 48.0);
                long lo = (long)f0, hi = (long)f1;
                if (hi <= lo) continue;
                long long nc = 0, nl = 0;
                for (long i = lo; i < hi; i++) if (GET(AD, i)) { nc++; nl += GET(L, i); }
                if (nc) printf(" %ld,%ld,%lld,%lld", lo, hi, nc, nl);
            }
            printf("\n");
        }

        long long pairs = 0, ok1 = 0, found = 0, runstart = 0, maximal = 0;
        long fa = 0, fd = 0; int nshow = 0;
        for (long m = 1; m <= mmax; m++) {
            long d = base * m, hi = X - (long)(n - 1) * d;
            if (hi < 1) continue;
            pairs += adm_count(hi);
            ok1 += pc_masked(L, hi);

            /* addition chain: C = AND_{k<have} (L >> k d), valid on [0, X-(have-1)d] */
            long nw2 = ((X - d) >> 6) + 10; if (nw2 > W) nw2 = W;
            shand(B, L, L, d, nw2);
            long have = 2;
            memcpy(C, B, (size_t)(nw2 + 10) * sizeof(u64));
            while (have * 2 <= n) {
                long nwn = ((X - (2 * have - 1) * d) >> 6) + 10; if (nwn > W) nwn = W;
                shand(D, C, C, have * d, nwn);
                u64 *t = C; C = D; D = t; have *= 2;
            }
            while (have < n) {
                long add = (n - have >= 2) ? 2 : 1;
                long nwn = ((X - (have + add - 1) * d) >> 6) + 10; if (nwn > W) nwn = W;
                if (add == 2) shand(D, C, B, have * d, nwn);
                else          shand(D, C, L, have * d, nwn);
                u64 *t = C; C = D; D = t; have += add;
            }
            long long f = pc_masked(C, hi);
            found += f;
            if (f) {
                for (long w = 0; w <= (hi >> 6); w++) {
                    u64 v = C[w] & AD[w];
                    if (w == (hi >> 6)) { long b = hi & 63; if (b < 63) v &= (1ULL << (b + 1)) - 1; }
                    while (v) {
                        long a = w * 64 + __builtin_ctzll(v); v &= v - 1;
                        int left = (a - d >= 1 && GET(L, a - d));
                        long nxt = a + (long)n * d;
                        int right = (nxt <= X && GET(L, nxt));
                        if (!left) runstart++;
                        if (!left && !right) maximal++;
                        if (!fa) { fa = a; fd = d; }
                        if (dolist && nshow < 200) {
                            printf("#   AP a=%ld d=%ld m=%ld last=%ld left=%d right=%d\n",
                                   a, d, m, a + (long)(n - 1) * d, left, right);
                            nshow++;
                        }
                    }
                }
            }
        }
        printf("n=%d base=%ld nd=%ld pairs=%lld rho1=%.6f COUNT=%lld RUNSTART=%lld MAXIMAL=%lld ex=%ld,%ld t=%.1fs\n",
               n, base, mmax, pairs, (double)ok1 / (double)pairs, found, runstart, maximal, fa, fd,
               omp_get_wtime() - t1);
        fflush(stdout);
    }
    return 0;
}

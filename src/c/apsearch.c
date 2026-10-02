/* apsearch -- Wroblewski-style two-stage search for long Loeschian APs.
 *
 * Build: cc -O3 -march=native -std=c11 -o build/apsearch src/c/apsearch.c
 *
 * STRUCTURE (see GOAL.md 2, 2b; technique after J. Wroblewski, "How to search
 * for 26 primes in arithmetic progression?", adapted from primes to Loeschian
 * numbers: "prime" -> "Loeschian", "all primes" -> "bad primes p = 2 mod 3").
 *
 *   d = K * D0,  D0 = 3 * 2*5*11*17*23*29*41*47*53 = 382160924970
 *
 * Bad primes split into three tiers:
 *
 *   TIER A (p | D0): no term is divisible by p at all, we only need p doesn't
 *       divide a. Handled in stage 2 as "a mod p != 0".
 *       Plus 3: every Loeschian number is 0 or 1 mod 3 and the primitive form
 *       has 3 doesn't divide a, so a = 1 mod 3. Handled in stage 1.
 *
 *   TIER B (the smallest bad primes q > n): at most one index is hit, so a can
 *       be chosen mod q to push the hit outside [0, n-1]. Good residues are
 *       a = -k*d (mod q) for k = n .. q-1, i.e. q-n of them, an arithmetic
 *       progression in a with common difference -d. So the CRT-admissible set
 *       over all of tier B is enumerated DIRECTLY by nested additive loops --
 *       an inadmissible a is never materialized. This is the whole trick: it
 *       buys prod(q/(q-n)) ~ 3e4 over sieving a.
 *
 *   TIER C (bad primes from the next one up to --b2): handled by 64-bit
 *       bitmask words. One residue R stands for 64 candidates
 *       a = R + (b+shift)*MOD, b = 0..63, so one modulo + one AND decides 64
 *       candidates at once. AND-chain short-circuits as soon as it hits zero.
 *
 * Stage 3 confirms survivors with an exact Loeschian test (factorization,
 * parity of bad-prime exponents). src/verify.py remains the only thing that
 * mints a record.
 *
 * Every (K, shift) pair is an independent unit of work; parallel workers must
 * be given disjoint --kmin/--kmax windows.
 */

#include <inttypes.h>
#include <stdbool.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef int64_t i64;
typedef __uint128_t u128;

#define D0_DEFAULT 382160924970ULL   /* 3 * 2*5*11*17*23*29*41*47*53 */
#define MAXTIERB 12
#define MAXTIERC 512

static volatile sig_atomic_t stop_requested = 0;
static void on_signal(int s) { (void)s; stop_requested = 1; }

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* ------------------------------------------------------------ modular ----- */

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }

static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}

static u64 inv_mod(u64 a, u64 m) {
    i64 t = 0, newt = 1, r = (i64)m, newr = (i64)(a % m);
    while (newr) { i64 q = r / newr;
        i64 tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp; }
    if (r != 1) return 0;
    return (u64)(t < 0 ? t + (i64)m : t);
}

static bool is_prime_u64(u64 n) {
    if (n < 2) return false;
    for (u64 p = 2; p < 38; p++) { if (p * p > n) return true;
        if (n % p == 0) return n == p; }
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    static const u64 wit[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) {
        u64 x = powmod(wit[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int j = 1; j < s; j++) { x = mulmod(x, x, n);
            if (x == n - 1) { ok = 1; break; } }
        if (!ok) return false;
    }
    return true;
}

static u64 gcd_u64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

/* Pollard-Brent with a BATCHED gcd. The textbook Floyd version takes one gcd
 * per iteration, which dominated everything: at ~1e4 iterations for a 1e16
 * cofactor that is ~1 ms per factorization and it was costing this search
 * roughly 240 of its 277 ns per residue. Accumulating the product of the
 * differences mod n and taking a single gcd every 128 steps removes almost all
 * of that. */
static u64 pollard(u64 n) {
    if (!(n & 1)) return 2;
    for (u64 c = 1;; c++) {
        u64 y = 2, m = 128, g = 1, r = 1, q = 1, x = 0, ys = 0;
        while (g == 1) {
            x = y;
            for (u64 i = 0; i < r; i++) y = (mulmod(y, y, n) + c) % n;
            for (u64 k = 0; k < r && g == 1; k += m) {
                ys = y;
                u64 lim = (m < r - k) ? m : r - k;
                for (u64 i = 0; i < lim; i++) {
                    y = (mulmod(y, y, n) + c) % n;
                    q = mulmod(q, x > y ? x - y : y - x, n);
                }
                g = gcd_u64(q, n);
            }
            r *= 2;
        }
        if (g == n) {                     /* back off one step at a time */
            g = 1;
            y = ys;
            while (g == 1) {
                y = (mulmod(y, y, n) + c) % n;
                g = gcd_u64(x > y ? x - y : y - x, n);
            }
        }
        if (g != n) return g;
    }
}

/* ------------------------------------------- exact Loeschian test (stage 3) */

static u32 *sp = NULL;        /* small primes for trial division */
static u64 n_sp = 0, sp_max = 0;

static void build_small_primes(u64 limit) {
    char *c = calloc(limit + 1, 1);
    for (u64 i = 2; i * i <= limit; i++)
        if (!c[i]) for (u64 j = i * i; j <= limit; j += i) c[j] = 1;
    u64 cnt = 0;
    for (u64 i = 2; i <= limit; i++) if (!c[i]) cnt++;
    sp = malloc(cnt * sizeof(u32));
    for (u64 i = 2; i <= limit; i++) if (!c[i]) sp[n_sp++] = (u32)i;
    sp_max = limit;
    free(c);
}

static u64 isqrt_u64(u64 n) {
    if (n == 0) return 0;
    u64 r = (u64)__builtin_sqrtl((long double)n);
    while (r > 0 && r > n / r) r--;
    while ((r + 1) <= n / (r + 1)) r++;
    return r;
}

/* Full factorization of the large cofactor into primes with multiplicity.
 * Every remaining factor exceeds sp_max, so there are at most a handful. */
static int nfac;
static u64 fac[64];

static void factor_rec(u64 n) {
    if (n == 1) return;
    if (is_prime_u64(n)) { fac[nfac++] = n; return; }
    u64 r = isqrt_u64(n);
    if (r * r == n) { factor_rec(r); factor_rec(r); return; }
    /* perfect cube (pollard-rho is unreliable on prime powers) */
    u64 c = (u64)__builtin_cbrtl((long double)n);
    for (u64 cc = (c > 1 ? c - 1 : 1); cc <= c + 1; cc++)
        if (cc > 1 && cc * cc * cc == n) {
            factor_rec(cc); factor_rec(cc); factor_rec(cc); return;
        }
    u64 f = pollard(n);
    factor_rec(f);
    factor_rec(n / f);
}

/* True iff t is Loeschian: every prime p = 2 (mod 3) divides t to an even
 * power. Exact for t < 2^64 -- this is the searcher's own test; src/verify.py
 * re-derives any hit independently in Python ints. */
static bool is_loeschian(u64 t) {
    if (t == 0) return true;
    for (u64 i = 0; i < n_sp; i++) {
        u32 p = sp[i];
        if ((u64)p * p > t) {
            /* what is left is 1 or a single prime */
            return (t == 1) ? true : (t % 3 != 2);
        }
        if (t % p) continue;
        int e = 0;
        while (t % p == 0) { t /= p; e++; }
        if (p % 3 == 2 && (e & 1)) return false;
    }
    if (t == 1) return true;
    /* Cofactor: all prime factors exceed sp_max. Factor it and check the
     * exponent parity of each BAD prime separately -- counting bad factors in
     * total would wrongly accept p*q for two distinct bad primes. */
    nfac = 0;
    factor_rec(t);
    for (int i = 0; i < nfac; i++) {
        if (fac[i] % 3 != 2) continue;
        int e = 0;
        for (int j = 0; j < nfac; j++) if (fac[j] == fac[i]) e++;
        if (e & 1) return false;
    }
    return true;
}

/* -------------------------------------------------------- Barrett modulo --- */

/* Modulo by a small prime, done in 32 bits.
 *
 * R < MOD < 2^51, so split R into three 17-bit limbs and fold with the
 * precomputed constants c1 = 2^17 mod m, c2 = 2^34 mod m:
 *     R mod m = (l0 + l1*c1 + l2*c2) mod m
 * The fold is at most 2^17 + 2*2^17*m < 2^32 for m <= 2000, so the division is
 * a 32-bit udiv rather than a 64-bit one. (A full 64-bit Barrett was tried
 * first and was wrong -- its magic truncated a 128-bit quotient to 64 bits.) */
typedef struct { u32 m, c1, c2; } bar;
static bar bar_make(u64 m) {
    bar b; b.m = (u32)m;
    b.c1 = (u32)(((u64)1 << 17) % m);
    b.c2 = (u32)(((u64)1 << 34) % m);
    return b;
}
static inline u32 bar_mod(u64 x, bar b) {
    u32 l0 = (u32)(x & 0x1FFFF), l1 = (u32)((x >> 17) & 0x1FFFF), l2 = (u32)(x >> 34);
    return (l0 + l1 * b.c1 + l2 * b.c2) % b.m;
}

/* ------------------------------------------------------------------ main --- */

static void usage(const char *p) {
    fprintf(stderr,
      "usage: %s --kmin K --kmax K [--nterms 58] [--shifts S] [--b2 2000]\n"
      "          [--report 50] [--out f.jsonl] [--D0 N] [--selftest]\n"
      "\n"
      "d = K*D0. Each (K, shift) is an independent work unit; give parallel\n"
      "workers disjoint --kmin/--kmax windows. --shifts S searches first terms\n"
      "a in [0, S*64*MOD).\n", p);
    exit(2);
}

int main(int argc, char **argv) {
    u64 kmin = 1, kmax = 0, nterms = 58, shifts = 16, b2 = 2000, D0 = D0_DEFAULT;
    u64 modcap = 6000000000000000ULL;   /* cap on MOD = 3*prod(tier B) */
    int report = 50;
    const char *out = NULL;
    bool selftest = false, isl_mode = false;

    for (int i = 1; i < argc; i++) {
        const char *k = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : (usage(argv[0]), ""))
        if (!strcmp(k, "--kmin")) kmin = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--kmax")) kmax = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--nterms")) nterms = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--shifts")) shifts = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--b2")) b2 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--report")) report = atoi(NEXT());
        else if (!strcmp(k, "--out")) out = NEXT();
        else if (!strcmp(k, "--D0")) D0 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--modcap")) modcap = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--selftest")) selftest = true;
        else if (!strcmp(k, "--isl")) isl_mode = true;
        else usage(argv[0]);
        #undef NEXT
    }
    if (!kmax) kmax = kmin;
    if (D0 % 3) { fprintf(stderr, "FATAL: 3 must divide D0\n"); return 2; }

    signal(SIGINT, on_signal); signal(SIGTERM, on_signal);
    build_small_primes(10000);

    if (isl_mode) {          /* read decimal values on stdin, print verdicts */
        u64 v;
        while (scanf("%" SCNu64, &v) == 1) printf("%" PRIu64 " %d\n", v, (int)is_loeschian(v));
        return 0;
    }

    if (selftest) {
        /* 1) known prefix of A003136 */
        static const int pre[] = {0,1,3,4,7,9,12,13,16,19,21,25,27,28,31,36,37,
                                  39,43,48,49,52,57,61,63,64,67,73,75,76,79,81,
                                  84,91,93,97,100};
        int j = 0, fail = 0;
        for (int n = 0; n <= 100; n++) {
            bool want = (j < 37 && pre[j] == n); if (want) j++;
            if (is_loeschian((u64)n) != want) {
                fprintf(stderr, "SELFTEST FAIL at n=%d\n", n); fail = 1; }
        }
        /* 2) bad primes and their powers */
        u64 bad[] = {2,5,11,17,23,29,41,47,53,59,71,83,89,101,1000037,1000121};
        for (unsigned i = 0; i < sizeof(bad)/sizeof(bad[0]); i++) {
            u64 p = bad[i];
            if (p % 3 != 2) continue;
            if (is_loeschian(p)) { fprintf(stderr, "SELFTEST FAIL: %" PRIu64
                                   " marked Loeschian\n", p); fail = 1; }
            if (!is_loeschian(p*p)) { fprintf(stderr, "SELFTEST FAIL: %" PRIu64
                                   "^2 not Loeschian\n", p); fail = 1; }
            if (is_loeschian(p*p*p)) { fprintf(stderr, "SELFTEST FAIL: %" PRIu64
                                   "^3 marked Loeschian\n", p); fail = 1; }
        }
        /* 3) Choudhry's parametric 9-term Loeschian AP (INTEGERS 25 (2025) #A14)
         *    a=2432167, d=1405152 -- an independent published test vector */
        for (int k = 0; k < 9; k++)
            if (!is_loeschian(2432167ULL + (u64)k * 1405152ULL)) {
                fprintf(stderr, "SELFTEST FAIL: Choudhry 9-term AP term %d\n", k);
                fail = 1; }
        if (is_loeschian(2432167ULL + 9ULL * 1405152ULL))
            fprintf(stderr, "note: Choudhry AP extends past 9 terms\n");
        /* 4) The case that a naive "parity of the total count of bad prime
         * factors" test gets WRONG: p*q with p, q distinct bad primes both
         * beyond the trial-division bound. Count parity is even, but each
         * exponent is odd, so it is NOT Loeschian. 1000037 and 1000121 are
         * both prime and both 2 mod 3. */
        if (is_loeschian(1000037ULL * 1000121ULL)) {
            fprintf(stderr, "SELFTEST FAIL: bad*bad semiprime\n"); fail = 1; }
        if (!is_loeschian(1000037ULL * 1000037ULL)) {
            fprintf(stderr, "SELFTEST FAIL: bad^2 large\n"); fail = 1; }
        if (is_loeschian(1000037ULL * 1000037ULL * 1000037ULL)) {
            fprintf(stderr, "SELFTEST FAIL: bad^3 large\n"); fail = 1; }
        /* good*bad^2 and bad^2*bad'^2 are Loeschian */
        if (!is_loeschian(1000003ULL * 1000037ULL * 1000037ULL)) {
            fprintf(stderr, "SELFTEST FAIL: good*bad^2\n"); fail = 1; }
        /* 5) agreement with a brute-force x^2+xy+y^2 enumeration up to 20000 */
        for (u64 t = 0; t <= 20000; t++) {
            bool bf = false;
            for (u64 x = 0; x * x <= t && !bf; x++)
                for (u64 y = 0; x * x + x * y + y * y <= t; y++)
                    if (x * x + x * y + y * y == t) { bf = true; break; }
            if (bf != is_loeschian(t)) {
                fprintf(stderr, "SELFTEST FAIL: brute force disagrees at %"
                        PRIu64 "\n", t); fail = 1; break; }
        }
        fprintf(stderr, "selftest: %s\n", fail ? "FAILED" : "OK");
        return fail;
    }

    FILE *of = out ? fopen(out, "a") : NULL;
    if (out && !of) { perror("open --out"); return 2; }

    int global_best = 0;
    u64 gb_a = 0, gb_d = 0;
    double t0 = now_s();
    u64 units = 0;
    double covered = 0;       /* raw a-range covered */

    for (u64 K = kmin; K <= kmax && !stop_requested; K++) {
        u64 d = K * D0;
        if (d / D0 != K) { fprintf(stderr, "K overflow at %" PRIu64 "\n", K); break; }

        /* ---- stage 1 components. Each is (modulus, start, step, count): the
         * good residues of a mod m are start + i*step (mod m), i < count, so
         * the CRT-admissible set is a product of arithmetic progressions and
         * can be walked by additions alone.
         *
         *   m = 3 : a = 1        (Loeschian numbers are never 2 mod 3)
         *   m = 2 : a = 1        (2 | d, so a must be odd)
         *   m = 5 : a != 0       (5 | d; all four nonzero residues are good)
         *   m = q : a = j*d for j = 1 .. q-n, for the smallest bad q > n
         *           (then q | a+i*d only at i = q-j >= n, outside the window)
         *
         * Putting 2 and 5 here rather than in the tier-C bitmask is worth
         * ~1.8x: candidates with a even or a = 0 mod 5 are then never
         * constructed at all, instead of being built and then rejected. */
        u64 cm[MAXTIERB + 4], cs[MAXTIERB + 4], ct[MAXTIERB + 4], cc[MAXTIERB + 4];
        int nc = 0;
        u64 MOD = 1;
        #define ADDC(m_, s_, t_, c_) do { cm[nc] = (m_); cs[nc] = (s_); \
            ct[nc] = (t_); cc[nc] = (c_); MOD *= (m_); nc++; } while (0)
        ADDC(3, 1, 0, 1);
        if (D0 % 2 == 0) ADDC(2, 1, 0, 1);
        if (D0 % 5 == 0) ADDC(5, 1, 1, 4);
        int ntb = 0;
        u64 tbq[MAXTIERB];
        for (u64 i = 0; i < n_sp && nc < MAXTIERB + 3; i++) {
            u64 q = sp[i];
            if (q <= nterms || q % 3 != 2 || d % q == 0) continue;
            if (MOD > (u64)4e18 / q) break;          /* keep MOD in u64 */
            if (MOD * q > modcap) break;             /* --modcap: unit size */
            u64 dq = d % q;
            tbq[ntb++] = q;
            ADDC(q, dq, dq, q - nterms);             /* a = j*d, j = 1..q-n */
        }
        #undef ADDC
        if (ntb < 4) { fprintf(stderr, "K=%" PRIu64 ": too few tier-B primes, "
                               "skipping\n", K); continue; }

        /* ---- CRT: idempotents, base residue R0, additive steps ---- */
        u64 R0 = 0, sstep[MAXTIERB + 4], subcyc[MAXTIERB + 4];
        for (int i = 0; i < nc; i++) {
            u64 mi = cm[i], co = MOD / mi;
            u64 e = mulmod(co % MOD, inv_mod(co % mi, mi), MOD); /* idempotent */
            R0 = (R0 + mulmod(cs[i], e, MOD)) % MOD;
            sstep[i] = mulmod(ct[i], e, MOD);
        }
        /* one full cycle of component i, to be undone when it wraps */
        for (int i = 0; i < nc; i++) subcyc[i] = mulmod(cc[i], sstep[i], MOD);

        /* Verify the CRT construction rather than trusting it: R0 must hit the
         * first good residue of every component, and sstep[i] must move ONLY
         * component i. A silent error here produces residues that are not
         * admissible at all, which looks like "the search found nothing". */
        for (int i = 0; i < nc; i++) {
            if (R0 % cm[i] != cs[i] % cm[i]) {
                fprintf(stderr, "FATAL CRT: R0=%" PRIu64 " mod %" PRIu64 " = %"
                        PRIu64 ", want %" PRIu64 "\n",
                        R0, cm[i], R0 % cm[i], cs[i] % cm[i]); return 2; }
            for (int j = 0; j < nc; j++) {
                u64 got = sstep[i] % cm[j], want2 = (i == j) ? ct[i] % cm[j] : 0;
                if (got != want2) {
                    fprintf(stderr, "FATAL CRT: sstep[%d] mod %" PRIu64 " = %"
                            PRIu64 ", want %" PRIu64 "\n", i, cm[j], got, want2);
                    return 2; }
            }
        }

        /* ---- tier C: bitmask primes ---- */
        u64 tcp[MAXTIERC]; bar tcb[MAXTIERC]; u64 *tcw[MAXTIERC]; int ntc = 0;
        for (u64 i = 0; i < n_sp && ntc < MAXTIERC; i++) {
            u64 r = sp[i];
            if (r > b2) break;
            bool inB = false;
            for (int t = 0; t < nc; t++) if (cm[t] == r) inB = true;
            if (inB) continue;                       /* already pinned in stage 1 */
            bool divD = (D0 % r == 0);
            if (!divD && (r % 3 != 2 || d % r == 0 || r <= nterms)) continue;
            tcp[ntc] = r;
            tcw[ntc] = malloc(r * sizeof(u64));
            ntc++;
        }

        /* Order tier C by how much of the candidate space each prime kills, so
         * the AND chain hits zero in a handful of steps. A bad prime q kills
         * nterms/q of residues; a prime p | D0 kills only 1/p. Ascending prime
         * order would put the feeble ones (11, 17, 23, ...) first and plough
         * through the whole chain almost every time. */
        for (int i = 0; i < ntc; i++) {
            int bestj = i;
            double bk = -1;
            for (int j = i; j < ntc; j++) {
                double kill = (D0 % tcp[j] == 0) ? 1.0 / (double)tcp[j]
                                                 : (double)nterms / (double)tcp[j];
                if (kill > bk) { bk = kill; bestj = j; }
            }
            u64 t1 = tcp[i]; tcp[i] = tcp[bestj]; tcp[bestj] = t1;
            u64 *t2 = tcw[i]; tcw[i] = tcw[bestj]; tcw[bestj] = t2;
        }
        for (int i = 0; i < ntc; i++) tcb[i] = bar_make(tcp[i]);
        if (K == kmin)
            fprintf(stderr, "tierC order: %" PRIu64 " %" PRIu64 " %" PRIu64
                    " %" PRIu64 " %" PRIu64 " ... (%d primes)\n",
                    tcp[0], tcp[1], tcp[2], tcp[3], tcp[4], ntc);

        for (u64 shift = 0; shift < shifts && !stop_requested; shift++) {
            u64 boff = shift * 64;

            /* build OKOK words for this (K, shift) */
            for (int t = 0; t < ntc; t++) {
                u64 r = tcp[t];
                char *ok = calloc(r, 1);
                if (D0 % r == 0) {
                    for (u64 x = 0; x < r; x++) ok[x] = (x != 0);
                } else {
                    for (u64 x = 0; x < r; x++) ok[x] = 1;
                    /* bad residues: a = -k*d (mod r), k = 0..nterms-1 */
                    u64 v = 0, dm = d % r;
                    for (u64 k = 0; k < nterms; k++) {
                        ok[v] = 0;
                        v = (v + r - dm) % r;
                    }
                }
                u64 modr = MOD % r;
                for (u64 x = 0; x < r; x++) {
                    u64 w = 0, y = (x + mulmod(boff % r, modr, r)) % r;
                    for (int b = 0; b < 64; b++) {
                        if (ok[y]) w |= (u64)1 << b;
                        y += modr; if (y >= r) y -= r;
                    }
                    tcw[t][x] = w;
                }
                free(ok);
            }

            /* ---- stage 1: nested additive enumeration of admissible R ---- */
            /* counts: component 0 is mod 3 (single residue); tier B follows.
             * The innermost loop gets the largest count. */
            u64 idx[MAXTIERB + 1];
            for (int i = 0; i < nc; i++) idx[i] = 0;
            u64 R = R0;
            /* (component counts drive the odometer directly) */


            u64 nres = 0, nsurv = 0, nconf = 0;
            double ts = now_s();
            {   u64 tot = 1;
                for (int i = 0; i < nc; i++) tot *= cc[i];
                fprintf(stderr, "unit K=%" PRIu64 " shift=%" PRIu64 " MOD=%"
                        PRIu64 " stage1=%d(", K, shift, MOD, nc);
                for (int i = 0; i < nc; i++)
                    fprintf(stderr, "%s%" PRIu64, i ? "," : "", cm[i]);
                fprintf(stderr, ") tierC=%d residues=%.4g range=%.4g\n",
                        ntc, (double)tot, (double)MOD * 64.0);
            }

            /* odometer over components 1..nc-1 */
            for (;;) {
                if (stop_requested) break;
                nres++;
                /* periodic progress: a work unit can be billions of residues,
                 * and a run with nothing to say for an hour is unreadable. */
                if ((nres & 0x3FFFFFF) == 0) {
                    double el2 = now_s() - ts;
                    fprintf(stderr, "  .. K=%" PRIu64 " shift=%" PRIu64 " R=%.4g"
                            " surv=%" PRIu64 " conf=%" PRIu64 " best=%d"
                            " %.0fs (%.3g res/s)\n", K, shift, (double)nres,
                            nsurv, nconf, global_best, el2, nres / el2);
                    if (of) {
                        fprintf(of, "{\"progress\":true,\"K\":%" PRIu64
                                ",\"shift\":%" PRIu64 ",\"res\":%" PRIu64
                                ",\"surv\":%" PRIu64 ",\"conf\":%" PRIu64
                                ",\"best\":%d,\"secs\":%.1f}\n",
                                K, shift, nres, nsurv, nconf, global_best, el2);
                        fflush(of);
                    }
                }

                /* ---- stage 2: 64 candidates at once ---- */
                u64 sito = ~(u64)0;
                for (int t = 0; t < ntc; t++) {
                    sito &= tcw[t][bar_mod(R, tcb[t])];
                    if (!sito) break;
                }
                while (sito) {
                    int b = __builtin_ctzll(sito);
                    sito &= sito - 1;
                    nsurv++;
                    u64 a0 = R + (boff + (u64)b) * MOD;
                    /* ---- stage 3: exact test of the whole admissible window.
                     * The LONGEST run need not start at index 0, so test all
                     * nterms and take the best run; then push outward past the
                     * window ends, which is where a 58 can hide next to a 56.
                     * Stage 3 fires ~1e-5 of residues, so this is free. ---- */
                    u64 best_run = 0, best_start = 0, cur = 0, cur_start = 0;
                    u64 want = (report > 1) ? (u64)report : 1;
                    for (u64 k = 0; k < nterms; k++) {
                        /* Give up on this window as soon as even a perfect
                         * tail cannot reach the reporting threshold. Costs
                         * ~19 term tests instead of 58 and recovers most of
                         * the throughput that full-window scanning spent. */
                        if (cur + (nterms - k) < want) break;
                        u64 t2 = a0 + k * d;
                        if (t2 < a0) break;              /* u64 overflow guard */
                        if (is_loeschian(t2)) {
                            if (cur == 0) cur_start = k;
                            cur++;
                            if (cur > best_run) { best_run = cur; best_start = cur_start; }
                        } else cur = 0;
                    }
                    if (best_run == 0) continue;
                    nconf++;
                    /* extend outward from the best run */
                    u64 a = a0 + best_start * d;
                    u64 run = best_run;
                    while (a >= d && is_loeschian(a - d)) { a -= d; run++; }
                    for (;;) {
                        u64 t2 = a + run * d;
                        if (t2 < a || !is_loeschian(t2)) break;
                        run++;
                    }
                    if ((int)run >= report || (int)run > global_best) {
                        if ((int)run > global_best) {
                            global_best = (int)run; gb_a = a; gb_d = d;
                            fprintf(stderr, "*** n=%d a=%" PRIu64 " d=%" PRIu64
                                    " (K=%" PRIu64 " shift=%" PRIu64 ")\n",
                                    global_best, a, d, K, shift);
                        }
                        if (of && (int)run >= report) {
                            fprintf(of, "{\"hit\":true,\"n\":%" PRIu64 ",\"a\":%"
                                    PRIu64 ",\"d\":%" PRIu64 ",\"K\":%" PRIu64
                                    "}\n", run, a, d, K);
                            fflush(of);
                        }
                    }
                }

                /* Advance the odometer. A component's cycle does NOT close
                 * (tbc[i]*sstep[i] is not 0 mod MOD), so on wrap we subtract
                 * the whole cycle explicitly to reset that component. */
                int i = nc - 1;
                bool exhausted = false;
                for (;;) {
                    while (i >= 0 && cc[i] <= 1) i--;   /* fixed components */
                    if (i < 0) { exhausted = true; break; }
                    R += sstep[i]; if (R >= MOD) R -= MOD;
                    idx[i]++;
                    if (idx[i] < cc[i]) break;
                    R = (R + MOD - subcyc[i]) % MOD;    /* undo cc[i] steps */
                    idx[i] = 0;
                    i--;
                }
                if (exhausted) break;
            }

            covered += (double)MOD * 64.0;
            units++;
            double el = now_s() - ts;
            fprintf(stderr, "K=%" PRIu64 " shift=%" PRIu64 " MOD=%" PRIu64
                    " tierB=%d tierC=%d R=%" PRIu64 " surv=%" PRIu64
                    " conf=%" PRIu64 " best=%d %.1fs (%.3g raw/s)\n",
                    K, shift, MOD, ntb, ntc, nres, nsurv, nconf, global_best,
                    el, (double)MOD * 64.0 / el);
            if (of) {
                fprintf(of, "{\"K\":%" PRIu64 ",\"shift\":%" PRIu64 ",\"MOD\":%"
                        PRIu64 ",\"res\":%" PRIu64 ",\"surv\":%" PRIu64
                        ",\"conf\":%" PRIu64 ",\"best\":%d,\"secs\":%.2f,"
                        "\"covered\":%.6g}\n",
                        K, shift, MOD, nres, nsurv, nconf, global_best, el, covered);
                fflush(of);
            }
        }
        for (int t = 0; t < ntc; t++) free(tcw[t]);
    }

    if (of) fclose(of);
    printf("BEST n=%d a=%" PRIu64 " d=%" PRIu64 "  units=%" PRIu64
           " covered=%.4g raw a  %.1fs%s\n", global_best, gb_a, gb_d, units,
           covered, now_s() - t0, stop_requested ? " INTERRUPTED" : "");
    return 0;
}

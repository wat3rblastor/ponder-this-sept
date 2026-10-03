/* loesch_core.h -- shared Loeschian arithmetic for the searchers.
 *
 * Extracted from apsearch.c so the CPU searcher (apsearch.c) and the Metal
 * searcher (apsearch_gpu.m) use ONE copy of the exact Loeschian test and the
 * modular helpers. Divergence between two copies of is_loeschian() is exactly
 * the kind of bug that would invalidate a record silently.
 */
#ifndef LOESCH_CORE_H
#define LOESCH_CORE_H

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

/* Montgomery arithmetic for odd n < 2^64. The u128 % in mulmod() is a library
 * call (__umodti3); stage 3 became the bottleneck of the GPU engine on hosts
 * with slower cores, and Miller-Rabin / Pollard are nothing but modular
 * multiplications. mg_mul is two widening multiplies and no division. */
typedef struct { u64 n, ni, one, r2; } mg_t;      /* ni = n^-1 mod 2^64 */
static inline mg_t mg_make(u64 n) {
    mg_t M; M.n = n;
    u64 x = n;                                    /* Newton: 3 correct bits -> 96 */
    x *= 2 - n * x; x *= 2 - n * x; x *= 2 - n * x; x *= 2 - n * x; x *= 2 - n * x;
    M.ni = x;
    M.one = (0 - n) % n;                          /* 2^64 mod n */
    M.r2 = (u64)((u128)M.one * M.one % n);        /* 2^128 mod n */
    return M;
}
static inline u64 mg_mul(u64 a, u64 b, const mg_t *M) {
    u128 t = (u128)a * b;
    u64 q = (u64)t * M->ni;
    u64 hi = (u64)(t >> 64), mh = (u64)(((u128)q * M->n) >> 64);
    return hi >= mh ? hi - mh : hi - mh + M->n;
}

/* Deterministic Miller-Rabin for odd n > 37 below 2^64 (7-base set of Sinclair). */
static bool mg_is_prime(u64 n, const mg_t *M) {
    static const u64 wit[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    const u64 one = M->one, mone = n - one;
    for (int i = 0; i < 7; i++) {
        u64 a = wit[i] % n;
        if (!a) continue;
        u64 b = mg_mul(a, M->r2, M), x = one, e = d;
        while (e) { if (e & 1) x = mg_mul(x, b, M); b = mg_mul(b, b, M); e >>= 1; }
        if (x == one || x == mone) continue;
        int ok = 0;
        for (int j = 1; j < s; j++) { x = mg_mul(x, x, M);
            if (x == mone) { ok = 1; break; } }
        if (!ok) return false;
    }
    return true;
}

/* Pollard-Brent with a batched gcd, in Montgomery form (n odd, composite).
 * x -> x^2 + c on Montgomery representatives is the same kind of random map;
 * differences and their product keep their gcd with n because 2^64 is a unit. */
static u64 mg_pollard(u64 n, const mg_t *M) {
    for (u64 c = 1;; c++) {
        u64 y = 2, m = 128, g = 1, r = 1, q = M->one, x = 0, ys = 0;
        #define MG_STEP(y_) do { y_ = mg_mul(y_, y_, M) + c; if (y_ < c || y_ >= n) y_ -= n; } while (0)
        while (g == 1) {
            x = y;
            for (u64 i = 0; i < r; i++) MG_STEP(y);
            for (u64 k = 0; k < r && g == 1; k += m) {
                ys = y;
                u64 lim = (m < r - k) ? m : r - k;
                for (u64 i = 0; i < lim; i++) {
                    MG_STEP(y);
                    q = mg_mul(q, x > y ? x - y : y - x, M);
                }
                g = gcd_u64(q, n);
            }
            r *= 2;
        }
        if (g == n) {                     /* back off one step at a time */
            g = 1;
            y = ys;
            while (g == 1) {
                MG_STEP(y);
                g = gcd_u64(x > y ? x - y : y - x, n);
            }
        }
        #undef MG_STEP
        if (g != n) return g;
    }
}

static u32 *sp = NULL;        /* small primes for trial division */
static u64 *sp_inv = NULL;    /* p^-1 mod 2^64 (odd p): t divisible by p  <=>  t * inv <= sp_lim */
static u64 *sp_lim = NULL;    /* floor((2^64 - 1) / p) */
static u64 n_sp = 0, sp_max = 0;

static void build_small_primes(u64 limit) {
    char *c = calloc(limit + 1, 1);
    for (u64 i = 2; i * i <= limit; i++)
        if (!c[i]) for (u64 j = i * i; j <= limit; j += i) c[j] = 1;
    u64 cnt = 0;
    for (u64 i = 2; i <= limit; i++) if (!c[i]) cnt++;
    sp = malloc(cnt * sizeof(u32));
    sp_inv = malloc(cnt * sizeof(u64));
    sp_lim = malloc(cnt * sizeof(u64));
    for (u64 i = 2; i <= limit; i++) if (!c[i]) {
        u64 x = i;
        x *= 2 - i * x; x *= 2 - i * x; x *= 2 - i * x; x *= 2 - i * x; x *= 2 - i * x;
        sp_inv[n_sp] = x; sp_lim[n_sp] = ~(u64)0 / i;
        sp[n_sp++] = (u32)i;
    }
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
static __thread int nfac;      /* thread-local: stage 3 runs under OpenMP in the CUDA engine */
static __thread u64 fac[64];

static void factor_rec(u64 n) {     /* n odd, every prime factor > 37 */
    if (n == 1) return;
    mg_t M = mg_make(n);
    if (mg_is_prime(n, &M)) { fac[nfac++] = n; return; }
    u64 r = isqrt_u64(n);
    if (r * r == n) { factor_rec(r); factor_rec(r); return; }
    /* perfect cube (pollard-rho is unreliable on prime powers) */
    u64 c = (u64)__builtin_cbrtl((long double)n);
    for (u64 cc = (c > 1 ? c - 1 : 1); cc <= c + 1; cc++)
        if (cc > 1 && cc * cc * cc == n) {
            factor_rec(cc); factor_rec(cc); factor_rec(cc); return;
        }
    u64 f = mg_pollard(n, &M);
    factor_rec(f);
    factor_rec(n / f);
}

/* True iff t is Loeschian: every prime p = 2 (mod 3) divides t to an even
 * power. Exact for t < 2^64 -- this is the searcher's own test; src/verify.py
 * re-derives any hit independently in Python ints.
 * Requires build_small_primes(limit) with limit >= 37. */
static bool is_loeschian(u64 t) {
    if (t == 0) return true;
    {   int e = __builtin_ctzll(t);
        if (e & 1) return false;
        t >>= e; }
    for (u64 i = 1; i < n_sp; i++) {             /* sp[0] = 2 is done */
        u32 p = sp[i];
        if ((u64)p * p > t) {
            /* what is left is 1 or a single prime */
            return (t == 1) ? true : (t % 3 != 2);
        }
        u64 q = t * sp_inv[i];
        if (q > sp_lim[i]) continue;             /* p does not divide t */
        int e = 0;
        do { t = q; e++; q = t * sp_inv[i]; } while (q <= sp_lim[i]);
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


#endif /* LOESCH_CORE_H */

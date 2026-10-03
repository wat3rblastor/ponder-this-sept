/* Equivalence + speed test: new is_loeschian (src/c/loesch_core.h) vs the reference.
 * Build: gcc -O3 -march=native -fopenmp -std=gnu11 -I src/c -o build/equiv \
 *        experiments/2026-10-03-profiling/stage3/equiv.c experiments/2026-10-03-profiling/stage3/ref.c -lm
 * Run:   build/equiv [nrandom]     exits 1 on the first mismatch */
#define _POSIX_C_SOURCE 200809L
#include <signal.h>
#include <time.h>
#include <omp.h>
#include "loesch_core.h"
void ref_init(u64 l); bool ref_is(u64 t);
static u64 rng(u64 *s) { u64 z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
static long bad = 0;
static void chk(u64 t) { if (is_loeschian(t) != ref_is(t)) {
    #pragma omp critical
    { bad++; fprintf(stderr, "MISMATCH t=%llu new=%d ref=%d\n", (unsigned long long)t, is_loeschian(t), ref_is(t)); } } }
int main(int argc, char **argv) {
    u64 N = argc > 1 ? strtoull(argv[1], 0, 10) : 2000000;
    (void)on_signal; (void)now_s; (void)inv_mod; (void)pollard; (void)is_prime_u64;
    build_small_primes(10000); ref_init(10000);
    /* 1. every t up to 3e6 */
    #pragma omp parallel for schedule(dynamic, 4096)
    for (u64 t = 0; t <= 3000000; t++) chk(t);
    printf("small range 0..3e6: mismatches so far %ld\n", bad);
    /* 2. structured: products / powers of primes just above the trial bound and near 2^32, edges of 2^64 */
    u64 ps[4000]; int np = 0;
    for (u64 x = 10001; np < 1500; x += 2) if (is_prime_u64(x)) ps[np++] = x;
    for (u64 x = 4294967295ULL; np < 2500; x -= 2) if (is_prime_u64(x)) ps[np++] = x;
    for (u64 x = 2642245; np < 3000; x -= 2) if (is_prime_u64(x)) ps[np++] = x;   /* cube root of 2^64 */
    for (u64 x = 65521; np < 3400; x -= 2) if (is_prime_u64(x)) ps[np++] = x;      /* fourth root */
    #pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < np; i++) {
        u64 p = ps[i]; u128 v;
        v = (u128)p * p; if (v >> 64 == 0) chk((u64)v);
        v = (u128)p * p * p; if (v >> 64 == 0) chk((u64)v);
        v = (u128)p * p * p * p; if (v >> 64 == 0) chk((u64)v);
        for (int j = i + 1; j < np && j < i + 60; j++) {
            v = (u128)p * ps[j]; if (v >> 64 == 0) { chk((u64)v); if (((u128)(u64)v * 4) >> 64 == 0) chk((u64)v * 4); }
            v = (u128)p * p * ps[j]; if (v >> 64 == 0) chk((u64)v);
            v = (u128)p * ps[j] * ps[(j * 7 + 3) % np]; if (v >> 64 == 0) chk((u64)v);
        }
    }
    for (u64 k = 0; k < 20000; k++) { chk(~(u64)0 - k); chk(((u64)1 << 63) + k - 10000); }
    printf("structured: mismatches so far %ld\n", bad);
    /* 3. random, the search's shape: t = 1 mod 3 across 1e15..2^64, and unrestricted */
    double tn = 0, tr = 0;
    #pragma omp parallel reduction(+:tn, tr)
    {
        u64 s = 12345 + 977 * omp_get_thread_num();
        #pragma omp for schedule(dynamic, 256)
        for (u64 i = 0; i < N; i++) {
            u64 r = rng(&s); int sh = rng(&s) % 15;
            u64 t = r >> sh; if (i & 1) t = t - t % 3 + 1;
            struct timespec a, b, c;
            clock_gettime(CLOCK_THREAD_CPUTIME_ID, &a); bool x = is_loeschian(t);
            clock_gettime(CLOCK_THREAD_CPUTIME_ID, &b); bool y = ref_is(t);
            clock_gettime(CLOCK_THREAD_CPUTIME_ID, &c);
            tn += (b.tv_sec - a.tv_sec) + 1e-9 * (b.tv_nsec - a.tv_nsec);
            tr += (c.tv_sec - b.tv_sec) + 1e-9 * (c.tv_nsec - b.tv_nsec);
            if (x != y) {
                #pragma omp critical
                { bad++; fprintf(stderr, "MISMATCH t=%llu new=%d ref=%d\n", (unsigned long long)t, x, y); }
            }
        }
    }
    printf("random %llu: new %.2f us/call, ref %.2f us/call, speedup %.2fx\n", (unsigned long long)N,
           1e6 * tn / N, 1e6 * tr / N, tr / tn);
    printf("TOTAL MISMATCHES %ld\n", bad);
    return bad ? 1 : 0;
}

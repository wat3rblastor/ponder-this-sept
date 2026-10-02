/* loeschsearch -- exhaustive search for long arithmetic progressions of
 * Loeschian numbers inside a sieved window [0, N].
 *
 * Build:  cc -O3 -march=native -std=c11 -o build/loeschsearch src/c/loeschsearch.c
 *
 * Method
 * ------
 * 1. Sieve a bitmap L[0..N] by direct enumeration of x^2 + x*y + y^2 <= N for
 *    x, y >= 0 (O(N) pairs, monotone writes per x, so it streams).
 * 2. For each step d in the requested family, compute for EVERY a in [0, N)
 *    the length of the Loeschian run a, a+d, a+2d, ... while terms stay <= N.
 *    Done in one descending pass with a circular buffer of d bytes:
 *        run(a) = L[a] ? run(a+d) + 1 : 0
 *    so one d costs O(N) time and O(d) memory, and covers all a at once.
 * 3. Report the best (a, d, n) and the run-length histogram (the histogram is
 *    the calibration data: it measures P(run >= k) empirically).
 *
 * This is an EXHAUSTIVE statement: when it finishes d, no a in [0, N) with
 * all terms <= N was missed for that d. Runs are truncated at N, so reported
 * lengths are lower bounds that only ever understate.
 *
 * All terms considered here are <= N <= 2^34, far below any overflow cliff;
 * the authoritative check is still src/verify.py on Python ints.
 */

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>

typedef uint64_t u64;
typedef uint8_t u8;

static u64 *L = NULL;       /* bitmap, bit n set iff n is Loeschian */
static u64 Nmax = 0;        /* highest index represented */

static inline bool getbit(u64 n) { return (L[n >> 6] >> (n & 63)) & 1u; }
static inline void setbit(u64 n) { L[n >> 6] |= (u64)1 << (n & 63); }

static volatile sig_atomic_t stop_requested = 0;
static void on_signal(int s) { (void)s; stop_requested = 1; }

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* ---------------------------------------------------------------- sieve --- */

static void build_sieve(u64 N) {
    u64 words = (N >> 6) + 1;
    L = calloc(words, sizeof(u64));
    if (!L) { fprintf(stderr, "out of memory for %" PRIu64 " words\n", words); exit(2); }
    Nmax = N;
    for (u64 x = 0; x * x <= N; x++) {
        u64 xx = x * x;
        for (u64 y = 0;; y++) {
            u64 v = xx + x * y + y * y;
            if (v > N) break;
            setbit(v);
        }
    }
}

/* Self-test against the known prefix of A003136 (GOAL.md section 2b). */
static void sieve_selftest(void) {
    static const int prefix[] = {0, 1, 3, 4, 7, 9, 12, 13, 16, 19, 21, 25, 27,
                                 28, 31, 36, 37, 39, 43, 48, 49, 52, 57, 61,
                                 63, 64, 67, 73, 75, 76, 79, 81, 84, 91, 93,
                                 97, 100};
    const int np = (int)(sizeof(prefix) / sizeof(prefix[0]));
    if (Nmax < 100) { fprintf(stderr, "selftest needs N >= 100\n"); exit(2); }
    int j = 0;
    for (int n = 0; n <= 100; n++) {
        bool want = (j < np && prefix[j] == n);
        if (want) j++;
        if (getbit((u64)n) != want) {
            fprintf(stderr, "SIEVE SELFTEST FAILED at n=%d (got %d want %d)\n",
                    n, (int)getbit((u64)n), (int)want);
            exit(2);
        }
    }
    /* bad primes must be absent, their squares present */
    static const int bad[] = {2, 5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83};
    for (unsigned i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
        u64 p = (u64)bad[i];
        if (p <= Nmax && getbit(p)) {
            fprintf(stderr, "SIEVE SELFTEST FAILED: %" PRIu64 " marked Loeschian\n", p);
            exit(2);
        }
        if (p * p <= Nmax && !getbit(p * p)) {
            fprintf(stderr, "SIEVE SELFTEST FAILED: %" PRIu64 "^2 not marked\n", p);
            exit(2);
        }
    }
    fprintf(stderr, "sieve selftest: OK\n");
}

/* ------------------------------------------------------------- scanning --- */

#define MAXHIST 256

typedef struct {
    int best_n;
    u64 best_a;
    u64 min_a;           /* smallest a whose run reaches nmin_g (0 = none) */
    int found_min;
    u64 hist[MAXHIST];   /* hist[k] = #a whose run is exactly k */
} scan_result;

static int nmin_g = 0;   /* G1 objective: minimise a + (nmin_g-1)*d */

/* Longest Loeschian run for step d over all a in [0, N], terms <= N. */
static scan_result scan_d(u64 d) {
    scan_result r;
    r.best_n = 0;
    r.best_a = 0;
    r.min_a = 0;
    r.found_min = 0;
    memset(r.hist, 0, sizeof(r.hist));

    u8 *buf = calloc(d, 1);
    if (!buf) { fprintf(stderr, "out of memory for buf of %" PRIu64 "\n", d); exit(2); }

    /* descending a; slot index (a mod d) is also ((a+d) mod d), so the slot
     * holds run(a+d) when we arrive and run(a) when we leave. */
    u64 a = Nmax;
    u64 i = a % d;
    for (;;) {
        u8 prev = (a + d <= Nmax) ? buf[i] : 0;
        u8 cur = getbit(a) ? (u8)(prev < 255 ? prev + 1 : 255) : 0;
        buf[i] = cur;
        if (cur > r.best_n) { r.best_n = cur; r.best_a = a; }
        if (nmin_g && cur >= nmin_g) { r.min_a = a; r.found_min = 1; }
        r.hist[cur]++;
        if (a == 0) break;
        a--;
        i = (i == 0) ? d - 1 : i - 1;
    }
    free(buf);
    return r;
}

/* ----------------------------------------------------------------- main --- */

static void usage(const char *p) {
    fprintf(stderr,
      "usage: %s --N <sieve-limit> --dbase <D0> --mmin <m> --mmax <m>\n"
      "          [--mstep <s>] [--out <file.jsonl>] [--resume] [--hist-min <k>]\n"
      "\n"
      "Scans every step d = D0*m for m in [mmin, mmax], exhaustively over all a\n"
      "with every term <= N. Appends one JSON line per d to --out; --resume skips\n"
      "d values already present there.\n", p);
    exit(2);
}

int main(int argc, char **argv) {
    u64 N = 0, D0 = 0, mmin = 1, mmax = 1, mstep = 1;
    const char *out = NULL;
    bool resume = false;
    int hist_min = 10;

    for (int i = 1; i < argc; i++) {
        const char *k = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : (usage(argv[0]), ""))
        if (!strcmp(k, "--N")) N = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--dbase")) D0 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--mmin")) mmin = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--mmax")) mmax = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--mstep")) mstep = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--hist-min")) hist_min = atoi(NEXT());
        else if (!strcmp(k, "--nmin")) nmin_g = atoi(NEXT());
        else if (!strcmp(k, "--out")) out = NEXT();
        else if (!strcmp(k, "--resume")) resume = true;
        else usage(argv[0]);
        #undef NEXT
    }
    if (!N || !D0 || mmin < 1 || mmax < mmin || mstep < 1) usage(argv[0]);

    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);

    /* --resume: collect m values already done */
    u8 *done = calloc(mmax + 1, 1);
    if (!done) { fprintf(stderr, "oom\n"); return 2; }
    if (resume && out) {
        FILE *f = fopen(out, "r");
        if (f) {
            char line[4096];
            while (fgets(line, sizeof(line), f)) {
                const char *p = strstr(line, "\"m\":");
                if (p) {
                    u64 m = strtoull(p + 4, NULL, 10);
                    if (m <= mmax) done[m] = 1;
                }
            }
            fclose(f);
            fprintf(stderr, "resume: skipping already-done m values\n");
        }
    }

    double t0 = now_s();
    fprintf(stderr, "sieving Loeschian numbers up to N=%" PRIu64
            " (%.1f MiB bitmap)...\n", N, (double)(N >> 3) / (1 << 20));
    build_sieve(N);
    sieve_selftest();
    fprintf(stderr, "sieve built in %.1fs\n", now_s() - t0);

    FILE *of = NULL;
    if (out) {
        of = fopen(out, "a");
        if (!of) { perror("open --out"); return 2; }
    }

    int global_best = 0;
    u64 gb_a = 0, gb_d = 0;
    u64 best_last = 0;          /* smallest a+(nmin-1)d seen so far */
    u64 bl_a = 0, bl_d = 0;
    for (u64 m = mmin; m <= mmax && !stop_requested; m += mstep) {
        if (done[m]) continue;
        u64 d = D0 * m;
        if (d > N) { fprintf(stderr, "d=%" PRIu64 " exceeds N, stopping\n", d); break; }
        double ts = now_s();
        scan_result r = scan_d(d);
        double el = now_s() - ts;

        /* suffix-sum the histogram into "number of a with run >= k" */
        u64 ge[MAXHIST];
        ge[MAXHIST - 1] = r.hist[MAXHIST - 1];
        for (int k = MAXHIST - 2; k >= 0; k--) ge[k] = ge[k + 1] + r.hist[k];

        if (r.best_n > global_best) {
            global_best = r.best_n;
            gb_a = r.best_a;
            gb_d = d;
            fprintf(stderr, "*** new best n=%d  a=%" PRIu64 "  d=%" PRIu64 "\n",
                    global_best, gb_a, gb_d);
        }
        if (nmin_g && r.found_min) {
            u64 last = r.min_a + (u64)(nmin_g - 1) * d;
            if (!best_last || last < best_last) {
                best_last = last; bl_a = r.min_a; bl_d = d;
                fprintf(stderr, "*** G1 n=%d last=%" PRIu64 " a=%" PRIu64
                        " d=%" PRIu64 "\n", nmin_g, best_last, bl_a, bl_d);
                if (of) { fprintf(of, "{\"g1\":true,\"n\":%d,\"last\":%"
                          PRIu64 ",\"a\":%" PRIu64 ",\"d\":%" PRIu64 "}\n",
                          nmin_g, best_last, bl_a, bl_d); fflush(of); }
            }
        }
        fprintf(stderr, "m=%" PRIu64 " d=%" PRIu64 " best_n=%d a=%" PRIu64
                " (%.1fs)\n", m, d, r.best_n, r.best_a, el);

        if (of) {
            fprintf(of, "{\"m\":%" PRIu64 ",\"d\":%" PRIu64 ",\"N\":%" PRIu64
                    ",\"best_n\":%d,\"best_a\":%" PRIu64 ",\"secs\":%.2f,\"ge\":{",
                    m, d, N, r.best_n, r.best_a, el);
            bool first = true;
            for (int k = hist_min; k < MAXHIST; k++) {
                if (!ge[k]) break;
                fprintf(of, "%s\"%d\":%" PRIu64, first ? "" : ",", k, ge[k]);
                first = false;
            }
            fprintf(of, "}}\n");
            fflush(of);
        }
    }

    if (of) fclose(of);
    free(L);
    free(done);
    if (nmin_g)
        printf("G1 BEST last=%" PRIu64 " a=%" PRIu64 " d=%" PRIu64 " (n=%d)\n",
               best_last, bl_a, bl_d, nmin_g);
    printf("BEST n=%d a=%" PRIu64 " d=%" PRIu64 "  (total %.1fs%s)\n",
           global_best, gb_a, gb_d, now_s() - t0,
           stop_requested ? ", INTERRUPTED" : "");
    return 0;
}

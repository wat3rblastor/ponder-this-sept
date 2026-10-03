#define _POSIX_C_SOURCE 200809L
/* loeschclass -- long Loeschian APs with steps far too large to sieve [0,N].
 *
 * Build: cc -O3 -march=native -std=c11 -o build/loeschclass src/c/loeschclass.c
 *
 * WHY THIS EXISTS
 * ---------------
 * For an n-term AP every bad prime p < n must divide d (GOAL.md 2, 2b), so for
 * n = 58 we need d = 0 mod 3*2*5*11*17*23*29*41*47*53 = 382160924970 and the
 * terms exceed 2e13. A flat bitmap over [0, N] cannot go there. But since
 * d = 0 (mod D), EVERY term of the AP lies in the single class a (mod D), so we
 * only ever need Loeschian-ness along the progression
 *      T_j = r + j*D,   j = 0 .. J-1
 * for one admissible residue r. Steps are d = m*D, so in j-coordinates the AP
 * is just j, j+m, j+2m, ... and one row of J bits serves every m at once.
 *
 * THE KEY FACT (no factoring anywhere)
 * ------------------------------------
 * For t = 1 (mod 3):  t is Loeschian  <=>  every bad prime p <= sqrt(t) divides
 * t to an even power.
 *   Proof: t has at most one prime factor > sqrt(t), and it occurs to the first
 *   power. For t coprime to 3, t = (-1)^(number of bad prime factors, with
 *   multiplicity) (mod 3), so t = 1 (mod 3) forces that count to be even. If
 *   every bad p <= sqrt(t) has even exponent, the count they contribute is
 *   even, hence the lone large prime factor (if any) cannot be bad. []
 * So an EXACT Loeschian bitmap for the row costs one parity sieve, ~1.6*J
 * operations, with no factoring and no cofactor bookkeeping.
 *
 * Parity by inclusion with alternating signs: for each bad prime p and each
 * e = 1, 2, 3, ... add (-1)^(e+1) at the j with p^e | T_j. The sum telescopes
 * to 1 when v_p(T_j) is odd and 0 when it is even, so cnt[j] ends up as the
 * number of bad primes of odd valuation and T_j is Loeschian iff cnt[j] == 0.
 *
 * Correctness is cross-checked against src/verify.py (--dump).
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
typedef uint8_t u8;
typedef __uint128_t u128;

static volatile sig_atomic_t stop_requested = 0;
static void on_signal(int s) { (void)s; stop_requested = 1; }

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static u64 isqrt_u64(u64 n) {
    if (n == 0) return 0;
    u64 x = (u64)1 << ((64 - __builtin_clzll(n) + 1) / 2);
    for (;;) {
        u64 y = (x + n / x) / 2;
        if (y >= x) break;
        x = y;
    }
    while (x > 0 && x > n / x) x--;
    while ((x + 1) <= n / (x + 1)) x++;
    return x;
}

/* modular inverse of a mod m (m need not be prime; gcd(a,m) must be 1) */
static u64 inv_mod(u64 a, u64 m) {
    int64_t t = 0, newt = 1;
    int64_t r = (int64_t)m, newr = (int64_t)(a % m);
    while (newr != 0) {
        int64_t q = r / newr;
        int64_t tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (r != 1) return 0;            /* not invertible */
    if (t < 0) t += (int64_t)m;
    return (u64)t;
}

/* ------------------------------------------------------- prime machinery --- */

static u32 *bad_p = NULL;      /* bad primes (p % 3 == 2) not dividing D */
static u32 *bad_inv = NULL;    /* D^{-1} mod p, precomputed (independent of r) */
static u64 n_bad = 0;

static void build_primes(u64 pmax, u64 D) {
    u64 nb = pmax + 1;
    u8 *comp = calloc(nb, 1);
    if (!comp) { fprintf(stderr, "oom primes (%" PRIu64 ")\n", nb); exit(2); }
    for (u64 i = 2; i * i <= pmax; i++)
        if (!comp[i]) for (u64 j = i * i; j <= pmax; j += i) comp[j] = 1;

    u64 cap = 0;
    for (u64 p = 2; p <= pmax; p++) if (!comp[p] && p % 3 == 2 && D % p) cap++;
    bad_p = malloc(cap * sizeof(u32));
    bad_inv = malloc(cap * sizeof(u32));
    if (!bad_p || !bad_inv) { fprintf(stderr, "oom prime tables\n"); exit(2); }
    for (u64 p = 2; p <= pmax; p++) {
        if (comp[p] || p % 3 != 2 || D % p == 0) continue;
        bad_p[n_bad] = (u32)p;
        bad_inv[n_bad] = (u32)inv_mod(D % p, p);
        n_bad++;
    }
    free(comp);
    fprintf(stderr, "primes: %" PRIu64 " bad primes <= %" PRIu64
            " not dividing D\n", n_bad, pmax);
}

/* small primes dividing D, used to pick admissible residues r */
static u32 dp[64];
static int n_dp = 0;
static void factor_D_small(u64 D) {
    u64 x = D;
    for (u64 p = 2; p * p <= x; p++)
        if (x % p == 0) { dp[n_dp++] = (u32)p; while (x % p == 0) x /= p; }
    if (x > 1) dp[n_dp++] = (u32)x;
}

/* r is admissible iff r = 1 (mod 3) and r is coprime to every prime of D.
 * (3 | D, so "coprime to D's primes" already includes 3 ∤ r; we additionally
 *  need the residue 1 rather than 2 mod 3 -- Loeschian numbers are never 2.) */
static bool admissible_r(u64 r, u64 D) {
    (void)D;
    if (r % 3 != 1) return false;
    for (int i = 0; i < n_dp; i++) if (r % dp[i] == 0) return false;
    return true;
}

/* ------------------------------------------------------------- one row ---- */

static u8 *cnt = NULL;         /* odd-valuation counter, length J */
static u64 *bits = NULL;       /* Loeschian bitmap over j, J bits */
static u8 *runbuf = NULL;      /* circular run buffer, length Mmax */

static inline bool getbit(u64 j) { return (bits[j >> 6] >> (j & 63)) & 1u; }

/* Exact Loeschian bitmap for T_j = r + j*D, j < J. */
static void sieve_row(u64 r, u64 D, u64 J, u64 Tmax) {
    memset(cnt, 0, J);
    for (u64 i = 0; i < n_bad; i++) {
        u64 p = bad_p[i];
        /* e = 1: j = -r * D^{-1} (mod p) */
        u64 j0 = (u64)((u128)((p - r % p) % p) * bad_inv[i] % p);
        for (u64 j = j0; j < J; j += p) cnt[j]++;
        if (j0 >= J && p >= J) continue;     /* no hit at all: no higher power */
        /* e >= 2, alternating sign; only while p^e <= Tmax */
        u64 pe = p;
        int sign = -1;
        for (;;) {
            if (pe > Tmax / p) break;
            pe *= p;
            u64 iv = inv_mod(D % pe, pe);
            if (iv == 0) break;              /* cannot happen: gcd(D, p) = 1 */
            u64 je = (u64)((u128)((pe - r % pe) % pe) * iv % pe);
            if (je >= J) { if (pe > J) break; }
            for (u64 j = je; j < J; j += pe) cnt[j] = (u8)(cnt[j] + sign);
            sign = -sign;
        }
    }
    u64 words = (J >> 6) + 1;
    memset(bits, 0, words * sizeof(u64));
    for (u64 j = 0; j < J; j++)
        if (cnt[j] == 0) bits[j >> 6] |= (u64)1 << (j & 63);
}

/* Longest run for step m over all start j, staying inside [0, J). */
static int scan_m(u64 J, u64 m, u64 *best_j) {
    int best = 0;
    u64 i = (J - 1) % m;
    for (u64 j = J - 1;; j--) {
        u8 prev = (j + m < J) ? runbuf[i] : 0;
        u8 cur = getbit(j) ? (u8)(prev < 255 ? prev + 1 : 255) : 0;
        runbuf[i] = cur;
        if (cur > best) { best = cur; *best_j = j; }
        if (j == 0) break;
        i = (i == 0) ? m - 1 : i - 1;
    }
    return best;
}

/* ---------------------------------------------------------------- main ---- */

static void usage(const char *p) {
    fprintf(stderr,
      "usage: %s --D <step base> --J <row length> --rstart <r> --rows <count>\n"
      "          [--nmin <report threshold>] [--out <file.jsonl>] [--resume]\n"
      "          [--dump <j>:<m>]\n"
      "\n"
      "Searches a = r + j*D, d = m*D. D must be divisible by 3 and by every bad\n"
      "prime p < target n. Rows (= admissible residues r, scanned upward from\n"
      "--rstart) are disjoint search spaces, so parallel workers must be given\n"
      "non-overlapping --rstart/--rows windows.\n", p);
    exit(2);
}

int main(int argc, char **argv) {
    u64 D = 0, J = 0, rstart = 1, rows = 1;
    int nmin = 40;
    const char *out = NULL;
    bool resume = false;
    long dump_j = -1, dump_m = -1;

    for (int i = 1; i < argc; i++) {
        const char *k = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : (usage(argv[0]), ""))
        if (!strcmp(k, "--D")) D = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--J")) J = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--rstart")) rstart = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--rows")) rows = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--nmin")) nmin = atoi(NEXT());
        else if (!strcmp(k, "--out")) out = NEXT();
        else if (!strcmp(k, "--resume")) resume = true;
        else if (!strcmp(k, "--dump")) { const char *s = NEXT();
            dump_j = atol(s); const char *c = strchr(s, ':');
            dump_m = c ? atol(c + 1) : 1; }
        else usage(argv[0]);
        #undef NEXT
    }
    if (!D || !J || J < 64) usage(argv[0]);
    if (D % 3) { fprintf(stderr, "FATAL: 3 must divide D (Loeschian numbers are "
                         "never 2 mod 3)\n"); return 2; }
    if (nmin < 2) nmin = 2;

    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);

    factor_D_small(D);
    u64 Mmax = (J - 1) / (u64)(nmin - 1);
    if (Mmax < 1) { fprintf(stderr, "FATAL: J too small for nmin=%d\n", nmin); return 2; }

    /* Tmax over the whole run of rows */
    u64 rlast = rstart + 64 * rows + 64;   /* generous upper bound on r used */
    u64 Tmax = rlast + (J - 1) * D;
    u64 pmax = isqrt_u64(Tmax);
    fprintf(stderr, "D=%" PRIu64 " J=%" PRIu64 " Mmax=%" PRIu64
            " Tmax~%.4g sqrt(Tmax)=%" PRIu64 "\n",
            D, J, Mmax, (double)Tmax, pmax);

    double t0 = now_s();
    build_primes(pmax, D);
    fprintf(stderr, "prime setup %.2fs\n", now_s() - t0);

    cnt = malloc(J);
    bits = malloc(((J >> 6) + 1) * sizeof(u64));
    runbuf = malloc(Mmax + 1);
    if (!cnt || !bits || !runbuf) { fprintf(stderr, "oom row buffers\n"); return 2; }

    /* --resume: skip r values already recorded */
    u64 resume_from = 0;
    if (resume && out) {
        FILE *f = fopen(out, "r");
        if (f) {
            char line[4096];
            while (fgets(line, sizeof(line), f)) {
                const char *p = strstr(line, "\"r\":");
                if (p) { u64 rr = strtoull(p + 4, NULL, 10);
                         if (rr >= resume_from) resume_from = rr + 1; }
            }
            fclose(f);
            if (resume_from > rstart) {
                fprintf(stderr, "resume: restarting at r=%" PRIu64 "\n", resume_from);
                rstart = resume_from;
            }
        }
    }

    FILE *of = NULL;
    if (out) { of = fopen(out, "a"); if (!of) { perror("open --out"); return 2; } }

    /* --dump -1:<count>: print the verdict for EVERY j, for a full cross-check */
    if (dump_j == -1 && dump_m > 0) {
        sieve_row(rstart, D, J, Tmax);
        u64 lim = ((u64)dump_m < J) ? (u64)dump_m : J;
        printf("DUMPALL r=%" PRIu64 " D=%" PRIu64 " n=%" PRIu64 "\n", rstart, D, lim);
        for (u64 j = 0; j < lim; j++)
            printf("%" PRIu64 " %d\n", rstart + j * D, (int)getbit(j));
        return 0;
    }

    /* --dump: print the first terms of one AP so Python can cross-check */
    if (dump_j >= 0) {
        sieve_row(rstart, D, J, Tmax);
        printf("DUMP r=%" PRIu64 " D=%" PRIu64 " j=%ld m=%ld\n", rstart, D, dump_j, dump_m);
        u64 run = 0;
        for (u64 k = 0; (u64)dump_j + k * (u64)dump_m < J; k++) {
            u64 j = (u64)dump_j + k * (u64)dump_m;
            bool L = getbit(j);
            printf("%" PRIu64 " %d\n", rstart + j * D, (int)L);
            if (!L) break;
            run++;
            if (run > 200) break;
        }
        printf("RUN %" PRIu64 "\n", run);
        return 0;
    }

    int global_best = 0;
    u64 gb_a = 0, gb_d = 0;
    u64 done_rows = 0, cand = 0;
    double tstart = now_s();

    for (u64 r = rstart; done_rows < rows && !stop_requested; r++) {
        if (!admissible_r(r, D)) continue;
        done_rows++;
        sieve_row(r, D, J, Tmax);

        int row_best = 0; u64 rb_j = 0, rb_m = 0;
        for (u64 m = 1; m <= Mmax && !stop_requested; m++) {
            u64 bj = 0;
            int b = scan_m(J, m, &bj);
            cand += J;
            if (b > row_best) { row_best = b; rb_j = bj; rb_m = m; }
            if (b > global_best) {
                global_best = b; gb_a = r + bj * D; gb_d = m * D;
                fprintf(stderr, "*** n=%d  a=%" PRIu64 "  d=%" PRIu64
                        "  (r=%" PRIu64 " j=%" PRIu64 " m=%" PRIu64 ")\n",
                        global_best, gb_a, gb_d, r, bj, m);
                if (of) { fprintf(of, "{\"best\":true,\"n\":%d,\"a\":%" PRIu64
                                  ",\"d\":%" PRIu64 ",\"r\":%" PRIu64 "}\n",
                                  global_best, gb_a, gb_d, r); fflush(of); }
            }
        }
        if (of) {
            fprintf(of, "{\"r\":%" PRIu64 ",\"D\":%" PRIu64 ",\"J\":%" PRIu64
                    ",\"row_best_n\":%d,\"a\":%" PRIu64 ",\"d\":%" PRIu64
                    ",\"cand\":%" PRIu64 ",\"secs\":%.2f}\n",
                    r, D, J, row_best, r + rb_j * D, rb_m * D, cand, now_s() - tstart);
            fflush(of);
        }
        if (done_rows % 16 == 0 || done_rows == rows)
            fprintf(stderr, "rows=%" PRIu64 "/%" PRIu64 " r=%" PRIu64
                    " best=%d cand=%.4g (%.1f Mcand/s)\n",
                    done_rows, rows, r, global_best, (double)cand,
                    cand / (now_s() - tstart) / 1e6);
    }

    if (of) fclose(of);
    printf("BEST n=%d a=%" PRIu64 " d=%" PRIu64 " cand=%.4g rows=%" PRIu64
           " %.1fs%s\n", global_best, gb_a, gb_d, (double)cand, done_rows,
           now_s() - tstart, stop_requested ? " INTERRUPTED" : "");
    return 0;
}

/* Ground truth for the counting model.
 *
 * Loeschian bitmap for 0..X with X < 4472^2, using
 *   v_q(t) mod 2 = XOR over k>=1 of [q^k | t]   for bad q <= sqrt(X)
 * plus the single remaining cofactor (1 or one prime > sqrt(X)).
 *
 * Then, for the family  d = base(n)*m (base = 3*prod(bad q <= n/2)),
 * a = 1 mod 3, a != 0 mod (bad q <= n/2), last term <= X:
 * EXACT number of n-term Loeschian APs, and the exact per-term rate.
 *
 * cc -O2 -o joint joint.c && ./joint
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static long X = 19000000L;
static unsigned char *loe;

static int isbad(long p) { return p % 3 == 2; }

static void build(void) {
    long n = X + 1, S = 1; while ((S + 1) * (S + 1) <= X) S++;
    unsigned char *par = calloc(n, 1);
    long *res = malloc(n * sizeof(long));
    for (long i = 0; i < n; i++) res[i] = i;
    char *sv = malloc(S + 1); memset(sv, 1, S + 1); sv[0] = sv[1] = 0;
    for (long i = 2; i * i <= S; i++) if (sv[i]) for (long j = i * i; j <= S; j += i) sv[j] = 0;
    for (long p = 2; p <= S; p++) {
        if (!sv[p]) continue;
        /* bit 1 = scratch: parity of v_p; bit 0 = "some bad prime has odd v" */
        for (long pk = p; pk <= X; pk *= p) {
            for (long m = pk; m <= X; m += pk) {
                res[m] /= p;
                if (isbad(p)) par[m] ^= 2;
            }
            if (pk > X / p) break;
        }
        if (isbad(p))
            for (long m = p; m <= X; m += p) {
                if (par[m] & 2) par[m] |= 1;
                par[m] &= ~2;
            }
    }
    loe = malloc(n);
    for (long i = 0; i < n; i++) {
        unsigned char q = par[i] & 1;
        if (res[i] > 1 && isbad(res[i])) q |= 1;
        loe[i] = (q == 0);
    }
    loe[0] = 1;
    free(par); free(res); free(sv);
}

static long bad_list[16], nbad;
static long base_of(int n) {
    long b = 3; nbad = 0;
    for (long q = 2; 2 * q <= n; q++) {
        int pr = q > 1; for (long i = 2; i * i <= q; i++) if (q % i == 0) pr = 0;
        if (pr && isbad(q)) { b *= q; bad_list[nbad++] = q; }
    }
    return b;
}

int main(int argc, char **argv) {
    if (argc > 1) X = atol(argv[1]);
    build();
    long cnt = 0; for (long i = 1; i <= X; i++) cnt += loe[i];
    printf("first 50 Loeschian:"); for (long i=1;i<50;i++) if (loe[i]) printf(" %ld", i); printf("\n");
    printf("X=%ld  Loeschian density = %.5f\n", X, (double)cnt / X);
    int nlo = argc > 2 ? atoi(argv[2]) : 20, nhi = argc > 3 ? atoi(argv[3]) : 26;
    for (int n = nlo; n <= nhi; n += 2) {
        long base = base_of(n);
        long mmax = (X - 1) / ((long)(n - 1) * base);
        long long pairs = 0, terms = 0, ok = 0, found = 0;
        long fa = 0, fd = 0;
        for (long m = 1; m <= mmax; m++) {
            long d = base * m, hi = X - (long)(n - 1) * d;
            for (long a = 1; a <= hi; a += 3) {
                int skip = 0;
                for (long i = 0; i < nbad; i++) if (a % bad_list[i] == 0) { skip = 1; break; }
                if (skip) continue;
                pairs++; terms++; ok += loe[a];
                if (!loe[a]) continue;
                int k = 1; while (k < n && loe[a + (long)k * d]) k++;
                if (k == n) { found++; if (!fa) { fa = a; fd = d; } }
            }
        }
        printf("n=%2d base=%7ld #d=%7ld pairs=%12lld per-term rho=%.5f  EXACT=%lld",
               n, base, mmax, pairs, (double)ok / terms, found);
        if (fa) printf("  example a=%ld d=%ld", fa, fd);
        printf("\n"); fflush(stdout);
    }
    return 0;
}

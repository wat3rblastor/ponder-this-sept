/* Deep search of the *symmetric* quadratic case:
 *      P(k) = alpha*(2k-h)^2 + C      (h = 0..114, so vertex at h/2 in [0,57])
 * stored as integer coefficients; hit at k iff P(k) is a perfect square.
 * Reports (hits, Tlow = 2*max_{k in A} sqrt(P(k))).
 *
 * alpha up to AMAX, y up to YMAX, h over all 0..114 (h even => true symmetry).
 * For each (alpha,h) we tally C = y^2 - alpha*(2k-h)^2 over all (k,y).
 */
#include <stdio.h>
#include <stdlib.h>

#define NK 58

static unsigned char *cnt;
static int *maxy;
static long *touched;

int main(int argc, char **argv) {
    long AMAX = atol(argv[1]);
    long YMAX = atol(argv[2]);
    int REPORT = atoi(argv[3]);
    long lo = -(AMAX * 114L * 114L) - 1;
    long hi = YMAX * YMAX + 1;
    long span = hi - lo + 1;
    cnt = calloc(span, 1);
    maxy = malloc(span * sizeof(int));
    touched = malloc((long)NK * (YMAX + 1) * sizeof(long));
    if (!cnt || !maxy || !touched) { fprintf(stderr, "oom span=%ld\n", span); return 1; }
    int best = 0;
    for (long al = 1; al <= AMAX; al++) {
        for (long h = 0; h <= 114; h++) {
            long nt = 0;
            for (long k = 0; k < NK; k++) {
                long w = 2 * k - h;
                long q = al * w * w;
                for (long y = 0; y <= YMAX; y++) {
                    long g = y * y - q;
                    if (g < lo || g > hi) continue;
                    long idx = g - lo;
                    if (cnt[idx] == 0) { touched[nt++] = idx; maxy[idx] = y; }
                    if (cnt[idx] < 250) cnt[idx]++;
                    if (y > maxy[idx]) maxy[idx] = y;
                }
            }
            for (long i = 0; i < nt; i++) {
                long idx = touched[i];
                long C = idx + lo;
                int c = cnt[idx];
                cnt[idx] = 0;
                /* degenerate: P a perfect square polynomial <=> C==0 */
                if (C == 0) continue;
                if (c >= REPORT)
                    printf("hits=%d alpha=%ld h=%ld C=%ld Tlow=%d\n",
                           c, al, h, C, 2 * maxy[idx]);
                if (c > best) best = c;
            }
        }
        if (al % 100 == 0) fprintf(stderr, "alpha=%ld best=%d\n", al, best);
    }
    printf("BEST=%d\n", best);
    return 0;
}

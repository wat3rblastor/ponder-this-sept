/* Exhaustive search over quadratics P(k) = a*k^2 + b*k + g, counting
 *   hits(P) = #{ k in [0,57] : P(k) is a perfect square }.
 *
 * Why this is the right object: see README.  For any family in which the
 * terms t_k = a + k*d are *identically* of the form L1^2 + 3 L2^2 (the only
 * way a term can be Loeschian for free), a and d are binary quadratic forms
 * in the family parameters and
 *      P(k) = -disc(a + k d)/12
 * is a quadratic polynomial in k whose square values are exactly the
 * candidate automatic indices.  Leading coefficient alpha = -disc(d)/12 >= 0.
 *
 * Search is exhaustive over 0 <= alpha <= AMAX, |beta| <= BMAX, and all
 * gamma such that every hit has y = sqrt(P(k)) <= YMAX  (equivalently the
 * term-size lower bound T >= 2*max_{k in A} y is <= 2*YMAX).
 *
 * For each (alpha,beta) we tally gamma = y^2 - alpha k^2 - beta k over all
 * (k,y), so the count for every gamma is obtained in one pass.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NK 58

int main(int argc, char **argv) {
    long ALO = atol(argv[1]); long AHI = atol(argv[2]);
    long BMAX = atol(argv[3]);
    long YMAX = atol(argv[4]);
    long REPORT = atol(argv[5]);

    long AM = (AHI>-ALO)?AHI:-ALO;
    long lo = -(AM * 57 * 57 + BMAX * 57) - 1;
    long hi = YMAX * YMAX + AM*57*57 + BMAX*57 + 1;
    long span = hi - lo + 1;
    unsigned char *cnt = calloc(span, 1);
    long *maxy = malloc(span * sizeof(long));
    long *touched = malloc((long)NK * (YMAX + 1) * sizeof(long));
    if (!cnt || !maxy || !touched) { fprintf(stderr, "oom\n"); return 1; }

    int best = 0;
    for (long al = ALO; al <= AHI; al++) {
        for (long be = -BMAX; be <= BMAX; be++) {
            long nt = 0;
            for (long k = 0; k < NK; k++) {
                long q = al * k * k + be * k;
                for (long y = 0; y <= YMAX; y++) {
                    long g = y * y - q;
                    if (g < lo || g > hi) continue;
                    long idx = g - lo;
                    if (cnt[idx] == 0) touched[nt++] = idx;
                    if (cnt[idx] < 250) cnt[idx]++;
                    if (cnt[idx] == 1 || y > maxy[idx]) maxy[idx] = y;
                }
            }
            for (long i = 0; i < nt; i++) {
                long idx = touched[i];
                long g = idx + lo;
                /* skip P = (e k + f)^2 identically: that is the degenerate
                   pencil d = lambda*a, i.e. a rescaling of a solution. */
                if (be * be == 4 * al * g) { cnt[idx] = 0; continue; }
                if (cnt[idx] >= REPORT) {
                    printf("hits=%d alpha=%ld beta=%ld gamma=%ld Tmin=%ld\n",
                           cnt[idx], al, be, idx + lo, 2 * maxy[idx]);
                    fflush(stdout);
                }
                if (cnt[idx] > best) best = cnt[idx];
                cnt[idx] = 0;
            }
        }
        fprintf(stderr, "alpha=%ld done best=%d\n", al, best);
    }
    printf("BEST=%d\n", best);
    return 0;
}

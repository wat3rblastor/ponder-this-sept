/* mc_q -- measure q(t) = P(t Loeschian | t divisible by no bad prime <= B),
 * which is exactly the per-term pass probability inside an apsearch admissible
 * window (every term of such a window is free of every bad prime <= b2, and the
 * bad primes q > b2 each divide at most one of the 58 indices, so the 58 term
 * events are independent given the window).
 *
 * Build: cc -O3 -march=native -std=gnu11 -I src/c -o mc_q mc_q.c -lm
 * Usage: mc_q <B> <lo> <hi> <nsamples_target_accepted> <seed>
 * Prints: B lo hi accepted loeschian q sigma_q
 */
#include <signal.h>
#include <time.h>
#include <math.h>
#include "loesch_core.h"

static u64 rs;
static inline u64 rnd(void) {           /* splitmix64 */
    u64 z = (rs += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: mc_q B lo hi naccept seed\n"); return 2; }
    u64 B = strtoull(argv[1], NULL, 10);
    double lo = strtod(argv[2], NULL), hi = strtod(argv[3], NULL);
    u64 want = strtoull(argv[4], NULL, 10);
    rs = strtoull(argv[5], NULL, 10);
    build_small_primes(10000);

    /* bad primes <= B */
    u64 bad[512]; int nbad = 0;
    for (u64 i = 0; i < n_sp; i++) {
        if (sp[i] > B) break;
        if (sp[i] % 3 == 2) bad[nbad++] = sp[i];
    }

    u64 span = (u64)(hi - lo);
    u64 tried = 0, acc = 0, loe = 0;
    /* per-prefix-length run statistic over synthetic independent windows is
     * computed in python from q; here we only need q itself. */
    while (acc < want) {
        u64 t = (u64)lo + rnd() % span;
        tried++;
        int ok = 1;
        for (int i = 0; i < nbad; i++) if (t % bad[i] == 0) { ok = 0; break; }
        if (!ok) continue;
        acc++;
        if (is_loeschian(t)) loe++;
    }
    double q = (double)loe / (double)acc;
    printf("B=%" PRIu64 " lo=%.4g hi=%.4g tried=%" PRIu64 " accepted=%" PRIu64
           " loeschian=%" PRIu64 " q=%.6f sigma=%.6f acceptrate=%.6f\n",
           B, lo, hi, tried, acc, loe, q, sqrt(q * (1 - q) / (double)acc),
           (double)acc / (double)tried);
    return 0;
}

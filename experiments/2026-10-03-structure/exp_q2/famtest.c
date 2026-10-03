#define _POSIX_C_SOURCE 200809L
/* famtest -- empirical comparison of AP families for Loeschian APs.
 * usage: famtest FAMILY N MODE SECONDS [LMIN] [SIZEMULT]
 *   FAMILY: string of flags:  base | no53 | no47 | no41 | no29 | sq (square option for pins 59..113)
 *           combined with '+', e.g. "no53+sq", "no41+no47+no53"
 *   MODE: runs | full
 * Candidates are sampled UNIFORMLY from the family's admissible (a, K) with
 * a < Xa, d = K*Dbase < Xd, where the box is scaled so that every family has the
 * same number of admissible candidates as the baseline box (a < 64*M0, K <= 50).
 */
#include <signal.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "/workspace/ponder-this-sept/src/c/loesch_core.h"

static const int BADS[] = {2,5,11,17,23,29,41,47,53,59,71,83,89,101,107,113};
#define NB 16

static inline u64 rng(u64 *s) { u64 x = *s; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; *s = x; return x * 2685821657736338717ULL; }


typedef struct { u64 n, ni, r1, r2; } mg;
static inline mg mg_make(u64 n){ mg m; m.n=n; u64 x=n; for(int i=0;i<6;i++) x*=2-n*x; m.ni=-x; m.r1=(u64)((((u128)1)<<64)%n); m.r2=(u64)((u128)m.r1*m.r1%n); return m; }
static inline u64 mg_mul(u64 a,u64 b,const mg*m){ u128 T=(u128)a*b; u64 q=(u64)T*m->ni; u64 t=(u64)((T+(u128)q*m->n)>>64); return t>=m->n?t-m->n:t; }
static u64 pollard_m(u64 n) {
    if (!(n & 1)) return 2;
    mg M = mg_make(n);
    for (u64 c = 1;; c++) {
        u64 y = 2, m = 128, g = 1, r = 1, q = M.r1, x = 0, ys = 0;
        while (g == 1) {
            x = y;
            for (u64 i = 0; i < r; i++) { y = mg_mul(y, y, &M) + c; if (y >= n) y -= n; }
            for (u64 k = 0; k < r && g == 1; k += m) {
                ys = y;
                u64 lim = (m < r - k) ? m : r - k;
                for (u64 i = 0; i < lim; i++) {
                    y = mg_mul(y, y, &M) + c; if (y >= n) y -= n;
                    q = mg_mul(q, x > y ? x - y : y - x, &M);
                }
                g = gcd_u64(q, n);
            }
            r *= 2;
        }
        if (g == n) {
            g = 1; y = ys;
            while (g == 1) { y = mg_mul(y, y, &M) + c; if (y >= n) y -= n; g = gcd_u64(x > y ? x - y : y - x, n); }
        }
        if (g != n) return g;
    }
}
static bool prime7(u64 n) {
    if (n < 2) return false;
    static const u64 smallp[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) { if (n % smallp[i] == 0) return n == smallp[i]; }
    mg M = mg_make(n);
    u64 d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
    u64 one = M.r1, mone = n - M.r1;
    static const u64 wit[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022ULL};
    for (int i = 0; i < 7; i++) {
        u64 a = wit[i] % n; if (!a) continue;
        u64 b = mg_mul(a, M.r2, &M), x = one, e = d;
        while (e) { if (e & 1) x = mg_mul(x, b, &M); b = mg_mul(b, b, &M); e >>= 1; }
        if (x == one || x == mone) continue;
        int ok = 0;
        for (int j = 1; j < s; j++) { x = mg_mul(x, x, &M); if (x == mone) { ok = 1; break; } }
        if (!ok) return false;
    }
    return true;
}
static __thread int nf2; static __thread u64 f2[64];
static void frec(u64 n) {
    if (n == 1) return;
    if (prime7(n)) { f2[nf2++] = n; return; }
    u64 r = isqrt_u64(n);
    if (r * r == n) { frec(r); frec(r); return; }
    u64 c = (u64)__builtin_cbrtl((long double)n);
    for (u64 cc = (c > 1 ? c - 1 : 1); cc <= c + 1; cc++)
        if (cc > 1 && cc * cc * cc == n) { frec(cc); frec(cc); frec(cc); return; }
    u64 f = pollard_m(n); frec(f); frec(n / f);
}
/* exact Loeschian test */
static bool loe(u64 t) {
    if (t == 0) return true;
    for (u64 i = 0; i < n_sp; i++) {
        u32 p = sp[i];
        if ((u64)p * p > t) return (t == 1) ? true : (t % 3 != 2);
        if (t % p) continue;
        int e = 0; while (t % p == 0) { t /= p; e++; }
        if (p % 3 == 2 && (e & 1)) return false;
    }
    if (t == 1) return true;
    if (t % 3 == 2) return false;          /* odd number of bad prime factors => some odd exponent */
    if (prime7(t)) return true;
    nf2 = 0; frec(t);
    for (int i = 0; i < nf2; i++) {
        if (f2[i] % 3 != 2) continue;
        int e = 0; for (int j = 0; j < nf2; j++) if (f2[j] == f2[i]) e++;
        if (e & 1) return false;
    }
    return true;
}
static double tcpu(void){ struct timespec ts; clock_gettime(CLOCK_THREAD_CPUTIME_ID,&ts); return ts.tv_sec+1e-9*ts.tv_nsec; }
static inline int maxrun(u64 x) { int c = 0; while (x) { x &= x >> 1; c++; } return c; }

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage\n"); return 2; }
    const char *fam = argv[1]; int n = atoi(argv[2]); const char *mode = argv[3];
    double secs = atof(argv[4]); int LMIN = argc > 5 ? atoi(argv[5]) : 14;
    double sizemult = argc > 6 ? atof(argv[6]) : 1.0;
    bool full = !strcmp(mode, "full");
    bool ind[NB]; for (int i = 0; i < NB; i++) ind[i] = BADS[i] < 58;
    bool sq = strstr(fam, "sq") != NULL;
    if (strstr(fam, "no53")) ind[8] = false;
    if (strstr(fam, "no47")) ind[7] = false;
    if (strstr(fam, "no41")) ind[6] = false;
    if (strstr(fam, "no29")) ind[5] = false;
    build_small_primes(2000);
    /* densities */
    double rho = 1, rho0 = 1; u64 Dbase = 3;
    for (int i = 0; i < NB; i++) {
        int q = BADS[i];
        /* baseline */
        if (q < 58) rho0 *= (1.0 - 1.0 / q) / q; else rho0 *= (double)(q - n) / q;
        if (ind[i]) { rho *= (1.0 - 1.0 / q) / q; Dbase *= q; }
        else {
            double wa = q > n ? (double)(q - n) * q : 0;
            int lo = n - q > 0 ? n - q : 0, hi = (n < q ? n : q) - 1;
            double ws = (sq || q < n) ? (hi - lo + 1 > 0 ? hi - lo + 1 : 0) : 0;
            if (wa + ws <= 0) { fprintf(stderr, "family impossible at q=%d\n", q); return 1; }
            rho *= (wa + ws) / ((double)q * q);
        }
    }
    double rrel = rho / rho0;
    double M0 = 3.0 * 59 * 71 * 83 * 89 * 101 * 107 * 113;
    double Xa_d = 64 * M0 / sqrt(rrel) * sizemult, Xd_d = 50.0 * 382160924970.0 / sqrt(rrel) * sizemult;
    u64 Xa = (u64)Xa_d; u64 Kmax = (u64)(Xd_d / Dbase); if (Kmax < 1) Kmax = 1;
    /* sieve primes: bad primes < 2000 not dividing Dbase */
    int sv[400], nsv = 0;
    for (u64 i = 0; i < n_sp; i++) { u32 p = sp[i]; if (p % 3 == 2 && Dbase % p) sv[nsv++] = p; }
    int pins[NB], npin = 0, dps[NB], ndp = 0;
    for (int i = 0; i < NB; i++) { if (ind[i]) dps[ndp++] = BADS[i]; else pins[npin++] = BADS[i]; }
    fprintf(stderr, "fam=%s n=%d Dbase=%" PRIu64 " rho_rel=%.4f Xa=%.3e Kmax=%" PRIu64 " (n*dmax=%.3e)\n",
            fam, n, Dbase, rrel, Xa_d, Kmax, (double)n * Kmax * Dbase);

    u64 tot_c = 0, hist[64] = {0}, tot_tests = 0;
    u64 pass_idx[64] = {0}, cnt_full = 0, hit_n = 0, hit_pass = 0, non_n = 0, non_pass = 0, alive_sum = 0;
    double sumlog = 0; u64 nlog = 0;
    double t0 = now_s(); double cpu = 0;
    #pragma omp parallel num_threads(10)
    {
        u64 seed = 0x9E3779B97F4A7C15ULL * (omp_get_thread_num() + 1) ^ (u64)(strlen(fam) * 1315423911u + n * 77 + full);
        for (int i = 0; i < 20; i++) rng(&seed);
        u64 l_c = 0, l_hist[64] = {0}, l_tests = 0;
        u64 l_pass[64] = {0}, l_full = 0, l_hn = 0, l_hp = 0, l_nn = 0, l_np = 0, l_alive = 0;
        double l_sumlog = 0; u64 l_nlog = 0;
        u32 dinv[400];
        double c0 = tcpu();
        while (tcpu() - c0 < secs) {
            u64 K = 1 + rng(&seed) % Kmax;
            bool okK = true; for (int i = 0; i < npin; i++) if (K % pins[i] == 0) okK = false;
            if (!okK) continue;
            u64 d = K * Dbase;
            for (int i = 0; i < nsv; i++) { u64 q = sv[i]; dinv[i] = (u32)inv_mod(d % q, q); }
            for (int batch = 0; batch < 2000; batch++) {
                /* sample admissible a */
                u128 a = 1, M = 3; u64 hitmask = 0;
                for (int i = 0; i < npin; i++) {
                    u64 q = pins[i];
                    u64 wa = q > (u64)n ? (q - n) * q : 0;
                    int lo = n > (int)q ? n - (int)q : 0, hi = (n < (int)q ? n : (int)q) - 1;
                    u64 ws = (sq || q < (u64)n) ? (u64)(hi - lo + 1) : 0;
                    u64 r = rng(&seed) % (wa + ws), m, c;
                    if (r < wa) { u64 k = n + rng(&seed) % (q - n); m = q; c = (q - (k * (d % q)) % q) % q; }
                    else { u64 j = lo + rng(&seed) % ws; m = q * q; c = (m - (j * (d % m)) % m) % m; hitmask |= 1ULL << j; }
                    u64 am = (u64)(a % m), Mm = (u64)(M % m);
                    u64 x = mulmod((c + m - am) % m, inv_mod(Mm, m), m);
                    a += M * x; M *= m;
                }
                u64 a0;
                if (M >= Xa) { if (a >= Xa) continue; a0 = (u64)a; }
                else { u64 L = (u64)(Xa / M); a0 = (u64)(a + M * (rng(&seed) % L)); }
                bool ok = true; for (int i = 0; i < ndp; i++) if (a0 % dps[i] == 0) { ok = false; break; }
                if (!ok) continue;
                l_c++;
                u64 allm = (n == 64) ? ~0ULL : ((1ULL << n) - 1), dead = 0;
                for (int i = 0; i < nsv; i++) {
                    u64 q = sv[i], r = a0 % q;
                    if (!dinv[i]) { if (r == 0) { dead = allm; break; } continue; }
                    u64 j = ((q - r) % q) * dinv[i] % q;
                    for (; j < (u64)n; j += q) {
                        u64 t = a0 + j * d; int e = 0;
                        while (t % q == 0) { t /= q; e++; }
                        if (e & 1) dead |= 1ULL << j;
                    }
                    if (!full && (i & 7) == 7 && maxrun(allm & ~dead) < LMIN) break;
                }
                u64 alive = allm & ~dead;
                if (full) {
                    l_full++; l_alive += __builtin_popcountll(alive);
                    int np = 0;
                    for (int j = 0; j < n; j++) {
                        bool p = false;
                        if (alive >> j & 1) { p = loe(a0 + (u64)j * d); l_tests++; }
                        if (p) { l_pass[j]++; np++; }
                        if (hitmask >> j & 1) { l_hn++; l_hp += p; } else { l_nn++; l_np += p; }
                    }
                    (void)np;
                    continue;
                }
                if (maxrun(alive) < LMIN) { l_hist[0]++; continue; }
                u64 tested = 0;
                for (;;) {
                    /* find longest alive run containing an untested index */
                    int best = -1, bestlen = 0;
                    for (int s = 0; s < n;) {
                        if (!(alive >> s & 1)) { s++; continue; }
                        int e = s; while (e < n && (alive >> e & 1)) e++;
                        int len = e - s;
                        if (len >= LMIN && len > bestlen) {
                            /* untested index nearest the middle */
                            int mid = (s + e - 1) / 2, pick = -1;
                            for (int off = 0; off < len; off++) {
                                int c1 = mid - off, c2 = mid + off;
                                if (c1 >= s && !(tested >> c1 & 1)) { pick = c1; break; }
                                if (c2 < e && !(tested >> c2 & 1)) { pick = c2; break; }
                            }
                            if (pick >= 0) { best = pick; bestlen = len; }
                        }
                        s = e;
                    }
                    if (best < 0) break;
                    tested |= 1ULL << best; l_tests++;
                    if (!loe(a0 + (u64)best * d)) alive &= ~(1ULL << best);
                }
                int mr = maxrun(alive);
                if (mr < LMIN) mr = 0;
                l_hist[mr]++;
                if (mr >= 40) {
                    int s = 0; u64 x = alive; for (int i = 1; i < mr; i++) x &= x >> 1; s = __builtin_ctzll(x);
                    #pragma omp critical
                    { printf("LONGRUN %d a=%" PRIu64 " d=%" PRIu64 " start=%d first=%" PRIu64 "\n", mr, a0, d, s, a0 + (u64)s * d); fflush(stdout); }
                }
            }
        }
        #pragma omp critical
        {
            cpu += tcpu() - c0; tot_c += l_c; tot_tests += l_tests; for (int i = 0; i < 64; i++) { hist[i] += l_hist[i]; pass_idx[i] += l_pass[i]; }
            cnt_full += l_full; hit_n += l_hn; hit_pass += l_hp; non_n += l_nn; non_pass += l_np; alive_sum += l_alive;
            sumlog += l_sumlog; nlog += l_nlog;
        }
    }
    (void)t0;
    if (full) {
        double sl = 0; for (int j = 0; j < n; j++) sl += log((double)pass_idx[j] / cnt_full);
        printf("FULL fam=%s n=%d cands=%" PRIu64 " cpu_s=%.0f sieve_alive=%.4f p_nonhit=%.5f (n=%" PRIu64 ") p_hit=%.5f (n=%" PRIu64 ") sum_ln_p=%.4f geo_p=%.5f rho_rel=%.4f Xa=%.3e\n",
               fam, n, cnt_full, cpu, (double)alive_sum / cnt_full / n, (double)non_pass / non_n, non_n,
               hit_n ? (double)hit_pass / hit_n : 0.0, hit_n, sl, exp(sl / n), rrel, Xa_d);
    } else {
        printf("RUNS fam=%s n=%d LMIN=%d cands=%" PRIu64 " cpu_s=%.0f cand_per_cpu_s=%.0f tests=%" PRIu64 " rho_rel=%.4f Xa=%.3e\n",
               fam, n, LMIN, tot_c, cpu, tot_c / cpu, tot_tests, rrel, Xa_d);
        u64 cum = 0; printf("  cum(run>=L):");
        u64 cumv[64];
        for (int L = 63; L >= LMIN; L--) { cum += hist[L]; cumv[L] = cum; }
        for (int L = LMIN; L < 46; L++) if (cumv[L]) printf(" %d:%" PRIu64, L, cumv[L]);
        printf("\n");
    }
    return 0;
}

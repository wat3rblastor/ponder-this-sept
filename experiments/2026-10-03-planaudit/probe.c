#define _POSIX_C_SOURCE 200809L
/* probe -- per-unit yield measurement for the plan-ranking audit.
 *
 * Reproduces the stage-1 (tier B additive odometer) and stage-2 (tier C
 * bitmask) enumeration of src/c/apsearch.c EXACTLY (code copied, not
 * modified), but replaces stage 3 by a statistics pass: for every surviving
 * 58-term window it tests all 58 terms and accumulates
 *
 *   nwin_ok   windows whose first term is not killed by a bad prime q | d
 *             (such windows have ALL terms divisible by q, run length 0; they
 *             are a separate factor prod_{q|d}(1-1/q), not a per-term rate)
 *   nterm[k]  windows counted at position k, nloe[k] of which were Loeschian
 *   hist[L]   best run length per window (diagnostic)
 *
 * Measured log yield per stage-1 residue:
 *     M = log(nwin_ok/res) + sum_k log(nloe[k]/nterm[k])
 * which is exactly the quantity tools/plan_units.py's score v models
 * (up to a unit-independent constant).
 *
 * Build: cc -O3 -march=native -std=gnu11 -o build/probe experiments/2026-10-03-planaudit/probe.c -lm
 * Usage: probe --kmin K [--shift S] [--b2 200] [--modcap M] [--maxres N] [--secs T]
 */
#include <signal.h>
#include <math.h>
#include <time.h>
#include "../../src/c/loesch_core.h"

int main(int argc, char **argv) {
    u64 K = 0, shift = 0, nterms = 58, b2 = 200, D0 = D0_DEFAULT;
    u64 modcap = 20000000000000000ULL, maxres = 0;
    double secs = 0;
    for (int i = 1; i < argc; i++) {
        const char *k = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : "")
        if (!strcmp(k, "--kmin") || !strcmp(k, "--K")) K = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--shift")) shift = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--nterms")) nterms = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--b2")) b2 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--modcap")) modcap = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--maxres")) maxres = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--secs")) secs = atof(NEXT());
        else { fprintf(stderr, "bad arg %s\n", k); return 2; }
        #undef NEXT
    }
    if (!K) { fprintf(stderr, "need --kmin\n"); return 2; }
    signal(SIGINT, on_signal); signal(SIGTERM, on_signal);
    build_small_primes(10000);

    u64 d = K * D0;
    if (d / D0 != K) { fprintf(stderr, "overflow\n"); return 2; }

    /* ---- stage 1 components (verbatim from apsearch.c) ---- */
    u64 cm[MAXTIERB + 4], cs[MAXTIERB + 4], ct[MAXTIERB + 4], cc[MAXTIERB + 4];
    int nc = 0;
    u64 MOD = 1;
    #define ADDC(m_, s_, t_, c_) do { cm[nc] = (m_); cs[nc] = (s_); \
        ct[nc] = (t_); cc[nc] = (c_); MOD *= (m_); nc++; } while (0)
    ADDC(3, 1, 0, 1);
    if (D0 % 2 == 0) ADDC(2, 1, 0, 1);
    if (D0 % 5 == 0) ADDC(5, 1, 1, 4);
    int ntb = 0;
    for (u64 i = 0; i < n_sp && nc < MAXTIERB + 3; i++) {
        u64 q = sp[i];
        if (q <= nterms || q % 3 != 2 || d % q == 0) continue;
        if (MOD > (u64)4e18 / q) break;
        if (MOD * q > modcap) break;
        u64 dq = d % q;
        ntb++;
        ADDC(q, dq, dq, q - nterms);
    }
    #undef ADDC
    if (ntb < 4) { fprintf(stderr, "too few tier-B\n"); return 3; }
    u64 total_res = 1;
    for (int i = 0; i < nc; i++) total_res *= cc[i];

    u64 R0 = 0, sstep[MAXTIERB + 4], subcyc[MAXTIERB + 4];
    for (int i = 0; i < nc; i++) {
        u64 mi = cm[i], co = MOD / mi;
        u64 e = mulmod(co % MOD, inv_mod(co % mi, mi), MOD);
        R0 = (R0 + mulmod(cs[i], e, MOD)) % MOD;
        sstep[i] = mulmod(ct[i], e, MOD);
    }
    for (int i = 0; i < nc; i++) subcyc[i] = mulmod(cc[i], sstep[i], MOD);

    /* ---- tier C ---- */
    u64 tcp[MAXTIERC]; bar tcb[MAXTIERC]; u64 *tcw[MAXTIERC]; int ntc = 0;
    for (u64 i = 0; i < n_sp && ntc < MAXTIERC; i++) {
        u64 r = sp[i];
        if (r > b2) break;
        bool inB = false;
        for (int t = 0; t < nc; t++) if (cm[t] == r) inB = true;
        if (inB) continue;
        bool divD = (D0 % r == 0);
        if (!divD && (r % 3 != 2 || d % r == 0 || r <= nterms)) continue;
        tcp[ntc] = r; tcw[ntc] = malloc(r * sizeof(u64)); ntc++;
    }
    for (int i = 0; i < ntc; i++) {
        int bestj = i; double bk = -1;
        for (int j = i; j < ntc; j++) {
            double kill = (D0 % tcp[j] == 0) ? 1.0 / (double)tcp[j]
                                             : (double)nterms / (double)tcp[j];
            if (kill > bk) { bk = kill; bestj = j; }
        }
        u64 t1 = tcp[i]; tcp[i] = tcp[bestj]; tcp[bestj] = t1;
        u64 *t2 = tcw[i]; tcw[i] = tcw[bestj]; tcw[bestj] = t2;
    }
    for (int i = 0; i < ntc; i++) tcb[i] = bar_make(tcp[i]);

    u64 boff = shift * 64;
    for (int t = 0; t < ntc; t++) {
        u64 r = tcp[t];
        char *ok = calloc(r, 1);
        if (D0 % r == 0) { for (u64 x = 0; x < r; x++) ok[x] = (x != 0); }
        else {
            for (u64 x = 0; x < r; x++) ok[x] = 1;
            u64 v = 0, dm = d % r;
            for (u64 kk = 0; kk < nterms; kk++) { ok[v] = 0; v = (v + r - dm) % r; }
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

    /* bad primes q > nterms dividing d: not filtered anywhere, q | a kills the
     * whole window. Collect them so those windows can be accounted separately. */
    u64 dvq[64]; int ndv = 0;
    for (u64 i = 0; i < n_sp && ndv < 64; i++) {
        u64 q = sp[i];
        if (q <= nterms || q % 3 != 2) continue;
        if (d % q == 0) dvq[ndv++] = q;
    }

    /* ---- enumerate ---- */
    u64 idx[MAXTIERB + 4]; for (int i = 0; i < nc; i++) idx[i] = 0;
    u64 R = R0;
    u64 nres = 0, nsurv = 0, nwin_ok = 0, nkill = 0;
    u64 nterm_k[64] = {0}, nloe_k[64] = {0}, hist[64] = {0};
    double t0 = now_s();
    for (;;) {
        if (stop_requested) break;
        if (maxres && nres >= maxres) break;
        if (secs > 0 && (nres & 0x3F) == 0 && now_s() - t0 > secs) break;
        nres++;
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
            bool killed = false;
            for (int j = 0; j < ndv; j++) if (a0 % dvq[j] == 0) killed = true;
            if (killed) { nkill++; continue; }
            nwin_ok++;
            u64 cur = 0, best = 0;
            for (u64 k = 0; k < nterms; k++) {
                u64 t2 = a0 + k * d;
                if (t2 < a0) break;
                nterm_k[k]++;
                if (is_loeschian(t2)) { nloe_k[k]++; cur++; if (cur > best) best = cur; }
                else cur = 0;
            }
            hist[best]++;
        }
        int i = nc - 1; bool exhausted = false;
        for (;;) {
            while (i >= 0 && cc[i] <= 1) i--;
            if (i < 0) { exhausted = true; break; }
            R += sstep[i]; if (R >= MOD) R -= MOD;
            idx[i]++;
            if (idx[i] < cc[i]) break;
            R = (R + MOD - subcyc[i]) % MOD;
            idx[i] = 0; i--;
        }
        if (exhausted) break;
    }

    double sumlog = 0; u64 tt = 0, tl = 0;
    for (u64 k = 0; k < nterms; k++) {
        if (!nterm_k[k]) { sumlog = NAN; break; }
        sumlog += log((double)(nloe_k[k] + 0.5) / (double)(nterm_k[k] + 1.0));
        tt += nterm_k[k]; tl += nloe_k[k];
    }
    printf("{\"K\":%llu,\"shift\":%llu,\"MOD\":%llu,\"total_res\":%.6g,"
           "\"res\":%llu,\"surv\":%llu,\"win_ok\":%llu,\"kill\":%llu,"
           "\"nterm\":%llu,\"nloe\":%llu,\"sumlogp\":%.6f,\"secs\":%.2f,"
           "\"hist\":[", (unsigned long long)K, (unsigned long long)shift,
           (unsigned long long)MOD, (double)total_res,
           (unsigned long long)nres, (unsigned long long)nsurv,
           (unsigned long long)nwin_ok, (unsigned long long)nkill,
           (unsigned long long)tt, (unsigned long long)tl, sumlog, now_s() - t0);
    for (int L = 0; L < 40; L++) printf("%s%llu", L ? "," : "", (unsigned long long)hist[L]);
    printf("]}\n");
    for (int t = 0; t < ntc; t++) free(tcw[t]);
    return 0;
}

#define _POSIX_C_SOURCE 200809L
/* apsearch_sq -- apsearch.c + the SQUARE OPTION.
 *
 * Build: cc -O3 -march=native -std=gnu11 -I src/c -o build/apsearch_sq \
 *            experiments/2026-10-03-sqmode/apsearch_sq.c -lm
 *
 * Difference from src/c/apsearch.c (which this file must not change):
 *
 *   A bad prime q with n/2 < q < n is NOT forced to divide d. Either q | d
 *   ("in-d mode", what the production engine always does), or q divides
 *   exactly one term t_{i0} of the window and q^2 | t_{i0} ("square mode").
 *   In square mode the admissible a form an arithmetic progression of
 *   residues modulo q^2:
 *        a = -i0*d (mod q^2),   i0 = max(0,n-q) .. min(n-1,q-1)
 *   which is exactly the (modulus, start, step, count) shape stage 1 already
 *   enumerates by additions, with
 *        modulus = q^2, start = -lo*d, step = -d, count = hi-lo+1.
 *   For q < n that count is 2q-n.
 *
 * The unit of work is now (K, shift, mode).  The mode is
 *     d    = K * D0 * prod(--ind primes)        [in-d primes]
 *     plus a q^2 component for every --square prime
 * and it is written into every output record ("ind" and "sq" fields).
 *
 * Primes in square mode are excluded from tier B and from tier C: tier C's
 * filter is "q misses the whole window", which is the exact negation of what
 * square mode generates, so applying it would reject every candidate we want.
 *
 * Extra instrumentation for measurement: a run-length histogram per unit
 * (cum counts at >= 30/35/40/45) and the smallest last term seen among runs
 * of length >= --histmin.
 */

#include <signal.h>
#include <time.h>

#include "loesch_core.h"

#define MAXC 40
#define MAXMODE 8

/* ------------------------------------------------------------------ main --- */

static void usage(const char *p) {
    fprintf(stderr,
      "usage: %s --kmin K --kmax K [--nterms 58] [--shifts S] [--b2 2000]\n"
      "          [--report 50] [--out f.jsonl] [--D0 N] [--selftest]\n"
      "          [--ind q,q,..] [--square q,q,..] [--mintierb 4] [--histmin 30]\n"
      "\n"
      "d = K*D0*prod(ind). Each (K, shift, mode) is an independent work unit.\n"
      "--square lists bad primes q put in SQUARE mode (q must not divide d;\n"
      "q^2 then divides exactly one term of the window).\n", p);
    exit(2);
}

static int parse_list(const char *s, u64 *out, int cap) {
    int n = 0;
    while (*s && n < cap) {
        char *e;
        u64 v = strtoull(s, &e, 10);
        if (e == s) break;
        out[n++] = v;
        s = e;
        while (*s == ',' || *s == ' ') s++;
    }
    return n;
}

static void list_str(char *buf, size_t bl, const u64 *v, int n) {
    size_t o = 0;
    buf[0] = 0;
    for (int i = 0; i < n && o + 24 < bl; i++)
        o += (size_t)snprintf(buf + o, bl - o, "%s%" PRIu64, i ? "," : "", v[i]);
}

int main(int argc, char **argv) {
    u64 kmin = 1, kmax = 0, nterms = 58, shifts = 16, b2 = 2000, D0 = D0_DEFAULT;
    u64 modcap = 6000000000000000ULL;   /* cap on MOD = 3*prod(tier B) */
    int report = 50, mintierb = 4, histmin = 30;
    const char *out = NULL;
    bool selftest = false, isl_mode = false;
    u64 ind[MAXMODE], sq[MAXMODE];
    int nind = 0, nsq = 0;
    double tlimit = 0;
    u64 atest = 0;            /* --atest A: is this exact a admissible? */

    for (int i = 1; i < argc; i++) {
        const char *k = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : (usage(argv[0]), ""))
        if (!strcmp(k, "--kmin")) kmin = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--kmax")) kmax = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--nterms")) nterms = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--shifts")) shifts = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--b2")) b2 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--report")) report = atoi(NEXT());
        else if (!strcmp(k, "--out")) out = NEXT();
        else if (!strcmp(k, "--D0")) D0 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--modcap")) modcap = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--mintierb")) mintierb = atoi(NEXT());
        else if (!strcmp(k, "--histmin")) histmin = atoi(NEXT());
        else if (!strcmp(k, "--tlimit")) tlimit = atof(NEXT());
        else if (!strcmp(k, "--atest")) atest = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--ind")) nind = parse_list(NEXT(), ind, MAXMODE);
        else if (!strcmp(k, "--square")) nsq = parse_list(NEXT(), sq, MAXMODE);
        else if (!strcmp(k, "--selftest")) selftest = true;
        else if (!strcmp(k, "--isl")) isl_mode = true;
        else usage(argv[0]);
        #undef NEXT
    }
    if (!kmax) kmax = kmin;
    if (D0 % 3) { fprintf(stderr, "FATAL: 3 must divide D0\n"); return 2; }

    signal(SIGINT, on_signal); signal(SIGTERM, on_signal);
    build_small_primes(10000);

    if (isl_mode) {          /* read decimal values on stdin, print verdicts */
        u64 v;
        while (scanf("%" SCNu64, &v) == 1) printf("%" PRIu64 " %d\n", v, (int)is_loeschian(v));
        return 0;
    }

    if (selftest) {
        static const int pre[] = {0,1,3,4,7,9,12,13,16,19,21,25,27,28,31,36,37,
                                  39,43,48,49,52,57,61,63,64,67,73,75,76,79,81,
                                  84,91,93,97,100};
        int j = 0, fail = 0;
        for (int n = 0; n <= 100; n++) {
            bool want = (j < 37 && pre[j] == n); if (want) j++;
            if (is_loeschian((u64)n) != want) {
                fprintf(stderr, "SELFTEST FAIL at n=%d\n", n); fail = 1; }
        }
        u64 bad[] = {2,5,11,17,23,29,41,47,53,59,71,83,89,101,1000037,1000121};
        for (unsigned i = 0; i < sizeof(bad)/sizeof(bad[0]); i++) {
            u64 p = bad[i];
            if (p % 3 != 2) continue;
            if (is_loeschian(p)) { fprintf(stderr, "SELFTEST FAIL: %" PRIu64
                                   " marked Loeschian\n", p); fail = 1; }
            if (!is_loeschian(p*p)) { fprintf(stderr, "SELFTEST FAIL: %" PRIu64
                                   "^2 not Loeschian\n", p); fail = 1; }
            if (is_loeschian(p*p*p)) { fprintf(stderr, "SELFTEST FAIL: %" PRIu64
                                   "^3 marked Loeschian\n", p); fail = 1; }
        }
        for (int k = 0; k < 9; k++)
            if (!is_loeschian(2432167ULL + (u64)k * 1405152ULL)) {
                fprintf(stderr, "SELFTEST FAIL: Choudhry 9-term AP term %d\n", k);
                fail = 1; }
        if (is_loeschian(1000037ULL * 1000121ULL)) {
            fprintf(stderr, "SELFTEST FAIL: bad*bad semiprime\n"); fail = 1; }
        if (!is_loeschian(1000037ULL * 1000037ULL)) {
            fprintf(stderr, "SELFTEST FAIL: bad^2 large\n"); fail = 1; }
        if (is_loeschian(1000037ULL * 1000037ULL * 1000037ULL)) {
            fprintf(stderr, "SELFTEST FAIL: bad^3 large\n"); fail = 1; }
        if (!is_loeschian(1000003ULL * 1000037ULL * 1000037ULL)) {
            fprintf(stderr, "SELFTEST FAIL: good*bad^2\n"); fail = 1; }
        for (u64 t = 0; t <= 20000; t++) {
            bool bf = false;
            for (u64 x = 0; x * x <= t && !bf; x++)
                for (u64 y = 0; x * x + x * y + y * y <= t; y++)
                    if (x * x + x * y + y * y == t) { bf = true; break; }
            if (bf != is_loeschian(t)) {
                fprintf(stderr, "SELFTEST FAIL: brute force disagrees at %"
                        PRIu64 "\n", t); fail = 1; break; }
        }
        fprintf(stderr, "selftest: %s\n", fail ? "FAILED" : "OK");
        return fail;
    }

    /* ---- mode validation ---- */
    u64 Dfull = D0;
    for (int i = 0; i < nind; i++) {
        u64 q = ind[i];
        if (!is_prime_u64(q) || q % 3 != 2) {
            fprintf(stderr, "FATAL --ind %" PRIu64 ": not a bad prime\n", q);
            return 2; }
        if (Dfull > (u64)1e18 / q) { fprintf(stderr, "FATAL: D overflow\n"); return 2; }
        Dfull *= q;
    }
    for (int i = 0; i < nsq; i++) {
        u64 q = sq[i];
        if (!is_prime_u64(q) || q % 3 != 2) {
            fprintf(stderr, "FATAL --square %" PRIu64 ": not a bad prime\n", q);
            return 2; }
        if (Dfull % q == 0) {
            fprintf(stderr, "FATAL --square %" PRIu64 ": divides D (mode "
                    "contradiction)\n", q); return 2; }
        if (2 * q <= nterms) {
            fprintf(stderr, "FATAL --square %" PRIu64 ": 2q <= n, q MUST divide "
                    "d (two forced hits)\n", q); return 2; }
    }
    char indstr[128], sqstr[128];
    list_str(indstr, sizeof indstr, ind, nind);
    list_str(sqstr, sizeof sqstr, sq, nsq);

    FILE *of = out ? fopen(out, "a") : NULL;
    if (out && !of) { perror("open --out"); return 2; }

    int global_best = 0;
    u64 gb_a = 0, gb_d = 0;
    double t0 = now_s();
    u64 units = 0;
    double covered = 0;       /* raw a-range covered */
    u64 ghist[256]; for (int i = 0; i < 256; i++) ghist[i] = 0;
    u64 gmin_last = 0;        /* smallest last term among runs >= histmin */
    u64 gres = 0, gsurv = 0, gconf = 0;

    for (u64 K = kmin; K <= kmax && !stop_requested; K++) {
        u64 d = K * Dfull;
        if (d / Dfull != K) { fprintf(stderr, "K overflow at %" PRIu64 "\n", K); break; }
        bool modebad = false;
        for (int i = 0; i < nsq; i++) if (d % sq[i] == 0) {
            fprintf(stderr, "K=%" PRIu64 ": square prime %" PRIu64 " divides d "
                    "(K multiple), skipping\n", K, sq[i]); modebad = true; }
        if (modebad) continue;

        /* ---- stage 1 components. Each is (modulus, start, step, count) ---- */
        u64 cm[MAXC], cs[MAXC], ct[MAXC], cc[MAXC];
        int nc = 0;
        u64 MOD = 1;
        #define ADDC(m_, s_, t_, c_) do { if (nc >= MAXC) { \
            fprintf(stderr, "FATAL: too many components\n"); return 2; } \
            cm[nc] = (m_); cs[nc] = (s_); \
            ct[nc] = (t_); cc[nc] = (c_); MOD *= (m_); nc++; } while (0)
        ADDC(3, 1, 0, 1);
        if (d % 2 == 0) ADDC(2, 1, 0, 1);
        if (d % 5 == 0) ADDC(5, 1, 1, 4);

        /* ---- SQUARE components: a = -i0*d (mod q^2), i0 = lo..hi ---- */
        for (int i = 0; i < nsq; i++) {
            u64 q = sq[i], m = q * q;
            u64 lo = (nterms > q) ? nterms - q : 0;
            u64 hi = (nterms - 1 < q - 1) ? nterms - 1 : q - 1;
            if (hi < lo) { fprintf(stderr, "FATAL square %" PRIu64 ": empty "
                                   "index range\n", q); return 2; }
            u64 dm = d % m;
            u64 start = (m - mulmod(lo % m, dm, m)) % m;
            u64 step = (m - dm) % m;
            ADDC(m, start, step, hi - lo + 1);
        }

        int ntb = 0;
        u64 tbq[MAXTIERB];
        for (u64 i = 0; i < n_sp && nc < MAXC - 1; i++) {
            u64 q = sp[i];
            if (q <= nterms || q % 3 != 2 || d % q == 0) continue;
            bool issq = false;
            for (int t = 0; t < nsq; t++) if (sq[t] == q) issq = true;
            if (issq) continue;                      /* handled mod q^2 above */
            if (MOD > (u64)4e18 / q) break;          /* keep MOD in u64 */
            if (MOD * q > modcap) break;             /* --modcap: unit size */
            if (ntb >= MAXTIERB) break;
            u64 dq = d % q;
            tbq[ntb++] = q;
            ADDC(q, dq, dq, q - nterms);             /* a = j*d, j = 1..q-n */
        }
        #undef ADDC
        if (ntb < mintierb) { fprintf(stderr, "K=%" PRIu64 ": too few tier-B "
                              "primes (%d), skipping\n", K, ntb); continue; }

        /* ---- CRT: idempotents, base residue R0, additive steps ---- */
        u64 R0 = 0, sstep[MAXC], subcyc[MAXC];
        for (int i = 0; i < nc; i++) {
            u64 mi = cm[i], co = MOD / mi;
            u64 e = mulmod(co % MOD, inv_mod(co % mi, mi), MOD); /* idempotent */
            R0 = (R0 + mulmod(cs[i], e, MOD)) % MOD;
            sstep[i] = mulmod(ct[i], e, MOD);
        }
        for (int i = 0; i < nc; i++) subcyc[i] = mulmod(cc[i], sstep[i], MOD);

        /* CRT self-check: R0 must hit the first good residue of every
         * component (including the new q^2 ones) and sstep[i] must move ONLY
         * component i. */
        for (int i = 0; i < nc; i++) {
            if (R0 % cm[i] != cs[i] % cm[i]) {
                fprintf(stderr, "FATAL CRT: R0=%" PRIu64 " mod %" PRIu64 " = %"
                        PRIu64 ", want %" PRIu64 "\n",
                        R0, cm[i], R0 % cm[i], cs[i] % cm[i]); return 2; }
            for (int j = 0; j < nc; j++) {
                u64 got = sstep[i] % cm[j], want2 = (i == j) ? ct[i] % cm[j] : 0;
                if (got != want2) {
                    fprintf(stderr, "FATAL CRT: sstep[%d] mod %" PRIu64 " = %"
                            PRIu64 ", want %" PRIu64 "\n", i, cm[j], got, want2);
                    return 2; }
            }
        }
        /* square-mode semantic self-check: every residue of the q^2 component
         * must make q^2 divide exactly one in-window term, at an index that
         * leaves no second hit. */
        for (int i = 0; i < nsq; i++) {
            u64 q = sq[i], m = q * q;
            u64 lo = (nterms > q) ? nterms - q : 0;
            u64 hi = (nterms - 1 < q - 1) ? nterms - 1 : q - 1;
            for (u64 i0 = lo; i0 <= hi; i0++) {
                u64 a = (m - mulmod(i0 % m, d % m, m)) % m;   /* a = -i0*d */
                int hits = 0, sqhits = 0;
                for (u64 k = 0; k < nterms; k++) {
                    u64 t2 = (a + mulmod(k % m, d % m, m)) % m;
                    if (t2 % q == 0) hits++;
                    if (t2 == 0) sqhits++;
                }
                if (hits != 1 || sqhits != 1) {
                    fprintf(stderr, "FATAL square self-check q=%" PRIu64
                            " i0=%" PRIu64 ": hits=%d sqhits=%d\n",
                            q, i0, hits, sqhits); return 2; }
            }
        }

        /* ---- tier C: bitmask primes ---- */
        u64 tcp[MAXTIERC]; bar tcb[MAXTIERC]; u64 *tcw[MAXTIERC]; int ntc = 0;
        for (u64 i = 0; i < n_sp && ntc < MAXTIERC; i++) {
            u64 r = sp[i];
            if (r > b2) break;
            bool inB = false;
            for (int t = 0; t < nc; t++) if (cm[t] == r) inB = true;
            /* CRITICAL: a square-mode prime must NEVER get tier C's
             * "q misses the window" mask -- that mask is the complement of the
             * square family and would reject every candidate we generate. */
            for (int t = 0; t < nsq; t++) if (sq[t] == r) inB = true;
            if (inB) continue;
            bool divD = (Dfull % r == 0);
            if (!divD && (r % 3 != 2 || d % r == 0 || r <= nterms)) continue;
            tcp[ntc] = r;
            tcw[ntc] = malloc(r * sizeof(u64));
            ntc++;
        }

        /* ---- --atest: would the engine's filters ACCEPT this exact a? ----
         * Checks stage 1 (a mod m is one of start+j*step), stage 2 / tier C
         * (a mod r is a good residue), then stage 3 (the run length). This is
         * the cheap form of the record regression: the full enumeration of a
         * unit this size is ~1e11 residues. */
        if (atest) {
            int bad = 0;
            for (int i = 0; i < nc; i++) {
                u64 r = atest % cm[i];
                int ok = 0;
                u64 v = cs[i] % cm[i];
                for (u64 j = 0; j < cc[i]; j++) {
                    if (v == r) { ok = 1; break; }
                    v = (v + ct[i]) % cm[i];
                }
                printf("  stage1 m=%-6" PRIu64 " a mod m=%-6" PRIu64 " %s\n",
                       cm[i], r, ok ? "OK" : "REJECT");
                if (!ok) bad = 1;
            }
            for (int t = 0; t < ntc; t++) {
                u64 r = tcp[t], ar = atest % r;
                int ok;
                if (Dfull % r == 0) ok = (ar != 0);
                else { ok = 1;
                    u64 v = 0, dm = d % r;
                    for (u64 k = 0; k < nterms; k++) {
                        if (ar == v) ok = 0;
                        v = (v + r - dm) % r; } }
                if (!ok) { printf("  stage2 r=%" PRIu64 " a mod r=%" PRIu64
                                  " REJECT\n", r, ar); bad = 1; }
            }
            u64 run = 0;
            while (run < 200 && is_loeschian(atest + run * d)) run++;
            printf("ATEST a=%" PRIu64 " d=%" PRIu64 " ind=[%s] sq=[%s]: "
                   "stage1+2 %s, loeschian run from a = %" PRIu64 "\n",
                   atest, d, indstr, sqstr, bad ? "REJECTED" : "ADMISSIBLE", run);
            for (int t = 0; t < ntc; t++) free(tcw[t]);
            continue;
        }

        for (int i = 0; i < ntc; i++) {
            int bestj = i;
            double bk = -1;
            for (int j = i; j < ntc; j++) {
                double kill = (Dfull % tcp[j] == 0) ? 1.0 / (double)tcp[j]
                                                 : (double)nterms / (double)tcp[j];
                if (kill > bk) { bk = kill; bestj = j; }
            }
            u64 t1 = tcp[i]; tcp[i] = tcp[bestj]; tcp[bestj] = t1;
            u64 *t2 = tcw[i]; tcw[i] = tcw[bestj]; tcw[bestj] = t2;
        }
        for (int i = 0; i < ntc; i++) tcb[i] = bar_make(tcp[i]);
        if (K == kmin)
            fprintf(stderr, "tierC order: %" PRIu64 " %" PRIu64 " %" PRIu64
                    " %" PRIu64 " %" PRIu64 " ... (%d primes)\n",
                    tcp[0], tcp[1], tcp[2], tcp[3], tcp[4], ntc);

        for (u64 shift = 0; shift < shifts && !stop_requested; shift++) {
            u64 boff = shift * 64;

            for (int t = 0; t < ntc; t++) {
                u64 r = tcp[t];
                char *ok = calloc(r, 1);
                if (Dfull % r == 0) {
                    for (u64 x = 0; x < r; x++) ok[x] = (x != 0);
                } else {
                    for (u64 x = 0; x < r; x++) ok[x] = 1;
                    u64 v = 0, dm = d % r;
                    for (u64 k = 0; k < nterms; k++) {
                        ok[v] = 0;
                        v = (v + r - dm) % r;
                    }
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

            /* ---- stage 1: nested additive enumeration of admissible R ---- */
            u64 idx[MAXC];
            for (int i = 0; i < nc; i++) idx[i] = 0;
            u64 R = R0;

            u64 nres = 0, nsurv = 0, nconf = 0;
            u64 hist[256]; for (int i = 0; i < 256; i++) hist[i] = 0;
            u64 min_last = 0;
            double ts = now_s();
            {   double tot = 1;
                for (int i = 0; i < nc; i++) tot *= (double)cc[i];
                fprintf(stderr, "unit K=%" PRIu64 " shift=%" PRIu64 " ind=[%s] "
                        "sq=[%s] d=%" PRIu64 " MOD=%" PRIu64 " stage1=%d(",
                        K, shift, indstr, sqstr, d, MOD, nc);
                for (int i = 0; i < nc; i++)
                    fprintf(stderr, "%s%" PRIu64 "x%" PRIu64, i ? "," : "", cm[i], cc[i]);
                fprintf(stderr, ") tierC=%d residues=%.4g range=%.4g\n",
                        ntc, tot, (double)MOD * 64.0);
            }

            for (;;) {
                if (stop_requested) break;
                nres++;
                if ((nres & 0x3FFFFFF) == 0) {
                    double el2 = now_s() - ts;
                    fprintf(stderr, "  .. K=%" PRIu64 " shift=%" PRIu64 " R=%.4g"
                            " surv=%" PRIu64 " conf=%" PRIu64 " best=%d"
                            " %.0fs (%.3g res/s)\n", K, shift, (double)nres,
                            nsurv, nconf, global_best, el2, nres / el2);
                    if (tlimit > 0 && now_s() - t0 > tlimit) stop_requested = 1;
                }

                /* ---- stage 2: 64 candidates at once ---- */
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
                    u64 best_run = 0, best_start = 0, cur = 0, cur_start = 0;
                    u64 want = (report > 1) ? (u64)report : 1;
                    for (u64 k = 0; k < nterms; k++) {
                        if (cur + (nterms - k) < want) break;
                        u64 t2 = a0 + k * d;
                        if (t2 < a0) break;              /* u64 overflow guard */
                        if (is_loeschian(t2)) {
                            if (cur == 0) cur_start = k;
                            cur++;
                            if (cur > best_run) { best_run = cur; best_start = cur_start; }
                        } else cur = 0;
                    }
                    if (best_run == 0) continue;
                    nconf++;
                    u64 a = a0 + best_start * d;
                    u64 run = best_run;
                    while (a >= d && is_loeschian(a - d)) { a -= d; run++; }
                    for (;;) {
                        u64 t2 = a + run * d;
                        if (t2 < a || !is_loeschian(t2)) break;
                        run++;
                    }
                    hist[run < 255 ? run : 255]++;
                    if ((int)run >= histmin) {
                        u64 last = a + (run - 1) * d;
                        if (!min_last || last < min_last) min_last = last;
                        if (!gmin_last || last < gmin_last) gmin_last = last;
                    }
                    if ((int)run >= report || (int)run > global_best) {
                        if ((int)run > global_best) {
                            global_best = (int)run; gb_a = a; gb_d = d;
                            fprintf(stderr, "*** n=%d a=%" PRIu64 " d=%" PRIu64
                                    " (K=%" PRIu64 " shift=%" PRIu64 " ind=[%s]"
                                    " sq=[%s])\n",
                                    global_best, a, d, K, shift, indstr, sqstr);
                        }
                        if (of && (int)run >= report) {
                            fprintf(of, "{\"hit\":true,\"n\":%" PRIu64 ",\"a\":%"
                                    PRIu64 ",\"d\":%" PRIu64 ",\"K\":%" PRIu64
                                    ",\"ind\":\"%s\",\"sq\":\"%s\"}\n",
                                    run, a, d, K, indstr, sqstr);
                            fflush(of);
                        }
                    }
                }

                /* Advance the odometer. */
                int i = nc - 1;
                bool exhausted = false;
                for (;;) {
                    while (i >= 0 && cc[i] <= 1) i--;   /* fixed components */
                    if (i < 0) { exhausted = true; break; }
                    R += sstep[i]; if (R >= MOD) R -= MOD;
                    idx[i]++;
                    if (idx[i] < cc[i]) break;
                    R = (R + MOD - subcyc[i]) % MOD;    /* undo cc[i] steps */
                    idx[i] = 0;
                    i--;
                }
                if (exhausted) break;
            }

            covered += (double)MOD * 64.0;
            units++;
            gres += nres; gsurv += nsurv; gconf += nconf;
            for (int i = 0; i < 256; i++) ghist[i] += hist[i];
            double el = now_s() - ts;
            /* cumulative counts */
            u64 c20 = 0, c25 = 0, c30 = 0, c35 = 0, c40 = 0;
            for (int i = 20; i < 256; i++) { c20 += hist[i];
                if (i >= 25) c25 += hist[i]; if (i >= 30) c30 += hist[i];
                if (i >= 35) c35 += hist[i]; if (i >= 40) c40 += hist[i]; }
            fprintf(stderr, "K=%" PRIu64 " shift=%" PRIu64 " MOD=%" PRIu64
                    " tierB=%d tierC=%d R=%" PRIu64 " surv=%" PRIu64
                    " conf=%" PRIu64 " ge20=%" PRIu64 " ge25=%" PRIu64
                    " ge30=%" PRIu64 " ge35=%" PRIu64 " ge40=%" PRIu64
                    " minlast=%" PRIu64 " best=%d %.1fs (%.3g raw/s)\n",
                    K, shift, MOD, ntb, ntc, nres, nsurv, nconf,
                    c20, c25, c30, c35, c40, min_last, global_best,
                    el, (double)MOD * 64.0 / el);
            if (of) {
                fprintf(of, "{\"K\":%" PRIu64 ",\"shift\":%" PRIu64 ",\"ind\":\"%s\""
                        ",\"sq\":\"%s\",\"MOD\":%" PRIu64 ",\"d\":%" PRIu64
                        ",\"res\":%" PRIu64 ",\"surv\":%" PRIu64
                        ",\"conf\":%" PRIu64 ",\"ge20\":%" PRIu64
                        ",\"ge25\":%" PRIu64 ",\"ge30\":%" PRIu64
                        ",\"ge35\":%" PRIu64 ",\"ge40\":%" PRIu64
                        ",\"minlast\":%" PRIu64 ",\"best\":%d,\"secs\":%.2f,"
                        "\"covered\":%.6g}\n",
                        K, shift, indstr, sqstr, MOD, d, nres, nsurv, nconf,
                        c20, c25, c30, c35, c40, min_last, global_best, el, covered);
                fflush(of);
            }
            if (tlimit > 0 && now_s() - t0 > tlimit) stop_requested = 1;
        }
        for (int t = 0; t < ntc; t++) free(tcw[t]);
    }

    if (of) fclose(of);
    double tot_s = now_s() - t0;
    u64 G20 = 0, G25 = 0, G30 = 0, G35 = 0, G40 = 0;
    for (int i = 20; i < 256; i++) { G20 += ghist[i];
        if (i >= 25) G25 += ghist[i]; if (i >= 30) G30 += ghist[i];
        if (i >= 35) G35 += ghist[i]; if (i >= 40) G40 += ghist[i]; }
    printf("BEST n=%d a=%" PRIu64 " d=%" PRIu64 "  units=%" PRIu64
           " covered=%.4g raw a  %.1fs%s\n", global_best, gb_a, gb_d, units,
           covered, tot_s, stop_requested ? " INTERRUPTED" : "");
    printf("TOTALS ind=[%s] sq=[%s] secs=%.2f res=%" PRIu64 " surv=%" PRIu64
           " conf=%" PRIu64 " ge20=%" PRIu64 " ge25=%" PRIu64 " ge30=%" PRIu64
           " ge35=%" PRIu64 " ge40=%" PRIu64 " minlast=%" PRIu64 "\n",
           indstr, sqstr, tot_s, gres, gsurv, gconf, G20, G25, G30, G35, G40,
           gmin_last);
    printf("HIST");
    for (int i = 1; i < 256; i++) if (ghist[i]) printf(" %d:%" PRIu64, i, ghist[i]);
    printf("\n");
    return 0;
}

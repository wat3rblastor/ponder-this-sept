/* apsearch_cuda -- the apsearch.c engine with stage 1 + stage 2 on an NVIDIA GPU.
 *
 * Build (see Makefile target `cuda`):
 *   gcc  -O3 -march=native -std=gnu11 -c -o build/loesch_core_api.o src/c/loesch_core_api.c
 *   nvcc -O3 -arch=sm_121 -Xcompiler -fopenmp -I src/c -o build/apsearch_cuda \
 *        src/c/apsearch_cuda.cu build/loesch_core_api.o -lm
 *
 * This is the CUDA twin of apsearch_gpu.m (the Metal engine). Same algorithm,
 * same work-unit definition (K, shift), same output format, so coverage logs
 * from the two are interchangeable:
 *
 *   d = K * D0.  Stage 1 pins a mod {3,2,5} and mod the smallest bad primes
 *   q > n (tier B) by direct CRT enumeration of admissible residues R < MOD.
 *   Stage 2 ANDs 64-bit OKOK words for every remaining bad prime r <= b2
 *   (tier C): one residue R stands for the 64 candidates a = R + (b+64*shift)*MOD.
 *   Stage 3 (CPU, OpenMP) runs the exact Loeschian test on survivors, which is
 *   the ONE copy in loesch_core.h shared with the CPU engine. src/verify.py
 *   still mints every record.
 *
 * GPU specifics:
 *   - one thread per prefix of the additive loop nest; each thread walks the
 *     two innermost (largest-count) components.
 *   - R < MOD < 2^64 is split into four 16-bit limbs and folded with
 *     precomputed 2^16, 2^32, 2^48 residues, so the per-prime index needs one
 *     32-bit modulo; that modulo is a magic-multiply (umulhi) plus at most
 *     two conditional subtractions, never a hardware divide. The fold is
 *     < 2^16 * (1 + 3*b2) < 2^32 for b2 <= 21845, asserted at startup, and
 *     the magic reduction is checked against % on the host for every prime.
 *   - per-prime constants live in __constant__ memory (warp-uniform reads).
 *   - two CUDA streams: the kernel for unit i runs while the CPU does stage 3
 *     for unit i-1 and builds the tables for unit i+1.
 *
 * Checkpointing: every finished unit is a line in --out; --resume re-reads
 * that file and skips the (K, shift) pairs already recorded there.
 */

#include <cuda_runtime.h>
#include <omp.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <set>
#include <string>
#include <algorithm>

#include "loesch_core_api.h"

typedef uint64_t u64;
typedef uint32_t u32;
typedef __uint128_t u128;

#define D0_DEFAULT 382160924970ULL   /* 3 * 2*5*11*17*23*29*41*47*53 */
#define MAXC 16
#define MAXTC 2560                   /* tier-C primes (bad primes <= b2) */

static volatile sig_atomic_t stop_requested = 0;
static void on_signal(int s) { (void)s; stop_requested = 1; }

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
    exit(2); } } while (0)

/* ------------------------------------------------------------- the kernel -- */

struct Params {
    u64 MOD;
    u64 s1, s2;      /* additive steps of the two innermost components */
    u32 c1, c2;      /* their trip counts */
    u32 ntc;         /* number of tier-C tables */
    u32 cap;         /* survivor buffer capacity */
    u32 nthreads;
};

/* per-prime constants, warp-uniform: offset into words, prime, magic, limb folds */
__constant__ u32 c_off[MAXTC];
__constant__ u32 c_rp[MAXTC];
__constant__ u32 c_mag[MAXTC];
__constant__ u32 c_k1[MAXTC];
__constant__ u32 c_k2[MAXTC];
__constant__ u32 c_k3[MAXTC];

__global__ void __launch_bounds__(256)
sieve(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
      Params P, u32 *cnt, u64 *hits)
{
    u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= P.nthreads) return;
    u64 R = Rpre[gid];
    for (u32 i1 = 0; i1 < P.c1; ++i1) {
        u64 Rb = R;
        for (u32 i2 = 0; i2 < P.c2; ++i2) {
            u32 l0 = (u32)(Rb & 0xFFFF);
            u32 l1 = (u32)((Rb >> 16) & 0xFFFF);
            u32 l2 = (u32)((Rb >> 32) & 0xFFFF);
            u32 l3 = (u32)(Rb >> 48);
            u64 sito = ~0ULL;
            for (u32 t = 0; t < P.ntc; ++t) {
                u32 x = l0 + l1 * c_k1[t] + l2 * c_k2[t] + l3 * c_k3[t];
                u32 p = c_rp[t];
                u32 q = __umulhi(x, c_mag[t]);
                u32 r = x - q * p;
                if (r >= p) r -= p;
                if (r >= p) r -= p;
                sito &= words[c_off[t] + r];
                if (sito == 0ULL) break;
            }
            if (sito != 0ULL) {
                u32 s = atomicAdd(cnt, 1u);
                if (s < P.cap) { hits[2 * s] = Rb; hits[2 * s + 1] = sito; }
            }
            Rb += P.s2; if (Rb >= P.MOD) Rb -= P.MOD;
        }
        R += P.s1; if (R >= P.MOD) R -= P.MOD;
    }
}

/* ------------------------------------------------------------------ host --- */

struct comp { u64 m, s, t, c; };

struct Unit {                 /* one (K, shift) work unit in flight */
    u64 K, shift, d, MOD, total, nthreads;
    int ntc;
    bool valid;
    double t_launch;
    cudaEvent_t done;
};

struct Hit { u64 n, a, d, K; };

static std::string unit_key(u64 K, u64 shift) {
    char b[64]; snprintf(b, sizeof b, "%llu/%llu", (unsigned long long)K,
                         (unsigned long long)shift); return b;
}

int main(int argc, char **argv) {
    u64 kmin = 1, kmax = 0, nterms = 58, shifts = 1, b2 = 10000, D0 = D0_DEFAULT;
    u64 shift0 = 0;
    u64 modcap = 2000000000000000ULL;
    int report = 44;
    const char *out = NULL;
    bool verify_mode = false, resume = false;
    int nomp = 0;

    for (int i = 1; i < argc; i++) {
        const char *k = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : "")
        if (!strcmp(k, "--kmin")) kmin = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--kmax")) kmax = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--nterms")) nterms = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--shifts")) shifts = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--shift0")) shift0 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--b2")) b2 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--modcap")) modcap = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--report")) report = atoi(NEXT());
        else if (!strcmp(k, "--out")) out = NEXT();
        else if (!strcmp(k, "--D0")) D0 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--threads")) nomp = atoi(NEXT());
        else if (!strcmp(k, "--verify")) verify_mode = true;
        else if (!strcmp(k, "--resume")) resume = true;
        else { fprintf(stderr, "unknown arg %s\n", k); return 2; }
        #undef NEXT
    }
    if (!kmax) kmax = kmin;
    if (D0 % 3) { fprintf(stderr, "FATAL: 3 must divide D0\n"); return 2; }
    if (b2 > 21845) { fprintf(stderr, "FATAL: --b2 > 21845 would overflow the "
                              "32-bit limb fold at 16-bit limbs\n"); return 2; }
    if (nomp > 0) omp_set_num_threads(nomp);
    signal(SIGINT, on_signal); signal(SIGTERM, on_signal);
    lc_build_small_primes(b2 > 10000 ? b2 : 10000);
    const u64 n_sp = lc_n_sp();

    /* ---- resume: skip units already recorded in --out ---- */
    std::set<std::string> done_units;
    if (resume && out) {
        FILE *rf = fopen(out, "r");
        if (rf) {
            char line[4096];
            while (fgets(line, sizeof line, rf)) {
                if (!strstr(line, "\"covered\"")) continue;
                unsigned long long K, sh;
                const char *pk = strstr(line, "\"K\":"), *ps = strstr(line, "\"shift\":");
                if (pk && ps && sscanf(pk, "\"K\":%llu", &K) == 1 &&
                    sscanf(ps, "\"shift\":%llu", &sh) == 1)
                    done_units.insert(unit_key(K, sh));
            }
            fclose(rf);
            fprintf(stderr, "resume: %zu units already recorded in %s\n",
                    done_units.size(), out);
        }
    }

    /* ---- CUDA setup ---- */
    int dev = 0; cudaDeviceProp prop; CK(cudaGetDeviceProperties(&prop, dev));
    fprintf(stderr, "GPU: %s (SMs=%d, cc %d.%d) omp_threads=%d\n", prop.name,
            prop.multiProcessorCount, prop.major, prop.minor, omp_get_max_threads());
    cudaStream_t stream[2]; CK(cudaStreamCreate(&stream[0])); CK(cudaStreamCreate(&stream[1]));
    const u32 CAP = 1u << 18;
    u64 *d_R[2] = {NULL, NULL}, *d_W[2] = {NULL, NULL}, *d_H[2];
    u32 *d_cnt[2];
    size_t cap_R[2] = {0, 0}, cap_W[2] = {0, 0};
    for (int b = 0; b < 2; b++) {
        CK(cudaMalloc(&d_H[b], (size_t)CAP * 16));
        CK(cudaMalloc(&d_cnt[b], 4));
    }
    std::vector<u64> h_hits((size_t)CAP * 2);
    Unit inflight[2]; inflight[0].valid = inflight[1].valid = false;
    for (int b = 0; b < 2; b++) CK(cudaEventCreate(&inflight[b].done));

    FILE *of = out ? fopen(out, "a") : NULL;
    if (out && !of) { perror("open --out"); return 2; }

    int global_best = 0;
    double t0 = lc_now_s(), covered = 0;
    u64 units = 0, total_conf = 0, total_surv = 0;

    /* host-side scratch, reused across units */
    std::vector<u64> Rpre, words;
    std::vector<u32> offs, rps, mags, k1s, k2s, k3s;
    u64 tcp[MAXTC];

    /* ---- finish a unit: wait for its kernel, pull survivors, stage 3 ---- */
    auto finish_unit = [&](int b) {
        Unit &U = inflight[b];
        if (!U.valid) return;
        CK(cudaEventSynchronize(U.done));
        double gpu_s = lc_now_s() - U.t_launch;
        u32 nsurv = 0;
        CK(cudaMemcpyAsync(&nsurv, d_cnt[b], 4, cudaMemcpyDeviceToHost, stream[b]));
        CK(cudaStreamSynchronize(stream[b]));
        bool overflow = nsurv > CAP;
        if (overflow) {
            fprintf(stderr, "WARNING: survivor buffer overflow (%u > %u); "
                    "results for K=%llu shift=%llu are INCOMPLETE\n",
                    nsurv, CAP, (unsigned long long)U.K, (unsigned long long)U.shift);
            nsurv = CAP;
        }
        if (nsurv)
            CK(cudaMemcpy(h_hits.data(), d_H[b], (size_t)nsurv * 16, cudaMemcpyDeviceToHost));

        /* ---- stage 3 on the CPU: exact test of each surviving window ---- */
        double ts3 = lc_now_s();
        u64 nconf = 0, nbits = 0;
        u64 boff = U.shift * 64, d = U.d, MOD = U.MOD;
        std::vector<Hit> hits;
        #pragma omp parallel for schedule(dynamic, 4) reduction(+:nconf, nbits)
        for (u32 s = 0; s < nsurv; s++) {
            u64 Rb = h_hits[2 * s], sito = h_hits[2 * s + 1];
            while (sito) {
                int bb = __builtin_ctzll(sito); sito &= sito - 1;
                nbits++;
                u64 a0 = Rb + (boff + (u64)bb) * MOD;
                u64 best_run = 0, best_start = 0, cur = 0, cur_start = 0;
                u64 want = (report > 1) ? (u64)report : 1;
                for (u64 kk = 0; kk < nterms; kk++) {
                    if (cur + (nterms - kk) < want) break;
                    u64 t2 = a0 + kk * d;
                    if (t2 < a0) break;
                    if (lc_is_loeschian(t2)) {
                        if (cur == 0) cur_start = kk;
                        cur++;
                        if (cur > best_run) { best_run = cur; best_start = cur_start; }
                    } else cur = 0;
                }
                if (!best_run) continue;
                nconf++;
                u64 a = a0 + best_start * d, run = best_run;
                while (a >= d && lc_is_loeschian(a - d)) { a -= d; run++; }
                for (;;) { u64 t2 = a + run * d;
                    if (t2 < a || !lc_is_loeschian(t2)) break; run++; }
                #pragma omp critical
                {
                    if ((int)run > global_best) {
                        global_best = (int)run;
                        fprintf(stderr, "*** n=%d a=%llu d=%llu (K=%llu shift=%llu)\n",
                                global_best, (unsigned long long)a, (unsigned long long)d,
                                (unsigned long long)U.K, (unsigned long long)U.shift);
                    }
                    if ((int)run >= report) hits.push_back(Hit{run, a, d, U.K});
                }
            }
        }
        double cpu_s = lc_now_s() - ts3;
        if (of) {
            for (const Hit &h : hits) {
                fprintf(of, "{\"hit\":true,\"n\":%llu,\"a\":%llu,\"d\":%llu,\"K\":%llu}\n",
                        (unsigned long long)h.n, (unsigned long long)h.a,
                        (unsigned long long)h.d, (unsigned long long)h.K);
            }
        }
        covered += (double)MOD * 64.0;
        units++;
        total_conf += nconf;
        total_surv += nbits;
        fprintf(stderr, "K=%llu sh=%llu MOD=%llu thr=%llu res=%.4g tierC=%d "
                "surv=%u bits=%llu conf=%llu best=%d gpu=%.2fs cpu=%.2fs (%.3g res/s)\n",
                (unsigned long long)U.K, (unsigned long long)U.shift,
                (unsigned long long)MOD, (unsigned long long)U.nthreads,
                (double)U.total, U.ntc, nsurv, (unsigned long long)nbits,
                (unsigned long long)nconf, global_best, gpu_s, cpu_s,
                (double)U.total / gpu_s);
        if (of) {
            fprintf(of, "{\"K\":%llu,\"shift\":%llu,\"res\":%.6g,\"surv\":%u,"
                    "\"conf\":%llu,\"best\":%d,\"gpu_s\":%.2f,\"cpu_s\":%.2f,"
                    "\"covered\":%.6g,\"overflow\":%s}\n",
                    (unsigned long long)U.K, (unsigned long long)U.shift,
                    (double)U.total, nsurv, (unsigned long long)nconf, global_best,
                    gpu_s, cpu_s, covered, overflow ? "true" : "false");
            fflush(of);
        }
        if (verify_mode)
            printf("VERIFY K=%llu shift=%llu residues=%.6g survivors=%u bits=%llu "
                   "confirmed=%llu\n", (unsigned long long)U.K,
                   (unsigned long long)U.shift, (double)U.total, nsurv,
                   (unsigned long long)nbits, (unsigned long long)nconf);
        U.valid = false;
    };

    int buf = 0;
    for (u64 K = kmin; K <= kmax && !stop_requested; K++) {
        u64 d = K * D0;
        if (d / D0 != K) { fprintf(stderr, "K overflow at %llu\n", (unsigned long long)K); break; }

        /* ---- stage-1 components (identical rules to apsearch.c) ---- */
        comp C[MAXC]; int nc = 0; u64 MOD = 1;
        #define ADDC(m_, s_, t_, c_) do { C[nc].m=(m_); C[nc].s=(s_); \
            C[nc].t=(t_); C[nc].c=(c_); MOD *= (m_); nc++; } while (0)
        ADDC(3, 1, 0, 1);
        if (D0 % 2 == 0) ADDC(2, 1, 0, 1);
        if (D0 % 5 == 0) ADDC(5, 1, 1, 4);
        int ntb = 0;
        for (u64 i = 0; i < n_sp && nc < MAXC - 1; i++) {
            u64 q = lc_sp(i);
            if (q <= nterms || q % 3 != 2 || d % q == 0) continue;
            if (MOD > (u64)4e18 / q || MOD * q > modcap) break;
            u64 dq = d % q;
            ADDC(q, dq, dq, q - nterms);
            ntb++;
        }
        #undef ADDC
        if (ntb < 4) { fprintf(stderr, "K=%llu: too few tier-B primes\n", (unsigned long long)K); continue; }

        /* ---- CRT ---- */
        u64 R0 = 0, sstep[MAXC], subcyc[MAXC];
        for (int i = 0; i < nc; i++) {
            u64 mi = C[i].m, co = MOD / mi;
            u64 e = lc_mulmod(co % MOD, lc_inv_mod(co % mi, mi), MOD);
            R0 = (R0 + lc_mulmod(C[i].s, e, MOD)) % MOD;
            sstep[i] = lc_mulmod(C[i].t, e, MOD);
        }
        for (int i = 0; i < nc; i++) subcyc[i] = lc_mulmod(C[i].c, sstep[i], MOD);
        for (int i = 0; i < nc; i++) {                     /* CRT self-check */
            if (R0 % C[i].m != C[i].s % C[i].m) {
                fprintf(stderr, "FATAL CRT R0 mod %llu\n", (unsigned long long)C[i].m); return 2; }
            for (int j = 0; j < nc; j++)
                if (sstep[i] % C[j].m != ((i == j) ? C[i].t % C[j].m : 0)) {
                    fprintf(stderr, "FATAL CRT sstep\n"); return 2; }
        }

        /* ---- the two innermost (largest-count) components run on-thread ---- */
        int in1 = -1, in2 = -1;
        for (int i = 0; i < nc; i++) {
            if (C[i].c <= 1) continue;
            if (in2 < 0 || C[i].c > C[in2].c) { in1 = in2; in2 = i; }
            else if (in1 < 0 || C[i].c > C[in1].c) { in1 = i; }
        }
        if (in1 < 0 || in2 < 0) { fprintf(stderr, "K=%llu: no inner loops\n", (unsigned long long)K); continue; }

        u64 nthreads = 1, total = 1;
        for (int i = 0; i < nc; i++) {
            total *= C[i].c;
            if (i != in1 && i != in2) nthreads *= C[i].c;
        }

        /* ---- per-thread start residues (additive odometer over the prefix) ---- */
        Rpre.resize(nthreads);
        {
            int pidx[MAXC]; memset(pidx, 0, sizeof(pidx));
            u64 R = R0, n = 0;
            for (;;) {
                Rpre[n++] = R;
                if (n >= nthreads) break;
                int i = nc - 1; bool done = false;
                for (;;) {
                    while (i >= 0 && (C[i].c <= 1 || i == in1 || i == in2)) i--;
                    if (i < 0) { done = true; break; }
                    R += sstep[i]; if (R >= MOD) R -= MOD;
                    if ((u64)(++pidx[i]) < C[i].c) break;
                    R = (R + MOD - subcyc[i]) % MOD;
                    pidx[i] = 0; i--;
                }
                if (done) break;
            }
            if (n != nthreads) {
                fprintf(stderr, "FATAL: prefix count %llu != %llu\n",
                        (unsigned long long)n, (unsigned long long)nthreads);
                return 2;
            }
        }

        /* ---- tier C tables, ordered by kill power ---- */
        int ntc = 0;
        for (u64 i = 0; i < n_sp && ntc < MAXTC; i++) {
            u64 r = lc_sp(i);
            if (r > b2) break;
            bool pinned = false;
            for (int t = 0; t < nc; t++) if (C[t].m == r) pinned = true;
            if (pinned) continue;
            bool divD = (D0 % r == 0);
            if (!divD && (r % 3 != 2 || d % r == 0 || r <= nterms)) continue;
            tcp[ntc++] = r;
        }
        for (int i = 0; i < ntc; i++) {
            int bj = i; double bk = -1;
            for (int j = i; j < ntc; j++) {
                double kill = (D0 % tcp[j] == 0) ? 1.0 / (double)tcp[j]
                                                 : (double)nterms / (double)tcp[j];
                if (kill > bk) { bk = kill; bj = j; }
            }
            u64 tt = tcp[i]; tcp[i] = tcp[bj]; tcp[bj] = tt;
        }
        u64 tot_words = 0;
        for (int i = 0; i < ntc; i++) tot_words += tcp[i];
        offs.resize(ntc); rps.resize(ntc); mags.resize(ntc);
        k1s.resize(ntc); k2s.resize(ntc); k3s.resize(ntc);
        words.resize(tot_words);
        {
            u64 acc = 0;
            for (int t = 0; t < ntc; t++) {
                u64 r = tcp[t];
                offs[t] = (u32)acc; rps[t] = (u32)r;
                mags[t] = (u32)(((u64)1 << 32) / r);
                k1s[t] = (u32)(((u64)1 << 16) % r);
                k2s[t] = (u32)(((u64)1 << 32) % r);
                k3s[t] = (u32)(((u64)1 << 48) % r);
                acc += r;
            }
            if (K == kmin) {
                /* magic-modulo self-check against %, for every prime, on a
                 * spread of x including the extremes of the fold range */
                for (int t = 0; t < ntc; t++) {
                    u32 p = rps[t], M = mags[t];
                    for (u64 xx = 0; xx < 200000; xx++) {
                        u32 x = (u32)((xx * 2654435761ULL) ^ (xx << 7));
                        if (xx < 8) x = (u32)(0xFFFFFFFFu - xx);
                        u32 q = (u32)(((u64)x * M) >> 32), rr = x - q * p;
                        if (rr >= p) rr -= p;
                        if (rr >= p) rr -= p;
                        if (rr != x % p) { fprintf(stderr, "FATAL magic mod p=%u x=%u\n", p, x); return 2; }
                    }
                }
                fprintf(stderr, "tierC order: %llu %llu %llu %llu %llu ... (%d primes, "
                        "%llu words) magic-mod OK\n", (unsigned long long)tcp[0],
                        (unsigned long long)tcp[1], (unsigned long long)tcp[2],
                        (unsigned long long)tcp[3], (unsigned long long)tcp[4], ntc,
                        (unsigned long long)tot_words);
            }
        }

        for (u64 shift = shift0; shift < shift0 + shifts && !stop_requested; shift++) {
            if (!done_units.empty() && done_units.count(unit_key(K, shift))) continue;
            u64 boff = shift * 64;
            #pragma omp parallel for schedule(dynamic, 1)
            for (int t = 0; t < ntc; t++) {
                u64 r = tcp[t], acc = offs[t];
                std::vector<char> ok(r, 0);
                if (D0 % r == 0) { for (u64 x = 0; x < r; x++) ok[x] = (x != 0); }
                else {
                    for (u64 x = 0; x < r; x++) ok[x] = 1;
                    u64 v = 0, dm = d % r;
                    for (u64 kk = 0; kk < nterms; kk++) { ok[v] = 0; v = (v + r - dm) % r; }
                }
                u64 modr = MOD % r;
                u64 b0 = lc_mulmod(boff % r, modr, r);
                for (u64 x = 0; x < r; x++) {
                    u64 w = 0, y = (x + b0) % r;
                    for (int b = 0; b < 64; b++) {
                        if (ok[y]) w |= (u64)1 << b;
                        y += modr; if (y >= r) y -= r;
                    }
                    words[acc + x] = w;
                }
            }

            /* ---- buffer set `buf` is free once the unit before last is finished ---- */
            finish_unit(buf);
            if (cap_R[buf] < nthreads) {
                if (d_R[buf]) CK(cudaFree(d_R[buf]));
                CK(cudaMalloc(&d_R[buf], nthreads * 8)); cap_R[buf] = nthreads;
            }
            if (cap_W[buf] < tot_words) {
                if (d_W[buf]) CK(cudaFree(d_W[buf]));
                CK(cudaMalloc(&d_W[buf], tot_words * 8)); cap_W[buf] = tot_words;
            }
            /* constants are shared by both streams: make sure the other
             * stream's kernel is done before overwriting them */
            if (inflight[buf ^ 1].valid) CK(cudaEventSynchronize(inflight[buf ^ 1].done));
            CK(cudaMemcpyToSymbol(c_off, offs.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_rp, rps.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_mag, mags.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_k1, k1s.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_k2, k2s.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_k3, k3s.data(), ntc * 4));
            CK(cudaMemcpyAsync(d_R[buf], Rpre.data(), nthreads * 8, cudaMemcpyHostToDevice, stream[buf]));
            CK(cudaMemcpyAsync(d_W[buf], words.data(), tot_words * 8, cudaMemcpyHostToDevice, stream[buf]));
            CK(cudaMemsetAsync(d_cnt[buf], 0, 4, stream[buf]));

            Params P;
            P.MOD = MOD; P.s1 = sstep[in1]; P.s2 = sstep[in2];
            P.c1 = (u32)C[in1].c; P.c2 = (u32)C[in2].c; P.ntc = (u32)ntc;
            P.cap = CAP; P.nthreads = (u32)nthreads;
            Unit &U = inflight[buf];
            U.K = K; U.shift = shift; U.d = d; U.MOD = MOD; U.total = total;
            U.nthreads = nthreads; U.ntc = ntc; U.valid = true;
            U.t_launch = lc_now_s();
            u32 blocks = (u32)((nthreads + 255) / 256);
            sieve<<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
            CK(cudaGetLastError());
            CK(cudaEventRecord(U.done, stream[buf]));
            buf ^= 1;
        }
    }
    finish_unit(buf);
    finish_unit(buf ^ 1);

    if (of) fclose(of);
    printf("GPU BEST n=%d  units=%llu covered=%.4g raw a  survivors=%llu confirmed=%llu  %.1fs%s\n",
           global_best, (unsigned long long)units, covered,
           (unsigned long long)total_surv, (unsigned long long)total_conf,
           lc_now_s() - t0, stop_requested ? " INTERRUPTED" : "");
    return 0;
}

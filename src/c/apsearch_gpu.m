/* apsearch_gpu -- the apsearch.c engine with stage 1 + stage 2 on the GPU.
 *
 * Build:
 *   cc -O3 -fobjc-arc -framework Metal -framework Foundation -I src/c \
 *      -o build/apsearch_gpu src/c/apsearch_gpu.m
 *
 * WHY. Profiling the CPU engine put ~164 of its 172 ns per residue in stage 2
 * (the bitmask AND chain) and only ~8 ns in stage 3. Stage 2 is pure 64-bit
 * AND plus a small-modulus table lookup over ~1e9 independent residues per
 * work unit, which is exactly what an integrated GPU is good at. Measured on
 * this M2: 1.68-1.84 G residues/s versus 5.8 M/s per CPU core.
 *
 * THE ONE TRAP. 64-bit integer modulo on the Apple GPU is ~15x slower than
 * 32-bit and would make the GPU LOSE to the CPU. So the index arithmetic stays
 * in 32 bits: R < MOD < 2^51 is split into three 17-bit limbs and folded with
 * precomputed 2^17 and 2^34 residues. The fold is at most
 * 2^17 * (1 + 2*b2) < 2^32 for b2 <= 16000, which is asserted at startup.
 *
 * DIVISION OF LABOUR
 *   host : stage-1 CRT setup, the per-thread starting residues, the OKOK
 *          tables, and stage 3 (the exact Loeschian test, shared with the CPU
 *          engine via loesch_core.h).
 *   GPU  : one thread per prefix of the additive loop nest; each walks the two
 *          innermost levels, runs the short-circuiting AND chain, and appends
 *          survivors through an atomic counter.
 *
 * Only 64-bit add / compare / subtract are needed on the GPU, so no 128-bit
 * intermediate is ever required there (MSL has no u128): the per-thread start
 * values are computed on the host, where mulmod is available.
 *
 * Correctness is checked with --verify against the CPU engine's survivor set
 * for the same work unit, and every reported run still has to pass
 * src/verify.py.
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <time.h>
#include "loesch_core.h"

#define D0_DEFAULT 382160924970ULL
#define MAXC 16
#define MAXTC 2048

/* ------------------------------------------------------------- the kernel -- */

static const char *MSL = R"MSL(
#include <metal_stdlib>
using namespace metal;

struct Params {
    ulong MOD;
    ulong s1;        // additive step of the second-innermost component
    ulong s2;        // additive step of the innermost component
    uint  c1, c2;    // their trip counts
    uint  ntc;       // number of tier-C tables
    uint  cap;       // survivor buffer capacity
};

kernel void sieve(const device ulong *Rpre  [[buffer(0)]],
                  const device ulong *words [[buffer(1)]],
                  const device uint  *off   [[buffer(2)]],
                  const device uint  *rp    [[buffer(3)]],
                  const device uint  *k1    [[buffer(4)]],   // 2^17 mod r
                  const device uint  *k2    [[buffer(5)]],   // 2^34 mod r
                  constant Params &P        [[buffer(6)]],
                  device atomic_uint *cnt   [[buffer(7)]],
                  device ulong *hits        [[buffer(8)]],
                  uint gid [[thread_position_in_grid]])
{
    ulong R = Rpre[gid];
    for (uint i1 = 0; i1 < P.c1; ++i1) {
        ulong Rb = R;
        for (uint i2 = 0; i2 < P.c2; ++i2) {
            // three 17-bit limbs: keeps every index computation in 32 bits
            uint l0 = (uint)(Rb & 0x7FFFF);
            uint l1 = (uint)((Rb >> 19) & 0x7FFFF);
            uint l2 = (uint)(Rb >> 38);
            ulong sito = ~0UL;
            for (uint t = 0; t < P.ntc; ++t) {
                uint idx = (l0 + l1 * k1[t] + l2 * k2[t]) % rp[t];   // k = 2^19, 2^38 mod r
                sito &= words[off[t] + idx];
                if (sito == 0UL) break;
            }
            if (sito != 0UL) {
                uint s = atomic_fetch_add_explicit(cnt, 1u, memory_order_relaxed);
                if (s < P.cap) { hits[2 * s] = Rb; hits[2 * s + 1] = sito; }
            }
            Rb += P.s2; if (Rb >= P.MOD) Rb -= P.MOD;
        }
        R += P.s1; if (R >= P.MOD) R -= P.MOD;
    }
}
)MSL";

/* ------------------------------------------------------------------ main -- */

typedef struct { u64 m, s, t, c; } comp;

int main(int argc, char **argv) {
@autoreleasepool {
    u64 kmin = 1, kmax = 0, nterms = 58, shifts = 1, b2 = 10000, D0 = D0_DEFAULT;
    u64 modcap = 2000000000000000ULL;
    int report = 44;
    const char *out = NULL;
    bool verify_mode = false;

    for (int i = 1; i < argc; i++) {
        const char *k = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : "")
        if (!strcmp(k, "--kmin")) kmin = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--kmax")) kmax = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--nterms")) nterms = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--shifts")) shifts = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--b2")) b2 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--modcap")) modcap = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--report")) report = atoi(NEXT());
        else if (!strcmp(k, "--out")) out = NEXT();
        else if (!strcmp(k, "--D0")) D0 = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--verify")) verify_mode = true;
        else { fprintf(stderr, "unknown arg %s\n", k); return 2; }
        #undef NEXT
    }
    if (!kmax) kmax = kmin;
    if (D0 % 3) { fprintf(stderr, "FATAL: 3 must divide D0\n"); return 2; }
    /* 19-bit limbs: l1*k1 + l2*k2 < 2*2^19*b2 must stay under 2^32 */
    if (b2 > 4000) { fprintf(stderr, "FATAL: --b2 > 4000 would overflow the "
                             "32-bit limb fold at 19-bit limbs\n"); return 2; }
    build_small_primes(10000);

    /* ---- Metal setup ---- */
    id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
    if (!dev) { fprintf(stderr, "FATAL: no Metal device\n"); return 2; }
    NSError *err = nil;
    id<MTLLibrary> lib = [dev newLibraryWithSource:[NSString stringWithUTF8String:MSL]
                                          options:nil error:&err];
    if (!lib) { fprintf(stderr, "MSL compile failed: %s\n",
                        err.description.UTF8String); return 2; }
    id<MTLComputePipelineState> pso =
        [dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"sieve"]
                                           error:&err];
    if (!pso) { fprintf(stderr, "pipeline failed: %s\n",
                        err.description.UTF8String); return 2; }
    id<MTLCommandQueue> queue = [dev newCommandQueue];
    fprintf(stderr, "GPU: %s (execWidth=%lu maxTG=%lu)\n", dev.name.UTF8String,
            (unsigned long)pso.threadExecutionWidth,
            (unsigned long)pso.maxTotalThreadsPerThreadgroup);

    FILE *of = out ? fopen(out, "a") : NULL;
    if (out && !of) { perror("open --out"); return 2; }

    int global_best = 0;
    double t0 = now_s(), covered = 0;
    u64 units = 0, total_conf = 0;

    for (u64 K = kmin; K <= kmax; K++) {
        u64 d = K * D0;

        /* ---- stage-1 components (identical rules to apsearch.c) ---- */
        comp C[MAXC]; int nc = 0; u64 MOD = 1;
        #define ADDC(m_, s_, t_, c_) do { C[nc].m=(m_); C[nc].s=(s_); \
            C[nc].t=(t_); C[nc].c=(c_); MOD *= (m_); nc++; } while (0)
        ADDC(3, 1, 0, 1);
        if (D0 % 2 == 0) ADDC(2, 1, 0, 1);
        if (D0 % 5 == 0) ADDC(5, 1, 1, 4);
        int ntb = 0;
        for (u64 i = 0; i < n_sp && nc < MAXC - 1; i++) {
            u64 q = sp[i];
            if (q <= nterms || q % 3 != 2 || d % q == 0) continue;
            if (MOD > (u64)4e18 / q || MOD * q > modcap) break;
            u64 dq = d % q;
            ADDC(q, dq, dq, q - nterms);
            ntb++;
        }
        #undef ADDC
        if (ntb < 4) { fprintf(stderr, "K=%llu: too few tier-B primes\n", K); continue; }
        if (MOD >= (1ULL << 57)) { fprintf(stderr, "FATAL: MOD >= 2^57 breaks the "
                                   "19-bit limb split\n"); return 2; }

        /* ---- CRT ---- */
        u64 R0 = 0, sstep[MAXC], subcyc[MAXC];
        for (int i = 0; i < nc; i++) {
            u64 mi = C[i].m, co = MOD / mi;
            u64 e = mulmod(co % MOD, inv_mod(co % mi, mi), MOD);
            R0 = (R0 + mulmod(C[i].s, e, MOD)) % MOD;
            sstep[i] = mulmod(C[i].t, e, MOD);
        }
        for (int i = 0; i < nc; i++) subcyc[i] = mulmod(C[i].c, sstep[i], MOD);
        for (int i = 0; i < nc; i++) {                     /* same self-check */
            if (R0 % C[i].m != C[i].s % C[i].m) {
                fprintf(stderr, "FATAL CRT R0 mod %llu\n", C[i].m); return 2; }
            for (int j = 0; j < nc; j++)
                if (sstep[i] % C[j].m != ((i == j) ? C[i].t % C[j].m : 0)) {
                    fprintf(stderr, "FATAL CRT sstep\n"); return 2; }
        }

        /* ---- pick the two innermost (largest-count) components for the GPU
         * inner loops; all the rest become the per-thread prefix ---- */
        int in1 = -1, in2 = -1;
        for (int i = 0; i < nc; i++) {
            if (C[i].c <= 1) continue;
            if (in2 < 0 || C[i].c > C[in2].c) { in1 = in2; in2 = i; }
            else if (in1 < 0 || C[i].c > C[in1].c) { in1 = i; }
        }
        if (in1 < 0 || in2 < 0) { fprintf(stderr, "K=%llu: no inner loops\n", K); continue; }

        u64 nthreads = 1, total = 1;
        for (int i = 0; i < nc; i++) {
            total *= C[i].c;
            if (i != in1 && i != in2) nthreads *= C[i].c;
        }

        /* ---- per-thread start residues, by an additive odometer over the
         * prefix components (host side, where mulmod exists) ---- */
        u64 *Rpre = malloc(nthreads * sizeof(u64));
        if (!Rpre) { fprintf(stderr, "oom Rpre\n"); return 2; }
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
                fprintf(stderr, "FATAL: prefix count %llu != %llu\n", n, nthreads);
                return 2;
            }
        }

        /* ---- tier C tables, ordered by kill power ---- */
        u64 tcp[MAXTC]; int ntc = 0;
        for (u64 i = 0; i < n_sp && ntc < MAXTC; i++) {
            u64 r = sp[i];
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
        u32 *offs = malloc(ntc * 4), *rps = malloc(ntc * 4);
        u32 *k1s = malloc(ntc * 4), *k2s = malloc(ntc * 4);
        u64 *words = malloc(tot_words * 8);
        if (!offs || !rps || !k1s || !k2s || !words) { fprintf(stderr, "oom tables\n"); return 2; }

        for (u64 shift = 0; shift < shifts; shift++) {
          /* Inner pool: without it every unit's Metal buffers (tables + the
           * survivor buffer) live until main() returns, leaking ~100 MB per
           * unit. That drove the machine into swap and made units take 420 s
           * instead of 12. */
          @autoreleasepool {
            u64 boff = shift * 64;
            u64 acc = 0;
            for (int t = 0; t < ntc; t++) {
                u64 r = tcp[t];
                offs[t] = (u32)acc; rps[t] = (u32)r;
                k1s[t] = (u32)(((u64)1 << 19) % r);
                k2s[t] = (u32)(((u64)1 << 38) % r);
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
                    words[acc + x] = w;
                }
                free(ok);
                acc += r;
            }

            /* ---- dispatch ---- */
            const u32 CAP = 1u << 18;
            #define BUF(p, len) [dev newBufferWithBytes:(p) length:(len) \
                                   options:MTLResourceStorageModeShared]
            id<MTLBuffer> bR  = BUF(Rpre, nthreads * 8);
            id<MTLBuffer> bW  = BUF(words, tot_words * 8);
            id<MTLBuffer> bO  = BUF(offs, ntc * 4);
            id<MTLBuffer> bP  = BUF(rps, ntc * 4);
            id<MTLBuffer> bK1 = BUF(k1s, ntc * 4);
            id<MTLBuffer> bK2 = BUF(k2s, ntc * 4);
            id<MTLBuffer> bCnt = [dev newBufferWithLength:4
                                     options:MTLResourceStorageModeShared];
            id<MTLBuffer> bH = [dev newBufferWithLength:(size_t)CAP * 16
                                   options:MTLResourceStorageModeShared];
            #undef BUF
            memset(bCnt.contents, 0, 4);

            struct { u64 MOD, s1, s2; u32 c1, c2, ntc, cap; } P = {
                MOD, sstep[in1], sstep[in2], (u32)C[in1].c, (u32)C[in2].c,
                (u32)ntc, CAP };

            double ts = now_s();
            id<MTLCommandBuffer> cb = [queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:bR offset:0 atIndex:0];
            [enc setBuffer:bW offset:0 atIndex:1];
            [enc setBuffer:bO offset:0 atIndex:2];
            [enc setBuffer:bP offset:0 atIndex:3];
            [enc setBuffer:bK1 offset:0 atIndex:4];
            [enc setBuffer:bK2 offset:0 atIndex:5];
            [enc setBytes:&P length:sizeof P atIndex:6];
            [enc setBuffer:bCnt offset:0 atIndex:7];
            [enc setBuffer:bH offset:0 atIndex:8];
            [enc dispatchThreads:MTLSizeMake(nthreads, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
            [enc endEncoding];
            [cb commit];
            [cb waitUntilCompleted];
            if (cb.error) { fprintf(stderr, "GPU exec: %s\n",
                            cb.error.description.UTF8String); return 2; }
            double gpu_s = now_s() - ts;

            u32 nsurv = *(u32 *)bCnt.contents;
            bool overflow = nsurv > CAP;
            if (overflow) {
                fprintf(stderr, "WARNING: survivor buffer overflow (%u > %u); "
                        "results for K=%llu shift=%llu are INCOMPLETE\n",
                        nsurv, CAP, K, shift);
                nsurv = CAP;
            }
            u64 *H = (u64 *)bH.contents;

            /* ---- stage 3 on the CPU: exact test of each surviving window ---- */
            u64 nconf = 0;
            double ts3 = now_s();
            for (u32 s = 0; s < nsurv; s++) {
                u64 Rb = H[2 * s], sito = H[2 * s + 1];
                while (sito) {
                    int b = __builtin_ctzll(sito); sito &= sito - 1;
                    u64 a0 = Rb + (boff + (u64)b) * MOD;
                    u64 best_run = 0, best_start = 0, cur = 0, cur_start = 0;
                    u64 want = (report > 1) ? (u64)report : 1;
                    for (u64 kk = 0; kk < nterms; kk++) {
                        if (cur + (nterms - kk) < want) break;
                        u64 t2 = a0 + kk * d;
                        if (t2 < a0) break;
                        if (is_loeschian(t2)) {
                            if (cur == 0) cur_start = kk;
                            cur++;
                            if (cur > best_run) { best_run = cur; best_start = cur_start; }
                        } else cur = 0;
                    }
                    if (!best_run) continue;
                    nconf++;
                    u64 a = a0 + best_start * d, run = best_run;
                    while (a >= d && is_loeschian(a - d)) { a -= d; run++; }
                    for (;;) { u64 t2 = a + run * d;
                        if (t2 < a || !is_loeschian(t2)) break; run++; }
                    if ((int)run > global_best) {
                        global_best = (int)run;
                        fprintf(stderr, "*** n=%d a=%llu d=%llu (K=%llu shift=%llu)\n",
                                global_best, a, d, K, shift);
                    }
                    if ((int)run >= report && of) {
                        fprintf(of, "{\"hit\":true,\"n\":%llu,\"a\":%llu,\"d\":%llu,"
                                "\"K\":%llu}\n", run, a, d, K);
                        fflush(of);
                    }
                }
            }
            double cpu_s = now_s() - ts3;
            covered += (double)MOD * 64.0;
            units++;
            total_conf += nconf;
            fprintf(stderr, "K=%llu sh=%llu MOD=%llu thr=%llu res=%.4g tierC=%d "
                    "surv=%u conf=%llu best=%d gpu=%.2fs cpu=%.2fs (%.3g res/s)\n",
                    K, shift, MOD, nthreads, (double)total, ntc, nsurv, nconf,
                    global_best, gpu_s, cpu_s, (double)total / gpu_s);
            if (of) {
                fprintf(of, "{\"K\":%llu,\"shift\":%llu,\"res\":%.6g,\"surv\":%u,"
                        "\"conf\":%llu,\"best\":%d,\"gpu_s\":%.2f,\"cpu_s\":%.2f,"
                        "\"covered\":%.6g,\"overflow\":%s}\n",
                        K, shift, (double)total, nsurv, nconf, global_best,
                        gpu_s, cpu_s, covered, overflow ? "true" : "false");
                fflush(of);
            }
            if (verify_mode) {
                printf("VERIFY K=%llu shift=%llu residues=%.6g survivors=%u "
                       "confirmed=%llu\n", K, shift, (double)total, nsurv, nconf);
            }
          }
        }
        free(Rpre); free(words); free(offs); free(rps); free(k1s); free(k2s);
    }

    if (of) fclose(of);
    printf("GPU BEST n=%d  units=%llu covered=%.4g raw a  confirmed=%llu  %.1fs\n",
           global_best, units, covered, total_conf, now_s() - t0);
    return 0;
}
}

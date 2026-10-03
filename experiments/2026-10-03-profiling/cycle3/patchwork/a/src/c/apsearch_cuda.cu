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
#include <deque>
#include <thread>
#include <memory>
#include <unistd.h>
#include <mutex>
#include <condition_variable>
#include <set>
#include <string>
#include <algorithm>

#include "loesch_core_api.h"

typedef uint64_t u64;
typedef uint32_t u32;
typedef __uint128_t u128;

#define D0_DEFAULT 382160924970ULL   /* 3 * 2*5*11*17*23*29*41*47*53 */
#define MAXC 16
#define MAXTC 2048                   /* tier-C primes; 7 constant arrays must fit 64 KB */

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


/* Variant 1: FLAT chain. In the nested version a warp runs until its slowest
 * lane's AND chain dies, so with an average chain of ~15 primes and a
 * max-of-32 of ~40 most lanes idle most of the time. Here every lane walks
 * its own residues and its own chain position t independently: when its
 * word dies it immediately moves to its next residue and restarts at t = 0.
 * The warp only synchronises on "all lanes finished". The price is that t is
 * no longer warp-uniform, so per-prime constants come from shared memory
 * (one 16-byte load per step) instead of the constant cache. */
#define MAXTC_SH 1024
__global__ void __launch_bounds__(256)
sieve_flat(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
           Params P, u32 *cnt, u64 *hits)
{
    __shared__ uint4 s_c[MAXTC_SH];     /* {rp | k1<<16, k2 | k3<<16, mag, off} */
    for (u32 t = threadIdx.x; t < P.ntc; t += blockDim.x)
        s_c[t] = make_uint4(c_rp[t] | (c_k1[t] << 16), c_k2[t] | (c_k3[t] << 16),
                            c_mag[t], c_off[t]);
    __syncthreads();
    u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    bool active = gid < P.nthreads;
    u64 R = active ? Rpre[gid] : 0, Rb = R;
    u32 i1 = 0, i2 = 0, t = 0;
    u64 sito = ~0ULL;
    u32 l0 = (u32)(Rb & 0xFFFF), l1 = (u32)((Rb >> 16) & 0xFFFF);
    u32 l2 = (u32)((Rb >> 32) & 0xFFFF), l3 = (u32)(Rb >> 48);
    const u32 ntc = P.ntc;
    for (;;) {
        if (active) {
            uint4 c = s_c[t];
            u32 p = c.x & 0xFFFF;
            u32 x = l0 + l1 * (c.x >> 16) + l2 * (c.y & 0xFFFF) + l3 * (c.y >> 16);
            u32 q = __umulhi(x, c.z);
            u32 r = x - q * p;
            if (r >= p) r -= p;
            if (r >= p) r -= p;
            sito &= words[c.w + r];
            t++;
            if (sito == 0ULL || t == ntc) {
                if (sito != 0ULL) {
                    u32 s = atomicAdd(cnt, 1u);
                    if (s < P.cap) { hits[2 * s] = Rb; hits[2 * s + 1] = sito; }
                }
                /* next residue */
                if (++i2 < P.c2) { Rb += P.s2; if (Rb >= P.MOD) Rb -= P.MOD; }
                else {
                    i2 = 0;
                    if (++i1 < P.c1) { R += P.s1; if (R >= P.MOD) R -= P.MOD; Rb = R; }
                    else active = false;
                }
                t = 0; sito = ~0ULL;
                l0 = (u32)(Rb & 0xFFFF); l1 = (u32)((Rb >> 16) & 0xFFFF);
                l2 = (u32)((Rb >> 32) & 0xFFFF); l3 = (u32)(Rb >> 48);
            }
        }
        if (__all_sync(0xFFFFFFFFu, !active)) break;
    }
}

/* Variant 2: nested chain, but the tables of the first `nsh` tier-C primes
 * (in kill order, i.e. the ones nearly every chain step touches) live in
 * shared memory. A warp's 32 random 8-byte gathers into a 1-3 KB table cost
 * ~14 L1 wavefronts from global/L1 but only a few bank-conflict cycles from
 * shared memory. */
__global__ void __launch_bounds__(256)
sieve_sh(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
         Params P, u32 nsh, u32 nshw, u32 *cnt, u64 *hits)
{
    extern __shared__ u64 s_w[];
    for (u32 i = threadIdx.x; i < nshw; i += blockDim.x) s_w[i] = words[i];
    __syncthreads();
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
            u32 t = 0;
            for (; t < nsh; ++t) {
                u32 x = l0 + l1 * c_k1[t] + l2 * c_k2[t] + l3 * c_k3[t];
                u32 p = c_rp[t];
                u32 q = __umulhi(x, c_mag[t]);
                u32 r = x - q * p;
                if (r >= p) r -= p;
                if (r >= p) r -= p;
                sito &= s_w[c_off[t] + r];
                if (sito == 0ULL) break;
            }
            if (sito != 0ULL) {
                for (t = nsh; t < P.ntc; ++t) {
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
            }
            Rb += P.s2; if (Rb >= P.MOD) Rb -= P.MOD;
        }
        R += P.s1; if (R >= P.MOD) R -= P.MOD;
    }
}

/* Variant 3: two independent chains per thread (residues i2 and i2+1 of the
 * innermost loop walk together). Each chain step is a dependent L1 load, so
 * a single chain per thread leaves the SM waiting on latency; two chains
 * double the loads in flight per warp at the cost of ~12 registers. */
__global__ void __launch_bounds__(256)
sieve_ilp2(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
           Params P, u32 *cnt, u64 *hits)
{
    u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= P.nthreads) return;
    u64 R = Rpre[gid];
    const u64 s2x2 = (P.s2 * 2 >= P.MOD) ? P.s2 * 2 - P.MOD : P.s2 * 2;
    for (u32 i1 = 0; i1 < P.c1; ++i1) {
        u64 Ra = R;
        u64 Rb = Ra + P.s2; if (Rb >= P.MOD) Rb -= P.MOD;
        for (u32 i2 = 0; i2 < P.c2; i2 += 2) {
            const bool two = (i2 + 1 < P.c2);
            u32 a0 = (u32)(Ra & 0xFFFF), a1 = (u32)((Ra >> 16) & 0xFFFF);
            u32 a2 = (u32)((Ra >> 32) & 0xFFFF), a3 = (u32)(Ra >> 48);
            u32 b0 = (u32)(Rb & 0xFFFF), b1 = (u32)((Rb >> 16) & 0xFFFF);
            u32 b2 = (u32)((Rb >> 32) & 0xFFFF), b3 = (u32)(Rb >> 48);
            u64 sa = ~0ULL, sb = two ? ~0ULL : 0ULL;
            for (u32 t = 0; t < P.ntc; ++t) {
                u32 k1 = c_k1[t], k2 = c_k2[t], k3 = c_k3[t], p = c_rp[t], m = c_mag[t], o = c_off[t];
                u32 xa = a0 + a1 * k1 + a2 * k2 + a3 * k3;
                u32 xb = b0 + b1 * k1 + b2 * k2 + b3 * k3;
                u32 ra = xa - __umulhi(xa, m) * p;
                u32 rb = xb - __umulhi(xb, m) * p;
                if (ra >= p) ra -= p;
                if (ra >= p) ra -= p;
                if (rb >= p) rb -= p;
                if (rb >= p) rb -= p;
                u64 wa = words[o + ra], wb = words[o + rb];
                sa &= wa; sb &= wb;
                if ((sa | sb) == 0ULL) break;
            }
            if (sa != 0ULL) {
                u32 s = atomicAdd(cnt, 1u);
                if (s < P.cap) { hits[2 * s] = Ra; hits[2 * s + 1] = sa; }
            }
            if (sb != 0ULL) {
                u32 s = atomicAdd(cnt, 1u);
                if (s < P.cap) { hits[2 * s] = Rb; hits[2 * s + 1] = sb; }
            }
            Ra += s2x2; if (Ra >= P.MOD) Ra -= P.MOD;
            Rb += s2x2; if (Rb >= P.MOD) Rb -= P.MOD;
        }
        R += P.s1; if (R >= P.MOD) R -= P.MOD;
    }
}

/* Variant 4+: NCH independent chains per thread (template), the generalisation
 * of sieve_ilp2. Chains are residues i2, i2+1, ..., i2+NCH-1 of the innermost
 * loop; the AND-chain loop runs until every chain's word is zero. */
template <int NCH>
__global__ void __launch_bounds__(256)
sieve_ilpN(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
           Params P, u32 *cnt, u64 *hits)
{
    u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= P.nthreads) return;
    u64 R = Rpre[gid];
    u64 sN = 0;                                  /* NCH * s2 mod MOD */
    for (int c = 0; c < NCH; c++) { sN += P.s2; if (sN >= P.MOD) sN -= P.MOD; }
    for (u32 i1 = 0; i1 < P.c1; ++i1) {
        u64 Rc[NCH];
        Rc[0] = R;
        #pragma unroll
        for (int c = 1; c < NCH; c++) { Rc[c] = Rc[c - 1] + P.s2; if (Rc[c] >= P.MOD) Rc[c] -= P.MOD; }
        for (u32 i2 = 0; i2 < P.c2; i2 += NCH) {
            u32 L0[NCH], L1[NCH], L2[NCH], L3[NCH];
            u64 sito[NCH];
            #pragma unroll
            for (int c = 0; c < NCH; c++) {
                L0[c] = (u32)(Rc[c] & 0xFFFF); L1[c] = (u32)((Rc[c] >> 16) & 0xFFFF);
                L2[c] = (u32)((Rc[c] >> 32) & 0xFFFF); L3[c] = (u32)(Rc[c] >> 48);
                sito[c] = (i2 + c < P.c2) ? ~0ULL : 0ULL;
            }
            for (u32 t = 0; t < P.ntc; ++t) {
                u32 k1 = c_k1[t], k2 = c_k2[t], k3 = c_k3[t], p = c_rp[t], m = c_mag[t], o = c_off[t];
                u64 any = 0;
                #pragma unroll
                for (int c = 0; c < NCH; c++) {
                    u32 x = L0[c] + L1[c] * k1 + L2[c] * k2 + L3[c] * k3;
                    u32 r = x - __umulhi(x, m) * p;
                    if (r >= p) r -= p;
                    if (r >= p) r -= p;
                    sito[c] &= words[o + r];
                    any |= sito[c];
                }
                if (any == 0ULL) break;
            }
            #pragma unroll
            for (int c = 0; c < NCH; c++) {
                if (sito[c] != 0ULL) {
                    u32 s = atomicAdd(cnt, 1u);
                    if (s < P.cap) { hits[2 * s] = Rc[c]; hits[2 * s + 1] = sito[c]; }
                }
                Rc[c] += sN; if (Rc[c] >= P.MOD) Rc[c] -= P.MOD;
            }
        }
        R += P.s1; if (R >= P.MOD) R -= P.MOD;
    }
}

/* Variant 10+N: FLAT x N. N chains per thread, and each chain refills with the
 * thread's next residue the moment its word dies (own chain position t per
 * chain), so no chain waits for the slowest of its group or of its warp.
 * Per-prime constants come from shared memory since t is not uniform.
 * Uses three 17-bit limbs: needs MOD < 2^51 and fold = 2^17*(1+2*b2) < 2^32. */
template <int NCH>
__global__ void __launch_bounds__(256)
sieve_flatN(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
            Params P, u32 *cnt, u64 *hits)
{
    __shared__ uint4 s_c[MAXTC_SH];     /* {rp | k1<<16, k2, mag, off}, k = 2^17, 2^34 mod r */
    for (u32 t = threadIdx.x; t < P.ntc; t += blockDim.x) {
        u32 r = c_rp[t];
        u32 k1 = (u32)((1ULL << 17) % r), k2 = (u32)((1ULL << 34) % r);
        s_c[t] = make_uint4(r | (k1 << 16), k2, c_mag[t], c_off[t]);
    }
    __syncthreads();
    u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= P.nthreads) return;
    /* residue dispenser */
    u64 R = Rpre[gid], Rn = R;
    u32 i1 = 0, i2 = 0;
    bool more = true;
    u64 Rc[NCH], sito[NCH];
    u32 tc[NCH];
    int nact = 0;
    #pragma unroll
    for (int c = 0; c < NCH; c++) {
        Rc[c] = Rn; tc[c] = 0;
        if (more) {
            sito[c] = ~0ULL; nact++;
            if (++i2 < P.c2) { Rn += P.s2; if (Rn >= P.MOD) Rn -= P.MOD; }
            else { i2 = 0;
                if (++i1 < P.c1) { R += P.s1; if (R >= P.MOD) R -= P.MOD; Rn = R; }
                else more = false; }
        } else sito[c] = 0ULL;
    }
    const u32 ntc = P.ntc;
    while (nact > 0) {
        #pragma unroll
        for (int c = 0; c < NCH; c++) {
            if (sito[c] != 0ULL) {
                uint4 k = s_c[tc[c]];
                u32 p = k.x & 0xFFFF;
                u64 rr = Rc[c];
                u32 x = (u32)(rr & 0x1FFFF) + (u32)((rr >> 17) & 0x1FFFF) * (k.x >> 16)
                      + (u32)(rr >> 34) * k.y;
                u32 r = x - __umulhi(x, k.z) * p;
                if (r >= p) r -= p;
                if (r >= p) r -= p;
                u64 w = sito[c] & words[k.w + r];
                u32 t = tc[c] + 1;
                if (w == 0ULL || t == ntc) {
                    if (w != 0ULL) {
                        u32 sidx = atomicAdd(cnt, 1u);
                        if (sidx < P.cap) { hits[2 * sidx] = rr; hits[2 * sidx + 1] = w; }
                    }
                    if (more) {
                        Rc[c] = Rn; w = ~0ULL; t = 0;
                        if (++i2 < P.c2) { Rn += P.s2; if (Rn >= P.MOD) Rn -= P.MOD; }
                        else { i2 = 0;
                            if (++i1 < P.c1) { R += P.s1; if (R >= P.MOD) R -= P.MOD; Rn = R; }
                            else more = false; }
                    } else { w = 0ULL; nact--; }
                }
                sito[c] = w; tc[c] = t;
            }
        }
    }
}

/* Variant 20: STRIDED groups. The innermost loop is walked WITHOUT reducing
 * mod MOD: A = R + i1*s1 + i2*s2 as a plain integer (< ~105*MOD < 2^57). An
 * unreduced A = Rred + w*MOD stands for the candidates a = A + (b+64*shift)*MOD,
 * i.e. the same residue class with its 64-candidate window moved up by w, so
 * every candidate is still a legitimate, distinct one; a unit (K, shift) then
 * covers, per class, the multiples [64*shift + w, 64*shift + 64 + w), w < c1+c2.
 *
 * The gain: A mod r now advances by exactly s2 mod r per i2 step, so with each
 * table stored pre-rotated, W'[y] = W[(y * s2) mod r], NCH consecutive i2
 * values read NCH ADJACENT words: one index computation and ~one cache line
 * per prime for the whole group, instead of one of each per residue. The
 * rotation index y0 = (A mod r) * (s2 mod r)^-1 comes straight out of the limb
 * fold by pre-multiplying the fold constants with the inverse. Tables carry
 * NCH-1 wrap-around words so y0 + j needs no reduction. */
__constant__ u32 c_k0[MAXTC];
template <int NCH, int U>
__global__ void __launch_bounds__(256)
sieve_strided(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
              Params P, u32 *cnt, u64 *hits)
{
    u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= P.nthreads) return;
    u64 A1 = Rpre[gid];
    const u64 sG = P.s2 * NCH;
    for (u32 i1 = 0; i1 < P.c1; ++i1) {
        u64 Ag = A1;
        for (u32 i2 = 0; i2 < P.c2; i2 += NCH) {
            u32 l0 = (u32)(Ag & 0xFFFF), l1 = (u32)((Ag >> 16) & 0xFFFF);
            u32 l2 = (u32)((Ag >> 32) & 0xFFFF), l3 = (u32)(Ag >> 48);
            u64 sito[NCH];
            #pragma unroll
            for (int c = 0; c < NCH; c++) sito[c] = (i2 + c < P.c2) ? ~0ULL : 0ULL;
            /* U primes per exit test: the loads of U table rows are issued back
             * to back and the branch waits once per U rows instead of once per
             * row. ntc is padded by the host to a multiple of U with all-ones
             * dummy rows (prime 3, never kills). */
            for (u32 t = 0; t < P.ntc; t += U) {
                const u64 *w[U];
                #pragma unroll
                for (int u = 0; u < U; u++) {
                    u32 p = c_rp[t + u];
                    u32 x = l0 * c_k0[t + u] + l1 * c_k1[t + u] + l2 * c_k2[t + u] + l3 * c_k3[t + u];
                    u32 y = x - __umulhi(x, c_mag[t + u]) * p;
                    if (y >= p) y -= p;
                    if (y >= p) y -= p;
                    w[u] = words + c_off[t + u] + y;
                }
                u64 any = 0;
                #pragma unroll
                for (int c = 0; c < NCH; c++) {
                    u64 v = w[0][c];
                    #pragma unroll
                    for (int u = 1; u < U; u++) v &= w[u][c];
                    sito[c] &= v; any |= sito[c];
                }
                if (any == 0ULL) break;
            }
            #pragma unroll
            for (int c = 0; c < NCH; c++) {
                if (sito[c] != 0ULL) {
                    u32 s = atomicAdd(cnt, 1u);
                    if (s < P.cap) { hits[2 * s] = Ag + (u64)c * P.s2; hits[2 * s + 1] = sito[c]; }
                }
            }
            Ag += sG;
        }
        A1 += P.s1;
    }
}

/* Variant 21: FLAT STRIDED. Same strided groups as variant 20, but each lane
 * moves on to its next group the moment its own group dies, instead of the
 * whole warp fetching rows until its slowest lane is done (about 2x the rows
 * an average lane needs). The row index t is then lane-specific, so the
 * per-row constants cannot come from the constant cache: they are packed
 * into 32-byte Row records in device memory, U of them (one 128-byte cache
 * line for U = 4) per batch. The batch itself is branch-free, so its
 * U * NCH word loads are all in flight together. */
struct Row { u32 rp, mag, k0, k1, k2, k3, off, pad; };
template <int NCH, int U>
__global__ void __launch_bounds__(256)
sieve_fs(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
         const Row *__restrict__ rows, Params P, u32 *cnt, u64 *hits)
{
    u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= P.nthreads) return;
    u64 A1 = Rpre[gid], Ag = A1;
    const u64 sG = P.s2 * NCH;
    u32 i1 = 0, i2 = 0, t = 0;
    u32 l0 = (u32)(Ag & 0xFFFF), l1 = (u32)((Ag >> 16) & 0xFFFF);
    u32 l2 = (u32)((Ag >> 32) & 0xFFFF), l3 = (u32)(Ag >> 48);
    u64 sito[NCH];
    #pragma unroll
    for (int c = 0; c < NCH; c++) sito[c] = ((u32)c < P.c2) ? ~0ULL : 0ULL;
    const u32 ntc = P.ntc;
    for (;;) {
        const u64 *w[U];
        #pragma unroll
        for (int u = 0; u < U; u++) {
            const Row rw = rows[t + u];
            u32 x = l0 * rw.k0 + l1 * rw.k1 + l2 * rw.k2 + l3 * rw.k3;
            u32 y = x - __umulhi(x, rw.mag) * rw.rp;
            if (y >= rw.rp) y -= rw.rp;
            if (y >= rw.rp) y -= rw.rp;
            w[u] = words + rw.off + y;
        }
        u64 any = 0;
        #pragma unroll
        for (int c = 0; c < NCH; c++) {
            u64 v = w[0][c];
            #pragma unroll
            for (int u = 1; u < U; u++) v &= w[u][c];
            sito[c] &= v; any |= sito[c];
        }
        t += U;
        if (any == 0ULL || t >= ntc) {
            if (any != 0ULL) {
                #pragma unroll
                for (int c = 0; c < NCH; c++) {
                    if (sito[c] != 0ULL) {
                        u32 sidx = atomicAdd(cnt, 1u);
                        if (sidx < P.cap) { hits[2 * sidx] = Ag + (u64)c * P.s2; hits[2 * sidx + 1] = sito[c]; }
                    }
                }
            }
            i2 += NCH;
            if (i2 < P.c2) Ag += sG;
            else {
                i2 = 0;
                if (++i1 >= P.c1) break;
                A1 += P.s1; Ag = A1;
            }
            t = 0;
            l0 = (u32)(Ag & 0xFFFF); l1 = (u32)((Ag >> 16) & 0xFFFF);
            l2 = (u32)((Ag >> 32) & 0xFFFF); l3 = (u32)(Ag >> 48);
            #pragma unroll
            for (int c = 0; c < NCH; c++) sito[c] = (i2 + c < P.c2) ? ~0ULL : 0ULL;
        }
    }
}

#include "kernels_cmp.cuh"

/* DIAGNOSTIC (KFIX env): every lane runs exactly P.ntc rows per group, no exit test,
 * no hit output; measures cost per lane-row at full lane utilisation. */
template <int NCH, int U>
__global__ void __launch_bounds__(256)
sieve_diag(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
              Params P, u32 *cnt, u64 *hits)
{
    u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= P.nthreads) return;
    u64 A1 = Rpre[gid];
    const u64 sG = P.s2 * NCH;
    u64 acc = 0;
    for (u32 i1 = 0; i1 < P.c1; ++i1) {
        u64 Ag = A1;
        for (u32 i2 = 0; i2 < P.c2; i2 += NCH) {
            u32 l0 = (u32)(Ag & 0xFFFF), l1 = (u32)((Ag >> 16) & 0xFFFF);
            u32 l2 = (u32)((Ag >> 32) & 0xFFFF), l3 = (u32)(Ag >> 48);
            u64 sito[NCH];
            #pragma unroll
            for (int c = 0; c < NCH; c++) sito[c] = ~0ULL;
            for (u32 t = 0; t < P.ntc; t += U) {
                const u64 *w[U];
                #pragma unroll
                for (int u = 0; u < U; u++) {
                    u32 p = c_rp[t + u];
                    u32 x = l0 * c_k0[t + u] + l1 * c_k1[t + u] + l2 * c_k2[t + u] + l3 * c_k3[t + u];
                    u32 y = x - __umulhi(x, c_mag[t + u]) * p;
                    if (y >= p) y -= p;
                    if (y >= p) y -= p;
                    w[u] = words + c_off[t + u] + (y & P.cap);
                }
                #pragma unroll
                for (int c = 0; c < NCH; c++) {
                    u64 v = w[0][c];
                    #pragma unroll
                    for (int u = 1; u < U; u++) v &= w[u][c];
                    sito[c] ^= v;
                }
            }
            #pragma unroll
            for (int c = 0; c < NCH; c++) acc += sito[c];
            Ag += sG;
        }
        A1 += P.s1;
    }
    if (acc == 0x123456789ULL) atomicAdd(cnt, 1u);
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
    const char *out = NULL, *units_file = NULL;
    u64 slice_i = 0, slice_n = 1;    /* --slice i n: take plan entries with index % n == i */
    bool verify_mode = false, resume = false, hist_mode = false;
    int tqdepth = 64;                /* stage-3 queue depth (units waiting for exact tests) */
    int nprep = 12;                  /* host-prep threads (units prepared ahead of the GPU) */
    int nomp = 0, kernel = 20;       /* strided groups; see the kernel comments for the measured ladder */
    u64 shbytes = 65536;
    int unr = 4;                     /* primes per exit test in the strided kernel */
    int ct0 = 28, tv = 0, cmp = 1;   /* --kernel 30/31 */
    int nch = 7;                     /* group size for the strided kernel (--kernel 20) */

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
        else if (!strcmp(k, "--kernel")) kernel = atoi(NEXT());
        else if (!strcmp(k, "--prep")) nprep = atoi(NEXT());
        else if (!strcmp(k, "--tq")) tqdepth = atoi(NEXT());
        else if (!strcmp(k, "--units")) units_file = NEXT();
        else if (!strcmp(k, "--slice")) { slice_i = strtoull(NEXT(), NULL, 10); slice_n = strtoull(NEXT(), NULL, 10); if (!slice_n || slice_i >= slice_n) { fprintf(stderr, "bad --slice\n"); return 2; } }
        else if (!strcmp(k, "--nch")) nch = atoi(NEXT());
        else if (!strcmp(k, "--unr")) unr = atoi(NEXT());
        else if (!strcmp(k, "--t0")) ct0 = atoi(NEXT());
        else if (!strcmp(k, "--tv")) tv = atoi(NEXT());
        else if (!strcmp(k, "--cmp")) cmp = atoi(NEXT());
        else if (!strcmp(k, "--shbytes")) shbytes = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(k, "--verify")) verify_mode = true;
        else if (!strcmp(k, "--resume")) resume = true;
        else if (!strcmp(k, "--hist")) hist_mode = true;
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
    {
        int maxsh = 0; CK(cudaDeviceGetAttribute(&maxsh, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev));
        if (shbytes > (u64)maxsh) shbytes = (u64)maxsh;
        CK(cudaFuncSetAttribute(sieve_sh, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shbytes));
        fprintf(stderr, "shared memory per block: %d max, using %llu\n", maxsh, (unsigned long long)shbytes);
    }
    cudaStream_t stream[2]; CK(cudaStreamCreate(&stream[0])); CK(cudaStreamCreate(&stream[1]));
    const u32 CAP = 1u << 23;        /* survivor words per unit (128 MB per buffer set) */
    u64 *d_R[2] = {NULL, NULL}, *d_W[2] = {NULL, NULL}, *d_H[2];
    Row *d_C[2] = {NULL, NULL};
    std::vector<Row> h_rows(MAXTC);
    u32 *d_cnt[2];
    size_t cap_R[2] = {0, 0}, cap_W[2] = {0, 0};
    for (int b = 0; b < 2; b++) {
        CK(cudaMalloc(&d_H[b], (size_t)CAP * 16));
        CK(cudaMalloc(&d_cnt[b], 4));
        CK(cudaMalloc(&d_C[b], sizeof(Row) * MAXTC));
    }
    std::vector<u64> h_hits((size_t)CAP * 2);
    Unit inflight[2]; inflight[0].valid = inflight[1].valid = false;
    for (int b = 0; b < 2; b++) CK(cudaEventCreate(&inflight[b].done));

    FILE *of = out ? fopen(out, "a") : NULL;
    if (out && !of) { perror("open --out"); return 2; }

    int global_best = 0;
    const bool kbench = getenv("KBENCH") != NULL; double kb_total = 0;
    const u32 TPB = getenv("KTPB") ? (u32)atoi(getenv("KTPB")) : 256u;
    double t0 = lc_now_s(), covered = 0;
    u64 units = 0, total_conf = 0, total_surv = 0;

    /* ---- host prep is done ahead of the GPU by a pool of threads (one unit per
     * thread, serial inside), delivered to the launch loop in plan order ---- */
    /* page-locked, grow-only host buffer: uploads from it are true async DMA copies
     * (pageable uploads from the prep threads' buffers were measured 10-40x slower) */
    struct PinBuf { u64 *p = NULL; size_t n = 0, cap = 0;
        void resize(size_t m) {
            if (m > cap) { if (p) cudaFreeHost(p);
                if (cudaHostAlloc((void **)&p, m * 8, cudaHostAllocDefault) != cudaSuccess) {
                    fprintf(stderr, "FATAL: cudaHostAlloc %zu bytes\n", m * 8); fflush(NULL); _exit(2); }
                cap = m; }
            n = m; }
        u64 *data() { return p; }
        u64 &operator[](size_t i) { return p[i]; } };
    struct Prep { int status; u64 K, shift, d, MOD, s1, s2, nthreads, total, tot_words;
                  u32 c1, c2; int ntc; double ksetup_s, prep_s;
                  PinBuf Rpre, words;
                  std::vector<u32> offs, rps, mags, k1s, k2s, k3s, k0s; u64 tcp5[5]; };
    double prep_s = 0, up_s = 0, wait_s = 0, ksetup_s = 0, fin_s = 0, launch_s = 0;

    /* ---- finish a unit: wait for its kernel, pull survivors, stage 3 ---- */
    /* ---- stage 3 runs on its own worker thread, fed through a small queue, so
     * the GPU never waits for the CPU: a large unit's exact tests can take tens
     * of seconds while the next kernels are short ---- */
    struct Task { Unit U; std::vector<u64> hv; u32 nsurv; bool overflow;
                  double gpu_s, prep_s, up_s, wait_s; };
    std::deque<Task> tq;
    std::mutex tq_m;
    std::condition_variable tq_cv, tq_space;
    bool tq_done = false;
    auto stage3 = [&](Task &T) {
        Unit &U = T.U;
        const std::vector<u64> &h_hits = T.hv;
        const u32 nsurv = T.nsurv; const bool overflow = T.overflow;
        const double gpu_s = T.gpu_s, prep_s = T.prep_s, up_s = T.up_s, wait_s = T.wait_s;
        /* ---- stage 3 on the CPU: exact test of each surviving window ---- */
        double ts3 = lc_now_s();
        u64 nconf = 0, nbits = 0;
        u64 boff = U.shift * 64, d = U.d, MOD = U.MOD;
        std::vector<Hit> hits;
        u64 hist[64]; memset(hist, 0, sizeof hist);   /* best run per window, after extension */
        #pragma omp parallel for schedule(dynamic, 4) reduction(+:nconf, nbits) reduction(+:hist[:64])
        for (u32 s = 0; s < nsurv; s++) {
            u64 Rb = h_hits[2 * s], sito = h_hits[2 * s + 1];
            while (sito) {
                int bb = __builtin_ctzll(sito); sito &= sito - 1;
                nbits++;
                u64 a0 = Rb + (boff + (u64)bb) * MOD;
                u64 best_run = 0, best_start = 0, cur = 0, cur_start = 0;
                const u64 want = hist_mode ? 1 : ((report > 1) ? (u64)report : 1);
                for (u64 kk = 0; kk < nterms; kk++) {
                    /* stop as soon as even a perfect tail cannot reach the
                     * reporting threshold (full scan only with --hist) */
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
                hist[run < 63 ? run : 63]++;
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
                "surv=%u bits=%llu conf=%llu best=%d gpu=%.2fs cpu=%.2fs prep=%.2fs up=%.3fs wait=%.3fs (%.3g res/s)\n",
                (unsigned long long)U.K, (unsigned long long)U.shift,
                (unsigned long long)MOD, (unsigned long long)U.nthreads,
                (double)U.total, U.ntc, nsurv, (unsigned long long)nbits,
                (unsigned long long)nconf, global_best, gpu_s, cpu_s, prep_s, up_s, wait_s,
                (double)U.total / gpu_s);
        if (of) {
            fprintf(of, "{\"K\":%llu,\"shift\":%llu,\"res\":%.6g,\"surv\":%u,"
                    "\"conf\":%llu,\"best\":%d,\"gpu_s\":%.2f,\"cpu_s\":%.2f,"
                    "\"covered\":%.6g,\"overflow\":%s}\n",
                    (unsigned long long)U.K, (unsigned long long)U.shift,
                    (double)U.total, nsurv, (unsigned long long)nconf, global_best,
                    gpu_s, cpu_s, covered, overflow ? "true" : "false");
            if (hist_mode) {
            fprintf(of, "{\"hist\":true,\"K\":%llu,\"shift\":%llu,\"runs\":[",
                    (unsigned long long)U.K, (unsigned long long)U.shift);
            for (int i = 0; i < 64; i++) fprintf(of, "%s%llu", i ? "," : "", (unsigned long long)hist[i]);
            fprintf(of, "]}\n");
            }
            fflush(of);
        }
        if (verify_mode)
            printf("VERIFY K=%llu shift=%llu residues=%.6g survivors=%u bits=%llu "
                   "confirmed=%llu\n", (unsigned long long)U.K,
                   (unsigned long long)U.shift, (double)U.total, nsurv,
                   (unsigned long long)nbits, (unsigned long long)nconf);
    };
    std::thread worker([&]() {
        for (;;) {
            Task T;
            {
                std::unique_lock<std::mutex> lk(tq_m);
                tq_cv.wait(lk, [&] { return tq_done || !tq.empty(); });
                if (tq.empty()) return;
                T = std::move(tq.front()); tq.pop_front();
            }
            tq_space.notify_all();
            stage3(T);
        }
    });

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
        if (nsurv) {   /* on this unit's stream: a plain cudaMemcpy uses the legacy default
                        * stream, which waits for the OTHER stream's running kernel */
            CK(cudaMemcpyAsync(h_hits.data(), d_H[b], (size_t)nsurv * 16, cudaMemcpyDeviceToHost, stream[b]));
            CK(cudaStreamSynchronize(stream[b]));
        }

        Task T;
        T.U = U; T.hv.assign(h_hits.begin(), h_hits.begin() + (size_t)nsurv * 2);
        T.nsurv = nsurv; T.overflow = overflow; T.gpu_s = gpu_s;
        T.prep_s = prep_s; T.up_s = up_s; T.wait_s = wait_s;
        {
            std::unique_lock<std::mutex> lk(tq_m);
            tq_space.wait(lk, [&] { return tq.size() < (size_t)tqdepth; });
            tq.push_back(std::move(T));
        }
        tq_cv.notify_one();
        U.valid = false;
    };

    /* work list: either the K range x shifts, or an explicit "K shift" file
     * (tools/plan_units.py) consumed in order */
    std::vector<std::pair<u64, u64> > ulist;
    if (units_file) {
        FILE *uf = fopen(units_file, "r");
        if (!uf) { perror("open --units"); return 2; }
        unsigned long long a_, b_;
        while (fscanf(uf, "%llu %llu", &a_, &b_) == 2) ulist.push_back(std::make_pair((u64)a_, (u64)b_));
        fclose(uf);
        fprintf(stderr, "units: %zu from %s\n", ulist.size(), units_file);
    } else {
        for (u64 K = kmin; K <= kmax; K++)
            for (u64 sh = shift0; sh < shift0 + shifts; sh++) ulist.push_back(std::make_pair(K, sh));
    }
    auto prepare = [&](const u64 K, const u64 shift, Prep &PP) -> int {
        PinBuf &Rpre = PP.Rpre, &words = PP.words;
        std::vector<u32> &offs = PP.offs, &rps = PP.rps, &mags = PP.mags, &k1s = PP.k1s,
                         &k2s = PP.k2s, &k3s = PP.k3s, &k0s = PP.k0s;
        std::vector<u64> tcpv(MAXTC); u64 *tcp = tcpv.data();
        double t_ks = lc_now_s();
        u64 d = K * D0;
        if (d / D0 != K) { fprintf(stderr, "K overflow at %llu\n", (unsigned long long)K); return 3; }

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
        if (ntb < 4) { fprintf(stderr, "K=%llu: too few tier-B primes\n", (unsigned long long)K); return 1; }

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
        if (in1 < 0 || in2 < 0) { fprintf(stderr, "K=%llu: no inner loops\n", (unsigned long long)K); return 1; }

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
            /* strided kernel: a prime with s2 = 0 mod r has no rotation; drop it
             * from the GPU filter for this K (stage 3 is exact, so nothing is lost) */
            if (kernel >= 20 && sstep[in2] % r == 0) continue;
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
        const int ntc_real = ntc;
        if (kernel >= 20) while (ntc % unr) tcp[ntc++] = 3;      /* all-ones dummy rows */
        const u64 pad = (kernel >= 20) ? (u64)(nch - 1) : 0;
        if (kernel >= 20) {
            u64 amax = MOD + (u64)C[in1].c * sstep[in1] + ((u64)C[in2].c + (u64)nch) * sstep[in2];
            (void)amax;
            /* unreduced walk: A < MOD * (1 + c1 + c2 + nch) must fit in 64 bits */
            if (b2 > 16383 || MOD >= (1ULL << 56) || C[in1].c + C[in2].c + (u64)nch >= 255) {
                fprintf(stderr, "K=%llu: strided kernel needs b2 <= 16383, MOD < 2^56; skipping\n",
                        (unsigned long long)K); return 1; }
        }
        u64 tot_words = 0;
        for (int i = 0; i < ntc; i++) tot_words += tcp[i] + pad;
        offs.resize(ntc); rps.resize(ntc); mags.resize(ntc);
        k1s.resize(ntc); k2s.resize(ntc); k3s.resize(ntc); k0s.resize(ntc);
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
                k0s[t] = 1;
                if (kernel >= 20 && t >= ntc_real) { k0s[t] = k1s[t] = k2s[t] = k3s[t] = 0; }
                else if (kernel >= 20) {     /* pre-multiply the fold by (s2 mod r)^-1 */
                    u64 inv = lc_inv_mod(sstep[in2] % r, r);
                    k0s[t] = (u32)inv;
                    k1s[t] = (u32)((u64)k1s[t] * inv % r);
                    k2s[t] = (u32)((u64)k2s[t] * inv % r);
                    k3s[t] = (u32)((u64)k3s[t] * inv % r);
                }
                acc += r + pad;
            }
        }
        {
            u64 boff = shift * 64;
            double t_prep = lc_now_s();
            PP.ksetup_s = t_prep - t_ks;
            std::vector<u64> tmp;
            for (int t = 0; t < ntc; t++) {
                u64 r = tcp[t], acc = offs[t];
                if (t >= ntc_real) { for (u64 y = 0; y < r + pad; y++) words[acc + y] = ~0ULL; continue; }
                /* bit b of W[x] is ok[(x + b0 + b*modr) mod r]; ok has only a few zeros
                 * (the nterms forbidden residues, or just 0 when r | D0), so start from
                 * all-ones and clear, for each zero z and bit b, x = z - b0 - b*modr:
                 * O(zeros * 64) per prime instead of O(r * 64), same table */
                u64 modr = MOD % r;
                u64 b0 = lc_mulmod(boff % r, modr, r);
                tmp.resize(r);
                u64 *dst = &words[acc];
                if (kernel >= 20) dst = tmp.data();
                for (u64 x = 0; x < r; x++) dst[x] = ~0ULL;
                {
                    const u64 nz = (D0 % r == 0) ? 1 : nterms, dm = d % r;
                    u64 v = 0;
                    for (u64 kk = 0; kk < nz; kk++) {
                        u64 x = (v + r - b0) % r;
                        for (int b = 0; b < 64; b++) {
                            dst[x] &= ~((u64)1 << b);
                            x = (x >= modr) ? x - modr : x + r - modr;
                        }
                        v = (v + r - dm) % r;
                    }
                }
                if (kernel >= 20) {          /* W'[y] = W[(y * s2) mod r], plus wrap pad */
                    u64 s2r = sstep[in2] % r, x = 0;
                    for (u64 y = 0; y < r + pad; y++) {
                        words[acc + y] = tmp[x];
                        x += s2r; if (x >= r) x -= r;
                    }
                }
            }

            PP.prep_s = lc_now_s() - t_prep;
        }
        PP.K = K; PP.shift = shift; PP.d = d; PP.MOD = MOD; PP.s1 = sstep[in1]; PP.s2 = sstep[in2];
        PP.c1 = (u32)C[in1].c; PP.c2 = (u32)C[in2].c; PP.nthreads = nthreads; PP.total = total;
        PP.ntc = ntc; PP.tot_words = tot_words;
        for (int i = 0; i < 5; i++) PP.tcp5[i] = tcp[i];
        return 0;
    };
    std::vector<std::pair<u64, u64> > jobs;
    for (size_t ui = 0; ui < ulist.size(); ui++) {
        if (ui % slice_n != slice_i) continue;
        if (!done_units.empty() && done_units.count(unit_key(ulist[ui].first, ulist[ui].second))) continue;
        jobs.push_back(ulist[ui]);
    }
    std::vector<Prep *> slot(jobs.size(), (Prep *)NULL);
    std::mutex pq_m; std::condition_variable pq_cv;
    size_t pq_next = 0, pq_consumed = 0; bool pq_quit = false;
    const size_t pq_window = (size_t)nprep * 2;
    std::vector<Prep *> pq_free;
    std::vector<std::thread> preppers;
    for (int w = 0; w < nprep; w++) preppers.emplace_back([&]() {
        for (;;) {
            size_t j;
            {
                std::unique_lock<std::mutex> lk(pq_m);
                pq_cv.wait(lk, [&] { return pq_quit || (pq_next < jobs.size() && pq_next < pq_consumed + pq_window); });
                if (pq_quit) return;
                j = pq_next++;
            }
            Prep *pp = NULL;
            {
                std::lock_guard<std::mutex> lk(pq_m);
                if (!pq_free.empty()) { pp = pq_free.back(); pq_free.pop_back(); }
            }
            if (!pp) pp = new Prep();     /* recycled: the big vectors keep their capacity */
            pp->status = prepare(jobs[j].first, jobs[j].second, *pp);
            { std::lock_guard<std::mutex> lk(pq_m); slot[j] = pp; }
            pq_cv.notify_all();
        }
    });
    auto prep_shutdown = [&]() {
        { std::lock_guard<std::mutex> lk(pq_m); pq_quit = true; }
        pq_cv.notify_all();
        for (auto &t : preppers) t.join();
    };
    bool first_unit = true;
    int buf = 0;
    for (size_t ji = 0; ji < jobs.size() && !stop_requested; ji++) {
        Prep *ppp;
        {
            std::unique_lock<std::mutex> lk(pq_m);
            pq_cv.wait(lk, [&] { return slot[ji] != NULL; });
            ppp = slot[ji]; slot[ji] = NULL; pq_consumed = ji + 1;
        }
        pq_cv.notify_all();
        auto pp_recycle = [&](Prep *q) { std::lock_guard<std::mutex> lk(pq_m); pq_free.push_back(q); };
        std::unique_ptr<Prep, decltype(pp_recycle)> pp_owner(ppp, pp_recycle);
        Prep &PP = *ppp;
        if (PP.status == 1) continue;
        if (PP.status == 3) break;
        if (PP.status == 2) { fflush(NULL); _exit(2); }
        const u64 K = PP.K, shift = PP.shift, d = PP.d, MOD = PP.MOD, nthreads = PP.nthreads,
                  total = PP.total, tot_words = PP.tot_words;
        const int ntc = PP.ntc;
        PinBuf &Rpre = PP.Rpre, &words = PP.words;
        std::vector<u32> &offs = PP.offs, &rps = PP.rps, &mags = PP.mags, &k1s = PP.k1s,
                         &k2s = PP.k2s, &k3s = PP.k3s, &k0s = PP.k0s;
        const u64 *tcp = PP.tcp5;
        ksetup_s = PP.ksetup_s; prep_s = PP.prep_s;
        {
            {
            if (first_unit) {
                first_unit = false;
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
                        if (rr != x % p) { fprintf(stderr, "FATAL magic mod p=%u x=%u\n", p, x); fflush(NULL); _exit(2); }
                    }
                }
                fprintf(stderr, "tierC order: %llu %llu %llu %llu %llu ... (%d primes, "
                        "%llu words) magic-mod OK\n", (unsigned long long)tcp[0],
                        (unsigned long long)tcp[1], (unsigned long long)tcp[2],
                        (unsigned long long)tcp[3], (unsigned long long)tcp[4], ntc,
                        (unsigned long long)tot_words);
            }
            }
        }
        {
            /* ---- buffer set `buf` is free once the unit before last is finished ---- */
            double t_fin = lc_now_s();
            finish_unit(buf);
            fin_s = lc_now_s() - t_fin;
            if (cap_R[buf] < nthreads) {
                if (d_R[buf]) CK(cudaFree(d_R[buf]));
                CK(cudaMalloc(&d_R[buf], nthreads * 8)); cap_R[buf] = nthreads;
            }
            if (cap_W[buf] < tot_words) {
                if (d_W[buf]) CK(cudaFree(d_W[buf]));
                CK(cudaMalloc(&d_W[buf], tot_words * 8)); cap_W[buf] = tot_words;
            }
            /* uploads go first, on this buffer set's own stream, so they overlap
             * the other stream's kernel; only the shared constants need it done */
            double t_up = lc_now_s();
            CK(cudaMemcpyAsync(d_R[buf], Rpre.data(), nthreads * 8, cudaMemcpyHostToDevice, stream[buf]));
            CK(cudaMemcpyAsync(d_W[buf], words.data(), tot_words * 8, cudaMemcpyHostToDevice, stream[buf]));
            CK(cudaMemsetAsync(d_cnt[buf], 0, 4, stream[buf]));
            if (kernel == 21) {
                for (int t = 0; t < ntc; t++) {
                    Row r_; r_.rp = rps[t]; r_.mag = mags[t]; r_.k0 = k0s[t]; r_.k1 = k1s[t];
                    r_.k2 = k2s[t]; r_.k3 = k3s[t]; r_.off = offs[t]; r_.pad = 0; h_rows[t] = r_;
                }
                CK(cudaMemcpyAsync(d_C[buf], h_rows.data(), sizeof(Row) * ntc, cudaMemcpyHostToDevice, stream[buf]));
            }
            CK(cudaStreamSynchronize(stream[buf]));
            up_s = lc_now_s() - t_up;
            double t_wait = lc_now_s();
            if (inflight[buf ^ 1].valid) CK(cudaEventSynchronize(inflight[buf ^ 1].done));
            wait_s = lc_now_s() - t_wait;
            CK(cudaMemcpyToSymbol(c_off, offs.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_rp, rps.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_mag, mags.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_k1, k1s.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_k2, k2s.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_k3, k3s.data(), ntc * 4));
            CK(cudaMemcpyToSymbol(c_k0, k0s.data(), ntc * 4));

            Params P;
            P.MOD = MOD; P.s1 = PP.s1; P.s2 = PP.s2;
            P.c1 = PP.c1; P.c2 = PP.c2; P.ntc = (u32)ntc;
            P.cap = CAP; P.nthreads = (u32)nthreads;
            Unit &U = inflight[buf];
            U.K = K; U.shift = shift; U.d = d; U.MOD = MOD; U.total = total;
            U.nthreads = nthreads; U.ntc = ntc; U.valid = true;
            U.t_launch = lc_now_s();
            u32 blocks = (u32)((nthreads + TPB - 1) / TPB);
            if (kernel == 1) {
                if (ntc > MAXTC_SH) { fprintf(stderr, "FATAL: ntc > MAXTC_SH\n"); return 2; }
                sieve_flat<<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
            } else if (kernel == 21) {
                #define LF2(N, UU) sieve_fs<N, UU><<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], d_C[buf], P, d_cnt[buf], d_H[buf])
                switch (nch * 10 + unr) {
                    case 74: LF2(7, 4); break; case 54: LF2(5, 4); break; case 44: LF2(4, 4); break;
                    case 64: LF2(6, 4); break; case 84: LF2(8, 4); break; case 72: LF2(7, 2); break;
                    case 42: LF2(4, 2); break; case 34: LF2(3, 4); break; case 78: LF2(7, 8); break;
                    case 48: LF2(4, 8); break; case 114: LF2(11, 4); break; case 24: LF2(2, 4); break;
                    case 32: LF2(3, 2); break; case 22: LF2(2, 2); break; case 52: LF2(5, 2); break;
                    case 82: LF2(8, 2); break; case 112: LF2(11, 2); break; case 142: LF2(14, 2); break;
                    case 144: LF2(14, 4); break; case 71: LF2(7, 1); break; case 41: LF2(4, 1); break;
                    default: fprintf(stderr, "bad --nch/--unr\n"); return 2;
                }
                #undef LF2
            } else if (kernel == 30 || kernel == 31) {
                if (ntc < ct0 || (ntc % 4) || (ct0 % 4)) { fprintf(stderr, "FATAL: kernel 30/31 needs ntc >= t0, multiples of 4\n"); return 2; }
                if (tv != 0) { fprintf(stderr, "FATAL: --tv (dual-copy rows) not ported; it was slower\n"); return 2; }
                if (TPB != 256) { fprintf(stderr, "FATAL: kernel 30/31 needs 256 threads/block\n"); return 2; }
                if (kernel == 30) {
                #define LC(T0_, TV_, CMP_) else if (nch == 7 && ct0 == T0_ && tv == TV_ && cmp == CMP_) \
                    sieve_cmp<7, T0_, TV_, CMP_><<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
                if (0) {} CMP_COMBOS(LC)
                else { fprintf(stderr, "--kernel 30: (t0, tv, cmp) = (%d, %d, %d) not instantiated\n", ct0, tv, cmp); return 2; }
                #undef LC
                } else {
                #define LW(N_, T0_) else if (nch == N_ && ct0 == T0_) \
                    sieve_cmpw<N_, T0_><<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
                if (0) {} CMPW_COMBOS(LW)
                else { fprintf(stderr, "--kernel 31: (nch, t0) = (%d, %d) not instantiated\n", nch, ct0); return 2; }
                #undef LW
                }
            } else if (kernel == 29) {
                P.ntc = (u32)atoi(getenv("KFIX")); P.ntc -= P.ntc % 4; P.cap = getenv("KYMASK") ? (u32)strtoul(getenv("KYMASK"), NULL, 0) : 0xFFFFFFFFu;
                if (getenv("KROW0")) { fprintf(stderr, "KROW0 unsupported\n"); return 2; }
                sieve_diag<7, 4><<<blocks, TPB, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
            } else if (kernel >= 20) {
                #define LS(N, UU) sieve_strided<N, UU><<<blocks, TPB, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf])
                switch (nch * 10 + unr) {
                    case 75: LS(7, 5); break; case 76: LS(7, 6); break; case 78: LS(7, 8); break;
                    case 64: LS(6, 4); break; case 54: LS(5, 4); break; case 94: LS(9, 4); break;
                    case 66: LS(6, 6); break; case 56: LS(5, 6); break; case 58: LS(5, 8); break; case 68: LS(6, 8); break;
                    case 41: LS(4, 1); break; case 81: LS(8, 1); break; case 111: LS(11, 1); break;
                    case 141: LS(14, 1); break; case 191: LS(19, 1); break; case 281: LS(28, 1); break;
                    case 82: LS(8, 2); break; case 112: LS(11, 2); break; case 142: LS(14, 2); break;
                    case 72: LS(7, 2); break; case 73: LS(7, 3); break; case 74: LS(7, 4); break;
                    case 83: LS(8, 3); break; case 113: LS(11, 3); break; case 143: LS(14, 3); break;
                    case 84: LS(8, 4); break; case 114: LS(11, 4); break; case 144: LS(14, 4); break;
                    case 86: LS(8, 6); break; case 116: LS(11, 6); break; case 118: LS(11, 8); break;
                    default: fprintf(stderr, "bad --nch/--unr\n"); return 2;
                }
                #undef LS
            } else if (kernel >= 10) {
                if (ntc > MAXTC_SH || MOD >= (1ULL << 51) || b2 > 16000) {
                    fprintf(stderr, "FATAL: flatN needs ntc<=%d, MOD<2^51, b2<=16000\n", MAXTC_SH); return 2; }
                #define LF(N) sieve_flatN<N><<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf])
                switch (kernel - 10) {
                    case 1: LF(1); break; case 2: LF(2); break; case 3: LF(3); break;
                    case 4: LF(4); break; case 5: LF(5); break; case 6: LF(6); break;
                    case 8: LF(8); break; case 10: LF(10); break; case 12: LF(12); break;
                    default: fprintf(stderr, "bad --kernel\n"); return 2;
                }
                #undef LF
            } else if (kernel == 4) {
                sieve_ilpN<4><<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
            } else if (kernel == 5) {
                sieve_ilpN<8><<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
            } else if (kernel == 6) {
                sieve_ilpN<3><<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
            } else if (kernel == 7) {
                sieve_ilpN<6><<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
            } else if (kernel == 3) {
                sieve_ilp2<<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
            } else if (kernel == 2) {
                u32 nsh = 0, nshw = 0;
                while (nsh < (u32)ntc && (nshw + rps[nsh]) * 8 <= shbytes) { nshw += rps[nsh]; nsh++; }
                sieve_sh<<<blocks, 256, nshw * 8, stream[buf]>>>(d_R[buf], d_W[buf], P, nsh, nshw, d_cnt[buf], d_H[buf]);
            } else
                sieve<<<blocks, 256, 0, stream[buf]>>>(d_R[buf], d_W[buf], P, d_cnt[buf], d_H[buf]);
            CK(cudaGetLastError());
            CK(cudaEventRecord(U.done, stream[buf]));
            if (kbench) { CK(cudaEventSynchronize(U.done));
                double ks = lc_now_s() - U.t_launch; kb_total += ks;
                fprintf(stderr, "ksolo K=%llu sh=%llu res=%.4g ntc=%d %.4fs %.3e res/s\n", (unsigned long long)K,
                        (unsigned long long)shift, (double)total * 1.0, ntc, ks, (double)total / ks); }
            launch_s = lc_now_s() - t_wait - wait_s;
            if (verify_mode) fprintf(stderr, "  timing K=%llu ksetup=%.3f prep=%.3f fin=%.3f up=%.3f wait=%.3f launch=%.3f\n",
                (unsigned long long)K, ksetup_s, prep_s, fin_s, up_s, wait_s, launch_s);
            buf ^= 1;
        }
    }
    prep_shutdown();
    finish_unit(buf);
    finish_unit(buf ^ 1);
    { std::lock_guard<std::mutex> lk(tq_m); tq_done = true; }
    tq_cv.notify_all();
    worker.join();

    if (of) fclose(of);
    if (kbench) fprintf(stderr, "KBENCH total kernel %.4fs\n", kb_total);
    printf("GPU BEST n=%d  units=%llu covered=%.4g raw a  survivors=%llu confirmed=%llu  %.1fs%s\n",
           global_best, (unsigned long long)units, covered,
           (unsigned long long)total_surv, (unsigned long long)total_conf,
           lc_now_s() - t0, stop_requested ? " INTERRUPTED" : "");
    return 0;
}

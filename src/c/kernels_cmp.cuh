/* Candidate kernels (research copy). Included into apsearch_cuda.cu right before
 * the host section (needs c_* constants, c_k0, Params, u64/u32).
 *
 * sieve_cmp<NCH, T0, TV, CMP>  (--kernel 30, --t0 T0 --tv TV --cmp 0|1)
 *
 * Phase 1: every lane runs rows [0, T0) of its next FRESH group, fully unrolled, so
 *   per-row constants are compile-time constant-bank operands (no LDC, no address
 *   arithmetic on t). Loads are still skipped per lane once its group is dead
 *   (checked every 4 rows), so L1 traffic stays what it is in sieve_strided.
 * TV: rows t < TV use the DUAL-COPY layout (row_ptr): 3 x LDG.128 + 1 x LDG.64 per
 *   row instead of 7 x LDG.64, i.e. 4 instead of 7 L1 sector reads per lane-row.
 *   Footprint of those rows doubles, so keep TV small enough that rows [0, ~32)
 *   still fit L1 (TV = 12: 39 KB + 61 KB for rows 12..31).
 * CMP = 1: after phase 1, surviving groups are packed (ballot + byte permutation in
 *   shared memory + 8 u64 shuffles) into a per-lane PENDING slot; when 32 are pending
 *   the warp runs rows [T0, ntc) on them (LDC + exit test every U = 4, like
 *   sieve_strided). All 32 lanes are at the same row, so constants stay
 *   warp-uniform while the lanes are dense: the warp no longer idles waiting for
 *   the max-of-32 chain (model: 60.7 -> 42 warp-rows per 32 groups at T0 = 28).
 * CMP = 0: no compaction; each lane continues its own group (A/B test of TV alone).
 * Host must pad ntc to a multiple of 4 and ntc >= T0; T0, TV multiples of 4.
 */

#define CMP_U 4
#ifndef CMP_GUARD0
#define CMP_GUARD0 12   /* no group dies before row ~12 (model): no load guard there */
#endif
#ifndef CMP_MINB
#define CMP_MINB 4
#endif

__device__ __forceinline__ const u64 *row_ptr(const u64 *__restrict__ words, u32 off,
                                              u32 y, u32 p, bool vec)
{   /* dual copy: copy 0 at even off (r+8 words), copy 1 at off + r + 8 (odd) */
    if (vec) return words + (off + y + (y & 1u) * (p + 8u));
    return words + (off + y);
}

/* live = false: the lane's group is dead, skip the loads (predicated, no branch) */
template <int NCH>
__device__ __forceinline__ void and_words(const u64 *w, bool vec, bool live, u64 (&s)[NCH])
{
    if (vec) {
        const ulonglong2 *w2 = reinterpret_cast<const ulonglong2 *>(w);
        #pragma unroll
        for (int c = 0; c + 1 < NCH; c += 2) {
            ulonglong2 q = make_ulonglong2(0, 0);
            if (live) q = __ldg(w2 + c / 2);
            s[c] &= q.x; s[c + 1] &= q.y; }
        if (NCH & 1) { u64 q = 0; if (live) q = __ldg(w + NCH - 1); s[NCH - 1] &= q; }
    } else {
        #pragma unroll
        for (int c = 0; c < NCH; c++) { u64 q = 0; if (live) q = __ldg(w + c); s[c] &= q; }
    }
}

template <int NCH>
__device__ __forceinline__ void row_and(const u64 *__restrict__ words, u32 t, bool vec, bool live,
        u32 l0, u32 l1, u32 l2, u32 l3, u64 (&s)[NCH])
{
    u32 p = c_rp[t];
    u32 x = l0 * c_k0[t] + l1 * c_k1[t] + l2 * c_k2[t] + l3 * c_k3[t];
    u32 y = x - __umulhi(x, c_mag[t]) * p;
    if (y >= p) y -= p;
    if (y >= p) y -= p;
    and_words<NCH>(row_ptr(words, c_off[t], y, p, vec), vec, live, s);
}

/* rows [t0, ntc), scalar layout, exit test every CMP_U rows, warp-uniform t */
template <int NCH>
__device__ __forceinline__ void tail_rows(const u64 *__restrict__ words, u32 t0, u32 ntc,
        u64 A, u64 (&s)[NCH])
{
    u32 l0 = (u32)(A & 0xFFFF), l1 = (u32)((A >> 16) & 0xFFFF);
    u32 l2 = (u32)((A >> 32) & 0xFFFF), l3 = (u32)(A >> 48);
    u64 any = 0;
    #pragma unroll
    for (int c = 0; c < NCH; c++) any |= s[c];
    if (any == 0ULL) return;
    for (u32 t = t0; t < ntc; t += CMP_U) {
        #pragma unroll
        for (int u = 0; u < CMP_U; u++) row_and<NCH>(words, t + u, false, true, l0, l1, l2, l3, s);
        any = 0;
        #pragma unroll
        for (int c = 0; c < NCH; c++) any |= s[c];
        if (any == 0ULL) break;
    }
}

template <int NCH>
__device__ __forceinline__ void emit(u64 A, const u64 (&s)[NCH], u64 s2, u32 cap,
                                     u32 *cnt, u64 *hits)
{
    #pragma unroll
    for (int c = 0; c < NCH; c++) {
        if (s[c] != 0ULL) {
            u32 k = atomicAdd(cnt, 1u);
            if (k < cap) { hits[2 * k] = A + (u64)c * s2; hits[2 * k + 1] = s[c]; }
        }
    }
}

template <int NCH, int T0, int TV, int CMP>
__global__ void __launch_bounds__(256, CMP_MINB)   /* 4: 64 registers, same 32 warps/SM as sieve_strided */
sieve_cmp(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
          Params P, u32 *cnt, u64 *hits)
{
    __shared__ unsigned char perm[256];             /* per warp: dest lane -> source lane */
    const u32 FULL = 0xFFFFFFFFu;
    const u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    const u32 lane = threadIdx.x & 31;
    unsigned char *wperm = perm + (threadIdx.x & ~31u);
    bool have = gid < P.nthreads;                   /* no early return: warp-wide shuffles */
    if (!CMP && !have) return;
    u64 A1 = have ? Rpre[gid] : 0, Ag = A1;
    u32 i1 = 0, i2 = 0;
    const u64 sG = P.s2 * NCH;
    u64 pA = 0, ps[NCH];
    #pragma unroll
    for (int c = 0; c < NCH; c++) ps[c] = 0;
    u32 npend = 0;                                  /* warp-uniform */
    while (CMP ? __any_sync(FULL, have) : have) {
        u64 s[NCH];
        #pragma unroll
        for (int c = 0; c < NCH; c++) s[c] = (have && i2 + c < P.c2) ? ~0ULL : 0ULL;
        {
            const u32 l0 = (u32)(Ag & 0xFFFF), l1 = (u32)((Ag >> 16) & 0xFFFF);
            const u32 l2 = (u32)((Ag >> 32) & 0xFFFF), l3 = (u32)(Ag >> 48);
            #pragma unroll
            for (int tb = 0; tb < T0; tb += CMP_U) {
                u64 o = 0;
                #pragma unroll
                for (int c = 0; c < NCH; c++) o |= s[c];
                const bool live = (tb < CMP_GUARD0) || o != 0ULL;   /* dead lanes stop loading */
                #pragma unroll
                for (int u = 0; u < CMP_U; u++)
                    row_and<NCH>(words, tb + u, (tb + u) < TV, live, l0, l1, l2, l3, s);
            }
        }
        const u64 Acur = Ag;
        if (have) {                                 /* advance this lane's fresh-group iterator */
            i2 += NCH;
            if (i2 < P.c2) Ag += sG;
            else { i2 = 0; if (++i1 < P.c1) { A1 += P.s1; Ag = A1; } else have = false; }
        }
        if (!CMP) {                                 /* no compaction: finish own group */
            tail_rows<NCH>(words, T0, P.ntc, Acur, s);
            emit<NCH>(Acur, s, P.s2, P.cap, cnt, hits);
            continue;
        }
        u64 o = 0;
        #pragma unroll
        for (int c = 0; c < NCH; c++) o |= s[c];
        u32 m = __ballot_sync(FULL, o != 0ULL);
        while (m) {                                 /* at most two passes */
            const u32 k = __popc(m), take = min(k, 32u - npend);
            const u32 rank = __popc(m & ((1u << lane) - 1u));
            if (((m >> lane) & 1u) && rank < take) wperm[npend + rank] = (unsigned char)lane;
            __syncwarp();
            const bool recv = lane >= npend && lane < npend + take;
            const u32 src = recv ? wperm[lane] : lane;
            __syncwarp();
            const u64 vA = __shfl_sync(FULL, Acur, src);
            if (recv) pA = vA;
            #pragma unroll
            for (int c = 0; c < NCH; c++) {
                const u64 v = __shfl_sync(FULL, s[c], src);
                if (recv) ps[c] = v;
            }
            npend += take;
            if (take == k) m = 0;                   /* the `take` lowest set bits are pending now */
            else for (u32 j = 0; j < take; j++) m &= m - 1;
            if (npend == 32) {
                tail_rows<NCH>(words, T0, P.ntc, pA, ps);
                emit<NCH>(pA, ps, P.s2, P.cap, cnt, hits);
                npend = 0;
            }
        }
    }
    if (CMP && npend) {                             /* flush the partial pending batch */
        if (lane >= npend) {
            #pragma unroll
            for (int c = 0; c < NCH; c++) ps[c] = 0;
        }
        tail_rows<NCH>(words, T0, P.ntc, pA, ps);
        emit<NCH>(pA, ps, P.s2, P.cap, cnt, hits);
    }
}

/* (T0, TV, CMP) combinations instantiated for --kernel 30 */
#define CMP_COMBOS(X) \
    X(28, 0, 1) X(24, 0, 1) X(20, 0, 1) X(32, 0, 1) \
    X(16, 12, 0) X(16, 16, 0) X(12, 12, 0) X(20, 20, 0) \
    X(28, 12, 1) X(28, 16, 1) X(24, 12, 1) X(24, 16, 1) X(32, 12, 1)


/* sieve_cmpw<NCH, T0>  (--kernel 31 --t0 T0): like sieve_cmp<NCH, T0, 0, 1>, but the
 * compaction packs live WORDS (one residue's 64 candidates), not whole groups: after
 * T0 rows a surviving group nearly always has one live word, so the tail does one
 * index + one load per lane-row instead of NCH loads, and its footprint in the
 * (L2-resident) late rows is one sector instead of a whole line. */
#ifndef CMPW_U
#define CMPW_U 2
#endif
__device__ __forceinline__ void tail_word(const u64 *__restrict__ words, u32 t0, u32 ntc,
        u64 A, u64 &s)
{
    if (s == 0ULL) return;
    const u32 l0 = (u32)(A & 0xFFFF), l1 = (u32)((A >> 16) & 0xFFFF);
    const u32 l2 = (u32)((A >> 32) & 0xFFFF), l3 = (u32)(A >> 48);
    for (u32 t = t0; t < ntc; t += CMPW_U) {
        #pragma unroll
        for (int u = 0; u < CMPW_U; u++) {
            u32 p = c_rp[t + u];
            u32 x = l0 * c_k0[t + u] + l1 * c_k1[t + u] + l2 * c_k2[t + u] + l3 * c_k3[t + u];
            u32 y = x - __umulhi(x, c_mag[t + u]) * p;
            if (y >= p) y -= p;
            if (y >= p) y -= p;
            s &= __ldg(words + c_off[t + u] + y);
        }
        if (s == 0ULL) break;
    }
}
__device__ __forceinline__ void emit_word(u64 A, u64 s, u32 cap, u32 *cnt, u64 *hits)
{
    if (s != 0ULL) {
        u32 k = atomicAdd(cnt, 1u);
        if (k < cap) { hits[2 * k] = A; hits[2 * k + 1] = s; }
    }
}

template <int NCH, int T0>
__global__ void __launch_bounds__(256, CMP_MINB)
sieve_cmpw(const u64 *__restrict__ Rpre, const u64 *__restrict__ words,
           Params P, u32 *cnt, u64 *hits)
{
    __shared__ unsigned char perm[256];
    const u32 FULL = 0xFFFFFFFFu;
    const u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
    const u32 lane = threadIdx.x & 31;
    unsigned char *wperm = perm + (threadIdx.x & ~31u);
    bool have = gid < P.nthreads;
    u64 A1 = have ? Rpre[gid] : 0, Ag = A1;
    u32 i1 = 0, i2 = 0;
    const u64 sG = P.s2 * NCH;
    u64 pA = 0, ps = 0;
    u32 npend = 0;                                  /* warp-uniform */
    while (__any_sync(FULL, have)) {
        u64 s[NCH];
        #pragma unroll
        for (int c = 0; c < NCH; c++) s[c] = (have && i2 + c < P.c2) ? ~0ULL : 0ULL;
        {
            const u32 l0 = (u32)(Ag & 0xFFFF), l1 = (u32)((Ag >> 16) & 0xFFFF);
            const u32 l2 = (u32)((Ag >> 32) & 0xFFFF), l3 = (u32)(Ag >> 48);
            #pragma unroll
            for (int tb = 0; tb < T0; tb += CMP_U) {
                u64 o = 0;
                #pragma unroll
                for (int c = 0; c < NCH; c++) o |= s[c];
                const bool live = (tb < CMP_GUARD0) || o != 0ULL;
                #pragma unroll
                for (int u = 0; u < CMP_U; u++)
                    row_and<NCH>(words, tb + u, false, live, l0, l1, l2, l3, s);
            }
        }
        const u64 Acur = Ag;
        if (have) {
            i2 += NCH;
            if (i2 < P.c2) Ag += sG;
            else { i2 = 0; if (++i1 < P.c1) { A1 += P.s1; Ag = A1; } else have = false; }
        }
        u32 lw = 0;                                 /* this lane's live words */
        #pragma unroll
        for (int c = 0; c < NCH; c++) lw |= (s[c] != 0ULL) ? (1u << c) : 0u;
        u32 m = __ballot_sync(FULL, lw != 0u);
        while (m) {
            const u32 c0 = lw ? (u32)(__ffs((int)lw) - 1) : 0u;
            u64 cs = 0;
            #pragma unroll
            for (int c = 0; c < NCH; c++) if ((u32)c == c0) cs = s[c];
            const u64 cA = Acur + (u64)c0 * P.s2;
            const u32 k = __popc(m), take = min(k, 32u - npend);
            const u32 rank = __popc(m & ((1u << lane) - 1u));
            const bool sent = ((m >> lane) & 1u) && rank < take;
            if (sent) wperm[npend + rank] = (unsigned char)lane;
            __syncwarp();
            const bool recv = lane >= npend && lane < npend + take;
            const u32 src = recv ? wperm[lane] : lane;
            __syncwarp();
            const u64 vA = __shfl_sync(FULL, cA, src);
            const u64 vs = __shfl_sync(FULL, cs, src);
            if (recv) { pA = vA; ps = vs; }
            if (sent) lw &= lw - 1u;
            npend += take;
            if (npend == 32) {
                tail_word(words, T0, P.ntc, pA, ps);
                emit_word(pA, ps, P.cap, cnt, hits);
                npend = 0;
            }
            m = __ballot_sync(FULL, lw != 0u);
        }
    }
    if (npend) {
        if (lane >= npend) ps = 0;
        tail_word(words, T0, P.ntc, pA, ps);
        emit_word(pA, ps, P.cap, cnt, hits);
    }
}
#define CMPW_COMBOS(X) X(7, 16) X(7, 20) X(7, 24) X(7, 28) X(7, 32) X(7, 36) X(6, 24) X(6, 28) X(8, 28) X(8, 24) X(5, 24) X(6, 20)

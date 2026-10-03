/* loesch_core_api.c -- C linkage shim so the CUDA engine (C++) uses the SAME
 * loesch_core.h as the CPU engine. The header is C with static functions, so
 * this translation unit instantiates it once and exports thin wrappers.
 *
 * Build: gcc -O3 -march=native -std=gnu11 -c -o build/loesch_core_api.o src/c/loesch_core_api.c
 */
#define _POSIX_C_SOURCE 200809L
#include <signal.h>
#include <time.h>
#include "loesch_core.h"

void lc_build_small_primes(u64 limit) { build_small_primes(limit); }
u64  lc_n_sp(void) { return n_sp; }
u32  lc_sp(u64 i) { return sp[i]; }
bool lc_is_loeschian(u64 t) { return is_loeschian(t); }
u64  lc_mulmod(u64 a, u64 b, u64 m) { return mulmod(a, b, m); }
u64  lc_inv_mod(u64 a, u64 m) { return inv_mod(a, m); }
double lc_now_s(void) { return now_s(); }

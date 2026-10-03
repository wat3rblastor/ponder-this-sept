#ifndef LOESCH_CORE_API_H
#define LOESCH_CORE_API_H
#include <stdint.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif
void     lc_build_small_primes(uint64_t limit);
uint64_t lc_n_sp(void);
uint32_t lc_sp(uint64_t i);
bool     lc_is_loeschian(uint64_t t);
uint64_t lc_mulmod(uint64_t a, uint64_t b, uint64_t m);
uint64_t lc_inv_mod(uint64_t a, uint64_t m);
double   lc_now_s(void);
#ifdef __cplusplus
}
#endif
#endif

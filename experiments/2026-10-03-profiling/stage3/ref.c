/* reference = the pre-Montgomery loesch_core.h (git 7f26d8f), kept for the equivalence test */
#define _POSIX_C_SOURCE 200809L
#include <signal.h>
#include <time.h>
#include "loesch_core_ref.h"
void ref_init(u64 l) { build_small_primes(l); }
bool ref_is(u64 t) { return is_loeschian(t); }

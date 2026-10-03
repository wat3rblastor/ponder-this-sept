#define _POSIX_C_SOURCE 200809L
#include <signal.h>
#include <time.h>
#include <stdio.h>
#include "loesch_core.h"
static double tcpu(void){struct timespec t;clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static u64 rs=88172645463325252ULL; static u64 rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
/* passes the GPU sieve's per-term condition: every bad prime <= 1e4 to an even power */
static int sieved_ok(u64 t){ int e=__builtin_ctzll(t); if(e&1) return 0; t>>=e;
  for(u64 i=1;i<n_sp;i++){u32 p=sp[i]; if(p%3!=2) continue; int k=0; while(t%p==0){t/=p;k++;} if(k&1) return 0;} return 1;}
int main(int argc,char**argv){
  build_small_primes(10000);
  double lo=atof(argv[1]), hi=atof(argv[2]); int N=atoi(argv[3]);
  u64 *v=malloc(N*8); int n=0;
  while(n<N){u64 t=(u64)(lo+(hi-lo)*(rnd()>>11)*(1.0/9007199254740992.0)); if(sieved_ok(t)) v[n++]=t;}
  int pass=0; double t0=tcpu(); for(int i=0;i<N;i++) pass+=is_loeschian(v[i]); double dt=tcpu()-t0;
  printf("range %.1e..%.1e N=%d pass=%.3f  %.2f us/test (thread cpu)\n",lo,hi,N,(double)pass/N,1e6*dt/N);
}

/* density: empirical density of Loeschian numbers that are ==1 mod 3 and
 * coprime to a given list of primes, measured at a ladder of cutoffs.
 * usage: ./density <N> [p1 p2 ...] */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
typedef uint64_t u64;
static u64 *L; static u64 Nmax;
static inline int getb(u64 n){return (L[n>>6]>>(n&63))&1u;}
static inline void setb(u64 n){L[n>>6]|=(u64)1<<(n&63);}
int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: density N [primes...]\n");return 2;}
  u64 N=strtoull(argv[1],NULL,10); Nmax=N;
  int np=argc-2; u64 pr[32]; for(int i=0;i<np;i++) pr[i]=strtoull(argv[i+2],NULL,10);
  L=calloc((N>>6)+1,8); if(!L){fprintf(stderr,"oom\n");return 2;}
  for(u64 x=0;x*x<=N;x++){u64 xx=x*x; for(u64 y=0;;y++){u64 v=xx+x*y+y*y; if(v>N)break; setb(v);} }
  /* sanity */
  if(getb(2)||getb(5)||!getb(4)||!getb(7)){fprintf(stderr,"SIEVE BAD\n");return 2;}
  u64 den=0,num=0; u64 next=10;
  printf("T,eligible,loeschian,density\n");
  for(u64 n=1;n<=N;n++){
    if(n%3!=1) goto tick;
    for(int i=0;i<np;i++) if(n%pr[i]==0) goto tick;
    den++; if(getb(n)) num++;
tick:
    if(n==next||n==N){printf("%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.8f\n",n,den,num,den?(double)num/den:0.0);
      if(n==next) next*=10;}
  }
  return 0;
}

#define _GNU_SOURCE
#include <signal.h>
#include <time.h>
#include "/workspace/ponder-this-sept/src/c/loesch_core.h"
#include <omp.h>
#include <math.h>
static inline u64 rng(u64 *s){ u64 x=*s; x^=x<<13; x^=x>>7; x^=x<<17; return *s=x; }
int main(int argc,char**argv){
  build_small_primes(10000);
  int nb=0; static u32 bad[2000];
  for(u64 i=0;i<n_sp;i++) if(sp[i]%3==2) bad[nb++]=sp[i];
  long N=atol(argv[1]);
  for(int ai=2;ai<argc;ai++){
    double T=atof(argv[ai]);
    u64 lo=(u64)(T*0.8), span=(u64)(T*0.4);
    long tot=0,ok=0; long c9[3]={0},o9[3]={0}; long c7[2]={0},o7[2]={0};
    #pragma omp parallel num_threads(24) reduction(+:tot,ok,c9[:3],o9[:3],c7[:2],o7[:2])
    { u64 s=0x9E3779B97F4A7C15ULL*(omp_get_thread_num()+1)+ai;
      long mine=N/24;
      while(mine>0){
        u64 t=lo+ (u64)(((u128)rng(&s)*span)>>64);
        t-= t%3; t+=1;
        int good=1; for(int i=0;i<nb;i++) if(t%bad[i]==0){good=0;break;}
        if(!good) continue;
        mine--; tot++; int L=is_loeschian(t); ok+=L;
        int r=(t%9)/3; c9[r]++; o9[r]+=L;
        int g=(t%7==0); c7[g]++; o7[g]+=L;
      }
    }
    double r=(double)ok/tot;
    printf("T=%.3g lnT=%.3f rho=%.5f +-%.5f | mod9: %.4f %.4f %.4f | 7∤t %.4f 7|t %.4f\n",T,log(T),r,sqrt(r*(1-r)/tot),
      (double)o9[0]/c9[0],(double)o9[1]/c9[1],(double)o9[2]/c9[2],(double)o7[0]/c7[0],(double)o7[1]/c7[1]);
    fflush(stdout);
  }
}

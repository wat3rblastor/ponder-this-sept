/* mc_win -- per-term pass probability measured on REAL apsearch AP terms.
 * For random a (satisfying the tier-A conditions) and d=K*D0, look at the 58
 * terms a+k*d, keep those free of every bad prime <= b2 (exactly what the
 * tier-C sieve guarantees for an admissible window), and measure how often
 * such a term is Loeschian.  Also records, per sampled window, how many of
 * the 58 terms were free, so the window-level and term-level effects separate.
 * Usage: mc_win K b2 nwindows seed */
#include <signal.h>
#include <time.h>
#include <math.h>
#include "loesch_core.h"
static u64 rs;
static inline u64 rnd(void){u64 z=(rs+=0x9E3779B97F4A7C15ULL);
 z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31);}
int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage: mc_win K b2 nwin seed\n");return 2;}
  u64 K=strtoull(argv[1],NULL,10), b2=strtoull(argv[2],NULL,10);
  u64 nwin=strtoull(argv[3],NULL,10); rs=strtoull(argv[4],NULL,10);
  build_small_primes(10000);
  u64 d=K*D0_DEFAULT;
  u64 bad[512]; int nbad=0;
  for(u64 i=0;i<n_sp;i++){ if(sp[i]>b2) break; if(sp[i]%3==2) bad[nbad++]=sp[i]; }
  u64 free_terms=0, loe=0, tot_terms=0;
  for(u64 w=0;w<nwin;w++){
    u64 a;
    for(;;){
      a=(u64)1e17 + rnd()%(u64)3e17;
      if(a%3!=1) continue;
      int ok=1;
      for(int i=0;i<nbad;i++) if(d%bad[i]==0 && a%bad[i]==0){ok=0;break;}
      if(ok) break;
    }
    for(int k=0;k<58;k++){
      u64 t=a+(u64)k*d; tot_terms++;
      int fr=1;
      for(int i=0;i<nbad;i++) if(t%bad[i]==0){fr=0;break;}
      if(!fr) continue;
      free_terms++;
      if(is_loeschian(t)) loe++;
    }
  }
  printf("K=%llu d=%llu b2=%llu windows=%llu terms=%llu free=%llu loe=%llu "
         "freerate=%.5f q_ap=%.5f sigma=%.5f\n",
    (unsigned long long)K,(unsigned long long)d,(unsigned long long)b2,
    (unsigned long long)nwin,(unsigned long long)tot_terms,
    (unsigned long long)free_terms,(unsigned long long)loe,
    (double)free_terms/tot_terms,(double)loe/free_terms,
    sqrt(((double)loe/free_terms)*(1-(double)loe/free_terms)/free_terms));
  return 0;}

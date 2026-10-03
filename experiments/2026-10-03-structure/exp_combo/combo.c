// combo.c: unified run-length measurement for AP families of Loeschian numbers (u64 terms).
// modes: G d Amax seed secs      general a (a=1 mod 3 random < Amax), fixed d
//        A m Xmax seed secs      a=x^2, d=3m^2, x random < Xmax
//        A2 M0 jmax              a=x^2, d=3m^2, m=M0*j, and x^2+3*k1*m^2 = z^2 (enumerated)
#include <time.h>
#include <signal.h>
#include "/workspace/ponder-this-sept/src/c/loesch_core.h"
#define N 57
#define TH 16
static u64 rs;
static u64 rnd(void){ u64 z=(rs+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31);}
static u32 qs[400], qinv[400]; static int nq;
static u64 hist[N+2], sqh[N+2], ncand, nexact;
static int longest(const unsigned char *f, int *st){ int best=0,cur=0; for(int i=0;i<N;i++){ if(!f[i]){cur++; if(cur>best){best=cur;*st=i-cur+1;}} else cur=0;} return best;}
static void setup(u64 d){
  nq=0; for(u32 q=2;q<1000;q++){ if(q%3!=2) continue; int pr=1; for(u32 e=2;e*e<=q;e++) if(q%e==0) pr=0; if(!pr) continue;
    if(d%q==0) continue; qs[nq]=q; qinv[nq]=(u32)inv_mod(d%q,q); nq++; }
}
static u64 badd[32]; static int nbadd;
static int test(u64 a,u64 d,int report){
  ncand++;
  unsigned char f[N]; memset(f,0,N); int sq=0, st=0;
  for(int i=0;i<nq;i++){ u32 q=qs[i]; u32 h=(u32)(((u64)(q-(u32)(a%q))%q)*qinv[i]%q); int ch=0;
    for(u32 k=h;k<N;k+=q){ u64 t=a+(u64)k*d; int e=0; while(t%q==0){t/=q;e++;} if(e&1){ if(!f[k]){f[k]=1;ch=1;} } else sq=1; }
    if(ch && longest(f,&st)<TH) return 0; }
  nexact++;
  for(int k=0;k<N;k++) if(!f[k]) { if(!is_loeschian(a+(u64)k*d)){ f[k]=1; if(longest(f,&st)<TH) return 0;} }
  int L=longest(f,&st);
  if(L>=TH){ hist[L]++; if(sq) sqh[L]++; }
  if(L>=38 && report){ printf("RUN %d a=%llu d=%llu start=%d (a_start=%llu)\n",L,(unsigned long long)a,(unsigned long long)d,st,(unsigned long long)(a+(u64)st*d)); fflush(stdout);}
  return L;
}
static void summary(const char*tag){
  double cpu=(double)clock()/CLOCKS_PER_SEC; u64 c=0,s=0; int mx=0;
  printf("%s cand=%llu exact=%llu cpu=%.1f\n",tag,(unsigned long long)ncand,(unsigned long long)nexact,cpu);
  for(int L=N;L>=TH;L--){ c+=hist[L]; s+=sqh[L]; if(hist[L]&&!mx) mx=L; if(L==16||L==20||L==25||L==30||L==35||L==40||L==45) printf("  >=%d: %llu (sq %llu)  per_cand=%.3e per_cpuh=%.3e\n",L,(unsigned long long)c,(unsigned long long)s,(double)c/ncand,(double)c/cpu*3600); }
  printf("  longest=%d\n",mx);
}
int main(int argc,char**argv){
  build_small_primes(10000);
  if(!strcmp(argv[1],"G")){ u64 d=strtoull(argv[2],0,10), Amax=strtoull(argv[3],0,10); rs=strtoull(argv[4],0,10); double secs=atof(argv[5]);
    setup(d);
    while((double)clock()/CLOCKS_PER_SEC<secs) for(int i=0;i<100000;i++){ u64 a=rnd()%Amax; a=a-a%3+1; test(a,d,1);} 
    summary(argv[1]);
  } else if(!strcmp(argv[1],"A")){ u64 m=strtoull(argv[2],0,10), Xmax=strtoull(argv[3],0,10); rs=strtoull(argv[4],0,10); double secs=atof(argv[5]);
    u64 d=3*m*m; setup(d);
    while((double)clock()/CLOCKS_PER_SEC<secs) for(int i=0;i<100000;i++){ u64 x=rnd()%Xmax+1; if(x%3==0||gcd_u64(x,m)!=1) continue; test(x*x,d,1);} 
    summary(argv[1]);
  } else if(!strcmp(argv[1],"A2")){ u64 M0=strtoull(argv[2],0,10); int jmax=atoi(argv[3]); int kmax=argc>4?atoi(argv[4]):56;
    for(int j=1;j<=jmax;j++){ u64 m=M0*j; if((long double)3*m*m*56>1.7e19L) break; u64 d=3*m*m; setup(d);
      for(int k1=1;k1<=kmax;k1++){ u64 NN=3*(u64)k1*m*m; // enumerate divisors
        u64 pp[32]; int ee[32], np=0; u64 t=NN; for(u64 p=2;p*p<=t;p++) if(t%p==0){ pp[np]=p; ee[np]=0; while(t%p==0){t/=p;ee[np]++;} np++; } if(t>1){pp[np]=t;ee[np]=1;np++;}
        static u64 dv[2000000]; int nd=1; dv[0]=1;
        for(int i=0;i<np;i++){ int old=nd; u64 pw=1; for(int e=1;e<=ee[i];e++){ pw*=pp[i]; for(int r=0;r<old;r++) dv[nd++]=dv[r]*pw; } }
        for(int r=0;r<nd;r++){ u64 u=dv[r], v=NN/u; if(u>=v||((u^v)&1)) continue; u64 x=(v-u)/2; if(x>4200000000ULL) continue; if(x%3==0||gcd_u64(x,m)!=1) continue;
          if((long double)x*x+56.0L*d>1.8e19L) continue; test(x*x,d,1); }
      } }
    summary("A2");
  }
  return 0;
}

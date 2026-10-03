// independent sampler: enumerate sieve-clean 58-windows for (K, shift) [reduced planes: multiples 64s..64s+63]
// and record exact per-position Loeschian results. usage: samp K shift [modcap]
#define _GNU_SOURCE
#include <signal.h>
#include <time.h>
#include "/workspace/ponder-this-sept/src/c/loesch_core.h"
#include <omp.h>
#include <math.h>
#define NT 58
#define D0 382160924970ULL
typedef struct { u64 m, s, t, c; } comp;
int main(int argc, char **argv) {
  u64 K = strtoull(argv[1],0,10), shift = strtoull(argv[2],0,10);
  double modcap = argc>3 ? atof(argv[3]) : 2e16;
  build_small_primes(10000);
  u64 d = K * D0;
  comp C[16]; int nc=0; u64 MOD=1;
  #define ADDC(m_,s_,t_,c_) do{C[nc].m=(m_);C[nc].s=(s_);C[nc].t=(t_);C[nc].c=(c_);MOD*=(m_);nc++;}while(0)
  ADDC(3,1,0,1); ADDC(2,1,0,1); ADDC(5,1,1,4);
  for (u64 i=0;i<n_sp;i++){ u64 q=sp[i]; if(q<=NT||q%3!=2||d%q==0) continue;
    if ((double)MOD*q>modcap) break; u64 dq=d%q; ADDC(q,dq,dq,q-NT); }
  u64 R0=0, ss[16];
  for (int i=0;i<nc;i++){ u64 mi=C[i].m, co=MOD/mi; u64 e=mulmod(co%MOD, inv_mod(co%mi,mi), MOD);
    R0=(R0+mulmod(C[i].s,e,MOD))%MOD; ss[i]=mulmod(C[i].t,e,MOD); }
  // tier C
  static u32 tc[2000]; int ntc=0;
  for (u64 i=0;i<n_sp;i++){ u64 r=sp[i]; int pinned=0; for(int t=0;t<nc;t++) if(C[t].m==r) pinned=1;
    if(pinned) continue; int divD=(D0%r==0); if(!divD&&(r%3!=2||d%r==0||r<=NT)) continue; tc[ntc++]=r; }
  // order by kill
  for (int i=0;i<ntc;i++){int bj=i;double bk=-1;for(int j=i;j<ntc;j++){double k=(D0%tc[j]==0)?1.0/tc[j]:(double)NT/tc[j]; if(k>bk){bk=k;bj=j;}} u32 t=tc[i];tc[i]=tc[bj];tc[bj]=t;}
  u64 **W = malloc(ntc*sizeof(u64*));
  for (int t=0;t<ntc;t++){ u32 r=tc[t]; W[t]=malloc(r*8); int divD=(D0%r==0);
    for (u32 x=0;x<r;x++){ u64 w=0; for(int b=0;b<64;b++){ u64 a=(x+(u64)((shift*64+b)%r)*(MOD%r))%r; int ok;
        if(divD) ok=(a!=0); else { ok=1; u64 dr=d%r; u64 v=a; for(int k=0;k<NT;k++){ if(v==0){ok=0;break;} v+=dr; if(v>=r) v-=r; } }
        if(ok) w|=1ULL<<b; } W[t][x]=w; } }
  // enumerate: outer loop over first non-trivial comps flattened
  int idx[16]; int act[16], na=0; for(int i=0;i<nc;i++) if(C[i].c>1) act[na++]=i;
  u64 total=1; for(int i=0;i<na;i++) total*=C[act[i]].c;
  u64 inner = C[act[na-1]].c * C[act[na-2]].c; u64 nouter = total/inner;
  fprintf(stderr,"K=%llu shift=%llu MOD=%llu res=%llu ntc=%d\n",(unsigned long long)K,(unsigned long long)shift,(unsigned long long)MOD,(unsigned long long)total,ntc);
  long pos_ok[NT]={0}, nwin=0, ndead=0; long runh[NT+1]={0};
  double sumlog[NT]={0};
  #pragma omp parallel for schedule(dynamic,64) num_threads(28) reduction(+:pos_ok[:NT],nwin,ndead,runh[:NT+1])
  for (u64 o=0;o<nouter;o++){
    u64 R=R0, oo=o;
    for (int i=0;i<na-2;i++){ u64 j=oo%C[act[i]].c; oo/=C[act[i]].c; R=(R+mulmod(j,ss[act[i]],MOD))%MOD; }
    u64 s1=ss[act[na-2]], s2=ss[act[na-1]];
    u64 R1=R;
    for (u64 i1=0;i1<C[act[na-2]].c;i1++){
      u64 R2=R1;
      for (u64 i2=0;i2<C[act[na-1]].c;i2++){
        u64 w=~0ULL;
        for (int t=0;t<ntc && w;t++) w&=W[t][R2%tc[t]];
        while (w){ int b=__builtin_ctzll(w); w&=w-1;
          u64 a0=R2+(shift*64+b)*MOD;
          if (K%7==0 && a0%7==0) continue;
          nwin++; int cur=0,best=0;
          int Ls[NT]; int tot=0; for (int k=0;k<NT;k++){ int L=is_loeschian(a0+(u64)k*d); Ls[k]=L; tot+=L; if(L){cur++; if(cur>best)best=cur;} else cur=0; }
          runh[best]++; if (tot<=3) { nwin--; ndead++; } else for (int k=0;k<NT;k++) pos_ok[k]+=Ls[k];
        }
        R2+=s2; if(R2>=MOD) R2-=MOD;
      }
      R1+=s1; if(R1>=MOD) R1-=MOD;
    }
  }
  // model rho at window-mean a
  double amid=(shift*64+32.0)*MOD;
  double ok=0, ex=0;
  printf("K=%llu shift=%llu alive=%ld dead=%ld Tlast~%.3g\n",(unsigned long long)K,(unsigned long long)shift,nwin,ndead,amid+57.0*d);
  for (int k=0;k<NT;k++){ double T=amid+(double)k*d; double rho=0.71992*pow(log(T)/39.144,-0.485); ok+=pos_ok[k]; ex+=rho*nwin; }
  printf("per-term pass: obs %.5f  model %.5f  ratio %.4f  (+-%.4f)\n", ok/(NT*nwin), ex/(NT*nwin), ok/ex, sqrt(0.21/(NT*(double)nwin))/0.7);
  printf("pos:"); for(int k=0;k<NT;k+=3) printf(" %.3f", (double)pos_ok[k]/nwin); printf("\n");
  printf("runhist:"); for(int n=0;n<=NT;n++) printf(" %ld", runh[n]); printf("\n");
}

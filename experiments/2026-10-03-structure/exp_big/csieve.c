#define main fam_main
#include "fam.c"
#undef main
#include <math.h>
// usage: csieve t n Q R seconds mode   (mode 0: C family x-sieve survivors; 1: per-term-conditioned control: nonsquare a same d; 2: 18-digit baseline control)
static u64 powm(u64 b,u64 e,u64 q){ u64 r=1; b%=q; while(e){ if(e&1) r=r*b%q; b=b*b%q; e>>=1; } return r; }
int main(int argc,char**argv){ u64 t=atoll(argv[1]); int n=atoi(argv[2]); int Q=atoi(argv[3]); double R=atof(argv[4]); double T=atof(argv[5]); int mode=atoi(argv[6]); sieve(); BUD128=1L<<20;
  u128 P=1; for(int i=0;i<npr&&pr[i]<=53;i++) if(isbad[i]) P*=pr[i]; u128 m=P*t, d=3*m*m; char sd[50]; p128(d,sd);
  int isfree[64]={0}; int nfree=0; for(int j=0;j*j<n;j++){ isfree[j*j]=1; nfree++; }
  static u32 qs[2000]; static unsigned char *ok[2000]; int nq=0; double dens=1,densb=1;
  for(int i=0;i<npr&&pr[i]<=(u32)Q;i++){ u32 q=pr[i]; if(!isbad[i]||q<=53||m%q==0) continue; qs[nq]=q; ok[nq]=malloc(q); u64 mi=powm((u64)(m%q),q-2,q), i3=powm(3,q-2,q); int cnt=0;
    for(u32 r=0;r<q;r++){ u64 z=r*mi%q; u64 k0=(q-z*z%q*i3%q)%q; int good=(r==0)||(k0>=(u64)n); ok[nq][r]=good; cnt+=good; } dens*=(double)cnt/q; densb*=(double)(q-n)/q; nq++; }
  double Xhi=R*(double)m;
  printf("t=%llu n=%d Q=%d d=%s (%d digits) free=%d nonfree=%d | sieve density C=%.3e (survivors per 1e9 x: %.3g) baseline-analog prod(1-n/q)=%.3e ratio=%.1f | x<=%.3g -> supply %.3g\n",(unsigned long long)t,n,Q,sd,digits(d),nfree,n-nfree,dens,dens*1e9,densb,dens/densb,Xhi,Xhi*dens*2/3);
  if(mode==0){ int nc=6; u128 M=1; for(int i=0;i<nc;i++) M*=qs[i]; 
    long surv=0,tp=0,tt=0,tu=0,hist[64]={0},histp[64]={0},freebad=0; long ppos[64]={0}; double t0=omp_get_wtime(); long walked=0;
    #pragma omp parallel
    { u64 s=smix(omp_get_thread_num()*7919+t*31+n); 
      while(omp_get_wtime()-t0<T){ u128 x0=0,MM=1; for(int i=0;i<nc;i++){ u32 q=qs[i]; u32 r; do{ s=smix(s); r=s%q; }while(!ok[i][r]); u64 inv=powm((u64)(MM%q),q-2,q); u64 cur=(u64)(x0%q); u64 k=(u64)((r+q-cur)%q)*inv%q; x0+=MM*k; MM*=q; }
        long lw=0; for(u128 x=x0;(double)x<Xhi;x+=M){ lw++; int good=1; for(int i=nc;i<nq;i++) if(!ok[i][(u64)(x%qs[i])]){ good=0; break; } if(!good) continue; if(x%3==0||gcd(x,m)!=1) continue;
          int res[64]; int run=-1,runp=-1; long lp=0,lu=0,lt=0,fb=0; for(int k=0;k<n;k++){ res[k]=loesch(x*x+(u128)k*d); if(isfree[k]){ if(res[k]==0) fb++; if(res[k]==2) res[k]=1; continue; } lt++; if(res[k]) lp++; if(res[k]==2) lu++; }
          for(int k=0;k<n;k++){ if(res[k]==0&&run<0) run=k; if(res[k]!=1&&runp<0) runp=k; } if(run<0)run=n; if(runp<0)runp=n;
          #pragma omp critical
          { surv++; tp+=lp; tu+=lu; tt+=lt; hist[run]++; histp[runp]++; freebad+=fb; for(int k=0;k<n;k++) if(res[k]) ppos[k]++; if(run>=30){ char sx[50]; p128(x,sx); printf("  LONG run=%d pes=%d x=%s a=x^2 d=%s\n",run,runp,sx,sd); } }
        }
        #pragma omp atomic
        walked+=lw; } }
    double po=(double)tp/tt, pp=(double)(tp-tu)/tt; int nf=n-nfree;
    printf("  survivors tested=%ld (cpu-s %.0f) nonfree per-term pass: optimistic %.4f proven %.4f (unknown %.2f%%) freebad=%ld\n",surv,T*omp_get_max_threads(),po,pp,100.0*tu/tt,freebad);
    printf("  P(full %d-run | survivor) ~ p^%d: opt %.3g  pes %.3g ; per 1e12 sieve-walked x: opt %.3g pes %.3g\n",n,nf,pow(po,nf),pow(pp,nf),pow(po,nf)*dens*1e12*2/3,pow(pp,nf)*dens*1e12*2/3);
    printf("  run-length from k=0 (opt) counts >=5,10,15,20,25,30,35,40: "); for(int L=5;L<=40;L+=5){ long c=0; for(int k=L;k<=n;k++) c+=hist[k]; printf("%ld ",c); } int mx=0,mxp=0; for(int k=0;k<=n;k++){ if(hist[k])mx=k; if(histp[k])mxp=k; } printf(" max %d (proven %d)\n",mx,mxp);
    printf("  per-position pass (opt), first 30: "); for(int k=0;k<30&&k<n;k++) printf("%.2f ",(double)ppos[k]/surv); printf("\n");
  } else { long tt=0,tp=0,tu=0; double t0=omp_get_wtime();
    #pragma omp parallel
    { u64 s=smix(omp_get_thread_num()*104729+5); while(omp_get_wtime()-t0<T){ s=smix(s); u128 a,dd; if(mode==1){ dd=d; a=((u128)(s>>1)<<20)%((u128)(R*(double)m)*(u128)(R*(double)m)); a|=1; } else { dd=3*P*(s%1000+1); u64 r=smix(s+9)>>4; r|=1ULL<<59; a=r; } a-=a%3; a+=1; if(gcd(a,dd)!=1) continue; u128 v=a+(u128)(s%n)*dd; int clean=1; for(int i=0;i<npr&&pr[i]<=(u32)Q;i++) if(isbad[i]&&v%pr[i]==0){ clean=0; break; } if(!clean) continue; int r=loesch(v);
        #pragma omp critical
        { tt++; if(r)tp++; if(r==2)tu++; } } }
    printf("  control mode %d: sieve-clean terms tested=%ld pass optimistic %.4f proven %.4f\n",mode,tt,(double)tp/tt,(double)(tp-tu)/tt); }
  return 0; }

// Exact (budgeted) Loeschian AP family measurement. gcc -O2 -fopenmp fam.c -o fam
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
typedef unsigned __int128 u128; typedef uint64_t u64; typedef uint32_t u32;
#define PB 65536
static u32 pr[7000]; static int npr; static unsigned char isbad[7000];
static void sieve(void){ static char c[PB+1]; for(int i=2;i<=PB;i++){ if(!c[i]){ pr[npr]=i; isbad[npr]=(i%3==2); npr++; for(long j=(long)i*i;j<=PB;j+=i)c[j]=1; } } }
static void mulfull(u128 a,u128 b,u128*hi,u128*lo){ u64 a0=a,a1=a>>64,b0=b,b1=b>>64; u128 p00=(u128)a0*b0,p01=(u128)a0*b1,p10=(u128)a1*b0,p11=(u128)a1*b1; u128 mid=(p00>>64)+(u64)p01+(u64)p10; *lo=(u64)p00|(mid<<64); *hi=p11+(p01>>64)+(p10>>64)+(mid>>64); }
typedef struct { u128 n,ninv,r1,r2; } mont;
static u128 mm(const mont*M,u128 a,u128 b){ u128 hi,lo,mh,ml; mulfull(a,b,&hi,&lo); u128 m=lo*M->ninv; mulfull(m,M->n,&mh,&ml); u128 t=hi+mh+(lo!=0); if(t>=M->n)t-=M->n; return t; }
static u128 slowmul(u128 a,u128 b,u128 n){ u128 r=0; a%=n; while(b){ if(b&1){ r+=a; if(r>=n)r-=n; } a+=a; if(a>=n)a-=n; b>>=1; } return r; }
static void minit(mont*M,u128 n){ M->n=n; u128 x=n; for(int i=0;i<7;i++) x*=2-n*x; M->ninv=(u128)0-x; M->r1=((u128)0-n)%n; M->r2=slowmul(M->r1,M->r1,n); }
static int sprp(const mont*M,u64 base){ u128 n=M->n,d=n-1; int s=0; while(!(d&1)){d>>=1;s++;} u128 b=mm(M,(u128)base%n,M->r2); if(b==0) return 1; u128 x=M->r1; u128 e=d; u128 bb=b; while(e){ if(e&1)x=mm(M,x,bb); bb=mm(M,bb,bb); e>>=1; } u128 one=M->r1, mone=n-M->r1; if(x==one||x==mone) return 1; for(int i=1;i<s;i++){ x=mm(M,x,x); if(x==mone) return 1; if(x==one) return 0; } return 0; }
static int isprime(u128 n){ if(n<4) return n>=2; if(!(n&1)) return 0; mont M; minit(&M,n); static const u64 bs[]={2,3,5,7,11,13,17,19,23,29,31,37,41}; for(int i=0;i<13;i++) if(!sprp(&M,bs[i])) return 0; return 1; }
static u128 gcd(u128 a,u128 b){ while(b){ u128 t=a%b; a=b; b=t; } return a; }
static u128 isqrt(u128 n){ if(n==0) return 0; u128 x=(u128)1<<64; if(n>>126) x=((u128)1<<64)-1; else { int b=0; u128 t=n; while(t){b++;t>>=1;} x=(u128)1<<((b+1)/2); }
  while(1){ u128 y=(x+n/x)>>1; if(y>=x) break; x=y; } return x; }
// rho brent, budgeted. returns nontrivial factor or 0
static u128 rho(u128 n,long budget,u64 seed){ mont M; minit(&M,n); for(u64 c0=1+seed; ; c0++){ u128 c=mm(&M,c0,M.r2); u128 y=mm(&M,2,M.r2),x=y,q=M.r1,g=1,ys=y; long r=1,used=0; do{ x=y; for(long i=0;i<r;i++){ y=mm(&M,y,y)+c; if(y>=n)y-=n; } long k=0; while(k<r&&g==1){ ys=y; long lim=(r-k<128)?r-k:128; for(long i=0;i<lim;i++){ y=mm(&M,y,y)+c; if(y>=n)y-=n; u128 df=x>y?x-y:y-x; q=mm(&M,q,df); } g=gcd(q,n); k+=lim; used+=lim; if(used>budget&&g==1) return 0; } r*=2; }while(g==1);
    if(g==n){ g=1; long cnt=0; while(g==1){ ys=mm(&M,ys,ys)+c; if(ys>=n)ys-=n; u128 df=x>ys?x-ys:ys-x; g=gcd(df,n); if(++cnt>budget*2) return 0; } }
    if(g!=n&&g!=1) return g; budget/=2; if(budget<64) return 0; } }
static long BUD128=1<<14;
// factor c (all prime factors > PB) into list; returns number of entries; flags composite unknowns
static int factor(u128 c,u128*f,int*comp,int n){ if(c==1) return n; if(c<((u128)1<<32)||isprime(c)){ f[n]=c; comp[n]=0; return n+1; }
  u128 s=isqrt(c); if(s*s==c){ int m=factor(s,f,comp,n); int k=m-n; for(int i=0;i<k;i++){ f[m+i]=f[n+i]; comp[m+i]=comp[n+i]; } return m+k; }
  long bud=(c>>64)?BUD128:(1L<<22); u128 g=rho(c,bud,0); if(!g){ f[n]=c; comp[n]=1; return n+1; }
  n=factor(g,f,comp,n); return factor(c/g,f,comp,n); }
// 1 pass, 0 fail, 2 unknown
static int loesch(u128 n){ if(n==0) return 1; while(n%3==0) n/=3;
  for(int i=0;i<npr;i++){ u32 p=pr[i]; if(p==3) continue; if(!(n>>64)){ u64 m=n; if((u64)p*p>m){ return (m==1||m%3==1)?1:0; } if(m%p==0){ int e=0; do{ m/=p; e++; }while(m%p==0); if(isbad[i]&&(e&1)) return 0; n=m; } }
    else if(n%p==0){ int e=0; do{ n/=p; e++; }while(n%p==0); if(isbad[i]&&(e&1)) return 0; } }
  if(n==1) return 1; if(n%3==2) return 0; if(n<((u128)1<<32)) return 1;
  u128 f[40]; int comp[40]; int k=factor(n,f,comp,0); int unk=0; for(int i=0;i<k;i++) if(comp[i]) unk=1;
  // parity of known bad primes
  for(int i=0;i<k;i++){ if(comp[i]||f[i]%3!=2) continue; int e=0; for(int j=0;j<k;j++) if(!comp[j]&&f[j]==f[i]) e++; if(e&1){ int div=0; for(int j=0;j<k;j++) if(comp[j]&&f[j]%f[i]==0) div=1; if(!div) return 0; unk=1; } }
  return unk?2:1; }
static u64 smix(u64 x){ x+=0x9E3779B97F4A7C15ULL; x=(x^(x>>30))*0xBF58476D1CE4E5B9ULL; x=(x^(x>>27))*0x94D049BB133111EBULL; return x^(x>>31); }
static void p128(u128 x,char*s){ char t[50]; int i=0; if(!x)t[i++]='0'; while(x){ t[i++]='0'+(int)(x%10); x/=10; } int j=0; while(i) s[j++]=t[--i]; s[j]=0; }
static int digits(u128 x){ int n=0; while(x){n++;x/=10;} return n; }
#define W 60
typedef struct { const char*name; int type; u128 d; u128 a0; u64 seed; } cfg;
// types: 0 seq a=a0+3i ; 1 K-family: K=i%1000+1, a=1+6*(i/1000) ; 2 a=x^2,x=a0+i ; 4 a prime seq ; 5 random 60-bit a ; 6 baseline: d*K (K=1..1000), random 60-bit a ; 7 a=1, d*(i+1)
static int gen(const cfg*c,u64 i,u128*a,u128*d){ *d=c->d; switch(c->type){
  case 0: *a=c->a0+3*(u128)i; break;
  case 1: *d=c->d*(i%1000+1); *a=1+6*(u128)(i/1000); break;
  case 2: { u128 x=c->a0+i; if(x%3==0) return 0; *a=x*x; break; }
  case 4: *a=c->a0+3*(u128)i; if(!isprime(*a)) return 0; break;
  case 5: { u64 r=smix(i*2654435761ULL+c->seed)>>4; r|=1ULL<<59; r-=r%3; r+=1; *a=r; break; }
  case 6: { u64 r=smix(i*2654435761ULL+c->seed)>>4; r|=1ULL<<59; r-=r%3; r+=1; *a=r; *d=c->d*(smix(i+77)%1000+1); break; }
  case 7: *a=1; *d=c->d*(u128)(i+1); break; }
  if((*d>>119)||(*a>>124)) return 0;
  return gcd(*a,*d)==1; }
int main(int argc,char**argv){ double T=argc>1?atof(argv[1]):8; int only=argc>2?atoi(argv[2]):-1; sieve();
  u128 Pb[400]; u128 Pall[400]; // product of bad primes <= y ; all primes <= y
  { u128 pb=1,pa=1; int ob=0,oa=0; int j=0; for(int y=0;y<400;y++){ while(j<npr&&pr[j]<=(u32)y){ if(isbad[j]){ if(pb>>100)ob=1; else pb*=pr[j]; } if(pa>>100)oa=1; else pa*=pr[j]; j++; } Pb[y]=ob?0:pb; Pall[y]=oa?0:pa; } }
  cfg C[200]; int nc=0; static char names[200][64];
#define ADD(...) do{ snprintf(names[nc],64,__VA_ARGS__); C[nc].name=names[nc]; }while(0)
  u128 D0=3*Pb[53]; u64 R60=((smix(12345)>>4)|(1ULL<<59)); R60-=R60%3; R60+=1;
  ADD("BASE D0*K rand a60"); C[nc].type=6; C[nc].d=D0; C[nc].a0=0; C[nc].seed=1; nc++;
  ADD("BASE D0 rand a60"); C[nc].type=5; C[nc].d=D0; C[nc].a0=0; C[nc].seed=2; nc++;
  int ys[]={29,53,59,71,89,113,150};
  for(int t=0;t<7;t++){ int y=ys[t]; ADD("A(%d) a small seq",y); C[nc].type=0; C[nc].d=3*Pb[y]; C[nc].a0=1; nc++;
    ADD("A(%d) a seq from 2^60",y); C[nc].type=0; C[nc].d=3*Pb[y]; C[nc].a0=R60; nc++; }
  int yb[]={53,71,89,113};
  for(int t=0;t<4;t++){ int y=yb[t]; ADD("B(%d) d*K a=1,7,13..",y); C[nc].type=1; C[nc].d=3*Pb[y]; nc++; }
  int yc[]={29,41,47,53};
  for(int t=0;t<4;t++){ int y=yc[t]; u128 d=3*Pb[y]*Pb[y]; ADD("C(%d) d=3P^2 a=x^2 small",y); C[nc].type=2; C[nc].d=d; C[nc].a0=1; nc++;
    ADD("C(%d) d=3P^2 a=x^2 x~2^30",y); C[nc].type=2; C[nc].d=d; C[nc].a0=(1u<<30)+12345; nc++;
    ADD("C(%d) d=3P^2 a nonsq ctrl",y); C[nc].type=0; C[nc].d=d; C[nc].a0=R60; nc++; }
  int ye[]={53,59,71,83};
  for(int t=0;t<4;t++){ int y=ye[t]; ADD("E1 d=K*primorial(%d) a=1",y); C[nc].type=7; C[nc].d=Pall[y]; nc++;
    ADD("E1b d=primorial(%d) a seq 2^60",y); C[nc].type=0; C[nc].d=Pall[y]; C[nc].a0=R60; nc++; }
  ADD("E2 D0 a prime seq 2^60"); C[nc].type=4; C[nc].d=D0; C[nc].a0=R60; nc++;
  ADD("E3 D0*7*13*19*31*37 rand a60"); C[nc].type=5; C[nc].d=D0*7*13*19*31*37; C[nc].seed=3; nc++;
  ADD("E3b D0*7*13*19*31*37 a small"); C[nc].type=0; C[nc].d=D0*7*13*19*31*37; C[nc].a0=1; nc++;
  ADD("E4 D0^2/3 (bad primes squared) a60"); C[nc].type=5; C[nc].d=D0*Pb[53]; C[nc].seed=4; nc++;
  ADD("E5 D0*K, a=1 (1+kd by hand)"); C[nc].type=7; C[nc].d=D0; nc++;
  printf("# T=%.1fs wall x %d threads per config. p = conditional pass prob given run reached position (opt: unknown=pass / pes: unknown=fail); u = unconditional pass prob from full-window subsample\n",T,omp_get_max_threads());
  printf("%-34s %3s %9s %7s | %-11s %-11s %-11s | %-6s %-6s | %8s %8s %8s %8s | %8s %8s %8s | %4s %4s | %6s\n","family","dig","cands","c/cpus","p0-9 o/p","p1-9 o/p","p10-29 o/p","u1-9","u10-29",">=5/M",">=10/M",">=15/M",">=20/M","ex>=20/M","ex>=25/M","ex>=30/M","max","maxA","unk%");
  for(int ci=0;ci<nc;ci++){ if(only>=0&&ci!=only) continue; cfg*c=&C[ci]; if(c->d==0||(c->d>>119)){ printf("%-34s overflow, skipped\n",c->name); continue; }
    u64 next=0; double t0=omp_get_wtime(); long hist[W+1]={0},histp[W+1]={0},tested[W]={0},passed[W]={0},unkn[W]={0},ut[W]={0},up[W]={0}; long ncand=0,ntests=0; int maxany=0;
    #pragma omp parallel
    { long h[W+1]={0},hp[W+1]={0},te[W]={0},pa[W]={0},un[W]={0},lut[W]={0},lup[W]={0}; long nca=0,nt=0; int ma=0;
      while(omp_get_wtime()-t0<T){ u64 i0;
        #pragma omp atomic capture
        { i0=next; next+=64; }
        for(u64 i=i0;i<i0+64;i++){ u128 a,d; if(!gen(c,i,&a,&d)) continue; nca++; int run=0,runp=-1; int res[W]; int k;
          for(k=0;k<W;k++){ int r=loesch(a+(u128)k*d); nt++; res[k]=r; te[k]++; if(r==0) break; pa[k]++; if(r==2){ un[k]++; if(runp<0) runp=k; } }
          run=k; if(runp<0) runp=run; h[run]++; hp[runp]++;
          if(run>=30){ char sa[50],sd[50]; p128(a,sa); p128(d,sd);
            #pragma omp critical
            { printf("  LONG run=%d pes=%d a=%s d=%s\n",run,runp,sa,sd); fflush(stdout);} }
          if(nca%40==0){ int cur=0,best=0; for(int kk=0;kk<W;kk++){ int r=(kk<run)?res[kk]:(kk==run?0:loesch(a+(u128)kk*d)); lut[kk]++; if(r){ lup[kk]++; cur++; if(cur>best)best=cur; } else cur=0; } if(best>ma) ma=best; }
        } }
      #pragma omp critical
      { for(int k=0;k<=W;k++){ hist[k]+=h[k]; histp[k]+=hp[k]; } for(int k=0;k<W;k++){ tested[k]+=te[k]; passed[k]+=pa[k]; unkn[k]+=un[k]; ut[k]+=lut[k]; up[k]+=lup[k]; } ncand+=nca; ntests+=nt; if(ma>maxany)maxany=ma; }
    }
    double el=omp_get_wtime()-t0; double cpus=el*omp_get_max_threads();
    #define RNG(arr,lo,hi) ({ long s_=0; for(int k_=lo;k_<=hi;k_++) s_+=arr[k_]; (double)s_; })
    double p09=RNG(passed,0,9)/RNG(tested,0,9), p19=RNG(passed,1,9)/RNG(tested,1,9), p1029=RNG(tested,10,29)>0?RNG(passed,10,29)/RNG(tested,10,29):-1;
    double q09=(RNG(passed,0,9)-RNG(unkn,0,9))/RNG(tested,0,9), q19=(RNG(passed,1,9)-RNG(unkn,1,9))/RNG(tested,1,9), q1029=RNG(tested,10,29)>0?(RNG(passed,10,29)-RNG(unkn,10,29))/RNG(tested,10,29):-1;
    double u19=RNG(ut,1,9)>0?RNG(up,1,9)/RNG(ut,1,9):-1, u1029=RNG(ut,10,29)>0?RNG(up,10,29)/RNG(ut,10,29):-1;
    long ge[W+2]; ge[W+1]=0; for(int k=W;k>=0;k--) ge[k]=ge[k+1]+hist[k]; int mx=0; for(int k=0;k<=W;k++) if(hist[k]) mx=k;
    double M=1e6/(double)(ncand?ncand:1);
    // extrapolation: P(run>=n) = prod of measured hazards for k<10 individually, then pooled p10-29 (or p1-9 if none)
    double pp=p1029>0&&RNG(tested,10,29)>200?p1029:p19; double s10=1; for(int k=0;k<10;k++) s10*= tested[k]? (double)passed[k]/tested[k]:pp;
    double e20=s10,e25,e30; for(int k=10;k<20;k++) e20*=pp; e25=e20; for(int k=0;k<5;k++) e25*=pp; e30=e25; for(int k=0;k<5;k++) e30*=pp;
    printf("%-34s %3d %9ld %7.0f | %.3f/%.3f %.3f/%.3f %.3f/%.3f | %.3f  %.3f  | %8.0f %8.1f %8.2f %8.2f | %8.3g %8.3g %8.3g | %4d %4d | %6.2f\n",c->name,digits(c->type==7||c->type==1||c->type==6?c->d*500:c->d),ncand,ncand/cpus,p09,q09,p19,q19,p1029,q1029,u19,u1029,ge[5]*M,ge[10]*M,ge[15]*M,ge[20]*M,e20*1e6,e25*1e6,e30*1e6,mx,maxany,100.0*RNG(unkn,0,W-1)/(double)(ntests?ntests:1));
    fflush(stdout);
  }
  return 0; }

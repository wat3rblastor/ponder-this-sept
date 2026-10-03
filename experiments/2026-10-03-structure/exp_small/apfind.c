// apfind: bitmap of Loeschian numbers == 1 mod 6 up to N. Modes:
//  find N n D0 T      : all maximal runs of length>=n, d multiple of D0 (D0 multiple of 6), terms<=N
//  hist N D T Kmax Lcap excl : survival histogram for d=D*K, K<=Kmax, gcd(K,excl)=1, a in [N/4,N/2], a==1 mod 6, gcd(a,D)=1
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
typedef uint64_t u64;
static u64 *B; static u64 NB;
static inline int bit(u64 i){ return (B[i>>6]>>(i&63))&1; }
static u64 gcd(u64 a,u64 b){ while(b){u64 t=a%b;a=b;b=t;} return a; }
int main(int argc,char**argv){
  u64 N=strtoull(argv[2],0,10);
  NB=N/6+2; B=calloc(NB/64+4,8);
  u64 xm=1; while((xm+1)*(xm+1)<=N) xm++;
  #pragma omp parallel for schedule(dynamic,64)
  for(u64 x=1;x<=xm;x++){
    u64 xx=x*x;
    for(u64 y=0;y<=x;y++){ u64 v=xx+x*y+y*y; if(v>N)break; if(v%6==1) __atomic_fetch_or(&B[((v-1)/6)>>6],1ULL<<(((v-1)/6)&63),__ATOMIC_RELAXED); }
  }
  fprintf(stderr,"sieve done\n");
  if(!strcmp(argv[1],"find")){
    int n=atoi(argv[3]); u64 D0=strtoull(argv[4],0,10); u64 s0=D0/6;
    u64 Kmax=(N-1)/((u64)(n-1)*D0);
    #pragma omp parallel for schedule(dynamic,4)
    for(u64 K=1;K<=Kmax;K++){
      u64 s=s0*K; u64 imax=(N-1)/6; if((u64)(n-1)*s>imax) continue;
      u64 lim=imax-(u64)(n-1)*s;
      for(u64 w=0;w*64<=lim;w++){
        u64 m=B[w]; if(!m) continue;
        for(int k=n-1;k>=1 && m;k--){ u64 o=w*64+(u64)k*s; u64 q=o>>6; int r=o&63; u64 v=B[q]>>r; if(r) v|=B[q+1]<<(64-r); m&=v; }
        while(m){ int b=__builtin_ctzll(m); m&=m-1; u64 i=w*64+b; if(i>lim) break;
        if(i>=s && bit(i-s)) continue;
        int L=n; while(i+(u64)L*s<=imax && bit(i+(u64)L*s)) L++;
        #pragma omp critical
        printf("%llu %llu %d\n",(unsigned long long)(6*i+1),(unsigned long long)(6*s),L);
        }
      }
    }
  } else {
    u64 D=strtoull(argv[3],0,10); u64 Kmax=strtoull(argv[4],0,10); int Lcap=atoi(argv[5]); u64 excl=strtoull(argv[6],0,10);
    u64 astep=argc>7?strtoull(argv[7],0,10):1;
    u64 hist[128]={0}; u64 cand=0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:hist[:128]) reduction(+:cand)
    for(u64 K=1;K<=Kmax;K++){
      if(gcd(K,excl)!=1) continue;
      u64 d=D*K, s=d/6;
      for(u64 i=N/24;i<N/12;i+=astep){
        if(gcd(6*i+1,D)!=1) continue;
        cand++;
        int L=0; while(L<Lcap && bit(i+(u64)L*s)) L++;
        hist[L]++;
      }
    }
    printf("D %llu cand %llu\n",(unsigned long long)D,(unsigned long long)cand);
    for(int L=0;L<=Lcap;L++) printf("%d %llu\n",L,(unsigned long long)hist[L]);
  }
  return 0;
}

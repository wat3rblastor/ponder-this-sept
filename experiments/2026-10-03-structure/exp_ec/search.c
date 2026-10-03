// find coprime (X,m) with X^2 + 8 e m^2 square for >=2 values e in 1..R
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
typedef unsigned long long u64;
static inline int issq(u64 v, u64 *r){ static const int ok[64]={1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0};
  if(!ok[v&63]) return 0; u64 s=(u64)sqrtl((long double)v); while(s*s>v)s--; while((s+1)*(s+1)<=v)s++; *r=s; return s*s==v;}
static u64 gcd(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
int main(int argc,char**argv){ long B=atol(argv[1]); int R=atoi(argv[2]); int MINC=argc>3?atoi(argv[3]):2;
  omp_set_num_threads(16);
  #pragma omp parallel for schedule(dynamic,8)
  for(long m=1;m<=B;m++){ char buf[4096];
    for(long X=1;X<=B;X++){ if(gcd(X,m)!=1) continue; u64 x2=(u64)X*X, m8=8ULL*m*m; int c=0; int es[128]; u64 r;
      u64 v=x2; for(int e=1;e<=R;e++){ v+=m8; if(issq(v,&r)) es[c++]=e; }
      if(c>=MINC){ int n=sprintf(buf,"%ld %ld",X,m); for(int i=0;i<c;i++) n+=sprintf(buf+n," %d",es[i]); 
        #pragma omp critical
        puts(buf); }
    }}
  return 0;}

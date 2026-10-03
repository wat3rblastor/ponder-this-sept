/* Exhaustive scan for  max #{k in [0,57] : P(k) is a perfect square}
 * over quadratics P(k) = -n k^2 + beta k + gamma  (n >= 1), in the RIGHT
 * coordinates.
 *
 * Completing the square: with  X_k = beta - 2 n k  and  y_k^2 = P(k),
 *      X_k^2 + 4 n y_k^2 = Delta,    Delta := beta^2 + 4 n gamma.
 * So a configuration is exactly: an integer n >= 1, an integer Delta, and a
 * choice of 58 consecutive terms of the arithmetic progression X = beta - 2nk;
 * the hits are the representations Delta = X^2 + 4 n y^2 whose X lies in that
 * window (automatically X = beta mod 2n).
 *
 * This enumerates (n, Delta) directly, so it is complete for the box
 * 1 <= n <= NMAX, 0 < Delta <= DMAX -- which contains every configuration
 * whose hits satisfy  max_k 4 n y_k^2 + X_k^2 <= DMAX.  The previous
 * experiment looped over beta and so could not reach the vertex-centred
 * configurations for n >~ 15 ( beta ~ 57 n ).
 *
 * usage: ./scan NMAX DMAX REPORT
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static long isqrtl(long n){ if(n<0) return -1; long r=(long)sqrtl((long double)n);
    while(r>0 && r*r>n) r--; while((r+1)*(r+1)<=n) r++; return r; }

int main(int argc, char **argv){
    long NMAX = atol(argv[1]);
    long DMAX = atol(argv[2]);
    int  REPORT = atoi(argv[3]);
    unsigned char *cnt = calloc(DMAX+1,1);
    if(!cnt){fprintf(stderr,"oom\n");return 1;}
    long *X = malloc(sizeof(long)*4096);
    int best=0;
    for(long n=1;n<=NMAX;n++){
        long ymax = isqrtl(DMAX/(4*n));
        /* pass 1: tally number of (X,y) representations per Delta, recording
           the Deltas that reach REPORT (so pass 2 need not sweep the array) */
        long ncand=0, candcap=1<<16; long *cand=malloc(sizeof(long)*candcap);
        for(long y=0;y<=ymax;y++){
            long base = 4*n*y*y;
            long xm = isqrtl(DMAX-base);
            /* X=0 counts once, |X|>0 counts twice */
            if(cnt[base] < 250){ cnt[base]++;
                if(cnt[base]==REPORT){ if(ncand==candcap){candcap*=2;cand=realloc(cand,sizeof(long)*candcap);} cand[ncand++]=base; } }
            for(long x=1;x<=xm;x++){
                long D = base + x*x;
                if(cnt[D] < 249){ cnt[D]+=2;
                    if(cnt[D]>=REPORT && cnt[D]<REPORT+2){ if(ncand==candcap){candcap*=2;cand=realloc(cand,sizeof(long)*candcap);} cand[ncand++]=D; } }
            }
        }
        /* clear by replaying the same loop */
        for(long y=0;y<=ymax;y++){
            long base = 4*n*y*y; long xm = isqrtl(DMAX-base);
            cnt[base]=0;
            for(long x=1;x<=xm;x++) cnt[base+x*x]=0;
        }
        /* pass 2: for Delta with enough representations, do the window test */
        for(long ci=0;ci<ncand;ci++){
            long D=cand[ci];
            if(D<1) continue;
            /* collect all X with (D - X^2)/(4n) a perfect square */
            long nx=0; long xm = isqrtl(D);
            for(long x=0;x<=xm;x++){
                long r = D - x*x;
                if(r % (4*n)) continue;
                long q = r/(4*n); long s=isqrtl(q);
                if(s*s!=q) continue;
                if(nx<4090){ X[nx++]=x; if(x) X[nx++]=-x; }
            }
            if(nx<REPORT) continue;
            /* group by residue of X mod 2n, window = 58 consecutive AP terms */
            /* sort */
            for(long i=1;i<nx;i++){long v=X[i],j=i-1;while(j>=0&&X[j]>v){X[j+1]=X[j];j--;}X[j+1]=v;}
            long mod = 2*n;
            for(long i=0;i<nx;i++){
                /* window anchored so that X[i] is the k=0 end (beta = X[i]) */
                long beta = X[i]; int c=0; long maxy=0;
                for(long j=i;j>=0;j--){
                    long diff = beta - X[j];
                    if(diff < 0) continue;
                    if(diff % mod) continue;
                    long k = diff/mod;
                    if(k>57) break;
                    c++;
                    long q=(D-X[j]*X[j])/(4*n); long s=isqrtl(q); if(s>maxy) maxy=s;
                }
                if(c>best) best=c;
                if(c>=REPORT){
                    /* gamma = (Delta - beta^2)/(4n) */
                    long num = D - beta*beta;
                    if(num % (4*n)) continue;
                    printf("hits=%d n=%ld beta=%ld gamma=%ld Delta=%ld maxy=%ld\n",
                           c,n,beta,num/(4*n),D,maxy);
                    fflush(stdout);
                }
            }
        }
        free(cand);
        if(n%25==0) fprintf(stderr,"n=%ld done best=%d\n",n,best);
    }
    printf("BEST=%d  (box: 1<=n<=%ld, Delta<=%ld)\n",best,NMAX,DMAX);
    return 0;
}

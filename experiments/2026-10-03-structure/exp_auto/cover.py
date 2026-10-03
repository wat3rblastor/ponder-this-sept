# Which index sets are automatic?  t_k = 3x^2 + m^2(24k+1), d=24m^2.
# k automatic if k = i + p, x^2+8 i m^2 = square (i in A), 24p+1 = square (p gen. pentagonal).
# Pure combinatorics first: best A (|A|=2,3,4 incl 0) for coverage of [0,57].
import itertools, math
def issq(n): 
    if n<0: return False
    r=math.isqrt(n); return r*r==n
P=[p for p in range(0,400) if issq(24*p+1)]
def cover(A): return {i+p for i in A for p in P if 0<=i+p<=57}
print("P",[p for p in P if p<=57],len([p for p in P if p<=57]))
R=range(-150,58)
best=sorted(((len(cover((0,i))),i) for i in R if i),reverse=True)[:8]; print("|A|=2",best)
b3=sorted(((len(cover((0,i,j))),i,j) for i in R for j in R if i<j and i and j),reverse=True)[:8]; print("|A|=3",b3)
import random
top=[(c,i,j) for c,i,j in sorted(((len(cover((0,i,j))),i,j) for i in R for j in R if i<j and i and j),reverse=True)[:300]]
b4=sorted(((len(cover((0,i,j,l))),i,j,l) for c,i,j in top for l in R if l not in(0,i,j)),reverse=True)[:5]; print("|A|=4 (greedy)",b4)

from sympy import factorint
import random
def lo(n):
    return n>=0 and all(e%2==0 for p,e in factorint(n).items() if p%3==2)
random.seed(1)
S1=set();S2=set();first=True
for _ in range(40):
    u=random.randrange(1,4000,2); v=random.randrange(1,4000,2)
    A=(3*u*u+v*v)//2; a=A*A; d=12*u*u*v*v
    ok={k for k in range(58) if lo(a+k*d)}
    S2 = ok if first else S2&ok
    A=random.randrange(1,10**6); B=random.randrange(1,10**6)
    ok={k for k in range(58) if lo(A*A+3*k*B*B)}
    S1 = ok if first else S1&ok
    first=False
print("1-anchor always-Loeschian k:",sorted(S1))
print("2-anchor always-Loeschian k:",sorted(S2))

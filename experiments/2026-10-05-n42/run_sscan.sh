#!/bin/sh
# The s-scan: exact counts in "omit-one-bad-prime" sub-families.
#
# A family d = B*m is exactly countable whenever base(n) | B.  Choosing B to
# contain every forced bad prime EXCEPT one, q_min, puts the family's
# structural parameter  s = n / (2 q_min)  anywhere in (0,1] at a MODERATE n,
# because s is set by q_min, not by n.  That is what makes the regime
# s = 0.65..0.97 -- which n=58 (s=0.707) actually sits in, and which no
# previous exact count has reached with good statistics -- affordable here:
# cost is  ~7 X^2 / (64 B (n-1))  word ops, inverse in B.
#
#   q_min=5   B = 2*3*11*17*23*29 = 748374   n<=9    s = n/10
#   q_min=11  B = 2*3*5*17*23*29  = 340170   n<=21   s = n/22
#   q_min=17  B = 2*3*5*11*23*29  = 220110   n<=33   s = n/34
#   q_min=23  B = 2*3*5*11*17*29  = 162690   n<=45   s = n/46
#   q_min=29  B = 129030                     n<=57   s = n/58
#   q_min=41  B = 3741870 = base(58)         n<=81   s = n/82   <- the real one
#
# n > 2*q_min is empty by the forced-divisibility theorem, so s <= 1 always.
cd "$(dirname "$0")"
X=${X:-2000000000}
T=${T:-3}
export OMP_NUM_THREADS=$T
N=${N:-18}
set -x
nice -n $N ./exact $X 1 1 1 --dens --base 748374  > dens_748374.txt  2>> log.txt
nice -n $N ./exact $X 1 1 1 --dens --base 340170  > dens_340170.txt  2>> log.txt
nice -n $N ./exact $X 1 1 1 --dens --base 220110  > dens_220110.txt  2>> log.txt
nice -n $N ./exact $X 1 1 1 --dens --base 162690  > dens_162690.txt  2>> log.txt
nice -n $N ./exact $X 1 1 1 --dens --base 129030  > dens_129030.txt  2>> log.txt
nice -n $N ./exact $X 1 1 1 --dens --base 3741870 > dens_3741870.txt 2>> log.txt
nice -n $N ./exact $X  4  9 1 --base 748374  > q5.txt   2>> log.txt
nice -n $N ./exact $X 10 21 1 --base 340170  > q11.txt  2>> log.txt
nice -n $N ./exact $X 16 33 1 --base 220110  > q17.txt  2>> log.txt
nice -n $N ./exact $X 24 45 1 --base 162690  > q23.txt  2>> log.txt
nice -n $N ./exact $X 24 42 1 --base 3741870 > q41.txt  2>> log.txt
nice -n $N ./exact $X 24 40 2 --base 129030  > q29.txt  2>> log.txt

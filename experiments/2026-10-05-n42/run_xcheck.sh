#!/bin/sh
# X-independence check: re-measure the two families that straddle s = 0.707
# with the most relevant q_min at a SECOND, independent X.  If the residual is
# a function of s alone it must not move with X.
cd "$(dirname "$0")"
X=${X:-2000000000}
export OMP_NUM_THREADS=${T:-3}
N=${N:-18}
set -x
nice -n $N ./exact $X 1 1 1 --dens --base 220110 > dens_220110_2e9.txt 2>> log2.txt
nice -n $N ./exact $X 1 1 1 --dens --base 162690 > dens_162690_2e9.txt 2>> log2.txt
nice -n $N ./exact $X 20 33 1 --base 220110 > q17b.txt 2>> log2.txt
nice -n $N ./exact $X 28 45 1 --base 162690 > q23b.txt 2>> log2.txt

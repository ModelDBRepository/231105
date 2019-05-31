#!/bin/sh
# Plotting the f-g_el curve
#

fla=a4
flb=a41
fma=b2
dir=$HOME/shares/docu/work_ms/tlcx/cvsync/

#../../simulation_program/tc.ex $fla
#../../simulation_program/tc.ex $flb
#../../simulation_program/tc.ex $fma

avfcc.py './' $fla 1.0 -1.0 -1.0 -1.0
avfcc.py './' $flb 1.0 -1.0 -1.0 -1.0

python ${dir}/raspl.py $fma

cat > line.xx <<EOF
2750.0 0.5
2800.0 0.5
EOF

xmgrace -graph 0 avfcc.sfr.$fla avfcc.sfr.$flb \
	-graph 1 avfcc.sch.$fla avfcc.sch.$flb \
        -graph 2 tc.ras.$fma.pl.P \
	-graph 3 line.xx \
	-hdevice EPS -p fel2.gr -printfile fel2.eps

/bin/rm fr.$fla.xx chi.$fla.xx
/bin/rm fr.$flb.xx chi.$flb.xx
/bin/rm tc.ras.$fma.pl.E tc.ras.$fma.pl.P line.xx

# cv.$fla.xx 
# cv.$flb.xx 

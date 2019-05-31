#!/bin/sh
# Plotting Fig. 1.

fla=c1
flb=c2
fld=e11

fma=b1
fmb=b2
fmc=b3
fmd=d11

#../../simulation_program/tc.ex $fla
#../../simulation_program/tc.ex $flb
#../../simulation_program/tc.ex $fld
#../../simulation_program/tc.ex $fma
#../../simulation_program/tc.ex $fmb
#../../simulation_program/tc.ex $fmc
#../../simulation_program/tc.ex $fmd

avfcc.py './' $fla 1.0 2.2 3.0 3.0
avfcc.py './' $flb 1.0 2.2 3.0 3.0
avfcc.py './' $fld 1.0 2.2 3.0 3.0

avfcc.py './' $fma 1.0 2.2 3.0 3.0
avfcc.py './' $fmb 1.0 2.2 3.0 3.0
avfcc.py './' $fmc 1.0 7.59 7.69 7.69
avfcc.py './' $fmd 1.0 7.59 7.69 7.69

cat >> chi.independent.xx <<EOF
2.4 0.08
42.0 0.08
EOF

cat > line.bs.xx <<EOF
0.0 -6.1
30.0 55.3
EOF

cat > arrow_atc.xx <<EOF
7.60178 0.0
7.60178 1.0

7.60178 1.0
EOF

cat > arrow_att.xx <<EOF
4.6 0.0
4.6 1.0

4.6 1.0
EOF

xmgrace -graph 0 avfcc.sfr.$fla avfcc.sfr.$flb \
	         f_ar.dat avfcc.sfr.$fld \
	-graph 1 avfcc.scv.$fla avfcc.scv.$flb avfcc.scv.$fld \
	-graph 2 avfcc.sch.$fla avfcc.sch.$flb avfcc.sch.$fld \
	-graph 3 avfcc.sfr.$fma avfcc.sfr.$fmb avfcc.sfr.$fmc \
	         line.bs.xx avfcc.sfr.$fmd \
	-graph 4 avfcc.scv.$fma avfcc.scv.$fmb avfcc.scv.$fmc avfcc.scv.$fmd \
	-graph 5 avfcc.sch.$fma avfcc.sch.$fmb avfcc.sch.$fmc avfcc.sch.$fmd \
	         chi.independent.xx \
        -graph 6 arrow_atc.xx \
        -graph 7 arrow_att.xx \
		 -hdevice EPS -p fatd.gr -printfile fatd.eps

/bin/rm chi.independent.xx line.bs.xx
/bin/rm arrow_atc.xx arrow_att.xx
/bin/rm avfcc.???.*


#!/bin/sh
# Plotting the f-g_el curve
#

fma=b2b

raspl.py $fma

cat > line.xx <<EOF
2100.0 0.5
2200.0 0.5
EOF

xmgrace -graph 0 tc.ras.$fma.pl.P \
	-graph 1 tc.ras.$fma.pl.T \
	-graph 2 line.xx \
	-hdevice EPS -p rast.gr -printfile rast.${fma}.eps

/bin/rm tc.ras.$fma.pl.E tc.ras.$fma.pl.P line.xx

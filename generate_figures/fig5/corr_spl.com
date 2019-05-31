#!/bin/sh
# Plotting cross-correlations.

fl=b2
	      
xmgrace -graph 0 crc.cor.$fl \
        -graph 1 crc.cor.$fl \
        -hdevice EPS -p corr_spl.gr -printfile corr_spl.$fl.eps


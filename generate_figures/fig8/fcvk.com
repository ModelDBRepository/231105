#!/bin/sh
# Plotting the f-CVKin_fact curve
#

fla=a3

avfcc.py $fla 1.0 -1.0 -1.0 -1.0

xmgrace -graph 0 -settype xydy avfcc.snf.P.$fla \
        -graph 1 -settype xydy avfcc.sfr.P.$fla \
        -graph 2 -settype xydy avfcc.sfs.P.$fla \
	-graph 3 -settype xydy avfcc.snf.E.$fla \
        -graph 4 -settype xydy avfcc.sfr.E.$fla \
        -graph 5 -settype xydy avfcc.sfs.E.$fla \
        -hdevice EPS -p fcvka.gr -printfile fcvka.eps

xmgrace -graph 0  -settype xydy avfcc.scv.P.$fla \
	-graph 1  -settype xydy avfcc.sch.P.$fla \
	-graph 2  -settype xydy avfcc.scv.E.$fla \
	-graph 3  -settype xydy avfcc.sch.E.$fla \
        -hdevice EPS -p fcvkb.gr -printfile fcvkb.eps

/bin/rm avfcc.???.E.$fla avfcc.???.P.$fla

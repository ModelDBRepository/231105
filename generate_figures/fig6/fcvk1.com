#!/bin/sh
# Plotting the f-CVKin_fact curve
#

fla=a3

#../../simulation_program/tc.ex a3

avfcc.py $fla 1.0 -1.0 -1.0 -1.0

awk '{print $1, $9}' tc.avr.$fla > fr.$fla.xx
#awk '{print $1, $9}' tc.avr.$flb > fr.$flb.xx

awk '{print $1, $10}' tc.avr.$fla > cv.$fla.xx
#awk '{print $1, $10}' tc.avr.$flb > cv.$flb.xx

awk '{print $1, $14}' tc.avr.$fla > chi.$fla.xx
#awk '{print $1, $14}' tc.avr.$flb > chi.$flb.xx

awk '{print $1, $17}' tc.avr.$fla > nofr.$fla.xx

awk '{print $1, $18}' tc.avr.$fla > frcv.$fla.xx

xmgrace -graph 0 -settype xydy avfcc.snf.$fla \
        -graph 1 -settype xydy avfcc.sfr.$fla \
        -graph 2 -settype xydy avfcc.sfs.$fla \
        -hdevice EPS -p fcvk1a.gr -printfile fnc1a.eps

xmgrace -graph 0  -settype xydy avfcc.scv.$fla \
	-graph 1  -settype xydy avfcc.sch.$fla \
        -hdevice EPS -p fcvk1b.gr -printfile fnc1b.eps

/bin/rm fr.$fla.xx cv.$fla.xx chi.$fla.xx nofr.$fla.xx frcv.$fla.xx
#/bin/rm fr.$flb.xx cv.$flb.xx chi.$flb.xx

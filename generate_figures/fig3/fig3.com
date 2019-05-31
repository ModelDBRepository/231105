#!/bin/sh
# Plotting Fig. 3.
#

fgr=fcc_e

fla=p001
flb=p051
flc=p101
fld=a001
fle=a051
flf=a101
flg=q001
flh=q051
fli=q101

#../../simulation_program/tc.ex $fla
#../../simulation_program/tc.ex $flb
#../../simulation_program/tc.ex $flc
#../../simulation_program/tc.ex $fld
#../../simulation_program/tc.ex $fleb
#../../simulation_program/tc.ex $flec
#../../simulation_program/tc.ex $flfb
#../../simulation_program/tc.ex $flfc
#../../simulation_program/tc.ex $flg
#../../simulation_program/tc.ex $flh
#../../simulation_program/tc.ex $fli

avfcc.py './' $fla 25.0 -1.0 -1.0 -1.0
avfcc.py './' $flb 25.0 -1.0 -1.0 -1.0
avfcc.py './' $flc 25.0 -1.0 -1.0 -1.0
avfcc.py './' $fld 25.0 -1.0 -1.0 -1.0
avfcc.py './' $fle 25.0 -1.0 -1.0 -1.0
#avfcc.py './' $flec 25.0 -1.0 -1.0 -1.0
avfcc.py './' $flf 25.0 -1.0 -1.0 -1.0
#avfcc.py './' $flfc 25.0 -1.0 -1.0 -1.0
avfcc.py './' $flg 25.0 -1.0 -1.0 -1.0
avfcc.py './' $flh 25.0 -1.0 -1.0 -1.0
avfcc.py './' $fli 25.0 -1.0 -1.0 -1.0

xmgrace -graph 0 avfcc.sfr.$fla avfcc.sfr.$flb avfcc.sfr.$flc \
	-graph 1 avfcc.sch.$fla avfcc.sch.$flb avfcc.sch.$flc \
	-graph 2 avfcc.sch.$fld avfcc.sch.$fle avfcc.sch.$flf \
	-graph 3 avfcc.sch.$flg avfcc.sch.$flh avfcc.sch.$fli \
	-graph 4 -nxy Data_gii_0.9_20_realiz.txt \
	-hdevice EPS -p ${fgr}.gr -printfile ${fgr}.eps

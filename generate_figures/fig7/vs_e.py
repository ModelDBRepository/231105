#!/usr/bin/python
# This program plots the firing rate nu_E vs. Isyn_E-I_syn_I .

import sys
import os
import math

def read_fr(frI, cvI, ffrIx):
    for line_read in ffrIx:
        line_list = line_read.split()
        if len(line_list) >= 1:
            if (int(line_list[2]) == 1):
                fr_val = float(line_list[1])
                cv_val = float(line_list[7])
                frE.append(fr_val)
                cvE.append(cv_val)
        
def read_synvar(synvarET, synvarEE, synvarEI, fzmpx):
    for line_read in fzmpx:
        line_list = line_read.split()
        if len(line_list) >= 1:
            if (int(line_list[0]) == 1):
                VsET_val = -float(line_list[6])
                synvarET.append(VsET_val)
                VsEE_val = -float(line_list[7])
                synvarEE.append(VsEE_val)
                VsEI_val = -float(line_list[8])
                synvarEI.append(VsEI_val)
        
# main
suffix = 'a2'
#'a1'

ffri = open('tc.fri.' + suffix, 'r')
fzmp = open('tc.zmp.' + suffix, 'r')
ffis = open('vs_e.fis.' + suffix, 'w')
ffcv = open('vs_e.fcv.' + suffix, 'w')
fout = open('vs_e.out.' + suffix, 'w')

Isnic = 0.7328
VsynE = 0.0
tsynE = 2.0
VsynI = -85.0
tsynI = 3.0

KET = 50.0
gET = 0.15 * 4

KEE = 200.0
gEE = 0.2 * 4.0

KEI = 25.0
gEI = 0.7 * 4

gsynET = gET / (math.sqrt(KET) * tsynE)
gsynEE = gEE / (math.sqrt(KEE) * tsynE)
gsynEI = gEI / (math.sqrt(KEI) * tsynI)
JET = gsynET
JEE = gsynEE
JEI = gsynEI

print 'JET=', JET, ' JEE=', JEE, 'JEI=', JEI,

frE = []
cvE = []
read_fr(frE, cvE, ffri)

synvarET = []
synvarEE = []
synvarEI = []
read_synvar(synvarET, synvarEE, synvarEI, fzmp)

non = len(frE)

for ion in range(0, non):
    Jsetx = JET * synvarET[ion]
    ffis.write('{0:g} {1:g}\n'.format(Jsetx, frE[ion]))

ffis.write(' \n')
for ion in range(0, non):
    Jseex = JEE * synvarEE[ion]
    ffis.write('{0:g} {1:g}\n'.format(Jseex, frE[ion]))

ffis.write(' \n')
for ion in range(0, non):
    Jseix = JEI * synvarEI[ion]
    ffis.write('{0:g} {1:g}\n'.format(Jseix, frE[ion]))

ffis.write(' \n')
for ion in range(0, non):
    Jsetx = JET * synvarET[ion]
    Jseex = JEE * synvarEE[ion]
    Jseix = JEI * synvarEI[ion]

    ffis.write('{0:g} {1:g} {2:g} {3:g} {4:g} {5:g} {6:g} {7:g} {8:d}\n'.\
    format(Jsetx + Jseex + Jseix, frE[ion], Jsetx, Jseex, Jseix, synvarET[ion],\
    synvarEE[ion], synvarEI[ion], ion+1))

    ffcv.write('{0:g} {1:g}\n'.format(frE[ion], cvE[ion]))
                                                                     
ffri.close()
fzmp.close()
ffis.close()
ffcv.close()
fout.close()

#!/usr/bin/python
# This program plots the firing rate nu_I vs. Isyn_E-I_syn_I .

import sys
import os
import math

def read_fr(frI, cvI, ffrIx):
    for line_read in ffrIx:
        line_list = line_read.split()
        if len(line_list) >= 1:
            if (int(line_list[2]) == 2):
                fr_val = float(line_list[1])
                cv_val = float(line_list[7])
                frI.append(fr_val)
                cvI.append(cv_val)
        
def read_synvar(synvarIT, synvarIE, synvarII, fzmpx):
    for line_read in fzmpx:
        line_list = line_read.split()
        if len(line_list) >= 1:
            if (int(line_list[0]) == 2):
                VsIT_val = -float(line_list[6])
                synvarIT.append(VsIT_val)
                VsIE_val = -float(line_list[7])
                synvarIE.append(VsIE_val)
                VsII_val = -float(line_list[8])
                synvarII.append(VsII_val)
        
# main
suffix = 'a2'
'a1'

ffri = open('tc.fri.' + suffix, 'r')
fzmp = open('tc.zmp.' + suffix, 'r')
ffis = open('vs.fis.' + suffix, 'w')
ffcv = open('vs.fcv.' + suffix, 'w')
fout = open('vs.out.' + suffix, 'w')

Isnic = 0.7328
VsynE = 0.0
tsynE = 2.0
VsynI = -85.0
tsynI = 3.0

KIT = 75.0
gIT = 0.2 * 4

KIE = 400.0
gIE = 0.6 * 4.0

KII = 25.0
gII = 0.5 * 4

gsynIT = gIT / (math.sqrt(KIT) * tsynE)
gsynIE = gIE / (math.sqrt(KIE) * tsynE)
gsynII = gII / (math.sqrt(KII) * tsynI)
JIT = gsynIT
JIE = gsynIE
JII = gsynII

print 'JIT=', JIT, ' JIE=', JIE, 'JII=', JII,

frI = []
cvI = []
read_fr(frI, cvI, ffri)

synvarIT = []
synvarIE = []
synvarII = []
read_synvar(synvarIT, synvarIE, synvarII, fzmp)

non = len(frI)

for ion in range(0, non):
    Jsitx = JIT * synvarIT[ion]
    ffis.write('{0:g} {1:g}\n'.format(Jsitx, frI[ion]))

ffis.write(' \n')
for ion in range(0, non):
    Jsiex = JIE * synvarIE[ion]
    ffis.write('{0:g} {1:g}\n'.format(Jsiex, frI[ion]))

ffis.write(' \n')
for ion in range(0, non):
    Jsiix = JII * synvarII[ion]
    ffis.write('{0:g} {1:g}\n'.format(Jsiix, frI[ion]))

ffis.write(' \n')
for ion in range(0, non):
    Jsitx = JIT * synvarIT[ion]
    Jsiex = JIE * synvarIE[ion]
    Jsiix = JII * synvarII[ion]

    ffis.write('{0:g} {1:g} {2:g} {3:g} {4:g} {5:g} {6:g} {7:g} {8:d}\n'.\
    format(Jsitx + Jsiex + Jsiix, frI[ion], Jsitx, Jsiex, Jsiix, synvarIT[ion],\
    synvarIE[ion], synvarII[ion], ion+1))

    ffcv.write('{0:g} {1:g}\n'.format(frI[ion], cvI[ion]))
                                                                     
ffri.close()
fzmp.close()
ffis.close()
ffcv.close()
fout.close()

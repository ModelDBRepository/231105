#!/usr/bin/python
# This program plots the firing rate nu_I vs. Isyn_E-I_syn_I .

import sys
import os
import math

def read_fr(frI, cvI, ffrIx):
    for line_read in ffrIx:
        line_list = line_read.split()
        fr_val = float(line_list[1])
        cv_val = float(line_list[7])
        frI.append(fr_val)
        cvI.append(cv_val)
        
def read_synvar(synvarIT, synvarII, fzmpx):
    for line_read in fzmpx:
        line_list = line_read.split()
        VsIT_val = -float(line_list[5])
        synvarIT.append(VsIT_val)
        VsII_val = -float(line_list[6])
        synvarII.append(VsII_val)
        
# main
#suffix = str(sys.argv[1]);
#print 'suffix=',suffix
suffix = 'b2b'
#suffix = 'b2e'
#'b2ra'

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

KII = 25.0
gII = 0.55 * 4

gsynIT =  gIT / (math.sqrt(KIT) * tsynE)
gsynII = gII / (math.sqrt(KII) * tsynI)
JIT = gsynIT
JII = gsynII

print 'gsynIT=', gsynIT, ' JIT=', JIT, 'gsynII=', gsynII, ' JII=', JII

frI = []
cvI = []
read_fr(frI, cvI, ffri)

synvarIT = []
synvarII = []
read_synvar(synvarIT, synvarII, fzmp)

non = len(frI)
for ion in range(0, non):
    Jsitx = JIT * synvarIT[ion]
    Jsiix = JII * synvarII[ion]

    ffis.write('{0:g} {1:g} {2:g} {3:g} {4:g} {5:g} {6:d}\n'.format( \
    Jsitx + Jsiix, frI[ion], Jsitx, Jsiix, synvarIT[ion], synvarII[ion], ion+1))

    ffcv.write('{0:g} {1:g}\n'.format(frI[ion], cvI[ion]))
                                                                     
ffis.write(' \n')
for ion in range(0, non):
    Jsitx = JIT * synvarIT[ion]
    ffis.write('{0:g} {1:g}\n'.format(Jsitx, frI[ion]))

ffis.write(' \n')
for ion in range(0, non):
    Jsiix = JII * synvarII[ion]
    ffis.write('{0:g} {1:g}\n'.format(Jsiix, frI[ion]))


ffri.close()
fzmp.close()
ffis.close()
ffcv.close()
fout.close()

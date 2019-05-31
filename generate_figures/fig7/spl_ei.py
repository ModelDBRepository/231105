#!/usr/bin/python
# This program plots Isyn_E and I_syn_I .

import sys
import os
import math

#suffix = str(sys.argv[1]);
#print 'suffix=',suffix
suffix = 'a2'
#'a1'

Isnic = 0.7328
VsynE = 0.0
DelVE = 65.0
tsynE = 2.0
VsynI = -85.0
tsynI = 3.0

KIT = 75.0
gIT = 0.2 * 4.0

KIE = 400.0
gIE = 0.6 * 4.0

KII = 25.0
gII = 0.5 * 4.0

gsynIT = gIT / (math.sqrt(KIT) * tsynE)
gsynIE = gIE / (math.sqrt(KIE) * tsynE)
gsynII = gII / (math.sqrt(KII) * tsynI)
JIT = gsynIT
JIE = gsynIE
JII = gsynII

print 'JIT=', JIT, ' JIE=', JIE, 'JII=', JII

print 'gsynIT=', gsynIT, 'gsynIE=', gsynIE, 'gsynII=', gsynII

fcol = open('tc.col.' + suffix, 'r')
fcol1 = open('tc.col.' + suffix + '1', 'r')
fvv = open('vv.' + suffix + '.xx', 'w')
fsit = open('sit.' + suffix + '.xx', 'w')
fsie = open('sie.' + suffix + '.xx', 'w')
fsii = open('sii.' + suffix + '.xx', 'w')
fsall = open('sall.' + suffix + '.xx', 'w')
fline = open('line.' + suffix + '.xx', 'w')
fras = open('tc.ras.' + suffix, 'r')
fspk = open('spk.' + suffix + '.xx', 'w')
fout = open('spl.out.' + suffix, 'w')

tbeg = 500.0002
ntstat = 0
sumIIT = 0.0
sumIIE = 0.0
sumIII = 0.0
sumJsitx = 0.0
sumJsiex = 0.0
sumJsiix = 0.0

for line in fcol1:
    line_list = line.split()
    time = float(line_list[0])
    VV = float(line_list[37])
    fvv.write('{0:g} {1:g}\n'.format(time, VV))

for line in fcol:
    line_list = line.split()
    time = float(line_list[0])
    VV = float(line_list[37])
    sit = float(line_list[38])
    sie = float(line_list[39])
    sii = float(line_list[40])
    Jsitx = -JIT * sit * (VV - VsynE)
    Jsiex = -JIE * sie * (VV - VsynE)
    Jsiix = -JII * sii * (VV - VsynI)
    fsit.write('{0:g} {1:g}\n'.format(time, Jsitx))
    fsie.write('{0:g} {1:g}\n'.format(time, Jsiex))
    fsii.write('{0:g} {1:g}\n'.format(time, Jsiix))
    fsall.write('{0:g} {1:g}\n'.format(time, Jsitx + Jsiex + Jsiix))

    if time > tbeg:
        ntstat += 1
        sumIIT += sit
        sumIIE += sie
        sumIII += sii
        sumJsitx += Jsitx
        sumJsiex += Jsiex
        sumJsiix += Jsiix
    
fline.write('{0:s} {1:g}\n'.format('0.0 ', Isnic))
fline.write('{0:s} {1:g}\n'.format('6000.0 ', Isnic))

sumIIT /= ntstat
sumIIE /= ntstat
sumIII /= ntstat
sumJsitx /= ntstat
sumJsiex /= ntstat
sumJsiix /= ntstat

fout.write('ntstat={0:d} sumIIT={1:g} sumIIE={2:g} sumIII={3:g}\n'.format( \
ntstat, sumIIT, sumIIE, sumIII))
fout.write('sumJsitx={0:g} sumJsiex={1:g} sumJsiix={2:g}\n'.format( \
sumJsitx, sumJsiex, sumJsiix))

for line in fras:
    line_list = line.split()
    if int(line_list[1]) == 1 and int(line_list[2]) == 1:
        tspk = float(line_list[0])
        if tspk >= 3700.0 and tspk <= 3000.0:
            fspk.write('{0:g} {1:g}\n'.format(tspk, 0.2))
            fspk.write('{0:g} {1:g}\n'.format(tspk, 0.8))
            fspk.write(' \n')

        
fcol.close()
fcol1.close()
fvv.close()
fsit.close()
fsie.close()
fsii.close()
fsall.close()
fline.close()

fras.close()
fspk.close()
fout.close()

flin = open('line.xx', 'w')
flin.write('2750.0 0.5\n')
flin.write('2800.0 0.5\n')
flin.close()

xm_com = 'xmgrace' + \
' -graph 0 spk.' + suffix + '.xx' + \
' -graph 1 sit.' + suffix + '.xx' + \
         ' sii.' + suffix + '.xx' + \
         ' sall.' + suffix + '.xx' + \
         ' line.' + suffix + '.xx' + \
' -graph 2 vv.' + suffix + '.xx' + \
' -graph 3 sit.' + suffix + '.xx' + \
         ' sie.' + suffix + '.xx' + \
         ' sii.' + suffix + '.xx' + \
' -graph 4 line.xx' + \
' -hdevice EPS -p spl_ei.gr -printfile spl_ei.' + suffix + '.eps'

os.system(xm_com)

xm_com = 'xmgrace' + \
' -graph 0 vs.fis.' + suffix + ' fi.dat' +  \
' -graph 1 vs.fcv.' + suffix + \
' -hdevice EPS -p spl_ei_dot.gr -printfile spl_ei_dot.' + suffix + '.eps'

os.system(xm_com)

del_com = '/bin/rm sit.' + suffix + '.xx' + ' sii.' + suffix + '.xx' + \
          ' sall.' + suffix + '.xx' + ' line.' + suffix + '.xx' \
          ' vv.' + suffix + '.xx' + ' spk.' + suffix + '.xx'

os.system(del_com)
os.system('/bin/rm line.xx')

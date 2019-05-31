#!/usr/bin/python
# This program plots Isyn_E and I_syn_I .

import sys
import os
import math

#suffix = str(sys.argv[1]);
#print 'suffix=',suffix
suffix = 'a2'
#'a1'

Isnic = 0.30261
VsynE = 0.0
DelVE = 65.0
tsynE = 2.0
VsynI = -85.0
tsynI = 3.0

KET = 50.0
gET = 0.15 * 4.0

KEE = 200.0
gEE = 0.2 * 4.0

KEI = 25.0
gEI = 0.7 * 4.0

gsynET = gET / (math.sqrt(KET) * tsynE)
gsynEE = gEE / (math.sqrt(KEE) * tsynE)
gsynEI = gEI / (math.sqrt(KEI) * tsynI)
JET = gsynET
JEE = gsynEE
JEI = gsynEI

print 'JET=', JET, ' JEE=', JEE, 'JEI=', JEI

print 'gsynET=', gsynET, 'gsynEE=', gsynEE, 'gsynEI=', gsynEI

fcol = open('tc.col.' + suffix, 'r')
fcol1 = open('tc.col.' + suffix + '1', 'r')
fvv = open('vv.' + suffix + '.xx', 'w')
fset = open('set.' + suffix + '.xx', 'w')
fsee = open('see.' + suffix + '.xx', 'w')
fsei = open('sei.' + suffix + '.xx', 'w')
fsall = open('sall.' + suffix + '.xx', 'w')
fline = open('line.' + suffix + '.xx', 'w')
fras = open('tc.ras.' + suffix, 'r')
fspk = open('spk.' + suffix + '.xx', 'w')
fout = open('spl.out.' + suffix, 'w')

tbeg = 500.0002
ntstat = 0
sumIET = 0.0
sumIEE = 0.0
sumIEI = 0.0
sumJsetx = 0.0
sumJseex = 0.0
sumJseix = 0.0

for line in fcol1:
    line_list = line.split()
    time = float(line_list[0])
    VV = float(line_list[3])
    fvv.write('{0:g} {1:g}\n'.format(time, VV))

for line in fcol:
    line_list = line.split()
    time = float(line_list[0])
    VV = float(line_list[3])
    set = float(line_list[4])
    see = float(line_list[5])
    sei = float(line_list[6])
    Jsetx = -JET * set * (VV - VsynE)
    Jseex = -JEE * see * (VV - VsynE)
    Jseix = -JEI * sei * (VV - VsynI)
    fset.write('{0:g} {1:g}\n'.format(time, Jsetx))
    fsee.write('{0:g} {1:g}\n'.format(time, Jseex))
    fsei.write('{0:g} {1:g}\n'.format(time, Jseix))
    fsall.write('{0:g} {1:g}\n'.format(time, Jsetx + Jseex + Jseix))

    if time > tbeg:
        ntstat += 1
        sumIET += set
        sumIEE += see
        sumIEI += sei
        sumJsetx += Jsetx
        sumJseex += Jseex
        sumJseix += Jseix
    
fline.write('{0:s} {1:g}\n'.format('0.0 ', Isnic))
fline.write('{0:s} {1:g}\n'.format('6000.0 ', Isnic))

sumIET /= ntstat
sumIEE /= ntstat
sumIEI /= ntstat
sumJsetx /= ntstat
sumJseex /= ntstat
sumJseix /= ntstat

fout.write('ntstat={0:d} sumIET={1:g} sumIEE={2:g} sumIEI={3:g}\n'.format( \
ntstat, sumIET, sumIEE, sumIEI))
fout.write('sumJsetx={0:g} sumJseex={1:g} sumJseix={2:g}\n'.format( \
sumJsetx, sumJseex, sumJseix))

for line in fras:
    line_list = line.split()
    if int(line_list[1]) == 1 and int(line_list[2]) == 1:
        tspk = float(line_list[0])
        if tspk >= 3700.0 and tspk <= 3000.0:
            fspk.write('{0:g} {1:g}\n'.format(tspk, 0.2))
            fspk.write('{0:g} {1:g}\n'.format(tspk, 0.8))
            fspk.write(' \n')

        
fcol.close()
fvv.close()
fset.close()
fsee.close()
fsei.close()
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
' -graph 1 set.' + suffix + '.xx' + \
         ' sei.' + suffix + '.xx' + \
         ' sall.' + suffix + '.xx' + \
         ' line.' + suffix + '.xx' + \
' -graph 2 vv.' + suffix + '.xx' + \
' -graph 3 set.' + suffix + '.xx' + \
         ' see.' + suffix + '.xx' + \
         ' sei.' + suffix + '.xx' + \
' -graph 4 line.xx' + \
' -hdevice EPS -p spl_ei_e.gr -printfile spl_ei_e.' + suffix + '.eps'

os.system(xm_com)

xm_com = 'xmgrace' + \
' -graph 0 vs_e.fis.' + suffix + ' fi.e.dat' +  \
' -graph 1 vs_e.fcv.' + suffix + \
' -hdevice EPS -p spl_ei_dot_e.gr -printfile spl_ei_dot_e.' + suffix + '.eps'

os.system(xm_com)

del_com = '/bin/rm set.' + suffix + '.xx' + ' sei.' + suffix + '.xx' + \
          ' sall.' + suffix + '.xx' + ' line.' + suffix + '.xx' \
          ' vv.' + suffix + '.xx' + ' spk.' + suffix + '.xx'

os.system(del_com)
os.system('/bin/rm line.xx')

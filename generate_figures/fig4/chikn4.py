#!/usr/bin/python
#  chi vs. Meff
import sys
import math
import collections
import os
import numpy as np

# This function processes the data for one parameter set.
# def read_col(suffix_avr, KKI0, nonI,  GII, DelVI, KKT0, nonT, GIT, DelVT)
def read_col(suffix_avr, Tc, Ic, IT, II):
    fsfr = open('avfcc.sfr.' + suffix_avr, 'r')
    fscv = open('avfcc.scv.' + suffix_avr, 'r')
    fsch = open('avfcc.sch.' + suffix_avr, 'r')
# K
    fres = open('chi.res.' + suffix_avr + '.xx', 'w')
# K_eff
    feff = open('chi.eff.' + suffix_avr + '.xx', 'w')
# frI
    ffri = open('chi.fri.' + suffix_avr + '.xx', 'w')
# cvI
    fcvi = open('chi.cvi.' + suffix_avr + '.xx', 'w')

    col_read = [];
    for line in fsch:
        var_list = line.split()
        KKT = IT['KK0'] * float(var_list[0])
        KKI = II['KK0'] * float(var_list[0])
        sfr_line = fsfr.readline()
        sfr_var_list = sfr_line.split()
        chiI = float(var_list[1])
        frI = float(sfr_var_list[1])
        cvI = float(var_list[1])
        lone = []
        lone.append(KKI)
        lone.append(KKT)
        lone.append(chiI)
        lone.append(frI)
        lone.append(cvI)
        col_read.append(lone)
        
    for kk_chi in col_read:
        fres.write('{0:g} {1:g}\n'.format(kk_chi[0], kk_chi[2]))

    factorII = (II['GG'] * Ic['DelV'] * Ic['nu'])**2
    factorIT = (IT['GG'] * Tc['DelV'] * Tc['nu'])**2

    if Ic['non'] > epsilon:
        for kk_chi in col_read:
            KKI = kk_chi[0]
            KKT = kk_chi[1]
            if math.fabs(KKI - Ic['non']) > epsilon and KKI > epsilon:
                denomI =  factorII * KKI * ((1.0 / KKI) - (1.0 / Ic['non']))
                if math.fabs(KKT - Tc['non']) > epsilon and KKT > epsilon:
                    denomT = factorIT * KKT * ((1.0 / KKT) - (1.0 / Tc['non']))
                else:
                    denomT = 0
                numerI = factorII * KKI
                Keff1 = numerI / (denomI + denomT)
                Keff  = 1.0 / ((1.0 / KKI) - (1.0 / Ic['non']))
#                feff.write('{0:g} {1:g} {2:g} {3:g} {4:g} {5:g} {6:g} {7:g} {8:g}\n'.format(Keff, kk_chi[2], denomI, denomT, numerI, KKI, Ic['non'], KKT, Tc['non']))
                feff.write('{0:g} {1:g}\n'.format(Keff, kk_chi[2]))
                ffri.write('{0:g} {1:g}\n'.format(Keff, kk_chi[3]))
                fcvi.write('{0:g} {1:g}\n'.format(Keff, kk_chi[4]))

    fsfr.close()
    fscv.close()
    fsch.close()
    fres.close()
    feff.close()
    ffri.close()
    fcvi.close()

#---------------------------------------------------------------
#main

epsilon = 1.0e-10

Tc = dict()
Ic = dict()
IT = dict()
II = dict()

Tc['DelV'] = 65.0
Tc['nu']   = 16.0

Ic['DelV'] = -20.0
Ic['nu']   = 15.5

IT['GG']  = 0.3
IT['KK0'] = 75.0

II['GG']  = 0.85 
II['KK0'] = 25.0

#fl='ksp'

fl = ''
fla = fl + 'p001' #'1'
flb = fl + 'p051' #'2'
flc = fl + 'p101' #'3'
fld = fl + 'p002' #'4'
fle = fl + 'p052' #'5'
flf = fl + 'p102' #'6'
flg = fl + 'p003' #'7'
flh = fl + 'p053' #'8'
fli = fl + 'p103' #'9'

os.system('avfcc.py ./ ' + fla + ' 1.0 -1.0 -1.0 -1.0')
os.system('avfcc.py ./ ' + flb + ' 1.0 -1.0 -1.0 -1.0')
os.system('avfcc.py ./ ' + flc + ' 1.0 -1.0 -1.0 -1.0')
os.system('avfcc.py ./ ' + fle + ' 1.0 -1.0 -1.0 -1.0')
os.system('avfcc.py ./ ' + flg + ' 1.0 -1.0 -1.0 -1.0')
os.system('avfcc.py ./ ' + flh + ' 1.0 -1.0 -1.0 -1.0')
os.system('avfcc.py ./ ' + fli + ' 1.0 -1.0 -1.0 -1.0')

Tc['non'] =   200.0
Ic['non'] =   150.0
read_col(fla, Tc, Ic, IT, II)

Tc['non'] =   200.0
Ic['non'] =   150.0
read_col(flb, Tc, Ic, IT, II)

Tc['non'] =   200.0
Ic['non'] =   150.0
read_col(flc, Tc, Ic, IT, II)

Tc['non'] =  2000.0
Ic['non'] =  1500.0
#read_col(fld, Tc, Ic, IT, II)

Tc['non'] =  2000.0
Ic['non'] =  1500.0
read_col(fle, Tc, Ic, IT, II)

Tc['non'] =  2000.0
Ic['non'] =  1500.0
#read_col(flf, Tc, Ic, IT, II)

Tc['non'] = 20000.0
Ic['non'] = 15000.0
read_col(flg, Tc, Ic, IT, II)

Tc['non'] = 20000.0
Ic['non'] = 15000.0
read_col(flh, Tc, Ic, IT, II)

Tc['non'] = 20000.0
Ic['non'] = 15000.0
read_col(fli, Tc, Ic, IT, II)


xm_command = 'xmgrace -graph 0 chi.fri.' + flb + '.xx ' + \
                              'chi.fri.' + fle + '.xx  '+ \
                              'chi.fri.' + flh + '.xx  '+ \
                     '-graph 1 chi.eff.' + flb + '.xx ' + \
                              'chi.eff.' + fle + '.xx  '+ \
                              'chi.eff.' + flh + '.xx  '+ \
                     '-graph 2 chi.eff.' + flg + '.xx ' + \
                              'chi.eff.' + flh + '.xx  '+ \
                              'chi.eff.' + fli + '.xx  '+ \
                     '-hdevice EPS -p chikn4.gr ' + \
                     '-printfile chikn4.' + fla + '.eps'

os.system(xm_command)

os.system('/bin/rm chi.res.????.xx chi.eff.????.xx')
os.system('/bin/rm chi.fri.????.xx chi.cvi.????.xx')


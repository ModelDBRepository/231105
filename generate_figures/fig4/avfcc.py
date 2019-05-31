#!/usr/bin/python
# This program computes statistics and standard deviation of FR, Zmd and Zphi 
# values from tc.avr.?? .

import sys
import math
import collections
import numpy as np

# This function calculate the average and the standard deviation, given the
# sum, the sum of squares, and the number of data items.
def avr_sd_cal(sumvalZ, sumvalZsq, nparZ):
    if (nparZ <= 0):
        avvalZ = -999.9
        avvalZsq = -999.8
        diffZ = -999.7
        sigvalZ = -999.6
    else:
        avvalZ = sumvalZ / nparZ
        avvalZsq = sumvalZsq / nparZ
        diffZ = avvalZsq - avvalZ * avvalZ
        if (diffZ >= 0.0):
            sigvalZ = math.sqrt(diffZ)
        else:
            sigvalZ = -999.7

#    print 'nparZ=', nparZ, 'avvalZ=', avvalZ, 'sigvalZ=', sigvalZ

    return avvalZ, sigvalZ

# This function reads a specific column and writes the average and the
# standard deviation.
def one_column_calculate_statistics(flread, flwrite, col_num, xval_min):
    par_old = -999.999
    npar_all = 0
    nparZ = 0
    sumvalZ = 0.0
    sumvalZsq = 0.0
    
    flread.seek(0)

    for line in flread:
        val_list = line.split()
        par = float(val_list[0])
        valZ = float(val_list[col_num])
        npar_all += 1
    
        if par != par_old:
            if (npar_all > 1):
                avvalZ, sigvalZ = avr_sd_cal(sumvalZ, sumvalZsq, nparZ)
#               flwrite.write('{0:g} {1:g} {2:g} {3:d}\n'.format( \
#               par_old * factor, sumvalZ, sumvalZsq, nparZ)) 
                if (avvalZ >= 0.0 and par_old * factor >= xval_min + 1.0e-10):
                    flwrite.write('{0:g} {1:g} {2:g}\n'.format( \
                    par_old * factor, avvalZ, sigvalZ))
    
            nparZ = 1
            sumvalZ = valZ
            sumvalZsq = valZ * valZ
            par_old = par
        else:
            nparZ += 1
            sumvalZ += valZ
            sumvalZsq += valZ * valZ
    
    
    if (npar_all > 0):
        avvalZ, sigvalZ = avr_sd_cal(sumvalZ, sumvalZsq, nparZ)
#       flwrite.write('{0:g} {1:g} {2:g} {3:d}\n'.format(par_old * factor, sumvalZ, sumvalZsq, nparZ)
        if (avvalZ >= 0.0):            
            flwrite.write('{0:g} {1:g} {2:g}\n'.format(par_old * factor, avvalZ, sigvalZ))

# This function calculates the average and the correlation matrix, given a
# vector of sums, a matrix of sums of squares, and the number of data items.
def avr_CC_cal(sumZcs, sumZZcs, nparZ):
    avZcs = np.array([-999.9, -999.8])
    CCZcs = np.array([[-998.9, -998.8], [-998.7, -998.7]])

    if (nparZ > 0):
        avZcs[0] = sumZcs[0] / nparZ
        avZcs[1] = sumZcs[1] / nparZ
        CCZcs[0,0] = sumZZcs[0,0] / nparZ - avZcs[0] * avZcs[0]
        CCZcs[0,1] = sumZZcs[0,1] / nparZ - avZcs[0] * avZcs[1]
        CCZcs[1,0] = sumZZcs[1,0] / nparZ - avZcs[1] * avZcs[0]
        CCZcs[1,1] = sumZZcs[1,1] / nparZ - avZcs[1] * avZcs[1]

    return avZcs, CCZcs

# This function computes the correlation matrix of Zmd and Zphi given the
# correnation matrix of Zcos and Zsin. Then, it calculates the standard
# deviations of Zmd and Zphi.
def compute_error_bars_md_phi(avZcs, CCZcs, avZmd):
    avZmdsq = avZmd * avZmd
    Amat = np.array([[ avZcs[0]/avZmd,   avZcs[1]/avZmd], \
                     [-avZcs[1]/avZmdsq, avZcs[0]/avZmdsq]])
    AmatT = Amat.T
    CCZmp = (Amat.dot(CCZcs)).dot(AmatT)
#   print 'Amat=', Amat, '\nAmatT=', AmatT, '\nCCZcs=', CCZcs, '\nCCZmp=', CCZmp

    if (CCZmp[0,0] >= 0.0):
        sigmd =  math.sqrt(CCZmp[0,0])
    else:
        sigmd = -997.9

    if (CCZmp[1,1] >= 0.0):
        sigphi =  math.sqrt(CCZmp[1,1])
    else:
        sigphi = -997.8

    return sigmd, sigphi


# This function reads a specific two column for Zmd and Zphi and writes their
# averages and standard deviations.
def two_column_calculate_statistics(flread, flwritea, flwriteb, col_num):
    par_old = -999.999
    npar_all = 0
    nparZ = 0
    sumZcs = np.array([0.0, 0.0])
    sumZZcs = np.array([[0.0, 0.0], [0.0, 0.0]])

    flread.seek(0)

    for line in flread:
        val_list = line.split()
        par = float(val_list[0])
        Zmd = float(val_list[col_num])
        Zphi = float(val_list[col_num+1])
        if (Zphi >= -math.pi):
            Zcos = Zmd * math.cos(Zphi)
            Zsin = Zmd * math.sin(Zphi)
            ZZ = np.array([Zcos, Zsin])

            npar_all += 1
    
            if par != par_old:
                if (npar_all > 1):
                    avZcs, CCZcs = avr_CC_cal(sumZcs, sumZZcs, nparZ)
#                   print 'avZcs=', avZcs[0], avZcs[1]
#                   print 'CCZcs=', CCZcs[0,0], CCZcs[0,1], CCZcs[1,0], CCZcs[1,1]
                    avZmd = math.sqrt(avZcs[0] * avZcs[0] + avZcs[1] * avZcs[1])
                    avZphi = math.atan2(avZcs[1], avZcs[0])
                    sigmd, sigphi = compute_error_bars_md_phi(avZcs, CCZcs, \
                    avZmd)
                    flwritea.write('{0:g} {1:g} {2:g}\n'.format( \
                    par_old * factor, avZmd, sigmd))
                    flwriteb.write('{0:g} {1:g} {2:g}\n'.format( \
                    par_old * factor, avZphi, sigphi))
    
                nparZ = 1
                sumZcs = ZZ
                sumZZcs = np.array([[ZZ[0] * ZZ[0], ZZ[0] * ZZ[1]], \
                                    [ZZ[1] * ZZ[0], ZZ[1] * ZZ[1]]])
                par_old = par

            else:
                nparZ += 1
                sumZcs += ZZ
                sumZZcs += np.array([[ZZ[0] * ZZ[0], ZZ[0] * ZZ[1]], \
                                     [ZZ[1] * ZZ[0], ZZ[1] * ZZ[1]]])
    
    if (npar_all > 0):
        avZcs, CCZcs = avr_CC_cal(sumZcs, sumZZcs, nparZ)
        avZmd = math.sqrt(avZcs[0] * avZcs[0] + avZcs[1] * avZcs[1])
        avZphi = math.atan2(avZcs[1], avZcs[0])
        sigmd, sigphi = compute_error_bars_md_phi(avZcs, CCZcs, avZmd)
        flwritea.write('{0:g} {1:g} {2:g}\n'.format(par_old * factor, avZmd, \
        sigmd))
        flwriteb.write('{0:g} {1:g} {2:g}\n'.format(par_old * factor, avZphi, \
        sigphi))
 

#        avZcos, sigZcos = avr_sd_cal(sumZcos, sumZcossq, nparZ)
#        avZsin, sigZsin = avr_sd_cal(sumZsin, sumZsinsq, nparZ)
#        avZmd = math.sqrt(avZcos * avZcos + avZsin * avZsin)
#        avZphi = math.atan2(avZsin, avZcos)
#        flwritea.write('{0:g} {1:g} 0.1\n'.format(par_old, avZmd))
#        flwriteb.write('{0:g} {1:g} 0.1\n'.format(par_old, avZphi))

#main
directory = str(sys.argv[1])
suffix = str(sys.argv[2])

factor = float(sys.argv[3])
xval_min_fr = float(sys.argv[4])
xval_min_cv = float(sys.argv[5])
xval_min_ch = float(sys.argv[6])
print 'suffix=', suffix, 'factor=', factor, ' xval_min=', xval_min_fr, xval_min_cv, xval_min_ch 

fl_avr_name = directory + 'tc.avr.' + suffix
print 'dir=', directory, ' fl_avr_name', fl_avr_name
favr = open(fl_avr_name, 'r')
fsfr = open('avfcc.sfr.' + suffix, 'w')
fscv = open('avfcc.scv.' + suffix, 'w')
fsch = open('avfcc.sch.' + suffix, 'w')

Nval_before = 8
one_column_calculate_statistics(favr, fsfr, Nval_before, xval_min_fr)

Nval_before = 9
one_column_calculate_statistics(favr, fscv, Nval_before, xval_min_cv)

Nval_before = 13
one_column_calculate_statistics(favr, fsch, Nval_before, xval_min_ch)

favr.close()
fsfr.close()
fscv.close()
fsch.close()


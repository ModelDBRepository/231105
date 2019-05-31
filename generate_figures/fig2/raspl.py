#!/usr/bin/python

import sys

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

suffix = str(sys.argv[1]);
print 'suffix=',suffix

fras = open('tc.ras.'+ suffix, 'r')
frplE = open('tc.ras.' + suffix + '.pl.E', 'w')
frplP = open('tc.ras.' + suffix + '.pl.P', 'w')
frplT = open('tc.ras.' + suffix + '.pl.T', 'w')

for line in fras:
    line_list = line.split()
    if line_list[1] == '2':
        frplE.write('%s %s \n' % (line_list[0], line_list[2]))
    elif line_list[1] == '1':
        frplP.write('%s %s \n' % (line_list[0], line_list[2]))
    elif line_list[1] == '0':
        frplT.write('%s %s \n' % (line_list[0], line_list[2]))

fras.close()
frplE.close()
frplP.close()
frplT.close()

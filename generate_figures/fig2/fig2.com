#!/bin/sh
# Plotting Fig. 2.

#../../simulation_program/tc.ex b2b
	      
rast.com
vs.py b2b
spl.py 	      
../../simulation_program/crc.ex b2b
corr_spl.com

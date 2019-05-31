#!/bin/sh
# Plotting Fig. 4.
#

#../../simulation_program/tc.ex p001
#../../simulation_program/tc.ex $051
#../../simulation_program/tc.ex $101

#../../simulation_program/tc.ex p002
#../../simulation_program/tc.ex $052
#../../simulation_program/tc.ex $102

#../../simulation_program/tc.ex p003
#../../simulation_program/tc.ex $053
#../../simulation_program/tc.ex $103

python chikn4.py

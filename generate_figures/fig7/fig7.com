#!/bin/sh
# Plotting Fig. 7.

fla=a2

#../../simulation_program/tc.ex $fla
#../../simulation_program/tc.ex ${fla}1

python vs.py $fla
python spl_ei.py $fla
python vs_e.py $fla
python spl_ei_e.py $fla

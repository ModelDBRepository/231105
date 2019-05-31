#!/bin/sh
# Plotting Fig. 7.

fla=a3

#../../simulation_program/tc.ex ${fla}
#../../simulation_program/tc.ex ${fla}1
#../../simulation_program/tc.ex ${fla}2
#../../simulation_program/tc.ex ${fla}3

cat tc.avr.${fla}1 >> cat tc.avr.${fla}
cat tc.avr.${fla}2 >> cat tc.avr.${fla}
cat tc.avr.${fla}3 >> cat tc.avr.${fla}

fcvk.com

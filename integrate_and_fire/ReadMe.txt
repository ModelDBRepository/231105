this is a temporary readme version.

the Neuron1 folder contains all the c++ files necessary for the program along with the makefile to compile them and an example net.ex resulting executable. 

net.ex is the executable; it needs to be run in the folder Neuron1 with the appropriate run_name (argv[1]) for each of the three runs: iafiKSweep2gNorm, iafiKSweep2tau0gNorm, iafiKSweep2tau1gNorm.

the matching parameter files are inside the folder NeuronIO.

the new results of the program should appear in the folder Neuron1.


the text files in the enclosing folder are the old outputs and the results averaged over the 20 realizations for each K value (in the files anding with _mean).

averaging was performed with th matlab script make_mean.m  .

the resultig files have 4 columns, which are: K, nu_I, chi, CV. (we only used K with chi, this is another matlab script, not included here).
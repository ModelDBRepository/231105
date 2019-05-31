#pragma once
#include "PoissonGen.h"

class NeuronMatrix
{
	//V 0
	//H 1
	//N 2
	//Z 3

public:
	double *** VarVals; //post_populations[2]XnumNeurons[Nt,Ne,Ni]XvariableValues[4]
	double **** spike_buffer;  //post_populations[2]Xpre_pop[3]XnumNeurons[Nt,Ne,Ni]XbufferSpotsForSpikes[32], to be updated in a cyclic manner!
	int buffer_ind;
	double** ready_to_fire;
	double*** traces_sum; //post_popXpre_popXpost_neur
						
	double * thalamic_times_buffer;
	double * save_last_spiketimes;
	double *current_th_spikes;
	Random * rand_gen;
	SystemConstants* sc;
	PoissonGen * pGen;
	//int dbcount = 0;
	//ofstream debug;

	NeuronMatrix();
	NeuronMatrix( PoissonGen *poissGen);
	void reset_neuronMatrix();
	void Update_thalamic_spikes(double* thalamic_spikes_buffer, double*current_th_spikes, double timestep, double time);
	void p_advance_spiketime(double *stepnum_ind, double time_in_run);
	~NeuronMatrix();
};


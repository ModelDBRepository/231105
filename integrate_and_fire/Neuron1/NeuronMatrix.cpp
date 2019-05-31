#include "NeuronMatrix.h"


NeuronMatrix::NeuronMatrix() {
	throw "don't use this constructor";
}

NeuronMatrix::NeuronMatrix(PoissonGen *poissGen)
{

	//debug.open("..\\NeuronIO\\debug.txt");

	rand_gen = poissGen->randGen;// randMain;
	sc = poissGen->sc;
	pGen = poissGen;
	VarVals = new double**[pop_num];
	spike_buffer = new double***[pop_num];
	current_th_spikes = new double[sc->N[0]];
	thalamic_times_buffer = new double[sc->N[0]];
	save_last_spiketimes = new double[sc->N[0]];


	ready_to_fire = new double*[pop_num];
	traces_sum = new double**[pop_num];
	//cout << "nMat1" << endl;
	
	for (int i_post_pop = 0; i_post_pop < pop_num; i_post_pop++) {
		ready_to_fire[i_post_pop] = new double[sc->N[i_post_pop + 1]];
		traces_sum[i_post_pop] = new double*[pop_num+1];

	
		for (int i_pre_p = 0; i_pre_p < pop_num + 1; i_pre_p++) {
			traces_sum[i_post_pop][i_pre_p] = new double[sc->N[i_post_pop + 1]];
			for (int j_post_n = 0; j_post_n < sc->N[i_post_pop + 1]; j_post_n++) {
			}
		}
	}
	//cout << "nMat11" << endl;//******
	//cout << "N0,N1,N2: " << sc->N[0] << ", " << sc->N[1] << ", " << sc->N[2] << endl;
	for (int i_post_pop = 0; i_post_pop < pop_num; i_post_pop++) {
		//cout << "nMat111" << endl;
		VarVals[i_post_pop] = new double*[sc->N[i_post_pop+1]];
		spike_buffer[i_post_pop] = new double**[pop_num+1];
		//cout << "nMat112" << endl;
		for (int i_pre_pop = 0; i_pre_pop < pop_num + 1; i_pre_pop++) {
			spike_buffer[i_post_pop][i_pre_pop] = new double*[sc->N[i_post_pop+1]];
			//cout << "nMat113" << endl;
			
			for (int i_neur = 0; i_neur < sc->N[i_post_pop+1]; i_neur++) {
				/*if (i_pre_pop == 1) {
					cout << "indices: " << i_post_pop << ", " << i_pre_pop << ", N[i_post_pop+1]: " << sc->N[i_post_pop + 1] << endl;
				}*/
				spike_buffer[i_post_pop][i_pre_pop][i_neur] = new double[sc->buffer_size];
				//cout << "nMat113" << endl;
			}
		}
			
	}
	//cout << "nMat2" << endl;//*******
	reset_neuronMatrix();

}

void NeuronMatrix::reset_neuronMatrix(){
	buffer_ind = 0;
	for (int i_t = 0; i_t < sc->N[0]; i_t++) {
		thalamic_times_buffer[i_t] = 0;
		p_advance_spiketime(&thalamic_times_buffer[i_t],0);
		save_last_spiketimes[i_t] = thalamic_times_buffer[i_t];
	}

	for (int n = 0; n < sc->N[0]; n++) {
		current_th_spikes[n] = 0;
	}

	for (int i_post_pop = 0; i_post_pop < pop_num; i_post_pop++) {
		for (int j_neur = 0; j_neur < sc->N[i_post_pop + 1]; j_neur++) {
			ready_to_fire[i_post_pop][j_neur] = 1;
		}
		for (int i_pre_p = 0; i_pre_p < pop_num + 1; i_pre_p++) {
			for (int j_post_n = 0; j_post_n < sc->N[i_post_pop + 1]; j_post_n++) {
				traces_sum[i_post_pop][i_pre_p][j_post_n] = 0;
			}
		}
	}
	//cout << "nMat3" << endl;

	for (int i_post_pop = 0; i_post_pop < pop_num; i_post_pop++) {
		VecDoub infs;
		for (int j_neur = 0; j_neur < sc->N[i_post_pop + 1]; j_neur++) {
			VarVals[i_post_pop][j_neur] = new double[sc->eq_num];
			// ------------initializations:-----------

			//VarVals[i_post_pop][j_neur][0] = -70 + 5 * (double)(j_neur) / (sc->N[i_post_pop + 1] - 1); //2.11.16 equalize
			VarVals[i_post_pop][j_neur][0] = VL;//for iaf check 27.6
			if (sc->eq_num == 4) {
				/*V*/ infs = sc->calcInf(VarVals[i_post_pop][j_neur][0]);
				/*H*/ VarVals[i_post_pop][j_neur][1] = infs[0];//0.985858945;//0.8;
				/*N*/ VarVals[i_post_pop][j_neur][2] = infs[1];// 0.0552263204;//0.2;
				/*z*/ VarVals[i_post_pop][j_neur][3] = infs[2];// 0.0;//0.1; 
			}
		}

		for (int i_pre_pop = 0; i_pre_pop < pop_num + 1; i_pre_pop++) {
			for (int i_neur = 0; i_neur < sc->N[i_post_pop + 1]; i_neur++) {
				for (int k_buf = 0; k_buf <sc->buffer_size; k_buf++) spike_buffer[i_post_pop][i_pre_pop][i_neur][k_buf] = 0;
			}
		}

	}
	//cout << "nMat4" << endl;

}



void NeuronMatrix::Update_thalamic_spikes(double* thalamic_times_buffer, double*current_th_spikes, double timestep, double time) {
	//dbcount++;
	//if ((dbcount % 10) == 0) cout << current_th_spikes[70] << endl;
	for (int i_n = 0; i_n < sc->N[0]; i_n++) {
		current_th_spikes[i_n] = 0;
		int ind = (int)(thalamic_times_buffer[i_n] / timestep);
		if (ind > 0) {
			thalamic_times_buffer[i_n] -= timestep;
		}
		else {
			save_last_spiketimes[i_n] = thalamic_times_buffer[i_n];
			p_advance_spiketime(&thalamic_times_buffer[i_n],time);
			current_th_spikes[i_n] ++;		
		}
	}
}


void NeuronMatrix::p_advance_spiketime(double *spike_time, double time_in_run) {
	double new_spike_time = *spike_time;
	double x_rand;
	double x2_rand;
	double r_max = pGen->Lambda_max_value();
	double r_t;
	if (r_max < 0.000001) return; //check that we actually have spikes and prevent division by zero
	for (int flag = 1,check=1; flag != 0;check++) {
		x_rand = rand_gen->RandomUniform0_to_1();
		//debug << x_rand << endl;
		new_spike_time = new_spike_time - log(x_rand) / r_max;
		x2_rand = rand_gen->RandomUniform0_to_1();
		r_t = pGen->Lambda_rate_function(time_in_run + new_spike_time);
		//debug << time_in_run<<" "<<r_t << " " << r_max <<" At: "<<sc->A_T<< " Bt: " << sc->B_T << " Ct: " << sc->C_T << endl;
		if (r_t / r_max>x2_rand) {
			flag = 0;//this spiketime is chosen
		}
		
		if (check > 1000) {
			printf("possible error in line %d file %s: over 1000 failures to allocate spike", __LINE__, __FILE__);
			throw "error";
		}
	}
	*spike_time = new_spike_time;
	//debug << time_in_run << " " << new_spike_time << endl;
}




NeuronMatrix::~NeuronMatrix()
{ 
	//debug.close();

	for (int i_post_pop = 0; i_post_pop < pop_num; i_post_pop++) {
		for (int i_pre_pop = 0; i_pre_pop < pop_num + 1; i_pre_pop++) {
			for (int i_neur = 0; i_neur < sc->N[i_post_pop+1]; i_neur++) {
				delete[] spike_buffer[i_post_pop][i_pre_pop][i_neur] ;				
			}
			delete[] spike_buffer[i_post_pop][i_pre_pop];
		}
		delete[] spike_buffer[i_post_pop] ;

	}
	delete[] spike_buffer;

	for (int i_post_pop = 0; i_post_pop < pop_num; i_post_pop++) {
		for (int j_neur = 0; j_neur < sc->N[i_post_pop+1]; j_neur++) {
			delete[] VarVals[i_post_pop][j_neur] ;

		}
		delete[] VarVals[i_post_pop];
		delete[] ready_to_fire[i_post_pop];
		for (int i_pre_p = 0; i_pre_p < pop_num + 1; i_pre_p++) {
			delete[] traces_sum[i_post_pop][i_pre_p];
		}
		delete[] traces_sum[i_post_pop];
	
	}

	
	
	delete[] VarVals;

	delete[] ready_to_fire;
	delete[] traces_sum;
	delete[] thalamic_times_buffer;
	delete[] save_last_spiketimes;
	delete[] current_th_spikes;
}

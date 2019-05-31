  #include "SystemConstants.h"



SystemConstants::SystemConstants()
{
	
	gL = new double[num_cell_types];
	gKz = new double[num_cell_types];
	tau_decay = new double[pop_num + 1];

  
	N = new int[pop_num + 1];

	mat_prob = new double*[pop_num]; 
	tau_delay = new double*[pop_num];
	d_tau = new double*[pop_num];
	K_ab = new double*[pop_num];
	K_norm = new double*[pop_num];
	G_ab = new double*[pop_num];
	spikeFromPop_val = new double*[pop_num];
	sig_ki = new double*[pop_num];

	for (int ex_neur_num = 0; ex_neur_num < 10; ex_neur_num++) {
		ex_neurons[ex_neur_num] = new double[2];
	}
	n_ex_neurons = 0;

	for (int i_pop = 0; i_pop < pop_num; i_pop++) {
		mat_prob[i_pop] = new double[pop_num + 1];
		tau_delay[i_pop] = new double[pop_num + 1];
		d_tau[i_pop] = new double[pop_num + 1];
		K_ab[i_pop] = new double[pop_num + 1];
		K_norm[i_pop] = new double[pop_num + 1];
		G_ab[i_pop] = new double[pop_num + 1];
		spikeFromPop_val[i_pop] = new double[pop_num + 1];
		sig_ki[i_pop] = new double[pop_num + 1];
		for (int pre_pop = 0; pre_pop < pop_num + 1; pre_pop++) {
			sig_ki[i_pop][pre_pop] = 0; //initialize sig_ki to default value (no heterogeneity)
		}
	}
	
}



void SystemConstants::setBuffer_size_andDelay() {
	max_delay = tau_delay[0][0];
	max_d_tau = d_tau[0][0];
	for (int i = 0; i < pop_num; i++) {
		for (int j = 0; j < pop_num + 1; j++) {
			max_delay = max_delay > tau_delay[i][j] ? max_delay : tau_delay[i][j];
			max_d_tau = max_d_tau > d_tau[i][j] ? max_d_tau : d_tau[i][j];
		}
	}
	buffer_size = (int)(round((max_delay + max_d_tau) / timestep) + 2);
}

void SystemConstants::setSpikevals_and_matprob() {
	for (int i_post_pop = 0; i_post_pop < pop_num; i_post_pop++) {
		for (int i_pre_p = 0; i_pre_p < pop_num + 1; i_pre_p++) {
			spikeFromPop_val[i_post_pop][i_pre_p] = (G_ab[i_post_pop][i_pre_p] / sqrt(K_ab[i_post_pop][i_pre_p]))
				*(1 / tau_decay[i_pre_p]);
			//cout << "spikeVal: " << spikeFromPop_val[i_post_pop][i_pre_p] << endl << "G_ab[i_post_pop][i_pre_p]: " << G_ab[i_post_pop][i_pre_p] <<
				//"K_ab[i_post_pop][i_pre_p]: " << K_ab[i_post_pop][i_pre_p] << ", 1 / tau_decay[i_pre_p]: " << 1 / tau_decay[i_pre_p] << endl;

			//spikeFromPop_val[i_post_pop][i_pre_p] = (G_ab[i_post_pop][i_pre_p] / sqrt(K_ab[i_post_pop][i_pre_p]*
				//K_ab[i_post_pop][i_pre_p]/K_norm[i_post_pop][i_pre_p])) *(1 / tau_decay[i_pre_p]);
			

      //cout<< "spikeFromPop_val["<<i_post_pop<<"]["<<i_pre_p<<"]: " << spikeFromPop_val[i_post_pop][i_pre_p]<<endl;


			mat_prob[i_post_pop][i_pre_p] = K_ab[i_post_pop][i_pre_p] / N[i_pre_p];
		}
	}
}

VecDoub SystemConstants::calcInf(Doub V) {
	Doub alphah, betah, alphan, betan;
	VecDoub infs = VecDoub(3);
	alphah = 0.7*exp(-(V + 44) / 20);
	betah = 10.0 / (1 + exp(-(V + 14) / 10));
	infs[0] = alphah / (alphah + betah);//Hinf

	

	alphan = 0.1*(V + 34) / (1 - exp(-(V + 34) / 10.0));
	betan = 1.25*exp(-(V + 44) / 80.0);
	infs[1] = alphan / (alphan + betan);//Ninf
	

	infs[2] = 1.0 / (1 + exp(-0.7*(V + 30)));//Zinf
	return infs;

	//----possible extension of the function for use in ODE_file:
	//infs[1] = 1 / (alphah + betah);//tauH
	//infs[3] = 1.0 / (alphan + betan);//tauN
	//-----------------------------------
}

void SystemConstants::setAt(double A_t) {
	A_T = A_t;
}

SystemConstants::~SystemConstants()
{
	delete[] gL ;
	delete[] gKz ;
	delete[] tau_decay;
	delete[] N;
	for (int i_pop = 0; i_pop < pop_num; i_pop++) {
		delete[] mat_prob[i_pop];
		delete[] tau_delay[i_pop] ;
		delete[] d_tau[i_pop];
		delete[] K_ab[i_pop] ;
		delete[] K_norm[i_pop];
		delete[] G_ab[i_pop];
		delete[] spikeFromPop_val[i_pop];
		delete[] sig_ki[i_pop];
	}
	for (int ex_neur_num = 0; ex_neur_num < 10; ex_neur_num++) {
		delete[] ex_neurons[ex_neur_num];
	}
	delete[] mat_prob;
	delete[] spikeFromPop_val;
	delete[] G_ab;
	delete[] K_ab;
	delete[] K_norm;
	delete[] tau_delay;
	delete[] d_tau;
	delete[] sig_ki;

}

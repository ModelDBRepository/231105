#include "PoissonGen.h"

PoissonGen::PoissonGen()//warning: this leaves randomGen uninitiated
{	
	throw "don't use this constructor";
}
PoissonGen::PoissonGen(Random* rand_Gen, SystemConstants* SC)
{
	randGen = rand_Gen;
	sc = SC;
	thalamicTrains = new VecDoub[sc->N[0]];
}

PoissonGen::~PoissonGen() {
	//if (thalamicTrains != nullptr) {
	delete[] thalamicTrains;
	//}
	
}




Doub PoissonGen::Lambda_rate_function(Doub time) {
	//-----------------
	//Lambda(t)=A_T(1+B_T*cos(2pi*t/T+phi)) +   C_T/tau_c *Heavyside(t-t_c)*Heavyside(t_c+tau_c-t)
	//-----------------
	Doub lambda_t;
	Int time_in_cycle;
	time_in_cycle = (int)floor(time);
	time_in_cycle = time_in_cycle % 100; //it is easier to do this with Int- use the modulus for refering to the cycle
	if ((time_in_cycle>=sc->t_c) && (time_in_cycle < sc->t_c + sc->tau_c)) {
		lambda_t = sc->A_T*(1 + sc->B_T*cos(2 * pi*time / sc->T + sc->phi)) + sc->C_T / sc->tau_c;
	}
	else {
		lambda_t = sc->A_T*(1 + sc->B_T*cos(2 * pi*time / sc->T + sc->phi));
	
	}
	return lambda_t;
}
Doub PoissonGen::Lambda_max_value() {
	return  sc->A_T*(1 + sc->B_T) + sc->C_T / sc->tau_c;
}

void PoissonGen::rate_check(VecDoub spikeTrain, Doub time_end, Int bins_per_cycle) {
	Doub Dt = sc->T / bins_per_cycle;
	Int spike_num = 0;
	VecInt spikes_in_bin = VecInt(bins_per_cycle);
	for (Int i_bin = 0; i_bin < bins_per_cycle; i_bin++) {
		spikes_in_bin[i_bin] = 0;
	}
	for (Int cycle_ind = 0; cycle_ind < time_end / sc->T; cycle_ind++) {
		for (Int i_bin = cycle_ind*bins_per_cycle; i_bin < (cycle_ind + 1)* bins_per_cycle; i_bin++) {
			Int bin_in_cycle = i_bin - cycle_ind*bins_per_cycle;
			while (spikeTrain[spike_num] < (i_bin*Dt)) {
				spikes_in_bin[bin_in_cycle]++;
				spike_num++;
			}


		}
	}
	for (Int i_bin = 0; i_bin < bins_per_cycle; i_bin++) {
		std::cout << "number of spikes in bin number " << i_bin << ": " << spikes_in_bin[i_bin] << "\n";
	}

}

VecDoub PoissonGen::PoissonTrain(Doub time_end) {

	int initial = 0;
	VecDoub train_buffer = VecDoub((round(time_end) + 5) * 20/*1.5*/, initial);
	Doub x_rand;
	Int flag_repeat = 0;
	Doub w = 0.5;
	Doub r_max = Lambda_max_value();
	Int spike_num;
	for (spike_num = 1; train_buffer[spike_num - 1] <time_end; spike_num++) {

		x_rand = randGen->RandomUniform0_to_1();
		if (flag_repeat == 0) {//no repeat
			train_buffer[spike_num] = train_buffer[spike_num - 1] - log(x_rand) / r_max;
		}
		else {//repeat
			train_buffer[spike_num] = train_buffer[spike_num] - log(x_rand) / r_max;
			//in repeat, we discard the previous spike and continue from its time to the next time

		}

		x_rand = randGen->RandomUniform0_to_1();

		if (x_rand > (Lambda_rate_function(train_buffer[spike_num])) / r_max) {
			flag_repeat = 1;
			spike_num--;
		}
		else
		{
			flag_repeat = 0;
		}

	}
	VecDoub shortened_train = VecDoub(spike_num - 1);
	for (int ind_spike = 0; ind_spike < spike_num - 1; ind_spike++) {
		shortened_train[ind_spike] = train_buffer[ind_spike];
	}
	return shortened_train;
}

void PoissonGen::CreateTrains(double time_end) {
	
	//if (thalamicTrains != nullptr) { //problematic *** unix
	//	delete[] thalamicTrains;
	//}

	int max_train_length = (round(time_end) + 5) * 100/*1.5*/;
	//thalamicTrains = new VecDoub[sc->N[0]]; //problematic *** unix
	for (int t_ind = 0; t_ind < sc->N[0]; t_ind++)
		thalamicTrains[t_ind] = VecDoub(max_train_length, -1);

	for (int train_ind = 0; train_ind < sc->N[0]; train_ind++) {
		VecDoub p_train = PoissonTrain(time_end);
		for (int spike_ind = 0; spike_ind < p_train.size(); spike_ind++) {

			thalamicTrains[train_ind][spike_ind] = p_train[spike_ind];
		}
	}

}



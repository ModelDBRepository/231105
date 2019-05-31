#include "IntegrationMethod.h"

IntegrationMethod::IntegrationMethod()
{
	throw "don't use this constructor";

}
IntegrationMethod::IntegrationMethod(NeuronMatrix* nMat, Connections *con_mat, ODE_file * ODE_in)
{
	p_n_Mat = nMat;
	p_rand_Main = nMat->rand_gen;
	sc = nMat->sc;
	p_conmat = con_mat;
	ODE = ODE_in;
	k1 = new double[sc->eq_num];
	k2 = new double[sc->eq_num];
	k3 = new double[sc->eq_num];
	k4 = new double[sc->eq_num];
	values_1loop = new double[sc->eq_num];
	rk_midvalues = new double[sc->eq_num];
	derivatives = new double[sc->eq_num];
	mean_inputs = new double**[pop_num];
	for (int i_pop = 0; i_pop < pop_num; i_pop++) {
		mean_inputs[i_pop] = new double *[sc->N[i_pop + 1]];
		for (int j_neur = 0; j_neur < sc->N[i_pop + 1]; j_neur++) {
			mean_inputs[i_pop][j_neur] = new double [pop_num+1];//including thalamic input and sum (last one is sum)
		}
	}


	flag_repeatHist=0;

	logHist = new double*[pop_num+1];
	for (int i_pop = 0; i_pop < pop_num+1; i_pop++) {
		logHist[i_pop] = new double[26];
	}

	for (int i_pop = 0; i_pop < pop_num + 1; i_pop++) {
		fired_spikes[i_pop] = new int[sc->N[i_pop]];
	}
	//------for  chi^2, CV:------:
	for (int i_pop = 0; i_pop < pop_num; i_pop++) {
		v_ij_sums[i_pop]= new double[sc->N[i_pop+1]];
		squared_v_ij_sums[i_pop]=new double[sc->N[i_pop+1]];
		squared_sigmas_ij[i_pop]= new double[sc->N[i_pop + 1]];//chi^2-----
	
		time_last_spike[i_pop]= new double[sc->N[i_pop + 1]];
		interspike_interval[i_pop]= new double[sc->N[i_pop + 1]];
		mean_interval[i_pop]= new double[sc->N[i_pop + 1]];
		mean_squared_interval[i_pop] = new double[sc->N[i_pop + 1]];
		CV_ij[i_pop] = new double[sc->N[i_pop + 1]];//CV-------
		logHist[i_pop] = new double[26];

	}
	//---------------------------
Igap_last=0;

}

IntegrationMethod::~IntegrationMethod()
{
	
	
	for (int i_pop = 0; i_pop < pop_num + 1; i_pop++) {
		delete[] fired_spikes[i_pop] ;
	}

	for (int i_pop = 0; i_pop < pop_num; i_pop++) {
		delete[] v_ij_sums[i_pop];
		delete[] squared_v_ij_sums[i_pop];
		delete[] squared_sigmas_ij[i_pop];//chi^2

		delete[] time_last_spike[i_pop];
		delete[] interspike_interval[i_pop];
		delete[] mean_interval[i_pop];
		delete[] mean_squared_interval[i_pop];
		delete[] CV_ij[i_pop];//CV-------
	}

	for (int i_pop = 0; i_pop < pop_num; i_pop++) {		
		for (int j_neur = 0; j_neur < sc->N[i_pop + 1]; j_neur++) {
			delete[] mean_inputs[i_pop][j_neur];//including thalamic input and sum (last one is sum)
		}
		delete[] mean_inputs[i_pop];
	}
	delete[] mean_inputs;

	for (int i_pop = 0; i_pop < pop_num + 1; i_pop++) {
		delete[] logHist[i_pop];
	}
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	delete[] values_1loop;
	delete[] rk_midvalues;
	delete[] derivatives;

	delete[] logHist;
}


void IntegrationMethod::Time_Evolution(string dirName) {
	double deriv4[4];
	double timeStart = 0;
	double timefinish = sc->runtime; 
	double*** varvals = p_n_Mat->VarVals;
	gap_coupling = sc->Gel / sqrt(sc->Kel);
	//cout << "IM1" << endl;
	if ((sc->method != 'e') && (sc->method != 'r')) throw "please choose 'r' or 'e' as method";
	if ((sc->thalamic_input != 'm') && (sc->thalamic_input != 's')) {
		printf("please choose 'm' (mean) or 's' (stochastic)as method, line %d file %s", __LINE__, __FILE__);
		throw "error";
	}
	if ((sc->p_spikes_per_touch_fr == 1) && (int(sc->runtime - sc->transient) % 100 != 0)) {
		cout << "when p_spikes_per_touch_fr == 1 spike count period must be a whole multiple of the whisking and touch cycle\n";
		throw "";
	}
	
	//ofstream debug;
	//debug.open("..\\NeuronIO\\debugChi.txt");

	//const char* dir_sep = "\\";
	string sep(dir_sep);/// *** unix ***********************


	/*string timestamp = to_string((long double)(((int)
	) % 86400));
	string dirName = ".."+sep+"NeuronIO"+sep + filename + "_" + string(1, sc->method) + "_" +
		string(1, sc->thalamic_input) + "_" +  to_string((long double)sc->runtime) + "_" +
		 to_string((long double)sc->timestep) + "_" + timestamp;// problematic *** unix 14.12 */
	
	
	
	

	ofstream first_e_volt, last_e_volt, first_i_volt, last_i_volt;
	ofstream spikesth0, spikese0 ;//, raster_e, raster_i, synaptic_EI1600;
	ofstream raster_cort[pop_num], raster_th;
	ofstream Synap10,*trace, mean_trace, inputcurrs, inputsum;

	if (sc->p_trace_firstrun == 1) {
		trace = new ofstream[sc->n_ex_neurons];
		for (int ex_n_num = 0; ex_n_num < sc->n_ex_neurons; ex_n_num++) {
			trace[ex_n_num].open(sc->outfilepath + "trace"+to_string((long double)ex_n_num)+".txt");
		}
	}
	if (sc->p_mean_trace_fr == 1) {
		mean_trace.open(sc->outfilepath + "mean_trace.txt");
		//inputsum.open(sc->outfilepath + "inputsum.txt");
	}
	if (sc->p_inputcurrs_firstrun == 1) {
		inputcurrs.open(sc->outfilepath + "inputcurrs.txt");
		//inputsum.open(sc->outfilepath + "inputsum.txt");
	}
	int exm_neur = 30, exm_pop=1;
	//Synap10 << "neuron: " << exm_neur << endl;
	
	/*if (sc->printflag == 1) {
		first_i_volt.open(dirName + sep+"first_i_volt.txt");
		last_i_volt.open(dirName + sep+"last_i_volt.txt");
		spikesth0.open(dirName + sep+"spikesth0.txt");
		raster_th.open(dirName + sep + "raster_th_full.txt");
		raster_cort[0].open(dirName + sep+"raster_e_full.txt");
		raster_cort[1].open(dirName + sep + "raster_i_full.txt");
	}*/
	if (sc->p_rasters_firstrun == 1) {
		for (int i_pop = 0; i_pop < pop_num; i_pop++) {
			raster_cort[i_pop].open(dirName + sep + "raster_pop="+to_string((long double)i_pop)+"_full.txt");
		}
		raster_th.open(dirName + sep + "raster_th_full.txt");
	}

	// calculate nsteps
	double ns = (timefinish - timeStart) / sc->timestep;
	ns = round(ns);		//rounds to nearest integer (returned as double)
	int nsteps = (int)ns;
	sc->timestep = (timefinish - timeStart) / ns;

	
	//for spike counts:---------------------------------------------------------------
	
	for (int i_pop = 0; i_pop < pop_num + 1; i_pop++) {
		//fired_spikes[j_neur] = 0;
		for (int j_neur = 0; j_neur < sc->N[i_pop]; j_neur++) {
			fired_spikes[i_pop][j_neur] = 0;
			//cout <<  ", i: " << i_pop << ", j_neur: " << j_neur << ", fired_spikes: " << fired_spikes[i_pop][j_neur] << endl;
			
		}
		pop_fired_spikes[i_pop] = 0;
		spikes_per_touch[i_pop] = 0;
	}
	int time_in_cycle;
	//---------------------------------------------------------------------------------
	
	//for chi^2:---------------------------------------------
	for (int i_pop = 0; i_pop < pop_num; i_pop++) {
		pop_mean_sums[i_pop] = 0;
		squared_pop_mean_sums[i_pop] = 0;
		for (int j_neur = 0; j_neur < sc->N[i_pop+1]; j_neur++) {
			v_ij_sums[i_pop][j_neur] = 0;
			squared_v_ij_sums[i_pop][j_neur] = 0;
		}
	}
	mean_v = 0;
	//-------------------------------------------------------
	//--------for CV:----------------------------------------
	for (int i_pop = 0; i_pop < pop_num; i_pop++) {
		pop_CV[i_pop]=0;
		for (int j_neur = 0; j_neur < sc->N[i_pop + 1]; j_neur++) {
			time_last_spike[i_pop][j_neur] = 0;
			interspike_interval[i_pop][j_neur] = 0;
			mean_interval[i_pop][j_neur] = 0;
			mean_squared_interval[i_pop][j_neur] = 0;
			CV_ij[i_pop][j_neur] = 0;
		}
	}
	//-------------------------------------------------------
	
	//--------------initialize mean currents----------------------
	for (int i_pop = 0; i_pop < pop_num; i_pop++) 
		for (int j_neur = 0; j_neur < sc->N[i_pop + 1]; j_neur++) 
			for (int input_curr = 0; input_curr < pop_num + 1; input_curr++)
				mean_inputs[i_pop][j_neur][input_curr] = 0;
	//---------------------------------------------------------------------		
			
		
	//cout << "IM2" << endl;
	//------------***************---------- iteration over run time----------------*************------------------
	for (int j = 0; j < nsteps + 1; j++) {
		t = timeStart + j*sc->timestep;
		//cout << "j: " << j << ", t:" << t << endl;//db 15.6
		/*if (sc->printflag == 1) {
			first_i_volt << setprecision(8) << fixed << t << " " << varvals[0][0][0] << "\n";
			last_i_volt << setprecision(8) << fixed << t << " " << varvals[0][sc->N[1] - 1][0] << "\n";
		}*/


		for (int i_pop = 0; i_pop < pop_num; i_pop++) {
#if pop_num==2 //Inhibition change
			ODE->SetNeurType(i_pop);
			//cout << "IM5" << endl;
#endif 
			//cout << "loop: i_pop: " << i_pop << endl;
			for (int j_neur = 0; j_neur < sc->N[i_pop + 1]; j_neur++) {
				/*if (((j == 1261) || (j == 1262)) && (i_pop == 1) && (j_neur == 116)) {//db 19.6.17
					cout << "j: " << j << ",varvals: " << varvals[i_pop][j_neur][0] << endl;//need to use shorter timestep
				}*/
			//	if (j_neur % 50 == 0) cout << "j_neur: " << j_neur<< endl;
				/*if ((i_pop == 1)&&(j_neur==152)) {
					cout << "IM6" <<", N: "<< sc->N[i_pop + 1]<<", j: "<<j_neur<< endl;
				}*/
				
				IAF_fired_spike = 0;
				//------------------------------update variables for calculation of chi^2---------------------
				if (t > sc->transient) {
					v_ij_sums[i_pop][j_neur] += varvals[i_pop][j_neur][0];
					squared_v_ij_sums[i_pop][j_neur] += varvals[i_pop][j_neur][0] * varvals[i_pop][j_neur][0];
					mean_v += varvals[i_pop][j_neur][0];
					
				}
				//--------------------------------------------------------------------------------------------
				
				/*double a0,b0;// f_n components for RK2 IAF, 18.5.17
				a0 = ODE->O_gL + ODE->g_AMPA + ODE->g_GABAa;
				b0 = ODE->O_gL*VL + ODE->g_AMPA*VAMPA + ODE->g_GABAa*VGABA - sc->Iapp + Igap_last;
				if ((j == 100) && (i_pop == 0) && (j_neur == 0)) {//db
					cout << "ODE->O_gL: " << ODE->O_gL << "***, ODE->g_AMPA: " << ODE->g_AMPA << ", ODE->g_GABAa: " << ODE->g_GABAa << endl;
					cout << "sc->Iapp: " << sc->Iapp << ", Igap_last:" << Igap_last << endl;
					cout << "a0: " << a0 << ", b0: " << b0 << endl;

				}*///this is wrong for 2 reasons: first, it gives data from previous neuron and not from previous timestep and 
					//second this would be a_(n-1), b_(n-1) rather than a_n, b_n


				//dealing with ***spikes***:
				int buffer_ind = p_n_Mat->buffer_ind;
				ODE->g_AMPA = 0;

				if (sc->thalamic_input == 's') {//3.11 
					p_n_Mat->traces_sum[i_pop][0][j_neur] += sc->spikeFromPop_val[i_pop][0] *
						p_n_Mat->spike_buffer[i_pop][0][j_neur][buffer_ind];
					p_n_Mat->spike_buffer[i_pop][0][j_neur][buffer_ind] = 0;
					ODE->g_AMPA += p_n_Mat->traces_sum[i_pop][0][j_neur];//E excitation

					
				}
				else {
					// constant through time and not integrated over the durateion of each spike decay.
			//therefore tau_decay does not influence the result (appears inversely inside spike value).
					ODE->g_AMPA += sc->spikeFromPop_val[i_pop][0]*sc->tau_decay[0] *sc->K_ab[i_pop][0]* sc->A_T;
					/*if ((j == 100) && (i_pop == 0) && (j_neur == 0)) {//db 20.6.19
						cout << "g_AMPA: " << ODE->g_AMPA << ",  spikeFromPop_val[i_pop][0]: " << sc->spikeFromPop_val[i_pop][0]<<endl;
						cout << "ODE->Isyn: " << ODE->Isyn <<", varvals[i_pop][j_neur][0]: "<< varvals[i_pop][j_neur][0]<< endl;
					}*/
				}


				//for (int pre_pop = 1; pre_pop < pop_num + 1; pre_pop++) {
				int pre_pop=1;
#if pop_num>1 // Inhibition change
					p_n_Mat->traces_sum[i_pop][pre_pop][j_neur] += sc->spikeFromPop_val[i_pop][pre_pop] *
						p_n_Mat->spike_buffer[i_pop][pre_pop][j_neur][buffer_ind];
					p_n_Mat->spike_buffer[i_pop][pre_pop][j_neur][buffer_ind] = 0;
					ODE->g_AMPA += p_n_Mat->traces_sum[i_pop][pre_pop][j_neur];//g_AMPA!

			   
				pre_pop = 2;
#endif 
					//double debtry = sc->spikeFromPop_val[j_neur][pre_pop];
					//debtry = p_n_Mat->traces_sum[j_neur][pre_pop][j_neur];
					//cout << " [j_neur][pre_pop][j_neur][buffer_ind]" << j_neur << " " << pre_pop<<" " << j_neur<<" " << buffer_ind << endl;
					//debtry = p_n_Mat->spike_buffer[j_neur][pre_pop][j_neur][buffer_ind];
				p_n_Mat->traces_sum[i_pop][pre_pop][j_neur] += sc->spikeFromPop_val[i_pop][pre_pop] *
					p_n_Mat->spike_buffer[i_pop][pre_pop][j_neur][buffer_ind];
				p_n_Mat->spike_buffer[i_pop][pre_pop][j_neur][buffer_ind] = 0;
				ODE->g_GABAa = p_n_Mat->traces_sum[i_pop][pre_pop][j_neur];//GABAa!
				
			   
				if (t > sc->transient) {
					if (pop_num == 1) {
						mean_inputs[i_pop][j_neur][0] += p_n_Mat->traces_sum[i_pop][0][j_neur] * (varvals[i_pop][j_neur][0] - VAMPA);//thalamic input
						mean_inputs[i_pop][j_neur][1] += p_n_Mat->traces_sum[i_pop][1][j_neur] * (varvals[i_pop][j_neur][0] - VGABA);//inh input
					}
					else if (pop_num == 2) {
						mean_inputs[i_pop][j_neur][0] += p_n_Mat->traces_sum[i_pop][0][j_neur] * (varvals[i_pop][j_neur][0] - VAMPA);//thalamic input
						mean_inputs[i_pop][j_neur][1] += p_n_Mat->traces_sum[i_pop][1][j_neur] * (varvals[i_pop][j_neur][0] - VAMPA);//exc input
						mean_inputs[i_pop][j_neur][2] += p_n_Mat->traces_sum[i_pop][2][j_neur] * (varvals[i_pop][j_neur][0] - VGABA);//inh input
					}
				}//need to make this compatible with inh only********





				//---------------initialize variable values----
				for (int i = 0; i < sc->eq_num; i++) {
					values_1loop[i] = varvals[i_pop][j_neur][i];

				}
				
				IsynT = ODE->g_AMPA*(values_1loop[0] - VAMPA);
				IsynI = ODE->g_GABAa*(values_1loop[0] - VGABA);
				Isyn = IsynT + IsynI;
				/*if ((j == 100) && (i_pop == 0) && (j_neur == 0)) {//db 20.6.17
					cout << "2: g_AMPA: " << ODE->g_AMPA << ",  spikeFromPop_val[i_pop][0]: " << sc->spikeFromPop_val[i_pop][0] << 
						", Isyn: "<<Isyn <<endl;
					cout << "ODE->Isyn: " << ODE->Isyn << ", varvals[i_pop][j_neur][0]: " << varvals[i_pop][j_neur][0] << endl;
				}*/
				//mean_IsynT[j_neur] += IsynT;//bug 29.5
				//mean_IsynI[j_neur] += IsynI;

				/*if (j_neur == ex_neur) {
					//Synap10 << ODE->g_AMPA << " "<< ODE->g_GABAa << " "<< values_1loop[0] <<" " << IsynT << " " << IsynI << " " << Isyn << endl;	
				}*/
				

				//---------------calculate gap current------
				double Igap=0;
				//Igap = 0;
				/*if (j_neur != (pop_num - 1)) {
					Igap = 0;
				}
				else {
					for (int other_neuron = 0; other_neuron < sc->N[pop_num];other_neuron++) {
						if (p_conmat->gap_prematrix[j_neur][other_neuron]) {
							Igap += (sc->Gel / sqrt(sc->Kel))*(values_1loop[0] - varvals[j_neur][other_neuron][0]);
						}
					}
				}*/
				//Igap = 0;
				
				if (sc->flag_el == 1) {
					if (i_pop == (pop_num - 1)) {
						for (int con_ind = 0; con_ind < p_conmat->gap_con_num[j_neur]; con_ind++) {
							Igap += (values_1loop[0] - varvals[i_pop][p_conmat->gap_mat[j_neur][con_ind]][0]);
						}
						Igap = Igap*gap_coupling;
					}
				}
				Igap_last = Igap;

				//-----------------Euler/Runge-Kutta---------------------------------------------

				if (strcmp(sc->modelType.c_str(), "WB")==0)// working with Wang-Buzsaki neurons, 15.5.17
				{
					ODE->derivWB(t, values_1loop, derivatives, Igap);
					for (int i = 0; i < sc->eq_num; i++) k1[i] = sc->timestep*derivatives[i];

					if (sc->method == 'r') {//4th order Runge-Kutta
						for (int i = 0; i < sc->eq_num; i++) rk_midvalues[i] = values_1loop[i] + k1[i] / 2;
						ODE->derivWB(t + sc->timestep / 2, rk_midvalues, derivatives, Igap);
						for (int i = 0; i < sc->eq_num; i++) k2[i] = sc->timestep*derivatives[i];

						for (int i = 0; i < sc->eq_num; i++) rk_midvalues[i] = values_1loop[i] + k2[i] / 2;
						ODE->derivWB(t + sc->timestep / 2 + 0.001, rk_midvalues, derivatives, Igap);
						for (int i = 0; i < sc->eq_num; i++) k3[i] = sc->timestep*derivatives[i];

						for (int i = 0; i < sc->eq_num; i++) rk_midvalues[i] = values_1loop[i] + k3[i];
						ODE->derivWB(t + sc->timestep - 0.001, rk_midvalues, derivatives, Igap);
						for (int i = 0; i < sc->eq_num; i++) k4[i] = sc->timestep*derivatives[i];
					}

					for (int i = 0; i < sc->eq_num; i++) {
						if (sc->method == 'r') {//4th order Runge-Kutta
							varvals[i_pop][j_neur][i] +=
								(k1[i] / 6 + k2[i] / 3 + k3[i] / 3 + k4[i] / 6);
						}
						else {// Euler method
							varvals[i_pop][j_neur][i] += k1[i];
						}

					}
				}
				else if ((sc->modelType.compare("IAFv") == 0) || (sc->modelType.compare("IAFi") == 0)|| (sc->modelType.compare("IAFampav") == 0)) 
					// working with integrate and fire neurons, 15.5.17
				{//*********************************&&&&&&&HERE@@@@@@@@@@
					
					double  dt_spike_IAF;
					if (values_1loop[0] > V_th) {
						cout << "voltage: " << values_1loop[0] << endl;
						cout << "j: " << j <<", i_pop: "<<i_pop<<", j_neur: "<<j_neur<< endl;
						printf("error: suprathreshold voltage outside step. line %d file %s", __LINE__, __FILE__);
						throw "error";
					}

					// ------define auxiliary variables    ---------- use a, b instead of alpha, beta
					double a1, b1, V_n_lb, k1_lb, k2_lb, k1_ub, k2_ub;//lb==lowbranch, ub==upperbranch
					double a0, b0;//I need to speak with David about this ... 20.6.19

					if (sc->modelType.compare("IAFv") == 0) {//synaptic input current is voltage-dependent 22.6.17
						a1 = ODE->O_gL + ODE->g_AMPA + ODE->g_GABAa; 
						b1 = ODE->O_gL*VL + ODE->g_AMPA*VAMPA + ODE->g_GABAa*VGABA - sc->Iapp + Igap;
					}
					else if (sc->modelType.compare("IAFi") == 0) {//synaptic input current is independent of neuron voltage 22.6.17
						a1 = ODE->O_gL;
						b1 = ODE->O_gL*VL + ODE->g_AMPA*(VAMPA-VL) + ODE->g_GABAa*(VGABA-VL) - sc->Iapp + Igap;
						if (t == 3.555) {
							cout << "varvals[i_pop][j_neur][0]: "<< varvals[i_pop][j_neur][0]<<"ODE->O_gL*VL: " << ODE->O_gL*VL << ", ODE->g_AMPA*(VAMPA-VL): " <<
								ODE->g_AMPA*(VAMPA - VL) << ", ODE->g_GABAa*(VGABA-VL): " <<
								ODE->g_GABAa*(VGABA - VL) <<			endl << "a1: " << a1 << ", b1: " << b1 << endl;
						}
					}
					else if (sc->modelType.compare("IAFampav") == 0) {//synaptic input current is independent of neuron voltage 22.6.17
						a1 = ODE->O_gL + ODE->g_AMPA;
						b1 = ODE->O_gL*VL + ODE->g_AMPA*VAMPA + ODE->g_GABAa*(VGABA - VL) - sc->Iapp + Igap;
						if (sc->thalamic_input != 'm') {
							cout << "possible error: IAFampav is meant for use with mean thalamic input" << endl;
							throw "error";
						}
					}
					a0 = a1; b0 = b1;//I need to speak with David about this ... 20.6.19
					// ------end of define auxiliary variables    ----------

					if (sc->method == 'r') {//2th order Runge-Kutta
						/*if ((j == 100) && (i_pop == 0) && (j_neur == 0)) {//db 20.6.17
							cout << "ODE->Isyn: " << ODE->Isyn << ", varvals[i_pop][j_neur][0]: " << varvals[i_pop][j_neur][0] <<" check1"<< endl;
							cout << "values_1loop[0]: " << values_1loop[0] << endl;
						}*/
					//make naive step
					//check if threshold has been reached
						//if threshold is reached, find dt_spike_IAF and interpolate with RK2. [make it so adding ref_period would be simple]
						//in RK2 we use derivIAF twice.
						
						
						
						/*if ((j == 100) && (i_pop == 1) && (j_neur == 0)) {//db 20.6.17
							cout << "Inh: ** ODE->O_gL: " << ODE->O_gL << ", ODE->g_AMPA: " << ODE->g_AMPA << ", ODE->g_GABAa: " << ODE->g_GABAa << endl;
							cout << "a1: " << a1 << ", b1: " << b1 << endl;
							
						}*/
						double V_new;
						k1_ub = -a0*values_1loop[0] + b0;
						k2_ub = -a1*(values_1loop[0] + k1_ub*sc->timestep) + b1;
						V_new = values_1loop[0] + sc->timestep*0.5*(k1_ub + k2_ub);
						//cout << "values_1loop[0]: " << values_1loop[0] << endl;
						
						/*if ((j == 1261) && (i_pop == 1) && (j_neur == 116)) {//db 19.6.17
							cout << "j: " << j << ",V_new_before thrsh check: " << V_new << endl;
						}*/
						if (V_new >= V_th) {
							IAF_fired_spike = 1;
							//cout << "V_new1: " <<setprecision(7)<< V_new << endl;
							dt_spike_IAF = ((V_th - values_1loop[0]) / (V_new - values_1loop[0]))*sc->timestep;
							V_n_lb = (V_reset - dt_spike_IAF*(b0 + b1 - a1*b0*sc->timestep) / 2) / (1 + dt_spike_IAF*(-a0 - a1 + a0*a1*sc->timestep) / 2);
							k1_lb = -a0*V_n_lb + b0;//this is an ODE step
							k2_lb = -a1*(V_n_lb + k1_lb*sc->timestep) + b1;//this is an ODE step
							V_new = V_n_lb + sc->timestep*0.5*(k1_lb + k2_lb);
							//cout << "V_new2: " << setprecision(7) << V_new <<", V_n_lb: "<<V_n_lb<<", k1_lb: "<<k1_lb<<", k2_lb: "<<k2_lb<< endl;
							//cout<< setprecision(7)<<"V_new3: "<< V_n_lb + sc->timestep*0.5*(k1_lb + k2_lb);
							/*if (V_new > V_th) {
								cout << "2voltage: " << V_new << endl;
								printf("error: suprathreshold voltage outside step. line %d file %s", __LINE__, __FILE__);
								throw "error";
							}*/
							/*if ((j == 1261) && (i_pop == 1) && (j_neur == 116)) {//db 19.6.17
								cout << "j: " << j << ",V_new_after thrsh check: " << V_new << endl;
							}*/

						}
						varvals[i_pop][j_neur][0] = V_new;
						/*if ((j == 100) && (i_pop == 0) && (j_neur == 0)) {//db 20.6.17
							cout << "ODE->Isyn: " << ODE->Isyn << ", varvals[i_pop][j_neur][0]: " << varvals[i_pop][j_neur][0] << " check2" << endl;
							cout << "values_1loop[0]: " << values_1loop[0] << endl;
						}*/
					}


					else {// Euler method
						values_1loop[0] += sc->timestep*(-a0*values_1loop[0] + b0);
						if (values_1loop[0] >= V_th) {
							dt_spike_IAF = ((V_th - varvals[i_pop][j_neur][0]) / (values_1loop[0] - varvals[i_pop][j_neur][0]))*sc->timestep;
							values_1loop[0] += V_reset - V_th;
							IAF_fired_spike = 1;
							
						}
						varvals[i_pop][j_neur][0] = values_1loop[0];
					}

					t_spike = t + dt_spike_IAF;
					/*
					if (values_1loop[0] >= V_th) {
					dt_spike_IAF = ((V_th - V_n) / (values_1loop[0] - V_n))*sc->timestep;
					V_n_lb = (V_reset - dt_spike_IAF*(b0 + b1 - a1*b0*sc->timestep) / 2) / (1 + dt_spike_IAF*(-a0 - a1 + a0*a1*sc->timestep) / 2);
					rk_midvalues[0] = V_n_lb;
					ODE->derivIAF(t, rk_midvalues, derivatives, Igap);
					k1_lb = sc->timestep*derivatives[0];
					rk_midvalues[0] = V_n_lb+k1_lb ;
					ODE->derivIAF(t+sc->timestep, rk_midvalues, derivatives, Igap);
					k2_lb = -a1*(V_n_lb + k1_lb*sc->timestep) + b1;
					values_1loop[0] = V_n_lb + sc->timestep*0.5*(k1_lb + k2_lb);
					}
					else {
					k1_ub = -a0*values_1loop[0] + b0;
					k2_ub = -a1*(values_1loop[0] + k1_lb*sc->timestep) + b1;
					values_1loop[0] = values_1loop[0] + sc->timestep*0.5*(k1_ub + k2_ub);
					}
					varvals[i_pop][j_neur][0] = values_1loop[0];*/
					

				}
				else
				{
					printf("error: unrecognized model type. line: %d, file: %s.\n", __LINE__, __FILE__);
					throw "error";
				}
				//--------------end of Euler/Runge-Kutta--------------------------------------


			//**************************non-thalamic spikes**************************************
				if (varvals[i_pop][j_neur][0] < -33) p_n_Mat->ready_to_fire[i_pop][j_neur] = 1;
				int ind_delay;
				int summed_cyclic_ind;

				//-------check spike threshold--------
				if ((((varvals[i_pop][j_neur][0] >= spike_thresholdWB) && (values_1loop[0] < spike_thresholdWB))
					&& (p_n_Mat->ready_to_fire[i_pop][j_neur] == 1))||(IAF_fired_spike==1)) {//there is a spike!!

					p_n_Mat->ready_to_fire[i_pop][j_neur] = 0;
					firespike = 1;

					double dt_spike_WB = sc->timestep/*because varvals is for next step*/ - ((varvals[i_pop][j_neur][0] - spike_thresholdWB)*sc->timestep)
						/ (varvals[i_pop][j_neur][0] - values_1loop[0]);

					if (strcmp(sc->modelType.c_str(), "WB") == 0) t_spike = t + dt_spike_WB;

					time_in_cycle = ((int)floor(t_spike)) % 100;//for spikes_per_touch

					
																
																
					//-----------------------for CV--------------------------------
					if (fired_spikes[i_pop + 1][j_neur] > 0) {
						interspike_interval[i_pop][j_neur] = t_spike - time_last_spike[i_pop][j_neur];
						mean_interval[i_pop][j_neur] += interspike_interval[i_pop][j_neur];
						mean_squared_interval[i_pop][j_neur] += interspike_interval[i_pop][j_neur] *
							interspike_interval[i_pop][j_neur];
					}
					//cout << "j: "<<j<<", i: "<<i_pop<<", j_neur: "<<j_neur<<", fired_spikes: " << fired_spikes[i_pop + 1][j_neur] << endl;
					time_last_spike[i_pop][j_neur] = t_spike;
					//-------------------------------------------------------------
					//-----gather data for mean spikes_per_touch and mean firing rates----
					if (t > sc->transient) {
						fired_spikes[i_pop + 1][j_neur]++;
						pop_fired_spikes[i_pop + 1]++;
						if (sc->p_spikes_per_touch_fr == 1) {
							if ((25 <= time_in_cycle) && (time_in_cycle < 50)) spikes_per_touch[i_pop + 1]--;
							if ((50 <= time_in_cycle) && (time_in_cycle < 75)) spikes_per_touch[i_pop + 1]++;
						}
					}
					
					//-------------------------------------------------------------------
					

					//----------print rasters and spikes to files if in print mode-------
					if (sc->p_rasters_firstrun == 1) {
						raster_cort[i_pop] << t_spike << " " << j_neur + 1 << " " << "\n";////db this one seems to be ok
						//if (j_neur == 0) spikese0 << dt_spike_IAF << "\n";
						
					}
					//-------------------------------------------------------------------
					
					//NOTE 23.5.17: it is not right that crucial algorithm parts  are freely intermingled with various prints. this should be fixed...

					// ------******-loop for inserting spikes to the appropriate destinations with delays---*****
					for (int i_receive_pop = 0; i_receive_pop < pop_num; i_receive_pop++) {
						int con_num = p_conmat->num_connections(i_receive_pop, i_pop + 1, j_neur);
						int j_receive_n;
						for (int receive_ind = 0; receive_ind < con_num; receive_ind++) {

							j_receive_n = p_conmat->connections[i_receive_pop][i_pop + 1][j_neur][receive_ind];
							ind_delay = (int)(((t_spike - t) + sc->tau_delay[i_receive_pop][i_pop + 1] +
								p_conmat->delay_mat[i_receive_pop][i_pop + 1][j_neur][receive_ind]) / sc->timestep);


							summed_cyclic_ind = (p_n_Mat->buffer_ind + 1 + ind_delay) % sc->buffer_size;
							/*if ((j == 1254) && (i_pop == 1) && (j_neur == 97) && (i_receive_pop == 1) && (receive_ind == 22)) {
								cout << "check524" << endl;
								cout << "t_spike: " << t_spike << ", t: " << t << endl <<
									", sc->tau_delay[i_receive_pop][i_pop + 1]: " << sc->tau_delay[i_receive_pop][i_pop + 1] <<
									", delay_mat: " 
									<< p_conmat->delay_mat[i_receive_pop][i_pop + 1][j_neur][receive_ind]  <<endl;
								cout << "p_n_Mat->buffer_ind: " << p_n_Mat->buffer_ind << ", ind_delay: " << ind_delay << endl;
								cout << "j_receive_n: " << j_receive_n << ", " << "summed_cyclic_ind: " << summed_cyclic_ind << endl;
							}*/ //db 19.6.17
							p_n_Mat->spike_buffer[i_receive_pop][i_pop + 1][j_receive_n][summed_cyclic_ind]++;
							
						}
					}//-------***------end of spike distribution loop----*****-
				}//--- end of over-spike-threshold if---------------------------------------------
				else {
					firespike = 0;
				}
				//*********************end of dealing with non-thalamic spikes**************
				
				//------print GABAa values for a specific neron--------*****
					/*if ((j_neur == 0) && (j_neur == 1599) && (sc->printflag == 1)) {
						synaptic_EI1600 << setprecision(8) << fixed << t << " "<< ODE->g_GABAa * 2 / 0.038103174 << "\n";
					}*/
					//-------endof print condition---------********************

				//decrease traces_sum for different decay times
				for (int pre_pop = 0; pre_pop < pop_num + 1; pre_pop++) {

					p_n_Mat->traces_sum[i_pop][pre_pop][j_neur] = exp(-sc->timestep / sc->tau_decay[pre_pop])*
						p_n_Mat->traces_sum[i_pop][pre_pop][j_neur];
				}


				/*if ((j == 100) && (i_pop == 0) && (j_neur == 0)) {//db 20.6.17
					cout << "3: g_AMPA: " << ODE->g_AMPA << ",  spikeFromPop_val[i_pop][0]: " << sc->spikeFromPop_val[i_pop][0] <<
						", Isyn: " << Isyn << endl;
					cout << "ODE->Isyn: " << ODE->Isyn << ", varvals[i_pop][j_neur][0]: " << varvals[i_pop][j_neur][0] << endl;
				}*/
				for (int i_ex_neur = 0; i_ex_neur < sc->n_ex_neurons; i_ex_neur++) {//loop over example neurons
					if ((sc->p_trace_firstrun == 1)&&(i_pop == sc->ex_neurons[i_ex_neur][0]) && (j_neur == sc->ex_neurons[i_ex_neur][1])) {
						trace[i_ex_neur] << t << " " << Isyn << " " << values_1loop[0] << " " << firespike << endl;//firespike is only relevant for debugging, not for graphs 13.6.17
																										 ////db				
					}
					
					if (sc->p_inputcurrs_firstrun == 1 && (i_pop == sc->ex_neurons[i_ex_neur][0]) && (j_neur == sc->ex_neurons[i_ex_neur][1])) {						
									//ODE->exm_g_th = p_n_Mat->traces_sum[i_pop][0][j_neur];
						//printing thalamic input currents
						inputcurrs << t << " " << p_n_Mat->traces_sum[i_pop][0][j_neur] *(varvals[i_pop][j_neur][0] - VAMPA) << " ";
						//printing exc input currents
						inputcurrs << p_n_Mat->traces_sum[i_pop][1][j_neur]*(varvals[i_pop][j_neur][0] - VAMPA) << " ";
						//printing of inh input currents
						inputcurrs << ODE->g_GABAa*(varvals[i_pop][j_neur][0] - VGABA) << endl;
						//print sum:
						inputsum << t << " "<< ODE->g_AMPA*(varvals[i_pop][j_neur][0] - VAMPA) + ODE->g_GABAa*(varvals[i_pop][j_neur][0] - VGABA) << endl;
						
					}////db
					

				}
			}//---------------------end of loop over neurons in population-------------------

			 //------------------------------update variables for calculation of chi^2---------------------
			if (t > sc->transient) {
				mean_v = mean_v / sc->N[i_pop + 1];
				if (sc->p_mean_trace_fr == 1) {
					mean_trace << t << " " << mean_v << endl;
				}
				pop_mean_sums[i_pop] += mean_v;
				squared_pop_mean_sums[i_pop] += mean_v*mean_v;
				double sigma_V_squared_db[pop_num];
				sigma_V_squared_db[i_pop] = squared_pop_mean_sums[i_pop]/(j+1) - pop_mean_sums[i_pop] * pop_mean_sums[i_pop]/(j+1)/(j+1);
				/*if ((t > 5999.95)&&(i_pop==0)) {//db 19.6.`7
					cout << "j: " << j << ", t: " << t << "sigma_V_squared_db[" << i_pop << "] " << sigma_V_squared_db[i_pop] << endl;
					cout << "squared_pop_mean_sums[i_pop]/(j+1): " << squared_pop_mean_sums[i_pop] / (j + 1) << endl;
					cout << "pop_mean_sums[i_pop] * pop_mean_sums[i_pop]/(j+1)/(j+1): "
						<< pop_mean_sums[i_pop] * pop_mean_sums[i_pop] / (j + 1) / (j + 1) << endl << endl;
				}
				if (sigma_V_squared_db[i_pop] < 0) {//db
					cout << "j: " << j << ", t: " << t << "sigma_V_squared_db[" << i_pop << "]" << sigma_V_squared_db[i_pop] << endl;
					throw "error";
				}*/

				mean_v = 0;
			}

			//--------------------------------------------------------------------------------------------

		}//----------------------------------end of loop over populations----------------------------------------


		//************************************************thalamic spikes**************************************************

		if (sc->thalamic_input == 's') {

			for (int i_t = 0; i_t < sc->N[0]; i_t++) {
				if (p_n_Mat->current_th_spikes[i_t] != 0) {//thalamic spike if

					double t_spike = t - sc->timestep + p_n_Mat->save_last_spiketimes[i_t];
					//-------------prints----------------------
					//if ((i_t == 0)&&(sc->printflag==1)) spikesth0 << setprecision(5) << fixed << dt_spike_IAF << endl;			
					if (sc->p_rasters_firstrun == 1) {
						raster_th << setprecision(5) << fixed << t_spike << " " << i_t+1 << endl;
						//cout << "rasterdb" << endl;
					}
					if (t_spike > sc->transient) {
						pop_fired_spikes[0]++;
						fired_spikes[0][i_t]++;
						time_in_cycle = ((int)floor(t_spike)) % 100;
						if (sc->p_spikes_per_touch_fr == 1) {
							if ((25 <= time_in_cycle) && (time_in_cycle < 50)) spikes_per_touch[0]--;
							if ((50 <= time_in_cycle) && (time_in_cycle < 75)) spikes_per_touch[0]++;
						}
					}
					//----end of prints---------------------

					// -----------loop for inserting spikes to the appropriate destinations with delays---
					for (int i_receive_pop = 0; i_receive_pop < pop_num; i_receive_pop++) {
						int num_con = p_conmat->num_connections(i_receive_pop, 0, i_t);
						for (int con_ind = 0; con_ind < num_con; con_ind++) {
							int j_neur = p_conmat->connections[i_receive_pop][0][i_t][con_ind];
							int ind_delay = (((t_spike - t) + sc->tau_delay[i_receive_pop][0] +//use ROUND?
								p_conmat->delay_mat[i_receive_pop][0][i_t][con_ind]) / sc->timestep);
							if (ind_delay < 0) ind_delay = 0;

							int summed_cyclic_ind = (p_n_Mat->buffer_ind + 1 + ind_delay) % sc->buffer_size;
							p_n_Mat->spike_buffer[i_receive_pop][0][j_neur][summed_cyclic_ind]++;

						}
					}//---end of spike distribution loop-----

				}//-----end of current spike if------

			}//----end of loop over thalamic neurons---

		}//end of 's' thalamic spikes if
		 //******************************************************************************************************



		//***spike buffer ind update:***
		p_n_Mat->buffer_ind = (p_n_Mat->buffer_ind + 1) % sc->buffer_size;
		if (sc->thalamic_input == 's') {
			p_n_Mat->Update_thalamic_spikes(p_n_Mat->thalamic_times_buffer, p_n_Mat->current_th_spikes, sc->timestep, t);
			//if (p_n_Mat->current_th_spikes[0] == 1) dbcountspikes++;
			//if ((j % 100) == 0) cout << "t: " << t << ", countspikes: " << dbcountspikes << endl;
		}//*************************
	}// end of timestep loop---------------------------------------------------------------------------------------



	// ******************* print to files: mean spikes per touch, mean firing rate********************* 
	
	double spike_count_period = sc->runtime - sc->transient;
	/*ofstream mean_rate_writer;
	mean_rate_writer.open(".."+sep+"NeuronIO" + sep + "firing_rates.txt", std::ios_base::app);// problematic *** unix 14.12 
	if (sc->thalamic_input == 's') {
		mean_rate_writer << "s: " << setprecision(5) << sc->C_T << " ";
		for (int j_neur = 0; j_neur < pop_num + 1; j_neur++) {
			mean_rate_writer << setprecision(5) << double(pop_fired_spikes[j_neur]) / sc->N[j_neur] / spike_count_period
				<< " ";
		}
		mean_rate_writer << endl;
	}
	else {
		mean_rate_writer <<setprecision(5) << sc->C_T << " ";
		for (int j_neur = 0; j_neur < pop_num; j_neur++) {
			mean_rate_writer << setprecision(5) << double(pop_fired_spikes[j_neur+1]) / sc->N[j_neur] / spike_count_period
				<< " ";
		}
		mean_rate_writer << endl;
	}
	mean_rate_writer.close();*/

	if (sc->p_spikes_per_touch_fr == 1) {
		ofstream spikes_per_touch_wr;
		spikes_per_touch_wr.open(".." + sep + "NeuronIO" + sep + "spikes_per_touch.txt", std::ios_base::app);
		for (int j_neur = 0; j_neur < pop_num + 1; j_neur++) {
			spikes_per_touch_wr << setprecision(5) << sc->C_T << " " <<
				spikes_per_touch[0] * 100 / sc->N[0] / spike_count_period << " ";
		}
		spikes_per_touch_wr << endl;
		spikes_per_touch_wr.close();
	}
	//-----------------------------------------------------------------------------------------------
	//--------------------------------calculate and print chi^2, CV-------------------------------------------------
	
	

	firing_rate[0]= double(pop_fired_spikes[0]) / sc->N[0] / spike_count_period;
	for (int i_pop = 0; i_pop < pop_num; i_pop++) {
		firing_rate[i_pop + 1] = double(pop_fired_spikes[i_pop + 1]) / sc->N[i_pop+1] / spike_count_period;
		mean_sigma_v_ij_squared[i_pop] = 0;
		for (int j_neur = 0; j_neur < sc->N[i_pop + 1]; j_neur++) {
			v_ij_sums[i_pop][j_neur] = v_ij_sums[i_pop][j_neur] * sc->timestep / spike_count_period;
			squared_v_ij_sums[i_pop][j_neur] = squared_v_ij_sums[i_pop][j_neur] * sc->timestep / spike_count_period;
			squared_sigmas_ij[i_pop][j_neur] = squared_v_ij_sums[i_pop][j_neur] - v_ij_sums[i_pop][j_neur] *v_ij_sums[i_pop][j_neur];
			mean_sigma_v_ij_squared[i_pop] += squared_sigmas_ij[i_pop][j_neur];
		}
		mean_sigma_v_ij_squared[i_pop] = mean_sigma_v_ij_squared[i_pop] / sc->N[i_pop + 1];
	}
	
	for (int i_pop = 0; i_pop < pop_num; i_pop++) {
		
		squared_pop_mean_sums[i_pop] = squared_pop_mean_sums[i_pop] * sc->timestep / spike_count_period;
		pop_mean_sums[i_pop]  = pop_mean_sums[i_pop] * sc->timestep / spike_count_period;
		sigma_V_squared[i_pop] = squared_pop_mean_sums[i_pop] -pop_mean_sums[i_pop]*pop_mean_sums[i_pop];
		//cout << "pop_mean_sums: " << pop_mean_sums[j_neur] << ", squared_pop_mean_sums: " 
			//<< squared_pop_mean_sums[j_neur] << endl;
		if ((-0.000001 < mean_sigma_v_ij_squared[i_pop]) && (mean_sigma_v_ij_squared[i_pop] < 0.000001)) {
			chi[i_pop] = -1;
		}
		else {
			chi_squared[i_pop] = sigma_V_squared[i_pop] / mean_sigma_v_ij_squared[i_pop];//chi^2
			if ((chi_squared[i_pop] < 0) && (chi_squared[i_pop] != -1)) {
				cout << "squared_pop_mean_sums[i_pop]: " << squared_pop_mean_sums[i_pop] << endl;
				cout << "pop_mean_sums[i_pop]*pop_mean_sums[i_pop]: " << pop_mean_sums[i_pop] * pop_mean_sums[i_pop] << endl;
				cout << "mean_sigma_v_ij_squared[i_pop]: " << mean_sigma_v_ij_squared[i_pop] << endl;
				cout << "sigma_V_squared[i_pop]: " << sigma_V_squared[i_pop] << endl << endl;
				cout << "*** chi_squared[" << i_pop << "]: " << chi_squared[i_pop] << endl;//db up to here
				printf("error: chi_squared is negative. line: %d, file: %s.\n", __LINE__, __FILE__);
				throw "error";
			}
			chi[i_pop] = sqrt(chi_squared[i_pop]);//chi
		}

		
		

		//cout << "pop_CV[j_neur] initial" << pop_CV[j_neur] << endl;
		int count_participating_neurs = 0;
		for (int j_neur = 0; j_neur < sc->N[i_pop + 1]; j_neur++) {
			mean_interval[i_pop][j_neur] = mean_interval[i_pop][j_neur] / (fired_spikes[i_pop + 1][j_neur]-1 );//fix
		//	if (j_neur == 1082) {
			//	cout << "IM7: fired_spikes[i_pop + 1][j_neur]: " << fired_spikes[i_pop + 1][j_neur] <<
		//			", mean_squared_interval[i_pop][j_neur]: " << mean_squared_interval[i_pop][j_neur] << endl;
		//	}
			mean_squared_interval[i_pop][j_neur] = mean_squared_interval[i_pop][j_neur] / (fired_spikes[i_pop + 1][j_neur]-1 );
			//cout << " im3 fired_spikes[1][1082]: " << fired_spikes[1][1082] << endl;
			sigma_cvij_squared = mean_squared_interval[i_pop][j_neur] -
				mean_interval[i_pop][j_neur] * mean_interval[i_pop][j_neur];
			//cout << "i_pop: " << i_pop << ", j_neur: " << j_neur <<", fired_spikes: "<< fired_spikes[i_pop + 1][j_neur]<<
			//	", mean_squared_interval: " << mean_squared_interval[i_pop][j_neur] << ", mean_interval^2" <<
			//	mean_interval[i_pop][j_neur] * mean_interval[i_pop][j_neur] << endl;
			if (sigma_cvij_squared < -0.0001) {
				printf("possible error: sigma_cvij_squared is negative. line: %d, file: %s.\n", __LINE__,__FILE__);
				throw "error";
			}
			else if (0> sigma_cvij_squared) {
				sigma_cvij_squared = -sigma_cvij_squared;
			}

			CV_ij[i_pop][j_neur] = sqrt(sigma_cvij_squared)/ mean_interval[i_pop][j_neur];
			//cout << "sigma_cvij: " << sigma_cvij_squared << endl;
			if (mean_interval== 0) {
				printf("possible error: mean interval is zero. line: %d, file: %s.\n", __LINE__, __FILE__);
				throw "error";
			}
			if (fired_spikes[i_pop+1][j_neur]>2) {//***CV limit check
				//cout << "CV_ij["<<j_neur<<"]["<<j_neur<<"]: " << CV_ij[j_neur][j_neur] << endl;
				pop_CV[i_pop] += CV_ij[i_pop][j_neur];
				//cout << "pop_CV[j_neur] non-normalized in-loop" << pop_CV[j_neur] << endl;
				count_participating_neurs++;
			}
			//cout << "pop_CV[j_neur] non-normalized in-loop" << pop_CV[j_neur] << endl;
		}
		//pop_CV[j_neur] = (sqrt(mean_squared_interval[j_neur] - mean_interval[j_neur] * mean_interval[j_neur])) / mean_interval[j_neur];
		//cout << "pop_CV[j_neur] non-normalized" << pop_CV[j_neur] << endl;
		pop_CV[i_pop] = pop_CV[i_pop] / count_participating_neurs;
		if (pop_CV[i_pop] != pop_CV[i_pop]) {
			pop_CV[i_pop] = -1;
			//this is to switch NaN with -1 for easier file data read\writing.
		}
		//cout << "num_participating in CV: " << count_participating_neurs << endl;
		//cout << "CV: " << pop_CV[j_neur] << endl;
		//CV----
	}

	//-----create rate histogram--------
	
	for (int i_pop = 0; i_pop < pop_num + 1; i_pop++) {
		if (flag_repeatHist == 0) {
			for (int indHist = 0; indHist < 26; indHist++) logHist[i_pop][indHist] = 0;
		}
		for (int j_neur = 0; j_neur < sc->N[i_pop]; j_neur++) {
			if (fired_spikes[i_pop][j_neur] == 0) {
				logHist[i_pop][0]++;
			}
			for (int indHist = 1; indHist < 26; indHist++) {
				//cout << fired_spikes[j_neur][j_neur] / spike_count_period << ", " << pow(10, (-11 + indHist) / 5) << " ";
				if ((fired_spikes[i_pop][j_neur] / spike_count_period*1000>=pow(10, (-11 + indHist)*0.2)) &&
					(fired_spikes[i_pop][j_neur] / spike_count_period *1000< pow(10, (-10 + indHist) *0.2))) {
					logHist[i_pop][indHist]++;
				}
			} //SOME NEURONS ARE NOT "REGISTERED". exactly 10 ,1 or 100? nned to fix->switched > to >=
			//cout << endl;
			
		}
	}

	//Print_firingRates();
	//PrintSync();
	//Quiescent_vs_sig();

	//CV_vs_r();
	//PrintHist();
	//PrintRates(1);
	//----------------------------------
	
	//Synap10.close();
	if (sc->p_trace_firstrun == 1) {
		for (int ex_n_num = 0; ex_n_num < sc->n_ex_neurons; ex_n_num++) {
			trace[ex_n_num].close();
		}
		delete[] trace;
	}
	if (sc->p_mean_trace_fr == 1) {
		mean_trace.close();
	}
	if (sc->p_inputcurrs_firstrun == 1) {
		inputcurrs.close();
		inputsum.close();
	}
	/*for (int neurn = 0; neurn < 150; neurn++) {
		mean_IsynT[neurn] = mean_IsynT[neurn] / (sc->runtime / sc->timestep);
		mean_IsynI[neurn] = mean_IsynI[neurn] / (sc->runtime / sc->timestep);
	}
	Synap10.open(sc->outfilepath + to_string((long double)sc->A_T)+".txt" , std::ios_base::app);
	double pop_mean_IsynT=0, pop_mean_IsynI=0;
	for (int neurn = 0; neurn < 150; neurn++) {
		Synap10 <<mean_IsynT[neurn] << " " << mean_IsynI[neurn] << " " << (mean_IsynT[neurn] + mean_IsynI[neurn]) <<
			" " << fired_spikes[1][neurn] / (sc->runtime - sc->transient)*1000 << endl;
		pop_mean_IsynT += mean_IsynT[neurn];
		pop_mean_IsynI += mean_IsynI[neurn];
	}
	pop_mean_IsynT = pop_mean_IsynT/150;
	pop_mean_IsynI = pop_mean_IsynI/150;
	//Synap10 << "population mean: " << pop_mean_IsynT << " " << pop_mean_IsynI << " " << (pop_mean_IsynT + pop_mean_IsynI) <<endl;
	Synap10.close();*/

	//Print_Y_vs_X(sc->A_T * 1000, firing_rate[1]*1000, "I_rate_vs_T_rate");
	//Print_firingRates(2);
	//PrintSync(0);
	//PrintSync(1);
	//CV_vs_r(0);
	//CV_vs_r(1);
	//inputcurrs_vs_r(0);
	//inputcurrs_vs_r(1);
	//Quiescent_vs_sig(0);
	//Quiescent_vs_sig(1);
	//Print_firingRates(1);
	//Print_firingRates(2);
	if (sc->p_sync_fig1 == 1) {
		for (int i_pop = 0; i_pop < pop_num; i_pop++) {
			PrintSync(i_pop);
		}
	}

	

	//*************************************************************************************************
	
	spikesth0.close();
	if (sc->p_rasters_firstrun == 1) {
		for (int i_pop = 0; i_pop < pop_num; i_pop++) raster_cort[i_pop].close();
		raster_th.close();
	}
	//first_i_volt.close();
	//last_i_volt.close();

	/*//print to file-------rate vs indegrees--------------
	ofstream r_vs_in_wr;
	r_vs_in_wr.open(sc->outfilepath+"_r_vs_in.txt", std::ios_base::app);// problematic *** unix 14.12
	for (int neur = 0; neur<sc->N[1]; neur++) {
	r_vs_in_wr << p_conmat->inDegrees_K[0][0][neur] << " " << fired_spikes[1][neur] / spike_count_period * 1000 <<" "<<p_conmat->inDegree_prob[0][0][neur] <<endl;
	}
	r_vs_in_wr.close();
	//----------------------------*/
}
void IntegrationMethod::Print_Y_vs_X(double x, double y, string filename) {
	ofstream Y_vs_X_wr;
	Y_vs_X_wr.open(sc->outfilepath + filename+".txt", std::ios_base::app);
	Y_vs_X_wr << setprecision(6) << fixed << x << " " << y << endl;
	Y_vs_X_wr.close();
}

void IntegrationMethod::Print_firingRates(int pop) {
	ofstream nu_i_wr, nu_dev_wr;
	double sum = 0; 
	double sumsquares = 0;
	double sqr_dev, dev;
	nu_i_wr.open(sc->outfilepath + "_nu_i2pop="+to_string((long double)pop)+".txt", std::ios_base::app);
	nu_dev_wr.open(sc->outfilepath + "_nu_dev2pop=" + to_string((long double)pop) + ".txt", std::ios_base::app);
	//nu_i_wr.open(sc->outfilepath + "_nu_i.txt", std::ios_base::app);
	//nu_dev_wr.open(sc->outfilepath + "_nu_dev.txt", std::ios_base::app);
	for (int j_neur = 0; j_neur < sc->N[pop]; j_neur++) {
		nu_i_wr << setprecision(6) << j_neur<<" "<<setprecision(6) <<
			fired_spikes[pop][j_neur ]/(sc->runtime - sc->transient) << endl;
		sum += fired_spikes[pop][j_neur ] / (sc->runtime - sc->transient);
		sumsquares += (fired_spikes[pop][j_neur ] / (sc->runtime - sc->transient))*(fired_spikes[pop][j_neur ] / (sc->runtime - sc->transient));
	}
	sqr_dev = (sumsquares / sc->N[pop]) - ((sum / sc->N[pop])*(sum / sc->N[pop]));
	if (sqr_dev < 0) {
		printf("sqrt<0, error line %d file %s", __LINE__, __FILE__);
		throw "error";
	}
	nu_dev_wr << setprecision(6) <<fixed << sc->run_param_val << " " << sqrt(sqr_dev) << endl;
	
	nu_i_wr.close();
	nu_dev_wr.close();
}

void IntegrationMethod::PrintSync(int pop) {
	ofstream synchrony_wr;
	synchrony_wr.open(sc->outfilepath + "_sync_pop=" + to_string((long double)pop) + ".txt", std::ios_base::app);

	if (!synchrony_wr) {
		cout << "error opening synchrony file: "<< sc->outfilepath + "_sync.txt" << endl;
		throw "error";
	}
	synchrony_wr << setprecision(6) << fixed <<sc->run_param_val<<" "<< firing_rate[pop + 1] <<
		" " << chi[pop] << " " << pop_CV[pop] << endl;
	
	synchrony_wr.close();
}

void IntegrationMethod::Quiescent_vs_sig(int pop) {
	ofstream Q_vs_sig_wr;
	Q_vs_sig_wr.open(sc->outfilepath + "_Q_vs_sig_pop="+to_string((long double)pop)+".txt", std::ios_base::app);
	Q_vs_sig_wr << setprecision(6) << fixed << sc->sig_ki[0][0] << " " << logHist[pop+1][0] << endl;
	Q_vs_sig_wr.close(); 
}

void IntegrationMethod::CV_vs_r(int pop) {
	ofstream cv_vs_r_wr;
	cv_vs_r_wr.open(sc->outfilepath + "_cv_vs_r_pop="+to_string((long double)pop)+".txt");
	for (int post_neur = 0; post_neur<sc->N[pop+1]; post_neur++) {
		if (CV_ij[pop][post_neur] == CV_ij[pop][post_neur]) {
			cv_vs_r_wr << setprecision(6) << fixed << fired_spikes[pop+1][post_neur] / (sc->runtime - sc->transient) * 1000 << " " << CV_ij[pop][post_neur] << endl;
		}
	}
	cv_vs_r_wr.close();
}

void IntegrationMethod::inputcurrs_vs_r(int pop) {
	ofstream inpcur_vs_r_wr;
	inpcur_vs_r_wr.open(sc->outfilepath + "_inputcurrs_vs_r_pop=" + to_string((long double)pop) + ".txt");
	for (int j_neur = 0; j_neur < sc->N[pop + 1]; j_neur++) {
		inpcur_vs_r_wr <<j_neur<<" "<< setprecision(6) << fired_spikes[pop + 1][j_neur] / (sc->runtime - sc->transient) * 1000 << " ";
		for (int input_curr = 0; input_curr < pop_num + 1; input_curr++) {
			mean_inputs[pop][j_neur][input_curr] = mean_inputs[pop][j_neur][input_curr] / ((sc->runtime - sc->transient) / sc->timestep);
			inpcur_vs_r_wr << setprecision(6)<< mean_inputs[pop][j_neur][input_curr] << " ";
		}
		inpcur_vs_r_wr << endl;
	}


	
	inpcur_vs_r_wr.close();
}

void IntegrationMethod::PrintHist() {
	ofstream hist_wr;
	hist_wr.open(sc->outfilepath + "_hist.txt");// problematic *** unix 14.12 
	for (int histInd = 0; histInd<26; histInd++) {
		hist_wr << logHist[1][histInd] << endl;
	}
	hist_wr.close();
}

void IntegrationMethod::PrintRates(int pop) {
	ofstream rates;
	rates.open(sc->outfilepath +"_" + to_string((long double)pop) +"rates.txt");
	for (int neur = 0; neur<sc->N[pop]; neur++) {
		rates << setprecision(6) << fixed << fired_spikes[pop][neur] / (sc->runtime- sc->transient)* 1000 << endl;		
	}
	rates.close();
}


//was a pre-function for printing mean population rates to "firing_rates.txt" [search for the name], mean spike per touch ("spikes_per_touch.txt")
//and "elapsed_times.txt" inWB_Neuron.
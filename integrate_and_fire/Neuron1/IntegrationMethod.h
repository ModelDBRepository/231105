#pragma once
#include "ODE_file.h"
#include "Connections.h"
#include "NeuronReader.h"



class IntegrationMethod
{
public:

	NeuronMatrix * p_n_Mat;
	Random* p_rand_Main;
	SystemConstants* sc;
	Connections * p_conmat;
	ODE_file * ODE;


	//int dbcountspikes=0;
	double *k1;
	double *k2;
	double *k3;
	double *k4;
	double *values_1loop;
	double *rk_midvalues;
	double *derivatives;
	double t;
	double gap_coupling;
	double Igap_last;

	double IIFvoltage, IIFrefSteps;
	double derivativeIAF[1];
	double t_spike;

	int IAF_fired_spike;
	int* fired_spikes[pop_num + 1];
	int pop_fired_spikes[pop_num + 1];
	int spikes_per_touch[pop_num + 1];
	//------calculate chi^2------:
	double *v_ij_sums[pop_num];
	double *squared_v_ij_sums[pop_num];
	double *squared_sigmas_ij[pop_num];
	double pop_mean_sums[pop_num];
	double squared_pop_mean_sums[pop_num];
	double mean_v;
	//---------------------------
	//------calculate CV------:
	double *time_last_spike[pop_num];
	double *interspike_interval[pop_num];
	double *mean_interval[pop_num];
	double *mean_squared_interval[pop_num];
	double *CV_ij[pop_num];
	double pop_CV[pop_num];
	double sigma_cvij_squared;
	double **logHist;//first entry is for zero rate. the start from 10^-2 with 0.2 skip in exponent
	int flag_repeatHist;

	double sigma_V_squared[pop_num];
	double mean_sigma_v_ij_squared[pop_num];
	double chi_squared[pop_num];
	double chi[pop_num];
	double firing_rate[pop_num + 1];
	double *** mean_inputs;
	

	int firespike;
	double IsynT;
	double IsynI;
	double Isyn;
	//---------------------------

	
	IntegrationMethod();
	IntegrationMethod(NeuronMatrix *nMat, Connections * con_mat, ODE_file * ODE);
	
	~IntegrationMethod();

	void Print_Y_vs_X(double x, double y, string filename);
	void Time_Evolution(string dirName);
	void Print_firingRates(int pop);
	void PrintSync(int pop);
	void Quiescent_vs_sig(int pop);
	void CV_vs_r(int pop);
	void inputcurrs_vs_r(int pop);
	void PrintHist();
	void PrintRates(int pop);

	/*void Time_Evolution(ODE_file ODE, double timeStart, double*** varvals,
		double timefinish, double timeStep, char method, char spike_stream,
		int printflag, ifstream *paramfile, NeuronReader *nRead, string run_Name, string dirName);*/
	/*rk->Time_Evolution(*ode_file, x_0, nMat->VarVals, x_end, SC->timestep, SC->method, SC->thalamic_input
			  , SC->printflag, &paramFile, Reader, run_name+"f_K", dirName);*/
	
};

#pragma once
#include "nr3.h"

/*#ifdef _WIN32
#define mat_path ((const char*)   "C:\\Users\\tomerarg\\Documents\\MATLAB\\")
//const char* dir_sep = "\\";
#elif __unix__
#define mat_path ((const char*)   "./chi_sqr.txt")
//const char* dir_sep = "/";
#endif*/

#define pop_num 1
#define num_cell_types 2
//#define eq_num 1
#define pi 3.1415926535898

//for ODE:
#define VL  -65
#define gNa  100
#define VNa  55
#define gKdr 40
#define VK  -90
#define VAMPA  0
#define VGABA  -85
#define spike_thresholdWB -20

#define V_th -54
#define V_reset -65 //usually needs to be equal to VL 
#define ref_period 2 //not yet in the model!

//temporary, for runs on K:
#define input_pop 0

enum Print {rasters, trace, inputcurrs, spikes_per_touch};
//#define logbin 1.58489319246





//#define N_t_DEF  1
//#define N_e_DEF  1600
//#define N_i_DEF  150
//#define prob_et_DEF 0.1
//#define prob_it_DEF 0.1
//#define prob_ei_DEF 0.05
//#define prob_ei_DEF 0.05
//#define prob_ee_DEF 0.05
//#define prob_ii_DEF 0.05
	//#define post_pop  2//4E, 4I
	//#define pre_pop  3//4E, 4I, thalamus 

	//-----------------------------------
//#define traceAmp 1
//#define traceTimeConst 1
//#define trainsTime 10000
	//-------------------------------------

/*static const double tau =1;
static const double d_tau= 0.5;
static const double timestep = 0.05;
static const int buffer_size = (int)(round((tau + d_tau) / timestep) + 2);

static const double* mat_prob[] = new double*[2];
mat_prob[0] = new double[3];
mat_prob[1] = new double[3];*/


	class SystemConstants

	{
	public:
		char method;
		char thalamic_input;
		string modelType;
   
   double size_factor;

   int seed;
    //int * K_norm;//for only I in cortex, scaling 1/K
    //int K_norm[2];//for only I in cortex, scaling 1/K
   double run_param_val;
		int* N;
		int popn;
		double** mat_prob;
		double** spikeFromPop_val;
		double * tau_decay;
		double** G_ab;
		double Gel;
		double** K_ab;
		double**K_norm;
		double Kel;
		double** tau_delay;
		double** d_tau;
		double A_T;
		double transient, runtime;
		//int printflag;
		double max_delay, max_d_tau;
		int eq_num;

		// ------runnning parameters----------
		string running_parameter;
		double min_val, max_val, par_step, realiz_num;
		double *ex_neurons[10];//population and neuron number for up to 10 example neurons
		int p_rasters_firstrun, p_trace_firstrun, p_mean_trace_fr, p_inputcurrs_firstrun, p_spikes_per_touch_fr;
		int n_ex_neurons;
		int p_sync_fig1, p_mean_rates, p_mean_CV, p_chi, p_quiescent_neurs, p_dev_rates;
		// end of runnning parameters----------

		int flag_el;

		double **sig_ki;//make into array
	
		int buffer_size;
		double  timestep;

		double *gL;
		double * gKz;


		double T; 
		double phi;
		double tau_c;
		double B_T;//0.25
		double C_T;//0.6
		double t_c;

		double Iapp ;
		double phi2 ;

		string outfilepath;
		//string win_IO_path;
		//string linux_IO_path;
		string IO_path;
		string matlab_path;


		SystemConstants();
		void setBuffer_size_andDelay();
		void setSpikevals_and_matprob();

		VecDoub calcInf(Doub V);
		void setAt(double A_T);
		~SystemConstants();
	};

	//static SystemConstants sc; //this produces multiples - can't change, can't work with!


#pragma once
#include "PoissonGen.h"



class ODE_file //:public Derivnfunction
{
public:

	SystemConstants* sc;

	double alpham;
	double betam;
	double Minf;
	double mqh;
	double INa;

	double V;
	double h;
	double n;
	double z;

	double alphah;
	double betah;
	double Hinf;
	double tauH;
	double alphan;
	double betan;
	double Ninf;
	double tauN;
	double Zinf;
	double tauZ;

	double O_gL, O_gKz;
	//int pop_type;

	/*double g_AMPA_th;
	double g_AMPA_e;
	double g_GABAa;*/
	double g_AMPA;
	double g_GABAa;
	//double exm_g_th;
	double Isyn;
	

	ODE_file();
	ODE_file(SystemConstants* SC);
	~ODE_file();
	void derivWB(double time, double initialConditions[], double derivatives[], double Igap);
	void derivIAF(double time, double * initialVoltage, double derivative[], double Igap);

	void SetNeurType(int n_type);

};
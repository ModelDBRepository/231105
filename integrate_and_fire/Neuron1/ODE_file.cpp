#include "ODE_file.h"

ODE_file::ODE_file() {
	throw "don't use this constructor";
}

ODE_file::ODE_file(SystemConstants* SC) {
	sc = SC;
	//pop_type = 2; //inhibitory
	O_gL = sc->gL[1]; //Inhibitory. possibly change for pop_num==3
	O_gKz = sc->gKz[1];


	
	/*g_AMPA_th = 0;
	g_AMPA_e = 0;
	g_GABAa = 0;*/
	
	//exm_g_th = 0;
	g_AMPA = 0;
	g_GABAa = 0;
	Isyn = 0;
	
	

}
ODE_file::~ODE_file() {
}

void ODE_file::derivWB(double time, double  initialConditions[], double derivatives[], double Igap) {
	V = initialConditions[0];
	h = initialConditions[1];
	n = initialConditions[2];
	z = initialConditions[3];
	
	if ((V + 30 > 0.00000001)|| (V + 30<-0.00000001)) {//this was abs
		alpham = (0.1*(V + 30)) / (1 - exp(-0.1*(V + 30)));
	}else {
		alpham = 10;
	}

	betam = 4 * exp(-(V + 55) / 18);
	Minf = alpham / (alpham + betam);
	mqh = (Minf*Minf*Minf*h);
	INa = gNa*mqh*(V - VNa);


	if (V > -120) {
		alphah = 0.7*exp(-(V + 44) / 20);
		betah = 10.0 / (1 + exp(-(V + 14) / 10));
		Hinf = alphah / (alphah + betah);
		tauH = 1 / (alphah + betah);

		alphan = 0.1*(V + 34) / (1 - exp(-(V + 34) / 10.0));
		betan = 1.25*exp(-(V + 44) / 80.0);
		Ninf = alphan / (alphan + betan);
		tauN = 1.0 / (alphan + betan);

		Zinf = 1.0 / (1 + exp(-0.7*(V + 30)));
		/* -------in case we want to do this with an extended version of the inf function used in initialization: ---- 
		infs = calcInf(V);
		Hinf = infs[0];
		tauH = infs[1];

		
		Ninf = infs[2];
		tauN = infs[3];

		Zinf = infs[4];----------------------------------------------------------------------------------*/
	}else {
		Hinf = 1.0;
		tauH = 1.0;
		Ninf = 0.0;
		tauN = 1.0;
		Zinf = 0.0;
	}
	tauZ = 60;
	
	Isyn = g_AMPA*(V - VAMPA) + g_GABAa*(V - VGABA);
	//-----------------------------------------------------------------------------------------------------------------------
	//note: David does this differently; his program has different variables for excitatory and inhibitory Isyn and 
	// additionally Isyn due to spikes is only caculated once per time step, whereas Isyn due to constant thalamic input
	// is calculated inside the general formula (below) 4 times per timestep in rk. this should account for a small difference
	// in numerical results when spikes occur using rk.
	//-----------------------------------------------------------------------------------------------------------------------

	
	/* V' */ derivatives[0] = (-O_gL*(V - VL) - INa - (gKdr*(n*n*n*n) + O_gKz*z)*(V - VK) + sc->Iapp - Isyn)-Igap;
	/* h' */ derivatives[1] = sc->phi2*(Hinf - h) / tauH;
	/* n' */ derivatives[2] = sc->phi2*(Ninf - n) / tauN;
	/* z' */ derivatives[3] = (Zinf - z) / tauZ;


};

void ODE_file::derivIAF(double time, double * initialVoltage, double derivative[], double Igap) {//unused! 22.6.17
	V = initialVoltage[0];
	//this function should only be executed outside the refractory period. 1.5.17

	Isyn = g_AMPA*(V - VAMPA) + g_GABAa*(V - VGABA);
	//-----------------------------------------------------------------------------------------------------------------------
	//note: David does this differently; his program has different variables for excitatory and inhibitory Isyn and 
	// additionally Isyn due to spikes is only caculated once per time step, whereas Isyn due to constant thalamic input
	// is calculated inside the general formula (below) 4 times per timestep in rk. this should account for a small difference
	// in numerical results when spikes occur using rk.
	//-----------------------------------------------------------------------------------------------------------------------

	/* V' */ derivative[0] = (-O_gL*(V - VL)  + sc->Iapp - Isyn) - Igap;

	//note: this is very simple and barely justifies a function. it is only important here that *O_gL* changes with *SetNeurType*; 
	//this is important for when we have excitatory neurons too. 1.5.17

};



void ODE_file::SetNeurType(int n_type){
	//pop_type = n_type;
	O_gL = sc->gL[n_type];
	O_gKz = sc->gKz[n_type];
}



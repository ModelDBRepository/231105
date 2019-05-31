#pragma once
#include "Random.h"
#include "SystemConstants.h"
//#include "Run_parameters.h"



class PoissonGen
{
public:
	Random* randGen;
	SystemConstants* sc;
	VecDoub* thalamicTrains;


	PoissonGen();
	PoissonGen(Random * rand_Gen, SystemConstants* SC);
	~PoissonGen();

	Doub Lambda_rate_function(Doub time);
	Doub Lambda_max_value();
	void rate_check(VecDoub spikeTrain, Doub time_end, Int num_bins);

	VecDoub PoissonTrain(Doub time_end);
	void CreateTrains(double time_end);
};


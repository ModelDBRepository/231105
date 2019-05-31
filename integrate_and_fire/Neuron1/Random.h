#pragma once
#include "nr3.h"
#include "ran.h"
#include "deviates.h"


class Random
{
public:
	int count;
	
	  //std::default_random_engine generator;
  	//std::uniform_real_distribution<double> distribution;
    Ran rand;
	Normaldev_BM normaldist;
  
	Random();
 
  Random(int seed);
  //void resetRandom(int seed);
	Doub RandomUniform0_to_1(); 
	~Random();
};


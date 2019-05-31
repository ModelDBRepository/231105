#pragma once
#include "SystemConstants.h"
//#include "Run_parameters.h"

using namespace std;


class NeuronReader
{
public:
	SystemConstants *sc;
	//Run_parameters *runpar;

	NeuronReader();
	NeuronReader(SystemConstants* SC/*, Run_parameters *RP*/);
	/*void parse(ifstream & parametersFile, ofstream & parseFile);*/
	void initialize(ifstream& parametersFile, double size_fact);
	void copyParamFile(ifstream& fileToRead,ofstream& fileToWrite);
	~NeuronReader();
};


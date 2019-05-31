#pragma once
#include "IntegrationMethod.h"
#include "NeuronReader.h"
//#include "Run_parameters.h"

void parameterSweep(double min_val, double max_val, double par_step, int realiz_num, IntegrationMethod *IM, string run_name, string dirName);
void set_running_parameter(SystemConstants * sc, string running_parameter, double value);
void f_kSweep(int kii_min, int kii_max, int step, int num_realiz, IntegrationMethod *IM, string run_name, string dirName);
void AtSweep(double at_min, double at_max, double step, int num_realiz, IntegrationMethod *IM, string run_name, string dirName);
void GelSweep(double gel_min, double gel_max, double step, int num_realiz, IntegrationMethod *IM, string run_name, string dirName);
void SigSweepAll(double sig_min, double sig_max, double step, int num_realiz, IntegrationMethod *IM, string run_name, string dirName);
void full_reset_and_run(IntegrationMethod *IM, string dirName);
void createRepeatedHist(IntegrationMethod *IM, string run_name, string dirName, int num_realiz);
void set_outpath(SystemConstants * sc, string dirName, string addition);
//void printRunningParToSync(SystemConstants *sc, double runningParameter);
void printProbsAndInDegrees(SystemConstants * sc, Connections* con, string run_name);

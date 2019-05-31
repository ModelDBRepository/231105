#pragma once
#include "NeuronMatrix.h"

class Connections
{
public:

	Random * randGen;
	SystemConstants * sc;
	
	int ****connections;
	double ****delay_mat;
	double ***inDegree_prob;
	double ***inDegrees_K;
	int **gap_prematrix;
	int **gap_mat;
	int *gap_con_num;
	Normaldev_BM ki_dist;
	/*int*base_mat;
	int**pbase_mat;
	int***ppbase_mat;
	int sum_sizes;
	int total_mat_size;
	int num_tot_preNeurons;*/
	ofstream* debug_mat;
	//double****randVals;
	

	Connections();
	Connections(Random *randGen, SystemConstants* SC);
	
	void CreateConnectionMatrix(Random* randGen);
	//void ReCreateConnectionMatrixSamerand(Random* randGen);
	void Set_ki_Prob(double mean_prob, double relative_deviation, int post_pop, int pre_pop);
 void Set_ki_Prob_FULLCOR(double mean_prob_it, double mean_prob_ii, double sigit, double sigii);
	void partial_matrix_fill( int post_pop_num, int post_num_neurons, int pre_pop_num, 
		int pre_num_neurons, /*double connection_p,*/ Random * randGen);
	//void partial_matrix_REfill_samerand(int post_pop_ind, int post_num_neurons, int pre_pop_ind,
		//int pre_num_neurons/*, double connection_p */, Random *randGen);
	void clean_connection_matrix();
	int num_connections(int post_pop_num, int pre_pop_num, int pre_ind);
	void Print_inDegrees(string path);
	void CreateGapMatrix();

	~Connections();
	
};







#include "Connections.h"
// t 0
// e 1
// i 2

Connections::Connections()
{
	throw "don't use this constructor";
}

Connections::Connections(Random* rand_Gen, SystemConstants* SC)
{
	randGen = rand_Gen;
	sc = SC;
	//CreateConnectionMatrix(sc->N,sc->mat_prob, randGen);
	// ********************TREE ARRAY VERSION******************

	connections = new int ***[pop_num]; //post populations
	delay_mat = new double***[pop_num];//rand delay
	inDegree_prob = new double**[pop_num];
	inDegrees_K = new double**[pop_num];
	//randVals = new double***[pop_num];
	for (int post_pop = 0; post_pop <pop_num; post_pop++) {
		connections[post_pop] = new int **[pop_num + 1]; //pre-populations
		delay_mat[post_pop] = new double **[pop_num + 1];
		inDegree_prob[post_pop] = new double *[pop_num + 1];
		inDegrees_K[post_pop] = new double *[pop_num + 1];
		//randVals[post_pop] = new double **[pop_num + 1];
		for (int pre_pop = 0; pre_pop < pop_num + 1; pre_pop++) {
			connections[post_pop][pre_pop] = new int *[sc->N[pre_pop]];
			delay_mat[post_pop][pre_pop] = new double *[sc->N[pre_pop]];
			inDegree_prob[post_pop][pre_pop] = new double[sc->N[post_pop + 1]];
			inDegrees_K[post_pop][pre_pop] = new double[sc->N[post_pop + 1]];
			//randVals[post_pop][pre_pop] = new double*[sc->N[pre_pop ]];
			//for (int pre_n = 0; pre_n < sc->N[pre_pop]; pre_n++) {
				//randVals[post_pop][pre_pop][pre_n] = new double[sc->N[post_pop + 1]];
			//}
		}
	}
	gap_prematrix = new int *[sc->N[pop_num]];
	gap_mat = new int *[sc->N[pop_num]];
	gap_con_num = new int[sc->N[pop_num]];
	for (int i_n = 0; i_n < sc->N[pop_num]; i_n++) {
		gap_prematrix[i_n] = new int[sc->N[pop_num]];
	}

	// *********************************************************

	CreateConnectionMatrix(randGen);
}

void Connections::CreateConnectionMatrix (Random* randGen) {
	//ofstream connections;
	//connections.open(sc->IO_path);
	sc->setSpikevals_and_matprob();
	//cout << "conmat1" << endl;
	for (int ind_post_pop = 0; ind_post_pop < pop_num ; ind_post_pop++) {
    //Set_ki_Prob_FULLCOR(sc->mat_prob[0][0],sc->mat_prob[0][1],sc->sig_ki[0][0],sc->sig_ki[0][1]);

		//cout << "conmat2, postpop "<<ind_post_pop << endl;
		for (int ind_pre_pop = 0; ind_pre_pop < pop_num + 1; ind_pre_pop++) {
			//cout << "conmat3, postpop " << ind_post_pop <<", prepop "<<ind_pre_pop<< endl;
			//cout << "prob: " << sc->mat_prob[ind_post_pop][ind_pre_pop] << endl;
			Set_ki_Prob(sc->mat_prob[ind_post_pop][ind_pre_pop], sc->sig_ki[ind_post_pop][ind_pre_pop],ind_post_pop,ind_pre_pop);
			//cout << "conmat4, postpop " << ind_post_pop << ", prepop " << ind_pre_pop << endl;
			partial_matrix_fill(ind_post_pop, sc->N[ind_post_pop+1], ind_pre_pop,
				sc->N[ind_pre_pop]/*, sc->mat_prob[ind_post_pop][ind_pre_pop] */,randGen);
			//cout << "createMatrix " << ind_post_pop << ": " << randGen->count << endl;
		}
	}
	//cout << "conmat5" << endl;
	if(sc->flag_el==1) CreateGapMatrix();
}


void Connections::Set_ki_Prob(double mean_prob, double relative_deviation, int post_pop,int pre_pop) {
	//ki_dist = Normaldev_BM(mean_prob, mean_prob*relative_deviation,seed+time(NULL));
	for (int ki_ind=0; ki_ind < sc->N[post_pop + 1]; ki_ind++) {
		inDegree_prob[post_pop][pre_pop][ki_ind] = mean_prob+ mean_prob*relative_deviation*randGen->normaldist.dev();
		if ((inDegree_prob[post_pop][pre_pop][ki_ind] < 0)||(inDegree_prob[post_pop][pre_pop][ki_ind] > 1)){
			ki_ind--;//special case choice
		}
	}
}

void Connections::Set_ki_Prob_FULLCOR(double mean_prob_it, double mean_prob_ii, double sigit, double sigii) {
	if(pop_num!=1){
		printf("error: Set_ki_Prob_FULLCOR only works with pop_num==1, line %d file %s",__LINE__,__FILE__);
		throw "error";
	}
	if (1== 1) {
		printf("FULLCOR is deprecated (time(NULL) for seed) pop_num==1, line %d file %s", __LINE__, __FILE__);
		throw "error";
	}
  
	ki_dist = Normaldev_BM(mean_prob_it, mean_prob_it*sigit, time(NULL));///FIXXXXX SEEDDDDD
	for (int ki_ind=0; ki_ind < sc->N[1]; ki_ind++) {
		inDegree_prob[0][0][ki_ind] = ki_dist.dev();
		inDegree_prob[0][1][ki_ind]=inDegree_prob[0][0][ki_ind]*mean_prob_ii/mean_prob_it;
		if ((inDegree_prob[0][0][ki_ind] < 0)||(inDegree_prob[0][0][ki_ind] > 1)){
			ki_ind--;//special case choice
		}
	}
}



void Connections::partial_matrix_fill(int post_pop_ind, int post_num_neurons, int pre_pop_ind,
	int pre_num_neurons/*, double connection_p */, Random *randGen) {
	double rand_x;
	int* buffer = new int[post_num_neurons];
	for (int post_n_ind = 0; post_n_ind < post_num_neurons; post_n_ind++) {
		inDegrees_K[post_pop_ind][pre_pop_ind][post_n_ind] = 0;
	}
	for (int pre_ind = 0; pre_ind < pre_num_neurons; pre_ind++) {
		int count_connections = 0;
		for (int post_ind = 0; post_ind < post_num_neurons; post_ind++) {
			rand_x = randGen->RandomUniform0_to_1();
			//randVals[post_pop_ind][pre_pop_ind][pre_ind][post_ind] = rand_x;
			if (rand_x < inDegree_prob[post_pop_ind][pre_pop_ind][post_ind]) {			
				buffer[count_connections] = post_ind;
				count_connections++;
				inDegrees_K[post_pop_ind][pre_pop_ind][post_ind]++; //indegrees probability distribution
			}
		}
		
		connections[post_pop_ind ][pre_pop_ind][pre_ind] = new int[count_connections+1];
		delay_mat[post_pop_ind ][pre_pop_ind][pre_ind] = new double[count_connections];

		for (int e_ind_toMat = 0; e_ind_toMat < count_connections; e_ind_toMat++) {
			connections[post_pop_ind ][pre_pop_ind][pre_ind][e_ind_toMat] = buffer[e_ind_toMat];

			if (sc->d_tau[post_pop_ind][pre_pop_ind] != 0) {
				printf("need to uncomment code dealing with delay variability, line %d file %s", __LINE__, __FILE__);
			}
			//double rand_del = sc->d_tau[post_pop_ind][pre_pop_ind] * (2 * randGen->RandomUniform0_to_1() - 1); //for DELAY VARIABILITY
			delay_mat[post_pop_ind][pre_pop_ind][pre_ind][e_ind_toMat] = 0;// rand_del;
		}
		connections[post_pop_ind][pre_pop_ind][pre_ind][count_connections] = -1;
	}
		delete[] buffer;
}


void Connections::clean_connection_matrix() {
	for (int post_pop = pop_num - 1; post_pop >= 0; post_pop--) {
		for (int pre_pop = pop_num; pre_pop >= 0; pre_pop--) {
			for (int pre_indNeur = 0; pre_indNeur < sc->N[pre_pop]; pre_indNeur++) {
				delete[] connections[post_pop][pre_pop][pre_indNeur];
				delete[] delay_mat[post_pop][pre_pop][pre_indNeur];
			}
		}
	}
	if (sc->flag_el == 1) {
		for (int i_n = 0; i_n < sc->N[pop_num]; i_n++) {
			delete[] gap_mat[i_n];
		}
	}
}

int Connections::num_connections(int post_pop_num, int pre_pop_num, int pre_ind) {
	int count_con = 0;
	while (connections[post_pop_num][pre_pop_num][pre_ind][count_con] != -1) count_con++;
	return count_con;
}

void Connections::Print_inDegrees(string path) {
	ofstream inDegrees_wr;
	inDegrees_wr.open(path);
	for (int inh_n = 0; inh_n < sc->N[pop_num]; inh_n++) {
		inDegrees_wr << inh_n << " " << inDegrees_K[0][0][inh_n] << " " << inDegrees_K[0][1][inh_n]<<" "<< 
			inDegree_prob[0][0][inh_n] << " " <<inDegree_prob[0][1][inh_n] << endl;
	}


}

void Connections::CreateGapMatrix() {
	for (int i_n = 0; i_n < sc->N[pop_num]; i_n++) {
		for (int j_n = 0; j_n <= i_n; j_n++) {
			double x_rand = randGen->RandomUniform0_to_1();
			double p = double(sc->Kel) / sc->N[pop_num];
			if (x_rand < p) {
				gap_prematrix[i_n][j_n] = 1;
			}
			else {
				gap_prematrix[i_n][j_n] = 0;
			}
		}
	}
	for (int i_n = 0; i_n < sc->N[pop_num]; i_n++) {
		for (int j_n = i_n + 1; j_n <sc->N[pop_num]; j_n++) {
			gap_prematrix[i_n][j_n] = gap_prematrix[j_n][i_n];
		}
	}
	int con_count;
	int *buffer = new int[sc->N[pop_num]];
	for (int i_n = 0; i_n < sc->N[pop_num]; i_n++) {
		con_count = 0;	
		for (int j_n = 0; j_n <sc->N[pop_num]; j_n++) {
			if (i_n == j_n) continue;
			if (gap_prematrix[i_n][j_n] == 1) {
				buffer[con_count] = j_n;
				con_count++;
			}
		}
		gap_con_num[i_n] = con_count;
		//cout << "gap_con_num[i_n]: " << gap_con_num[i_n] << endl;
		gap_mat[i_n] = new int[con_count];
		for (int con_ind = 0; con_ind < con_count; con_ind++) {
			gap_mat[i_n][con_ind] = buffer[con_ind];
		}
	}
	delete[] buffer;
	
}


Connections::~Connections()
{
	if (connections != nullptr) {
		//**************************TREE ARRAY VERSION*******************
		
		for (int post_pop = pop_num-1; post_pop >= 0; post_pop--) {
			for (int pre_pop =pop_num; pre_pop >=0; pre_pop--) {
				for (int pre_indNeur = 0; pre_indNeur < sc->N[pre_pop]; pre_indNeur++) {
				
							delete[] connections[post_pop][pre_pop][pre_indNeur];
							delete[] delay_mat[post_pop][pre_pop][pre_indNeur];
							//delete[] randVals[post_pop][pre_pop][pre_indNeur];
				}
				delete[] connections[post_pop][pre_pop];
				delete[] delay_mat[post_pop][pre_pop];
				delete[] inDegree_prob[post_pop][pre_pop];
				delete[] inDegrees_K[post_pop][pre_pop];
				//delete[] randVals[post_pop][pre_pop];
			}
			delete[] connections[post_pop];
			delete[] delay_mat[post_pop];
			delete[] inDegree_prob[post_pop];
			delete[] inDegrees_K[post_pop];
			//delete[] randVals[post_pop];
		}
		delete[] connections;
		delete[] delay_mat;
		delete[] inDegree_prob;
		delete[] inDegrees_K;
		//delete[] randVals;

		for (int i_n = 0; i_n < sc->N[pop_num]; i_n++) {
			delete[] gap_prematrix[i_n];
		}
		delete[] gap_prematrix;
		delete[] gap_mat;
		delete[] gap_con_num;
		//**************************************************************

	}
		
	
}



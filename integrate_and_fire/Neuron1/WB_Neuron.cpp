#include "WB_Neuron.h"
using namespace std;

#ifdef __unix__
int main(int argc, char*argv[]) {
//throw "error";
#elif _WIN32
int main() {
#endif
	//string filename;
	//string filepath;
	string temp;
	string IO_dir_path;
	//string clearbuffer;
	//char filetype;
	string sep(dir_sep);/// *** unix ***********************
#ifdef __unix__

cout<<"main1"<<endl;
	string run_name(argv[1]);

	//double size_fact = atof(argv[2]);
	double size_fact=1;

#elif _WIN32
	string run_name;
	//int seed;
	
	cout << "pop_num: " << pop_num << endl;
	cout << "run_name: " << endl;
	cin >> run_name;
	//cout << "seed: " << endl;
	//cin >> seed;
	double size_fact = 1;
#endif
	//cout << "enter run name: \n";
	//cin >> run_name;







	//Run_parameters * runpar = new Run_parameters();
	SystemConstants *sc = new SystemConstants();
	NeuronReader *Reader = new NeuronReader(sc);
	sc->size_factor = size_fact;
	//filepath = ".."+ str_sep+"NeuronIO" + str_sep; //problematic *** unix
/*#ifdef _WIN32
	IO_dir_path = SC->win_IO_path;
#elif __unix__
	IO_dir_path = SC->linux_IO_path;
#endif*/


	ifstream paramFile;
	//ofstream parseFile;
	//ifstream parsedFile;

//----

	paramFile.open(".." + sep + "NeuronIO" + sep + run_name + "_in.txt");
	if (!paramFile.is_open()) {
		paramFile.open(".." + sep + "NeuronIO" + sep + run_name + ".in");//this should become a comment
		if (!paramFile.is_open()) {
			printf("file is not open (file %s line %d)\n", __FILE__, __LINE__);
			throw "error";
		}
	}
	/*parseFile.open(".." + sep + "NeuronIO" + sep + run_name + "_parsed.in.txt"); //paramFile replace with parseFile 12.6.17
	if (!parseFile.is_open()) {
		printf("file is not open (file %s line %d)\n", __FILE__, __LINE__);
		throw "error";
	}*/

	//Reader->parse(paramFile, parseFile);
	Reader->initialize(paramFile, size_fact); //***********SC is initialized!***********************************************************
	if (pop_num != sc->popn) {
		cout << "number of populations in file does not match the chosen program" << endl;
		throw "error";
	}
	cout << "running parameter: " << sc->running_parameter<<endl;

	if ((sc->modelType.compare( "IAFv")==0)|| (sc->modelType.compare("IAFi") == 0)|| (sc->modelType.compare("IAFampav") == 0)) {
		sc->eq_num = 1;
	}
	else if (strcmp(sc->modelType.c_str(), "WB") == 0) {
		sc->eq_num = 4;
	}
	else {
		printf("problem with modelType (file %s line %d)\n", __FILE__, __LINE__);
		throw "error";
	}
	IO_dir_path = sc->IO_path;
	cout << "main: directory - " << IO_dir_path << "\n";

	cout << "method: " << sc->method << "\n";
	cout << "thalamic_input: " << sc->thalamic_input << "\n";
	//cout << "A_T is " << SC->A_T << "\n";
	cout << "runtime: " << sc->runtime << ", transient: " << sc->transient << "\n";
	//cout << "main3" << endl;
	Random *randMain = new Random(sc->seed);
	//cout << "main31" << endl;
	PoissonGen* pMain = new PoissonGen(randMain, sc);
	//cout << "main32" << endl;
	NeuronMatrix* nMat = new NeuronMatrix(pMain);

	Connections* con_mat = new Connections(randMain, sc);//needs to be set with the right parameters in SC after file reading


	con_mat->Print_inDegrees(sc->matlab_path + sep + "inDegrees2.txt");
	ODE_file * ode_file = new ODE_file(sc);

	IntegrationMethod * IM = new IntegrationMethod(nMat, con_mat, ode_file);

	ofstream elapsed_times;
	//cout << "G_IT: " << SC->G_ab[0][0] << ", G_II: " << SC->G_ab[0][1] << endl;
	//double x_0 = 0;
	//double x_end =SC->runtime;

	//ofstream paramCopy;
	///take care of IO_directory---------------------------

	string dirName = IO_dir_path + sep + run_name +/*"F"+to_string((long double)SC->size_factor)+*/ "_output";
	set_outpath(sc, run_name, "");
	//cout << "main4" << endl;
	/*ofstream seedlog;
	seedlog.open(IO_dir_path + sep + "seedlog.txt", std::ios_base::app);
	seedlog << "run_name: " << run_name << ", time: " << time(NULL) << ", seed: " << seed << endl;
	seedlog.close();*/

	if ((sc->p_rasters_firstrun == 1) || (sc->p_trace_firstrun == 1)||(sc->p_inputcurrs_firstrun == 1)) {
#ifdef _WIN32
		_mkdir(dirName.c_str());
#elif __unix__
		mkdir(dirName.c_str(), 0700); // problematic *** unix 14.12 
#endif
		/*paramCopy.open(dirName + sep + run_name+sep+".in.txt");// problematic *** unix 14.12 
		Reader->copyParamFile(paramFile, paramCopy);
		paramCopy.close();*/
	}
	else {
		//paramCopy.open(".."+str_sep+"NeuronIO"+str_sep+"PARAM " + string(1, th_input) + " " +// problematic *** unix 14.12
		//string(1, spike_stream) + " " + to_string((long double)sc->runtime) + " " +
		//to_string((long double)sc->timestep) + " " + timestamp + ".txt");
		/*paramCopy.open(IO_dir_path + sep + run_name + ".in.txt");
		Reader->copyParamFile(paramFile, paramCopy);
		paramCopy.close();*/
	}
	//------------end of IO directory--------------------

	
	long time_bef = time(NULL);
 //**********************************************************************************************
	if (strcmp(sc->running_parameter.c_str(), "single_run")==0) {
		full_reset_and_run(IM, dirName);
	}
	else {
		parameterSweep(sc->min_val, sc->max_val, sc->par_step, sc->realiz_num, IM, run_name, dirName);
	}

	/*string sigstring;
	char th_input;
	cout << "enter sig value: " << endl;
	cin >> sigstring;
	getchar();
	cout << "enter th_input: " << endl;
	th_input=getchar();
	cout << "th_input: " << th_input << endl;

	string str_thinput = string(1, th_input);
	double sig;
	
	sig = stoi(sigstring, nullptr, 10);//wrong! sig is not an int...
	//char th_input = 's';
	sc->thalamic_input = th_input;
		sc->sig_ki[0][0] = sig;
		sc->sig_ki[0][1] = sig;
		createRepeatedHist(IM, "th_input_"+str_thinput+"sig_"+sigstring, dirName, 10);*/
	/*sc->runtime = 60000;
	sc->thalamic_input = 's';
	sc->sig_ki[0][0] = 0;
	sc->sig_ki[0][1] = 0;
	createRepeatedHist(IM, run_name+"60s_th_input_ssig_0", dirName, 1);

	sc->thalamic_input = 's';
	sc->sig_ki[0][0] = 0.2;
	sc->sig_ki[0][1] = 0.2;
	createRepeatedHist(IM, run_name + "60s_th_input_ssig_0.2", dirName, 1);

	sc->thalamic_input = 'm';
	sc->sig_ki[0][0] = 0;
	sc->sig_ki[0][1] = 0;
	createRepeatedHist(IM, run_name + "60s_th_input_msig_0", dirName, 1);

	sc->thalamic_input = 'm';
	sc->sig_ki[0][0] = 0.2;
	sc->sig_ki[0][1] = 0.2;
	createRepeatedHist(IM, run_name + "60s_th_input_msig_0.2", dirName, 1);*/

	
	//SigSweep(0, 0.3, 0.05, 1, 1, 10, IM, run_name, dirName);
	//GelSweep(0, 0.4, 0.05, 1, IM, run_name, dirName);

	//printProbsAndInDegrees(sc,  con_mat, run_name);
	
	//AtSweep(20, 20, 1, 1, IM, run_name+"rate_vs_input", dirName);
	//SigSweepAll(0, 0.3, 0.15, 5, IM, run_name + "8.6", dirName);
	//full_reset_and_run(IM, dirName);
	
  /*full_reset_and_run(IM, dirName);	
  sc->sig_ki[0][0] = 0.2;
  sc->sig_ki[0][1] = 0.2;
  set_outpath( sc,run_name , "sig_0.2");
  full_reset_and_run(IM, dirName);*/


 
	//f_kSweep(25, 25, 1, 1, rk, run_name, dirName);
	//GelSweep(0, 0.04, 0.002, 1, rk, run_name, dirName);
//**********************************************************************************************

	/*set_outpath(sc, run_name, "debug");

	randMain->count = 0;
	sc->sig_ki[0][0] = 0.2;
	sc->sig_ki[0][1] = 0.2;
	con_mat->clean_connection_matrix();
	cout << "1rand count: " << randMain->count << endl;
	con_mat->CreateConnectionMatrix(randMain);
	cout << "2rand count: " << randMain->count << endl;
	IM->Time_Evolution(dirName);
	
	con_mat->Print_inDegrees(sc->matlab_path + sep + "inDegrees2.txt");*/

	/*sc->sig_ki[0][0] = 0.2;
	sc->sig_ki[0][1] = 0.2;
	con_mat->clean_connection_matrix();
	con_mat->ReCreateConnectionMatrixSamerand(randMain);
	con_mat->Print_inDegrees(sc->matlab_path + sep + "inDegrees22.txt");*/
	
	
	
	//cout << "rand count: " << randMain->count << endl;
	/*con_mat->Set_ki_Prob(1 / 6.0, 0.2, 0, 1);
	for (int pre_ind = 0; pre_ind < sc->N[1]; pre_ind ++ ) {
		int*old_indegrees = con_mat->connections[0][1][pre_ind];
	}*/
	//full_reset_and_run(IM, dirName);
	//cout << to_string(long double(randMain->RandomUniform0_to_1())) << endl;
	
	//cout << "rand count: " << randMain->count << endl;


	
	
	paramFile.close();

	
	long time_af = time(NULL);
	long time_diff = time_af - time_bef;
	cout << "elapsed time: " << time_diff << endl;
	//elapsed_times.open(IO_dir_path+sep+"elapsed_times.txt" , std::ios_base::app);
	//elapsed_times << "time_diff: " << time_diff <<", time: " <<(time_af /3600.0) << "\n";
	//elapsed_times.close();
	
	
	delete ode_file;
	delete IM;
	delete con_mat;
	delete nMat; 
	delete pMain;
	delete sc;
	//delete runpar;
	delete randMain;

};

void parameterSweep(double min_val, double max_val, double par_step, int realiz_num, IntegrationMethod *IM, string run_name, string dirName) {
	SystemConstants *sc = IM->sc;
	set_outpath(IM->sc, run_name, sc->running_parameter.c_str());
	for (double value = min_val; value <= max_val; value += par_step) {
		for (int realization = 0; realization<realiz_num; realization++) {
			cout << sc->running_parameter << " = " << value << endl;
			set_running_parameter(sc, sc->running_parameter, value);
			full_reset_and_run(IM, dirName);

			if ((value == min_val) && (realization == 0)) {
				if (sc->p_rasters_firstrun == 1) cout << "rasters printed"<<endl;
				if (sc->p_trace_firstrun == 1) cout << "trace printed" << endl;
				if (sc->p_inputcurrs_firstrun == 1) cout << "inputcurrs printed" << endl;
				if (sc->p_spikes_per_touch_fr == 1) cout << "spikes_per_touch printed" << endl;
				sc->p_rasters_firstrun = 0;
				sc->p_trace_firstrun = 0;
				sc->p_inputcurrs_firstrun = 0;
				sc->p_spikes_per_touch_fr = 0;
				/*for (int printInt = rasters; printInt != spikes_per_touch+1; printInt++)
				{
					Print print = static_cast<Print>(printInt);
					cout << "enum trial!" << endl;
				}*/
			}
		}
	}
}

void set_running_parameter(SystemConstants * sc,string running_parameter, double value) {
	sc->run_param_val = value;

	if (strcmp(running_parameter.c_str(), "A_T") == 0) { sc->A_T = value/1000; }
	else if (strcmp(running_parameter.c_str(), "sig_all") == 0) {
		for (int post_pop = 0; post_pop < pop_num; post_pop++)
			for (int pre_pop = 0; pre_pop < pop_num + 1; pre_pop++)
				sc->sig_ki[post_pop][pre_pop] = value;
	}
	else if (strcmp(running_parameter.c_str(), "K_all") == 0) {
		for (int post_pop = 0; post_pop < pop_num; post_pop++) {
			for (int pre_pop = 0; pre_pop < pop_num + 1; pre_pop++) {
				sc->K_ab[post_pop][pre_pop] = (value / 25)*sc->K_norm[post_pop][pre_pop];
				cout << "K_ab[" << post_pop << "][" << pre_pop << "]: " << sc->K_ab[post_pop][pre_pop] << "  ";
				if (sc->K_ab[post_pop][pre_pop]>sc->N[pre_pop]) {
					cout << "error: K_ab is bigger than N_b" << endl;
					throw "error";
				}
			}
		}
		cout << endl;
	}
	else {
		cout << "run has not yet been defined for parameter " << running_parameter << endl;
		throw "error";
	}

}

void f_kSweep(int kii_min, int kii_max, int step, int num_realiz, IntegrationMethod *IM, string run_name, string dirName){
	SystemConstants *sc = IM->sc;
	set_outpath(IM->sc, run_name, "f_K");
	for (double kii = kii_min; kii <= kii_max; kii += step) {
		for (int realization = 0; realization<num_realiz; realization++) {			
			sc->K_ab[0][0] = kii * 3;
			sc->K_ab[0][1] = kii;
			cout << "f_K = " << double(kii)/25 << endl;
			full_reset_and_run(IM, dirName);
		}
	}	
}

void AtSweep(double at_min, double at_max, double step, int num_realiz, IntegrationMethod *IM, string run_name, string dirName) {
	SystemConstants *sc = IM->sc;
	set_outpath(IM->sc, run_name, "AtSweep");
	for (double at = at_min; at <= at_max; at += step) {
		for (int realization = 0; realization<num_realiz; realization++) {
			sc->A_T = at/1000;
			cout << "At = " << at<< endl;
			full_reset_and_run(IM, dirName);
		}
	}
}

void GelSweep(double gel_min, double gel_max, double step, int num_realiz, IntegrationMethod *IM, string run_name, string dirName) {
	SystemConstants *sc = IM->sc;
	set_outpath(sc, run_name, "Gel");
	int tempflag = sc->flag_el;
	//sc->flag_el = 1; //Not good. the matrix is only allocated value the beginning if flag==1.

	for (double gel = gel_min; gel <= gel_max; gel += step) {
		for (int realization = 0; realization<num_realiz; realization++) {
			sc->Gel = gel;
			cout << "gel is: " << gel << endl;
			full_reset_and_run(IM,  dirName);
		}
	}
	//sc->flag_el = tempflag;
}

void SigSweepAll(double sig_min, double sig_max, double step, int num_realiz, IntegrationMethod *IM, string run_name, string dirName) {
	SystemConstants *sc = IM->sc;
	
	for (double sig = sig_min; sig <= sig_max; sig += step) {
		for (int realization = 0; realization < num_realiz; realization++) {
			set_outpath(sc, run_name, "sig_all");
			for (int post_pop = 0; post_pop < pop_num; post_pop++) {
				for (int pre_pop = 0; pre_pop < pop_num + 1; pre_pop++){
					sc->sig_ki[post_pop][pre_pop] = sig;
				}
			}
			cout << "sig_all: " << sig << endl;
			/*if ((flag_it == 1) && (flag_ii == 0)) {
				set_outpath(sc, run_name, "sig_it");
				sc->sig_ki[0][0] = sig;
				cout << "sig_it: " << sig << endl;
			}
			else if ((flag_it == 0) && (flag_ii == 1)) {
				set_outpath(sc, run_name, "sig_ii");
				sc->sig_ki[0][1] = sig;
				cout << "sig_ii: " << sig << endl;
			}
			else if ((flag_it == 1) && (flag_ii == 1)) {
				set_outpath(sc, run_name, "sig_all");
				sc->sig_ki[0][0] = sig;
				sc->sig_ki[0][1] = sig;
				cout << "sig_all: " << sig << endl;
			}
			else {
				printf("error: value least 1 heterogeneity sig must be varied, line %d file %s", __LINE__, __FILE__);
				throw "error";
			}*/
			
			full_reset_and_run(IM, dirName);
		}
	}
}

void full_reset_and_run(IntegrationMethod *IM, string dirName) {
	Connections * con_mat = IM->p_conmat;
	Random * randMain = IM->p_rand_Main;
	NeuronMatrix * nMat = IM->p_n_Mat;
	//cout << "frar 1st rand count: " << randMain->count << endl;
	con_mat->clean_connection_matrix();
	con_mat->CreateConnectionMatrix(randMain);
	//cout << "frar 2st rand count: " << randMain->count << endl;
	IM->Time_Evolution( dirName);
	//cout << "frar 3st rand count: " << randMain->count << endl;
	nMat->reset_neuronMatrix();
}

void createRepeatedHist(IntegrationMethod *IM, string run_name, string dirName, int num_realiz) {
	
	SystemConstants *sc = IM->sc;
	set_outpath(IM->sc, run_name, "mHist");
	for (int realization = 0; realization<num_realiz; realization++) {
			full_reset_and_run(IM, dirName);
			if (realization == 0) IM->flag_repeatHist = 1;
	}
	for (int histInd = 0; histInd < 26; histInd++) {
		IM->logHist[1][histInd] = IM->logHist[1][histInd] / num_realiz;
	}

	IM->PrintHist();
	IM->flag_repeatHist = 0;

}

void set_outpath(SystemConstants * sc,string run_name, string addition) {
#ifdef _WIN32
	sc->outfilepath = sc->matlab_path + run_name + addition;
#elif __unix__
	sc->outfilepath = sc->IO_path + run_name + addition;
#endif
}

/*void printRunningParToSync(SystemConstants *sc,double runningParameter) {//this is not efficient. runningparameter should be an input to TimeEvolution or in sc/parametersfile
	ofstream synchrony_wr;
	synchrony_wr.open(sc->outfilepath + "_sync.txt", std::ios_base::app);
	synchrony_wr << setprecision(6) << fixed << runningParameter << " ";
	synchrony_wr.close();
}*/

void printProbsAndInDegrees(SystemConstants * sc,Connections* con, string run_name) {
	set_outpath(sc, run_name, "Indegrees");
	ofstream indegrees_wr;
	indegrees_wr.open(sc->outfilepath+ ".txt");
	for (int j_n = 0; j_n < sc->N[1]; j_n++) {
		indegrees_wr << setprecision(6) << fixed << con->inDegree_prob[0][0][j_n] <<
			" " << con->inDegree_prob[0][1][j_n] << " " << con->inDegrees_K[0][0][j_n] << " "
			<< con->inDegrees_K[0][1][j_n] << endl;
	}
	indegrees_wr.close();

}


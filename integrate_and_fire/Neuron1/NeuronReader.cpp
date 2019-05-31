#include "NeuronReader.h"



NeuronReader::NeuronReader()
{
	throw "don't use this constructor";
}

NeuronReader::NeuronReader(SystemConstants* SC/*, Run_parameters *RP*/)
{
	sc = SC;
	//runpar = RP;
}

/*void NeuronReader::parse(ifstream & parametersFile, ofstream & parseFile) {
	string line;
	while (getline(parametersFile, line)) {

		stringstream contents(line);
		
		if (contents.rdbuf()->in_avail() == 0) { //skip empty lines
			continue;
		}
		if (line.compare("\r") == 0) {//check for a unix conversion error
			printf("windows to unix conversion error: not converted. file %s line %d\n", __FILE__, __LINE__);
			throw "error";
		}
		if (line.at(0) == '#') { //skip comments
			continue;
		}

		if (line.at(0) == '$') { // execute parsing commands
			line = line.substr(1);
			continue;
		}

		if (!(contents >> param_name >> param_val)) {
			//printf( "error reading line from file.\nline content: %s\nfile %s line %d\n",line,__FILE__,__LINE__);
			printf("error reading line from file. position: file %s line %d\n", __FILE__, __LINE__);
			//problematic *** unix
			throw "error reading file";
		}
		//	cout << param_name <<" "<< param_val<<"\n";// - if I want to cout the parameters read from the file for debugging****


		if (param_name.compare("modelType") == 0) sc->modelType = param_val;
	}
}*/


void NeuronReader::initialize(ifstream & parametersFile, double size_fact)
{
  //cout<<"NRp2: size_fact = "<<size_fact<<endl;
	stringstream char_content;
	string line;
	string param_name;
	string temp;
	string temp_win_path, temp_lin_path;
	string endline;
	double param_val;
	int count = 1;
	char finalchar;

	getline(parametersFile, line);
	char_content.str(line);
	char_content >> temp >> sc->method;

	getline(parametersFile, line);
	char_content.str(line);
	char_content >> temp >> sc->thalamic_input;

	//strings:-----
	getline(parametersFile, line);
	char_content.str(line);
	char_content >> temp >> sc->modelType >> finalchar;

	getline(parametersFile, line);
	char_content.str(line);
	char_content >> temp >> temp_win_path>>finalchar;

	getline(parametersFile, line);
	char_content.str(line);
	char_content >> temp >> temp_lin_path>>finalchar;


	getline(parametersFile, line);
	char_content.str(line);
	char_content >> temp>> sc->matlab_path>>finalchar;

	getline(parametersFile, line);
	char_content.str(line);
	char_content >> temp >> sc->running_parameter >> finalchar;

#ifdef _WIN32
	sc->IO_path = temp_win_path;
#elif __unix__
	sc->IO_path = temp_lin_path;
#endif

	
	//cout << "NeuronReader: IO_path: " << sc->IO_path << endl << "matlab_path: " << sc->matlab_path << endl;


	while (getline(parametersFile, line)) {
  
		stringstream contents(line);
/*#ifdef _WIN32
	if (contents.rdbuf()->in_avail() == 0) { //problematic *** unix
#elif __unix__
  if (contents.rdbuf()->in_avail() == 0) { //problematic *** unix
	//if (line.compare("\r")==0) {//problematic *** unix
#endif*/

    if (contents.rdbuf()->in_avail() == 0) { //problematic *** unix
			//***cout << "--\n";// - if I want to cout the parameters read from the file for debugging****
			continue;
		}
   if (line.compare("\r")==0) {//problematic *** unix
     printf( "windows to unix conversion error: not converted. file %s line %d\n",__FILE__,__LINE__);
     throw "error";
   }
   //cout << "first char: " << line.at(0) << endl;
   //cout << "line: " << line << endl;
   if (line.at(0) == '#') { //skip comments
	   contents.str("");
	   continue;
   }
		if (!(contents >> param_name >> param_val)) {
			//printf( "error reading line from file.\nline content: %s\nfile %s line %d\n",line,__FILE__,__LINE__);
			printf( "error reading line from file. position: file %s line %d\n",__FILE__,__LINE__);
			//problematic *** unix
			throw "error reading file";
		}
		//cout << param_name <<" "<< param_val<<"\n";// - if I want to cout the parameters read from the file for debugging****
	

		if (param_name.compare("modelType") == 0) sc->modelType = param_val;// duplicate? 12.6.17 should not be doing anything

		if (param_name.compare("seed") == 0) sc->seed = param_val;
		
		if (param_name.compare("pop_num") == 0) sc->popn = param_val;
		if (param_name.compare("runtime") == 0) sc->runtime = param_val;
		if (param_name.compare("transient") == 0) sc->transient = param_val;
		//if (param_name.compare("printflag") == 0) sc->printflag = param_val;

		//--------------running parameters-----------
		if (param_name.compare("min_val") == 0) sc->min_val = param_val;
		if (param_name.compare("max_val") == 0) sc->max_val = param_val;
		if (param_name.compare("par_step") == 0) sc->par_step = param_val;
		if (param_name.compare("realiz_num") == 0) sc->realiz_num = param_val;

		//sc->n_ex_neurons = 0; //this happens in sc initialization
		for (int ex_neur_num = 0; ex_neur_num < 10; ex_neur_num++) {
			if (param_name.compare("ex" + to_string((long double)(ex_neur_num + 1)) + "pop") == 0) {
				sc->ex_neurons[ex_neur_num][0] = param_val;
				sc->n_ex_neurons++;//in case we need to use the number of example neurons
			}
			if (param_name.compare("ex" + to_string((long double)(ex_neur_num + 1)) + "neur") == 0) sc->ex_neurons[ex_neur_num][1] = param_val;
		}
		if (sc->n_ex_neurons > 10) {
			cout << "error: sc->n_ex_neurons>10" << endl;
			throw "error";
		}

		if (param_name.compare("p_rasters_firstrun") == 0) sc->p_rasters_firstrun = param_val;
		if (param_name.compare("p_trace_firstrun") == 0) sc->p_trace_firstrun = param_val;
		if (param_name.compare("p_mean_trace_fr") == 0) sc->p_mean_trace_fr = param_val;
		if (param_name.compare("p_inputcurrs_firstrun") == 0) sc->p_inputcurrs_firstrun = param_val; 
		if (param_name.compare("p_spikes_per_touch_fr") == 0) sc->p_spikes_per_touch_fr = param_val;

		
		if (param_name.compare("p_sync_fig1") == 0) sc->p_sync_fig1 = param_val;
		if (param_name.compare("p_mean_rates") == 0) sc->p_mean_rates = param_val;
		if (param_name.compare("p_mean_CV") == 0) sc->p_mean_CV = param_val;
		if (param_name.compare("p_chi") == 0) sc->p_chi = param_val;
		if (param_name.compare("p_quiescent_neurs") == 0) sc->p_quiescent_neurs = param_val;
		if (param_name.compare("p_dev_rates") == 0) sc->p_dev_rates = param_val;
		


		//-----------end of running parameters-------

		if (param_name.compare("timestep") == 0) sc->timestep = param_val;

		if (param_name.compare("A_T") == 0) sc->A_T = param_val;
		if (param_name.compare("B_T") == 0) sc->B_T = param_val;
		if (param_name.compare("C_T") == 0) sc->C_T = param_val;

		

		if (param_name.compare("gL_E") == 0) sc->gL[0] = param_val;
		if (param_name.compare("gL_I") == 0) sc->gL[1] = param_val;
		if (param_name.compare("gKz_E") == 0) sc->gKz[0] = param_val;
		if (param_name.compare("gKz_I") == 0) sc->gKz[1] = param_val;

		if (param_name.compare("T") == 0) sc->T = param_val;
		if (param_name.compare("tau_c") == 0) sc->tau_c = param_val;
		if (param_name.compare("t_c") == 0) sc->t_c = param_val;
		if (param_name.compare("phi") == 0) sc->phi = param_val;

		if (param_name.compare("Iapp") == 0) sc->Iapp = param_val;
		if (param_name.compare("phi2") == 0) sc->phi2 = param_val;

		if (param_name.compare("flag_el") == 0) sc->flag_el = param_val;
		if (param_name.compare("Kel") == 0) sc->Kel = param_val;
		if (param_name.compare("Gel") == 0) sc->Gel = param_val;

#ifdef _WIN32    
		if (param_name.compare("win_IO_path") == 0) sc->IO_path = param_val;
#elif __unix__
		if (param_name.compare("linux_IO_path") == 0) sc->IO_path = param_val;
#endif
		if (param_name.compare("matlab_path") == 0) sc->matlab_path = param_val;


#if pop_num==2 //inhibition change

		if (param_name.compare("sig_ki_ET") == 0) sc->sig_ki[0][0] = param_val;  
		if (param_name.compare("sig_ki_EE") == 0) sc->sig_ki[0][1] = param_val;
		if (param_name.compare("sig_ki_EI") == 0) sc->sig_ki[0][2] = param_val;  
		if (param_name.compare("sig_ki_IT") == 0) sc->sig_ki[1][0] = param_val;
		if (param_name.compare("sig_ki_IE") == 0) sc->sig_ki[1][1] = param_val;  
		if (param_name.compare("sig_ki_II") == 0) sc->sig_ki[1][2] = param_val;

		if (param_name.compare("N_T") == 0) sc->N[0] = param_val;
		if (param_name.compare("N_E") == 0) sc->N[1] = param_val;
		if (param_name.compare("N_I") == 0) sc->N[2] = param_val;

		if (param_name.compare("tau_decay_T") == 0) sc->tau_decay[0] = param_val;
		if (param_name.compare("tau_decay_E") == 0) sc->tau_decay[1] = param_val;
		if (param_name.compare("tau_decay_I") == 0) sc->tau_decay[2] = param_val;

		if (param_name.compare("K_ET") == 0) { sc->K_ab[0][0] = param_val; sc->K_norm[0][0] = param_val; }
		if (param_name.compare("K_EE") == 0) {sc->K_ab[0][1] = param_val; sc->K_norm[0][1] = param_val; }
		if (param_name.compare("K_EI") == 0) {sc->K_ab[0][2] = param_val; sc->K_norm[0][2] = param_val; }
		if (param_name.compare("K_IT") == 0) {sc->K_ab[1][0] = param_val; sc->K_norm[1][0] = param_val; }
		if (param_name.compare("K_IE") == 0) {sc->K_ab[1][1] = param_val; sc->K_norm[1][1] = param_val; }
		if (param_name.compare("K_II") == 0) {sc->K_ab[1][2] = param_val; sc->K_norm[1][2] = param_val; }
		

		if (param_name.compare("G_ET") == 0) sc->G_ab[0][0] = param_val;
		if (param_name.compare("G_EE") == 0) sc->G_ab[0][1] = param_val;
		if (param_name.compare("G_EI") == 0) sc->G_ab[0][2] = param_val;
		if (param_name.compare("G_IT") == 0) sc->G_ab[1][0] = param_val;
		if (param_name.compare("G_IE") == 0) sc->G_ab[1][1] = param_val;
		if (param_name.compare("G_II") == 0) sc->G_ab[1][2] = param_val;

		if (param_name.compare("tau_delay_ET") == 0) sc->tau_delay[0][0] = param_val;
		if (param_name.compare("tau_delay_EE") == 0) sc->tau_delay[0][1] = param_val;
		if (param_name.compare("tau_delay_EI") == 0) sc->tau_delay[0][2] = param_val;
		if (param_name.compare("tau_delay_IT") == 0) sc->tau_delay[1][0] = param_val;
		if (param_name.compare("tau_delay_IE") == 0) sc->tau_delay[1][1] = param_val;
		if (param_name.compare("tau_delay_II") == 0) sc->tau_delay[1][2] = param_val;

		if (param_name.compare("d_tau_ET") == 0) sc->d_tau[0][0] = param_val;
		if (param_name.compare("d_tau_EE") == 0) sc->d_tau[0][1] = param_val;
		if (param_name.compare("d_tau_EI") == 0) sc->d_tau[0][2] = param_val;
		if (param_name.compare("d_tau_IT") == 0) sc->d_tau[1][0] = param_val;
		if (param_name.compare("d_tau_IE") == 0) sc->d_tau[1][1] = param_val;
		if (param_name.compare("d_tau_II") == 0) sc->d_tau[1][2] = param_val;

#elif pop_num==1
		if (param_name.compare("sig_ki_IT") == 0) sc->sig_ki[0][0] = param_val;  //need to generalize to 2 populations?
		if (param_name.compare("sig_ki_II") == 0) sc->sig_ki[0][1] = param_val;

		if (param_name.compare("N_T") == 0) sc->N[0] = (int)(param_val*size_fact);
		if (param_name.compare("N_I") == 0) sc->N[1] = (int)(param_val*size_fact);

		if (param_name.compare("tau_decay_T") == 0) sc->tau_decay[0] = param_val;
		if (param_name.compare("tau_decay_I") == 0) sc->tau_decay[1] = param_val;


		if (param_name.compare("K_IT") == 0) { sc->K_ab[0][0] = (int)(param_val*size_fact); sc->K_norm[0][0] = (int)(param_val*size_fact); }
		if (param_name.compare("K_II") == 0) { sc->K_ab[0][1] = (int)(param_val*size_fact); sc->K_norm[0][1] = (int)(param_val*size_fact); }

		if (param_name.compare("G_IT") == 0) sc->G_ab[0][0] = param_val;
		if (param_name.compare("G_II") == 0) sc->G_ab[0][1] = param_val;

		if (param_name.compare("tau_delay_IT") == 0) sc->tau_delay[0][0] = param_val;
		if (param_name.compare("tau_delay_II") == 0) sc->tau_delay[0][1] = param_val;

		if (param_name.compare("d_tau_IT") == 0) sc->d_tau[0][0] = param_val;
		if (param_name.compare("d_tau_II") == 0) sc->d_tau[0][1] = param_val;
#endif

	}


	sc->setBuffer_size_andDelay();
	sc->setSpikevals_and_matprob();//possibly redundant: appears also in CreateConnectionMatrix. 29.12.16
	

	parametersFile.clear();
	parametersFile.seekg(0, ios::beg);
	
}

void NeuronReader::copyParamFile(ifstream& fileToRead, ofstream& fileToWrite) {
	string line;
	string temp;

	while (getline(fileToRead, line)) {
		fileToWrite << line << endl;
	}


	fileToRead.clear();
	fileToRead.seekg(0, ios::beg);
	
}



NeuronReader::~NeuronReader()
{
}

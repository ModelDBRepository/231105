/* This program simulates a network for one parameter set. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "tca.h"
#include "tcn.h"
#include "tcwb.h"
#include "ippn.h"
#include "fp.h"

/* This function simulates the system for one parameter set. */
void one_par(par_all *parall, avr_val *av, fl_st fl)
{
  func_cell_model *fcm;
  net_par netpar;
  run_par runpar;
  syn_str synstr;
  spk_str spkstr;
  double ***Varbar;

  read_npop(&netpar, av, fl);

  fcm = define_fcm_array(1, netpar.npop);
  read_input(fcm, &netpar, &runpar, parall, fl);
  define_variables_arrays(&netpar, &runpar, &Varbar, fl);
  compute_heterogeneous_intrinsic_parameters(&netpar, &runpar, parall, fl);
  compute_opto_conductance(&netpar, &runpar, fl);
  substitute_connectivity(&netpar, &runpar, parall, fl);
  in_con(fcm, Varbar, &netpar, &runpar, parall, fl);
  define_synaptic_variables(&synstr, &netpar, &runpar, fl);
  initiate_synaptic_strengths_and_variables(&synstr, &netpar, &runpar, fl);
  initiate_electrical_strengths(&synstr, &netpar, &runpar, fl);
  define_spike_struct(&spkstr, &netpar, &runpar, fl);
  determine_thal_Cv(&netpar, &runpar, parall, fl);
  if (netpar.T.thal_input == 'p')
    initialize_thalamic_variables_and_spikes(&spkstr, &netpar, &runpar, parall,
    fl);
  n_run(fcm, Varbar, &synstr, &spkstr, &netpar, &runpar, parall, av,
    fl);
  compute_spike_statistics(&spkstr, &netpar, &runpar, av, fl);
  compute_voltage_statistics(spkstr.Vav, &synstr, &netpar, &runpar, av, fl);

  free_fcm_array(fcm, 1, netpar.npop);
  free_nwritear_arrays(&netpar, &runpar, fl);
  free_variables_arrays(&netpar, &runpar, Varbar, fl);
  free_heterogeneous_intrinsic_parameters(&netpar, &runpar, fl);
  free_opto_arrays(&netpar, &runpar, fl);
  free_connectivity_arrays(&netpar, &runpar, fl);
  free_synaptic_variables(&synstr, &netpar, &runpar, fl);
  free_spike_struct(&spkstr, &netpar, &runpar, fl);
  free_dvector(runpar.deltat, 1, runpar.ndeltat);
  free_ivector(runpar.nt, 1, runpar.ndeltat);
  if (netpar.T.thal_input == 'p')
  {
    free_dvector(netpar.T.phi, 1, netpar.T.ntl);
  }
  free_dvector(netpar.T.Cv, 1, netpar.T.ntl);
}

/* This function reads the number of populations npop. */
void read_npop(net_par *netpar, avr_val *av, fl_st fl)
{
  int ipop;
  char fmt[5];
  
  rewind(fl.tmp);

  if (fscanf(fl.tmp, "npop=%d pop_name=", &netpar->npop) != 0)
    fprintf(fl.out, "npop=%d pop_name=", netpar->npop);
  else
    printf("npop\n");

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    if (ipop < netpar->npop)
      strcpy(fmt, "%c ");
    else
      strcpy(fmt, "%c\n");
    
    if (fscanf(fl.tmp, fmt, &netpar->pop_name[ipop]) != 0)
    fprintf(fl.out, fmt, netpar->pop_name[ipop]);
  else
    printf("npop\n");
  }

  av->npop = netpar->npop;
  strcpy(av->pop_name, netpar->pop_name);
}

/* allocates an fcm vector with subscript range v[nl..nh] */
func_cell_model *define_fcm_array(long nl, long nh)
{
	func_cell_model *v;

	v=(func_cell_model *)malloc((size_t) ((nh-nl+1+NR_END)*
          sizeof(func_cell_model)));
	if (!v) nrerror("allocation failure in func_cell_model()");
	return v-nl+NR_END;
}

/* Free a func_cell_model vector allocated with func_cell_model() */
void free_fcm_array(func_cell_model *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

/* Free nwritear */
void free_nwritear_arrays(net_par *netpar, run_par *runpar, fl_st fl)
{
  int ipop;

  for (ipop=0; ipop<=netpar->npop; ipop++)
    free_ivector(runpar->nwritear[ipop], 1, runpar->nwrite[ipop]);

  free_pivector(runpar->nwritear, 0, netpar->npop);
  free_ivector(runpar->nwrite, 0, netpar->npop);
}

/* allocate an int *vector with subscript range v[nl..nh] */
int **pivector(long nl, long nh)
{
	int **v;

	v=(int **)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int*)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

/* free an int *vector allocated with ivector() */
void free_pivector(int **v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

/* allocate an double *vector with subscript range v[nl..nh] */
double **pdvector(long nl, long nh)
{
	double **v;

	v=(double **)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double*)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

/* free an double *vector allocated with ivector() */
void free_pdvector(double **v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

/* This function reads the input parameters. */
void read_input(func_cell_model *fcm, net_par *netpar, run_par *runpar,
     par_all *parall, fl_st fl)
{
  double ftau;
  char line[Mlinea];
  int ipop, jpop, iwrite, ifscan, ideltat;

  runpar->epsilon = 1.0e-15;
  runpar->nwrite = ivector(0, netpar->npop);
  runpar->nwritear = pivector(0, netpar->npop);
  
  /* Reading the data */

  if (fscanf(fl.tmp, "noise=%lf\n", &netpar->noise) != 0)
    fprintf(fl.out, "noise=%lf\n", netpar->noise);
  else
    printf("noise\n");

  /* T cells */

  if (fgets(line, Mlinea, fl.tmp) != NULL)
    fprintf(fl.out, "%s", line);

  if (fscanf(fl.tmp, "Av=%lf Bv=%lf Tper=%lf phi_read=%lf Pi Bp=%lf tcfrac=%lf "
    "cycle tauc=%lf\n", &netpar->T.Av, &netpar->T.Bv, &netpar->T.Tper,
    &netpar->T.phi_read, &netpar->T.Bp, &netpar->T.tcfrac, &netpar->T.tauc)
    != 0)
  {
    netpar->T.Av /= 1000.0;
    fprintf(fl.out, "Av=%lf Bv=%lf Tper=%lf phi_read=%lf Bp=%lf\n",
    netpar->T.Av, netpar->T.Bv, netpar->T.Tper, netpar->T.phi_read,
    netpar->T.Bp);
    netpar->T.tc = netpar->T.tcfrac * netpar->T.Tper; 
    fprintf(fl.out, "tcfrac=%lf tc=%lf tauc=%lf\n",  netpar->T.tcfrac,
    netpar->T.tc, netpar->T.tauc);
  }
  else
    printf("Av=%lf\n", netpar->T.Av);

  if (fscanf(fl.tmp, "wd=%c Cvmin=%lf Cvmax=%lf ds_type=%c frac_only_p=%lf "
    "frac_only_r=%lf\n", &netpar->T.wd, &netpar->T.Cvmin, &netpar->T.Cvmax,
    &netpar->T.ds_type, &netpar->T.frac_only_p, &netpar->T.frac_only_r) != 0)
  {
    netpar->T.Cvmin /= 1000.0;
    netpar->T.Cvmax /= 1000.0;
    fprintf(fl.out, "wd=%c Cvmin=%lf Cvmax=%lf ds_type=%c frac_only_p=%lf "
    "frac_only_r=%lf\n", netpar->T.wd, netpar->T.Cvmin, netpar->T.Cvmax,
    netpar->T.ds_type, netpar->T.frac_only_p, netpar->T.frac_only_r);
  }
  else
    printf("wd=%c\n", netpar->T.wd);
    
  if (fscanf(fl.tmp, "Tall=%lf nspike_max=%d ntl=%d determine_phi=%c "
    "thal_input=%c\n", &netpar->T.Tall, &netpar->T.nspike_max, &netpar->T.ntl,
    &netpar->T.determine_phi, &netpar->T.thal_input) != 0)
  {
    netpar->T.Tall *= 1000.0;
    fprintf(fl.out, "Tall=%lf msec nspike_max=%d ntl=%d determine_phi=%c "
    "thal_input=%c\n", netpar->T.Tall, netpar->T.nspike_max, netpar->T.ntl,
    netpar->T.determine_phi, netpar->T.thal_input);
  }
  else
    printf("Tall\n");

  if ((netpar->T.thal_input != 'p') && (netpar->T.thal_input != 'a'))
  {
    printf("thal_input should be p or a!\n");
    exit(0);
  }

  netpar->T.AvBcTo2p = netpar->T.Av * netpar->T.Bv * netpar->T.Tper / 
    (2.0 * Pi);
  fprintf(fl.out, "AvBcTo2p=%lf\n", netpar->T.AvBcTo2p);

  netpar->T.pr = 0;

  netpar->T.Cv = dvector(1, netpar->T.ntl);
  
  if (fscanf(fl.tmp, "nw=%c Avnw=%lf Tnw=%lf Telev=%lf Tw=%lf\n",
    &netpar->T.nw, &netpar->T.Avnw, &netpar->T.Tnw, &netpar->T.Telev,
    &netpar->T.Tw) != 0)
  {
    netpar->T.Avnw /= 1000.0;
    fprintf(fl.out, "nw=%c Avnw=%lf Tnw=%lf Telev=%lf Tw=%lf\n",
    netpar->T.nw, netpar->T.Avnw, netpar->T.Tnw, netpar->T.Telev, netpar->T.Tw);
  }
  else
    printf("nw=%c\n", netpar->T.nw);

  /*
  if (fgets(line, Mlinea, fl.tmp) != NULL)
    fprintf(fl.out, "%s", line);
  */

  /* Reading data for one cell type */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    sprintf(line, "%c_CELL: %%c%%c model\n", netpar->pop_name[ipop]);

    fflush(fl.out);
    if (fscanf(fl.tmp, line, &netpar->C[ipop].model_type[0], 
    &netpar->C[ipop].model_type[1]) == 0)
    {
      printf("cannot read model_type!\n");
      exit(0);
    }
    netpar->C[ipop].model_type[2] = '\n';

    fprintf(fl.out, line, netpar->C[ipop].model_type[0], 
    netpar->C[ipop].model_type[1]);

    if (strncmp(netpar->C[ipop].model_type, "WB", 2) == 0)  
      /* Golomb-Hansel model for FS neurons */
    {
      substitute_function_WB(&fcm[ipop], fl);
    }
    /*
    else if (strncmp(netpar->C[ipop].model_type, "GY", 2) == 0)  
      (* Golomb-Yaari model for E neurons  *)
    {
      substitute_function_GY(&fcm[ipop], fl);
    }
    else if (strncmp(netpar->C[ipop].model_type, "GH", 2) == 0)  
      (* Golomb-Hansel model for FS neurons *)
    {
      substitute_function_GH(&fcm[ipop], fl);
    }
    */
    else
    {
      printf("wrong model_type=%s", netpar->C[ipop].model_type);
      exit(0);
    }

    (fcm[ipop]).read_cell_par(&netpar->C[ipop], fl);

    if (fscanf(fl.tmp, "rhd=%lf Vinc1=%lf Vinc2=%lf\n", 
      &netpar->C[ipop].rhd, &netpar->C[ipop].Vinc1, &netpar->C[ipop].Vinc2)
      != 0)
      fprintf(fl.out, "rhd=%lf Vinc1=%lf Vinc2=%lf\n", 
      netpar->C[ipop].rhd, netpar->C[ipop].Vinc1, netpar->C[ipop].Vinc2);
    else
      printf("rhd\n");

    if (fscanf(fl.tmp, "opto: amp=%lf sig=%lf freq=%lf\n",
      &netpar->C[ipop].opto_amp, &netpar->C[ipop].opto_sig,
      &netpar->C[ipop].opto_freq) != 0) 
      fprintf(fl.out, "opto: amp=%lf sig=%lf freq=%lf\n",
      netpar->C[ipop].opto_amp, netpar->C[ipop].opto_sig,
      netpar->C[ipop].opto_freq);
    else
      printf("opto\n");

    if (fscanf(fl.tmp, "inject_current=%c ion_inject=%d Iinject=%lf tinject=%lf"
      "\n", &netpar->C[ipop].inject_current, &netpar->C[ipop].ion_inject,
      &netpar->C[ipop].Iinject, &netpar->C[ipop].tinject) != 0)
      fprintf(fl.out, "inject_current=%c ion_inject=%d Iinject=%lf tinject=%lf"
      "\n", netpar->C[ipop].inject_current, netpar->C[ipop].ion_inject,
      netpar->C[ipop].Iinject, netpar->C[ipop].tinject);
    else
      printf("inject_current\n");
  }

  if (fgets(line, Mlinea, fl.tmp) != NULL)
    fprintf(fl.out, "%s", line);
  else
    printf("SYNAPSE\n");

  if (fscanf(fl.tmp, "scalingc=%c scaleK=%c Kfactor=%lf CVKin_fact=%lf\n",
    &netpar->scalingc, &netpar->scaleK, &netpar->Kfactor, &netpar->CVKin_fact)
    != 0)
    fprintf(fl.out, "scalingc=%c scaleK=%c Kfactor=%lf CVKin_fact=%lf\n",
    netpar->scalingc, netpar->scaleK, netpar->Kfactor, netpar->CVKin_fact);
  else
    printf("scalingc\n");

  if ((netpar->scalingc != 'v') && (netpar->scalingc != 'k') &&
      (netpar->scalingc != 'G'))
  {
    printf("scalingc=%c should be v or k or G!\n", netpar->scalingc);
    exit(0);
  }

  if ((netpar->scaleK != 's') && (netpar->scaleK != 'w') 
   && (netpar->scaleK != 'n'))
  {
    printf("scaleK=%c should be s or w or n!\n", netpar->scaleK);
    exit(0);
  }

  if (netpar->Kfactor < runpar->epsilon)
  {
    printf("Kfactor=%lf <= 0.0 (epsilon)!\n", netpar->Kfactor);
    exit(0);
  }

  if (fscanf(fl.tmp, "Length=%lf geom_dim=%d con_shape=%c rho_concur=%lf "
    "consider=%c\n", &netpar->Length, &netpar->inter.geom_dim,
    &netpar->inter.con_shape, &netpar->inter.rho_concur, 
    &netpar->inter.consider) != 0)
    fprintf(fl.out, "Length=%lf geom_dim=%d con_shape=%c rho_concur=%lf "
    "consider=%c\n", netpar->Length, netpar->inter.geom_dim,
    netpar->inter.con_shape, netpar->inter.rho_concur, 
    netpar->inter.consider);
  else
    printf("Length\n");

  if ((netpar->inter.consider != 'y') && (netpar->inter.consider != 'n'))
  {
    printf("consider=%c should be y or n!\n", netpar->inter.consider);
    exit(0);
  }

  /* substituing non for cortical populations */
  netpar->C[0].non = netpar->T.ntl;

  /* computing non and rho for cortical populations */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    if (netpar->C[ipop].rhd < runpar->epsilon)
    {
      printf("ipop=%d rhd=%lf should be positive!\n", ipop, 
      netpar->C[ipop].rhd );
      exit(0);
    }

    netpar->C[ipop].nod = (int) ((netpar->Length + runpar->epsilon) *
    netpar->C[ipop].rhd);

    if (netpar->inter.geom_dim <= 1)
    {
      netpar->C[ipop].non = netpar->C[ipop].nod;
      netpar->C[ipop].rho = netpar->C[ipop].rhd;
    }
    else if (netpar->inter.geom_dim == 2)
    {
      netpar->C[ipop].non = netpar->C[ipop].nod * netpar->C[ipop].nod;
      netpar->C[ipop].rho = netpar->C[ipop].rhd * netpar->C[ipop].rhd;
    }
    else
    {
      printf("wrong geom_dim=%d\n", netpar->inter.geom_dim);
      exit(0);
    }

    fprintf(fl.out, "ipop=%d nod=%d non=%d\n", ipop, netpar->C[ipop].nod,
    netpar->C[ipop].non);
  }

  /* Defining jpop_start */
  if (netpar->T.thal_input == 'p')
  {
    /* Thalamocortical synaptic connections - jpop=0 */
    netpar->jpop_start = 0;
  }
  else
  {
    /* Intracortical synaptic connections - jpop=1,...,npop */
    netpar->jpop_start = 1;
  }
  fprintf(fl.out, "jpop_start=%d\n", netpar->jpop_start);

  netpar->Volt_thresh = -20.0;
  fprintf(fl.out, "Volt_thresh=%lf\n", netpar->Volt_thresh);

  if (fscanf(fl.tmp, "AMPA: ths=%lf sigs=%lf tsynd=%lf Vrev=%lf\n",
    &netpar->P.ths, &netpar->P.sigs, &netpar->P.tsynd, &netpar->P.Vrev) != 0)
    {
      netpar->P.tsynr = 0;
      fprintf(fl.out, "AMPA: ths=%lf sigs=%lf tsynr=%lf tsynd=%lf Vrev=%lf\n",
      netpar->P.ths, netpar->P.sigs, netpar->P.tsynr, netpar->P.tsynd,
       netpar->P.Vrev);
    }
    else
    printf("AMPA!\n");

  if (fscanf(fl.tmp, "NMDA: tsynr=%lf tsynd=%lf Vrev=%lf\n",
    &netpar->N.tsynr, &netpar->N.tsynd, &netpar->N.Vrev) != 0)
    fprintf(fl.out, "NMDA: tsynr=%lf tsynd=%lf Vrev=%lf\n",
    netpar->N.tsynr, netpar->N.tsynd, netpar->N.Vrev);
  else
    printf("NMDA1\n");

  if (fscanf(fl.tmp, "ths=%lf sigs=%lf thetanp=%lf sigmanp=%lf\n",
    &netpar->N.ths, &netpar->N.sigs, &netpar->N.thetanp,
    &netpar->N.sigmanp) != 0)
    fprintf(fl.out, "ths=%lf sigs=%lf thetanp=%lf sigmanp=%lf\n",
    netpar->N.ths, netpar->N.sigs, netpar->N.thetanp,
    netpar->N.sigmanp);
  else
    printf("NMDA1\n");

  if (fscanf(fl.tmp, "GABAA: ths=%lf sigs=%lf tsynd=%lf Vrev=%lf\n",
    &netpar->A.ths, &netpar->A.sigs, &netpar->A.tsynd, &netpar->A.Vrev) != 0)
    {
      netpar->A.tsynr = 0.0;
      fprintf(fl.out, "GABAA: ths=%lf sigs=%lf tsynr=%lf tsynd=%lf Vrev=%lf\n",
      netpar->A.ths, netpar->A.sigs, netpar->A.tsynr, netpar->A.tsynd,
      netpar->A.Vrev);
    }
  else
    printf("GABAA\n");

  if (fscanf(fl.tmp, "GABAA_PP: DelVrev_O_DelIext=%lf\n",
    &netpar->A.DelVrev_O_DelIext) != 0)
    fprintf(fl.out, "GABAA_PP: DelVrev_O_DelIext=%lf\n",
    netpar->A.DelVrev_O_DelIext);
  else
    printf("GABAA_PP\n");

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=0; jpop<=netpar->npop; jpop++) /* jpop=0: T */
    {
      netpar->S[ipop][jpop].AMPA.g = 0.0;
      netpar->S[ipop][jpop].NMDA.g = 0.0;
      netpar->S[ipop][jpop].GABAA.g = 0.0;
    }
  }

  if (fscanf(fl.tmp, "factETPT=%lf process_time_delay=%c\n", &netpar->factETPT,
    &netpar-> process_time_delay) != 0)
  {
    fprintf(fl.out, "factETPT=%lf process_time_delay=%c\n", netpar->factETPT,
    netpar-> process_time_delay);
  }
  else
    printf("factETPT\n");

  for (jpop=0; jpop<=netpar->npop; jpop++)
  {
    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      synaptic_parameters_read_write(ipop, jpop, netpar, fl);
      
      if ((netpar->pop_name[ipop] == 'P') && (netpar->pop_name[jpop] == 'P') &&
	  (ipop == jpop))
      {
	if (fscanf(fl.tmp, "gel=%lf Gel=%lf Kel=%lf\n", &netpar->S[ipop][jpop].gel,
        &netpar->S[ipop][jpop].Gel, &netpar->S[ipop][jpop].Kel) !=0)
        {
          fprintf(fl.out, "gel=%lf Gel=%lf Kel=%lf\n", netpar->S[ipop][jpop].gel,
          netpar->S[ipop][jpop].Gel, netpar->S[ipop][jpop].Kel);
        }
        else
          printf("PP gel\n");
      }
    }
  }
 
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=1; jpop<=netpar->npop; jpop++)
    {
      if (netpar->S[ipop][jpop].Kin <= 0.0)
      {
        printf("ipop=%d jpop=%d Kin=%lf<0!\n", ipop, jpop, 
        netpar->S[ipop][jpop].Kin);
        exit(0);
      }
    }
  }

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    if (netpar->T.thal_input == 'a')
    {
      if ((netpar->scalingc =='k') || (netpar->scalingc =='G'))
      {
	if (netpar->scalingc == 'k')
	{
          netpar->C[ipop].gback_nor = netpar->S[ipop][0].AMPA.g * 
          netpar->S[ipop][0].UU * netpar->T.Av * sqrt(netpar->S[ipop][0].Kin);
	}
	else if (netpar->scalingc == 'G')
	{
          netpar->C[ipop].gback_nor = netpar->S[ipop][0].AMPA.G* 
          netpar->S[ipop][0].UU * netpar->T.Av * netpar->S[ipop][0].Kin;
	}       
        fprintf(fl.out, "ipop=%d gAMPA_CT=%lf GAMPA_CT=%lf\n", ipop,
	netpar->S[ipop][0].AMPA.g, netpar->S[ipop][0].AMPA.G);
        fprintf(fl.out, "UU=%lf Av=%lf Kin=%lf gback_nor=%lf\n",
	netpar->S[ipop][0].UU, netpar->T.Av, netpar->S[ipop][0].Kin,
	netpar->C[ipop].gback_nor);
        fprintf(fl.out, "Iext_equiv=%lf\n", netpar->C[ipop].gback_nor * 
        (netpar->P.Vrev - netpar->C[ipop].VL));
      }
      else if (netpar->scalingc =='v')
      {
        ftau = functau(netpar->C[ipop].gL, netpar->C[ipop].Cm, netpar->P.tsynd);

        netpar->C[ipop].g_one_psp = -netpar->S[ipop][0].AMPA.Vpsp *
        netpar->C[ipop].Cm / ((netpar->C[ipop].VL - netpar->P.Vrev) * ftau);

        netpar->C[ipop].gback_nor = netpar->C[ipop].g_one_psp *
        netpar->T.Av * netpar->S[ipop][0].Kin;

        fprintf(fl.out, "ipop=%d Vpsp_AMPA=%lf Av=%lf Kin=%lf DelV=%lf\n",
        ipop, netpar->S[ipop][0].AMPA.Vpsp, netpar->T.Av,
        netpar->S[ipop][0].Kin, netpar->C[ipop].VL - netpar->P.Vrev);
        fprintf(fl.out, "tAMPA=%lf ftau=%lf\n", netpar->P.tsynd, ftau);
        fprintf(fl.out, "g_one_psp=%lf gback_nor=%lf\n",
        netpar->C[ipop].g_one_psp, netpar->C[ipop].gback_nor);
      }
      else
      {
        printf("scalingc=%c should be either v or k\n", netpar->scalingc);
        exit(0); 
      }
    }
    else
    {
      netpar->C[ipop].gback_nor = 0.0;
    }
  }

  if (fscanf(fl.tmp, "%s\n", line) != 0)
    fprintf(fl.out, "GENERAL\n");

  if (fscanf(fl.tmp, "ndeltat=%d", &runpar->ndeltat) !=0)
    fprintf(fl.out, "ndeltat=%d", runpar->ndeltat);
  else
   printf("ndeltat\n");

  runpar->deltat = dvector(1, runpar->ndeltat);
  runpar->nt = ivector(1, runpar->ndeltat);
  runpar->ntall = 0;
  runpar->time_all = 0.0;

  for (ideltat=1; ideltat<=runpar->ndeltat; ideltat++)
  {
    if (fscanf(fl.tmp, " deltat=%lf nt=%d", &runpar->deltat[ideltat], 
      &runpar->nt[ideltat])!= 0)
    {
        fprintf(fl.out, " deltat=%lf nt=%d", runpar->deltat[ideltat],
        runpar->nt[ideltat]);
        runpar->ntall += runpar->nt[ideltat];
        runpar->time_all += runpar->nt[ideltat] * runpar->deltat[ideltat];
    }
    else
      printf("deltat, ideltat=%d\n", ideltat);
  }
  if (fscanf(fl.tmp, "\n") == 0) fprintf(fl.out, "\n");
  fprintf(fl.out, "ntall=%d time_all=%lf\n", runpar->ntall,
  runpar->time_all);

  if (fscanf(fl.tmp, "method=%c incond=%c fpcal=%c smforce=%c\n",
    &runpar->method, &runpar->incond, &runpar->fpcal, &runpar->smforce) != 0)
    fprintf(fl.out, "method=%c incond=%c fpcal=%c smforce=%c\n",
    runpar->method, runpar->incond, runpar->fpcal, runpar->smforce);
  else
    printf("method\n");

  /* Reading data about which cell information to print */
  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    sprintf(line, "%c: nwrite=%%d nwritear=", netpar->pop_name[ipop]);
    if (fscanf(fl.tmp, line, &runpar->nwrite[ipop]) != 0)
    {
      fprintf(fl.out, line, runpar->nwrite[ipop]);

      if (runpar->nwrite[ipop] >= 1)
      {
        runpar->nwritear[ipop] = ivector(1, runpar->nwrite[ipop]);

        if (fscanf(fl.tmp, "%d", &runpar->nwritear[ipop][1]) == 0)
	{
          printf("nwritear1\n");
          exit(0);
	}

        for (iwrite=2; iwrite<=runpar->nwrite[ipop]; iwrite++)
	{
          if (fscanf(fl.tmp, " %d", &runpar->nwritear[ipop][iwrite]) == 0)
	  {
            printf("nwritear%d\n", iwrite);
            exit(0);
	  }
	}

        ifscan = fscanf(fl.tmp, "\n");

        for (iwrite=1; iwrite<=runpar->nwrite[ipop]; iwrite++) 
        {
         if (runpar->nwritear[ipop][iwrite] > netpar->C[ipop].non)
           runpar->nwritear[ipop][iwrite] = netpar->C[ipop].non;
        }

        fprintf(fl.out, "%d", runpar->nwritear[ipop][1]);
        for (iwrite=2; iwrite<=runpar->nwrite[ipop]; iwrite++)
	{
          fprintf(fl.out, " %d", runpar->nwritear[ipop][iwrite]);
	}
        fprintf(fl.out, "\n");
	fflush(fl.out);
      }
    } 
  }

  if (fscanf(fl.tmp, "write_aux=%c twrite=%d tmcol=%lf tstat=%lf traster=%lf "
    "sp=%d\n", &runpar->write_aux, &runpar->twrite, &runpar->tmcol,
    &runpar->tstat, &runpar->traster, &runpar->sp) != 0)
    {
      if (runpar->tstat < runpar->epsilon) runpar->tstat = runpar->epsilon;
      if (runpar->tstat > runpar->time_all) runpar->tstat = runpar->time_all;

      fprintf(fl.out,"write_aux=%c twrite=%d tmcol=%lf tstat=%lf traster=%lf "
      "sp=%d\n", runpar->write_aux, runpar->twrite, runpar->tmcol,
      runpar->tstat, runpar->traster, runpar->sp);
    }
  else
    printf("twrite\n");

  if (fscanf(fl.tmp, "nhist=%d t_touch_interval=%lf\n", &runpar->nhist,
    &runpar->t_touch_interval) != 0)
  {
    fprintf(fl.out, "nhist=%d t_touch_interval=%lf\n", runpar->nhist,
    runpar->t_touch_interval);
  }
  else
    printf("nhist\n");
  
  /* printing flag */
  runpar->sm = parall->sm; /* 0 - print as if there is only one parameter set */
                           /* 1 - no such printing                            */
  if (runpar->smforce == 'p')      /* always print    */
    runpar->sm = 0;
  else if (runpar->smforce == 'n') /* always no print */
    runpar->sm = 1;
  else if (runpar->smforce != 'l') /* leave as is     */
  {
    printf("smforce should be either p or n or l or a !!! smforce=%c\n",
    runpar->smforce);
    exit(0);
  }
  fprintf(fl.out, "sm=%d sp=%d\n", runpar->sm, runpar->sp);

  fflush(fl.out);
}


/* This function reads parameters of the synaptic connections.          */
void synaptic_parameters_read_write(int ipop, int jpop, net_par *netpar,
     fl_st fl)
{
  char ab[3];

  ab[0] = netpar->pop_name[ipop];
  ab[1] = netpar->pop_name[jpop];
  ab[2] = '\n';

  ab_Kin_read_write(&ab[0], &netpar->S[ipop][jpop], fl);
  utxt_read_write(&ab[0], &netpar->S[ipop][jpop], fl);  
  if ((netpar->pop_name[jpop] == 'T') || (netpar->pop_name[jpop] == 'E'))
  {
    AMPA_read_write(&ab[0], &netpar->S[ipop][jpop], &netpar->C[jpop],
    netpar->scalingc, netpar->scaleK, netpar->Kfactor, netpar->CVKin_fact,
    netpar->factETPT, fl);
    NMDA_read_write(&ab[0], &netpar->S[ipop][jpop], &netpar->C[jpop],
    netpar->scalingc, netpar->scaleK, netpar->Kfactor, netpar->CVKin_fact,
    netpar->factETPT, fl);
  }
  else if (netpar->pop_name[jpop] == 'P')
  {
    GABAA_read_write(&ab[0], &netpar->S[ipop][jpop], &netpar->C[jpop],
    netpar->scalingc, netpar->scaleK, netpar->Kfactor, netpar->CVKin_fact, fl);
  }
  else
  {
    printf("wrong jpop=%d\n", jpop);
    exit(0);
  }
}

/* This function reads parameters of the architecture of the synaptic   */
/* connections.                                                         */
void ab_Kin_read_write(char *ab, syn_coup_par *Sij, fl_st fl)
{
  char line[Mlinea];

  if (fgets(line, Mlinea, fl.tmp) != NULL)
    fprintf(fl.out, "%s", line);
  else
    printf("%s\n", ab);

  if (fscanf(fl.tmp, "Kin=%lf CVKin=%lf lam=%lf\n", &Sij->Kin, &Sij->CVKin,
    &Sij->lam) != 0)
  {
    fprintf(fl.out, "Kin=%lf CVKin=%lf lam=%lf\n", Sij->Kin, Sij->CVKin,
    Sij->lam);
  }
  else
  {
    printf("Kin\n");
  }
}

/* This function reads synaptic parameters */
void utxt_read_write(char *ab, syn_coup_par *Sij, fl_st fl)
{
  char string_format[160];

  strcpy(string_format,
  "UU=%lf taur=%lf tauf=%lf xic=%lf tau_delay=%lf Del_tau_delay=%lf\n");
  
  if (fscanf(fl.tmp, string_format, &Sij->UU, &Sij->taur, &Sij->tauf, &Sij->xic,
    &Sij->tau_delay, &Sij->Del_tau_delay) != 0)
  {
    if (Sij->tau_delay < 0.0)
    {
      printf("tau_delay=%lf < 0!\n", Sij->tau_delay);
      exit(0);
    }
    if (Sij->Del_tau_delay < 0.0)
    {
      printf("Del_tau_delay=%lf < 0!\n", Sij->Del_tau_delay);
      exit(0);
    }
    if (Sij->Del_tau_delay > Sij->tau_delay)
    {
      printf("Del_tau_delay=%lf > tau_delay=%lf!\n", Sij->Del_tau_delay,
      Sij->tau_delay);
      exit(0);
    }
    
    strcpy(string_format,
    "UU=%lf taur=%lf tauf=%lf xic=%lf\ntau_delay=%lf Del_tau_delay=%lf\n");
    fprintf(fl.out, string_format, Sij->UU, Sij->taur, Sij->tauf, Sij->xic,
    Sij->tau_delay, Sij->Del_tau_delay);
  }
  else
    printf("%c%c UU\n", ab[0], ab[1]);
}

/* This function reads AMPA synaptic parameters */
void AMPA_read_write(char *ab, syn_coup_par *Sij, cell_par *Cj, char scalingc,
     char scaleK, double Kfactor, double CVKin_fact, double factETPT, fl_st fl)
{
  char string_format[160];
  char line[Mlinea];
  
  strcpy(string_format, "gAMPA=%lf GAMPA=%lf Vpsp_AMPA=%lf\n");

  if (fscanf(fl.tmp, string_format, &Sij->AMPA.g, &Sij->AMPA.G, &Sij->AMPA.Vpsp)
      != 0)
  {
    if ((scalingc == 'k') && (scaleK != 'n'))
    {
      Sij->Kin *= Kfactor;
      if (scaleK == 'w') Sij->AMPA.g /= sqrt(Kfactor);
    }

    if (Sij->Kin > Cj->non)
      Sij->Kin = 1.0 * Cj->non;

    Sij->CVKin *= CVKin_fact;

    Sij->AMPA.g *= factETPT;

    fprintf(fl.out, string_format, Sij->AMPA.g, Sij->AMPA.g, Sij->AMPA.Vpsp);
  }
  else
    printf("%c%c AMPA\n", ab[0], ab[1]);

  fprintf(fl.out, "Kin=%lf CVKin=%lf\n", Sij->Kin, Sij->CVKin);
}

/* This function reads NMDA synaptic parameters */
void NMDA_read_write(char *ab, syn_coup_par *Sij, cell_par *Cj, char scalingc,
     char scaleK, double Kfactor, double CVKin_fact, double factETPT, fl_st fl)
{
  char string_format[160];

  strcpy(string_format, "gNMDA=%lf GNMDA=%lf Vpsp_NMDA=%lf\n");
  
  if (fscanf(fl.tmp, string_format, &Sij->NMDA.g, &Sij->NMDA.G, &Sij->NMDA.Vpsp)
    != 0)
  {
    if ((scalingc == 'k') && (scaleK != 'n'))
    {
      if (scaleK == 'w') Sij->NMDA.g /= sqrt(Kfactor);
    }

    Sij->NMDA.g *= factETPT;

    fprintf(fl.out, string_format, Sij->NMDA.g, Sij->NMDA.G, Sij->NMDA.Vpsp);
  }
  else
    printf("%c%c NMDA\n", ab[0], ab[1]);
}

/* This function reads GABAA synaptic parameters */
void GABAA_read_write(char *ab, syn_coup_par *Sij, cell_par *Cj, char scalingc,
     char scaleK, double Kfactor,  double CVKin_fact, fl_st fl)
{
  char string_format[160];

  strcpy(string_format,
  "gGABAA=%lf GGABAA=%lf Vpsp_GABAA=%lf\n");
  
  if (fscanf(fl.tmp, string_format, &Sij->GABAA.g, &Sij->GABAA.G,
    &Sij->GABAA.Vpsp) != 0)
  {
    if ((scalingc == 'k') && (scaleK != 'n'))
    {
      Sij->Kin *= Kfactor;
      if (scaleK == 'w') Sij->GABAA.g /= sqrt(Kfactor);
    }

    if (Sij->Kin > Cj->non)
      Sij->Kin = 1.0 * Cj->non;

    Sij->CVKin *= CVKin_fact;

    fprintf(fl.out, string_format, Sij->GABAA.g, Sij->GABAA.G, Sij->GABAA.Vpsp);
  }
  else
    printf("%c%c GABAA\n", ab[0], ab[1]);

  fprintf(fl.out, "Kin=%lf CVKin=%lf\n", Sij->Kin, Sij->CVKin);
}

/* Defining Varbar */
void define_variables_arrays(net_par *netpar, run_par *runpar,
     double ****Varbar, fl_st fl)
{
  double ***Varb;
  int ipop, nl, nh;

  /* Varbar */

  nl = 1;
  nh = netpar->npop;

  Varb = (double ***)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int**)));
  if (!Varb) 
    nrerror(strcat("allocation failure in define_variables_arrays()"
                   , " while allocating Varb"));
  *Varbar = Varb-nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    (*Varbar)[ipop] = dmatrix(1, netpar->C[ipop].non, 1, netpar->C[ipop].neq); 
  }
}

/* Freeing Varbar */
void free_variables_arrays(net_par *netpar, run_par *runpar,
     double ***Varbar, fl_st fl)
{
  int ipop, nl, nh;

  /* Varbar */
  nl = 1;
  nh = netpar->npop;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    free_dmatrix(Varbar[ipop], 1, netpar->C[ipop].non, 1, netpar->C[ipop].neq); 
  }
  
  free((FREE_ARG) (Varbar+nl-NR_END));
}

/* Defining Varold */
void define_old_variables_arrays(net_par *netpar, run_par *runpar,
     double ****Varold, int ***after_max_vol, fl_st fl)
{
  double ***Varb;
  int **amv;
  int ipop, nl, nh;

  nl = 1;
  nh = netpar->npop;

  Varb = (double ***)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int**)));
  if (!Varb) 
    nrerror(strcat("allocation failure in define_old_variables_arrays()"
                   , " while allocating Varb"));
  *Varold = Varb-nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    (*Varold)[ipop] = dmatrix(1, netpar->C[ipop].non, 1, 2); 
  }

  amv = (int **)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int*)));
  if (!amv) 
    nrerror(strcat("allocation failure in define_old_variables_arrays()"
                   , " while allocating amv"));
  *after_max_vol = amv-nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    (*after_max_vol)[ipop] = ivector(1, netpar->C[ipop].non); 
  }
}

/* Freeing Varold */
void free_old_variables_arrays(net_par *netpar, run_par *runpar,
     double ***Varold, int **after_max_vol, fl_st fl)
{
  int ipop, nl, nh;

  nl = 1;
  nh = netpar->npop;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    free_dmatrix(Varold[ipop], 1, netpar->C[ipop].non, 1, 2); 
  }
  
  free((FREE_ARG) (Varold+nl-NR_END));

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    free_ivector(after_max_vol[ipop], 1, netpar->C[ipop].non); 
  }
  
  free((FREE_ARG) (after_max_vol+nl-NR_END));
}

/* This function computes neuronal intrinsic properties that are       */
/* heterogeneous.                                                      */
void compute_heterogeneous_intrinsic_parameters(net_par *netpar, 
     run_par *runpar, par_all *parall, fl_st fl)
{
  double rand_num;
  int ipop, ion;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    netpar->C[ipop].gLar = dvector(1, netpar->C[ipop].non);
    netpar->C[ipop].Iextar = dvector(1, netpar->C[ipop].non);
  }

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      rand_num = get_rn_dbl(parall->rng_ptr);
      netpar->C[ipop].gLar[ion] = netpar->C[ipop].gL + netpar->C[ipop].DelgL *
        (2.0 * rand_num - 1.0);

      rand_num = get_rn_dbl(parall->rng_ptr);
      if (ion <=
          (int) (netpar->C[ipop].fracIext * netpar->C[ipop].non + 
          runpar->epsilon))
      {
        netpar->C[ipop].Iextar[ion] = netpar->C[ipop].Iext + 
          netpar->C[ipop].DelIext * (2.0 * rand_num - 1.0);
      }
      else
      {
        netpar->C[ipop].Iextar[ion] = 0.0;
      }
    }

    if (!runpar->sm)
    {
      fprintf(fl.out, "\ngLar  Iextar  ipop=%d\n", ipop);
      for (ion=1; ion<=netpar->C[ipop].non; ion++)
      {
        fprintf(fl.out, "%d %lf %lf\n", ion, netpar->C[ipop].gLar[ion],
        netpar->C[ipop].Iextar[ion]);
      }
    }
  }
}

/* This program frees the arrays of neuronal intrinsic properties. */
void free_heterogeneous_intrinsic_parameters(net_par *netpar, 
     run_par *runpar, fl_st fl)
{
  int ipop;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    free_dvector(netpar->C[ipop].gLar, 1, netpar->C[ipop].non);
  }
}

/* This function computes the strength of optogenetic activation. */
void compute_opto_conductance(net_par *netpar, run_par *runpar, fl_st fl)
{
  double sqtpi, xcen_opt, ycen_opt, xpos, ypos, rdiff_sq, tsigsq, factsig;
  int ipop, ion, ionx, iony;

  /* Defining the arrays */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    netpar->C[ipop].opt_con = dvector(1, netpar->C[ipop].non);
  }

  sqtpi = sqrt(2.0 * Pi);
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    if (netpar->inter.geom_dim <= 1)
    {
      tsigsq = 2.0 * netpar->C[ipop].opto_sig * netpar->C[ipop].opto_sig;
      xcen_opt = netpar->Length / 2.0;

      for (ion=1; ion<=netpar->C[ipop].non; ion++)
      {
        xpos = 1.0 * (ion - 1) / netpar->C[ipop].rhd;
        netpar->C[ipop].opt_con[ion] = netpar->C[ipop].opto_amp *
          exp(-(xpos - xcen_opt) * (xpos - xcen_opt) / tsigsq)
          / (sqtpi * netpar->C[ipop].opto_sig);
      }

      if (!runpar->sm)
      {
        fprintf(fl.out, "ipop=%d opt_con=\n", ipop);
        for (ion=1; ion<=netpar->C[ipop].non; ion++)
	{
          fprintf(fl.out, "ion=%d opt_con=%lf\n", ion,
	  netpar->C[ipop].opt_con[ion]);
	}
      }

    }
    else if (netpar->inter.geom_dim == 2)
    {
      factsig = 2.0 * Pi * netpar->C[ipop].opto_sig * netpar->C[ipop].opto_sig;
      tsigsq = 2.0 * netpar->C[ipop].opto_sig * netpar->C[ipop].opto_sig;
      xcen_opt = netpar->Length / 2.0;
      ycen_opt = netpar->Length / 2.0;
      for (iony=1; iony<=netpar->C[ipop].nod; iony++)
      {
        for (ionx=1; ionx<=netpar->C[ipop].nod; ionx++)
        {
          ion = (iony - 1) * netpar->C[ipop].nod + ionx;
          xpos = 1.0 * (ionx - 1) / netpar->C[ipop].rhd;
          ypos = 1.0 * (iony - 1) / netpar->C[ipop].rhd;
          rdiff_sq = (xpos - xcen_opt) * (xpos - xcen_opt) +  
	             (ypos - ycen_opt) * (ypos - ycen_opt);
          netpar->C[ipop].opt_con[ion] = netpar->C[ipop].opto_amp *
	    exp(-( rdiff_sq /  tsigsq))  / factsig;
	  /*
	  fprintf(fl.out, "ionx=%d iony=%d ion=%d xpos=%lf ypos=%lf rd=%lf "
          "opt_con=%lf\n", ionx, iony, ion, xpos, ypos, rdiff_sq, 
          netpar->C[ipop].opt_con[ion]);
	  */
	}
      }
      if ((!runpar->sm) && (2 > 1))
      {
        fprintf(fl.out, "ipop=%d opt_con=\n", ipop);
        for (iony=1; iony<=netpar->C[ipop].nod; iony++)
	{
          for (ionx=1; ionx<=netpar->C[ipop].nod; ionx++)
	  {
            ion = (iony - 1) * netpar->C[ipop].nod + ionx;
            fprintf(fl.out, "ionx=%d iony=%d ion=%d opt_con=%lf\n", ionx, iony,
	    ion, netpar->C[ipop].opt_con[ion]);
	  }
	}
      }
    }
    else
    {
      printf("wrong geom_dim=%d\n", netpar->inter.geom_dim);
      exit(0);
    }
  }
}

/* This program frees the arrays of optogenetic conductances */
void free_opto_arrays(net_par *netpar, run_par *runpar, fl_st fl)
{
  int ipop;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    free_dvector(netpar->C[ipop].opt_con, 1, netpar->C[ipop].non);
  }
}

/* This function substitutes the connectivity matrices */
void substitute_connectivity(net_par *netpar, run_par *runpar, par_all *parall,
     fl_st fl)
{
  int ipop, jpop;

  /* Chemical coupling */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=netpar->jpop_start; jpop<=netpar->npop; jpop++)
    {
      fprintf(fl.out, "\ncon: ipop=%d jpop=%d Kin=%lf\n", ipop, jpop, 
      netpar->S[ipop][jpop].Kin);
      
      /* Computing the normalization factor */
      if (netpar->inter.geom_dim >= 1)
      {
        compute_normalization_factor(&(netpar->S[ipop][jpop]), &netpar->inter, 
        ipop, jpop, netpar->C[jpop].nod, netpar->C[jpop].non,
	netpar->C[jpop].rhd, netpar->C[jpop].rho, runpar, fl);
      }

      /* Finding the coupling coefficients w_ij */
      if (netpar->inter.geom_dim == 0)
      {
        find_coupling_matrix_zero_d(&(netpar->S[ipop][jpop]), 
        &netpar->C[ipop], &netpar->C[jpop], ipop, jpop, runpar, parall, fl);	
      }
      else if (netpar->inter.geom_dim == 1)
      {
        find_coupling_matrix_one_d(&(netpar->S[ipop][jpop]), 
        netpar->inter.con_shape, netpar->inter.geom_dim, netpar->Length, 
        &netpar->C[ipop], &netpar->C[jpop], ipop, jpop, runpar, parall, fl);
     }
      else if (netpar->inter.geom_dim == 2)
      {
       fprintf(fl.out, "y32 ipop=%d jpop=%d\n", ipop, jpop);
       find_coupling_matrix_two_d(&(netpar->S[ipop][jpop]), 
        netpar->inter.con_shape, netpar->Length, &netpar->C[ipop],
        &netpar->C[jpop], ipop, jpop, runpar, parall, fl);
      }
    }
  }

  /* Electrical coupling */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=netpar->jpop_start; jpop<=netpar->npop; jpop++)
    {
      if ((netpar->pop_name[ipop] == 'P') && (netpar->pop_name[jpop] == 'P'))
      {
        fprintf(fl.out, "\ncon el: ipop=%d jpop=%d Kin=%lf\n", ipop, jpop, 
        netpar->S[ipop][jpop].Kel);
      
        if (netpar->inter.geom_dim <= 1)
        {
          find_electrical_matrix_zero_d(&(netpar->S[ipop][jpop]), 
          netpar->inter.con_shape, netpar->inter.geom_dim, netpar->Length,
          &netpar->C[ipop], &netpar->C[jpop], ipop, jpop, runpar, parall, fl);
        }
        else if (netpar->inter.geom_dim >= 1)
        {
          printf("Electrical coupling for geom_dim=%d not implemented yet\n",
          netpar->inter.geom_dim);
          exit(0);
        }
      }
    }
  }
}

/* This function computes the normalization factor Anor for the probability */
/* of two neurons to be coupled.                                            */
void compute_normalization_factor(syn_coup_par *scp, syn_par_all *inter, 
     int ipop, int jpop, int nod_pre, int non_pre, double rhd_pre, 
     double rho_pre, run_par *runpar, fl_st fl)
{
  double Anor_anal;

  printf("in compute_normalization_factor\n");
  fprintf(fl.out, "ipop=%d jpop=%d Kin=%lf lam=%lf rho_pre=%lf non=%d\n", 
  ipop, jpop, scp->Kin, scp->lam, rho_pre, non_pre);

  if (inter->geom_dim == 0)
  {
    scp->Anor = scp->Kin / non_pre;
    fprintf(fl.out, "non_pre=%d Kin=%lf Anor=%lf\n", non_pre, scp->Kin,
    scp->Anor);
  }
  else if (inter->geom_dim == 1)
  {
    scp->Anor = numerical_normalization_factor_one_d(scp->Kin, scp->lam, 
    rho_pre, non_pre, inter->con_shape, runpar, fl);
  }
  else if (inter->geom_dim == 2)
  {
    scp->Anor = numerical_normalization_factor_two_d(scp->Kin, scp->lam, 
    rhd_pre, nod_pre, inter->con_shape, runpar, fl);
  }
  else
  {
    printf("wrong geom_dim=%d\n", inter->geom_dim);
    exit(0);
  }
  if (!runpar->sm) fprintf(fl.out, "Anor=%lf\n", scp->Anor);
}

/* This function computes the normalization factor Anor numericaly for a 1d  */
/* geometry.                                                               */
double numerical_normalization_factor_one_d(double Kin, double lam, 
       double rho_pre, int non_pre, char con_shape, run_par *runpar, fl_st fl)
{
  double abs_Delx, sum_func, Anor;
  double sqtpi, prob_zero;
  int ion;

  sqtpi = sqrt(2.0 * Pi);

  sum_func = 0.0;

  for (ion=0; ion<=non_pre-1; ion++)
  {
    abs_Delx = min(fabs(1.0 * ion / rho_pre),
		   fabs(1.0 * (non_pre - ion) / rho_pre));

    if (con_shape == 'e')
    {
      sum_func += exp(-abs_Delx / lam) / (2.0 * lam *rho_pre);
    }
    else if (con_shape == 'g')
    {
      sum_func += exp(-abs_Delx * abs_Delx / (2.0 * lam * lam)) /
      (sqtpi * lam * rho_pre);
    }
    else
    {
      printf("wrong con_shape=%c\n", con_shape);
      exit(0);
    }
  }

  Anor = Kin / sum_func;

  /* Cheking that the propability is not above 1 */

  if (con_shape == 'e')
  {
    prob_zero = Anor / (2.0 * lam *rho_pre);
  }
  else if (con_shape == 'g')
  {
    prob_zero =  Anor / (sqtpi * lam * rho_pre);
  }

  fprintf(fl.out, "sum_func=%lf Anor=%lf prob_zero=%lf\n", sum_func, Anor,
  prob_zero);

  if (prob_zero > 1.0)
  {
    printf("prob_zero=%lf > 0!\n", prob_zero);
    exit(0); 
  }
  
  return Anor;
}

/* This function computes the normalization factor Anor numericaly for a 2d  */
/* geometry.                                                               */
double numerical_normalization_factor_two_d(double Kin, double lam, 
       double rhd_pre, int nod_pre, char con_shape, run_par *runpar, fl_st fl)
{
  double abs_Delx, abs_Dely, Delr, sum_func, Anor;
  double factlr, prob_zero;
  int ion, jon;

  factlr = 2.0 * Pi * lam * lam * rhd_pre * rhd_pre;

  sum_func = 0.0;

  for (ion=0; ion<=nod_pre-1; ion++)
  {
    for (jon=0; jon<=nod_pre-1; jon++)
    {
      abs_Delx = min(fabs(1.0 * ion / rhd_pre),
		   fabs(1.0 * (nod_pre - ion) / rhd_pre));

      abs_Dely = min(fabs(1.0 * jon / rhd_pre),
		   fabs(1.0 * (nod_pre - jon) / rhd_pre));

      Delr = sqrt(abs_Delx * abs_Delx + abs_Dely * abs_Dely);

      if (con_shape == 'e')
      {
        sum_func += exp(-abs_Delx / lam) / factlr;
      }
      else if (con_shape == 'g')
      {
        sum_func += exp(-abs_Delx * abs_Delx / (2.0 * lam * lam)) / factlr;
      }
      else
      {
        printf("wrong con_shape=%c\n", con_shape);
        exit(0);
      }
    }
  }

  Anor = Kin / sum_func;

  /* Cheking that the propability is not above 1 */

  if (con_shape == 'e')
  {
    prob_zero = Anor / factlr;
  }
  else if (con_shape == 'g')
  {
    prob_zero =  Anor / factlr;
  }

  fprintf(fl.out, "sum_func=%lf Anor=%lf prob_zero=%lf\n", sum_func, Anor,
  prob_zero);

  if (prob_zero > 1.0)
  {
    printf("prob_zero=%lf > 0!\n", prob_zero);
    exit(0); 
  }
  return Anor;
}

/* This function finds the connectivity matrix for a 0d geometry.        */
void find_coupling_matrix_zero_d(syn_coup_par *scp, cell_par *cpost,
     cell_par *cpre, int ipop, int jpop, run_par *runpar, par_all *parall,
     fl_st fl)
{
  double rand_num, probx;
  double *td_cells;
  int nget_cells, *get_cells, ipre, ipost;
  int n_input_cells, i_input_cells;
  int iwrite;
  double sum;

  /* Finding connectivity matrix */
  scp->nwcoup = ivector(1, cpre->non);
  scp->wcoup = pivector(1, cpre->non);
  scp->tdcoup = pdvector(1, cpre->non);
  scp->nwprob = dvector(1, cpost->non);
  scp->npostw = ivector(1, cpost->non);
  scp->npostp = ivector(1, cpost->non);

  for (ipost=1; ipost<=cpost->non; ipost++) scp->npostw[ipost] = 0;

  compute_coupling_in_prob(scp, cpost, cpre, ipop, jpop, runpar, parall, fl);
  
  /* nget_cells = 10 * ((int) scp->Kin); */
  nget_cells = cpost->non;
  fprintf(fl.out, "ipop=%d jpop=%d nget_cells=%d\n", ipop, jpop, nget_cells);

  get_cells = ivector(1, nget_cells);
  td_cells = dvector(1, nget_cells);

  fprintf(fl.out, "non_pre=%d non_post=%d\n", cpre->non, cpost->non);
  for (ipre=1; ipre <=cpre->non; ipre++)
  {
    n_input_cells = 0;

    /* sum = 0.0; */
    for (ipost=1; ipost<=cpost->non; ipost++)
    {
      /* random_number */
      rand_num = get_rn_dbl(parall->rng_ptr);
 
      probx = scp->nwprob[ipost] / cpre->non;
      /* probx = scp->Kin / cpre->non; */

      if ((probx < 0.0) || (probx > 1.0))
      {
        printf("wrong probx=%lf\n", probx);
        exit(0);
      }

      if (rand_num < probx)
      {
        n_input_cells++;
	scp->npostw[ipost]++;

        if (n_input_cells > nget_cells)
	{
          printf("n_input_cells=%d > nget_cells=%d\n", n_input_cells,
	  nget_cells);
	  exit(0);
	}

        get_cells[n_input_cells] = ipost;

	/* Finding tau_delay for this particular synaptic connection. */
	rand_num = get_rn_dbl(parall->rng_ptr);
	td_cells[n_input_cells] = scp->tau_delay +
	scp->Del_tau_delay * (2.0 * rand_num - 1.0);	
      }
    }

    scp->nwcoup[ipre] = n_input_cells;
    scp->wcoup[ipre] = ivector(1, scp->nwcoup[ipre]);
    scp->tdcoup[ipre] = dvector(1, scp->nwcoup[ipre]);
    /* fprintf(fl.out, "define wcoup: ipre=%d scp->nwcoup[ipre]=%d\n", ipre, 
       scp->nwcoup[ipre]); */
    for (i_input_cells=1; i_input_cells<=n_input_cells; i_input_cells++)
    {
      scp->wcoup[ipre][i_input_cells] = get_cells[i_input_cells];
      scp->tdcoup[ipre][i_input_cells] = td_cells[i_input_cells];
    }

  }

  free_ivector(get_cells, 1, nget_cells);
  free_dvector(td_cells, 1, nget_cells);

  if ((!runpar->sm) && (!runpar->sp))
  {
    for (ipre=1; ipre <=cpre->non; ipre++)
    {
      fprintf(fl.out, "ipre=%d nwcoup=%d\n", ipre, scp->nwcoup[ipre]);
      fflush(fl.out);
    }

    fprintf(fl.out, "ipop=%d jpop=%d \n", ipop, jpop);
    for (ipre=1; ipre <=cpre->non; ipre++)
    {
      fprintf(fl.out, "ipre=%d nwcoup=%d\n", ipre, scp->nwcoup[ipre]);
      for (i_input_cells=1; i_input_cells<=scp->nwcoup[ipre]; i_input_cells++)
      {
        fprintf(fl.out, "%d %lf\n", scp->wcoup[ipre][i_input_cells],
	scp->tdcoup[ipre][i_input_cells]);
      }
      fprintf(fl.out, "\n");
    }
  }

  fprintf(fl.out, "nwcoup stat ipop=%d jpop=%d\n", ipop, jpop);
  compute_connectivity_statistics(scp->nwcoup, cpre->non, runpar, fl);
  fprintf(fl.out, "npostw stat ipop=%d jpop=%d\n", ipop, jpop);
  fprintf(fl.out, "npostw=%d %d %d %d non=%d\n", scp->npostw[1], scp->npostw[2],
	  scp->npostw[3], scp->npostw[4], cpost->non);
  compute_connectivity_statistics(scp->npostw, cpost->non, runpar, fl);
  fflush(fl.out);

  free_ivector(scp->npostw, 1, cpost->non);
  free_ivector(scp->npostp, 1, cpost->non);

}

/* This function computes the "in-degree" value for which the coupling      */
/* probability is computed.                                                 */
void compute_coupling_in_prob(syn_coup_par *scp, cell_par *cpost,
     cell_par *cpre, int ipop, int jpop, run_par *runpar, par_all *parall,
     fl_st fl)
{
  double sigKin, Kin_i, xnoise;
  int ipost;

  sigKin = scp->CVKin * scp->Kin;

  fprintf(fl.out, "ipop=%d jpop=%d Kin=%lf sigKin=%lf\n", ipop, jpop, scp->Kin,
  sigKin);
  
  for (ipost=1; ipost<=cpost->non; ipost++)
  {
    xnoise = gasdev(parall);
    Kin_i = scp->Kin + sigKin * xnoise;

    if (Kin_i < 0.0)
    {
      printf("ipop=%d jpop=%d ipost=%d Kin_i=%lf < 0!\n", ipop, jpop, ipost,
      Kin_i);
      Kin_i = 0.0;
    }
    else if (Kin_i > cpre->non)
    {
      printf("ipop=%d jpop=%d ipost=%d Kin_i=%lf < 0!\n", ipop, jpop, ipost,
      Kin_i);
      Kin_i = 1.0 * cpre->non;
    }

    scp->nwprob[ipost] = Kin_i; /* scp->Kin; */

    if ((!runpar->sm) && (!runpar->sp))
    {
      fprintf(fl.out, "xnoise=%lf nwprob=%lf Kin=%lf sigKin=%lf\n", xnoise, scp->nwprob[ipost], scp->Kin, sigKin); 
    }
  }
  
  compute_connectivity_d_statistics(scp->nwprob, cpost->non, runpar, fl);
}

/* This function finds the connectivity matrix for a 1d geometry.        */
void find_coupling_matrix_one_d(syn_coup_par *scp, char con_shape,
     int geom_dim, double Length, cell_par *cpost, cell_par *cpre, int ipop,
     int jpop, run_par *runpar, par_all *parall, fl_st fl)
{
  double sqtpi, x_pre, x_post, abs_Delx, rand_num, probx;
  double *td_cells;
  int nget_cells, *get_cells, ipre, ipost;
  int n_input_cells, i_input_cells;
  int iwrite;
  double sum;

  sqtpi = sqrt(2.0 * Pi);

  /* Finding connectivity matrix */
  scp->nwcoup = ivector(1, cpre->non);
  scp->wcoup = pivector(1, cpre->non);
  scp->tdcoup = pdvector(1, cpre->non);

  /* nget_cells = 10 * ((int) scp->Kin); */
  nget_cells = cpost->non;
  fprintf(fl.out, "ipop=%d jpop=%d nget_cells=%d\n", ipop, jpop, nget_cells);

  get_cells = ivector(1, nget_cells);
  td_cells = dvector(1, nget_cells);

  fprintf(fl.out, "non_pre=%d non_post=%d\n", cpre->non, cpost->non);
  for (ipre=1; ipre <=cpre->non; ipre++)
  {
    n_input_cells = 0;

    /* sum = 0.0; */
    for (ipost=1; ipost<=cpost->non; ipost++)
    {
      /* abs_Delx */
      x_pre = 1.0 * (ipre - 1) / cpre->rhd;
      x_post = 1.0 * (ipost - 1) / cpost->rhd;
      abs_Delx = min(fabs(x_pre - x_post), 
                     Length - fabs(x_pre - x_post));
      /* fprintf(fl.out, "ipre=%d ipost=%d x_pre=%lf x_post=%lf abs_Delx=%lf\n",
	 ipre, ipost, x_pre, x_post, abs_Delx); */

      /* random_number */
      rand_num = get_rn_dbl(parall->rng_ptr);
 
      if (geom_dim == 1)
      {
        if (con_shape == 'e')
        {
          probx = scp->Anor * exp(-abs_Delx / scp->lam) / 
	    (2.0 * scp->lam * cpre->rhd);
        }
        else if (con_shape == 'g')
        {
          probx = scp->Anor * exp(-abs_Delx * abs_Delx / 
            (2.0 * scp->lam * scp->lam)) / (sqtpi * scp->lam * cpre->rhd);
        }
      }
      else if (geom_dim == 0)
      {
        probx = scp->Anor;
      }
      else
      {
        printf("wrong geom_dim=%d\n", geom_dim);
        exit(0);
      }

      if ((probx < 0.0) || (probx > 1.0))
      {
        printf("wrong probx=%lf\n", probx);
        exit(0);
      }

      /* sum += probx; */

      /* probx = scp->Kin / cpre->non; */
      /* fprintf(fl.out, "rand_num=%lf probx=%lf\n", rand_num, probx); */

      if (rand_num < probx)
      {
        n_input_cells++;

        if (n_input_cells > nget_cells)
	{
          printf("n_input_cells=%d > nget_cells=%d\n", n_input_cells,
	  nget_cells);
	  exit(0);
	}

        get_cells[n_input_cells] = ipost;

	/* Finding tau_delay for this particular synaptic connection. */
	rand_num = get_rn_dbl(parall->rng_ptr);
	td_cells[n_input_cells] = scp->tau_delay +
	scp->Del_tau_delay * (2.0 * rand_num - 1.0);
	
      }
    }
    /* fprintf(fl.out, "ipre=%d sum=%lf\n", ipre, sum); */

    scp->nwcoup[ipre] = n_input_cells;
    scp->wcoup[ipre] = ivector(1, scp->nwcoup[ipre]);
    scp->tdcoup[ipre] = dvector(1, scp->nwcoup[ipre]);
    /* fprintf(fl.out, "define wcoup: ipre=%d scp->nwcoup[ipre]=%d\n", ipre, 
       scp->nwcoup[ipre]); */
    for (i_input_cells=1; i_input_cells<=n_input_cells; i_input_cells++)
    {
      scp->wcoup[ipre][i_input_cells] = get_cells[i_input_cells];
      scp->tdcoup[ipre][i_input_cells] = td_cells[i_input_cells];
    }

  }

  free_ivector(get_cells, 1, nget_cells);
  free_dvector(td_cells, 1, nget_cells);

  if ((!runpar->sm) && (!runpar->sp))
  {
    for (ipre=1; ipre <=cpre->non; ipre++)
    {
      fprintf(fl.out, "ipre=%d nwcoup=%d\n", ipre, scp->nwcoup[ipre]);
      fflush(fl.out);
    }

    fprintf(fl.out, "ipop=%d jpop=%d \n", ipop, jpop);
    for (ipre=1; ipre <=cpre->non; ipre++)
    {
      fprintf(fl.out, "ipre=%d nwcoup=%d\n", ipre, scp->nwcoup[ipre]);
      for (i_input_cells=1; i_input_cells<=scp->nwcoup[ipre]; i_input_cells++)
      {
        fprintf(fl.out, "%d %lf\n", scp->wcoup[ipre][i_input_cells],
	scp->tdcoup[ipre][i_input_cells]);
      }
      fprintf(fl.out, "\n");
    }
  }
  /*
  else
  {
    fprintf(fl.out, "ipop=%d jpop=%d nwrite[jpop]=%d\n", ipop, jpop, 
    runpar->nwrite[jpop]);
    for (iwrite=1; iwrite<=runpar->nwrite[jpop]; iwrite++)
    {
      ipre = runpar->nwritear[jpop][iwrite];
      fprintf(fl.out, "iwrite=%d ipre=%d nwcoup=%d\n", iwrite, ipre,
      scp->nwcoup[ipre]);

      for (i_input_cells=1; i_input_cells<=scp->nwcoup[ipre]; i_input_cells++)
      {
        fprintf(fl.out, "%d ", scp->wcoup[ipre][i_input_cells]);
      }
      fprintf(fl.out, "\n");
    }
  }
  */

  compute_connectivity_statistics(scp->nwcoup, cpre->non, runpar, fl);
  fflush(fl.out);
}

/* This function finds the connectivity matrix for a 1d geometry.        */
void find_coupling_matrix_two_d(syn_coup_par *scp, char con_shape,
     double Length, cell_par *cpost, cell_par *cpre, int ipop, int jpop, 
     run_par *runpar, par_all *parall, fl_st fl)
{
  double x_pre, y_pre, x_post, y_post, abs_Delx, abs_Dely, Delr;
  double factlr, rand_num, probx;
  double sum;
  int nget_cells, *get_cells, iprex, iprey, ipre, ipostx, iposty;
  int n_input_cells, i_input_cells;

  factlr = 2.0 * Pi * scp->lam * scp->lam * cpre->rhd * cpre->rhd;

  /* checking */
  sum = 0.0;
  for (iprex=1; iprex<=cpre->nod; iprex++)
  {
    for (iprey=1; iprey<=cpre->nod; iprey++)
    {
      x_pre = 1.0 * (iprex - 1) / cpre->rhd;
      x_post = 0.0;
      abs_Delx = min(fabs(x_pre - x_post), 
                     Length - fabs(x_pre - x_post));

      y_pre = 1.0 * (iprey - 1) / cpre->rhd;
      y_post = 0.0;
      abs_Dely = min(fabs(y_pre - y_post), 
                     Length - fabs(y_pre - y_post));

      Delr = sqrt(abs_Delx * abs_Delx + abs_Dely * abs_Dely);
      
      if (con_shape == 'e')
      {
        probx = scp->Anor * exp(-abs_Delx / scp->lam) / factlr;
      }
      else if (con_shape == 'g')
      {
        probx = scp->Anor * exp(-abs_Delx * abs_Delx /
                (2.0 * scp->lam * scp->lam)) / factlr;
      }
      sum += probx;
    }
  }
  fprintf(fl.out, "non=%d sum=%lf\n", cpre->non, sum);

  /* Finding connectivity matrix */
  scp->nwcoup = ivector(1, cpre->non);
  scp->wcoup = pivector(1, cpre->non);

  nget_cells = 10 * ((int) scp->Kin);
  fprintf(fl.out, "ipop=%d jpop=%d nget_cells=%d non_pre=%d non_post=%d\n", 
  ipop, jpop, nget_cells, cpre->non, cpost->non);

  get_cells = ivector(1, nget_cells);

  for (iprex=1; iprex <=cpre->nod; iprex++)
  {
    for (iprey=1; iprey <=cpre->nod; iprey++)
    {
      n_input_cells = 0;

      /* sum = 0.0; */
      for (ipostx=1; ipostx<=cpost->nod; ipostx++)
      {
        for (iposty=1; iposty<=cpost->nod; iposty++)
        {
          /* abs_Delx, abs_Dely, Delr */
          x_pre = 1.0 * (iprex - 1) / cpre->rhd;
          x_post = 1.0 * (ipostx - 1) / cpost->rhd;
          abs_Delx = min(fabs(x_pre - x_post), 
                     Length - fabs(x_pre - x_post));

          y_pre = 1.0 * (iprey - 1) / cpre->rhd;
          y_post = 1.0 * (iposty - 1) / cpost->rhd;
          abs_Dely = min(fabs(y_pre - y_post), 
                     Length - fabs(y_pre - y_post));

	  Delr = sqrt(abs_Delx * abs_Delx + abs_Dely * abs_Dely);

          /*
          fprintf(fl.out, "ipre=%d %d ipost=%d %d xy_pre=%lf %lf xy_post="
          "%lf %lf abs_Delxy=%lf %lf Delr=%lf\n", iprex, iprey, ipostx, iposty,
	  x_pre, y_pre, x_post, y_post, abs_Delx, abs_Dely, Delr);
	  */

          /* random_number */
          rand_num = get_rn_dbl(parall->rng_ptr);
 
          if (con_shape == 'e')
          {
            probx = scp->Anor * exp(-Delr / scp->lam) / factlr;
          }
          else if (con_shape == 'g')
          {
            probx = scp->Anor * exp(-Delr * Delr /
                    (2.0 * scp->lam * scp->lam)) / factlr;
          }

          if (rand_num < probx)
          {
            n_input_cells++;

            if (n_input_cells > nget_cells)
   	    {
              printf("n_input_cells=%d > nget_cells=%d\n", n_input_cells,
	      nget_cells);
	      exit(0);
   	    }

            get_cells[n_input_cells] = iposty * cpost->nod + ipostx;
          }
        }
      }

      ipre = (iprey - 1) * cpre->nod + iprex;
      if ((ipre < 1) || (ipre > cpre->non))
      {
        printf("ipre=%d is outside of range [1, %d]\n", ipre, cpre->non); 
      }
      scp->nwcoup[ipre] = n_input_cells;
      scp->wcoup[ipre] = ivector(1, scp->nwcoup[ipre]);
      for (i_input_cells=1; i_input_cells<=n_input_cells; i_input_cells++)
      {
        scp->wcoup[ipre][i_input_cells] = get_cells[i_input_cells];
      }
      /* fprintf(fl.out, "ipre=%d n_input_cells=%d\n", ipre, n_input_cells); */
    }
  }

  free_ivector(get_cells, 1, nget_cells);

  if ((!runpar->sm) && (!runpar->sp))
  {
    for (ipre=1; ipre <=cpre->non; ipre++)
    {
      fprintf(fl.out, "ipre=%d nwcoup=%d\n", ipre, scp->nwcoup[ipre]);
      for (i_input_cells=1; i_input_cells<=scp->nwcoup[ipre]; i_input_cells++)
      {
        fprintf(fl.out, "%d ", scp->wcoup[ipre][i_input_cells]);
      }
      fprintf(fl.out, "\n");
    }
  }

  compute_connectivity_statistics(scp->nwcoup, cpre->non, runpar, fl);
}

/* This function finds the connectivity matrix electrical coupling for    */
/* 0-d geometry.                                                          */
void find_electrical_matrix_zero_d(syn_coup_par *scp, char con_shape,
     int geom_dim, double Length, cell_par *cpost, cell_par *cpre, int ipop,
     int jpop, run_par *runpar, par_all *parall, fl_st fl)
{
  double rand_num, probx;
  int ipre, ipost;

  fprintf(fl.out, "el: ipop=%d jpop=%d\n", ipop, jpop);
  fprintf(fl.out, "el: non_pre=%d non_post=%d\n", cpre->non, cpost->non);

  /* Finding connectivity matrix */
  scp->nw_e_coup = ivector(1, cpost->non);
  for (ipost=1; ipost<=cpre->non; ipost++) scp->nw_e_coup[ipost] = 0;
  scp->w_e_coup = imatrix(1, cpost->non, 1, cpre->non);
    
  for (ipost=1; ipost<=cpost->non; ipost++)
  {
    for (ipre=1; ipre < ipost; ipre++)
    {
     /* random_number */
      rand_num = get_rn_dbl(parall->rng_ptr);

      probx = scp->Kel / cpre->non;

      if ((probx < 0.0) || (probx > 1.0))
      {
        printf("el: wrong probx=%lf\n", probx);
        exit(0);
      }

      if (rand_num < probx)
      {
        scp->nw_e_coup[ipost]++;
        scp->w_e_coup[ipost][scp->nw_e_coup[ipost]] = ipre;

        scp->nw_e_coup[ipre]++;
        scp->w_e_coup[ipre][scp->nw_e_coup[ipre]] = ipost;
      }
    }
  }

  if ((!runpar->sm) && (!runpar->sp))
  {
    for (ipost=1; ipost <=cpost->non; ipost++)
    {
      fprintf(fl.out, "el: ipost=%d nw_e_coup=%d\n", ipost, 
      scp->nw_e_coup[ipost]);
      fflush(fl.out);
    }

    fprintf(fl.out, "el: ipop=%d jpop=%d \n", ipop, jpop);
    for (ipost=1; ipost <=cpre->non; ipost++)
    {
      fprintf(fl.out, "xy ipost=%d nw_e_coup=%d\n", ipost,
      scp->nw_e_coup[ipost]);
      for (ipre=1; ipre<=scp->nw_e_coup[ipost]; ipre++)
      {
        fprintf(fl.out, "%d ", scp->w_e_coup[ipost][ipre]);
      }
      fprintf(fl.out, "\n");
      fflush(fl.out);
    }
  }

  compute_connectivity_statistics(scp->nw_e_coup, cpre->non, runpar, fl);
}

/* Free connectivity arrays */
void free_connectivity_arrays(net_par *netpar, run_par *runpar, fl_st fl)
{
  int ipop, jpop, ipre, ipost;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=1; jpop<=netpar->npop; jpop++)
    {
      for (ipre=1; ipre<=netpar->C[jpop].non; ipre++)
      {
        free_ivector(netpar->S[ipop][jpop].wcoup[ipre], 1, 
	netpar->S[ipop][jpop].nwcoup[ipre]);
	free_dvector(netpar->S[ipop][jpop].tdcoup[ipre], 1, 
	netpar->S[ipop][jpop].nwcoup[ipre]);
      }

      free_pivector(netpar->S[ipop][jpop].wcoup, 1, netpar->C[jpop].non);
      free_ivector(netpar->S[ipop][jpop].nwcoup, 1, netpar->C[jpop].non);
      free_pdvector(netpar->S[ipop][jpop].tdcoup, 1, netpar->C[jpop].non);
      free_dvector(netpar->S[ipop][jpop].nwprob, 1, netpar->C[ipop].non);
    }
  }

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=1; jpop<=netpar->npop; jpop++)
    {
      if ((netpar->pop_name[ipop] == 'P') && (netpar->pop_name[jpop] == 'P'))
      {
        free_imatrix(netpar->S[ipop][jpop].w_e_coup, 1, netpar->C[ipop].non,
        1, netpar->C[jpop].non);
        free_ivector(netpar->S[ipop][jpop].nw_e_coup, 1, netpar->C[ipop].non);
      }
    }
  }
}

/* This function computes the average and the standard deviation of the  */
/* number of neurons projecting from one pre-synaptic neuron.            */
void compute_connectivity_statistics(int *ncoup, int non_pre,
     run_par *runpar, fl_st fl)
{
  double sum_num, sum_num_2, diff_sum, av_num, sig_num;
  int ipre;

  if (!runpar->sm)
  {
    fprintf(fl.out, "nwcoup\n");
    for (ipre=1; ipre<=non_pre; ipre++)
    {
      fprintf(fl.out, "%d %d\n", ipre, ncoup[ipre]);
    }
  }
  
  sum_num = 0.0;
  sum_num_2 = 0.0;

  for (ipre=1; ipre<=non_pre; ipre++)
  {
    sum_num += ncoup[ipre];
    sum_num_2 += ncoup[ipre] * ncoup[ipre];
  }

  sum_num /= non_pre;
  sum_num_2 /= non_pre;
  
  av_num = sum_num;
  diff_sum = sum_num_2 - sum_num * sum_num;

  if (diff_sum >= 0.0)
  {
    sig_num = sqrt(diff_sum);
  }
  else if (diff_sum >= -runpar->epsilon)
  {
    sig_num = 0.0;
  }
  else
  {
    fprintf(fl.out, "diff_sum=%le sig_num=%le\n", diff_sum, sig_num);
    sig_num = -999.9;
  }

  if (!runpar->sm)
    fprintf(fl.out, "av_num=%lf sig_num=%lf\n", av_num, sig_num);
}

/* This function computes the average and the standard deviation of the  */
/* number of neurons projecting to one post-synaptic neuron.             */
void compute_connectivity_d_statistics(double *ncoup, int non_post,
     run_par *runpar, fl_st fl)
{
  double sum_num, sum_num_2, diff_sum, av_num, sig_num;
  int ipost;

  if (!runpar->sm)
  {
    fprintf(fl.out, "npost\n");
    for (ipost=1; ipost<=non_post; ipost++)
    {
      fprintf(fl.out, "%d %lf\n", ipost, ncoup[ipost]);
    }
  }
  
  sum_num = 0.0;
  sum_num_2 = 0.0;

  for (ipost=1; ipost<=non_post; ipost++)
  {
    sum_num += ncoup[ipost];
    sum_num_2 += ncoup[ipost] * ncoup[ipost];
  }

  sum_num /= non_post;
  sum_num_2 /= non_post;
  
  av_num = sum_num;
  diff_sum = sum_num_2 - sum_num * sum_num;

  if (diff_sum >= 0.0)
  {
    sig_num = sqrt(diff_sum);
  }
  else
  {
    sig_num = -999.9;
  }

  if (!runpar->sm)
    fprintf(fl.out, "av_num=%lf sig_num=%lf\n", av_num, sig_num);
}

/* This function substitutes the initial conditions.     */
void in_con(func_cell_model *fcm, double ***Varbar, net_par *netpar,
     run_par *runpar, par_all *parall, fl_st fl)
{
  int ipop;
  int ifscan;

  if (runpar->incond == 'r')
  {
    ifscan = fscanf(fl.tmp, "\n");
    ifscan = fscanf(fl.tmp, "INITIAL CONDITIONS\n");
  }

  if (!runpar->sm)
  {
    fprintf(fl.out, "\nInitial conditions\n");
  }

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    in_con_one_pop(&fcm[ipop], Varbar[ipop], &netpar->C[ipop], &netpar->inter,
    &netpar->P, ipop, runpar, parall, fl);
  }
}

/* This function substitutes the initial conditions for one population of  */
/* neurons.                                                                */
void in_con_one_pop(func_cell_model *fcm, double **Varbar, cell_par *cpar, 
     syn_par_all *inter, syn_receptor_par *AMPA, int ipop, run_par *runpar,
     par_all *parall, fl_st fl)
{
  double rand_num;
  int ion, ieq, ifscan;
  char line1[Mlinea], line2[Mlinea];

  fprintf(fl.out, "model=%c%c\n", cpar->model_type[0], cpar->model_type[1]); 
  if (runpar->incond == 'r')
  {
    if (fgets(line1, Mlinea, fl.tmp) == NULL) printf("line1=NULL\n");
    if (fgets(line2, Mlinea, fl.tmp) == NULL) printf("line1=NULL\n");
    /* E: V     h     n     b     z */
    /* I: V     h     n     a     b */

    for (ion=1; ion<=cpar->non; ion++)  
    {
      for (ieq=1; ieq<=cpar->nceq; ieq++)
      {
        if (fscanf(fl.tmp, "%lf", &(Varbar[ion][ieq])) == 0)
	{
          fprintf(fl.out, "Reading Varbar failed!\n");
          exit(0);
	}
      }
      ifscan = fscanf(fl.tmp, "\n");
    }
  }
  else if (runpar->incond == 'b')
  {
    for (ion=1; ion<=cpar->non; ion++)
    {
      rand_num = get_rn_dbl(parall->rng_ptr);
      Varbar[ion][1] = cpar->Vinc1 + (cpar->Vinc2 - cpar->Vinc1) * rand_num;
      fcm->steady_state_var(cpar, Varbar[ion], runpar, fl);
    }
  }
  else if (runpar->incond == 'e')
  {
    for (ion=1; ion<=cpar->non; ion++)
    {
      if (cpar->non >= 2)
      {
        Varbar[ion][1] = cpar->Vinc1 + (cpar->Vinc2 - cpar->Vinc1) *
        ((1.0 * ion) - 1.0) / ((1.0 * cpar->non) - 1.0);
      }
      else
      {
	Varbar[ion][1] = cpar->Vinc1;
      }
      fcm->steady_state_var(cpar, Varbar[ion], runpar, fl);
    }
  }
  else if (runpar->incond == 'c')
  {
    for (ion=1; ion<=cpar->non; ion++)
    {
      rand_num = get_rn_dbl(parall->rng_ptr);
      Varbar[ion][1] = cpar->Vinc1 + (cpar->Vinc2 - cpar->Vinc1) * rand_num;
      Varbar[ion][2] = 0.8;
      Varbar[ion][3] = 0.2;
      Varbar[ion][4] = 0.1;
    }    
  }
  else
  {
    printf("Wrong incond:%c\n", runpar->incond);
    exit(1);
  }

  if (runpar->fpcal == 'y')
  {
    for (ion=1; ion<=cpar->non; ion++)
    {
      find_fixed_point(fcm, cpar, runpar, inter, AMPA, Varbar[ion], ipop, ion,
      fl);
    }
  }
  else if (runpar->fpcal != 'n')
  {
    printf("fpcal should be y or n !!!\n");
    exit(0);
  }

  for (ion=1; ion<=cpar->non; ion++)  
  {
    for (ieq=cpar->nceq+1; ieq<=cpar->neq; ieq++)
    {
      Varbar[ion][ieq] = 0.0;
    }
  }

  if (!runpar->sm)
  {
    fprintf(fl.out, "Printing initial conditions\n");
    if (runpar->incond == 'r')
    {
      fprintf(fl.out, "%s", line1);
      fprintf(fl.out, "%s", line2);
    }

    for (ion=1; ion<=cpar->non; ion++)
    {
      fprintf(fl.out, "%5d", ion);
      fprintf(fl.out, " %10.5lf", Varbar[ion][1]);
      for (ieq=2; ieq<=cpar->neq; ieq++) 
        fprintf(fl.out, " %8.5lf", Varbar[ion][ieq]);
      fprintf(fl.out, "\n");
    }

  }
  fflush(fl.out);
}

/* Defining the arrays in synstr */
void define_synaptic_variables(syn_str *synstr, net_par *netpar, 
     run_par *runpar, fl_st fl)
{
  syn_par **synpar_var, *synpar_var_vec;
  double **Isyn_cont_var, *Isyn_cont_var_vec;
  double **Iel_cont_var, *Iel_cont_var_vec;
  double **el_coef_sum_var, *el_coef_sum_var_vec;
  double **el_coef_V_var, *el_coef_V_var_vec;
  double ***v3, **v2, *v1;
  char ***c3, **c2, *c1;
  int ipop, jpop, ion, jon, nl, nh;
  int isynvar;

  if (!runpar->sm) fprintf(fl.out, "\ndefine_synaptic_variables\n");

  synstr->nsynvar = ivector(1, netpar->npop);
  synstr->nmda_on_population = ivector(1, netpar->npop);

  /* defining csynvar */
  synstr->mc = 100;
  nl = 1;
  nh = netpar->npop;
  c3 = (char ***)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(char**)));
  if (!c3) nrerror(
    "allocation failure in define_synaptic_variables() while defining c3");
  synstr->csynvar = c3 -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    nh = synstr->mc;
    c2 = (char **)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(char*)));
    if (!c2) nrerror(
      "allocation failure in define_synaptic_variables() while defining c2");
    synstr->csynvar[ipop] = c2 -nl+NR_END;

    for (isynvar=1; isynvar<=synstr->mc; isynvar++)
    {
      nl = 0;
      nh = 1;
      c1 = (char *)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(char)));
      if (!c1) nrerror(
       "allocation failure in define_synaptic_variables() while defining c1");
      synstr->csynvar[ipop][isynvar] = c1 -nl+NR_END;
    }
  }

  /* counting to find nsynvar */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    if (!runpar->sm) fprintf(fl.out, "def ipop=%d\n", ipop);

    synstr->nsynvar[ipop] = 0;
    synstr->nmda_on_population[ipop] = 0;

    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      /* AMPA */
      substitute_condition_check_syn(ipop, jpop, 'P', 1,
      &netpar->S[ipop][jpop].condition_for_nonzero_synapse_AMPA,
      netpar->S[ipop][jpop].AMPA, netpar->scalingc, synstr, runpar, fl);

      /* NMDA */
      substitute_condition_check_syn(ipop, jpop, 'N', 2,
      &netpar->S[ipop][jpop].condition_for_nonzero_synapse_NMDA,
      netpar->S[ipop][jpop].NMDA, netpar->scalingc, synstr, runpar, fl);

      /* GABA_A */
      substitute_condition_check_syn(ipop, jpop, 'A', 1,
      &netpar->S[ipop][jpop].condition_for_nonzero_synapse_GABAA,
      netpar->S[ipop][jpop].GABAA, netpar->scalingc, synstr, runpar, fl);
    }

    if (!runpar->sm)
    {
      fprintf(fl.out, "ipop=%d nsynvar=%d nmda_on_population=%d\n", ipop, 
      synstr->nsynvar[ipop], synstr->nmda_on_population[ipop]);
    }
  }

  /* defining synapr */
  nl = 1;
  nh = netpar->npop;
  synpar_var = (syn_par **)malloc((size_t) ((nh-nl+1+NR_END) * 
               sizeof(syn_par*)));
  if (!synpar_var) nrerror(strcat(
     "allocation failure in define_synaptic_variables()",
     " while defining synpar_var"));
  synstr->synpar = synpar_var -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    nh = synstr->nsynvar[ipop];

    synpar_var_vec = (syn_par *)malloc((size_t) ((nh-nl+1+NR_END) * 
		     sizeof(syn_par)));
    if (!synpar_var_vec) nrerror(strcat(
      "allocation failure in define_synaptic_variables()",
      " while defining synpar_var_vec"));
    synstr->synpar[ipop] = synpar_var_vec -nl+NR_END;

    for (isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
    {
      synstr->synpar[ipop][isynvar].Vsyn = dvector(1, netpar->C[ipop].non);
    }
  }

  /* defining synvar */
  nl = 1;
  nh = netpar->npop;
  v3 = (double ***)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double**)));
  if (!v3) nrerror(
    "allocation failure in define_synaptic_variables() while defining v3");
  synstr->synvar = v3 -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    nh = netpar->C[ipop].non;
    v2 = (double **)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double*)));
    if (!v2) nrerror(
      "allocation failure in define_synaptic_variables() while defining v2");
    synstr->synvar[ipop] = v2 -nl+NR_END;

    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      nl = 1;
      nh = synstr->nsynvar[ipop];
      v1 = (double *)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double)));
      if (!v1) nrerror(
       "allocation failure in define_synaptic_variables() while defining v1");
      synstr->synvar[ipop][ion] = v1 -nl+NR_END;
    }
  }

  /* defining synvar_avt */
  nl = 1;
  nh = netpar->npop;
  v3 = (double ***)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double**)));
  if (!v3) nrerror(
    "allocation failure in define_synaptic_variables() while defining v3");
  synstr->synvar_avt = v3 -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    nh = netpar->C[ipop].non;
    v2 = (double **)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double*)));
    if (!v2) nrerror(
      "allocation failure in define_synaptic_variables() while defining v2");
    synstr->synvar_avt[ipop] = v2 -nl+NR_END;

    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      nl = 1;
      nh = synstr->nsynvar[ipop];
      v1 = (double *)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double)));
      if (!v1) nrerror(
       "allocation failure in define_synaptic_variables() while defining v1");
      synstr->synvar_avt[ipop][ion] = v1 -nl+NR_END;
    }
  }

  /* defining Isynvar_avt */
  nl = 1;
  nh = netpar->npop;
  v3 = (double ***)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double**)));
  if (!v3) nrerror(
    "allocation failure in define_synaptic_variables() while defining v3");
  synstr->Isynvar_avt = v3 -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    nh = netpar->C[ipop].non;
    v2 = (double **)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double*)));
    if (!v2) nrerror(
      "allocation failure in define_synaptic_variables() while defining v2");
    synstr->Isynvar_avt[ipop] = v2 -nl+NR_END;

    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      nl = 1;
      nh = synstr->nsynvar[ipop];
      v1 = (double *)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double)));
      if (!v1) nrerror(
       "allocation failure in define_synaptic_variables() while defining v1");
      synstr->Isynvar_avt[ipop][ion] = v1 -nl+NR_END;
    }
  }

  /* defining isyn_to_send */
  /*                              |---- post ---|  |---- pre ----|      */
  synstr->isyn_to_send = i3tensor(1, netpar->npop, 0, netpar->npop, 1, 2);

  /* defining nonzero_gsyn */
  /*                             |---- post ---|  |---- pre ----|*/
  synstr->nonzero_gsyn = imatrix(1, netpar->npop, 0, netpar->npop);

  /* defining Isyn_cont */
  nl = 1;
  nh = netpar->npop;
  Isyn_cont_var = (double **)malloc((size_t) ((nh-nl+1+NR_END) * 
               sizeof(double*)));
  if (!Isyn_cont_var) nrerror(strcat(
     "allocation failure in define_synaptic_variables()",
     " while defining Isyn_cont_var"));
  synstr->Isyn_cont = Isyn_cont_var -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    nh = netpar->C[ipop].non;

    Isyn_cont_var_vec = (double *)malloc((size_t) ((nh-nl+1+NR_END) * 
		     sizeof(double)));
    if (!Isyn_cont_var_vec) nrerror(strcat(
      "allocation failure in define_synaptic_variables()",
      " while defining synpar_var_vec"));
    synstr->Isyn_cont[ipop] = Isyn_cont_var_vec -nl+NR_END;
  }

  /* defining xvar */
  nl = 1;
  nh = netpar->npop;
  v3 = (double ***)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double**)));
  if (!v3) nrerror(
    "allocation failure in define_synaptic_variables() while defining v3");
  synstr->xvar = v3 -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 0;
    nh = netpar->npop;
    v2 = (double **)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double*)));
    if (!v2) nrerror(
      "allocation failure in define_synaptic_variables() while defining v2");
    synstr->xvar[ipop] = v2 -nl+NR_END;

    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      nl = 1;
      nh = netpar->C[jpop].non;
      v1 = (double *)malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double)));
      if (!v1) nrerror(
       "allocation failure in define_synaptic_variables() while defining v1");
      synstr->xvar[ipop][jpop] = v1 -nl+NR_END;

      for (jon=1; jon<=netpar->C[jpop].non; jon++)
      {
        synstr->xvar[ipop][jpop][jon] = netpar->S[ipop][jpop].xic;
      }
    }
  }

  if (!runpar->sm)
  {
    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      for (jpop=0; jpop<=netpar->npop; jpop++)
      {
	fprintf(fl.out, "ipop=%d jpop=%d xvar\n", ipop, jpop);
	for (jon=1; jon<=netpar->C[jpop].non; jon++)
	{
	  fprintf(fl.out, "jon=%d xvar=%lf\n", jon,
          synstr->xvar[ipop][jpop][jon]);
	}
      }
    }
  }

  /* defining Iel_cont */
  nl = 1;
  nh = netpar->npop;
  Iel_cont_var = (double **)malloc((size_t) ((nh-nl+1+NR_END) * 
               sizeof(double*)));
  if (!Iel_cont_var) nrerror(strcat(
     "allocation failure in define_synaptic_variables()",
     " while defining Iel_cont_var"));
  synstr->Iel_cont = Iel_cont_var -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    nh = netpar->C[ipop].non;

    Iel_cont_var_vec = (double *)malloc((size_t) ((nh-nl+1+NR_END) * 
		     sizeof(double)));
    if (!Iel_cont_var_vec) nrerror(strcat(
      "allocation failure in define_synaptic_variables()",
      " while defining synpar_var_vec"));
    synstr->Iel_cont[ipop] = Iel_cont_var_vec -nl+NR_END;
  }

  /* defining el_coef_sum */
  nl = 1;
  nh = netpar->npop;
  el_coef_sum_var = (double **)malloc((size_t) ((nh-nl+1+NR_END) * 
               sizeof(double*)));
  if (!el_coef_sum_var) nrerror(strcat(
     "allocation failure in define_synaptic_variables()",
     " while defining el_coef_sum_var"));
  synstr->el_coef_sum = el_coef_sum_var -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    nh = netpar->C[ipop].non;

    el_coef_sum_var_vec = (double *)malloc((size_t) ((nh-nl+1+NR_END) * 
		     sizeof(double)));
    if (!el_coef_sum_var_vec) nrerror(strcat(
      "allocation failure in define_synaptic_variables()",
      " while defining synpar_var_vec"));
    synstr->el_coef_sum[ipop] = el_coef_sum_var_vec -nl+NR_END;
  }

  /* defining el_coef_V */
  nl = 1;
  nh = netpar->npop;
  el_coef_V_var = (double **)malloc((size_t) ((nh-nl+1+NR_END) * 
               sizeof(double*)));
  if (!el_coef_V_var) nrerror(strcat(
     "allocation failure in define_synaptic_variables()",
     " while defining el_coef_V_var"));
  synstr->el_coef_V = el_coef_V_var -nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    nh = netpar->C[ipop].non;

    el_coef_V_var_vec = (double *)malloc((size_t) ((nh-nl+1+NR_END) * 
		     sizeof(double)));
    if (!el_coef_V_var_vec) nrerror(strcat(
      "allocation failure in define_synaptic_variables()",
      " while defining synpar_var_vec"));
    synstr->el_coef_V[ipop] = el_coef_V_var_vec -nl+NR_END;
  }

  fflush(fl.out);
}

/* This function suvstitutes synaptic parameters. */
void substitute_condition_check_syn(int ipop, int jpop, char cell_ch,
     int num_add, int *condition_for_nonzero_synapse, conductance cnd,
     char scalingc, syn_str *synstr, run_par *runpar, fl_st fl)
{
  if (scalingc == 'k')
  {
    /* scaling according to sqrt(k) */
    *condition_for_nonzero_synapse = fabs(cnd.g) >= runpar->epsilon;
  }
  else if (scalingc == 'G')
  {
    /* scaling according to G only */
    *condition_for_nonzero_synapse = fabs(cnd.G) >= runpar->epsilon;
  }
  else if (scalingc == 'v')
  {
    /* scaling according to Vpsp */
    *condition_for_nonzero_synapse = fabs(cnd.Vpsp) >= runpar->epsilon;
  }
  else
  {
    printf("scalingc=%c should be either k or v\n", scalingc);
    exit(0);
  }

  if (*condition_for_nonzero_synapse)
  {
    synstr->nsynvar[ipop] += num_add;
    update_check_nsynvar_mc(synstr, ipop, jpop, cell_ch, fl);

    if (cell_ch == 'N')
      synstr->nmda_on_population[ipop] = 1;
  }

  if (!runpar->sm)
    fprintf(fl.out, "%c ipop=%d jpop=%d cond=%d nsynvar=%d\n", cell_ch, ipop,
     jpop, *condition_for_nonzero_synapse, synstr->nsynvar[ipop]);
}

/* This function updates the labels in csynvar and checks whether isynvar  */
/* is indeed smaller than mc.                                              */
void update_check_nsynvar_mc(syn_str *synstr,  int ipop, int jpop, char label,
     fl_st fl)
{
  char str_pr[2];

  if (synstr->nsynvar[ipop] > synstr->mc)
  {
    printf("ipop=%d nsynvar=%d > mc=%d!\n", ipop, synstr->nsynvar[ipop],
    synstr->mc);
    exit(0); 
  }

  sprintf(str_pr, "%d", jpop);

  if ((label == 'P') || label == 'A')
  {
    synstr->csynvar[ipop][synstr->nsynvar[ipop]][0] = str_pr[0];
    synstr->csynvar[ipop][synstr->nsynvar[ipop]][1] = label;
  }
  else if (label == 'N')
  {
    synstr->csynvar[ipop][synstr->nsynvar[ipop]-1][0] = str_pr[0];
    synstr->csynvar[ipop][synstr->nsynvar[ipop]-1][1] = 'N';
    synstr->csynvar[ipop][synstr->nsynvar[ipop]][0] = str_pr[0];
    synstr->csynvar[ipop][synstr->nsynvar[ipop]][1] = 'M';
  }
  else
  {
    printf("wrong label=%c\n", label);
    exit(0);
  }
}

/* Freeing the arrays in synstr */
void free_synaptic_variables(syn_str *synstr, net_par *netpar, 
     run_par *runpar, fl_st fl)
{
  int nl, ipop, jpop, ion, isynvar;

  free_ivector(synstr->nsynvar, 1, netpar->npop);
  free_ivector(synstr->nmda_on_population, 1, netpar->npop);

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
    {
      free_dvector(synstr->synpar[ipop][isynvar].Vsyn, 1, netpar->C[ipop].non);
    }

    nl = 1;
    free((FREE_ARG) (synstr->synpar[ipop] +nl-NR_END));
  }

  nl  = 1;
  free((FREE_ARG) (synstr->synpar +nl-NR_END));

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      nl = 1;
      free((FREE_ARG) (synstr->synvar[ipop][ion] +nl-NR_END));
    }

    nl = 1;
    free((FREE_ARG) (synstr->synvar[ipop] +nl-NR_END));
  }
  nl = 1;
  free((FREE_ARG) (synstr->synvar +nl-NR_END));
  
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      nl = 1;
      free((FREE_ARG) (synstr->synvar_avt[ipop][ion] +nl-NR_END));
    }

    nl = 1;
    free((FREE_ARG) (synstr->synvar_avt[ipop] +nl-NR_END));
  }
  nl = 1;
  free((FREE_ARG) (synstr->synvar_avt +nl-NR_END));

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      nl = 1;
      free((FREE_ARG) (synstr->Isynvar_avt[ipop][ion] +nl-NR_END));
    }

    nl = 1;
    free((FREE_ARG) (synstr->Isynvar_avt[ipop] +nl-NR_END));
  }
  nl = 1;
  free((FREE_ARG) (synstr->Isynvar_avt +nl-NR_END));

  free_i3tensor(synstr->isyn_to_send, 1, netpar->npop, 0, netpar->npop, 1, 2);

  free_imatrix(synstr->nonzero_gsyn, 1, netpar->npop, 0, netpar->npop);

  /* Isyn_cont */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    free((FREE_ARG) (synstr->Isyn_cont[ipop] +nl-NR_END));
  }

  nl  = 1;
  free((FREE_ARG) (synstr->Isyn_cont +nl-NR_END));

  /* xvar */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      nl = 1;
      free((FREE_ARG) (synstr->xvar[ipop][jpop] +nl-NR_END));
    }

    nl = 0;
    free((FREE_ARG) (synstr->xvar[ipop] +nl-NR_END));
  }

  nl = 1;
  free((FREE_ARG) (synstr->xvar +nl-NR_END));

  /* Iel_cont */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    free((FREE_ARG) (synstr->Iel_cont[ipop] +nl-NR_END));
  }

  nl  = 1;
  free((FREE_ARG) (synstr->Iel_cont +nl-NR_END));

  /* el_coef_sum */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    free((FREE_ARG) (synstr->el_coef_sum[ipop] +nl-NR_END));
  }

  nl  = 1;
  free((FREE_ARG) (synstr->el_coef_sum +nl-NR_END));

  /* el_coef_V */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 1;
    free((FREE_ARG) (synstr->el_coef_V[ipop] +nl-NR_END));
  }

  nl  = 1;
  free((FREE_ARG) (synstr->el_coef_V +nl-NR_END));
}

/* This function initiates the structure of synaptic variables */
void initiate_synaptic_strengths_and_variables(syn_str *synstr, net_par *netpar,
     run_par *runpar, fl_st fl)
{
  double ftau, Vextr, g_one_psp, G_norm_bal, fNMDA;
  int ipop, jpop, ion, isynvar, first_input_to_ipop;

  /* counting to find the parameter values for syn_par and initiating  */
  /* isyn_to_send                                                      */

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    isynvar = 0;

    for (jpop=netpar->jpop_start; jpop<=netpar->npop; jpop++)
    {
      first_input_to_ipop = 0;
      synstr->isyn_to_send[ipop][jpop][1] = 0;

      /* AMPA */
      if (netpar->S[ipop][jpop].condition_for_nonzero_synapse_AMPA)
      {
	one_synaptic_type_strengths_variables(ipop, jpop, 'P', &netpar->P,
	netpar->S[ipop][jpop].AMPA, netpar->S[ipop][jpop].Kin,
	netpar->S[ipop][jpop].UU, &netpar->C[ipop], netpar->scalingc,
	&first_input_to_ipop, &isynvar, synstr, runpar, fl);
      }

      /* NMDA */
      if (netpar->S[ipop][jpop].condition_for_nonzero_synapse_NMDA)
      {
	one_synaptic_type_strengths_variables(ipop, jpop, 'N', &netpar->N,
	netpar->S[ipop][jpop].NMDA, netpar->S[ipop][jpop].Kin,
	netpar->S[ipop][jpop].UU, &netpar->C[ipop], netpar->scalingc,
	&first_input_to_ipop, &isynvar, synstr, runpar, fl);	
      }
      
      /* GABA_A */
      if (netpar->S[ipop][jpop].condition_for_nonzero_synapse_GABAA)
      {
	one_synaptic_type_strengths_variables(ipop, jpop, 'A', &netpar->A,
	netpar->S[ipop][jpop].GABAA, netpar->S[ipop][jpop].Kin,
	netpar->S[ipop][jpop].UU, &netpar->C[ipop], netpar->scalingc,
	&first_input_to_ipop, &isynvar, synstr, runpar, fl);	
      }
    
      if (synstr->isyn_to_send[ipop][jpop][1] == 0)
      {
        synstr->nonzero_gsyn[ipop][jpop] = 0;
        synstr->isyn_to_send[ipop][jpop][2] = 0;
      }
      else
      {
        synstr->nonzero_gsyn[ipop][jpop] = 1;
        synstr->isyn_to_send[ipop][jpop][2] = isynvar;
      }

    }
  }

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      synstr->Isyn_cont[ipop][ion] = 0.0;
    }
  }

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      synstr->Iel_cont[ipop][ion] = 0.0;
    }
  }

  if (!runpar->sm)
  {
    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      for (isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
      {
        fprintf(fl.out, "ipop=%d isynvar=%d gsyn=%lf Vsyn=%lf tsyn=%lf nmda=%d"
	  "\n", ipop, isynvar, synstr->synpar[ipop][isynvar].gsyn, 
          synstr->synpar[ipop][isynvar].Vsyn[1], 
          synstr->synpar[ipop][isynvar].tsyn, 
          synstr->synpar[ipop][isynvar].is_it_nmda);
	  /* Halorhodopsin */
  	  for (ion=1; ion<=netpar->C[ipop].non; ion++)
	  {
	    fprintf(fl.out, "ipop=%d isynvar=%d Vsyn=%lf\n", ipop, isynvar,
	    synstr->synpar[ipop][isynvar].Vsyn[ion]);
	  }
      }
    }
    fflush(fl.out);

    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      for (jpop=0; jpop<=netpar->npop; jpop++)
      {
        fprintf(fl.out, "ipop=%d jpop=%d nonzero_gsyn=%d isyn_to_send=[%d, %d]"
	"\n", ipop, jpop, synstr->nonzero_gsyn[ipop][jpop],
	synstr->isyn_to_send[ipop][jpop][1], synstr->isyn_to_send[ipop][jpop][2]
	);
      }
    }

    fprintf(fl.out, "\n");
  }

  /* initiating synvar */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      for (isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
      {
        synstr->synvar[ipop][ion][isynvar] = 0.0;
        synstr->synvar_avt[ipop][ion][isynvar] = 0.0;
        synstr->Isynvar_avt[ipop][ion][isynvar] = 0.0;
      }
    }
  }
}

/* This function initiates the structure of synaptic variables for one   */
/* type of synapses (AMPA, NMDA, GABAA).                                 */
void one_synaptic_type_strengths_variables(int ipop, int jpop, char cell_ch,
     syn_receptor_par *R, conductance cnd, double Kin, double UU, cell_par *CC,
     char scalingc, int *first_input_to_ipop, int *isynvar, syn_str *synstr,
     run_par *runpar, fl_st fl)
{
  double G_norm_bal, g_one_psp, ftau, Vextr, fNMDA;
  int ion;
  
  (*isynvar)++;
  for (ion=1; ion<=CC->non; ion++)
  {
    synstr->synpar[ipop][*isynvar].Vsyn[ion] = R->Vrev;

    /* Halorhodopsin */
    if ((cell_ch == 'A') && (ipop == 2) && (jpop == 2))
    {
      synstr->synpar[ipop][*isynvar].Vsyn[ion] += 
        CC->Iextar[ion] * R->DelVrev_O_DelIext;
    }
  }
  synstr->synpar[ipop][*isynvar].tsyn = R->tsynd;
  synstr->synpar[ipop][*isynvar].is_it_nmda = (cell_ch != 'N') ? 0 : 1;

  if (scalingc == 'k')
  {
    synstr->synpar[ipop][*isynvar].gsyn = cnd.g / (sqrt(Kin) *
      (R->tsynd - R->tsynr));
    G_norm_bal = cnd.g;
    g_one_psp  = G_norm_bal / sqrt(Kin);
  }
  else if (scalingc == 'G')
  {
    synstr->synpar[ipop][*isynvar].gsyn = cnd.G / (R->tsynd - R->tsynr);
    G_norm_bal = cnd.G * sqrt(Kin);
    g_one_psp  = cnd.G;
  }
 else if (scalingc == 'v')
  {
    if ((cell_ch == 'P') || (cell_ch == 'A'))
    {
      ftau = functau(CC->gL, CC->Cm, R->tsynd);
      fNMDA = 1.0;
    }
    else if (cell_ch == 'N')
    {
      ftau = functau_NMDA(CC->gL, CC->Cm, R->tsynr, R->tsynd, runpar, fl);
      fNMDA = Gammaf(CC->VL, R->thetanp, R->sigmanp);
      if (!runpar->sm)
	 fprintf(fl.out, "ipop=%d jpop=%d fNMDA=%lf\n", ipop, jpop, fNMDA);
    }
    else
    {
      printf("cell_ch=%c should be P, N or A!\n", cell_ch);
      exit(0);
    }

    g_one_psp = -cnd.Vpsp * CC->Cm / ((CC->VL - R->Vrev) * ftau) / fNMDA;
    synstr->synpar[ipop][*isynvar].gsyn = g_one_psp / (R->tsynd - R->tsynr);
    G_norm_bal = g_one_psp * sqrt(Kin);
  }

  (*first_input_to_ipop)++;
  if (*first_input_to_ipop == 1) 
    synstr->isyn_to_send[ipop][jpop][1] = *isynvar;
  
  if  ((cell_ch == 'P') || (cell_ch == 'A'))
  {
    ftau = functau(CC->gL, CC->Cm, synstr->synpar[ipop][*isynvar].tsyn);
  }
  else if (cell_ch == 'N')
  {
    (*isynvar)++;
    for (ion=1; ion<=CC->non; ion++)
    {
      synstr->synpar[ipop][*isynvar].Vsyn[ion] = R->Vrev;
    }
    
    synstr->synpar[ipop][*isynvar].tsyn = R->tsynr;
    synstr->synpar[ipop][*isynvar].is_it_nmda = 1;
    synstr->synpar[ipop][*isynvar].gsyn =
      -synstr->synpar[ipop][*isynvar-1].gsyn;
    ftau = functau_NMDA(CC->gL, CC->Cm, R->tsynr, R->tsynd, runpar, fl);
  }
    
  Vextr = -(synstr->synpar[ipop][*isynvar].gsyn * UU *
    (CC->VL - R->Vrev) / CC->Cm) * ftau * (R->tsynd - R->tsynr);
   
  fprintf(fl.out, "ipop=%d jpop=%d %c gsyn=%lf UU=%lf Vpsp=%lf\n",
    ipop, jpop, cell_ch, cnd.g, UU, cnd.Vpsp);
  fprintf(fl.out, "g_one_psp=%lf G_norm_bal=%lf\n", g_one_psp, G_norm_bal);
  fprintf(fl.out, "Kin=%lf tsynr=%lf tsynd=%lf DelV=%lf\n", Kin, R->tsynr,
   R->tsynd, CC->VL - R->Vrev);
  fprintf(fl.out, "gsyn_cal=%lf", synstr->synpar[ipop][*isynvar].gsyn);
  if (cell_ch == 'N')
    fprintf(fl.out, " %lf", synstr->synpar[ipop][*isynvar-1].gsyn);
  fprintf(fl.out, " ftau=%lf Vextr=%lf\n", ftau, Vextr); 
    }

/* This function computes the strength of single electrical synapse  */
void initiate_electrical_strengths(syn_str *synstr, net_par *netpar,
     run_par *runpar, fl_st fl)
{
  int ipop, jpop;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=netpar->jpop_start; jpop<=netpar->npop; jpop++)
    {
      if ((netpar->pop_name[ipop] == 'P') && (netpar->pop_name[jpop] == 'P'))
      {
        if (netpar->scalingc == 'k')
        {
          synstr->gel_one_syn = netpar->S[ipop][jpop].gel / 
            sqrt(netpar->S[ipop][jpop].Kel);
        }
        else if (netpar->scalingc == 'G')
        {
          synstr->gel_one_syn = netpar->S[ipop][jpop].Gel;
        }

        fprintf(fl.out, "ipop=%d jpop=%d gel=%lf Gl=%lf gel_one_syn=%lf\n",
	ipop, jpop, netpar->S[ipop][jpop].gel, netpar->S[ipop][jpop].Gel,
        synstr->gel_one_syn);
      }
    }
  }
}

/* This function defines the arrays in the spike structure spk_str */
void define_spike_struct(spk_str *spkstr, net_par *netpar, run_par *runpar, 
     fl_st fl)
{
  int ipop, ihist;

  spkstr->non_all = 0;
  for (ipop=0; ipop<=netpar->npop; ipop++)
    spkstr->non_all += netpar->C[ipop].non;
  fprintf(fl.out, "non_all=%d\n", spkstr->non_all);

  spkstr->tspike = dvector(1, spkstr->non_all);
  spkstr->jpop = ivector(1, spkstr->non_all);
  spkstr->jon = ivector(1, spkstr->non_all);

  spkstr->mspk = 10;
  define_pop_non_nsome(netpar, runpar, &spkstr->t_p_delay_spike, spkstr->mspk, 
  fl);

  define_pop_pop_non_spk(netpar, runpar, &spkstr->x_u_delay_spike, 
  spkstr->mspk, fl);
  define_int_pop_non(netpar, runpar, &spkstr->spike_exist, fl);
  spkstr->nspike_pop = ivector(0, netpar->npop);
  for (ipop=0; ipop<=netpar->npop; ipop++) spkstr->nspike_pop[ipop] = 0;

  if (netpar->process_time_delay == 's') /* spread of time delays */
    define_storage_spikes_td(netpar, runpar, &spkstr->time_spike_delay,
    &spkstr->ndelay, &spkstr->iptr_delay, &spkstr->sp_in_dt,
    spkstr->non_all, fl);
  
  define_pop_non(netpar, runpar, &spkstr->tfire_prev, fl);
  set_val_pop_non(-1.0, netpar, runpar, spkstr->tfire_prev, fl);
  define_pop_non(netpar, runpar, &spkstr->av_deltat_fire, fl);
  define_pop_non(netpar, runpar, &spkstr->fr, fl);
  define_pop_non(netpar, runpar, &spkstr->frinter, fl);
  define_pop_non(netpar, runpar, &spkstr->av_deltat_fire_sq, fl);
  define_pop_non(netpar, runpar, &spkstr->sig_deltat_fire, fl);
  define_pop_non(netpar, runpar, &spkstr->cv_deltat_fire, fl);
  define_int_pop_non(netpar, runpar, &spkstr->nfire, fl);
  define_pop_non(netpar, runpar, &spkstr->Z1cos, fl);
  define_pop_non(netpar, runpar, &spkstr->Z1sin, fl);
  define_pop_non(netpar, runpar, &spkstr->Z1md, fl);
  define_pop_non(netpar, runpar, &spkstr->Z1phi, fl);

  spkstr->av_deltat_fire_pop     = dvector(1, netpar->npop);
  spkstr->sig_av_deltat_fire_pop = dvector(1, netpar->npop);
  spkstr->fr_pop                 = dvector(1, netpar->npop);
  spkstr->fr_subpop              = dmatrix(1, netpar->npop, 1, 2);
  spkstr->fr_pop_sd               = dvector(1, netpar->npop);
  spkstr->frinter_pop            = dvector(1, netpar->npop);
  spkstr->cv_deltat_fire_pop     = dvector(1, netpar->npop);
  spkstr->sig_cv_deltat_fire_pop = dvector(1, netpar->npop);

  spkstr->fr_hist                = imatrix(0, netpar->npop+1, 0,
    runpar->nhist-1);
  define_pop_non(netpar, runpar, &spkstr->spk_touch, fl);
  define_pop_non(netpar, runpar, &spkstr->spk_before_touch, fl);

  spkstr->non_with_cv = ivector(0, netpar->npop);
  spkstr->non_no_firing = ivector(0, netpar->npop);
      
  for (ipop=0; ipop<=netpar->npop+1; ipop++)
    for (ihist=0; ihist<=runpar->nhist-1; ihist++)
      spkstr->fr_hist[ipop][ihist] = 0;

  define_Vav(netpar, runpar, &spkstr->Vav, fl);
}

/* This function frees the arrays in the spike structure spk_str */
void free_spike_struct(spk_str *spkstr, net_par *netpar, run_par *runpar, 
     fl_st fl)
{
  free_dvector(spkstr->tspike, 1, spkstr->non_all);
  free_ivector(spkstr->jpop, 1, spkstr->non_all);
  free_ivector(spkstr->jon, 1, spkstr->non_all);

  free_pop_non_nsome(netpar, runpar, spkstr->t_p_delay_spike, spkstr->mspk, fl);
  free_pop_pop_non_spk(netpar, runpar, spkstr->x_u_delay_spike, 
  spkstr->mspk, fl);
  free_int_pop_non(netpar, runpar, spkstr->spike_exist, fl);

  free_ivector(spkstr->nspike_pop, 0, netpar->npop);

  if (netpar->process_time_delay == 's')  /* spread of time delays */
    free_storage_spikes_td(netpar, runpar, spkstr->time_spike_delay,
    spkstr->ndelay, spkstr->iptr_delay, spkstr->sp_in_dt, spkstr->non_all, fl);

  free_pop_non(netpar, runpar, spkstr->tfire_prev, fl);
  free_pop_non(netpar, runpar, spkstr->av_deltat_fire, fl);
  free_pop_non(netpar, runpar, spkstr->av_deltat_fire_sq, fl);
  free_pop_non(netpar, runpar, spkstr->fr, fl);
  free_pop_non(netpar, runpar, spkstr->frinter, fl);
  free_pop_non(netpar, runpar, spkstr->sig_deltat_fire, fl);
  free_pop_non(netpar, runpar, spkstr->cv_deltat_fire, fl);
  free_int_pop_non(netpar, runpar, spkstr->nfire, fl);
  free_pop_non(netpar, runpar, spkstr->Z1cos, fl);
  free_pop_non(netpar, runpar, spkstr->Z1sin, fl);
  free_pop_non(netpar, runpar, spkstr->Z1md, fl);
  free_pop_non(netpar, runpar, spkstr->Z1phi, fl);
  free_ivector(spkstr->non_with_cv, 0, netpar->npop);
  free_ivector(spkstr->non_no_firing, 0, netpar->npop);


  free_dvector(spkstr->av_deltat_fire_pop    , 1, netpar->npop);
  free_dvector(spkstr->sig_av_deltat_fire_pop, 1, netpar->npop);
  free_dvector(spkstr->fr_pop                , 1, netpar->npop);
  free_dmatrix(spkstr->fr_subpop             , 1, netpar->npop, 1, 2);
  free_dvector(spkstr->fr_pop_sd             , 1, netpar->npop);
  free_dvector(spkstr->frinter_pop           , 1, netpar->npop);
  free_dvector(spkstr->cv_deltat_fire_pop    , 1, netpar->npop);
  free_dvector(spkstr->sig_cv_deltat_fire_pop, 1, netpar->npop);

  free_imatrix(spkstr->fr_hist,                0, netpar->npop+1, 0, 
    runpar->nhist-1);
  free_pop_non(netpar, runpar, spkstr->spk_touch, fl);
  free_pop_non(netpar, runpar, spkstr->spk_before_touch, fl);

  free_Vav(netpar, runpar, spkstr->Vav, fl);
}

/* This function defines the arrays time_spike_delay, ndelay, iptr_delay,  */
/* sp_in_dt .                                                      */
void define_storage_spikes_td(net_par *netpar, run_par *runpar,
     double *****time_spike_delay, int ***ndelay, int ***iptr_delay,
     spike_in_deltat **sp_in_dt, int non_all, fl_st fl)
{
  double deltat, max_delay_time;
  double ****tt4, ***tt3;
  int ipop, jpop, ion, jon, ndelta_in_delay, idelay, nl, nh;

  deltat = runpar->deltat[runpar->ndeltat];
  fprintf(fl.out, "\ndelay  deltat=%lf\n", deltat);

  /* Defining ndelay, the number of time steps in the delay stack. */
  *ndelay = imatrix(1, netpar->npop, 0, netpar->npop);
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      max_delay_time = netpar->S[ipop][jpop].tau_delay +
      netpar->S[ipop][jpop].Del_tau_delay;
      ndelta_in_delay = (int) ((max_delay_time / deltat) + runpar->epsilon) + 2;
      (*ndelay)[ipop][jpop] = ndelta_in_delay;
      fprintf(fl.out, "ipop=%d jpop=%d max_delay_time=%lf ndelay=%d\n", ipop,
      jpop, max_delay_time, (*ndelay)[ipop][jpop]);
    }
  }

  /* Defining **iptr_delay, the pointer indicating the present time step */
  /* in the stack.                                                       */
  *iptr_delay = imatrix(1, netpar->npop, 0, netpar->npop);
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      (*iptr_delay)[ipop][jpop] = (*ndelay)[ipop][jpop];
      fprintf(fl.out, "ipop=%d jpop=%d iptr_ndelay=%d\n", ipop, jpop,
      (*iptr_delay)[ipop][jpop]);
    }
  }

  /* Defining time_spike_delay [ipop][jpop][ion][ndelay] */
  nl = 1;
  nh = netpar->npop;
  tt4 = (double ****)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double***)));
  if (!tt4) 
    nrerror(strcat("allocation failure in define_pop_pop_non_spk()"
                   , " while allocating tt4"));
  *time_spike_delay = tt4-nl+NR_END;
  
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 0;
    nh = netpar->npop;
    tt3 = (double ***)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double**)));
    if (!tt3) 
    nrerror(strcat("allocation failure in define_pop_pop_non_spk()"
                   , " while allocating tt3"));
    (*time_spike_delay)[ipop] = tt3-nl+NR_END;
    
    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      
      (*time_spike_delay)[ipop][jpop] = pdvector(1, netpar->C[ipop].non);
      for (ion=1; ion<=netpar->C[ipop].non; ion++)
      {
        (*time_spike_delay)[ipop][jpop][ion] = dvector(1,
	(*ndelay)[ipop][jpop]);

        for (idelay=1; idelay<=(*ndelay)[ipop][jpop]; idelay++)
        {
          (*time_spike_delay)[ipop][jpop][ion][idelay] = 0.0;
        }
      }
    }
  }

  /* Defining sp_in_dt [0=<ispk<=sum_non][2: jpop, jon] */
  nl = 1;
  nh = non_all;
  *sp_in_dt = (spike_in_deltat *) malloc((size_t)((nh-nl+1+NR_END)*
    sizeof(spike_in_deltat)));
  if (!(*sp_in_dt))
    nrerror("allocation failure 1 in vector _spike_in_deltat()");
  *sp_in_dt += -nl+NR_END;
}

/* This function frees the arrays time_spike_delay, ndelay, iptr_delay, */
void free_storage_spikes_td(net_par *netpar, run_par *runpar,
     double ****time_spike_delay, int **ndelay, int **iptr_delay,
     spike_in_deltat *sp_in_dt, int non_all, fl_st fl)
{
  int ipop, jpop, ion, nl, nh;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
   for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      for (ion=1; ion<=netpar->C[ipop].non; ion++)
      {
        free_dvector(time_spike_delay[ipop][jpop][ion], 1, ndelay[ipop][jpop]);
      }

      free_pdvector(time_spike_delay[ipop][jpop], 1, netpar->C[ipop].non);
    }
  
    nl = 0;
    nh = netpar->npop;
    free((FREE_ARG) (time_spike_delay[ipop]+nl-NR_END));
  }

  nl = 1;
  nh = netpar->npop;
  free((FREE_ARG) (time_spike_delay+nl-NR_END));

  free_imatrix(ndelay, 1, netpar->npop, 0, netpar->npop);
  free_imatrix(iptr_delay, 1, netpar->npop, 0, netpar->npop);

  nl = 1;
  nh = non_all;
  free((FREE_ARG) (sp_in_dt+nl-NR_END));
}

/* Defining tsp[ipop][ion] */
void define_pop_non(net_par *netpar, run_par *runpar, double ***tsp, fl_st fl)
{
  double **ttt;
  int ipop, nl, nh, ion, non;

  nl = 0;
  nh = netpar->npop+1;

  ttt = (double **)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double*)));
  if (!ttt) 
    nrerror(strcat("allocation failure in define_pop_non()"
                   , " while allocating Varb"));
  *tsp = ttt-nl+NR_END;

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    if (ipop <= netpar->npop)
      non = netpar->C[ipop].non;
    else
      non = netpar->C[ipop-1].non;

    (*tsp)[ipop] = dvector(1, non);

    for (ion=1; ion<=non; ion++)
    {
      (*tsp)[ipop][ion] = 0.0;
    }
  }
}

/* This function sets the value val to all the elements of jagar. */
void set_val_pop_non(double val, net_par *netpar, run_par *runpar,
     double **jagar, fl_st fl)
{
  int ipop, ion;

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      jagar[ipop][ion] = val;
    }
  }
}

/* Freeing tsp[ipop][ion] */
void free_pop_non(net_par *netpar, run_par *runpar, double **tsp, fl_st fl)
{
  int ipop, nl, nh, non;

  nl = 0;
  nh = netpar->npop;

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    if (ipop <= netpar->npop)
      non = netpar->C[ipop].non;
    else
      non = netpar->C[ipop-1].non;

    free_dvector(tsp[ipop], 1, non);
  }
  
  free((FREE_ARG) (tsp+nl-NR_END));
}

/* Defining t_p_delay_spike[ipop][ion][kspk] */
void define_pop_non_nsome(net_par *netpar, run_par *runpar, double ****tsp, 
     int nsome, fl_st fl)
{
  double ***ttt;
  int ipop, nl, nh, ion, ispk;

  nl = 0;
  nh = netpar->npop;

  ttt = (double ***)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double**)));
  if (!ttt) 
    nrerror(strcat("allocation failure in define_pop_non()"
                   , " while allocating Varb"));
  *tsp = ttt-nl+NR_END;

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    (*tsp)[ipop] = pdvector(1, netpar->C[ipop].non);

    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      (*tsp)[ipop][ion] = dvector(1, nsome);

      for (ispk=1; ispk<=nsome; ispk++)
      {
        (*tsp)[ipop][ion][ispk] = 0.0;
      }
    }
  }
}

/* Freeing t_p_delay_spike[ipop][ion][kspk] */
void free_pop_non_nsome(net_par *netpar, run_par *runpar, double ***tsp,
     int nsome, fl_st fl)
{
  int ipop, nl, nh, ion;

  nl = 0;
  nh = netpar->npop;

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      free_dvector(tsp[ipop][ion], 1, nsome);
    }

    free_pdvector(tsp[ipop], 1, netpar->C[ipop].non);
  }
  
  free((FREE_ARG) (tsp+nl-NR_END));
}

/* Defining x_u_delay_spike[jpop][ipop][ion][kspk] */
void define_pop_pop_non_spk(net_par *netpar, run_par *runpar, double *****tsp, 
     int mspk, fl_st fl)
{
  double ****tt4, ***tt3;
  int ipop, jpop, nl, nh, jon, kspk;

  nl = 1;
  nh = netpar->npop;

  tt4 = (double ****)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double***)));
  if (!tt4) 
    nrerror(strcat("allocation failure in define_pop_pop_non_spk()"
                   , " while allocating tt4"));
  *tsp = tt4-nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    nl = 0;
    nh = netpar->npop;
    tt3 = (double ***)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double**)));
    if (!tt3) 
    nrerror(strcat("allocation failure in define_pop_pop_non_spk()"
                   , " while allocating tt3"));
    (*tsp)[ipop] = tt3-nl+NR_END;

    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      (*tsp)[ipop][jpop] = pdvector(1, netpar->C[jpop].non);

      for (jon=1; jon<=netpar->C[jpop].non; jon++)
      {
        (*tsp)[ipop][jpop][jon] = dvector(1, mspk);

        for (kspk=1; kspk<=mspk; kspk++)
        {
          (*tsp)[ipop][jpop][jon][kspk] = 0.0;
        }
      }
    }
  }
}

/* Freeing x_u_delay_spike[ipop][ion][kspk] */
void free_pop_pop_non_spk(net_par *netpar, run_par *runpar, double ****tsp, 
     int mspk, fl_st fl)
{
  int ipop, jpop, jon, nl, nh;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      for (jon=1; jon<=netpar->C[jpop].non; jon++)
      {
        free_dvector(tsp[ipop][jpop][jon], 1, mspk);
      }

      free_pdvector(tsp[ipop][jpop], 1, netpar->C[jpop].non);
    }
  
    nl = 0;
    nh = netpar->npop;
    free((FREE_ARG) (tsp[ipop]+nl-NR_END));
  }

  nl = 1;
  nh = netpar->npop;
  free((FREE_ARG) (tsp+nl-NR_END));
}

/* Defining nsp[ipop][ion] */
void define_int_pop_non(net_par *netpar, run_par *runpar, int ***nsp,
     fl_st fl)
{
  int **ttt;
  int ipop, nl, nh, ion;

  nl = 0;
  nh = netpar->npop;

  ttt = (int **)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int*)));
  if (!ttt) 
    nrerror(strcat("allocation failure in define_pop_non()"
                   , " while allocating Varb"));
  *nsp = ttt-nl+NR_END;

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    (*nsp)[ipop] = ivector(1, netpar->C[ipop].non);

    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      (*nsp)[ipop][ion] = 0;
    }
  }
}

/* Freeing nsp[ipop][ion] */
void free_int_pop_non(net_par *netpar, run_par *runpar, int **nsp, fl_st fl)
{
  int ipop, nl, nh;

  nl = 0;
  nh = netpar->npop;

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    free_ivector(nsp[ipop], 1, netpar->C[ipop].non);
  }
  
  free((FREE_ARG) (nsp+nl-NR_END));
}

/* This function defines the array Vav  */
void define_Vav(net_par *netpar, run_par *runpar, V_aver **Vav, fl_st fl)
{
  V_aver *Vav_def;
  int ipop, ion, nl, nh;

  /* Vav */

  nl = 1;
  nh = netpar->npop;

  Vav_def = (V_aver *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(V_aver)));
  if (!Vav_def) 
    nrerror(strcat("allocation failure in define_variables_arrays()"
                   , " while allocating Varb"));
  *Vav = Vav_def-nl+NR_END;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    (*Vav)[ipop].V_avt = dvector(1, netpar->C[ipop].non);
    (*Vav)[ipop].V_sq_avt = dvector(1, netpar->C[ipop].non);
  }

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      (*Vav)[ipop].V_avt[ion] = 0.0;
      (*Vav)[ipop].V_sq_avt[ion] = 0.0;
    }
    (*Vav)[ipop].Vpop_avt = 0.0;
    (*Vav)[ipop].Vpop_sq_avt = 0.0;
  }

  
 
}

/* This function freess the array Vav  */
void free_Vav(net_par *netpar, run_par *runpar, V_aver *Vav, fl_st fl)
{
  int ipop, nl, nh;

  /* Vav */

  nl = 1;
  nh = netpar->npop;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    free_dvector(Vav[ipop].V_avt, 1, netpar->C[ipop].non);
    free_dvector(Vav[ipop].V_sq_avt, 1, netpar->C[ipop].non);
  }
  
  free((FREE_ARG) (Vav+nl-NR_END));
}

/* This function determines the values of thalamic touch responses Cv during */
/* protraction and retraction.                                               */
void determine_thal_Cv(net_par *netpar, run_par *runpar, par_all *parall,
     fl_st fl)
{
  int itl_both_Cv, itl_p_Cv, itl_r_Cv, itl;
  
  if ((netpar->T.wd != 'p') && (netpar->T.wd != 'r'))
  {
    printf("wrong wd=%c\n", netpar->T.wd);
    exit(0);
  }

  if ((netpar->T.ds_type != 'a') && (netpar->T.ds_type != 'u') &&
      (netpar->T.ds_type != 'n'))
  {
    printf("wrong ds_type=%c\n", netpar->T.ds_type);
    exit(0);
  }

  itl_both_Cv = (int) ((1.0 - netpar->T.frac_only_p - netpar->T.frac_only_r) *
  netpar->T.ntl + runpar->epsilon);
  itl_p_Cv = (int) ((netpar->T.frac_only_p) * netpar->T.ntl + runpar->epsilon);
  itl_r_Cv = (int) ((netpar->T.frac_only_r) * netpar->T.ntl + runpar->epsilon);
  fprintf(fl.out, "itl_both_Cv=%d itl_p_Cv=%d itl_r_Cv=%d\n", itl_both_Cv,
	  itl_p_Cv, itl_r_Cv); 
  
  if (netpar->T.ds_type == 'n')    /* Cvmin */
  {
    for (itl=1; itl<=netpar->T.ntl; itl++)
    {
      netpar->T.Cv[itl] = netpar->T.Cvmin;
    }
  }
  if (netpar->T.ds_type == 'a')    /* average */
  {
    for (itl=1; itl<=netpar->T.ntl; itl++)
    {
      if (itl <= itl_both_Cv)
      {
	netpar->T.Cv[itl] = 0.5 * (netpar->T.Cvmin + netpar->T.Cvmax);
      }
      else if (itl <= itl_both_Cv + itl_p_Cv)
      {
	netpar->T.Cv[itl] = (netpar->T.wd == 'p') ?
	  0.5 * (netpar->T.Cvmin + netpar->T.Cvmax) :0.0;
      }
      else
      {
	netpar->T.Cv[itl] = (netpar->T.wd == 'r') ?
	  0.5 * (netpar->T.Cvmin + netpar->T.Cvmax) :0.0;
      }
    }
  }
  else if (netpar->T.ds_type == 'u')    /* uniform */
  {
    for (itl=1; itl<=netpar->T.ntl; itl++)
    {
      if (itl <= itl_both_Cv)
      {
	if (itl_both_Cv > 1)
          netpar->T.Cv[itl] = netpar->T.Cvmin +
	  (netpar->T.Cvmax - netpar->T.Cvmin) * (itl - 1) / (itl_both_Cv - 1);
	else
	  netpar->T.Cv[itl] = 0.5 * (netpar->T.Cvmin + netpar->T.Cvmax);
      }
      else if (itl <= itl_both_Cv + itl_p_Cv)
      {
	if (netpar->T.wd == 'p')
	{
	  if (itl_p_Cv > 1)
	    netpar->T.Cv[itl] = netpar->T.Cvmin +
	    (netpar->T.Cvmax - netpar->T.Cvmin) * (itl - itl_both_Cv - 1) /
	    (itl_p_Cv - 1);
	  else
	    netpar->T.Cv[itl] = 0.5 * (netpar->T.Cvmin + netpar->T.Cvmax);
	}
	else
	  netpar->T.Cv[itl] = 0.0;
      }
      else
      {
	if (netpar->T.wd == 'r')
	{
	  if (itl_r_Cv > 1)
	    netpar->T.Cv[itl] = netpar->T.Cvmin +
	    (netpar->T.Cvmax - netpar->T.Cvmin) *
	    (itl - itl_both_Cv - itl_r_Cv - 1) / (itl_r_Cv - 1);
	  else
	    netpar->T.Cv[itl] = 0.5 * (netpar->T.Cvmin + netpar->T.Cvmax);
	}
	else
	  netpar->T.Cv[itl] = 0.0;
      }
    }
  }

  fprintf(fl.out, "\n itl   Cv\n");
  for (itl=1; itl<=netpar->T.ntl; itl++)
  {
    fprintf(fl.out, "%d %lf\n", itl, netpar->T.Cv[itl]);
  }
}

/* This function initializes the thalamic variables. The thalamic input */
/* is modeled as an inhomogeneous Poisson process.                      */
void initialize_thalamic_variables_and_spikes(spk_str *spkstr, net_par *netpar,
     run_par *runpar, par_all *parall, fl_st fl)
{
  double rand_num, phi_cal;
  int itl, ipop;

  netpar->T.phi = dvector(1, netpar->T.ntl);

  for (itl=1; itl<=netpar->T.ntl; itl++)
  {
    spkstr->spike_exist[0][itl] = 1;
  }
  
  if (netpar->T.determine_phi == 'f')       /* fixed */
  {
    for (itl=1; itl<=netpar->T.ntl; itl++)
    {
      netpar->T.phi[itl] = netpar->T.phi_read;
    }
  }
  else if (netpar->T.determine_phi == 'u') /* uniformly random drawing of phi */
  {
    for (itl=1; itl<=netpar->T.ntl; itl++)
    {
      rand_num = get_rn_dbl(parall->rng_ptr);
      netpar->T.phi[itl] = 2.0 * Pi * rand_num;
    }
  }
  else if (netpar->T.determine_phi == 'p') /* uniformly random drawing of phi */
  {
    for (itl=1; itl<=netpar->T.ntl; itl++)
    {
      rand_num = get_rn_dbl(parall->rng_ptr);
      phi_cal = solve_phi_p_sin_phi_eq(rand_num, netpar->T.Bp, runpar, fl);
      netpar->T.phi[itl] = netpar->T.phi_read + phi_cal;
      if (netpar->T.phi[itl] > 2.0) netpar->T.phi[itl] -= 2.0;
      if (netpar->T.phi[itl] < 0.0) netpar->T.phi[itl] += 2.0;
    }
    if (!runpar->sm)
    {
      fprintf(fl.out, "\n itl phi\n");
      for (itl=1; itl<=netpar->T.ntl; itl++)
        fprintf(fl.out, "%d %lf\n", itl, netpar->T.phi[itl]);
    }
  }


  for (itl=1; itl<=netpar->T.ntl; itl++)
  {
    rand_num = get_rn_dbl(parall->rng_ptr);
    if (netpar->T.nw == 'w')
    {
      spkstr->t_p_delay_spike[0][itl][1] = find_time_next_spike_w(0.0, itl, 
      rand_num, &netpar->T, runpar, fl);
    }
    else if (netpar->T.nw == 'n')
    {
      spkstr->t_p_delay_spike[0][itl][1] = find_time_next_spike_n(0.0, itl, 
      rand_num, &netpar->T, runpar, fl);
    }
    else if (netpar->T.nw == 'l')
    {
      spkstr->t_p_delay_spike[0][itl][1] = find_time_next_spike_l(0.0, itl, 
      rand_num, &netpar->T, runpar, fl);
    }
    else
    {
      printf ("nw %c should be w or n or l!", netpar->T.nw);
      exit(0);
    }

    /* x_u
    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      spkstr->x_u_delay_spike[ipop][0][itl][1] = netpar->S[ipop][0].UU;
    }
    */
  }

  if (!runpar->sm)
  {
    fprintf(fl.out, "\n itl  phi    tl_spike\n");
    for (itl=1; itl<=netpar->T.ntl; itl++)
      fprintf(fl.out, "%d %lf %lf\n", itl, netpar->T.phi[itl], 
      spkstr->t_p_delay_spike[0][itl][1]);
    fflush(fl.out);
  }
}

/* This function solves the equation:                              */
/* phi + Bp sin(phi) = rand_num                                    */
double solve_phi_p_sin_phi_eq(double rand_num, double Bp, run_par *runpar,
       fl_st fl)
{
  transfer_to_func_solve ttfuncs;
  double tol;
  double phi_solve;
  void *ptr;
  
  tol=1.0e-12;

  ttfuncs.Bp = &Bp;
  ttfuncs.rand_num = &rand_num;
  ttfuncs.fl = &fl;
  ptr = (void*) (&ttfuncs);
  phi_solve = zbrent(Bp_phi_func, 0.0, 2.0, tol, ptr);

  fprintf(fl.out, "rand_num=%lf Bp=%lf phi_sol=%lf\n", rand_num, Bp, phi_solve);

  return (phi_solve);
}

/* This function is used by zbrent to compute the solution of the equation:  */
/* phi + Bp sin(phi) = rand_num                                              */
/* 0 <= phi <= 2                                                             */
double Bp_phi_func(double phi, void *ptr)
{
  transfer_to_func_solve *ttfuncs;
  double Bp, rand_num, val_cal;
  fl_st fl;
  
  ttfuncs = (transfer_to_func_solve*) ptr;
  Bp = *(ttfuncs->Bp);
  rand_num = *(ttfuncs->rand_num);
  fl = *(ttfuncs->fl);

  val_cal = phi + Bp * sin(Pi * phi) / Pi;

  return val_cal - 2.0 * rand_num;
}

/* This function solves the differential equations. */
void n_run(func_cell_model *fcm, double ***Varbar, syn_str *synstr, 
     spk_str *spkstr, net_par *netpar, run_par *runpar, par_all *parall,
     avr_val *av, fl_st fl)
{
  double ***k0, ***k1, ***k2, ***k3, ***k4, ***Varc, ***Varold;
  double time, deltat;
  double xnoise;
  int **after_max_vol;
  int it, ipop, ion, ieq, iold;

  define_variables_arrays(netpar, runpar, &k0, fl);
  define_variables_arrays(netpar, runpar, &k1, fl);
  define_variables_arrays(netpar, runpar, &k2, fl);
  define_variables_arrays(netpar, runpar, &k3, fl);
  define_variables_arrays(netpar, runpar, &k4, fl);
  define_variables_arrays(netpar, runpar, &Varc, fl);
  define_old_variables_arrays(netpar, runpar, &Varold, &after_max_vol, fl);

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)  
    { 
      for (iold=1; iold<=2; iold++)
      {
        Varold[ipop][ion][iold] = 0.0;
      }
    }
  }

  time = 0.0;
  it = 0;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    spkstr->Vav[ipop].Vpop = 0.0;
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      spkstr->Vav[ipop].Vpop += Varbar[ipop][ion][1];      
    }
    spkstr->Vav[ipop].Vpop /= netpar->C[ipop].non;
  }

  if ((!runpar->sm) && (runpar->tmcol + runpar->epsilon >= runpar->time_all))
    pr_fct(Varbar, synstr, spkstr, netpar, runpar, time, it, fl);
  

  for (runpar->ideltat=1; runpar->ideltat<=runpar->ndeltat; runpar->ideltat++)
  {
    it=0;
    deltat = runpar->deltat[runpar->ideltat];

    while ((++it) <= runpar->nt[runpar->ideltat])
    {
      time += deltat;                   /* runpar->deltat[runpar->ideltat]; */

      /*
      if (netpar->T.nw == 'y')
        determine_thalamic_input_values(time, netpar, runpar, fl);
      */
      /* update_i_Iapp(netpar, it, fl); */

      if (netpar->inter.consider == 'y')
        compute_total_synaptic_conductance_on_a_neuron(Varbar, synstr, time, it,
        deltat, netpar, runpar, fl);

      if (runpar->method == 'r')
                                                     /* Runge-Kutta-4 method */
      {
        one_integration_step(fcm, netpar, runpar, synstr, Varbar, k0, k1, Varc,
        0.0,        time, it, fl); 
        one_integration_step(fcm, netpar, runpar, synstr, Varbar, k1, k2, Varc,
        deltat/2.0, time, it, fl);
        one_integration_step(fcm, netpar, runpar, synstr, Varbar, k2, k3, Varc,
        deltat/2.0, time, it, fl);
        one_integration_step(fcm, netpar, runpar, synstr, Varbar, k3, k4, Varc,
        deltat,     time, it, fl);

        for (ipop=1; ipop<=netpar->npop; ipop++)
        {
          for (ion=1; ion<=netpar->C[ipop].non; ion++)
          {
            for (ieq=1; ieq<=netpar->C[ipop].neq; ieq++)
  	    {
              Varbar[ipop][ion][ieq] += (deltat/6.0) * 
              (k1[ipop][ion][ieq]  + 2.0 * (k2[ipop][ion][ieq]+
               k3[ipop][ion][ieq]) + k4[ipop][ion][ieq]);
	    }
          }
        }
      }
      else if (runpar->method =='t') 
                                                     /* Runge-Kutta-2 method */
      {
        one_integration_step(fcm, netpar, runpar, synstr, Varbar, k0, k1, Varc,
        0.0,        time, it, fl);
        one_integration_step(fcm, netpar, runpar, synstr, Varbar, k1, k2, Varc,
        deltat/2.0, time, it, fl);

        for (ipop=1; ipop<=netpar->npop; ipop++)
        {
          for (ion=1; ion<=netpar->C[ipop].non; ion++)
          {
            for (ieq=1; ieq<=netpar->C[ipop].neq; ieq++)
	    {
              Varbar[ipop][ion][ieq] +=  
              deltat * k2[ipop][ion][ieq];
	    }
          }
        }
      }
      else if (runpar->method =='e')                         /* Eulear method */
      {
        one_integration_step(fcm, netpar, runpar, synstr, Varbar, k0, k1, Varc,
        0.0,        time, it, fl);

        for (ipop=1; ipop<=netpar->npop; ipop++)
        {
          for (ion=1; ion<=netpar->C[ipop].non; ion++)
          {
            for (ieq=1; ieq<=netpar->C[ipop].neq; ieq++)
  	    {
              Varbar[ipop][ion][ieq] +=
              deltat * k1[ipop][ion][ieq];
	    }
          }
        }
      }
      else
      {
        printf("wrong method!\n");
        exit(0);
      }

      /* Checking for nan */
      for (ipop=1; ipop<=netpar->npop; ipop++)
      {
        for (ion=1; ion<=netpar->C[ipop].non; ion++)
        {
          if (isnan(Varbar[ipop][ion][1]))
	  {
            printf("nan: it=%d time=%lf ipop=%d ion=%d V=%lf", it, time,
            ipop, ion, Varbar[ipop][ion][1]);
	    for (ieq=2; ieq<=netpar->C[ipop].neq; ieq++)
	      printf(" %lf", Varbar[ipop][ion][ieq]);
	    printf("\n");
            exit(0);
	  }
        }
      }

      /* Noise */
      if (netpar->noise >= runpar->epsilon)
      {
        for (ipop=1; ipop<=netpar->npop; ipop++)
        {
          for (ion=1; ion<=netpar->C[ipop].non; ion++)
          {
            xnoise = gasdev(parall);
            Varbar[ipop][ion][1] += sqrt(2.0 * netpar->noise * deltat) * xnoise;
	  }
	}
      }

      spike_detect(Varbar, Varold, after_max_vol, it, time, synstr, spkstr, 
      netpar, runpar, av, fl);

      spkstr->nspike = 0;
      if (netpar->T.thal_input == 'p')
        update_thalamic_variables(it, time, deltat, synstr, spkstr, netpar,
          runpar, parall, av, fl);

      if (netpar->process_time_delay == 's')  /* spread of time delays */
        multiple_store_spikes_plus_td(it, time, deltat, synstr, spkstr, netpar,
  	  runpar, fl);

      update_delayed_cortical_spikes(it, time, deltat, spkstr, netpar, runpar, 
        av, fl);

      decay_post_synaptic_variables(synstr, time, it, netpar, runpar, fl);

      if (netpar->process_time_delay == 'a')
        update_post_synaptic_variables_for_pre_synaptic_spikes_a(synstr, spkstr,
        time, it, deltat, netpar, runpar, fl);
      else if (netpar->process_time_delay == 's')
        update_post_synaptic_variables_for_pre_synaptic_spikes_s(synstr, spkstr,
        time, it, deltat, netpar, runpar, fl);

      for (ipop=1; ipop<=netpar->npop; ipop++)
      {
        for (ion=1; ion<=netpar->C[ipop].non; ion++)  
        { 
          Varold[ipop][ion][2] = Varold[ipop][ion][1];
          Varold[ipop][ion][1] = Varbar[ipop][ion][1];
        }
      }

      if (time >= runpar->time_all - runpar->tstat + 0.5 * deltat)
      {
	update_Vav_arrays(Varbar, spkstr->Vav, synstr, it, time, deltat, netpar, 
	runpar, fl);
      }
      else
      {
        for (ipop=1; ipop<=netpar->npop; ipop++)
        {
          spkstr->Vav[ipop].Vpop = 0.0;
          for (ion=1; ion<=netpar->C[ipop].non; ion++)
          {
            spkstr->Vav[ipop].Vpop += Varbar[ipop][ion][1];      
	  }
          spkstr->Vav[ipop].Vpop /= netpar->C[ipop].non;
	}
      }

      if ((!runpar->sm)  && !(it%runpar->twrite) &&
          (time >= runpar->time_all - runpar->tmcol - runpar->epsilon))
        pr_fct(Varbar, synstr, spkstr, netpar, runpar, time, it, fl);
    }
  }

  free_variables_arrays(netpar, runpar, k0, fl);
  free_variables_arrays(netpar, runpar, k1, fl);
  free_variables_arrays(netpar, runpar, k2, fl);
  free_variables_arrays(netpar, runpar, k3, fl);
  free_variables_arrays(netpar, runpar, k4, fl);
  free_variables_arrays(netpar, runpar, Varc, fl);
  free_old_variables_arrays(netpar, runpar, Varold, after_max_vol, fl);
}

/* This function computes one integration step */
void one_integration_step(func_cell_model *fcm, net_par *netpar,
     run_par *runpar, syn_str *synstr, double ***Varbar, double ***kin,
     double ***kout, double ***Varc, double delt, double time, int it,
     fl_st fl)
{
  double ***Varo, Iapp_now;
  int ipop, ion, ieq;

  /* Runge-Kutta input variables */
  if (delt > runpar->epsilon)
  {
    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      for (ion=1; ion<=netpar->C[ipop].non; ion++)
      {
        for (ieq=1; ieq<=netpar->C[ipop].neq; ieq++)
	{
          Varc[ipop][ion][ieq] = Varbar[ipop][ion][ieq] + delt * 
                                 kin[ipop][ion][ieq];
	}
      }
    }
    Varo = Varc;
  }
  else
  {
    Varo = Varbar;
  }

  /*
  if (it > 200 && it < 400)
  {
    fprintf(fl.out, "it=%d time=%lf delt=%lf Varo=%lf %lf S=%lf ", it, time, delt, Varo[2][150][1], Varbar[2][150][1], synstr->Isyn_cont[2][150]);
  }
  */

  /* Updating cells */

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      if (netpar->C[ipop].inject_current == 'y')
      {
        if (ion == netpar->C[ipop].ion_inject)
        {
          if (time <= netpar->C[ipop].tinject + runpar->epsilon)
	  {
            Iapp_now = netpar->C[ipop].Iinject;
	  }
          else
	  {
            Iapp_now = 0.0;
	  }
        }
        else
        {
          Iapp_now = 0.0;
        }
      }
      else if (netpar->C[ipop].inject_current == 'n')
      {
        Iapp_now = 0.0;
      }
      else
      {
        printf("inject_current=%c should be y or n!\n", 
          netpar->C[ipop].inject_current);
        exit(0);
      }

      /* if ((ipop == 1) && (ion == 1) && (delt < runpar->epsilon))
	 fprintf(fl.col, "is1=%lf\n", synstr->Isyn_cont[ipop][ion]); */

      fcm[ipop].update_cell(Varo[ipop][ion], kout[ipop][ion], ipop, ion,
      &netpar->C[ipop], &netpar->inter, &netpar->P, Iapp_now, 
      synstr->Isyn_cont[ipop][ion], synstr->Iel_cont[ipop][ion], it, fl);
    }
  }
  /*
  if (it > 200 && it <400)
  {
    fprintf(fl.out, "k=%lf\n", kout[2][150][1]);
  }
  */
}

/* This function computes the synaptic conductances on a neuron. */
void compute_total_synaptic_conductance_on_a_neuron(double ***Varbar, 
     syn_str *synstr, double time, int it, double deltat, net_par *netpar,
     run_par *runpar, fl_st fl)
{
  double fNMDA, Vpost;
  int ipop, ion, isynvar;
  int jpop, jon, jpre;
  syn_coup_par *scp;

  /* Chemical synapses */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {  
    if ((!runpar->sm) && (!runpar->sp))
      fprintf(fl.out, "ipop=%d nm=%d\n", ipop, 
      synstr->nmda_on_population[ipop]);

    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      Vpost = Varbar[ipop][ion][1];

      /* if ((!runpar->sm) && (!runpar->sp))
	   fprintf(fl.out, "ion=%d Vb=%lf\n", ion, Vpost); */

      synstr->Isyn_cont[ipop][ion] = 0.0;

      for (isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
      {
        fNMDA = (synstr->synpar[ipop][isynvar].is_it_nmda) ?
 	          Gammaf(Vpost, netpar->N.thetanp, netpar->N.sigmanp) : 1.0;

        synstr->Isyn_cont[ipop][ion] -= synstr->synpar[ipop][isynvar].gsyn * 
          synstr->synvar[ipop][ion][isynvar] *
	    (netpar->inter.rho_concur * 
             (Vpost - synstr->synpar[ipop][isynvar].Vsyn[ion])+
	     (1.0 - netpar->inter.rho_concur) *
	     (netpar->C[ipop].VL - synstr->synpar[ipop][isynvar].Vsyn[ion])) * 
             fNMDA;

        if ((ipop == 2) && (ion == 150) && (1 < 2))  
	{
	  /* if ((!runpar->sm) && (runpar->sp)) */
	  if (it > 200 && it <400)
	  {
            fprintf(fl.out, "%lf is=%d gsyn=%lf s=%lf DelV=%lf Isc=%lf"
	      " Va=%lf\n", time,
	      isynvar, synstr->synpar[ipop][isynvar].gsyn,
	      synstr->synvar[ipop][ion][isynvar], 
	      Vpost - synstr->synpar[ipop][isynvar].Vsyn[ion],
            synstr->Isyn_cont[ipop][ion], Varbar[ipop][ion][1]);
	  }
        }
      }
    }
  }

  /* Electrical synapses */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      synstr->Iel_cont[ipop][ion] = 0.0;
    }

    if (netpar->pop_name[ipop] == 'P')
    {
      jpop = ipop;
      scp = &(netpar->S[ipop][jpop]);

      for (ion=1; ion<=netpar->C[ipop].non; ion++)
      {
       Vpost = Varbar[ipop][ion][1];

       synstr->el_coef_V[ipop][ion] = 0.0;

       for (jpre=1; jpre<=scp->nw_e_coup[ion]; jpre++)
       {
         jon = scp->w_e_coup[ion][jpre];
         synstr->el_coef_V[ipop][ion] += Vpost - Varbar[jpop][jon][1];
       }

       synstr->Iel_cont[ipop][ion] = -synstr->gel_one_syn * 
         synstr->el_coef_V[ipop][ion];
       /*       printf("ipop=%d ion=%d gel=%lf DelV=%lf Iel=%lf\n", ipop, ion, synstr->gel_one_syn, synstr->el_coef_V[ipop][ion], synstr->Iel_cont[ipop][ion]); */
      }
    }
  }
}

/* This function writes the data as functions of time.  */
void pr_fct(double ***Varbar, syn_str *synstr, spk_str *spkstr, net_par *netpar,
     run_par *runpar, double time, int it, fl_st fl)
{
  int ipop, ion, iwrite, ieq, isynvar;
  int icol_wrt, cond;

  cond =  (it == 0) ||
          (fabs(runpar->tmcol - runpar->time_all) <= runpar->epsilon);
  icol_wrt = 1;
  if (cond) fprintf(fl.out, "\n%d time\n", icol_wrt);
  fprintf(fl.col, "%12.5lf", time);

  /*
  fprintf(fl.col, " %d", spkstr->nspike_pop[0]);
  */

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    if (cond) fprintf(fl.out, "%d nspike_pop %d\n", ++icol_wrt, ipop);
    fprintf(fl.col, " %d", spkstr->nspike_pop[ipop]);
    if (cond) fprintf(fl.out, "%d Vav %d\n", ++icol_wrt, ipop);
    fprintf(fl.col, " %lf", spkstr->Vav[ipop].Vpop);

    for (iwrite=1; iwrite<=runpar->nwrite[ipop]; iwrite++)
    {
      ion = runpar->nwritear[ipop][iwrite];
      if (cond) fprintf(fl.out, "%d V(%d,%d)\n", ++icol_wrt, ipop, ion);
      fprintf(fl.col, " %12.5lf", Varbar[ipop][ion][1]);

      if (runpar->write_aux == 'y')
      {
        for (ieq=2; ieq<=netpar->C[ipop].nceq; ieq++)
	{
          if (cond) fprintf(fl.out, "%d V(%d,%d,%d)\n", ++icol_wrt, ipop,
            ion, ieq);
          fprintf(fl.col, " %12.5lf", Varbar[ipop][ion][ieq]);
	}
      }
      else if (runpar->write_aux != 'n')
      {
        printf("write_aux=%c should be y or n!\n", runpar->write_aux);
        exit(0);
      }

      if (synstr->nsynvar[ipop] >= 1)
      {
        if (cond)
          fprintf(fl.out, "%d S(%d,%d,1) %c%c\n", ++icol_wrt, ipop, ion,
	  synstr->csynvar[ipop][1][0], synstr->csynvar[ipop][1][1]);
        fprintf(fl.col, " %14.8lf", synstr->synvar[ipop][ion][1]);

        if ((2 > 1)) /* (runpar->write_aux == 'y') || */
        {
          for (isynvar = 2; isynvar<=synstr->nsynvar[ipop]; isynvar++)
	  {
            if (cond)
              fprintf(fl.out, "%d S(%d,%d,%d) %c%c\n", ++icol_wrt, ipop, ion,
	      isynvar, synstr->csynvar[ipop][isynvar][0],
              synstr->csynvar[ipop][isynvar][1]);
            fprintf(fl.col, " %12.5lf", synstr->synvar[ipop][ion][isynvar]);
	    if (synstr->synvar[ipop][ion][isynvar] > 10000.0)
	    {
              fprintf(fl.col, "large synvar! ipop=%d ion=%d isynvar=%d\n", ipop, ion, isynvar);
	      printf("large synvar! ipop=%d ion=%d isynvar=%d\n", ipop, ion, isynvar);
	      exit(0);
	    }
	  }
        }
      }
    }
  }
  if (cond) fprintf(fl.out, "  \n");
  fprintf(fl.col, "\n");
  fflush(fl.col);
}

/* This function finds the spikes fired by pre-synaptic neurons and update  */
/* the synaptic variables of the post-synaptic neurons.                     */
void spike_detect(double ***Varbar, double ***Varold, int **after_max_vol,
     int it, double time, syn_str *synstr, spk_str *spkstr, net_par *netpar,
     run_par *runpar, avr_val *av, fl_st fl)
{
  double after_min_vol=-33.0, xb, xc, tpeak, Vpeak;
  double V0, V1, V2;
  double diff_fire_time, tspk_relative_to_per, tt_mod_Tper;
  int ipop, jpop, jon, nspk, ihist;
  int spike_detected;
  char spike_detect_criterion;

  spike_detect_criterion = 't';
  /* spkstr->nspike = 0; */
  spkstr->nsp_in_dt = 0;
  
  for (jpop=1; jpop<=netpar->npop; jpop++)
    spkstr->nspike_pop[jpop] = 0;

  for (jpop=1; jpop<=netpar->npop; jpop++)
  {
    for (jon=1; jon<=netpar->C[jpop].non; jon++)  
    {
      V0 = Varbar[jpop][jon][1];
      if ((!after_max_vol[jpop][jon]) && it>=10)
      {
        V1 = Varold[jpop][jon][1];
 
        if (spike_detect_criterion == 'p')
	{
          V2 = Varold[jpop][jon][2];
          spike_detected = spike_detect_peak(V0, V1, V2, time, &tpeak, &Vpeak,
          netpar, runpar, fl);
	}
        else if (spike_detect_criterion == 't')
	{
          spike_detected = spike_detect_threshold(V0, V1, time, &tpeak, &Vpeak,
          netpar, runpar, fl);
	}
        else
	{
          printf("spike_detect_criterion=%c should be p or t\n",
            spike_detect_criterion);
          exit(0);
	}

	
        if (spike_detected)
	{
          after_max_vol[jpop][jon] = 1;
          if (netpar->process_time_delay == 'a')
	  {
	    /* only one value of tau_delay for all synapses */
	    single_td_storing_spikes(jpop, jon, tpeak, it, time, spkstr, netpar,
	    runpar, fl);
	  }
	  else if (netpar->process_time_delay == 's')
	  {
	    /* spread of tau_delay values among synapses */
	    multiple_td_storing_spikes(jpop, jon, tpeak, it, time, spkstr,
	    netpar, runpar, fl);
	  }
	  
	  /*
          for (ipop=1; ipop<=netpar->npop; ipop++)
          {
            spkstr->x_u_delay_spike[ipop][jpop][jon][nspk] = 
            netpar->S[ipop][jpop].UU * synstr->xvar[ipop][jpop][jon];
          }
	  */

          short_term_plasticity(it, time, synstr->xvar, 
          spkstr->tfire_prev[jpop][jon], tpeak, jpop, jon, netpar, 
          runpar, fl);

          /* store spike data */
          if (time >= runpar->time_all - runpar->tstat + runpar->epsilon)
	  {
            spkstr->nfire[jpop][jon]++;
	    tspk_relative_to_per = tpeak - 
              ((int) (tpeak / netpar->T.Tper) * netpar->T.Tper);
            spkstr->Z1cos[jpop][jon] += cos(2.0 * Pi * tspk_relative_to_per / 
	      netpar->T.Tper);
            spkstr->Z1sin[jpop][jon] += sin(2.0 * Pi * tspk_relative_to_per / 
	      netpar->T.Tper);

	    ihist = (int) (tspk_relative_to_per * runpar->nhist /
	    netpar->T.Tper);

  	    if (jpop <= 1)
	    {
	      spkstr->fr_hist[jpop][ihist]++;
	    }
	    else
	    {
	      if (jon <= 
                  (int) (netpar->C[jpop].fracIext * netpar->C[jpop].non + 
		  runpar->epsilon))
	      {
                spkstr->fr_hist[jpop][ihist]++;
	      }
	      else
	      {
                spkstr->fr_hist[jpop+1][ihist]++;
	      }
	    }

	    if ((tt_mod_Tper = t_mod(tspk_relative_to_per - netpar->T.tc,
	      netpar->T.Tper)) < runpar->t_touch_interval)
	    {
  	      if (jpop <= 1)
	      {
	        spkstr->spk_touch[jpop][jon]++;
	      }
	      else
	      {
	        if (jon <= 
                    (int) (netpar->C[jpop].fracIext * netpar->C[jpop].non + 
	  	    runpar->epsilon))
	        {
                  spkstr->spk_touch[jpop][jon]++;
	        }
	        else
	        {
                  spkstr->spk_touch[jpop][jon]++;
	        }
	      }	      
	    }
	    else if (t_mod(netpar->T.tc - tspk_relative_to_per,
	      netpar->T.Tper) < runpar->t_touch_interval)
	    {
  	      if (jpop <= 1)
	      {
	        spkstr->spk_before_touch[jpop][jon]++;
	      }
	      else
	      {
	        if (jon <= 
                    (int) (netpar->C[jpop].fracIext * netpar->C[jpop].non + 
	  	    runpar->epsilon))
	        {
                  spkstr->spk_before_touch[jpop][jon]++;
	        }
	        else
	        {
                  spkstr->spk_before_touch[jpop][jon]++;
	        }
	      }	      
	    }


            diff_fire_time = -1.0;

            if (spkstr->nfire[jpop][jon] > 1)
	    {
              diff_fire_time = tpeak - spkstr->tfire_prev[jpop][jon];
              spkstr->av_deltat_fire[jpop][jon] += diff_fire_time;
              spkstr->av_deltat_fire_sq[jpop][jon] += diff_fire_time * 
                diff_fire_time;
	    }
            spkstr->tfire_prev[jpop][jon] = tpeak;
	  }

          if (!runpar->sm)
  	  {
            fprintf(fl.ras, "%14.8e %d %d %lf %d\n", tpeak, jpop, jon, Vpeak, 
            spkstr->nspike_pop[jpop]);
	    /*
            fprintf(fl.ras, "%14.8e %d %d %lf %d\n", 
              spkstr->tspike[spkstr->nspike], spkstr->jpop[spkstr->nspike],
	      spkstr->jon[spkstr->nspike], Vpeak, spkstr->nspike);
	    */
            fflush(fl.ras);
	  }
	}
      }
      else
      {
        if (V0 < after_min_vol) after_max_vol[jpop][jon] = 0;
      }
    }
  }

  /*
  for (jpop=1; jpop<=netpar->npop; jpop++) 
    fprintf(fl.col, "jpop=%d nspike=%d\n", jpop, spkstr->nspike_pop[jpop]);
  */
}

/* This function detects a peak of the membrane potential. */
int spike_detect_peak(double V0, double V1, double V2, double time,
    double *tpeak, double *Vpeak, net_par *netpar, run_par *runpar, fl_st fl)
{
  double xb, xc;
  int spike_detected;

  if (V1 >= V0 && V1 >= V2 && V1 > netpar->Volt_thresh)  /* detect a spike */
  {
    spike_detected = 1;
    xb = V2 - V0;
    xc = V0 - 2.0 * V1 + V2;
    if (fabs(xc) < runpar->epsilon)
    {
      *tpeak = time - runpar->deltat[runpar->ideltat];
      *Vpeak = V1;
    }
    else        
    {
      *tpeak = time - runpar->deltat[runpar->ideltat] + 
               0.5 * (xb / xc) * runpar->deltat[runpar->ideltat];
      *Vpeak = V1 - 0.125 * xb * xb / xc;
    }
  }
  else
  {
    spike_detected = 0;
  }

  return(spike_detected);
}

/* This function detects a threshold crossing of the membrane potential. */
int spike_detect_threshold(double V0, double V1, double time, double *tpeak,
    double *Vpeak, net_par *netpar, run_par *runpar, fl_st fl)
{
  int spike_detected;

  if ((V0 >= netpar->Volt_thresh) && (V1 < netpar->Volt_thresh))  
    /* detect a spike */
  {
    spike_detected = 1;
    *tpeak = lininter(V1, V0, netpar->Volt_thresh, 
             time - runpar->deltat[runpar->ideltat], time);
    *Vpeak = netpar->Volt_thresh;
  }
  else
  {
    spike_detected = 0;
  }

  return(spike_detected);
}

/* This function stores newly-detected spikes in the aray spkstr->spike_exist */
/* It is used for cases in which all the synaptic delays tau_delay are equal. */
void single_td_storing_spikes(int jpop, int jon, double tpeak, int it,
     double time, spk_str *spkstr, net_par *netpar, run_par *runpar, fl_st fl)
{
  int nspk;
  
  spkstr->nspike_pop[jpop]++;

  spkstr->spike_exist[jpop][jon]++;
  if (spkstr->spike_exist[jpop][jon] > spkstr->mspk)
  {
    printf("jpop=%d jon=%d spike_exist=%d > mspk=%d\n", jpop, jon,
    spkstr->spike_exist[jpop][jon], spkstr->mspk);
    exit(0);
  }
  else
  {
    nspk = spkstr->spike_exist[jpop][jon];
    spkstr->t_p_delay_spike[jpop][jon][nspk] = tpeak + 
      netpar->S[1][jpop].tau_delay;  
  }

  /*
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    spkstr->x_u_delay_spike[ipop][jpop][jon][nspk] = 
      netpar->S[ipop][jpop].UU * synstr->xvar[ipop][jpop][jon];
  }
  */
}

/* This function stores newly-detected spikes. in the array spkstr->sp_in_dt. */
/* It is used for cases in which the values of synaptic delays tau_delay vary */
/* among synapses.                                                            */
void multiple_td_storing_spikes(int jpop, int jon, double tpeak, int it,
     double time, spk_str *spkstr, net_par *netpar, run_par *runpar, fl_st fl)
{
  spkstr->nspike_pop[jpop]++;

  spkstr->nsp_in_dt++;
  spkstr->sp_in_dt[spkstr->nsp_in_dt].tspike = tpeak;
  spkstr->sp_in_dt[spkstr->nsp_in_dt].jpop = jpop;
  spkstr->sp_in_dt[spkstr->nsp_in_dt].jon = jon;
}

/* This function finds if there are spikes fired by thalamiuc neurons   */
/* in this time interval deltat. If there are, the variable nspike      */
/* is updated and the variable ion is set to the number of the neuron   */
/* that fired.                                                          */
void update_thalamic_variables(int it, double time, double deltat, 
     syn_str *synstr, spk_str *spkstr, net_par *netpar, run_par *runpar,
     par_all *parall, avr_val *av, fl_st fl)
{
  double rand_num, tspk_relative_to_per;
  double diff_fire_time;
  int ipop, jtl, ihist;

  spkstr->nspike_pop[0] = 0;

  for (jtl=1; jtl<=netpar->T.ntl; jtl++)
  {
    if ((spkstr->t_p_delay_spike[0][jtl][1] >= time - deltat) && 
        (spkstr->t_p_delay_spike[0][jtl][1] < time))
    {
      spkstr->nspike_pop[0]++;

      spkstr->nspike++;
      spkstr->jpop[spkstr->nspike] = 0;
      spkstr->jon[spkstr->nspike] = jtl;
      spkstr->tspike[spkstr->nspike] = spkstr->t_p_delay_spike[0][jtl][1];

      if (netpar->process_time_delay == 's')
        multiple_td_storing_spikes(0, jtl, spkstr->t_p_delay_spike[0][jtl][1],
          it, time, spkstr, netpar, runpar, fl);

      /*
      for (ipop=1; ipop<=netpar->npop; ipop++)
      {
        spkstr->x_u_delay_spike[ipop][0][jtl][1] = 
          netpar->S[ipop][0].UU * synstr->xvar[ipop][0][jtl];
      }
      */
      
      short_term_plasticity(it, time, synstr->xvar, spkstr->tfire_prev[0][jtl],
      spkstr->t_p_delay_spike[0][jtl][1], 0, jtl, netpar, runpar, fl);

      /* store spike data */
      if (spkstr->t_p_delay_spike[0][jtl][1] >= 
          runpar->time_all - runpar->tstat + runpar->epsilon)
      {
        spkstr->nfire[0][jtl]++;
	tspk_relative_to_per = spkstr->t_p_delay_spike[0][jtl][1] - 
          ((int) (spkstr->t_p_delay_spike[0][jtl][1] / netpar->T.Tper) * 
          netpar->T.Tper);
        spkstr->Z1cos[0][jtl] += cos(2.0 * Pi * tspk_relative_to_per / 
	      netpar->T.Tper);
        spkstr->Z1sin[0][jtl] += sin(2.0 * Pi * tspk_relative_to_per / 
	      netpar->T.Tper);

        ihist = (int) (tspk_relative_to_per * runpar->nhist / netpar->T.Tper);
        spkstr->fr_hist[0][ihist]++;

        if (t_mod(tspk_relative_to_per - netpar->T.tc, netpar->T.Tper) <
	      runpar->t_touch_interval)
	{
	  spkstr->spk_touch[0][jtl]++;
	}
        else if (t_mod(netpar->T.tc - tspk_relative_to_per,
	      netpar->T.Tper) < runpar->t_touch_interval)
	{
	  spkstr->spk_before_touch[0][jtl]++;
	}

        if (spkstr->nfire[0][jtl] > 1)
        {
          diff_fire_time = spkstr->t_p_delay_spike[0][jtl][1] - 
            spkstr->tfire_prev[0][jtl];
          spkstr->av_deltat_fire[0][jtl] += diff_fire_time;
          spkstr->av_deltat_fire_sq[0][jtl] += diff_fire_time * 
            diff_fire_time;
        }
        spkstr->tfire_prev[0][jtl] = spkstr->t_p_delay_spike[0][jtl][1];
      }

      if (!runpar->sm)
      {
        fprintf(fl.ras, "%14.8e %d %d 65.0 %d\n", 
          spkstr->tspike[spkstr->nspike], spkstr->jpop[spkstr->nspike],
          spkstr->jon[spkstr->nspike], spkstr->nspike);
        fflush(fl.ras);
      }

      /* Generating a new thalamic spike */
      rand_num = get_rn_dbl(parall->rng_ptr);
      if (netpar->T.nw == 'w')
      {
        spkstr->t_p_delay_spike[0][jtl][1] = find_time_next_spike_w(
        spkstr->t_p_delay_spike[0][jtl][1], jtl, rand_num, &netpar->T, runpar,
        fl);
      }
      else if (netpar->T.nw == 'n')
      {
        spkstr->t_p_delay_spike[0][jtl][1] = find_time_next_spike_n(
        spkstr->t_p_delay_spike[0][jtl][1], jtl, rand_num, &netpar->T, runpar,
        fl);
      }
      else if (netpar->T.nw == 'l')
      {
        spkstr->t_p_delay_spike[0][jtl][1] = find_time_next_spike_l(
        spkstr->t_p_delay_spike[0][jtl][1], jtl, rand_num, &netpar->T, runpar,
        fl);
      }

      if (spkstr->t_p_delay_spike[0][jtl][1] < time + runpar->epsilon) 
	spkstr->t_p_delay_spike[0][jtl][1] = time + runpar->epsilon;
    }
  }  

  /* 15/7/2014
  if (!runpar->sm)
  {
    for (jtl = spkstr->nspike - spkstr->nspike_pop[0] + 1; 
         jtl<= spkstr->nspike; jtl++)
      fprintf(fl.out, "new tl spike %lf %d %d %d %lf\n", time, jtl, 
      spkstr->jon[jtl], spkstr->jpop[jtl], spkstr->tspike[jtl]);
  }
  */
}

/* This function takes all the spikes fired within the time step deltat,      */
/* adds the specific synaptic delay for each pair of pre- and post-synaptic   */
/* neurons, and stores the existence of the coming spike to the post-synaptic */
/* neuron, at the time that is the sum of spike time + tau_delay, in the      */
/* array time_spike_delay.                                                    */
void multiple_store_spikes_plus_td(int it, double time, double deltat, 
     syn_str *synstr, spk_str *spkstr, net_par *netpar, run_par *runpar,
     fl_st fl)
{
  double tspk, tspike_arrive;
  int ispk, ipop, jpop, ion, jon, iwcoup, idelay_a, idelay_b, idelay;

  /*
  if ((!runpar->sm) && (!runpar->sp))
  {
    for (ispk=1; ispk<=spkstr->nsp_in_dt; ispk++)
    {
      fprintf(fl.ras, "%lf %d %d --\n", spkstr->sp_in_dt[ispk].tspike,
      spkstr->sp_in_dt[ispk].jpop, spkstr->sp_in_dt[ispk].jon);
    }
  }
  */

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      spkstr->iptr_delay[ipop][jpop] = (spkstr->iptr_delay[ipop][jpop] + 1);
      while (spkstr->iptr_delay[ipop][jpop] > spkstr->ndelay[ipop][jpop])
      {
	spkstr->iptr_delay[ipop][jpop] -= spkstr->ndelay[ipop][jpop];
      }
      /*
      printf("ipop=%d jpop=%d iptr_delay=%d ndelay=%d\n", ipop, jpop, 
      spkstr->iptr_delay[ipop][jpop], spkstr->ndelay[ipop][jpop]);
      */
    }
  }

  /*
  if (spkstr->nsp_in_dt > 0)
     printf("it=%d nsp_in_dt=%d\n", it, spkstr->nsp_in_dt);
  */
  for (ispk=1; ispk<=spkstr->nsp_in_dt; ispk++)
  {
    tspk = spkstr->sp_in_dt[ispk].tspike;
    jpop = spkstr->sp_in_dt[ispk].jpop;
    jon = spkstr->sp_in_dt[ispk].jon;

    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      if (synstr->nonzero_gsyn[ipop][jpop])
      {
        for (iwcoup=1; iwcoup<=netpar->S[ipop][jpop].nwcoup[jon]; iwcoup++)
        {
          ion = netpar->S[ipop][jpop].wcoup[jon][iwcoup];
          tspike_arrive = tspk + netpar->S[ipop][jpop].tdcoup[jon][iwcoup];
	  idelay_a = (int) ((deltat + tspike_arrive - time + runpar->epsilon) /
	    deltat);
	  idelay_b = idelay_a + spkstr->iptr_delay[ipop][jpop];
          idelay = (idelay_b <= spkstr->ndelay[ipop][jpop])? idelay_b :
	    idelay_b - spkstr->ndelay[ipop][jpop];
	  spkstr->time_spike_delay[ipop][jpop][ion][idelay] +=
	    synstr->xvar[ipop][jpop][jon] * netpar->S[ipop][jpop].UU;

	  
	  /* 
          if ((ipop == 2) && (jpop == 1) && (jon == 485) && (ion == 40))
	  if ((ipop == 2) && (jpop == 0) && (jon == 196) && (ion == 86))
	  {
	    printf("it=%d time=%lf ipop=%d ion=%d nonzero_gsyn=%d nwcoup=%d\n",
	    it, time, ipop, ion, synstr->nonzero_gsyn[ipop][jpop],
	    netpar->S[ipop][jpop].nwcoup[jon]);
   	    printf("iwcoup=%d tspike_arrive=%lf t-t=%lf idelay a b 0=%d %d %d"
	    "\n", iwcoup, tspike_arrive, tspike_arrive - time, idelay_a,
	    idelay_b, idelay);
	    printf("iptr_delay=%d time_spike_delay=%lf\n",
	    spkstr->iptr_delay[ipop][jpop],
	    spkstr->time_spike_delay[ipop][jpop][ion][idelay]);
	    printf("ndelay=%d idelay_b=%d idelay=%d\n",
 	    spkstr->ndelay[ipop][jpop], idelay_b, idelay);
	  }
	  */
	}
      }
    }
  }
}

/* This function cpmputes the new values of the short-term plasticity    */
/* variables according to the Tsodyks-Markram dynamics.                  */
void short_term_plasticity(int it, double time, double ***xvar, double tprev,
     double tnow, int jpop, int jon, net_par *netpar, run_par *runpar, 
     fl_st fl)
{
  double xold, expDelt, onemU;
  int ipop;

  if (tprev >= 0.0)
  {
    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      xold = xvar[ipop][jpop][jon];
      expDelt = exp(-(tnow - tprev) / netpar->S[ipop][jpop].taur);
      onemU = 1.0 - netpar->S[ipop][jpop].UU;
      xvar[ipop][jpop][jon] = 1.0 -  expDelt + xold * onemU * expDelt;

      if ((!runpar->sm) && (!runpar->sp))
      {
	if ((netpar->pop_name[ipop] == 'P') && (netpar->pop_name[jpop] == 'P'))
	{
	  fprintf (fl.zmp, "%lf %d %d %d %lf %lf\n", time, ipop, jpop, jon,
	  xold, xvar[ipop][jpop][jon]);
	}
      }

      /*
      fprintf(fl.out, "ipop=%d jpop=%d jon=%d tprev=%lf tnow=%lf xold=%lf "
      "xnew=%lf\n", ipop, jpop, jon, tprev, tnow, xold, xvar[ipop][jpop][jon]);
      */
    }
  }
}

/* This function finds if there are spikes fired by cortical neurons   */
/* in this time interval deltat. If there are, the variable nspike      */
/* is updated and the variable ion is set to the number of the neuron   */
/* that fired.                                                          */
void update_delayed_cortical_spikes(int it, double time, double deltat, 
     spk_str *spkstr, net_par *netpar, run_par *runpar, avr_val *av, fl_st fl)
{
  int jpop, jon, ispk;

  for (jpop=0; jpop<=netpar->npop; jpop++)
  {
    for (jon=1; jon<=netpar->C[jpop].non; jon++)
    {
      if (spkstr->spike_exist[jpop][jon] >= 1)
      {
        if ((spkstr->t_p_delay_spike[jpop][jon][1] >= time - deltat) && 
          (spkstr->t_p_delay_spike[jpop][jon][1] < time))
        {
          spkstr->nspike++;
          spkstr->jpop[spkstr->nspike] = jpop;
          spkstr->jon[spkstr->nspike] = jon;
          spkstr->tspike[spkstr->nspike] = 
            spkstr->t_p_delay_spike[jpop][jon][1];

	  spkstr->spike_exist[jpop][jon] -= 1;

	  if (spkstr->spike_exist[jpop][jon] >= 1)
	  {
	    /*
            fprintf(fl.out, "\nstored spikes: it=%d jpop=%d jon=%d nspk=%d\n",
	      it, jpop, jon, spkstr->spike_exist[jpop][jon]);
            fprintf(fl.out, "before: ");
            for (ispk=1; ispk<=spkstr->spike_exist[jpop][jon]+1; ispk++)
              fprintf(fl.out, " %lf", spkstr->t_p_delay_spike[jpop][jon][ispk]);
            fprintf(fl.out, "\n");
	    */
            for (ispk=1; ispk<=spkstr->spike_exist[jpop][jon]; ispk++)
	    {
              spkstr->t_p_delay_spike[jpop][jon][ispk] = 
              spkstr->t_p_delay_spike[jpop][jon][ispk+1];
	    }
	    /*
            fprintf(fl.out, "after: ");
            for (ispk=1; ispk<=spkstr->spike_exist[jpop][jon]; ispk++)
              fprintf(fl.out, " %lf", spkstr->t_p_delay_spike[jpop][jon][ispk]);
            fprintf(fl.out, "\n");
	    */
	  }
	}
      }
    }
  }
}

/* This function compute the decay of synaptic variables.   */
void decay_post_synaptic_variables(syn_str *synstr, double time, int it,
     net_par *netpar, run_par *runpar, fl_st fl)
{
  int ipop, ion, isynvar;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      for(isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
      {
        synstr->synvar[ipop][ion][isynvar] *= 
          exp(-runpar->deltat[runpar->ideltat] /
          synstr->synpar[ipop][isynvar].tsyn);
      } 
    }
  }
}

/* This function updates the synaptic variables of post-synaptic neurons      */
/* in response to firing of pre-synaptic neurons, for process_time_delay='a'. */
void update_post_synaptic_variables_for_pre_synaptic_spikes_a(syn_str *synstr, 
     spk_str *spkstr, double time, int it, double deltat, net_par *netpar, 
     run_par *runpar, fl_st fl)
{
  double difft, tsyn;
  int ispike, ipop, jpop, ion, jon, iwcoup, isynvar;

  for (ispike=1; ispike<=spkstr->nspike; ispike++)
  {
    difft = time - spkstr->tspike[ispike];
    jpop = spkstr->jpop[ispike];              /* pre-synaptic population */
    jon = spkstr->jon[ispike];                /* pre-synaptic neuron     */
    if ((!runpar->sm) && (!runpar->sp)) 
      fprintf(fl.out, "ispike=%d jpop=%d jon=%d tspike=%lf time=%lf difft=%lf\n"
      , ispike, jpop, jon, spkstr->tspike[ispike], time, difft);

    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      /*
      printf("ipop=%d jpop=%d jon=%d\n", ipop, jpop, jon);
      printf("nonzero_gsyn=%d nwcoup=%d\n", synstr->nonzero_gsyn[ipop][jpop], netpar->S[ipop][jpop].nwcoup[jon]);
      */

      if (synstr->nonzero_gsyn[ipop][jpop])
      {
        for (iwcoup=1; iwcoup<=netpar->S[ipop][jpop].nwcoup[jon]; iwcoup++)
        {
          ion = netpar->S[ipop][jpop].wcoup[jon][iwcoup];
          for (isynvar= synstr->isyn_to_send[ipop][jpop][1]; 
               isynvar<=synstr->isyn_to_send[ipop][jpop][2]; isynvar++)
	  {
	    /*
            if ((!runpar->sm) && (!runpar->sp)) 
              fprintf(fl.out, "ipop=%d iw=%d ion=%d isv=%d sb=%lf\n", ipop, 
              iwcoup, ion, isynvar, synstr->synvar[ipop][ion][isynvar]);
	    */

            tsyn = synstr->synpar[ipop][isynvar].tsyn;
	    /* ----- */
	    /* synstr->synvar[ipop][ion][isynvar] += exp(-difft / tsyn); */
	    synstr->synvar[ipop][ion][isynvar] += 
              synstr->xvar[ipop][jpop][jon] * netpar->S[ipop][jpop].UU;

	    /*	    
	    if ((ipop == 2) && (ion == 2) && (isynvar == 1))
	      printf("it=%d ipop=%d ion=%d isynvar=%d synvar=%lf jpop=%d "
              "jon=%d xvar=%lf UU=%lf\n", it, ipop, ion, isynvar,
              synstr->synvar[ipop][ion][isynvar], jpop, jon,
              synstr->xvar[ipop][jpop][jon], netpar->S[ipop][jpop].UU);  
	    */

	    /*
            if ((!runpar->sm) && (!runpar->sp)) 
              fprintf(fl.out, "tsyn=%lf sa=%lf\n", tsyn, 
              synstr->synvar[ipop][ion][isynvar]);
	    */
	  }
	}
      }
    }
  }
}

/* This function updates the synaptic variables of post-synaptic neurons      */
/* in response to firing of pre-synaptic neurons, for process_time_delay='s'. */
void update_post_synaptic_variables_for_pre_synaptic_spikes_s(syn_str *synstr, 
     spk_str *spkstr, double time, int it, double deltat, net_par *netpar, 
     run_par *runpar, fl_st fl)
{
  int ipop, jpop, ion, idelay, isynvar;

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (jpop=0; jpop<=netpar->npop; jpop++)
    {
      idelay = spkstr->iptr_delay[ipop][jpop];
      /* printf("it=%d ipop=%d jpop=%d idelay=%d\n", it, ipop, jpop, idelay); */

      for (ion=1; ion<=netpar->C[ipop].non; ion++)
      {
	if (spkstr->time_spike_delay[ipop][jpop][ion][idelay] > runpar->epsilon)
        {
          /*
          if ((ipop == 1) && (jpop == 0) && (ion == 485))
	  {
            
            printf("it=%d tds=%lf idelay=%d ipop=%d ion=%d jpop=%d\n", it,
            spkstr->time_spike_delay[ipop][jpop][ion][idelay], idelay, ipop,
            ion, jpop);
	    printf("isyn_to_send=%d %d\n", synstr->isyn_to_send[ipop][jpop][1],
            synstr->isyn_to_send[ipop][jpop][2]);
            
          }
          */
          for (isynvar= synstr->isyn_to_send[ipop][jpop][1]; 
               isynvar<=synstr->isyn_to_send[ipop][jpop][2]; isynvar++)
	  {
	    synstr->synvar[ipop][ion][isynvar] += 
	      spkstr->time_spike_delay[ipop][jpop][ion][idelay];
	    /*
	    if ((ipop == 1) && (jpop == 1) && (it == 176) && (ion == 7))
  	    {
              printf("isynvar=%d synvar=%lf\n", isynvar,
              synstr->synvar[ipop][ion][isynvar]);
            }
	    */
	  }
	  spkstr->time_spike_delay[ipop][jpop][ion][idelay] = 0.0;
	}
      }
    }
  }
}

/* This function updates the arranys of Vav every time step.      */
/* The arrays are needed to compute the synchrony measure chi.    */
void update_Vav_arrays(double ***Varbar, V_aver *Vav, syn_str *synstr, int it,
     double time, double deltat, net_par *netpar, run_par *runpar, fl_st fl)
{
  double xx;
  int ipop, ion, isynvar; 

  /* updating Vav */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    Vav[ipop].Vpop = 0.0;

    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      Vav[ipop].V_avt[ion] += Varbar[ipop][ion][1] * deltat;
      Vav[ipop].V_sq_avt[ion] += pow(Varbar[ipop][ion][1], 2.0) * deltat;
      Vav[ipop].Vpop += Varbar[ipop][ion][1];      
    }

    Vav[ipop].Vpop /= netpar->C[ipop].non;

    Vav[ipop].Vpop_avt += Vav[ipop].Vpop  * deltat;
    Vav[ipop].Vpop_sq_avt += pow(Vav[ipop].Vpop, 2.0) * deltat;

    /*
    fprintf(fl.col, "Vpop=%lf Vpop_avt=%lf\n", Vav[ipop].Vpop, 
    Vav[ipop].Vpop_avt);
    */
  }

  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      for (isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
      {
        synstr->synvar_avt[ipop][ion][isynvar] +=
	  synstr->synvar[ipop][ion][isynvar] * deltat;
        synstr->Isynvar_avt[ipop][ion][isynvar] +=
	  synstr->synvar[ipop][ion][isynvar] * deltat *
	  (Varbar[ipop][ion][1] - synstr->synpar[ipop][isynvar].Vsyn[ion]);
      }
    }
  }
}

/* This function computes spike statistics.   */
void compute_spike_statistics(spk_str *spkstr, net_par *netpar, run_par *runpar,
     avr_val *av, fl_st fl)
{
  double var_deltat_fire, av_deltat_fire_pop_sq, cv_deltat_fire_pop_sq;
  double xpos;
  double *Z1cos_av, *Z1sin_av, *Z1md_av, *Z1phi_av, *nfire_av;
  double x_factor, fr_x_hist;
  double *fr_pop_sq, fr_pop_diff;
  int ipop, ion, non_with_cv, ihist;
  int non_subpop;

  nfire_av = dvector(0, netpar->npop);
  Z1cos_av = dvector(0, netpar->npop);
  Z1sin_av = dvector(0, netpar->npop);
  Z1md_av  = dvector(0, netpar->npop);
  Z1phi_av = dvector(0, netpar->npop);
  fr_pop_sq = dvector(0, netpar->npop+1);

  if (runpar->tstat < runpar->epsilon)
  {
    printf("tstat = 0!\n");
    exit(0);
  }

  if (!runpar->sm)
  {
    int sumn;

    fprintf(fl.out, "Tn\n");
    sumn = 0;
    for (ion=1; ion<=netpar->C[0].non; ion++)
    {
      sumn += spkstr->nfire[0][ion];
      fprintf(fl.out, "Tnfire %d %d\n",ion, spkstr->nfire[0][ion]);
    }
    fprintf(fl.out, "sumn=%d\n", sumn);
  }

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)  
    {
      spkstr->frinter[ipop][ion] = 1000.0 * spkstr->nfire[ipop][ion] /
        (runpar->tstat);
      spkstr->spk_touch[ipop][ion] *= netpar->T.Tper / runpar->tstat; 
      spkstr->spk_before_touch[ipop][ion] *= netpar->T.Tper / runpar->tstat;
     
      /* Computing the time-average firing rate and CV */
      if (spkstr->nfire[ipop][ion] >= 3)
      {
        spkstr->av_deltat_fire[ipop][ion] /= 
          spkstr->nfire[ipop][ion] - 1;
        spkstr->av_deltat_fire_sq[ipop][ion] /=
          spkstr->nfire[ipop][ion] - 1;
        var_deltat_fire = spkstr->av_deltat_fire_sq[ipop][ion] - 
          spkstr->av_deltat_fire[ipop][ion] * spkstr->av_deltat_fire[ipop][ion];

        if ((var_deltat_fire >= 0.0) && 
            (spkstr->av_deltat_fire[ipop][ion] > runpar->epsilon))
	{
          spkstr->sig_deltat_fire[ipop][ion] = sqrt(var_deltat_fire);
          if (spkstr->av_deltat_fire[ipop][ion] > runpar->epsilon)
	  {
            spkstr->cv_deltat_fire[ipop][ion] = 
              spkstr->sig_deltat_fire[ipop][ion] / 
              spkstr->av_deltat_fire[ipop][ion];
            spkstr->fr[ipop][ion] = 1000.0 / spkstr->av_deltat_fire[ipop][ion];
	  }
          else
	  {
            spkstr->fr[ipop][ion] = -999.7;
            spkstr->cv_deltat_fire[ipop][ion] = -999.7;
	  }
	}
        else
	{
          spkstr->fr[ipop][ion] = -999.8;
          spkstr->sig_deltat_fire[ipop][ion] = -999.8;
          spkstr->cv_deltat_fire[ipop][ion] = -999.8;
	}
      }
      else
      {
        spkstr->av_deltat_fire[ipop][ion] = -999.9;
        spkstr->sig_deltat_fire[ipop][ion] = -999.9;
        spkstr->cv_deltat_fire[ipop][ion] = -999.9;
      }
    }
  }

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)  
    {
      /* Computing the time-average Z1 */
      compute_Z_md_phi(&spkstr->Z1cos[ipop][ion], &spkstr->Z1sin[ipop][ion],
	spkstr->nfire[ipop][ion], &spkstr->Z1md[ipop][ion], 
        &spkstr->Z1phi[ipop][ion], runpar, fl);
    }
  }

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    Z1cos_av[ipop] = 0;
    Z1sin_av[ipop] = 0;
    for (ion=1; ion<=netpar->C[ipop].non; ion++)  
    {
      nfire_av[ipop] += spkstr->nfire[ipop][ion];
      Z1cos_av[ipop] += spkstr->Z1cos[ipop][ion];
      Z1sin_av[ipop] += spkstr->Z1sin[ipop][ion];

      if (!runpar->sm)
      {
        fprintf(fl.out, "ipop=%d ion=%d nfire_av=%lf Z1cos=%lf Z1cos_av=%lf "
        "Z1sin=%lf Z1sin_av=%lf\n", ipop, ion, nfire_av[ipop], 
        spkstr->Z1cos[ipop][ion], Z1cos_av[ipop], spkstr->Z1sin[ipop][ion],
        Z1sin_av[ipop]);
      }
    }
    fprintf(fl.out, "ipop=%d nfire_av=%lf\n", ipop, nfire_av[ipop]);
    compute_Z_md_phi(&Z1cos_av[ipop], &Z1sin_av[ipop], nfire_av[ipop],
      &Z1md_av[ipop], &Z1phi_av[ipop], runpar, fl);
    nfire_av[ipop] *= 1000.0 / (netpar->C[ipop].non * runpar->tstat);
  }

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    av->nfire_av[ipop] = nfire_av[ipop];
    av->Z1md_av[ipop] = Z1md_av[ipop];
    av->Z1phi_av[ipop] = Z1phi_av[ipop];
  }

  if ((!runpar->sm) && (0 > 1))
  {
    fprintf(fl.out, "\n");
    for (ipop=0; ipop<=netpar->npop; ipop++)
    {
      for (ion=1; ion<=netpar->C[ipop].non; ion++)  
      {
	fprintf(fl.out, "ipop=%d ion=%d nfire=%d ", ipop, ion,
	spkstr->nfire[ipop][ion]);
	fprintf(fl.out, "frinter=%lf cv=%lf ", spkstr->frinter[ipop][ion],
	spkstr->cv_deltat_fire[ipop][ion]);
        fprintf(fl.out, "Z1cos=%lf Z1sin=%lf Z1md=%lf Z1phi=%lf ",
        spkstr->Z1cos[ipop][ion], spkstr->Z1sin[ipop][ion],
	spkstr->Z1md[ipop][ion], spkstr->Z1phi[ipop][ion]);
        fprintf(fl.out, "spk_touch=%lf spk_before_touch=%lf\n",
	spkstr->spk_touch[ipop][ion], spkstr->spk_before_touch[ipop][ion]);
      }
      for (ion=1; ion<=netpar->C[ipop].non; ion++)  
      {
        fprintf(fl.zmp, "%d %d", ipop, ion);
        fprintf(fl.zmp, " %lf", 1000.0 * spkstr->nfire[ipop][ion] / 
          runpar->tstat);
        fprintf(fl.zmp, " %lf", spkstr->cv_deltat_fire[ipop][ion]);
        fprintf(fl.zmp, " %lf %lf\n", spkstr->Z1md[ipop][ion], 
          spkstr->Z1phi[ipop][ion]);
      }
    }

    fprintf(fl.out, "\n");
    for (ipop=0; ipop<=netpar->npop; ipop++)
    {

      fprintf(fl.out, "ipop=%d nfire_av=%lf Z1cos_av=%lf Z1sin_av=%lf "
	"Z1md_av=%lf Z1phi=%lf\n", ipop, nfire_av[ipop], Z1cos_av[ipop],
      Z1sin_av[ipop], Z1md_av[ipop], Z1phi_av[ipop]);
    }
  }

  if (netpar->inter.geom_dim <= 1)
  {
    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
      for (ion=1; ion<=netpar->C[ipop].non; ion++)  
      {
        xpos = 1.0 * (ion - 1) / netpar->C[ipop].rhd;
        fprintf(fl.fri, "%lf %lf %d %d ", xpos, spkstr->frinter[ipop][ion],
        ipop, ion);
        fprintf(fl.fri, "%lf %lf %lf ", netpar->C[ipop].Iextar[ion], 
	spkstr->spk_touch[ipop][ion], spkstr->spk_before_touch[ipop][ion]);
	fprintf(fl.fri, "%lf\n", spkstr->cv_deltat_fire[ipop][ion]);
      }
      if (ipop < netpar->npop) fprintf(fl.fri, "\n");
    }
  }
 
  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    non_with_cv = 0;
    spkstr->non_no_firing[ipop] = 0;
    spkstr->av_deltat_fire_pop[ipop] = 0.0;
    av_deltat_fire_pop_sq = 0.0;
    spkstr->fr_pop[ipop] = 0.0;
    fr_pop_sq[ipop] = 0.0;
    spkstr->cv_deltat_fire_pop[ipop] = 0.0;
    cv_deltat_fire_pop_sq = 0.0;

    if (ipop == 2)
    {
      spkstr->fr_subpop[ipop][1] = 0.0;
      spkstr->fr_subpop[ipop][2] = 0.0;
    }

    for (ion=1; ion<=netpar->C[ipop].non; ion++)  
    {
      spkstr->fr_pop[ipop] += spkstr->frinter[ipop][ion];
      fr_pop_sq[ipop] += spkstr->frinter[ipop][ion] *
	spkstr->frinter[ipop][ion];

      if (ipop == 2)
      {
	if (ion <= (int) (netpar->C[ipop].fracIext * netpar->C[ipop].non + 
	    1.0e5 * runpar->epsilon))
	{
          spkstr->fr_subpop[ipop][1] += spkstr->frinter[ipop][ion];
	}
        else
	{
          spkstr->fr_subpop[ipop][2] += spkstr->frinter[ipop][ion];
        }
      }

      if (!runpar->sm)
      {
        fprintf(fl.out, "ipop=%d ion=%d frinter=%lf fr_pop=%lf", ipop, ion, 
	spkstr->frinter[ipop][ion], spkstr->fr_pop[ipop]);
        if (ipop == 2)
	{
	  fprintf(fl.out, " fr_pop1=%lf fr_pop2=%lf\n", 
	    spkstr->fr_subpop[ipop][1], spkstr->fr_subpop[ipop][2]);
	}
        fprintf(fl.out, "\n");
      }

      spkstr->av_deltat_fire_pop[ipop] += spkstr->av_deltat_fire[ipop][ion];
      av_deltat_fire_pop_sq += pow(spkstr->av_deltat_fire[ipop][ion], 2.0);

     if ((spkstr->nfire[ipop][ion] >= 3) &&
          (spkstr->cv_deltat_fire[ipop][ion] > runpar->epsilon))
      {
        non_with_cv++;
        spkstr->cv_deltat_fire_pop[ipop] += spkstr->cv_deltat_fire[ipop][ion];
        cv_deltat_fire_pop_sq += pow(spkstr->cv_deltat_fire[ipop][ion], 2.0);
      }

      if (spkstr->nfire[ipop][ion] == 0)
      {
        spkstr->non_no_firing[ipop]++;
        /* printf("ipop=%d ion=%d non_no_firing=%d\n", ipop, ion,
	   spkstr->non_no_firing[ipop]); */
      }
    }
  
    spkstr->fr_pop[ipop] /= netpar->C[ipop].non;
    fr_pop_sq[ipop] /= netpar->C[ipop].non;
    fr_pop_diff = fr_pop_sq[ipop] - spkstr->fr_pop[ipop] * spkstr->fr_pop[ipop];
    if (fr_pop_diff >= 0.0)
    {
      spkstr->fr_pop_sd[ipop] = sqrt(fr_pop_diff);
    }
    else if (fr_pop_diff >= -runpar->epsilon)
    {
      spkstr->fr_pop_sd[ipop] = 0.0;
    }
    else
    {
      printf("ipop=%d fr_pop_diff=%lf < 0!\n", ipop, fr_pop_diff);
      spkstr->fr_pop_sd[ipop] = -999.9;
    }
    
    if ((netpar->pop_name[ipop] == 'P') && (ipop == 2))
    {
      fprintf(fl.out, "div: %lf %lf\n",  netpar->C[ipop].fracIext * 
	      netpar->C[ipop].non, (1.0 - netpar->C[ipop].fracIext) * 
	      netpar->C[ipop].non);
      if (netpar->C[ipop].fracIext * netpar->C[ipop].non > 
        1.0e5 * runpar->epsilon)
      {
        spkstr->fr_subpop[ipop][1] /= netpar->C[ipop].fracIext * 
          netpar->C[ipop].non;
      }
      else
      {
        spkstr->fr_subpop[ipop][1] = -999.9;
      }

      if ((1.0 - netpar->C[ipop].fracIext) * netpar->C[ipop].non > 
        1.0e5 * runpar->epsilon)
      {
        spkstr->fr_subpop[ipop][2] /= (1.0 - netpar->C[ipop].fracIext) * 
          netpar->C[ipop].non;
      }
      else
      {
        spkstr->fr_subpop[ipop][2] = -999.9;
      }
    }

    spkstr->av_deltat_fire_pop[ipop] /= netpar->C[ipop].non;
    av_deltat_fire_pop_sq /= netpar->C[ipop].non;
    spkstr->sig_av_deltat_fire_pop[ipop] = 
      compute_sd(spkstr->av_deltat_fire_pop[ipop], av_deltat_fire_pop_sq);

    if (non_with_cv >= 3)
    {
      spkstr->cv_deltat_fire_pop[ipop] /= non_with_cv;
      cv_deltat_fire_pop_sq /= non_with_cv;

      spkstr->sig_cv_deltat_fire_pop[ipop] = 
        compute_sd(spkstr->cv_deltat_fire_pop[ipop], cv_deltat_fire_pop_sq);
    }
    else
    {
      spkstr->sig_av_deltat_fire_pop[ipop] = -99.9;
      spkstr->sig_cv_deltat_fire_pop[ipop] = -99.9;
    }

    spkstr->non_with_cv[ipop] = non_with_cv;
  }

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    fprintf(fl.out, "ipop=%d av_fr=%lf sd_fr=%lf\n", ipop,
    spkstr->fr_pop[ipop], spkstr->fr_pop_sd[ipop]);
    fprintf(fl.out, "av_dtm=%lf sig_dtm=%lf av_cv=%lf sig_cv=%lf\n",
    spkstr->av_deltat_fire_pop[ipop], spkstr->sig_av_deltat_fire_pop[ipop], 
    spkstr->cv_deltat_fire_pop[ipop], spkstr->sig_cv_deltat_fire_pop[ipop]);
    fprintf(fl.out, "non_no_firing=%d non_with_cv=%d\n", 
    spkstr->non_no_firing[ipop], spkstr->non_with_cv[ipop]);

    av->fr_pop[ipop] = spkstr->fr_pop[ipop];
    av->fr_pop_sd[ipop] = spkstr->fr_pop_sd[ipop];
    av->av_cv[ipop]  = spkstr->cv_deltat_fire_pop[ipop];
    av->sig_cv[ipop] = spkstr->sig_cv_deltat_fire_pop[ipop];
    av->ratio_no_firing[ipop] = ((double) spkstr->non_no_firing[ipop]) /
    netpar->C[ipop].non;
  }

  for (ipop=0; ipop<=netpar->npop; ipop++)
  {
    if (netpar->pop_name[ipop] == 'P')
    {
      av->fr_subpop[ipop][1] = spkstr->fr_subpop[ipop][1];
      av->fr_subpop[ipop][2] = spkstr->fr_subpop[ipop][2];
      fprintf(fl.out, "ipop=%d av_subpop_fr=%lf %lf\n", ipop, 
        spkstr->fr_subpop[ipop][1], spkstr->fr_subpop[ipop][2]);
    }
  }
  
  free_dvector(nfire_av, 0, netpar->npop);
  free_dvector(Z1cos_av, 0, netpar->npop);
  free_dvector(Z1sin_av, 0, netpar->npop);
  free_dvector(Z1md_av , 0, netpar->npop);
  free_dvector(Z1phi_av, 0, netpar->npop);
  free_dvector(fr_pop_sq, 0, netpar->npop+1);

  for (ipop=0; ipop<=netpar->npop+1; ipop++)
  {
    if (ipop <= 1)
    {
      non_subpop = netpar->C[ipop].non;
    }
    else if (ipop == 2)
    {
      non_subpop = (int) (netpar->C[ipop].fracIext * netpar->C[ipop].non + 
        runpar->epsilon);
    }
    else if (ipop == 3)
    {
      non_subpop = (int) ((1.0 - netpar->C[ipop-1].fracIext) * 
        netpar->C[ipop-1].non + runpar->epsilon);
    }
    else
    {
      printf("wrong ipop=%d\n", ipop);
      exit(0);
    }

    fprintf(fl.out, "ipop=%d fracIext=%lf non_subpop=%d\n", ipop, 
      netpar->C[ipop].fracIext, non_subpop);

    if (non_subpop > 0)
    {
      x_factor = 1000.0 * runpar->nhist / (runpar->tstat * non_subpop);
    } 
   else
    {
      x_factor = 0.0;
    } 

    for (ihist=0; ihist<=runpar->nhist-1; ihist++)  
    {
      fr_x_hist = x_factor * spkstr->fr_hist[ipop][ihist];
      fprintf(fl.his, "%d %lf %d\n", ihist, fr_x_hist, ipop);
    }
    fprintf(fl.his, "  \n");
  }
}

/* This function updates the arranys of Vav every time step.      */
/* The arrays are needed to compute the synchrony measure chi.    */
void compute_voltage_statistics(V_aver *Vav, syn_str *synstr, net_par *netpar,
     run_par *runpar, avr_val *av, fl_st fl)
{
  double *sigma_V_sq, sigma_Vpop_sq, pop_av_sigma_V_sq;
  int ipop, ion, isynvar;

  /* Normalizing Vav */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    sigma_V_sq = dvector(1, netpar->C[ipop].non);

    Vav[ipop].Vpop_avt    /= runpar->tstat;
    Vav[ipop].Vpop_sq_avt /= runpar->tstat;
    sigma_Vpop_sq = Vav[ipop].Vpop_sq_avt - pow(Vav[ipop].Vpop_avt, 2.0);

    pop_av_sigma_V_sq = 0.0;
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      Vav[ipop].V_avt[ion]    /= runpar->tstat;
      Vav[ipop].V_sq_avt[ion] /= runpar->tstat;
      sigma_V_sq[ion] = Vav[ipop].V_sq_avt[ion] - pow(Vav[ipop].V_avt[ion], 
      2.0);
      pop_av_sigma_V_sq += sigma_V_sq[ion];
    }
    pop_av_sigma_V_sq /= netpar->C[ipop].non;

    if (pop_av_sigma_V_sq < 0.0)
      Vav[ipop].chi = -9999.0;
    else if (sigma_Vpop_sq < 0.0)
      Vav[ipop].chi = -9998.0;
    else
    {
      Vav[ipop].chi = sqrt(sigma_Vpop_sq / pop_av_sigma_V_sq);
    }

    fprintf(fl.out, "Vpop_avt=%lf Vpop_sq_avt=%lf sigma_Vpop_sq=%lf\n",
    Vav[ipop].Vpop_avt, Vav[ipop].Vpop_sq_avt, sigma_Vpop_sq);
    fprintf(fl.out, "pop_av_sigma_V_sq=%lf chi=%lf\n", pop_av_sigma_V_sq,
    Vav[ipop].chi);

    av->chi[ipop] = Vav[ipop].chi;

    free_dvector(sigma_V_sq, 1, netpar->C[ipop].non);
  }

  /* Normalizing synvar_avt */
  for (ipop=1; ipop<=netpar->npop; ipop++)
  {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
    {
      for (isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
      {
        synstr->synvar_avt[ipop][ion][isynvar] /= runpar->tstat;
        synstr->Isynvar_avt[ipop][ion][isynvar] /= runpar->tstat;
      }
    }
  }

  if ((!runpar->sm) && (1 > 0))
  {
    for (ipop=1; ipop<=netpar->npop; ipop++)
    {
    for (ion=1; ion<=netpar->C[ipop].non; ion++)
      {
        fprintf(fl.zmp, "%d %d", ipop, ion);
	fprintf(fl.zmp, " %lf", Vav[ipop].V_avt[ion]);
        for (isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
        {
          fprintf(fl.zmp, " %lf", synstr->synvar_avt[ipop][ion][isynvar]);
        }
        for (isynvar=1; isynvar<=synstr->nsynvar[ipop]; isynvar++)
        {
          fprintf(fl.zmp, " %lf", synstr->Isynvar_avt[ipop][ion][isynvar]);
        }
	fprintf(fl.zmp, "\n");
      }
    }
  }
}

/* This function compute the standard deviation.   */
double compute_sd(double av, double av_sq)
{
  double var, sd;
 
  var = av_sq - av * av;

  if (var > 0.0)
  {
    sd = sqrt(var);
  }
  else
  {
    sd = -998.9;
  }

  return(sd);
}

/* This function compute the function of tau0=Cm/gL and taus=tsyn that   */
/* yields the exremum PSP for an exponential PSC.                        */
double functau(double gL, double Cm, double tsyn)
{
  double tau0, taus, dft, rtt, tp, ftau;

  tau0 = Cm / gL;
  taus = tsyn;
  dft = tau0 - taus;
  rtt = taus / tau0;
  tp = -tau0 * taus * log(rtt) / dft;
  ftau = (tau0 / dft) * (pow(rtt, (taus / dft)) - pow(rtt, (tau0 / dft)));

  return ftau;
}

/* This function compute the maximim of the function of tau0=Cm/gL tausa  */
/* and tausb for extremem NMDA PSP.                                       */
double functau_NMDA(double gL, double Cm, double tausa, double tausb,
       run_par *runpar, fl_st fl)
{
  double tau0, taus, dft, rtt, tp, ftau;
  double coef0, coef1, coef2, tmax, tt, Vfun, Vfun_max;
  double t0, t1, t2, V0, V1, V2, xb, xc;
  int it, itmax, numt;

  numt = 1000;
  tau0 = Cm / gL;
  tmax = 10.0 * max(max(tau0, tausa), tausb);

  if (fabs(tau0 - tausa) < runpar->epsilon)
  {
    printf("tau0=%lf = tausa=%lf\n", tau0, tausa);
    exit(0);
  }
  else if (fabs(tau0 - tausb) < runpar->epsilon)
  {
    printf("tau0=%lf = tausb=%lf\n", tau0, tausb);
    exit(0);
  }
  else if (fabs(tausa - tausb) < runpar->epsilon)
  {
    printf("tausa=%lf = tausb=%lf\n", tausa, tausb);
    exit(0);
  }

  /*  	    
  fprintf(fl.out, "tau0=%lf tausa=%lf tausb=%lf\n", tau0, tausa, tausb);
  fprintf(fl.out, "numt=%d tmax=%lf tau0=%lf\n", numt, tmax, tau0);
  */
  coef0 = tau0 * tau0  / ((tausa - tau0) * (tausb - tau0));
  coef1 = tau0 * tausa / ((tausa - tau0) * (tausb - tausa));
  coef2 = tau0 * tausb / ((tausb - tau0) * (tausb - tausa));
  /*
  fprintf(fl.out, "coef0=%lf coef1=%lf coef2=%lf\n", coef0, coef1, coef2);
  */
  itmax = 0;
  Vfun_max = 0.0;

  for (it=1; it<=numt; it++)
  {
    tt = it * tmax / numt;
    Vfun = coef0 * exp(-tt / tau0) - coef1 * exp(-tt / tausa) + 
      coef2 * exp(-tt / tausb);
    if (Vfun > Vfun_max)
    {
      itmax = it;
      Vfun_max = Vfun;
    }
  }

  /*
  fprintf(fl.out, "itmax=%d ttmax=%lf Vfun_max=%lf\n", itmax,
  itmax * tmax / numt, Vfun_max);
  */
  if ((itmax >= 2) && (itmax < numt-1))
  {
    t0 = (itmax-1) * tmax / numt;
    t1 = itmax * tmax / numt;
    t2 = (itmax+1) * tmax / numt;
    V0 = coef0 * exp(-t0 / tau0) - coef1 * exp(-t0 / tausa) + 
      coef2 * exp(-t0 / tausb);
    V1 = coef0 * exp(-t1 / tau0) - coef1 * exp(-t1 / tausa) + 
      coef2 * exp(-t1 / tausb);
    V2 = coef0 * exp(-t2 / tau0) - coef1 * exp(-t2 / tausa) + 
      coef2 * exp(-t2 / tausb);
    /*
    fprintf(fl.out, "t0=%lf t1=%lf t2=%lf\n", t0, t1, t2);
    fprintf(fl.out, "V0=%lf V1=%lf V2=%lf\n", V0, V1, V2);
    */

    xb = V2 - V0;
    xc = V0 - 2.0 * V1 + V2;
    if (fabs(xc) < runpar->epsilon)
    {
      Vfun_max = V1;
    }
    else        
    {
       Vfun_max = V1 - 0.125 * xb * xb / xc;
    }
  }
  /*
  fprintf(fl.out, "itmax=%d ttmax=%lf Vfun_max=%lf\n", itmax,
  itmax * tmax / numt, Vfun_max);
  */
  ftau  = Vfun_max;
  return ftau;
}

/* This function computes Zmd and Zphi.     */
void compute_Z_md_phi(double *Zcos, double *Zsin, int nfire, double *Zmd, 
     double *Zphi, run_par *runpar, fl_st fl)
{
  if (nfire > 0)
  {
    *Zmd = 2.0 * sqrt(*Zcos * *Zcos + *Zsin * *Zsin) / (1.0 * nfire);

    if (*Zmd > runpar->epsilon)
    {
      *Zphi = atan2(*Zcos, *Zsin) / Pi; 
    }
    else
    {
      *Zphi = -999.9;
    }
  }
  else
  {
    *Zmd = 0.0;
    *Zphi = -999.8;
  }
}

/* This function generates Gaussian random numbers                        */
double gasdev(par_all *parall)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if (iset == 0) {
    do {
      v1=2.0*get_rn_dbl(parall->rng_ptr)-1.0;
      v2=2.0*get_rn_dbl(parall->rng_ptr)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq>=1.0 || rsq == 0.0);

    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  }
  else {
    iset=0;
    return gset;
  }
}


/* Linear interpolation */
double lininter(double x1, double x2, double xc, double y1, double y2)
{
  double linter;

  linter = ((xc-x1)*y2+(x2-xc)*y1) / (x2-x1) ;
  return(linter);
}

/* tt mod TT: result within [0..T)] */
double t_mod(double tt, double TT)
{
  double tcal, tmod;

  tcal = (tt >= 0) ? tt : tt + TT;
  
  if (tcal < 0.0)
  {
    printf("tcal=%lf < 0", tcal);
    exit(0);
  }

  tmod = fmod(tcal, TT);

  if (tmod > TT)
  {
    printf("tcal=%lf > TT=%lf", tcal, TT);
    exit(0);
  }
  
  return tmod;
}

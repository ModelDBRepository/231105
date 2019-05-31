#ifndef Mlinea
#define Mlinea 120
#endif

#ifndef MIapp
#define MIapp 10
#endif

typedef struct cell_par{
  double gNa, gNaP, gKdr, gKdr1, gKdr3, gA, gZ, gKz, gKd, gL, DelgL, *gLar;
  double Iext, DelIext, *Iextar, fracIext, gback_nor;
  double g_one_psp;
  double Cm, VNa, VK, VL, Vopto;
  double phi, gamma, two_pow_gamma;
  double thetam, sigmam, thetap, sigmap;
  double thetah, sigmah, thetan, sigman;
  double thetaa, sigmaa, thetab, sigmab, tauA, tauB;
  double thetaz, sigmaz, tauZ;
  double rhd, rho, Vinc1, Vinc2;
  double opto_amp, opto_sig, opto_freq, *opt_con, Iinject, tinject;
  /* double Iapp[MIapp], Iapp_a[MIapp];                */
  /* int n_Iapp, nin_Iapp_a, nt_Iapp[MIapp], i_Iapp;   */
  int nod, non, nceq, nseq, neq;
  int ion_inject;
  char model_type[3];
  char inject_current, gamma_for_syn;
} cell_par;

typedef struct tl_par{
  char nw;
  double Av, Bv, Tper, phi_read, Bp, *phi, tcfrac, tc, tauc;
  double Cvmin, Cvmax, *Cv, frac_only_p, frac_only_r; 
  double Avnw;
  double Tall, Tnw, Telev, Tw;
  double AvBcTo2p;
  int nspike_max;
  int ntl, itl, pr;
  char thal_input, determine_phi;
  char wd, ds_type;
} tl_par;

typedef struct syn_par_all{
  double rho_concur;
  int geom_dim;
  char con_shape, consider;
} syn_par_all;

typedef struct syn_receptor_par{
  double ths, sigs, tsynr, tsynd, thetanp, sigmanp;;
  double Vrev, DelVrev_O_DelIext;
} syn_receptor_par;

typedef struct conductance{
  double g, G, Vpsp;
} conductance;

typedef struct syn_coup_par{
  conductance AMPA, NMDA, GABAA;
  double UU, taur, tauf, xic;
  double tau_delay, Del_tau_delay;
  double Kin, CVKin, lam, rho_in;
  double gel, Gel, Kel, rho_el;
  double Anor;
  double **tdcoup, *nwprob;
  int *nwcoup, **wcoup;
  int *npostw, *npostp;
  int *nw_e_coup, **w_e_coup;
  int condition_for_nonzero_synapse_AMPA, condition_for_nonzero_synapse_NMDA;
  int condition_for_nonzero_synapse_GABAA;
} syn_coup_par;

typedef struct net_par{
  cell_par C[4];           /* cells: 1 - E, 2 - P, 3 - S */
  tl_par T;
  syn_par_all inter;       /* Synaptic interaction - general parameters */
  syn_receptor_par P, N, A;
  syn_coup_par S[4][4];    /* s[a][b]: b -> a connections */
                           /* b=0: T neurons              */
  double factETPT;

  double Length;
  double Volt_thresh;
  double noise, Kfactor, CVKin_fact;
  int npop;
  int jpop_start;   /* 0 - with thalamic input, 1 - without thalamic input. */
  char pop_name[10];
  char scalingc, scaleK;
  char process_time_delay;
} net_par;

typedef struct syn_par{
  double gsyn, *Vsyn, tsyn;
  int is_it_nmda;
} syn_par;

typedef struct syn_str{

  /* post: ipop */
  int *nsynvar, *nmda_on_population;  /* [ipop]               */
  int mc;
  char ***csynvar;                    /* [ipop][isynvar][0,1] */
  double gel_one_syn;

  /* post: ipop */
  syn_par **synpar;                   /* [ipop][isynvar] [ion] */
  double ***synvar;                   /* [ipop][ion][isynvar]  */
  double ***synvar_avt;               /* [ipop][ion][isynvar]  */
  double ***Isynvar_avt;              /* [ipop][ion][isynvar]  */

  /* post: ipop, ion; pre: jpop  */
  int ***isyn_to_send;                /* [ipop][jpop][2]      */
  int **nonzero_gsyn;                 /* [ipop][jpop]         */

  /* post: ipop */
  double **Isyn_cont;                 /* [ipop][ion]          */

  double ***xvar;                     /* [ipop][jpop][jon]    */

  double **Iel_cont;                  /* [ipop][ion]          */
  double **el_coef_sum, **el_coef_V;  /* [ipop][ion]          */
} syn_str;

typedef struct V_aver{
  double Vpop, chi;
  double Vpop_avt, Vpop_sq_avt, *V_avt, *V_sq_avt;
} V_aver;

typedef struct spike_in_deltat{
  double tspike;
  int jpop, jon;
} spike_in_deltat;

typedef struct spk_str{
  int mspk;

  /* double *tl_tspike; */                  /* [itl]                */

  double *tspike;
  int non_all, *jpop, *jon, nspike; 

  double ***t_p_delay_spike;               /* [jpop][jon][kspk] */
  double ****x_u_delay_spike;              /* [ipop][jpop][jon][kspk] */
  int  **spike_exist, *nspike_pop;
  double **tfire_prev;
  double **av_deltat_fire, **av_deltat_fire_sq, **fr, **frinter;
  double **sig_deltat_fire, **cv_deltat_fire;
  double *av_deltat_fire_pop, *sig_av_deltat_fire_pop, *fr_pop, *frinter_pop;
  double *fr_pop_sd;
  double **fr_subpop;
  double *cv_deltat_fire_pop, *sig_cv_deltat_fire_pop;
  double **Z1cos, **Z1sin, **Z1md, **Z1phi;
  double **spk_touch, **spk_before_touch;
  int **fr_hist;
  int **nfire;
  int *non_with_cv, *non_no_firing;

  double ****time_spike_delay;     /* [ipop][jpop][ion][ndelay] */
  int **ndelay, **iptr_delay;      /* [ipop][jpop] */
  double **spikes_in_deltat;       /* jpop, jon, tpeak */
  spike_in_deltat *sp_in_dt;
  int nsp_in_dt;

  V_aver *Vav;
} spk_str;

typedef struct run_par{
  
  double epsilon, *deltat;
  double time_all, tmcol, tstat, traster;
  double t_touch_interval;
  int ndeltat, ideltat, *nt, ntall;
  int *nwrite, **nwritear;
  int twrite;
  int sm, sp;
  int nhist;
  char method, incond, fpcal, smforce;
  char write_aux;
} run_par;

typedef struct func_cell_model
{
  void (*read_cell_par)(cell_par *cpar, fl_st fl);
  void (*write_cell_par)(cell_par *cpar, run_par runpar, fl_st fl);
  void (*steady_state_var)(cell_par *cpar, double *Varb, run_par *runpar, 
     fl_st fl);
  void (*update_cell)(double *Varc, double *kout, int ipop, int ion, 
     cell_par *cpar, syn_par_all *inter, syn_receptor_par *AMPA,
     double Iapp_now, double Isyn_cont, double Iel_cont, int it, fl_st fl);
} func_cell_model;

typedef struct transfer_to_func_solve{
  double *rand_num, *Bp;
  fl_st *fl;  
} transfer_to_func_solve;


/* Function Declaration */
void one_par(par_all *parall, avr_val *av, fl_st fl);
void read_npop(net_par *netpar, avr_val *av, fl_st fl);
func_cell_model *define_fcm_array(long nl, long nh);
void free_fcm_array(func_cell_model *v, long nl, long nh);
void free_nwritear_arrays(net_par *netpar, run_par *runpar, fl_st fl);
int **pivector(long nl, long nh);
void free_pivector(int **v, long nl, long nh);
double **pdvector(long nl, long nh);
void free_pdvector(double **v, long nl, long nh);
void read_input(func_cell_model *fcm, net_par *netpar, run_par *runpar,
     par_all *parall, fl_st fl);
void synaptic_parameters_read_write(int ipop, int jpop, net_par *netpar,
     fl_st fl);
void ab_Kin_read_write(char *ab, syn_coup_par *Sij, fl_st fl);
void AMPA_read_write(char *ab, syn_coup_par *Sij, cell_par *Cj, char scalingc,
     char scaleK, double Kfactor, double CVKin_fact, double factETPT, fl_st fl);
void NMDA_read_write(char *ab, syn_coup_par *Sij, cell_par *Cj, char scalingc,
     char scaleK, double Kfactor, double CVKin_fact, double factETPT, fl_st fl);
void GABAA_read_write(char *ab, syn_coup_par *Sij, cell_par *Cj, char scalingc,
     char scaleK, double Kfactor, double CVKin_fact, fl_st fl);
void utxt_read_write(char *ab, syn_coup_par *Sij, fl_st fl);
void define_variables_arrays(net_par *netpar, run_par *runpar,
     double ****Varbar, fl_st fl);
void initiate_synaptic_strengths_and_variables(syn_str *synstr, net_par *netpar,
     run_par *runpar, fl_st fl);
void one_synaptic_type_strengths_variables(int ipop, int jpop, char cell_ch,
     syn_receptor_par *R, conductance cnd, double Kin, double UU, cell_par *CC,
     char scalingc, int *first_input_to_ipop, int *isynvar, syn_str *synstr,
     run_par *runpar, fl_st fl);
void define_Vav(net_par *netpar, run_par *runpar, V_aver **Vav, fl_st fl);
void free_Vav(net_par *netpar, run_par *runpar, V_aver *Vav, fl_st fl);
void initiate_electrical_strengths(syn_str *synstr, net_par *netpar,
     run_par *runpar, fl_st fl);
void free_variables_arrays(net_par *netpar, run_par *runpar,
     double ***Varbar, fl_st fl);
void define_old_variables_arrays(net_par *netpar, run_par *runpar,
     double ****Varold, int ***after_max_vol, fl_st fl);
void free_old_variables_arrays(net_par *netpar, run_par *runpar,
     double ***Varold, int **after_max_vol, fl_st fl);
void compute_heterogeneous_intrinsic_parameters(net_par *netpar, 
     run_par *runpar, par_all *parall, fl_st fl);
void free_heterogeneous_intrinsic_parameters(net_par *netpar, 
     run_par *runpar, fl_st fl);
void compute_opto_conductance(net_par *netpar, run_par *runpar, fl_st fl);
void free_opto_arrays(net_par *netpar, run_par *runpar, fl_st fl);
void substitute_connectivity(net_par *netpar, run_par *runpar, par_all *parall,
     fl_st fl);
void compute_normalization_factor(syn_coup_par *scp, syn_par_all *inter, 
     int ipop, int jpop, int nod_pre, int non_pre, double rhd_pre, 
     double rho_pre, run_par *runpar, fl_st fl);
double numerical_normalization_factor_one_d(double Kin, double lam, 
       double rho_pre, int non_pre, char con_shape, run_par *runpar, fl_st fl);
double numerical_normalization_factor_two_d(double Kin, double lam, 
       double rhd_pre, int nod_pre, char con_shape, run_par *runpar, fl_st fl);
void find_coupling_matrix_zero_d(syn_coup_par *scp, cell_par *cpost,
     cell_par *cpre, int ipop, int jpop, run_par *runpar, par_all *parall,
     fl_st fl);
void compute_coupling_in_prob(syn_coup_par *scp, cell_par *cpost,
     cell_par *cpre, int ipop, int jpop, run_par *runpar, par_all *parall,
     fl_st fl);
void find_coupling_matrix_one_d(syn_coup_par *scp, char con_shape,
     int geom_dim, double Length, cell_par *cpost, cell_par *cpre, int ipop,
     int jpop, run_par *runpar, par_all *parall, fl_st fl);
void find_coupling_matrix_two_d(syn_coup_par *scp, char con_shape,
     double Length, cell_par *cpost, cell_par *cpre, int ipop, int jpop, 
     run_par *runpar, par_all *parall, fl_st fl);
void find_electrical_matrix_zero_d(syn_coup_par *scp, char con_shape,
     int geom_dim, double Length, cell_par *cpost, cell_par *cpre, int ipop,
     int jpop, run_par *runpar, par_all *parall, fl_st fl);
void free_connectivity_arrays(net_par *netpar, run_par *runpar, fl_st fl);
void compute_connectivity_statistics(int *ncoup, int non_pre,
     run_par *runpar, fl_st fl);
void compute_connectivity_d_statistics(double *ncoup, int non_post,
     run_par *runpar, fl_st fl);
void in_con(func_cell_model *fcm, double ***Varbar, net_par *netpar,
     run_par *runpar, par_all *parall, fl_st fl);
void in_con_one_pop(func_cell_model *fcm, double **Varbar, cell_par *cpar, 
     syn_par_all *inter, syn_receptor_par *AMPA, int ipop, run_par *runpar,
     par_all *parall, fl_st fl);
void define_synaptic_variables(syn_str *synstr, net_par *netpar, 
     run_par *runpar, fl_st fl);
void substitute_condition_check_syn(int ipop, int jpop, char cell_ch,
     int num_add, int *condition_for_nonzero_synapse, conductance cnd,
     char scalingc, syn_str *synstr, run_par *runpar, fl_st fl);
void update_check_nsynvar_mc(syn_str *synstr,  int ipop, int jpop, char label,
     fl_st fl);
void free_synaptic_variables(syn_str *synstr, net_par *netpar, 
     run_par *runpar, fl_st fl);
void define_spike_struct(spk_str *spkstr, net_par *netpar, run_par *runpar, 
     fl_st fl);
void free_spike_struct(spk_str *spkstr, net_par *netpar, run_par *runpar, 
     fl_st fl);
void define_storage_spikes_td(net_par *netpar, run_par *runpar,
     double *****time_spike_delay, int ***ndelay, int ***iptr_delay,
     spike_in_deltat **sp_in_dt, int non_all, fl_st fl);
void free_storage_spikes_td(net_par *netpar, run_par *runpar,
     double ****time_spike_delay, int **ndelay, int **iptr_delay,
     spike_in_deltat *sp_in_dt, int non_all, fl_st fl);
void define_pop_non(net_par *netpar, run_par *runpar, double ***tsp, fl_st fl);
void set_val_pop_non(double val, net_par *netpar, run_par *runpar,
     double **jagar, fl_st fl);
void free_pop_non(net_par *netpar, run_par *runpar, double **tsp, fl_st fl);
void define_pop_non_nsome(net_par *netpar, run_par *runpar, double ****tsp, 
     int nsome, fl_st fl);
void free_pop_non_nsome(net_par *netpar, run_par *runpar, double ***tsp,
     int nsome, fl_st fl);
void define_pop_pop_non_spk(net_par *netpar, run_par *runpar, double *****tsp, 
     int mspk, fl_st fl);
void free_pop_pop_non_spk(net_par *netpar, run_par *runpar, double ****tsp, 
     int mspk, fl_st fl);
void define_int_pop_non(net_par *netpar, run_par *runpar, int ***nsp,
     fl_st fl);
void free_int_pop_non(net_par *netpar, run_par *runpar, int **nsp, fl_st fl);
void determine_thal_Cv(net_par *netpar, run_par *runpar, par_all *parall,
     fl_st fl);
void initialize_thalamic_variables_and_spikes(spk_str *spkstr, net_par *netpar,
     run_par *runpar, par_all *parall, fl_st fl);
double solve_phi_p_sin_phi_eq(double rand_num, double Bp, run_par *runpar,
       fl_st fl);
double Bp_phi_func(double phi, void *ptr);
void n_run(func_cell_model *fcm, double ***Varbar, syn_str *synstr, 
     spk_str *spkstr, net_par *netpar, run_par *runpar, par_all *parall,
     avr_val *av, fl_st fl);
void one_integration_step(func_cell_model *fcm, net_par *netpar,
     run_par *runpar, syn_str *synstr, double ***Varbar, double ***kin,
     double ***kout, double ***Varc, double delt, double time, int it,
     fl_st fl);
void compute_total_synaptic_conductance_on_a_neuron(double ***Varbar, 
     syn_str *synstr, double time, int it, double deltat, net_par *netpar,
     run_par *runpar, fl_st fl);
void pr_fct(double ***Varbar, syn_str *synstr, spk_str *spkstr, net_par *netpar,
     run_par *runpar, double time, int it, fl_st fl);
void spike_detect(double ***Varbar, double ***Varold, int **after_max_vol,
     int it, double time, syn_str *synstr, spk_str *spkstr, net_par *netpar,
     run_par *runpar, avr_val *av, fl_st fl);
int spike_detect_peak(double V0, double V1, double V2, double time,
    double *tpeak, double *Vpeak, net_par *netpar, run_par *runpar, fl_st fl);
int spike_detect_threshold(double V0, double V1, double time, double *tpeak, 
    double *Vpeak, net_par *netpar, run_par *runpar, fl_st fl);
void single_td_storing_spikes(int jpop, int jon, double tpeak, int it,
     double time, spk_str *spkstr, net_par *netpar, run_par *runpar, fl_st fl);
void multiple_td_storing_spikes(int jpop, int jon, double tpeak, int it,
     double time, spk_str *spkstr, net_par *netpar, run_par *runpar, fl_st fl);
void update_thalamic_variables(int it, double time, double deltat, 
     syn_str *synstr, spk_str *spkstr, net_par *netpar, run_par *runpar,
     par_all *parall, avr_val *av, fl_st fl);
void multiple_store_spikes_plus_td(int it, double time, double deltat, 
     syn_str *synstr, spk_str *spkstr, net_par *netpar, run_par *runpar,
     fl_st fl);
void short_term_plasticity(int it, double time, double ***xvar, double tprev,
     double tnow, int jpop, int jon, net_par *netpar, run_par *runpar, 
     fl_st fl);
void update_delayed_cortical_spikes(int it, double time, double deltat, 
     spk_str *spkstr, net_par *netpar, run_par *runpar, avr_val *av, fl_st fl);
void decay_post_synaptic_variables(syn_str *synstr, double time, int it,
     net_par *netpar, run_par *runpar, fl_st fl);
void update_post_synaptic_variables_for_pre_synaptic_spikes_a(syn_str *synstr, 
     spk_str *spkstr, double time, int it, double deltat, net_par *netpar, 
     run_par *runpar, fl_st fl);
void update_post_synaptic_variables_for_pre_synaptic_spikes_s(syn_str *synstr, 
     spk_str *spkstr, double time, int it, double deltat, net_par *netpar, 
     run_par *runpar, fl_st fl);
void update_Vav_arrays(double ***Varbar, V_aver *Vav, syn_str *synstr, int it,
     double time, double deltat, net_par *netpar, run_par *runpar, fl_st fl);
void compute_spike_statistics(spk_str *spkstr, net_par *netpar, run_par *runpar,
     avr_val *av, fl_st fl);
void compute_voltage_statistics(V_aver *Vav, syn_str *synstr, net_par *netpar,
     run_par *runpar, avr_val *av, fl_st fl);
double compute_sd(double av, double av_sq);
double functau(double gL, double Cm, double tsyn);
double functau_NMDA(double gL, double Cm, double tausa, double tausb,
       run_par *runpar, fl_st fl);
void compute_Z_md_phi(double *Zcos, double *Zsin, int nfire, double *Zmd, 
     double *Zphi, run_par *runpar, fl_st fl);
double gasdev(par_all *parall);
double lininter(double x1, double x2, double xc, double y1, double y2);
double t_mod(double tt, double TT);

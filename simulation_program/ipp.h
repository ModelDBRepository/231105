/* Structure Declaration */

typedef struct run_par{
  double epsilon;
} run_par;

typedef struct tl_par{
  char nw;
  double Av, Bv, Cv, Tper, phi_read, *phi, tcfrac, tc, tauc;
  double Avnw;
  double Tall, Tnw, Tw;
  double AvBcTo2p;
  int nspike_max;
  int ntl, itl;
  char thal_input, determine_phi, pr;
} tl_par;


/* Function Declaration */
void read_thalamic_par(tl_par *tlpar, run_par *runpar, fl_st fl);
void compute_spike_train(tl_par *tlpar, run_par *runpar, int **rng_ptr,
     fl_st fl);


#ifndef min
#define min(x,y) ( ((x) < (y)) ? (x) : (y) )
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

#ifndef Pi
#define Pi 3.1415926535897931
#endif

#ifndef NR_END
#define NR_END 1
#endif

#ifndef FREE_ARG
#define FREE_ARG char*
#endif

#ifndef div_nz
#define div_nz(x,y) ( (fabs(y) < (runpar.epsilon)) ? (0) : (x/y) )
#endif

#ifndef Gammaf
#define Gammaf(VV, theta,sigma) ( 1.0/(1.0+exp(-(VV-(theta))/(sigma))) )
#endif

#ifndef LinGam
#define LinGam(VV, theta,sigma) ( (VV-(theta)) /(1.0+exp(-(VV-(theta))/(sigma))) )
#endif

#ifndef Meq
#define Meq 101
#endif


/* Structure Declaration */

typedef struct fl_st{
  FILE *in, *tmp, *avr, *out, *col, *ras, *fri, *zmp, *his;
  FILE *spk;
} fl_st;

typedef struct par_all{
  int *rng_ptr;
  int sm;
  char scan_type;
} par_all;

typedef struct avr_val{
  double par;
  int ipar, irepeat;
  double fr_pop[21], fr_pop_sd[21], av_cv[21], sig_cv[21];
  double nfire_av[21], Z1md_av[21], Z1phi_av[21];
  double chi[21];
  double fr_subpop[21][3];
  double ratio_no_firing[21];
  int npop;
  char pop_name[10];
}  avr_val;


/* Function Declaration */

/* lfibrng6a.c */
int **init_rng_s_dbl(int ngen, int length, int seed);
double get_rn_dbl(int *genptr);

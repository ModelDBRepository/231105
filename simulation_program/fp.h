typedef struct transfer_to_func_fp{
  func_cell_model *fcm;
  cell_par *cpar;
  syn_par_all *inter;
  syn_receptor_par *AMPA;
  run_par *runpar;
  fl_st *fl;
  int ipop, ion;
} transfer_to_func_fp;


/* Function Declaration */
void find_fixed_point(func_cell_model *fcm, cell_par *cpar, run_par *runpar,
     syn_par_all *inter, syn_receptor_par *AMPA, double *Varbar, int ipop,
     int ion, fl_st fl);
void usrfun(double *xvec, int nvec, double *fvec, double **fjac, void *ptr);
void vecfunc(int nvec, double *xvec, double *fvec, void *ptr);
/*
int stability_fixed_point(int nvec, double *xvec, void *ptr, fl_st fl);
*/

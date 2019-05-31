void substitute_function_WB(func_cell_model *fcm, fl_st fl);
void read_cell_par_WB(cell_par *cpar, fl_st fl);
void steady_state_var_WB(cell_par *cpar, double *Varb, run_par *runpar, 
     fl_st fl);
void update_cell_WB(double *Varc, double *kout, int ipop, int ion, 
     cell_par *cpar, syn_par_all *inter, syn_receptor_par *AMPA,
     double Iapp_now, double Isyn_cont, double Iel_cont, int it, fl_st fl);

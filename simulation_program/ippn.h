
/* Structure Declaration */

/*
typedef struct time_div{
  double t_in_cycle, tc_updated;
  int past_n_cycles;
} time_div;
*/

typedef struct transfer_to_func{
  tl_par *tlpar;
  double *t_in_cycle, *tc_updated, *cnc, *log_rand_num;
  fl_st *fl;
} transfer_to_func;


/* Function Declaration */

double find_time_next_spike_w(double time, int itl, double rand_num,
       tl_par *tlpar, run_par *runpar, fl_st fl);
double find_time_next_spike_n(double time, int itl, double rand_num,
       tl_par *tlpar, run_par *runpar, fl_st fl);
double find_time_next_spike_l(double time, int itl, double rand_num,
       tl_par *tlpar, run_par *runpar, fl_st fl);
int find_tns_or_update_in_this_episode_n(double *t_next_spike,
    double *log_rand_num, tl_par *tlpar, run_par *runpar, fl_st fl);
int find_tns_or_update_in_this_episode_l(double *t_next_spike,
    double *log_rand_num, tl_par *tlpar, run_par *runpar, fl_st fl);
int find_tns_or_update_lrn_nw(double *t_next_spike, double *log_rand_num, 
    double t_in_episode, double Tend, double AA, int pr, run_par *runpar, 
    fl_st fl);
int find_tns_or_update_lrn_w(double *t_next_spike, double *log_rand_num, 
    double t_in_episode, tl_par *tlpar, run_par *runpar, fl_st fl);
int find_tns_or_update_lrn_elev(double *t_next_spike, double *log_rand_num, 
    double t_in_episode, tl_par *tlpar, run_par *runpar, fl_st fl);
void determine_contact_next_change(double t_in_episode_run, double *cnc,
     double *tc_updated, tl_par *tlpar, run_par *runpar, fl_st fl);
double compute_LQ(double t_in_cycle, double tc_updated, double cnc, 
       tl_par *tlpar, fl_st fl);
double find_time_LQ(double t_in_cycle, double tc_updated, double cnc,
        double log_rand_num, tl_par *tlpar, fl_st fl);
double LQ_func(double tLQ, void *ptr);

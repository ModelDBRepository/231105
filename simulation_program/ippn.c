/* This file generates spike trains according to an inhomogeneous Poisson */
/* process.                                                                  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "tca.h"
#include "tcn.h"
#include "ippn.h"

/* This function computes the time of the next spike. */
/* The animal is constantly whising.                  */
double find_time_next_spike_w(double time, int itl, double rand_num,
       tl_par *tlpar, run_par *runpar, fl_st fl)
{
  double log_rand_num, tc_updated;
  double t_next_spike, t_from_last_contact_start;
  double delta_t_next_spike_new, t_next_spike_new;
  double cnc;
  double LQnc, LQc;
  double Tper_before_spike;
  int nTper_before_spike, tspike_range_located;

  tlpar->itl = itl;

  t_next_spike = time;
  log_rand_num = log(rand_num);
  tspike_range_located = 0;

  if (tlpar->pr)
  {
   fprintf(fl.out, "\nbeginning: t_next_spike=%lf\n", t_next_spike);
   fprintf(fl.out, " log_rand_num=%lf tpdiv=%lf\n",log_rand_num, 
   time - log_rand_num/0.02);
  }

  t_from_last_contact_start = t_next_spike - 
    ((int) (t_next_spike /  tlpar->Tper)) * tlpar->Tper - tlpar->tc;
  if (t_from_last_contact_start < 0) t_from_last_contact_start += tlpar->Tper;

  if (tlpar->pr)
    fprintf(fl.out, "t_from_last_contact_start=%lf\n", 
    t_from_last_contact_start);

  /* in_contact, cnc: 0 - nc, 1 - c */

  if (t_from_last_contact_start < tlpar->tauc)
  /* time is within a contact period */
  {
    if (tlpar->pr)
      fprintf(fl.out, "time is within a contact period\n");
   
    nTper_before_spike = (int) (t_next_spike / tlpar->Tper);

    tc_updated = t_next_spike + fmod( (nTper_before_spike + 1) * tlpar->Tper + 
                       tlpar->tc + tlpar->tauc - t_next_spike, tlpar->Tper);
		   
    if (tlpar->pr) fprintf(fl.out, "a tc_updated=%lf\n", tc_updated);

    cnc = 1.0;

    LQc = compute_LQ(t_next_spike, tc_updated, cnc, tlpar, fl);
    
    if (tlpar->pr)
     fprintf(fl.out, "t_next_spike=%lf tc_updated=%lf cnc=%lf LQc=%lf\n",
     t_next_spike, tc_updated, cnc, LQc);

    if (LQc < log_rand_num)
    {
      tspike_range_located = 1;
    }
    else
    {
      t_next_spike = tc_updated;
      tc_updated += tlpar->Tper - tlpar->tauc;
      log_rand_num -= LQc;
      cnc = 0.0;

      if (tlpar->pr)
      {
        fprintf(fl.out, "\nno spike during contact:\n");
        fprintf(fl.out, "t_next_spike=%lf tc_updated=%lf log_rand_num=%lf "
        "cnc=%lf\n", t_next_spike, tc_updated, log_rand_num, cnc);
      }
    }

    if (tlpar->pr)
      fprintf(fl.out, "trl=%d t_next_spike=%lf\n", tspike_range_located,
      t_next_spike);
  }
  else
  {
    /* time is not within a contact period */
    if (tlpar->pr)
      fprintf(fl.out, "time is not within a contact period\n");
 
    tc_updated = ((int) (t_next_spike / tlpar->Tper)) *  tlpar->Tper + 
    tlpar->tc;
    if (tc_updated < t_next_spike) tc_updated += tlpar->Tper;
    cnc = 0.0;
 
    if (tlpar->pr)
     fprintf(fl.out, "tc_updated=%lf cnc=%lf\n", tc_updated, cnc);
  }

  if (!tspike_range_located)
  /* time is not within a contact period, or next spike was not located */
  /* in the contact period.                                             */
  {
    if (tlpar->pr)
      fprintf(fl.out, "time is not within a contact period!\n");

    LQnc = compute_LQ(t_next_spike, tc_updated, cnc, tlpar, fl);

    if (tlpar->pr)
      fprintf(fl.out, "t_next_spike=%lf tc_updated=%lf LQnc=%lf\n", 
      t_next_spike, tc_updated, LQnc);

    if (LQnc < log_rand_num)
    {
      tspike_range_located = 1;

      if (tlpar->pr)
        fprintf(fl.out, "tspike_range_located=%d\n", tspike_range_located);
    }
    else
    {
      log_rand_num -= LQnc;
      t_next_spike = tc_updated;
      cnc = 1.0;

      if (tlpar->pr)
        fprintf(fl.out, "log_rand_num=%lf t_next_spike=%lf cnc=%lf "
        "tspike_range_located=%d\n", log_rand_num, t_next_spike,
        cnc, tspike_range_located);
    }
  }

  if (!tspike_range_located)
  /* No next spike before the next contact */
  {
    LQc = compute_LQ(t_next_spike, t_next_spike + tlpar->tauc, 
    1.0, tlpar, fl);
    LQnc = compute_LQ(t_next_spike + tlpar->tauc, t_next_spike + tlpar->Tper,
    0.0, tlpar, fl);
    nTper_before_spike = (int) (log_rand_num / (LQc + LQnc));

    t_next_spike += nTper_before_spike *  tlpar->Tper;
    log_rand_num -= nTper_before_spike * (LQc + LQnc);

    if (tlpar->pr)
    {
      fprintf(fl.out, "\nno spike during %d cycles: t_next_spike=%lf "
      "log_rand_num=%lf\n",  nTper_before_spike, t_next_spike, log_rand_num);
      fprintf(fl.out, "log_rand_num=%lf LQC=%lf LQnc=%lf sum=%lf "
      " nTper_before_spike=%d\n", log_rand_num, LQc, LQnc, LQc+LQnc,
      nTper_before_spike);
    }

    if (LQc < log_rand_num)
    {
      tc_updated = t_next_spike + tlpar->tauc;
      cnc = 1.0;
    }
    else
    {
      t_next_spike += tlpar->tauc;
      log_rand_num -= LQc;
      tc_updated = t_next_spike + (tlpar->Tper - tlpar->tauc);
      cnc = 0.0;
    }
  }

  if (tlpar->pr)
    fprintf(fl.out, "t_next_spike=%lf log_rand_num=%lf tc_updated=%lf cnc=%lf"
    "\n", t_next_spike, log_rand_num, tc_updated, cnc);

  nTper_before_spike = (int) (t_next_spike / tlpar->Tper);
  Tper_before_spike = nTper_before_spike * tlpar->Tper;

  delta_t_next_spike_new = find_time_LQ(t_next_spike - Tper_before_spike, 
  tc_updated - Tper_before_spike, cnc, log_rand_num, tlpar, fl);
  t_next_spike_new = Tper_before_spike + delta_t_next_spike_new;

  if (fabs(t_next_spike_new - t_next_spike) < runpar->epsilon)
  {
    fprintf(fl.out, "t_next_spike=%lf = t_next_spike_new=%lf\n", t_next_spike,
    t_next_spike_new);
    fprintf(fl.out, "cnc=%lf\n", cnc);
  }

  t_next_spike = t_next_spike_new;

  if (tlpar->pr)
    fprintf(fl.out, "t_next_spike=%lf\n", t_next_spike);

  return(t_next_spike);
}

/* This function computes the time of the next spike.                     */
/* The animal is switching from non-whisking to whisking with a specific  */
/* phase.                                                                 */
double find_time_next_spike_n(double time, int itl, double rand_num,
       tl_par *tlpar, run_par *runpar, fl_st fl)
{
  double log_rand_num, t_next_spike, LQn, LQw;
  int find_tns, n_episode_no_spikes;

  tlpar->itl = itl;

  t_next_spike = time;
  log_rand_num = log(rand_num);

  if (tlpar->pr)
  {
    fprintf(fl.out, "in find_time_next_spike_n\n");
    fprintf(fl.out, "t_next_spike=%lf log_rand_num=%lf\n", t_next_spike, 
    log_rand_num);
  }

  if (find_tns = find_tns_or_update_in_this_episode_n(&t_next_spike, 
    &log_rand_num, tlpar, runpar, fl))
  {
    if (tlpar->pr) fprintf(fl.out, "in one episode find_tns=%d t_next_spike="
      "%lf\n", find_tns, t_next_spike);
    return(t_next_spike);
  }

  LQn = tlpar->Avnw * tlpar->Tnw;
  LQw = (tlpar->Av + tlpar->Cv[tlpar->itl] * (1000.0 / tlpar->Tper)) *
  tlpar->Tw;

  if (fabs(LQn + LQw) < runpar->epsilon)
  {
    printf("LQn+LQw=%lf\n", LQn + LQw);
    exit(0);
  }

  if (tlpar->pr) 
  {
    fprintf(fl.out, "\nLQn=%lf LQw=%lf\n", LQn, LQw);
    fprintf(fl.out, "t_next_spike=%lf log_rand_num=%lf\n", t_next_spike, 
    log_rand_num);
  }

  n_episode_no_spikes = (int) (-log_rand_num / (LQn + LQw));
  log_rand_num += n_episode_no_spikes * (LQn + LQw);
  t_next_spike += n_episode_no_spikes * (tlpar->Tnw + tlpar->Tw);

  if (tlpar->pr) 
    fprintf(fl.out, "n_episode_no_spikes=%d log_rand_num=%lf t_next_spike="
    "%lf\n\n", n_episode_no_spikes, log_rand_num, t_next_spike);

  if (find_tns = find_tns_or_update_in_this_episode_n(&t_next_spike, 
    &log_rand_num, tlpar, runpar, fl))
  {
    if (tlpar->pr) fprintf(fl.out, "in one episode find_tns=%d t_next_spike="
      "%lf\n", find_tns, t_next_spike);
    return(t_next_spike);
  }
  else
  {
    printf("no spike is found in the second call to "
    "find_tns_or_update_in_this_episode_n!\n");
    exit(0);
  }
}

/* This function computes the time of the next spike.                     */
/* The animal is switching linearly from non-whisking to whisking         */
double find_time_next_spike_l(double time, int itl, double rand_num,
       tl_par *tlpar, run_par *runpar, fl_st fl)
{
  double log_rand_num, t_next_spike, LQn, LQe, LQw;
  int find_tns, n_episode_no_spikes;

  tlpar->itl = itl;

  t_next_spike = time;
  log_rand_num = log(rand_num);

  if (tlpar->pr)
  {
    fprintf(fl.out, "in find_time_next_spike_l\n");
    fprintf(fl.out, "t_next_spike=%lf log_rand_num=%lf\n", t_next_spike, 
    log_rand_num);
  }

  if (find_tns = find_tns_or_update_in_this_episode_l(&t_next_spike, 
    &log_rand_num, tlpar, runpar, fl))
  {
    if (tlpar->pr) fprintf(fl.out, "in one episode find_tns=%d t_next_spike="
      "%lf\n", find_tns, t_next_spike);
    return(t_next_spike);
  }

  LQn = tlpar->Avnw * tlpar->Tnw;
  LQe = 0.5 * (tlpar->Avnw + tlpar->Av) * tlpar->Telev;
  LQw = (tlpar->Av + tlpar->Cv[tlpar->itl] * (1000.0 / tlpar->Tper)) *
  tlpar->Tw;

  if (fabs(LQn + LQe + LQw) < runpar->epsilon)
  {
    printf("LQn+LQw=%lf\n", LQn + LQe + LQw);
    exit(0);
  }

  if (tlpar->pr) 
  {
    fprintf(fl.out, "\nLQn=%lf LQe=%lf LQw=%lf\n", LQn, LQe, LQw);
    fprintf(fl.out, "t_next_spike=%lf log_rand_num=%lf\n", t_next_spike, 
    log_rand_num);
  }

  n_episode_no_spikes = (int) (-log_rand_num / (LQn + LQe + LQw));
  log_rand_num += n_episode_no_spikes * (LQn + LQe + LQw);
  t_next_spike += n_episode_no_spikes * (tlpar->Tnw +tlpar->Telev + tlpar->Tw); 

  if (tlpar->pr) 
    fprintf(fl.out, "n_episode_no_spikes=%d log_rand_num=%lf t_next_spike="
    "%lf\n\n", n_episode_no_spikes, log_rand_num, t_next_spike);

  if (find_tns = find_tns_or_update_in_this_episode_l(&t_next_spike, 
    &log_rand_num, tlpar, runpar, fl))
  {
    if (tlpar->pr) fprintf(fl.out, "in one episode find_tns=%d t_next_spike="
      "%lf\n", find_tns, t_next_spike);
    return(t_next_spike);
  }
  else
  {
    printf("no spike is found in the second call to "
    "find_tns_or_update_in_this_episode!\n");
    exit(0);
  }
}

/* This function finds the time of the next spike or increases log_rand_num */
/* within the present episode.                                              */ 
/* The animal is switching from non-whisking to whisking with a specific    */
/* phase.                                                                   */
int find_tns_or_update_in_this_episode_n(double *t_next_spike,
    double *log_rand_num, tl_par *tlpar, run_par *runpar, fl_st fl)
{
  double t_in_episode;
  int find_tns;

  t_in_episode = fmod(*t_next_spike + runpar->epsilon, tlpar->Tnw + tlpar->Tw) -
  runpar->epsilon;

  if (t_in_episode - tlpar->Tnw < 0.0)
  {
    if (tlpar->pr) 
      fprintf(fl.out, "call nw t_next_spike=%lf t_in_episode=%lf\n", 
      *t_next_spike, t_in_episode);
    if (find_tns = find_tns_or_update_lrn_nw(t_next_spike, log_rand_num,
	t_in_episode, tlpar->Tnw, tlpar->Avnw, tlpar->pr, runpar, fl))
    {
      if (tlpar->pr) fprintf(fl.out, "find_tns=%d t_next_spike=%lf\n",
        find_tns, *t_next_spike);
      return(find_tns);
    }
  }

  if (tlpar->pr) 
    fprintf(fl.out, "end_of_nw t_in_episode=%lf t_in_episode-Tnw=%e\n",
    t_in_episode, t_in_episode-tlpar->Tnw);

  t_in_episode = fmod(*t_next_spike, tlpar->Tnw + tlpar->Tw);

  if (tlpar->pr) 
    fprintf(fl.out, "before_of_w t_in_episode=%lf t_in_episode-Tnw=%e\n",
    t_in_episode, t_in_episode-tlpar->Tnw);

  if (t_in_episode - tlpar->Tnw > -runpar->epsilon)
  {
    if (tlpar->pr)
      fprintf(fl.out, "call w t_next_spike=%lf t_in_episode=%lf\n",
      *t_next_spike, t_in_episode);
    if (find_tns = find_tns_or_update_lrn_w(t_next_spike, log_rand_num,
       t_in_episode, tlpar, runpar, fl))
    {
     if (tlpar->pr) fprintf(fl.out, "find_tns=%d\n", find_tns);
     return(find_tns);
    }
  }

  t_in_episode = fmod(*t_next_spike, tlpar->Tnw + tlpar->Tw);

  if (tlpar->pr) fprintf(fl.out, "end_of_episode t_in_episode=%lf \n",
  t_in_episode);

  return(find_tns);
 }

/* This function finds the time of the next spike or increases log_rand_num */
/* within the present episode.                                              */ 
/* The animal is switching linearly from non-whisking to whisking         */
int find_tns_or_update_in_this_episode_l(double *t_next_spike,
    double *log_rand_num, tl_par *tlpar, run_par *runpar, fl_st fl)
{
  double t_in_episode;
  int find_tns;

  if (tlpar->Telev < runpar->epsilon)
  {
    printf("Telev=%lf should be positive!\n", tlpar->Telev);
    exit(0);
  }

  /* During non-whisking */
  t_in_episode = fmod(*t_next_spike + runpar->epsilon, tlpar->Tnw + 
    tlpar->Telev + tlpar->Tw) - runpar->epsilon;

  if (t_in_episode - tlpar->Tnw < 0.0)
  {
    if (tlpar->pr) 
      fprintf(fl.out, "call nw t_next_spike=%lf t_in_episode=%lf\n", 
      *t_next_spike, t_in_episode);
    if (find_tns = find_tns_or_update_lrn_nw(t_next_spike, log_rand_num,
      t_in_episode, tlpar->Tnw, tlpar->Avnw, tlpar->pr, runpar, fl))
    {
      if (tlpar->pr) fprintf(fl.out, "find_tns=%d t_next_spike=%lf\n",
        find_tns, *t_next_spike);
      return(find_tns);
    }
  }

  if (tlpar->pr) 
    fprintf(fl.out, "end_of_nw t_in_episode=%lf t_in_episode-Tnw=%e\n",
    t_in_episode, t_in_episode-tlpar->Tnw);

  /* During elevation of firing rate */
  t_in_episode = fmod(*t_next_spike + runpar->epsilon, tlpar->Tnw + 
    tlpar->Telev + tlpar->Tw) - runpar->epsilon;

  if (tlpar->pr) 
    fprintf(fl.out, "in elev t_in_episode=%lf t_in_episode-Tnw=%e\n",
    t_in_episode, t_in_episode-tlpar->Tnw);

  if ((t_in_episode - tlpar->Tnw > -runpar->epsilon) &&
      (t_in_episode - tlpar->Tnw - tlpar->Telev < 0.0))
  {
    if (tlpar->pr)
      fprintf(fl.out, "call e t_next_spike=%lf t_in_episode=%lf\n",
      *t_next_spike, t_in_episode);
    if (find_tns = find_tns_or_update_lrn_elev(t_next_spike, log_rand_num,
	t_in_episode, tlpar, runpar, fl))
    {
     if (tlpar->pr) fprintf(fl.out, "find_tns=%d\n", find_tns);
     return(find_tns);
    }
  }

  if (tlpar->pr) 
    fprintf(fl.out, "end_of_elev t_in_episode=%lf t_in_episode-Tnw=%e\n",
    t_in_episode, t_in_episode-tlpar->Tnw);

  /* During whisking, constant firing rate  */
  t_in_episode = fmod(*t_next_spike + runpar->epsilon, tlpar->Tnw + 
    tlpar->Telev + tlpar->Tw) - runpar->epsilon;

  if (tlpar->pr) 
    fprintf(fl.out, "in_w t_in_episode=%lf t_in_episode-Tnw=%e\n",
    t_in_episode, t_in_episode-tlpar->Tnw);

  if (t_in_episode - tlpar->Tnw - tlpar->Telev >  -runpar->epsilon)
  {
    if (tlpar->pr)
      fprintf(fl.out, "call w t_next_spike=%lf t_in_episode=%lf\n",
      *t_next_spike, t_in_episode);
    if (find_tns = find_tns_or_update_lrn_nw(t_next_spike, log_rand_num,
      t_in_episode, tlpar->Tnw +tlpar->Telev + tlpar->Tw, tlpar->Av, tlpar->pr,
      runpar, fl))
    {
      if (tlpar->pr) fprintf(fl.out, "find_tns=%d t_next_spike=%lf\n",
        find_tns, *t_next_spike);
      return(find_tns);
    }
  }

  if (tlpar->pr)
    fprintf(fl.out, "no spike detected t_next_spike=%lf\n", *t_next_spike);

  return(find_tns);
 }

/* This function works if the present time is in the non-whisking episode.  */
/* it finds the time of the next spike or increases log_rand_num within     */
/* the time interval until the end of the non-whisking episode.             */ 
int find_tns_or_update_lrn_nw(double *t_next_spike, double *log_rand_num, 
    double t_in_episode, double Tend, double AA, int pr, run_par *runpar, 
    fl_st fl)
{
  double delta_t_next_spike_new, LQ_n_epi;
  int find_tns;

  LQ_n_epi = AA * (Tend - t_in_episode);

  if (*log_rand_num > -LQ_n_epi)
  {
    delta_t_next_spike_new = -(*log_rand_num) / AA;
    *t_next_spike += delta_t_next_spike_new;

    if (pr)
    {
      fprintf(fl.out, "> end LQ_n_epi=%lf\n", LQ_n_epi);
      fprintf(fl.out, "delta_t_next_spike_new=%lf t_next_spike=%lf\n",
      delta_t_next_spike_new, *t_next_spike);
    }
    find_tns = 1;
  }
  else
  {
    *t_next_spike += Tend - t_in_episode;
    *log_rand_num += LQ_n_epi;

    if (pr)
    {
      fprintf(fl.out, "< cont LQ_n_epi=%lf\n", LQ_n_epi);
      fprintf(fl.out, "log_rand_num=%lf t_next_spike=%lf\n",
      *log_rand_num, *t_next_spike);
    }
    find_tns = 0;
  }

  return(find_tns);
}

/* This function works if the present time is in the whisking episode.      */
/* it finds the time of the next spike or increases log_rand_num within     */
/* the time interval until the end of the whisking episode.                 */ 
int find_tns_or_update_lrn_w(double *t_next_spike, double *log_rand_num, 
    double t_in_episode, tl_par *tlpar, run_par *runpar, fl_st fl)
{
  /* in_contact, cnc: 0 - nc, 1 - c */
  double t_in_episode_run, cnc, tc_updated, LQc, delta_t_next_spike_new;
  int find_tns;

  t_in_episode_run = t_in_episode;

  while (t_in_episode_run < tlpar->Tnw + tlpar->Tw - runpar->epsilon)
  {
    determine_contact_next_change(t_in_episode_run, &cnc, &tc_updated, tlpar,
    runpar, fl);
    LQc = compute_LQ(t_in_episode_run, tc_updated, cnc, tlpar, fl);

    if (tlpar->pr)
    {
      fprintf(fl.out, "t_in_episode_run=%lf log_rand_num=%lf\n",
      t_in_episode_run, *log_rand_num);
      fprintf(fl.out, "tc_updated=%lf cnc=%lf LQc=%lf\n", tc_updated, cnc, LQc);
    }

    if (LQc <= *log_rand_num)
    {
      t_in_episode_run = find_time_LQ(t_in_episode_run, tc_updated, cnc,
       *log_rand_num, tlpar, fl);

      if (tlpar->pr)
      {
        fprintf(fl.out, "> end t_in_episode_run=%lf t_in_episode=%lf\n", 
        t_in_episode_run, t_in_episode);
      }

      if (t_in_episode_run <= tlpar->Tnw + tlpar->Tw + runpar->epsilon)
      {
        /* spike is located in this episode during whisking */
        *t_next_spike += t_in_episode_run - t_in_episode;
        find_tns = 1;
      
        if (tlpar->pr)
        {
          fprintf(fl.out, ">real t_next_spike=%lf\n", *t_next_spike);
	}

        return(find_tns);
      }
      else
      {
	/* no spike is located in this episode during whisking */
        printf("t_in_episode_run=%lf > Tnw=%lf + Tw=%lf\n", t_in_episode_run,
	tlpar->Tnw, tlpar->Tw);
        exit(0);
      }
    }
    else
    {
      *log_rand_num -= LQc;
      t_in_episode_run = tc_updated;
      if (tlpar->pr)
      {
        fprintf(fl.out, "<cont log_rand_num=%lf t_in_episode_run=%lf\n",
        *log_rand_num, t_in_episode_run);
      }
    }
  }

  *t_next_spike += t_in_episode_run - t_in_episode;
  find_tns = 0;
  
  if (tlpar->pr)
  {
    fprintf(fl.out, ">real t_next_spike=%lf\n", *t_next_spike);
  }

  return(find_tns);
}

/* This function works if the present time is in the elevation episode.     */
/* it finds the time of the next spike or increases log_rand_num within     */
/* the time interval until the end of the elevation episode.                */ 
int find_tns_or_update_lrn_elev(double *t_next_spike, double *log_rand_num, 
    double t_in_episode, tl_par *tlpar, run_par *runpar, fl_st fl)
{
  double delta_t_next_spike_new, LQ_e_epi;
  double area, slope, aux, disc, BB;
  int find_tns;

  LQ_e_epi = tlpar->Avnw * (tlpar->Tnw + tlpar->Telev - t_in_episode) + 0.5 * 
    (tlpar->Av - tlpar->Avnw) * (tlpar->Tnw + tlpar->Telev + t_in_episode - 
    2.0 * tlpar->Tnw) * (tlpar->Tnw + tlpar->Telev - t_in_episode) / 
    tlpar->Telev;

  if (*log_rand_num > -LQ_e_epi)
  {
    area = -(*log_rand_num);
    slope = (tlpar->Av - tlpar->Avnw) / tlpar->Telev;

    if (fabs(slope) > runpar->epsilon)
    {    
      aux = tlpar->Avnw + slope * (t_in_episode - tlpar->Tnw);
      disc = 2.0 * area * slope + aux * aux;
      if (disc < 0.0)
      {
        printf("disc=%lf should be positive!\n", disc);
        exit(0);
      }

      BB = tlpar->Avnw - slope * tlpar->Tnw;
      delta_t_next_spike_new = ((-BB + sqrt(disc)) / slope) - t_in_episode;
    }
    else
    {
      delta_t_next_spike_new = -(*log_rand_num) / tlpar->Avnw;
    }

    *t_next_spike += delta_t_next_spike_new;

    if (tlpar->pr)
    {
      fprintf(fl.out, "> end LQ_e_epi=%lf\n", LQ_e_epi);
      fprintf(fl.out, "delta_t_next_spike_new=%lf t_next_spike=%lf\n",
      delta_t_next_spike_new, *t_next_spike);
    }
    find_tns = 1;
  }
  else
  {
    *t_next_spike += tlpar->Tnw + tlpar->Telev - t_in_episode;
    *log_rand_num += LQ_e_epi;

    if (tlpar->pr)
    {
      fprintf(fl.out, "< cont LQ_e_epi=%lf\n", LQ_e_epi);
      fprintf(fl.out, "log_rand_num=%lf t_next_spike=%lf\n",
      *log_rand_num, *t_next_spike);
    }
    find_tns = 0;
  }

  return(find_tns);
}

void determine_contact_next_change(double t_in_episode_run, double *cnc,
     double *tc_updated, tl_par *tlpar, run_par *runpar, fl_st fl)
{
  double t_in_whisking_period, t_beg_per;

  t_in_whisking_period = fmod(t_in_episode_run - tlpar->Tnw, tlpar->Tper);

  if (t_in_whisking_period < tlpar->tc)
  {
    t_in_whisking_period += tlpar->Tper;
  }

  t_beg_per = ((int) ((t_in_episode_run - tlpar->Tnw) / tlpar->Tper)) * 
  tlpar->Tper;

  if ((t_in_whisking_period >= tlpar->tc) && 
      (t_in_whisking_period < tlpar->tc + tlpar->tauc))
  {
    *cnc = 1.0;
    *tc_updated = tlpar->Tnw + t_beg_per + tlpar->tc + tlpar->tauc; 
  }
  else
  {
    *cnc = 0.0;
    *tc_updated = tlpar->Tnw + t_beg_per + tlpar->tc;   
  }

  if (*tc_updated < t_in_episode_run)
    *tc_updated += tlpar->Tper;

  if (*tc_updated > tlpar->Tnw + tlpar->Tw)
    *tc_updated = tlpar->Tnw + tlpar->Tw;

  if (tlpar->pr)
  {
    fprintf(fl.out, "deter\nt_in_episode_run=%lf t_in_whisking_period=%lf "
    "t_beg_per=%lf\n", t_in_episode_run, t_in_whisking_period, t_beg_per);
    fprintf(fl.out, "cnc=%lf tc_updated=%lf\n", *cnc, *tc_updated);
    fflush(fl.out);
  }

  if (*tc_updated < t_in_episode_run)
  {
    printf("tc_updated=%lf < t_in_episode_run=%lf\n", *tc_updated,
    t_in_episode_run);
    exit(0);
  }
}

/* This function computes the log of Q.                               */
/* cnc=0: no contact. cnc=1: contact.                                 */ 
double compute_LQ(double t_in_cycle, double tc_updated, double cnc, 
       tl_par *tlpar, fl_st fl)
{
  double AvpCvtauc, costip, costcu, LQnc, phi, diff_cos;
  double tsum, tsum_mod, tdiff;

  if (tc_updated < t_in_cycle)
  {
    printf("tc_updated=%lf < t_in_cycle=%lf\n", tc_updated, t_in_cycle);
    exit(0);
  }

  AvpCvtauc = tlpar->Av + cnc * tlpar->Cv[ tlpar->itl] *
    (1000.0 / tlpar->Tper) / (tlpar->tauc / tlpar->Tper);
  phi = tlpar->phi[tlpar->itl];

  /*
  costip = cos((2.0 * Pi * t_in_cycle / tlpar->Tper) + phi * Pi);
  costcu = cos((2.0 * Pi * tc_updated / tlpar->Tper) + phi * Pi);
  diff_cos = costcu - costip;
  LQnc = -AvpCvtauc * (tc_updated - t_in_cycle) + tlpar->AvBcTo2p * diff_cos;
  fprintf(fl.out, "c diff_cos=%20.15lf LQnc=%lf\n", diff_cos, LQnc);
  */

  tsum = tc_updated + t_in_cycle;
  tsum_mod = tsum - (((int) (tsum / tlpar->Tper / 2.0)) * 2.0 * tlpar->Tper);
  tdiff =  tc_updated - t_in_cycle;
  diff_cos = -2.0 * sin((Pi * tsum_mod / tlpar->Tper) + phi * Pi) * 
    sin(Pi * tdiff / tlpar->Tper) ;
  LQnc = -AvpCvtauc * (tc_updated - t_in_cycle) + tlpar->AvBcTo2p * diff_cos;

  /*
  fprintf(fl.out, "tsum=%lf tsum_mod=%lf\n", tsum, tsum_mod);
  fprintf(fl.out, "s diff_cos=%20.15lf LQnc=%lf\n", diff_cos, LQnc);
  fprintf(fl.out, "t=%lf %lf AvpCvtauc=%lf costip=%lf costcu=%lf LQnc=%lf\n",
  t_in_cycle, tc_updated, AvpCvtauc, costip, costcu, LQnc);
  */

  return LQnc;
}

/* This function solves the equation:                              */
/* log(Q(t1, t2)) = log_rand_num                                   */
double find_time_LQ(double t_in_cycle, double tc_updated, double cnc,
        double log_rand_num, tl_par *tlpar, fl_st fl)
{
  transfer_to_func ttfunc;
  double tol;
  double time_LQ;
  void *ptr;

  /*
  fprintf(fl.out, "find_time: t_in_cycle=%lf tc_updated=%lf cnc=%lf "
  "\nlog_rand_num=%lf\n", t_in_cycle, tc_updated, cnc, log_rand_num);
  */

  tol=1.0e-12;

  ttfunc.tlpar = tlpar;
  ttfunc.t_in_cycle = &t_in_cycle;
  ttfunc.tc_updated = &tc_updated;
  ttfunc.cnc = &cnc;
  ttfunc.log_rand_num = &log_rand_num;
  ttfunc.fl = &fl;
  ptr = (void*) (&ttfunc);
  time_LQ = zbrent(LQ_func, t_in_cycle, tc_updated, tol, ptr);

  /*
  fprintf(fl.out, "time_LQ=%lf\n", time_LQ);
  */

  return(time_LQ);
}

/* This function is used by zbrent to compute the solution of the equation  */
/* LQ = log_rand_num .                                                      */
double LQ_func(double tLQ, void *ptr)
{
  transfer_to_func *ttfunc; 
  tl_par *tlpar;
  double t_in_cycle, cnc, log_rand_num, LQ_cal;
  fl_st fl;

  ttfunc = (transfer_to_func*) ptr;
  tlpar = ttfunc->tlpar;
  t_in_cycle = *(ttfunc->t_in_cycle);
  cnc = *(ttfunc->cnc);
  log_rand_num = *(ttfunc->log_rand_num);
  fl = *(ttfunc->fl);

  /*
  fprintf(fl.out, "LQ_func: t_in_cycle=%lf cnc=%lf "
  "log_rand_num=%lf\n", t_in_cycle, cnc, log_rand_num);
  fflush(fl.out);
  */

  LQ_cal = compute_LQ(t_in_cycle, tLQ, cnc, tlpar, fl);

  if (tlpar->pr)
    fprintf(fl.out, "z2 tLQ=%18.14lf log=%18.14lf LQ_cal=%18.14lf "
    "diff=%18.14lf\n", tLQ, log_rand_num, LQ_cal, LQ_cal - log_rand_num);

  return(LQ_cal - log_rand_num);
}

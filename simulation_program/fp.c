/* Finiding the fixed point and its stabiliy   */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "tca.h"
#include "tcn.h"
#include "tcwb.h"
#include "fp.h"

/* This function finds a fixed point of the differential equations and */
/* calculates its stability                                            */
void find_fixed_point(func_cell_model *fcm, cell_par *cpar, run_par *runpar,
     syn_par_all *inter, syn_receptor_par *AMPA, double *Varbar, int ipop, int ion,
     fl_st fl)
{
  transfer_to_func_fp ttfunc_fp, *ttf_fp_ptr;
  double *neuron_vec, *neuron_vec_old;
  double tolx, tolf;
  int ieq, ntrial, num_pos_eigval, singular;
  void *ptr;

  /* Finding the fixed point */

  neuron_vec = dvector(1, cpar->neq);
  neuron_vec_old = dvector(1, cpar->neq);

  for (ieq=1; ieq<=cpar->neq; ieq++) neuron_vec[ieq] = Varbar[ieq];

  ttfunc_fp.fcm = fcm;
  ttfunc_fp.cpar = cpar;
  ttfunc_fp.runpar = runpar;
  ttfunc_fp.inter = inter;
  ttfunc_fp.AMPA = AMPA;
  ttfunc_fp.ipop = ipop;
  ttfunc_fp.ion = ion;
  ttfunc_fp.fl = &fl;
  ttf_fp_ptr = &ttfunc_fp;
  ptr = (void *) ttf_fp_ptr;

  ntrial=500;
  tolx=1.0e-8;
  tolf=1.0e-8;

  /* Newton-Raphson method */
  mnewt(ntrial, neuron_vec, cpar->neq, tolx, tolf, &singular, usrfun, ptr);

  if (singular == 1)   /* mnewt returned a solution */
  {
    if ((neuron_vec[1] < cpar->VK) || (neuron_vec[1] > cpar->VNa))
    {
      singular = -1;
      fprintf(fl.out, "mnewt returns V out of range\n");

    }
    else
    {
      for (ieq=1; ieq<=cpar->neq; ieq++)
      {
        neuron_vec_old[ieq] = neuron_vec[ieq];
      }

      /*
      if (!runpar->sm)
      {
        fprintf(fl.out, "ipop=%d ion=%d FP:", ipop, ion);
        for (ieq=1; ieq<=cpar->neq; ieq++)
        {
          fprintf(fl.out, "%lf ", neuron_vec[ieq]);
        }
        if (!runpar->sm) fprintf(fl.out, "\n");
      }
      */
    }
  }
  else if (singular == 0)
  {
    fprintf(fl.out, "mnewt failed");
  }
  else
  {
    printf("singular=%d should be 0 or 1\n", singular);
    exit(0);
  }

  for (ieq=1; ieq<=cpar->neq; ieq++)
  {
     Varbar[ieq] = neuron_vec_old[ieq];
  }

  free_dvector(neuron_vec, 1, cpar->neq);
  free_dvector(neuron_vec_old, 1, cpar->neq);
}

/* User function for Newton-Raphson routine mnewt */
void usrfun(double *xvec, int nvec, double *fvec, double **fjac, void *ptr)
{
  fl_st *fl;

  fl = ((transfer_to_func_fp *) ptr)->fl;
  vecfunc(nvec, xvec, fvec, ptr);
  fdjac(nvec, xvec, fvec, fjac, vecfunc, ptr);
}

/* User function for Newton-Raphson routines mnewt and fdjac */
void vecfunc(int nvec, double *xvec, double *fvec, void *ptr)
{
  func_cell_model *fcm;
  cell_par *cpar;
  syn_par_all *inter;
  syn_receptor_par *AMPA;
  run_par *runpar;
  fl_st *fl;
  double *Varbar, *kout;
  double time;
  double Iapp_now, Isyn_cont, Iel_cont ;
  int ipop, ion, it;
  int ieq;
  
  fcm = ((transfer_to_func_fp *) ptr)->fcm;
  cpar = ((transfer_to_func_fp *) ptr)->cpar;
  runpar = ((transfer_to_func_fp *) ptr)->runpar;
  inter = ((transfer_to_func_fp *) ptr)->inter;
  AMPA = ((transfer_to_func_fp *) ptr)->AMPA;
  fl = ((transfer_to_func_fp *) ptr)->fl;
  ipop = ((transfer_to_func_fp *) ptr)->ipop;
  ion = ((transfer_to_func_fp *) ptr)->ion;

  Varbar = dvector(1, cpar->neq);
  kout   = dvector(1, cpar->neq);

  for (ieq=1; ieq<=cpar->neq; ieq++)
  {
      Varbar[ieq] = xvec[ieq];
  }

  time = 0.0;
  Iapp_now = 0.0;
  Isyn_cont = 0.0;
  Iel_cont = 0.0;
  it = 0;

  fcm->update_cell(Varbar, kout, ipop, ion, cpar, inter, AMPA,
    Iapp_now, Isyn_cont, Iel_cont, it, *fl);
  /*
  one_integration_step(fcm, cpar, runpar, Varbar, kin, kout, 0.0, 0, time,
  Varc, *fl);
  */
  for (ieq=1; ieq<=cpar->neq; ieq++)
  {
      fvec[ieq] = kout[ieq];
  }

  free_dvector(Varbar, 1, cpar->neq);
  free_dvector(kout  , 1, cpar->neq);
}

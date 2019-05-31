/* Functions for the E cells - GA model.                                  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "tca.h"
#include "tcn.h"
#include "tcwb.h"

/* This function substitutes the functions for the WB model */
void substitute_function_WB(func_cell_model *fcm, fl_st fl)
{
  fcm->read_cell_par = &read_cell_par_WB;
  fcm->steady_state_var = &steady_state_var_WB;
  fcm->update_cell = &update_cell_WB;
}

/* This function reads the data for the E cells */
void read_cell_par_WB(cell_par *cpar, fl_st fl)
{
  
  if (fscanf(fl.tmp, "gNa=%lf gKdr=%lf gKz=%lf gL=%lf DelgL=%lf Iext=%lf "
    "DelIext=%lf fracIext=%lf\n", &cpar->gNa, &cpar->gKdr, &cpar->gKz, 
     &cpar->gL, &cpar->DelgL, &cpar->Iext, &cpar->DelIext, &cpar->fracIext)
    != 0)
  {
    fprintf(fl.out, "gNa=%lf gKdr=%lf gKz=%lf gL=%lf DelgL=%lf Iext=%lf "
    "DelIext=%lf fracIext=%lf\n", cpar->gNa, cpar->gKdr, cpar->gKz, 
     cpar->gL, cpar->DelgL, cpar->Iext, cpar->DelIext, cpar->fracIext);

    if (cpar->gL <= 1.0e-10)
    {
      printf("gL should be positive\n");
      exit(0);
    }
  }
  else
    printf("gNa\n");

  /* I_app_read(cpar, fl); */

  if (fscanf(fl.tmp, "Cm=%lf VNa=%lf VK=%lf VL=%lf Vopto=%lf phi=%lf "
    "gamma=%lf gamma_for_syn=%c\n", &cpar->Cm, &cpar->VNa, &cpar->VK, &cpar->VL,
    &cpar->Vopto, &cpar->phi, &cpar->gamma, &cpar->gamma_for_syn) != 0)
    {
      cpar->two_pow_gamma = pow(2.0, cpar->gamma);
      fprintf(fl.out, "Cm=%lf VNa=%lf VK=%lf VL=%lf Vopto=%lf\n", cpar->Cm, 
      cpar->VNa, cpar->VK, cpar->VL, cpar->Vopto);
      fprintf(fl.out, "phi=%lf gamma=%lf two_pow_gamma=%lf "
      "gamma_for_syn=%c\n", cpar->phi, cpar->gamma,
      cpar->two_pow_gamma, cpar->gamma_for_syn);
    }
  else
    printf("Cm\n");

  /* Defining the structure of the cell's compartments and variables */
 
  cpar->nceq = 4;    /* number of cellular variables: V, h, n, z */
  cpar->nseq = 0;    /* number of synaptic variables: sP, sNa, sNb, sA */

  cpar->neq = cpar->nceq + cpar->nseq;
  fprintf(fl.out, "E: nceq=%d nseq=%d neq=%d \n", cpar->nceq, cpar->nseq,
  cpar->neq);
  fflush(fl.out);
}

/* This function calculates the steady-state variables for a specific V */
void steady_state_var_WB(cell_par *cpar, double *Varb, run_par *runpar, 
     fl_st fl)
{
  double alphahs, betahs, Hinfs, tauHs;
  double alphans, betans, Ninfs, tauNs;
  double Vc;
  int ieq;

  Vc = Varb[1];                                                /* V */

  alphahs = 0.7*exp(-(Vc+44.0)/20.0);
  betahs  = 10.0/(1.0+exp(-(Vc+14.0)/10.0));
  Varb[2] = alphahs/(alphahs+betahs);                          /* h */

  alphans = 0.1*(Vc+34.0)/(1.0-exp(-(Vc+34.0)/10.0));
  betans  = 1.25*exp(-(Vc+44.0)/80.0);
  Varb[3] = alphans/(alphans+betans);                          /* n */

  Varb[4] = 1.0 / (1.0 + exp(-0.7 * (Vc + 30.0)));             /* z */

  for (ieq=cpar->nceq+1; ieq<=cpar->neq; ieq++) Varb[ieq] = 0.0;
}

/* This function updates the variables of an excitatory cell   */
void update_cell_WB(double *Varc, double *kout, int ipop, int ion, 
     cell_par *cpar, syn_par_all *inter, syn_receptor_par *AMPA,
     double Iapp_now, double Isyn_cont, double Iel_cont, int it, fl_st fl)
{
  double Vc, hc, nc, zc, tc, sc, s_inf;
  double alphams, betams, Minfs, alphahs, betahs, Hinfs, tauHs;
  double alphans, betans, Ninfs, tauNs;
  double Zinfs, tauZs;
  double INa, IKdr, IKz, Iback, Iopto ;
  double I_sum_intrinsic, I_sum_synaptic;
  int ieq;
  
  Vc = Varc[1];
  hc = Varc[2];
  nc = Varc[3];
  zc = Varc[4];

  alphams = (fabs(Vc+30.0) > 10e-8) ? 
               0.1*(Vc+30.0)/(1.0-exp(-(Vc+30.0)/10.0)) : 
               10.0;
  betams  = 4.0*exp(-(Vc+55.0)/18.0);
  Minfs   = alphams/(alphams+betams);

  if (Vc > -120.0)
  {  
    alphahs = 0.7*exp(-(Vc+44.0)/20.0);
    betahs  = 10.0/(1.0+exp(-(Vc+14.0)/10.0));
    Hinfs   = alphahs/(alphahs+betahs);
    tauHs   = 1.0 / (cpar->phi * (alphahs+betahs));


    alphans = 0.1*(Vc+34.0)/(1.0-exp(-(Vc+34.0)/10.0));
    betans  = 1.25*exp(-(Vc+44.0)/80.0);
    Ninfs   = alphans/(alphans+betans);
    tauNs   = 1.0 / (cpar->phi * (alphans+betans));

    Zinfs   = 1.0 / (1.0 + exp(-0.7 * (Vc + 30.0)));
    tauZs   = 60.0;
  }
  else
  {
    Hinfs   = 1.0;
    tauHs   = 1.0;
    Ninfs   = 0.0;
    tauNs   = 1.0;
    Zinfs   = 0.0;
    tauZs   = 60.0;
  }

  INa  = cpar->gNa * Minfs * Minfs * Minfs * hc * (Vc - cpar->VNa);
  IKdr = cpar->gKdr * nc * nc * nc * nc * (Vc - cpar->VK);
  IKz  = cpar->gKz * zc * (Vc - cpar->VK);
  Iopto = cpar->opt_con[ion] * (Vc - cpar->Vopto);
  Iback = cpar->gback_nor * (inter->rho_concur * (Vc - AMPA->Vrev) +
	  (1.0 - inter->rho_concur) * (cpar->VL - AMPA->Vrev));

  I_sum_intrinsic = -cpar->gLar[ion] * (Vc - cpar->VL) - INa - IKdr - IKz - 
          Iopto + cpar->Iextar[ion] + Iapp_now;
  I_sum_synaptic = -Iback + Isyn_cont + Iel_cont;

  if (cpar->gamma_for_syn == 'y')
  {
    kout[1] = (I_sum_intrinsic + I_sum_synaptic) * cpar->two_pow_gamma / 
    cpar->Cm;                                                 /* V */
  }
  else if (cpar->gamma_for_syn == 'n')                        
  {
    kout[1] = (I_sum_intrinsic * cpar->two_pow_gamma + I_sum_synaptic) / 
    cpar->Cm;                                                 /* V */
  }
                              

  /*  if ((ipop == 2) && (ion == 150))
  {
        printf("V=%lf Vl=%lf gL=%lf k1=%lf Iback=%lf IL=%lf\n", Vc, cpar->VL, cpar->gLar[ion], kout[1], Iback, cpar->gLar[ion] * (Vc - cpar->VL));
	  printf("V=%lf gKdr=%lf nc=%lf VK=%lf IKdr=%lf\n", Vc, cpar->gKdr, nc, cpar->VK, IKdr);
  } */


  kout[2] = cpar->two_pow_gamma * (Hinfs - hc) / tauHs;       /* h */

  kout[3] = cpar->two_pow_gamma * (Ninfs - nc) / tauNs;       /* n */

  kout[4] = cpar->two_pow_gamma * (Zinfs - zc) / tauZs;       /* z */

  /*
  if ((ipop == 1) && (ion == 1))
    fprintf(fl.col, "Iback=%lf Iopto=%lf Isyn_cont=%lf\n", Iback, Iopto,
    Isyn_cont);
  */
  /*
  if ((ipop == 1) && (ion == 1)) fprintf(fl.col, "%lf ", Iback);
  fprintf(fl.col, "kout=");
  for (ieq=1; ieq<=4; ieq++) fprintf(fl.col, "%lf ", kout[ieq]);
  fprintf(fl.col, "Minfs=%lf mqh=%lf INa=%lf\n", Minfs, 
  Minfs * Minfs * Minfs * hc, INa);
  fprintf(fl.col, "\n");
  */
}

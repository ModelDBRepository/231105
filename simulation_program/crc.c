/* This program computes cross-correlations */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "crc.h"

int main(int argc, char*argv[])
{
  crc_par crcpar;
  fl_st fl;
  double **datar, *corar;
  int ic;
  char suffix[7]="a1", file_name[36], fnm[8];

  if (argc >= 2) strcpy(suffix,argv[1]);
  fl.in  = fopen("crc.n", "r");
  fl.dat = fopen(strcat(strcpy(file_name,"tc.col."),suffix), "r");
  fl.cor = fopen(strcat(strcpy(file_name,"crc.cor."),suffix), "w");
  fl.out = fopen(strcat(strcpy(file_name,"crc.out."),suffix), "w");
  fl.res = fopen(strcat(strcpy(file_name,"crc.res."),suffix), "w");

  read_input(&crcpar, fl);
  write_input(&crcpar, fl);

  datar  = dmatrix(1, 2, 1, crcpar.n_line_read);
  corar  = dvector(1 - crcpar.n_line_read, crcpar.n_line_read - 1);

  for (ic=1; ic<=crcpar.nicol2ar; ic++)
  {
    printf("ic=%d nicol2ar=%d\n", ic, crcpar.nicol2ar);
    crcpar.icol[2] = crcpar.icol2ar[ic];
    if ((ic > 1) && (!crcpar.smap)) fprintf(fl.cor, "  \n");
    one_cor_cal_quantify(&crcpar, datar, corar, fl);
  }

  free_dmatrix(datar, 1, 2, 1, crcpar.n_line_read);
  free_dvector(corar, 1 - crcpar.n_line_read, crcpar.n_line_read - 1);
  free_ivector(crcpar.icol2ar, 1, crcpar.nicol2ar);

  fclose(fl.in);
  fclose(fl.dat);
  fclose(fl.cor);
  fclose(fl.out);
  fclose(fl.res);
}

/* This function read the input parameters */
void read_input(crc_par *crcpar, fl_st fl)
{
  int ic;
  int fsc;
  crcpar->epsilon = 1.0e-10;

  fsc = fscanf(fl.in, "ncol=%d detcol=%c\n", &crcpar->ncol, &crcpar->detcol);

  if (crcpar->detcol == 'o')
  {
    fsc = fscanf(fl.in, "icol1=%d icol2=%d\n", &crcpar->icol[1],
    &crcpar->icol[2]);
    crcpar->nicol2ar = 1;
    crcpar->icol2ar = ivector(1, crcpar->nicol2ar);
    crcpar->icol2ar[1] = crcpar->icol[2];
  }
  else if (crcpar->detcol == 'l')
  {
    fsc = fscanf(fl.in, "icol1=%d nicol2ar=%d icol2ar=", &crcpar->icol[1],
    &crcpar->nicol2ar);
    crcpar->icol2ar = ivector(1, crcpar->nicol2ar);
    for (ic=1; ic<=crcpar->nicol2ar; ic++)
    {
      fsc = fscanf(fl.in, "%d", &crcpar->icol2ar[ic]);
      if (ic < crcpar->nicol2ar) fsc = fscanf(fl.in, " ");
    }
    fsc = fscanf(fl.in, "\n");
  }

  fsc = fscanf(fl.in, "n_line_skip=%d n_line_read=%d smap=%d\n", 
  &crcpar->n_line_skip, &crcpar->n_line_read, &crcpar->smap);
  fsc = fscanf(fl.in, "Tmax=%lf T_between_peaks=%lf\n", &crcpar->Tmax,
  &crcpar->T_between_peaks);
  fsc = fscanf(fl.in, "kpoints=%d base_intersect=%lf inter_points=%d\n",
  &crcpar->kpoints, &crcpar->base_intersect, &crcpar->inter_points);
}

/* This function writes the input parameters */
void write_input(crc_par *crcpar, fl_st fl)
{
  int ic;

  fprintf(fl.out, "ncol=%d detcol=%c\n", crcpar->ncol, crcpar->detcol);
  fprintf(fl.out, "icol1=%d nicol2ar=%d icol2ar=", crcpar->icol[1],
  crcpar->nicol2ar);
  
  for (ic=1; ic<=crcpar->nicol2ar; ic++)
  {
    fprintf(fl.out, "%d", crcpar->icol2ar[ic]);
    if (ic < crcpar->nicol2ar) fprintf(fl.out, " ");
  }
  fprintf(fl.out, "\n");

  fprintf(fl.out, "n_line_skip=%d n_line_read=%d smap=%d\n",
  crcpar->n_line_skip, crcpar->n_line_read, crcpar->smap);
  fprintf(fl.out, "Tmax=%lf T_between_peaks=%lf\n", crcpar->Tmax,
  crcpar->T_between_peaks);
  fprintf(fl.out, "kpoints=%d base_intersect=%lf inter_points=%d\n",
  crcpar->kpoints, crcpar->base_intersect, crcpar->inter_points);
  fprintf(fl.out, "\n");
  fflush(fl.out);
}

/* This function calculate the cross-correlation and quantify it  */
/* for two specific columns.                                      */
void one_cor_cal_quantify(crc_par *crcpar, double **datar, double *corar,
     fl_st fl)
{
  read_data(crcpar, datar, fl);
  calculate_correlation(crcpar, datar, 2, 2, corar, fl);
  /* amp_ph_cal(crcpar, 2, 2, corar, fl); */
  calculate_correlation(crcpar, datar, 1, 2, corar, fl);
  /* amp_ph_cal(crcpar, 1, 2, corar, fl); */
}

/* This function read the data file */
void read_data(crc_par *crcpar, double **datar, fl_st fl)
{
  double *vec, time1, time2;
  int iread, icol;
  char line[Mline], *p1;

  rewind(fl.dat);
  vec = dvector(1, crcpar->ncol);

  iread = 1;
  while (iread <= crcpar->n_line_skip) 
  {
    if (fgets(line, Mline, fl.dat) != NULL)
    {
      iread++;
    }
    else
    {
      break;
    }
  }

  if (iread != crcpar->n_line_skip + 1)
  {
    printf("iread=%d != n_line_skip=%d + 1\n", iread, crcpar->n_line_skip);
    exit(0);
  }

  iread  = 1;
  while (iread <= crcpar->n_line_read) 
  {
    if (fgets(line, Mline, fl.dat) != NULL)
    {
      p1= strtok(line, " ");
      if (iread == 1) sscanf(p1, "%lf\n", &time1);
      else if (iread == 2) 
      {
        sscanf(p1, "%lf\n", &time2);
        crcpar->deltat = time2 - time1;
        fprintf(fl.out, "deltat=%lf\n", crcpar->deltat);
      }

      sscanf(p1, "%lf\n", &vec[1]);
      for(icol=2; icol<=crcpar->ncol; icol++)
      {
        p1= strtok(NULL, " ");
        sscanf(p1, "%lf\n", &vec[icol]);
      }

      datar[1][iread] = vec[crcpar->icol[1]];
      datar[2][iread] = vec[crcpar->icol[2]];
      iread++;
    }
    else
    {
      printf("iread=%d not enough data!\n", iread);
      break;
    }
  }

  if (iread != crcpar->n_line_read + 1)
  {
    printf("iread=%d != n_line_read=%d + 1\n", iread, crcpar->n_line_read);
    exit(0);
  }

  /*
  for (iread=1; iread<=crcpar->n_line_read; iread++)
    fprintf(fl.out, "%d %lf %lf\n", iread, datar[1][iread], datar[2][iread]);
  fprintf(fl.out, "\n");
  */
  free_dvector(vec, 1, crcpar->ncol);
}

/* This function calculates the cross-correlation */
void calculate_correlation(crc_par *crcpar, double **datar, int idat1,
     int idat2, double *corar, fl_st fl)
{
  double **var, *ans, V1av, V2av, xcc;
  int npt, ipt, jpt;
  int ido;

  npt = crcpar->n_line_read;
  var = dmatrix(1, 2, 1, 2*npt);
  ans = dvector(1, 4*npt);

  /*  calculating the non-connected correlation <V_1 V_2> */

  for (ipt=1; ipt<=npt; ipt++) var[1][ipt] = datar[idat1][ipt];
  for (ipt=1; ipt<=npt; ipt++) var[1][npt+ipt] = 0.0;
  for (ipt=1; ipt<=npt; ipt++) var[2][ipt] = datar[idat2][ipt];
  for (ipt=1; ipt<=npt; ipt++) var[2][npt+ipt] = 0.0;


  fprintf(fl.out, "var\n");
  for (ipt=1; ipt<=2*npt; ipt++)
    fprintf(fl.out, "%d %lf %lf\n", ipt, var[1][ipt], var[2][ipt]);
  fprintf(fl.out, "\n");

  xcorrel(var[1], var[2], (unsigned long) (2*npt), ans);

  for (ipt=1-npt; ipt<=-1; ipt++)
      corar[ipt] = ans[2*npt+ipt+1]; 
  for (ipt=0; ipt<=npt-1; ipt++)
    corar[ipt] = ans[ipt+1];

  /* calculating the average <V> */

  V1av = 0.0;
  for (ipt=1; ipt<=npt; ipt++) V1av += datar[idat1][ipt];
  V1av /= npt;

  V2av = 0.0;
  for (ipt=1; ipt<=npt; ipt++) V2av += datar[idat2][ipt];
  V2av /= npt;

  fprintf(fl.out, "V1av=%lf V2av=%lf\n",  V1av, V2av);

  fprintf(fl.out, "cc npt=%d smap=%d\n", npt, crcpar->smap);
  for (ipt=1-npt; ipt<=npt-1; ipt++)
  {
    corar[ipt] /= npt - abs(ipt);
    corar[ipt] -= V1av * V2av; 
    if (!crcpar->smap)
    {
      fprintf(fl.cor, "%lf %lf\n", ipt * crcpar->deltat, corar[ipt]);
    }
  }
  fprintf(fl.cor, "  \n");

  free_dmatrix(var, 1, 2, 1, 2*npt);
  free_dvector(ans, 1, 4*npt);
}

/* ------------------------------------------------------------ */

/* This function computes the amplitudes and the phases of the  */
/* cross-correlation                                            */
/* Adapted from ampph.c                                         */
void amp_ph_cal (crc_par *crcpar, int idat1, int idat2, double *corar,
     fl_st fl)
{
  peak_str peakstr;

  crcpar->ntc = (crcpar->Tmax + crcpar->epsilon) / crcpar->deltat;
  fprintf(fl.out, "ntc=%d\n", crcpar->ntc);

  peakstr.mp = (int) (2.0 * crcpar->Tmax / crcpar->T_between_peaks);
  peakstr.iarmx = ivector(1, peakstr.mp);
  peakstr.np = 0;
  peakstr.iarmn = ivector(1, 4);
  peakstr.t_isect = dmatrix(1, 3, 1, 3);

  find_local_max(crcpar, &peakstr, corar, fl);
  find_peak_near_zero(crcpar, &peakstr, corar, fl);
  find_local_min(crcpar, &peakstr, corar, fl);
  find_intersections_with_zero(crcpar, &peakstr, corar, fl);

  free_ivector(peakstr.iarmx, 1, peakstr.mp);
  free_ivector(peakstr.iarmn, 1, 4);
  free_dmatrix(peakstr.t_isect, 1, 3, 1, 3);

  fprintf(fl.res, "%d %d %lf %lf %lf %lf\n", crcpar->icol[idat1],
  crcpar->icol[idat2], peakstr.amplitude, peakstr.t_period, peakstr.t_shift,
  peakstr.phase);
  fflush(fl.res);
}

/* This function calculates the local maxima, within time interval */
/* T_between_peaks                                                 */
void find_local_max(crc_par *crcpar, peak_str *peakstr, double *corar, 
     fl_st fl)
{
  int ii;
  int ipeak;

  for (ii=crcpar->kpoints - crcpar->ntc; ii<=crcpar->ntc - crcpar->kpoints; 
  ii++)
  {
    if (find_very_local_max(ii, crcpar, corar, fl))
    {
      if (!crcpar->smap)
        fprintf(fl.out, "%lf ii=%d\n", corar[ii], ii);
      if (peakstr->np == 0)
      {
        peakstr->np++;
        peakstr->iarmx[peakstr->np] = ii;
      }
      else 
      {
        if (ii * crcpar->deltat -  
            peakstr->iarmx[peakstr->np]* crcpar->deltat < 
            crcpar->T_between_peaks)
        {
          if (corar[ii] >  corar[peakstr->iarmx[peakstr->np]])
          { 
   	    peakstr->iarmx[peakstr->np] = ii;         
	  }
	}
        else
	{
          peakstr->np++;
          if (peakstr->np > peakstr->mp)
          {
            printf("np=%d > mp=%d\n", peakstr->np, peakstr->mp);
            exit(0);
	  }
          peakstr->iarmx[peakstr->np] = ii;
	}
      }
    }
  }

  if (!crcpar->smap)
  {
    fprintf(fl.out, "\n");
    for (ipeak=1; ipeak<=peakstr->np; ipeak++)
    {
      ii = peakstr->iarmx[ipeak];
      fprintf(fl.out, "%lf %lf %d %d\n", ii * crcpar->deltat, corar[ii], ii,
      ipeak);
    }
  }
}

/* This function calculates a very local maximum */
int find_very_local_max(int ii, crc_par *crcpar, double *corar, fl_st fl)
{
  int ll1, ll2, vlm, jj;

  ll1 = max(ii - crcpar->kpoints, -crcpar->ntc);
  ll2 = min(ii + crcpar->kpoints, crcpar->ntc);

  vlm = 1;

  for (jj=ll1; jj<=ii-1; jj++)
  {
    if (corar[jj] > corar[ii])
    {
      vlm=0;
      break;
    }
  }

  if (vlm)
  {
    for(jj=ii+1; jj<=ll2; jj++)
    {
      if (corar[jj] > corar[ii])
      {
        vlm=0;
        break;
      }
    }
  }

  return(vlm);
}

/* This function find the peak that is closest to tau=0.                */
void find_peak_near_zero(crc_par *crcpar, peak_str *peakstr, double *corar,
     fl_st fl)
{
  int ipeak, iold, inew;

  peakstr->zp = 1;

  for (ipeak=2; ipeak<=peakstr->np; ipeak++)
  {
    iold = peakstr->iarmx[peakstr->zp];
    inew = peakstr->iarmx[ipeak];
    if (fabs(inew) < fabs(iold))
    {
      peakstr->zp = ipeak;
    }
  }

  if (!crcpar->smap)
    fprintf(fl.out, "\nzp=%d ii=%d c=%lf %lf\n", peakstr->zp,
    peakstr->iarmx[peakstr->zp], peakstr->iarmx[peakstr->zp] * crcpar->deltat,
    corar[peakstr->iarmx[peakstr->zp]]);

  if (peakstr->zp < 3)
  {
    printf("zp=%d < 3\n", peakstr->zp);
    exit(0);
  }
}

/* This function calculates the four centeral local manima */
void find_local_min(crc_par *crcpar, peak_str *peakstr, double *corar,
     fl_st fl)
{
  int ipeak, jpeak, ileft, iright, ii;

  if (!crcpar->smap) fprintf(fl.out, "\n");
  for (ipeak=1; ipeak<=4; ipeak++)
  {
    jpeak = ipeak + peakstr->zp - 3;
    ileft = peakstr->iarmx[jpeak]+1;
    iright = peakstr->iarmx[jpeak+1]-1;
    
    if (!crcpar->smap)
      fprintf(fl.out, "ipeak=%d jpeak=%d ileft=%d iright=%d\n", ipeak, jpeak, 
      ileft, iright);

    peakstr->iarmn[ipeak] = ileft;
    for (ii=ileft+1; ii<=iright; ii++)
    {
      if (corar[ii] < corar[peakstr->iarmn[ipeak]])
      {
        peakstr->iarmn[ipeak] = ii;
      }
    }

    if (!crcpar->smap)
      fprintf(fl.out, "%lf %lf %d %d\n", peakstr->iarmn[ipeak]*crcpar->deltat,
      corar[peakstr->iarmn[ipeak]], ipeak, peakstr->iarmn[ipeak]);
  }

  if (!crcpar->smap)
  {
    fprintf(fl.out, "\nmax for amp\n");
    fprintf(fl.out, "%lf %lf %d\n", peakstr->iarmx[peakstr->zp-1] *
    crcpar->deltat, corar[peakstr->iarmx[peakstr->zp-1]], peakstr->zp-1);
    fprintf(fl.out, "%lf %lf %d\n", peakstr->iarmx[peakstr->zp+1] *
    crcpar->deltat, corar[peakstr->iarmx[peakstr->zp+1]], peakstr->zp-1);
  }

  peakstr->amplitude = 0.0;
  for (ipeak=-1; ipeak<=1; ipeak+=2)
  {
    peakstr->amplitude += 2.0 * corar[peakstr->iarmx[peakstr->zp+ipeak]];
  }
  for (ipeak=1; ipeak<=4; ipeak++)
  {
    peakstr->amplitude -= corar[peakstr->iarmn[ipeak]];
  }
  peakstr->amplitude /= 4.0;
  if (!crcpar->smap) fprintf(fl.out, "amplitude=%lf\n", peakstr->amplitude);
}

/* This function calculates intersections of the cross-correlation with */
/* zero                                                                 */
void find_intersections_with_zero(crc_par *crcpar, peak_str *peakstr, 
     double *corar, fl_st fl)
{
  int ipeak, pleft, pright;

  /* Finds the closest point to zero */
  if (!crcpar->smap) fprintf(fl.out, "\n");
  for (ipeak=1; ipeak<=3; ipeak++)
  {
    pleft = peakstr->iarmn[ipeak] + 1;
    pright = peakstr->iarmx[peakstr->zp + ipeak - 2] - 1;
    if (!crcpar->smap)
      fprintf(fl.out, "\npleft=%d nleft=%d\n", pleft, pright);
    find_one_intersection(pleft, pright, crcpar, corar,
    &peakstr->t_isect[ipeak][1],fl);
    if (!crcpar->smap) 
      fprintf(fl.out, "t_isect=%lf\n", peakstr->t_isect[ipeak][1]); 
  }

  for (ipeak=1; ipeak<=3; ipeak++)
  {
    pleft = peakstr->iarmx[peakstr->zp + ipeak - 2] + 1;
    pright = peakstr->iarmn[ipeak+1] - 1;
    if (!crcpar->smap)
      fprintf(fl.out, "\npleft=%d nleft=%d\n", pleft, pright);
    find_one_intersection(pleft, pright, crcpar, corar,
    &peakstr->t_isect[ipeak][2],fl);
    if (!crcpar->smap) 
      fprintf(fl.out, "t_isect=%lf\n", peakstr->t_isect[ipeak][2]);
  }

  if (!crcpar->smap) fprintf(fl.out, "\n");
  for (ipeak=1; ipeak<=3; ipeak++)
  {
    peakstr->t_isect[ipeak][3] = (peakstr->t_isect[ipeak][1] +
    peakstr->t_isect[ipeak][2]) / 2.0;
    if (!crcpar->smap)
      fprintf(fl.out, "%lf %lf %lf\n", peakstr->t_isect[ipeak][1], 
      peakstr->t_isect[ipeak][2], peakstr->t_isect[ipeak][3]);
  }

  peakstr->t_period = (peakstr->t_isect[3][3] - peakstr->t_isect[1][3]) / 2.0;
  peakstr->t_shift = (peakstr->t_isect[1][3] + peakstr->t_isect[2][3] +
    peakstr->t_isect[3][3]) / 3.0;
  peakstr->phase = peakstr->t_shift / peakstr->t_period;
  if (!crcpar->smap)
    fprintf(fl.out, "\nt_period=%lf t_shift=%lf phase=%lf\n",
    peakstr->t_period, peakstr->t_shift, peakstr->phase);
}

/* This function finds one intersection of the cross-correlation with zero. */
void find_one_intersection(int pleft, int pright, crc_par *crcpar, 
     double *corar, double *t_isect, fl_st fl)
{
  int ipl, ipr, find_isect_l, find_isect_r;

  /* Finds left intersection point */
  ipl = pleft;
  while (!(find_isect_l = cond_intersect(ipl, crcpar, corar, fl)) && 
         ipl < pright)
  {
   ipl++;
  }

  ipr = pright;
  while (!(find_isect_r = cond_intersect(ipr-1, crcpar, corar, fl)) &&
         ipr > pleft)
  {
   ipr--;
  }

  if (!crcpar->smap)
    fprintf(fl.out, "find_isect_l=%d ipl=%d find_isect_r=%d ipr=%d\n",
  find_isect_l, ipl, find_isect_r, ipr);

  if (find_isect_l && find_isect_r)
  {
    linear_regression(ipl, ipr, crcpar, corar, t_isect, fl);
  }

}

/* This function formulates the condition for intersection with zero  */
/* between two adjacent points.                                       */
int cond_intersect(int ip, crc_par *crcpar, double *corar, fl_st fl)
{
  int cond;

  cond = 0;

  if (fabs(corar[ip]) < crcpar->epsilon)
  {
    cond = 1;
  }
  else if (corar[ip] * corar[ip+1] < 0.0)
  {
    cond = 1;
  }

  return(cond);
}

void linear_regression(int ipl, int ipr, crc_par *crcpar, double *corar,
     double *t_isect, fl_st fl)
{
  double *xx, *yy;
  double sumx, sumxx, sumy, sumxy, denom, aa, bb;
  int lenv, ii;

  lenv = ipr - ipl + 2 * crcpar->inter_points - 1; 
  if (!crcpar->smap) fprintf(fl.out, "lenv=%d\n", lenv);

  xx = dvector(1, lenv);
  yy = dvector(1, lenv);

  for (ii=1; ii<=lenv; ii++)
  {
    xx[ii] = (ipl - crcpar->inter_points + ii) * crcpar->deltat;
    yy[ii] = corar[ipl - crcpar->inter_points + ii];
    if (!crcpar->smap)
      fprintf(fl.out, "ii=%d jj=%d xx=%lf yy=%lf\n", ii, 
      ipl - crcpar->inter_points + ii, xx[ii], yy[ii]);
  }

  sumx  = 0.0;
  sumxx = 0.0;
  sumy  = 0.0;
  sumxy = 0.0;

  for (ii=1; ii<=lenv; ii++)
  {
    sumx  += xx[ii];
    sumxx += xx[ii] * xx[ii];
    sumy  += yy[ii];
    sumxy += xx[ii] * yy[ii];
  }

  sumx  /= lenv;
  sumxx /= lenv;
  sumy  /= lenv;
  sumxy /= lenv;

  denom = sumxx - sumx * sumx;
  if (fabs(denom) < crcpar->epsilon)
  {
    *t_isect = -999.0;
  }
  else
  {
    aa = (-sumxy * sumx + sumxx * sumy) / denom;
    if (fabs(aa) < crcpar->epsilon)
    {
      *t_isect = -998.0;
    }
    else
    {
      bb = (sumxy - sumx * sumy) / denom;
      *t_isect = (crcpar->base_intersect - aa) / bb;
      if (!crcpar->smap)
        fprintf(fl.out, "denom=%lf aa=%lf bb=%lf ti=%lf\n", denom, aa, bb,
        *t_isect);
    }
  }

  free_dvector(xx, 1, lenv);
  free_dvector(yy, 1, lenv);
}

/* CORREL */
void xcorrel(double data1[], double data2[], unsigned long n, double ans[])
{
        unsigned long no2,i;
        double dum,*fft;

        fft=dvector(1,n<<2);
        xtwofft(data1,data2,fft,ans,n);

        no2=n>>1;
        for (i=2;i<=n+2;i+=2) {
            ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
            ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
	}
        ans[2]=ans[n+1];
        realft(ans,n,-1);
        free_dvector(fft,1,n<<2);
}

void xtwofft(double data1[], double data2[], double fft1[], double fft2[],
	unsigned long n)
{
	unsigned long nn3,nn2,jj,j;
	double rep,rem,aip,aim;

	nn3=1+(nn2=2+n+n);

	for (j=1,jj=2;j<=n;j++,jj+=2) {
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}

	xfour1(fft1,n,1);

	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;

	for (j=3;j<=n+1;j+=2) {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j]=aip;
		fft2[j+1] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	}
}

void xfour1(double data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi,temp;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;

	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;

		for (m=1;m<mmax;m+=2) {

			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}

			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}

		mmax=istep;
	}

}

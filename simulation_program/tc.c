/* This program simulates cortical networks composed of conductance-based   */
/* neurons.                                                                 */
/* In this file, the loop over the parameter is executed, and the           */
/* parameters of each run are determined.                                   */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "tca.h"
#include "tcn.h"
#include "tc.h"

int main(int argc, char*argv[])
{
  scan_val sval, svala;
  par_all parall;
  avr_val av;
  fl_st fl;
  int **ran_gen, ngen, length;
  char suffix[10]="a11", file_name[42];

  if (argc >= 2) strcpy(suffix,argv[1]);

  fl.in  = fopen(strcat(strcat(strcpy(file_name,"tc."),"n."   ),suffix), "r");
  fl.tmp = fopen(strcat(strcat(strcpy(file_name,"tc."),"tmp." ),suffix), "w+");
  fl.out = fopen(strcat(strcat(strcpy(file_name,"tc."),"out." ),suffix), "w");
  fl.col = fopen(strcat(strcat(strcpy(file_name,"tc."),"col." ),suffix), "w");
  fl.ras = fopen(strcat(strcat(strcpy(file_name,"tc."),"ras." ),suffix), "w");
  fl.fri = fopen(strcat(strcat(strcpy(file_name,"tc."),"fri." ),suffix), "w");
  fl.zmp = fopen(strcat(strcat(strcpy(file_name,"tc."),"zmp." ),suffix), "w");
  fl.his = fopen(strcat(strcat(strcpy(file_name,"tc."),"his." ),suffix), "w");

  if (fscanf(fl.in, "scan=%c\n", &sval.scan_type) == 0)
  {
    printf("cannot read scan_type\n");
    exit(0);
  }

  fprintf(fl.out, "scan=%c\n", sval.scan_type);
  parall.scan_type = sval.scan_type;

  if (sval.scan_type == 'n')
  {
    /* Simulating one parameter set */
    parall.sm = 0;
    if (fscanf(fl.in, "seed=%d\n", &sval.seed) == 0) exit(0);
    fprintf(fl.out, "seed=%d\n", sval.seed);
    fl.avr = fl.out;
    update_file_old(&sval, 2, fl);

    av.ipar = 0;
    av.irepeat = 0;
    av.par = 0.0;

    ngen = 2;
    length = 17;
    ran_gen = init_rng_s_dbl(ngen, length, sval.seed);
    parall.rng_ptr = *ran_gen;

    one_par(&parall, &av, fl);
    write_avr(&parall, &av, fl);
  }
  else if ((sval.scan_type == 'e') || (sval.scan_type == 'u'))
  {
    /* Simulating several parameter sets */
    parall.sm = 1;
    read_first_input_line(&sval, fl);
    if (fscanf(fl.in, "seed=%d\n", &sval.seed) == 0) exit(0);
    fprintf(fl.out, "seed=%d\n", sval.seed);
    fl.avr = fopen(strcat(strcat(strcpy(file_name,"tc."),"avr." ),suffix),
    "w");

    ngen = 5 * (sval.npar+2) * sval.nrepeat;
    length = 17;
    ran_gen = init_rng_s_dbl(ngen, length, sval.seed);

    for (sval.ipar=0; sval.ipar <= sval.npar; sval.ipar++)
    {
      av.ipar = sval.ipar;
      av.par = sval.par_ar[sval.ipar];
      fprintf(fl.his, "#p %d %lf %d\n", av.ipar, av.par, sval.nrepeat);

      for (sval.irepeat=1; sval.irepeat <= sval.nrepeat; sval.irepeat++)
      {
        printf("ipar=%d npar=%d irepeat=%d nrepeat=%d\n", sval.ipar, 
        sval.npar, sval.irepeat, sval.nrepeat);
        if (sval.ipar > 0) fprintf(fl.out, "\n");

        update_file_old(&sval, 3, fl);

        av.irepeat = sval.irepeat;

        parall.rng_ptr = *(ran_gen + (sval.ipar * sval.nrepeat) + 
                          sval.irepeat - 1);

        fprintf(fl.his, "#r %d\n", av.irepeat);
        one_par(&parall, &av, fl);
        write_avr(&parall, &av, fl);
        fflush(fl.his);
        fflush(fl.avr);
      }
    }

    fclose(fl.avr);
  }

  fclose(fl.in);
  fclose(fl.tmp);
  fclose(fl.out);
  fclose(fl.col);
  fclose(fl.ras);
  fclose(fl.fri);
  fclose(fl.zmp);
  fclose(fl.his);
}

/* This function reads and processes the first input line */
void read_first_input_line(scan_val *sval, fl_st fl)
{
  int ipar;
  char line[Mline], *p1, par2all[Mword], *p2, *p3;

  if (fgets(line, Mline, fl.in) == NULL)
  {
    printf("empty input file!!!\n");
    exit(0);
  }
  
  p1 = strtok(line, " ");
  sscanf(p1, " %s", sval->par1);
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2all);

  if (strstr(par2all, "#") == NULL)
  {
    strcpy(sval->par2, par2all);
    sval->npt = 1;
  }
  else
  {
    p3 = &sval->par2[0];
    for (p2 = &par2all[0]; *p2 != '#'; p2++) *p3++ = *p2;
    *p3 = '\0';
    p2++;
   sscanf(p2, "%d", &sval->npt);
  }

  fprintf(fl.out, "par1=%s par2=%s npt=%d", sval->par1, sval->par2,
  sval->npt);

  if (sval->scan_type == 'e')
  {
    p1 = strtok(NULL, " ");
    sscanf(p1, "parmin=%lf", &sval->parmin);
    p1 = strtok(NULL, " ");
    sscanf(p1, "parmax=%lf", &sval->parmax);
    p1 = strtok(NULL, " ");
    sscanf(p1, "npar=%d", &sval->npar);
    p1 = strtok(NULL, " ");
    sscanf(p1, "nrepeat=%d", &sval->nrepeat);
    fprintf(fl.out, " parmin=%lf parmax=%lf\nnpar=%d nrepeat=%d\n", 
    sval->parmin, sval->parmax, sval->npar, sval->nrepeat);

    if (sval->npar == 0)
    {
      sval->par_ar[0] = sval->parmin;
    }
    else
    {
      for (ipar=0; ipar<=sval->npar; ipar++)
      {
        sval->par_ar[ipar] = sval->parmin + (sval->parmax - sval->parmin) *
        ipar / sval->npar;
      } 
    }
  }
  else if (sval->scan_type == 'u')
  {
    ipar=-1;
    while((p1 = strtok(NULL, " ")) != NULL)
    {
      sscanf(p1, "%lf", &sval->par_ar[++ipar]);
    }
    sval->npar=ipar;
  }
  else
  {
    printf("wrong scan_type!!!\n");
    exit(0);
  }

  for (ipar=0; ipar<=sval->npar; ipar++)
    fprintf(fl.out, "ipar=%d par_ar=%lf\n", ipar, sval->par_ar[ipar]);

  return;
}

/* This function updates the input file and write the new parameter   */
/* value(s)                                                           */
void update_file_old(scan_val *sval, int skip_lines, fl_st fl)
{
  int nget, nchange;
  char line[Mline], iline;

  rewind(fl.in);
  rewind(fl.tmp);
  for (iline=1; iline<=skip_lines; iline++)
  {
    if (fgets(line, Mline, fl.in) == NULL)
    {
      printf("no line to read!\n");
      exit(0);
    }
  }

  /* no scanning */
  if (sval->scan_type == 'n')
  {
    while (fgets(line, Mline, fl.in) != NULL) fputs(line, fl.tmp);
    rewind(fl.tmp);
    return;
  }

  /* scanning - multiplying all the occurrences of the specific parameter */
  /* value                                                                */

  if (strcmp(sval->par1, "ALL") == 0)
  {
    nchange=0;
    while (nget = (fgets(line, Mline, fl.in)) != NULL)
    {
      if (process_line_old(sval, line, fl.tmp)) nchange++;
    }
    rewind(fl.tmp);
    return;
  }

  /* scanning - changing the specific parameter value */

  while (nget = (fgets(line, Mline, fl.in)) != NULL)
  {
    if (strncmp(line, sval->par1, strlen(sval->par1)) != 0)
    {
      fputs(line, fl.tmp);
    }
    else
      break;
  }
  fputs(line, fl.tmp);
  
  if (nget == 0) 
  {
    printf("par1=%s not found!!!\n", sval->par1);
    exit(0);
  }

  while (nget = (fgets(line, Mline, fl.in)) != NULL)
  {
    if (process_line_old(sval, line, fl.tmp))
    {
       break;
    }
  }

  /* checking for end of file */
  if (fgets(line, Mline, fl.in) != NULL) 
    fputs(line, fl.tmp);
  else
  {
    printf("match not found!!!\n");
   exit(0);
  }

  while (fgets(line, Mline, fl.in) != NULL) fputs(line, fl.tmp);
  rewind(fl.tmp);
}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2                                           */
int process_line_old(scan_val *sval, char line[], FILE *ftmp)
{
  int ret_val;

  if (strstr(sval->par2, ":") == NULL)
  {
    /* printf("par2=%s=NU\n", sval->par2); */
    ret_val = process_line_no_colon_old(sval, line, ftmp);
  }
  else
  {
    /* printf("par2=%s!=NU\n", sval->par2); */
    ret_val = process_line_yes_colon_old(sval, line, ftmp);
  }

  return(ret_val);
}


/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2 .                                         */
/* There is no colon (:) in sval->par2 .                                 */
int process_line_no_colon_old(scan_val *sval, char line[], FILE *ftmp)
{
  double par_ref;
  int  il, il_end, im, ipt, iw;
  int cond1, cond2, cond3;
  char *pline, word[Mword];

  il_end = -1;
  while (line[il_end+1] != '\n' && il_end < Mline-2) il_end++;
  pline = line;
  il = -1;
  
  while (++il <= il_end)
  {
    /* Condition: matched pattern, '=' at the end, ' ' or beginning of  */
    /* line at the beginning                                            */
    cond1 = strncmp(pline+il, sval->par2, strlen(sval->par2)) == 0;
    cond2 = line[il + strlen(sval->par2)] == '=';
    if (il == 0)
      cond3 = 1;
    else
      cond3 = line[il - 1] == ' ';
    if (cond1 && cond2 && cond3) break;
  }

  if (il >= il_end-1)
  /* par2 does not appear in line */
  {
    fputs(line, ftmp);
    return(0);
  }
  else
  /* par2 appears in line */
  {
    for (im=0; im<il; im++) fputc(line[im], ftmp);
    fprintf(ftmp, "%s=", sval->par2);

    while (line[il-1] != '=') il++;
    ipt=0;
    while (ipt < sval->npt-1)
    {
      putc(line[++il], ftmp);
      if (line[il-1] != ' ' && line[il] == ' ') ipt++;
    }

    while (line[il] == ' ') il++;
    iw=-1;
    while ((line[il] != ' ') && (line[il] != '\n'))
    {
      word[++iw] = line[il];
      il++;
    }
    word[++iw] = '\0';

    if (strcmp(sval->par1, "ALL") != 0)
      fprintf(ftmp, "%lf", sval->par_ar[sval->ipar]);
    else
    {
      sscanf(word, "%lf", &par_ref);
      fprintf(ftmp, "%lf", sval->par_ar[sval->ipar] * par_ref);
    }

    for (im=il; im<=il_end; im++) fputc(line[im], ftmp);
    fputc('\n', ftmp);

    return(1);
  }

}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2 .                                         */
/* There is a colon (:) in sval->par2 .                                  */
int process_line_yes_colon_old(scan_val *sval, char line[], FILE *ftmp)
{
  double par_ref;
  int  il, il_end, im, ipt, iw;
  int cond1, cond2, cond3;
  int len2an;
  char line_cp[Mword], par2a[Mword], par2an[Mword], par2b[Mword];
  char *pline, word[Mword], *p1;

  il_end = -1;
  while (line[il_end+1] != '\n' && il_end < Mline-2) il_end++;
  pline = line;

  strcpy(line_cp, sval->par2);
  p1 = strtok(line_cp, ":");
  sscanf(p1, " %s", par2a);
  strcat(strcpy(par2an, par2a),":");
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2b);
  len2an = strlen(par2an);

  /* Checking for par2an */

  if (strncmp(line, par2an, len2an))
  {
    fputs(line, ftmp);
    return(0);
  }

  il = len2an;

  while (++il <= il_end)
  {
    /* Condition: matched pattern, '=' at the end, ' ' or beginning of  */
    /* line at the beginning                                            */
    cond1 = strncmp(pline+il, par2b, strlen(par2b)) == 0;
    cond2 = line[il + strlen(par2b)] == '=';
    if (il == 0)
      cond3 = 1;
    else
      cond3 = line[il - 1] == ' ';
    if (cond1 && cond2 && cond3) break;
  }


  if (il >= il_end-1)
  /* par2 does not appear in line */
  {
    fputs(line, ftmp);
    return(0);
  }
  else
  /* par2 appears in line */
  {
    for (im=0; im<il; im++) fputc(line[im], ftmp);
    fprintf(ftmp, "%s=", par2b);

    while (line[il-1] != '=') il++;
    ipt=0;
    while (ipt < sval->npt-1)
    {
      putc(line[++il], ftmp);
      if (line[il-1] != ' ' && line[il] == ' ') ipt++;
    }

    while (line[il] == ' ') il++;
    iw=-1;
    while (line[il] != ' ')
    {
      word[++iw] = line[il];
      il++;
    }
    word[++iw] = '\0';

    if (strcmp(sval->par1, "ALL") != 0)
      fprintf(ftmp, "%lf", sval->par_ar[sval->ipar]);
    else
    {
      sscanf(word, "%lf", &par_ref);
      fprintf(ftmp, "%lf", sval->par_ar[sval->ipar] * par_ref);
    }

    for (im=il; im<=il_end; im++) fputc(line[im], ftmp);
    fputc('\n', ftmp);
    return(1);
  }
}

void write_avr(par_all *parall, avr_val *av, fl_st fl)
{
  int ipop;
    
  fprintf(fl.avr, "%lf", av->par);
  fprintf(fl.avr, " %d %d", av->ipar, av->irepeat);

  for (ipop=0; ipop<=av->npop; ipop++)
  {
    fprintf(fl.avr, " %lf %lf %lf", av->fr_pop[ipop], av->av_cv[ipop],
    av->sig_cv[ipop]);

    /* av->nfire_av[ipop] */
    fprintf(fl.avr, " %lf %lf", av->Z1md_av[ipop], av->Z1phi_av[ipop]);

    if (ipop >= 1)
    {
      fprintf(fl.avr, " %lf", av->chi[ipop]);
      fprintf(fl.avr, " %lf %lf", av->ratio_no_firing[ipop],
      av->fr_pop_sd[ipop]);
    }

    if (av->pop_name[ipop] == 'P')
    {
      fprintf(fl.avr, " %lf %lf", av->fr_subpop[ipop][1],
      av->fr_subpop[ipop][2]);
    }
  }

  fprintf(fl.avr, "\n");
}


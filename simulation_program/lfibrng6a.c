/*************************************************************************/
/*************************************************************************/
/* A MULTI-PROCESSOR LAGGED-FIBONACCI RANDOM NUMBER GENERATION SYSTEM    */
/*                                                                       */ 
/* Authors: Steven A. Cuccaro and Daniel V. Pryor,                       */
/*            IDA/Supercomputing Research Center (SRC)                   */
/* E-Mail: cuccaro@super.org      pryor@super.org                        */
/*                                                                       */ 
/* Copyright 1993 December 20, United States Government as Represented   */
/* by the Director, National Security Agency. All rights reserved.       */
/*                                                                       */
/* Disclaimer: SRC expressly disclaims any and all warranties, expressed */
/* or implied, concerning the enclosed software. This software was       */
/* developed at SRC for use in internal research. The intent in sharing  */
/* this software is to promote the productive interchange of ideas       */
/* throughout the research community. All software is furnished on an    */
/* "as is" basis. No further updates to this software should be          */
/* expected. Although this may occur, no commitment exists. The authors  */
/* certainly invite your comments as well as the reporting of any bugs.  */
/* SRC cannot commit that any or all bugs will be fixed.                 */
/*************************************************************************/
/*************************************************************************/

/*      This is version 6.3, created 23 January 1995                     */

#include <stdio.h>
#include <math.h>

/*      BITS_IN_INT_GEN is the log_2 of the modulus of the generator     */
/*      BITS_IN_INT_GEN must be = 32 to guarantee that everything works  */
/*      INT_MOD_MASK is used to perform modular arithmetic - specifying  */
/*           this value compensates for different sized words on         */
/*           different architectures                                     */
/*      INT_MASK is used to mask out the part of the generator which     */
/*           is not in the canonical form; it should be                  */
/*           2^{BITS_IN_INT_GEN-1}-1                                     */
/*      MAX_BIT_INT is the largest bit position allowed in the index     */
/*           of the node - it equals BITS_IN_INT_GEN - 2                 */
#define BITS_IN_INT_GEN 32
#if (BITS_IN_INT_GEN==32)
#define INT_MOD_MASK 0xffffffff
#else
#define INT_MOD_MASK ((unsigned)(1<<BITS_IN_INT_GEN)-1)
#endif
#define INT_MASK ((unsigned)(1<<(BITS_IN_INT_GEN-1))-1)
#define MAX_BIT_INT (BITS_IN_INT_GEN-2)

/*      BITS_IN_FLT_GEN is the log_2 of the modulus of the generator     */
/*      BITS_IN_FLT_GEN must be < or = 8*sizeof(int)                     */
/*      FLT_MASK is used to mask out the part of the generator which     */
/*            is not in the canonical form; it should be                 */
/*            2^{BITS_IN_FLT_GEN-1}-1                                    */
/*      MAX_BIT_FLT is the largest bit position allowed in the index     */
/*            of the node - it equals BITS_IN_FLT_GEN - 2                */
#define BITS_IN_FLT_GEN 23
#define FLT_MASK ((1<<(BITS_IN_FLT_GEN-1))-1)
#define MAX_BIT_FLT (BITS_IN_FLT_GEN-2)
 
/*      BITS_IN_DBL_GEN is the log_2 of the modulus of the generator     */
/*           45 is chosen since CRAYs have only 46 bits of precision     */
/*      DBL_MASK is used to mask out the part of the generator which     */
/*           is not in the canonical form; it should be                  */
/*           2^{BITS_IN_DBL_GEN-BITS_IN_INT_GEN}-1                       */
/*	DBL_CONV1 and DBL_CONV2 are used to divide the integer fills     */
/*           generated to convert to double precision                    */
#define BITS_IN_DBL_GEN 45
#define DBL_CONV1 (4.0*(1<<30))
#define DBL_CONV2 (1<<13)
#define DBL_MASK ((1<<13)-1)

/*      RUNUP, FRUNUP and DRUNUP keep certain generators from looking    */
/*           too similar in the first few words output                   */
#define RUNUP (2*BITS_IN_INT_GEN)
#define FRUNUP (2*BITS_IN_FLT_GEN)
#define DRUNUP (2*BITS_IN_DBL_GEN)

/*      GS0 gives a more "random" distribution of generators when the    */
/*      user uses small integers as seeds                                */
#define GS0 0x372f05ac
#define TOOMANY "generator has branched maximum number of times;\nindependence of generators no longer guaranteed"

/*************************************************************************/
/*************************************************************************/
/*                  STRUCTURES AND GLOBAL DATA                           */
/*************************************************************************/
/*************************************************************************/

struct lfgen {
      int *si;      /* sets next branch seed  */
      int *r;       /* pointer to the current fill  */
};
struct lfgenf {
        int *si;    /* sets next branch seed  */
        float *r;   /* pointer to the current fill  */
        int hptr;   /* integer pointer into fill  */
};
struct lfgend {
        int *si;    /* sets next branch seed  */
        double *r;  /* pointer to the current fill  */
        int hptr;   /* integer pointer into fill  */
};

struct vstruct {
      int L;
      int K;
      int LSBS;
};

/*      MAXPOS is maximum position of set bit (see valid[])              */
#define MAXPOS 21 

struct vstruct valid[] = { {2, 1, 3}, {3, 2, 1}, {5, 3, 6},
{7, 4, 2}, {10, 3, 2}, {11, 2, 1}, {15, 4, 2}, {17, 5, 1024},
{29, 2, 1}, {31, 6, 4}, {35, 2, 1}, {55, 24, 2048}, {63, 31, 16384},
{127, 97, 2097152}, {0,0,0} };

int gseed=0,lval=0,kval=0;

/*************************************************************************/
/*************************************************************************/
/*                    ERROR PRINTING FUNCTION                            */
/*************************************************************************/
/*************************************************************************/

#ifdef ANSI_C
int errprint(char *level, char *routine, char *error)
#else
int errprint(level, routine, error)
char *level,*routine,*error;
#endif
{
      printf("%s from %s: %s\n",level,routine,error);
}

/*************************************************************************/
/*************************************************************************/
/*            ROUTINES USED TO CREATE GENERATOR FILLS                    */
/*************************************************************************/
/*************************************************************************/

/**************************/
/* function bitcnt:       */
/**************************/
#ifdef ANSI_C
int bitcnt(int x)
#else
int bitcnt(x)
int x;
#endif
{
      unsigned i=0,y;

      for (y=x; y; y &= (y-1) ) i++;
      return(i);
}

/**************************/
/* function advance_reg:  */
/**************************/
#ifdef ANSI_C
void advance_reg(int *reg_fill)
#else
int advance_reg(reg_fill)
int *reg_fill;
#endif
{
/*      the register steps according to the primitive polynomial         */
/*           (64,4,3,1,0); each call steps register 64 times             */
/*      we use two words to represent the register to allow for integer  */
/*           size of 32 bits                                             */

#ifdef ANSI_C
      const int mask = 0x1b;
      const int adv_64[4][2];
#else
      int mask = 0x1b;
      int adv_64[4][2];
#endif

      int i,new_fill[2];
      unsigned temp;

      adv_64[0][0] = 0xb0000000;
      adv_64[0][1] = 0x1b;
      adv_64[1][0] = 0x60000000;
      adv_64[1][1] = 0x2d;
      adv_64[2][0] = 0xc0000000;
      adv_64[2][1] = 0x5a;
      adv_64[3][0] = 0x80000000;
      adv_64[3][1] = 0xaf;
      new_fill[1] = new_fill[0] = 0;
      temp = mask<<27;
      for (i=27;i>=0;i--) {
            new_fill[0] = (new_fill[0]<<1) | (1&bitcnt(reg_fill[0]&temp));
            new_fill[1] = (new_fill[1]<<1) | (1&bitcnt(reg_fill[1]&temp));
            temp >>= 1;
            }
      for (i=28;i<32;i++) {
            temp = bitcnt(reg_fill[0]&(mask<<i));
            temp ^= bitcnt(reg_fill[1]&(mask>>(32-i)));
            new_fill[0] |= (1&temp)<<i;
            temp = bitcnt(reg_fill[0]&adv_64[i-28][0]);
            temp ^= bitcnt(reg_fill[1]&adv_64[i-28][1]);
            new_fill[1] |= (1&temp)<<i;
            }
      reg_fill[0] = new_fill[0];
      reg_fill[1] = new_fill[1];
}

/**************************/
/* functions get_fill_?:  */
/**************************/

#ifdef ANSI_C
int get_fill_int(const int *n, int *r0)
#else
int get_fill_int(n,r0)
int *n,*r0;
#endif
{
      int i,j,temp[2];

/*      initialize the shift register with the node number XORed with    */
/*           the global seed                                             */
/*      fill the shift register with two copies of this number           */
/*           except when equal to zero                                   */
      temp[1] = temp[0] = n[0]^gseed;
      if (!temp[0]) temp[0] = GS0;

/*      advance the shift register some                                  */
      advance_reg(temp);
      advance_reg(temp);

/*      the generator is filled with the lower 32 bits of the shift      */
/*           register at each time                                       */
/*      the node number is XORed into the fill and the canonical form    */
/*           for the LSB is instituted                                   */
      r0[0] = (INT_MASK&n[0])<<1;
      for (i=1;i<lval-1;i++) {
            advance_reg(temp);
            r0[i] = (unsigned)(INT_MASK&(temp[0]^n[i]))<<1;
            }
      r0[lval-1] = 0;
      i = -1;
      while (valid[++i].L) if (valid[i].L==lval) break;
      temp[0] = (MAXPOS<lval-1) ? MAXPOS : lval-1;
      for (j=0;j<=temp[0];j++) r0[j] |= (valid[i].LSBS>>j)&1;
}

#ifdef ANSI_C
int get_fill_flt(const int *n, float *r0)
#else
int get_fill_flt(n,r0)
int *n;
float *r0;
#endif
{
        int i,j,temp[2],*arr;
 
/*      initialize the shift register with the node number XORed with    */
/*           the global seed                                             */
/*      fill the shift register with two copies of this number           */
/*           except when equal to zero                                   */
        arr = (int *) malloc(lval*sizeof(int));
        temp[1] = temp[0] = n[0]^gseed;
        if (!temp[0]) temp[0] = GS0;
 
/*      advance the shift register some                                  */
        advance_reg(temp);
        advance_reg(temp);
 
/*      the generator is filled with the lower BITS_IN_FLT_GEN bits of   */
/*           the shift register at each time                             */
/*      the node number is XORed into the fill and the canonical form    */
/*           for the LSB is instituted                                   */
        arr[0] = (FLT_MASK&n[0])<<1;
        for (i=1;i<lval-1;i++) {
                advance_reg(temp);
                arr[i] = (FLT_MASK&(temp[0]^n[i]))<<1;
                }
        arr[lval-1] = 0;
        i = -1;
        while (valid[++i].L)
                if (valid[i].L==lval) break;
        if (!valid[i].L) return(1);
        temp[0] = (MAXPOS<lval-1) ? MAXPOS : lval-1;
        for (j=0;j<=temp[0];j++) arr[j] |= (valid[i].LSBS>>j)&1;
        for (i=0;i<lval;i++) r0[i] = arr[i]/(float)(1<<BITS_IN_FLT_GEN);
        free(arr);
        return(0);
}

#ifdef ANSI_C
int get_fill_dbl(const int *n, double *r0)
#else
int get_fill_dbl(n,r0)
int *n;
double *r0;
#endif
{
        int i,j,k,temp[2];
        unsigned *arr;

/*      initialize the shift register with the node number XORed with    */
/*           the global seed                                             */
/*      fill the shift register with two copies of this number           */
/*           except when equal to zero                                   */
        arr = (unsigned *) malloc(lval*sizeof(unsigned));
        temp[1] = temp[0] = n[0]^gseed;
        if (!temp[0]) temp[0] = GS0;

/*      advance the shift register some                                  */
        advance_reg(temp);
        advance_reg(temp);

/*      the generator is filled with the lower BITS_IN_INT_GEN bits of   */
/*           the shift register at this time                             */
/*      the node number is XORed into the fill and the canonical form    */
/*           for the LSB is instituted                                   */
        arr[0] = (INT_MASK&n[0])<<1;
        for (i=1;i<lval-1;i++) {
                advance_reg(temp);
                arr[i] = (INT_MASK&(temp[0]^n[i]))<<1;
                }
        arr[lval-1] = 0;
        i = -1; 
        while (valid[++i].L)
                if (valid[i].L==lval) break;
        if (!valid[i].L) return(1);
        k = (MAXPOS<lval-1) ? MAXPOS : lval-1;
        for (j=0;j<=k;j++) arr[j] |= (valid[i].LSBS>>j)&1;
/*      the generator is converted to double precision, keeping the      */
/*           lowest significant bit at 2^{-BITS_IN_DBL_GEN}              */
/*      the remaining bits of initial fill are added now                 */
        for (i=0;i<lval-1;i++) {
                advance_reg(temp);
                r0[i] = arr[i]/(double)DBL_CONV1 + (temp[0]&DBL_MASK);
                r0[i] /= (double)DBL_CONV2;
                }
        r0[lval-1] = 0.0;
        free(arr);
        return(0);
}


/*************************************************************************/
/*************************************************************************/
/*            SI_DOUBLE: updates index for next spawning                 */
/*************************************************************************/
/*************************************************************************/

#ifdef ANSI_C
void si_double_int(int *a, const int *b)
#else
void si_double_int(a,b)
int *a,*b;
#endif
{
      int i;

      if (b[lval-2]&(1<<MAX_BIT_INT))
            errprint("WARNING","si_double_int",TOOMANY);
      a[lval-2] = INT_MASK&((unsigned)b[lval-2]<<1);
      for (i=lval-3;i>=0;i--) {
            if (b[i]&(1<<MAX_BIT_INT)) a[i+1]++;
            a[i] = INT_MASK&((unsigned)b[i]<<1);
            }
}

#ifdef ANSI_C
void si_double_flt(int *a, int *b)
#else
void si_double_flt(a,b)
int *a,*b;
#endif
{
        int i;
 
        if (b[lval-2]&(1<<MAX_BIT_FLT))
              errprint("WARNING","si_double_flt",TOOMANY);
        a[lval-2] = FLT_MASK&(b[lval-2]<<1);
        for (i=lval-3;i>=0;i--) {
                if (b[i]&(1<<MAX_BIT_FLT)) a[i+1]++;
                a[i] = FLT_MASK&(b[i]<<1);
                }
}

#ifdef ANSI_C
void si_double_dbl(int *a, int *b)
#else
void si_double_dbl(a,b)
int *a,*b;
#endif
{
        int i;
/*      NOTE: generator artificially limited to BITS_IN_INT_GEN*(lval-1) */
/*           branches                                                    */
        if (b[lval-2]&(1<<(BITS_IN_INT_GEN-1)))
              errprint("WARNING","si_double_dbl",TOOMANY);
        a[lval-2] = INT_MOD_MASK&(b[lval-2]<<1);
        for (i=lval-3;i>=0;i--) {
                if (b[i]&(1<<(BITS_IN_INT_GEN-1))) a[i+1]++;
                a[i] = INT_MOD_MASK&(b[i]<<1);
                }
}

/*************************************************************************/
/*************************************************************************/
/*            GET_RN: returns generated random number                    */
/*************************************************************************/
/*************************************************************************/

#ifdef ANSI_C
int get_rn_int(int *genptr)
#else
int get_rn_int(genptr)
int *genptr;
#endif
/*      returns value put into new position                              */
{
      unsigned new_val,lptr;
      int *r0,*hptr;

      r0 = ((struct lfgen *)genptr)->r;
      hptr = ((struct lfgen *)genptr)->si+lval-1;
      lptr = *hptr + kval;
      if (lptr>=lval) lptr -= lval;
/*    INT_MOD_MASK causes arithmetic to be modular when integer size is  */
/*         greater than 32 bits                                          */
      new_val = INT_MOD_MASK&(r0[*hptr] + r0[lptr]);
      r0[*hptr] = new_val;
      if (--*hptr < 0) *hptr = lval - 1;
/*    cast to unsigned before shifting ensures that a 0 is shifted in    */
      return ((unsigned int)new_val >>1);
}

#ifdef ANSI_C
float get_rn_flt(int *genptr)
#else
float get_rn_flt(genptr)
int *genptr;
#endif
/*      returns value put into new position                              */
{
        unsigned lptr;
        int *hptr = &((struct lfgenf *)genptr)->hptr;
        float *r0,new_val;
 
        r0 = ((struct lfgenf *)genptr)->r;
        lptr = *hptr + kval;
        if (lptr>=lval) lptr -= lval;
        new_val = r0[*hptr] + r0[lptr];
        if (new_val >= 1.0) new_val -= 1.0;
        r0[*hptr] = new_val;
        if (--*hptr < 0) *hptr = lval - 1;
        return(new_val);
} 

#ifdef ANSI_C
double get_rn_dbl(int *genptr)
#else
double get_rn_dbl(genptr)
int *genptr;
#endif
/*      returns value put into new position                              */
{
        unsigned lptr;
        int *hptr = &((struct lfgend *)genptr)->hptr;
        double *r0,new_val;

        r0 = ((struct lfgend *)genptr)->r;
        lptr = *hptr + kval;
        if (lptr>=lval) lptr -= lval;
        new_val = r0[*hptr] + r0[lptr];
        if (new_val >= 1.0) new_val -= 1.0;
        r0[*hptr] = new_val;
        if (--*hptr < 0) *hptr = lval - 1;
        return(new_val);
}

/*************************************************************************/
/*************************************************************************/
/*            INITIALIZE: starts the whole thing going                   */
/*************************************************************************/
/*************************************************************************/

#ifdef ANSI_C
int **initialize_int(int ngen, int length, int seed, const int *nstart)
#else
int **initialize_int(ngen,length,seed,nstart)
int ngen,length,seed,*nstart;
#endif
{
      int i,j,k,l,*nindex,*order;
      struct lfgen **q;

/*      allocate memory for node number and fill of each generator       */
      order = (int *)malloc(ngen*sizeof(int));
      q = (struct lfgen **) malloc(ngen*sizeof(struct lfgen *));
      if (q == NULL) return((int **)NULL);
      for (i=0;i<ngen;i++) {
            q[i] = (struct lfgen *) malloc(sizeof(struct lfgen));
            q[i]->si = (int *) malloc(length*sizeof(int));
            q[i]->r = (int *) malloc(length*sizeof(int));
            if (q[i]->r == NULL) return((int **)NULL);
            q[i]->si[length-1] = length - 1;
            }
/*      specify register fills and node number arrays                    */
/*      do fills in tree fashion so that all fills branch from index     */
/*           contained in nstart array                                   */
      get_fill_int(nstart,q[0]->r);
      si_double_int(q[0]->si,nstart);
      q[0]->si[0]++;
      i = 1;
      order[0] = 0;
      if (ngen>1) while (1) {
            l = i;
            for (k=0;k<l;k++) {
                  nindex = q[order[k]]->si;
                  get_fill_int(nindex,q[i]->r);
                  si_double_int(nindex,nindex);
                  for (j=0;j<length-1;j++) q[i]->si[j] = nindex[j];
                  q[i]->si[0]++;
                  if (ngen == ++i) break;
                  }
            if (ngen == i) break;
            for (k=l-1;k>0;k--) {
                  order[2*k+1] = l+k;
                  order[2*k] = order[k];
                  }
            order[1] = l;
            }
      for (i=ngen-1;i>=0;i--) {
            k = 0;
            for (j=1;j<lval-1;j++) if (q[i]->si[j]) k = 1;
            if (!k) break;
            for (j=0;j<length*RUNUP;j++) get_rn_int((int *)q[i]);
            }
      while (i>=0) {
            for (j=0;j<4*length;j++) get_rn_int((int *)q[i]);
            i--;
            }
      return((int **)q);
}

#ifdef ANSI_C
int **initialize_flt(int ngen, int length, int seed, int *nstart)
#else
int **initialize_flt(ngen,length,seed,nstart)
int ngen,length,seed,*nstart;
#endif
{
        int doexit=0,i,j,k,l,*nindex,*order;
        struct lfgenf **q;

/*      allocate memory for node number and fill of each generator       */
        order = (int *)malloc(ngen*sizeof(int));
        q = (struct lfgenf **) malloc(ngen*sizeof(struct lfgenf *));
        if (q == NULL) return((int **)NULL);
        for (i=0;i<ngen;i++) {
                q[i] = (struct lfgenf *) malloc(sizeof(struct lfgenf));
                q[i]->si = (int *) malloc((length-1)*sizeof(int));
                q[i]->r = (float *) malloc(length*sizeof(float));
                if (q[i]->r == NULL) return((int **)NULL);
                q[i]->hptr = length - 1;
                }
/*      specify register fills and node number arrays                    */
/*      do fills in tree fashion so that all fills branch from index     */
/*           contained in nstart array                                   */
        get_fill_flt(nstart,q[0]->r);
        si_double_flt(q[0]->si,nstart);
        q[0]->si[0]++;
        i = 1;
        order[0] = 0;
        if (ngen>1) while (1) {
                l = i;
                for (k=0;k<l;k++) {
                        nindex = q[order[k]]->si;
                        get_fill_flt(nindex,q[i]->r);
                        si_double_flt(nindex,nindex);
                        for (j=0;j<length-1;j++) q[i]->si[j] = nindex[j];
                        q[i]->si[0]++;
                        if (ngen == ++i) break;
                        }
                if (ngen == i) break;
                for (k=l-1;k>0;k--) {
                        order[2*k+1] = l+k;
                        order[2*k] = order[k];
                        }
                order[1] = l;
                }
        for (i=ngen-1;i>=0;i--) {
                k = 0;
                for (j=1;j<lval-1;j++) if (q[i]->si[j]) k = 1;
                if (!k) break;
                for (j=0;j<length*FRUNUP;j++) get_rn_flt((int *)q[i]);
                }
        while (i>=0) {
                for (j=0;j<4*length;j++) get_rn_flt((int *)q[i]);
                i--;
                }   
        return((int **)q);
}

#ifdef ANSI_C
int **initialize_dbl(int ngen, int length, int seed, int *nstart)
#else
int **initialize_dbl(ngen,length,seed,nstart)
int ngen,length,seed,*nstart;
#endif
{
        int doexit=0,i,j,k,l,*nindex,*order;
        struct lfgend **q;
 
/*      allocate memory for node number and fill of each generator       */
        order = (int *)malloc(ngen*sizeof(int));
        q = (struct lfgend **) malloc(ngen*sizeof(struct lfgend *));
        if (q == NULL) return((int **)NULL);
        for (i=0;i<ngen;i++) {
                q[i] = (struct lfgend *) malloc(sizeof(struct lfgend));
                q[i]->si = (int *) malloc((length-1)*sizeof(int));
                q[i]->r = (double *) malloc(length*sizeof(double));
                if (q[i]->r == NULL) return((int **)NULL);
                q[i]->hptr = length - 1;
                }
/*      specify register fills and node number arrays                    */
/*      do fills in tree fashion so that all fills branch from index     */
/*           contained in nstart array                                   */
        get_fill_dbl(nstart,q[0]->r);
        si_double_dbl(q[0]->si,nstart);
        q[0]->si[0]++;
        i = 1;
        order[0] = 0;
        if (ngen>1) while (1) {
                l = i;
                for (k=0;k<l;k++) {
                        nindex = q[order[k]]->si;
                        get_fill_dbl(nindex,q[i]->r);
                        si_double_dbl(nindex,nindex);
                        for (j=0;j<length-1;j++) q[i]->si[j] = nindex[j];
                        q[i]->si[0]++;
                        if (ngen == ++i) break;
                        }
                if (ngen == i) break;
                for (k=l-1;k>0;k--) {
                        order[2*k+1] = l+k;
                        order[2*k] = order[k];
                        }
                order[1] = l;
                }
        for (i=ngen-1;i>=0;i--) {
                k = 0;
                for (j=1;j<lval-1;j++) if (q[i]->si[j]) k = 1;
                if (!k) break;
                for (j=0;j<length*DRUNUP;j++) get_rn_dbl((int *)q[i]);
                }
        while (i>=0) {
                for (j=0;j<4*length;j++) get_rn_dbl((int *)q[i]);
                i--;
                }   
        return((int **)q);
}

/*************************************************************************/
/*************************************************************************/
/*            INIT_RNG's: user interface to start things off             */
/*************************************************************************/
/*************************************************************************/

#ifdef ANSI_C
int *init_rng_d_int(int gennum, int length, int seed, int total_gen)
#else
int *init_rng_d_int(gennum,length,seed,total_gen)
int gennum,length,seed,total_gen;
#endif
{
      int doexit=0,i,k,k1;
      int **p=NULL,*nstart=NULL,*si;

/*      gives back one generator (node gennum) with updated spawning     */
/*      info; should be called total_gen times, with different value     */
/*      of gennum in [0,total_gen) each call                             */

/*      check values of gennum and total_gen                             */
      if (total_gen <= 0) {
            total_gen = 1;
            errprint("NOTICE","init_rng_d_int","default value of 1 used for total_gen");
            }
      if (gennum < 0 || gennum >= total_gen) {
            errprint("ERROR","init_rng_d_int","gennum out of range"); 
            return((int *)NULL);
            }
/*      check whether generators have previously been defined            */
/*      guard against access while defining generator parameters for     */
/*            the 1st time                                               */
      if (!lval) {
/*            determine generator to be used                             */
            i = -1;
            while (valid[++i].L) if (length == valid[i].L) break;
            if (!(k=k1=valid[i].K)) {
                  length = 17;
                  k = 5;
                  }
/*        define parameters of & allocate structure for generator        */
            lval = length;
            kval = k;
            gseed = seed^GS0;
            }
      else {
/*      check values of parameters for consistency                       */
            if (!length) length = lval;
            else if (lval!=length) doexit++;
            if (seed && seed!=(gseed^GS0)) doexit += 2;
            }
      if (!k1) errprint("NOTICE","init_rng_d_int","default value of 17 used for length");
      if (doexit) {
            if (doexit&1) errprint("ERROR","init_rng_d_int","changing global L value!");
            if (doexit&2) errprint("ERROR","init_rng_d_int","changing global seed value!");
            return((int *)NULL);
            }
/*      define the starting vector for the initial node                  */
      nstart = (int *)malloc((length-1)*sizeof(int));
      if (nstart == NULL) {
            errprint("ERROR","init_rng_d_int","insufficient memory");
            return((int *)NULL);
            }
      nstart[0] = gennum;
      for (i=1;i<length-1;i++) nstart[i] = 0;
      p = initialize_int(1,lval,gseed,nstart);
      if (p==NULL) {
            errprint("ERROR","init_rng_d_int","insufficient memory");
            return((int *)NULL);
            }
      si = ((struct lfgen *)(p[0]))->si;
      while (si[0] < total_gen && !si[1]) si_double_int(si,si);
      free(nstart);
      return(*p);
}

#ifdef ANSI_C
int *init_rng_d_flt(int gennum, int length, int seed, int total_gen)
#else
int *init_rng_d_flt(gennum,length,seed,total_gen)
int gennum,length,seed,total_gen;
#endif
{
        int doexit=0,i,k,k1;
        int **p=NULL,*nstart=NULL,*si;
 
/*      gives back one generator (node gennum) with updated              */
/*      spawning info; should be called total_gen times, with            */
/*      different value of gennum in [0,total_gen) each call             */

/*      check values of gennum and total_gen                             */
        if (total_gen <= 0) {
                total_gen = 1;
                errprint("NOTICE","init_rng_d_flt","default value of 1 used for total_gen\n");
                }
        if (gennum < 0 || gennum >= total_gen) {
                errprint("ERROR","init_rng_d_flt","gennum out of range \n");
                return((int *)NULL);
                }
/*      check whether generators have previously been defined            */
/*      guard against access while defining generator parameters for     */
/*            the 1st time                                               */
        if (!lval) {
/*              determine generator to be used                           */
                i = -1;
                while (valid[++i].L) if (length == valid[i].L) break;
                if (!(k=k1=valid[i].K)) {
                        length = 17;
                        k = 5;
                        }
/*        define parameters of & allocate structure for generator        */
                lval = length;
                kval = k;
                gseed = seed^GS0;
                }
        else {
/*      check values of parameters for consistency                       */
                if (!length) length = lval;
                else if (lval!=length) doexit++;
                if (seed && seed!=gseed^GS0) doexit += 2;
                }
        if (!k1) errprint("NOTICE","init_rng_d_flt","default value of 17 used for length");
        if (doexit) {
                if (doexit&1) errprint("ERROR","init_rng_d_flt","changing global L value!");
                if (doexit&2) errprint("ERROR","init_rng_d_flt","changing global seed value!");
                return((int *)NULL);
                }
/*      define the starting vector for the initial node                  */
        nstart = (int *)malloc((length-1)*sizeof(int));
        if (nstart == NULL) {
                errprint("ERROR","init_rng_d_flt","insufficient memory");
                return((int *)NULL);
                }
        nstart[0] = gennum;
        for (i=1;i<length-1;i++) nstart[i] = 0;
        p = initialize_flt(1,lval,gseed,nstart);
        if (p==NULL) {
                errprint("ERROR","init_rng_d_flt","insufficient memory");
                return((int *)NULL);
                }
        si = ((struct lfgenf *)(p[0]))->si;
        while (si[0] < total_gen && !si[1]) si_double_flt(si,si);
        free(nstart);
        return(*p);
}

#ifdef ANSI_C
int *init_rng_d_dbl(int gennum, int length, int seed, int total_gen)
#else
int *init_rng_d_dbl(gennum,length,seed,total_gen)
int gennum,length,seed,total_gen;
#endif
{
        int doexit=0,i,k,k1;
        int **p=NULL,*nstart=NULL,*si;

/*      gives back one generator (node gennum) with updated              */
/*      spawning info; should be called total_gen times, with            */
/*      different value of gennum in [0,total_gen) each call             */

/*      check values of gennum and total_gen                             */
        if (total_gen <= 0) {
                total_gen = 1;
                errprint("NOTICE","init_rng_d_dbl","default value of 1 used for total_gen\n");
                }
        if (gennum < 0 || gennum >= total_gen) {
                errprint("ERROR","init_rng_d_dbl","gennum out of range \n");
                return((int *)NULL);
                }
/*      check whether generators have previously been defined            */
/*      guard against access while defining generator parameters for     */
/*            the 1st time                                               */
        if (!lval) {
/*              determine generator to be used                           */
                i = -1;
                while (valid[++i].L) if (length == valid[i].L) break;
                if (!(k=k1=valid[i].K)) {
                        length = 17;
                        k = 5;
                        }
/*        define parameters of & allocate structure for generator        */
                lval = length;
                kval = k;
                gseed = seed^GS0;
                }
        else {
/*      check values of parameters for consistency                       */
                if (!length) length = lval;
                else if (lval!=length) doexit++;
                if (seed && seed!=gseed^GS0) doexit += 2;
                }
        if (!k1) errprint("NOTICE","init_rng_d_dbl","default value of 17 used for length");
        if (doexit) {
                if (doexit&1) errprint("ERROR","init_rng_d_dbl","changing global L value!");
                if (doexit&2) errprint("ERROR","init_rng_d_dbl","changing global seed value!");
                return((int *)NULL);
                }
/*      define the starting vector for the initial node                  */
        nstart = (int *)malloc((length-1)*sizeof(int));
        if (nstart == NULL) {
                errprint("ERROR","init_rng_d_dbl","insufficient memory");
                return((int *)NULL);
                }
        nstart[0] = gennum;
        for (i=1;i<length-1;i++) nstart[i] = 0;
        p = initialize_dbl(1,lval,gseed,nstart);
        if (p==NULL) {
                errprint("ERROR","init_rng_d_dbl","insufficient memory");
                return((int *)NULL);
                }
        si = ((struct lfgend *)(p[0]))->si;
        while (si[0] < total_gen && !si[1]) si_double_dbl(si,si);
        free(nstart);
        return(*p);
}

#ifdef ANSI_C
int **init_rng_s_int(int ngen, int length, int seed)
#else
int **init_rng_s_int(ngen,length,seed)
int ngen,length,seed;
#endif
{
      int doexit=0,i,k,k1;
      int **p=NULL,*nstart=NULL;

/*      gives back vector of generators with nodes [0,ngen)              */
/*      should be called once only                                       */

/*      check value of ngen                                              */
      if (ngen <= 0 ) {
            ngen = 1;
            errprint("NOTICE","init_rng_s_int","default value of 1 used for ngen");
            }
/*      check whether generators have previously been defined            */
/*      guard against access while defining generator parameters for     */
/*            the 1st time                                               */
      if (!lval) {
/*            determine generator to be used                             */
            i = -1;
            while (valid[++i].L) if (length == valid[i].L) break;
            if (!(k=k1=valid[i].K)) {
                  length = 17;
                  k = 5;
                  }
/*        define parameters of & allocate structure for generator        */
            lval = length;
            kval = k;
            gseed = seed^GS0;
            }
      else {
/*      check values of parameters for consistency                       */
            if (!length) length = lval;
            else if (lval!=length) doexit++;
            if (seed && seed!=(gseed^GS0)) doexit += 2;
            }
      if (!k1) errprint("NOTICE","init_rng_s_int","default value of 17 used for length");
      if (doexit) {
            if (doexit&1) errprint("ERROR","init_rng_s_int","changing global L value!");
            if (doexit&2) errprint("ERROR","init_rng_s_int","changing global seed value!");
            return((int **)NULL);
            }
/*      define the starting vector(s) for the initial node(s)            */
      nstart = (int *)malloc((length-1)*sizeof(int));
      if (nstart == NULL) {
            errprint("ERROR","init_rng_s_int","insufficient memory");
            return((int **)NULL);
            }
      for (i=0;i<length-1;i++) nstart[i] = 0;
      p = initialize_int(ngen,lval,gseed,nstart);
      free(nstart);
      if (p == NULL) errprint("ERROR","init_rng_s_int","insufficient memory");
      return(p);
}

#ifdef ANSI_C
int **init_rng_s_flt(int ngen, int length, int seed)
#else
int **init_rng_s_flt(ngen,length,seed)
int ngen,length,seed;
#endif
{
        int doexit=0,i,k,k1;
        int **p=NULL,*nstart=NULL;

/*      gives back vector of generators with nodes [0,ngen)              */
/*      should be called once only                                       */

/*      check value of ngen                                              */
        if (ngen <= 0 ) {
                ngen = 1;
                errprint("NOTICE","init_rng_s_flt","default value of 1 used for ngen");
                }
/*      check whether generators have previously been defined            */
/*      guard against access while defining generator parameters for     */
/*            the 1st time                                               */
          if (!lval) {
/*              determine generator to be used                           */
                i = -1;
                while (valid[++i].L) if (length == valid[i].L) break;
                if (!(k = k1= valid[i].K)) {
                        length = 17;
                        k = 5;
                        }
/*              define parameters of & allocate structure for generator  */
                lval = length;
                kval = k;
                gseed = seed^GS0;
                }
          else {
/*              check values of parameters for consistency               */
                if (!length) length = lval;
                else if (lval!=length) doexit++;
                if (seed && seed!=gseed^GS0) doexit += 2;
                }
        if (!k1) errprint("NOTICE","init_rng_s_flt","default value of 17 used for length");
        if (doexit) {
                if (doexit&1) errprint("ERROR","init_rng_s_flt","changing global L value!");
                if (doexit&2) errprint("ERROR","init_rng_s_flt","changing global seed value!");
                return((int **)NULL);
                }
/*      define the starting vector(s) for the initial node(s)            */
        nstart = (int *)malloc((length-1)*sizeof(int));
        if (nstart == NULL) {
                errprint("ERROR","init_rng_s_flt","insufficient memory");
                return((int **)NULL);
                }
        for (i=0;i<length-1;i++) nstart[i] = 0;
        p = initialize_flt(ngen,lval,gseed,nstart);
        free(nstart);
        if (p == NULL) errprint("ERROR","init_rng_s_flt","insufficient memory");
        return(p);
}

#ifdef ANSI_C
int **init_rng_s_dbl(int ngen, int length, int seed)
#else
int **init_rng_s_dbl(ngen,length,seed)
int ngen,length,seed;
#endif
{
        int doexit=0,i,k,k1;
        int **p=NULL,*nstart=NULL;

/*      gives back vector of generators with nodes [0,ngen)              */
/*      should be called once only                                       */

/*      check value of ngen                                              */
        if (ngen <= 0 ) {
                ngen = 1;
                errprint("NOTICE","init_rng_s_dbl","default value of 1 used for ngen");
                }
/*      check whether generators have previously been defined            */
/*      guard against access while defining generator parameters for     */
/*            the 1st time                                               */
          if (!lval) {
/*              determine generator to be used                           */
                i = -1;
                while (valid[++i].L) if (length == valid[i].L) break;
                if (!(k = k1= valid[i].K)) {
                        length = 17;
                        k = 5;
                        }
/*              define parameters of & allocate structure for generator  */
                lval = length;
                kval = k;
                gseed = seed^GS0;
                }
          else {
/*              check values of parameters for consistency               */
                if (!length) length = lval;
                else if (lval!=length) doexit++;
                if (seed && seed!=gseed^GS0) doexit += 2;
                }
        if (!k1) errprint("NOTICE","init_rng_s_dbl","default value of 17 used for length");
        if (doexit) {
                if (doexit&1) errprint("ERROR","init_rng_s_dbl","changing global L value!");
                if (doexit&2) errprint("ERROR","init_rng_s_dbl","changing global seed value!");
                return((int **)NULL);
                }
/*      define the starting vector(s) for the initial node(s)            */
        nstart = (int *)malloc((length-1)*sizeof(int));
        if (nstart == NULL) {
                errprint("ERROR","init_rng_s_dbl","insufficient memory");
                return((int **)NULL);
                }
        for (i=0;i<length-1;i++) nstart[i] = 0;
        p = initialize_dbl(ngen,lval,gseed,nstart);
        free(nstart);
        if (p == NULL) errprint("ERROR","init_rng_s_dbl","insufficient memory");        return(p);
}

/*************************************************************************/
/*************************************************************************/
/*                  SPAWN_RNG: spawns new generators                     */
/*************************************************************************/
/*************************************************************************/

#ifdef ANSI_C
int **spawn_rng_int(int *genptr, int nspawned)
#else
int **spawn_rng_int(genptr,nspawned)
int *genptr,nspawned;
#endif
{
      int i,*p,**q=NULL;

      if (nspawned > 0) {
            p = ((struct lfgen *)genptr)->si;
            q = initialize_int(nspawned,lval,gseed,p);
            si_double_int(p,p);
            if (q == NULL)
                  errprint("ERROR","spawn_rng_int","insufficient memory");
            }
      return(q);
}

#ifdef ANSI_C
int **spawn_rng_flt(int *genptr, int nspawned)
#else   
int **spawn_rng_flt(genptr,nspawned)
int *genptr,nspawned;
#endif  
{
        int i,*p,**q=NULL;

        if (nspawned > 0) {
                p = ((struct lfgenf *)genptr)->si;
                q = initialize_flt(nspawned,lval,gseed,p);
                si_double_flt(p,p);
                if (q == NULL) errprint("ERROR","spawn_rng_flt","insufficient memory");
                }
        return(q);
}

#ifdef ANSI_C
int **spawn_rng_dbl(int *genptr, int nspawned)
#else
int **spawn_rng_dbl(genptr,nspawned)
int *genptr,nspawned;
#endif
{
        int i,*p,**q=NULL;

        if (nspawned > 0) {
                p = ((struct lfgend *)genptr)->si;
                q = initialize_dbl(nspawned,lval,gseed,p);
                si_double_dbl(p,p);
                if (q == NULL) errprint("ERROR","spawn_rng_dbl","insufficient memory");
                }
        return(q);
}


/*************************************************************************/
/*************************************************************************/
/*                  UTILITY ROUTINES                                     */
/*************************************************************************/
/*************************************************************************/

#ifdef ANSI_C
int get_llag_rng(void)
#else
int get_llag_rng()
#endif
{
        return(lval);
}

#ifdef ANSI_C
int get_klag_rng(void)
#else
int get_klag_rng()
#endif
{
        return(kval);
}

#ifdef ANSI_C
int get_seed_rng(void)
#else
int get_seed_rng()
#endif
{
        return(GS0^gseed);
}

#ifdef ANSI_C
int get_hptr_rng_int(const int *genptr)
#else
int get_hptr_rng_int(genptr)
int *genptr;
#endif
{
        return (((struct lfgen *)genptr)->si[lval-1]);
}

#ifdef ANSI_C
int get_hptr_rng_flt(const int *genptr)
#else
int get_hptr_rng_flt(genptr)
int *genptr;
#endif
{          
          return (((struct lfgenf *)genptr)->hptr);
}
 
#ifdef ANSI_C
int get_hptr_rng_dbl(const int *genptr)
#else
int get_hptr_rng_dbl(genptr)
int *genptr;
#endif
{
          return (((struct lfgend *)genptr)->hptr);
}

#ifdef ANSI_C
int *get_fill_rng_int(const int *genptr)
#else
int *get_fill_rng_int(genptr)
int *genptr;
#endif
{
      int i,*p;
      int *pp;

      pp = ((struct lfgen *)genptr)->r;
      p = (int *) malloc(lval*sizeof(int));
      for (i=0;i<lval;i++) p[i] = pp[i];
      return(p);
}

#ifdef ANSI_C
float *get_fill_rng_flt(const int *genptr)
#else
float *get_fill_rng_flt(genptr)
int *genptr;
#endif
{
        int i;
        float *p,*pp;
 
        pp = ((struct lfgenf *)genptr)->r;
        p = (float *) malloc(lval*sizeof(float));
        for (i=0;i<lval;i++) p[i] = pp[i];
        return(p);
}

#ifdef ANSI_C
double *get_fill_rng_dbl(const int *genptr)
#else
double *get_fill_rng_dbl(genptr)
int *genptr;
#endif
{
        int i;
        double *p,*pp;

        pp = ((struct lfgend *)genptr)->r;
        p = (double *) malloc(lval*sizeof(double));
        for (i=0;i<lval;i++) p[i] = pp[i];
        return(p);
}

#ifdef ANSI_C
int *get_next_index_rng_int(const int *genptr)
#else
int *get_next_index_rng_int(genptr)
int *genptr;
#endif
{
      int i,*p,*pp;

      pp = ((struct lfgen *)genptr)->si;
      p = (int *) malloc((lval-1)*sizeof(int));
      for (i=0;i<lval-1;i++) p[i] = pp[i];
      return(p);
}

#ifdef ANSI_C
int *get_next_index_rng_flt(const int *genptr)
#else
int *get_next_index_rng_flt(genptr)
int *genptr;
#endif
{
        int i;
        int *p,*pp;
 
        pp = ((struct lfgenf *)genptr)->si;
        p = (int *) malloc((lval-1)*sizeof(int));
        for (i=0;i<lval-1;i++) p[i] = pp[i];
        return(p);
}

#ifdef ANSI_C
int *get_next_index_rng_dbl(const int *genptr)
#else
int *get_next_index_rng_dbl(genptr)
int *genptr;
#endif
{
        int i;
        int *p,*pp;

        pp = ((struct lfgend *)genptr)->si;
        p = (int *) malloc((lval-1)*sizeof(int));
        for (i=0;i<lval-1;i++) p[i] = pp[i];
        return(p);
}
 
#ifdef ANSI_C
void si_halve_int(int *a)
#else
void si_halve_int(a)
int *a;
#endif
{
      int i;

      for (i=0;i<lval-2;i++) {
            a[i] >>= 1;
            if (a[i+1]&1) a[i] ^= 1<<MAX_BIT_INT;
            }
      a[lval-2] >>= 1;
}

#ifdef ANSI_C
void si_halve_flt(int *a)
#else
void si_halve_flt(a)
int *a;
#endif
{
        int i;
 
        for (i=0;i<lval-2;i++) {
                a[i] >>= 1;
                if (a[i+1]&1) a[i] ^= 1<<MAX_BIT_FLT;
                }
        a[lval-2] >>= 1;
}

#ifdef ANSI_C
void si_halve_dbl(int *a)
#else
void si_halve_dbl(a)
int *a;
#endif
{
/*      NOTE: generator artificially limited to BITS_IN_INT_GEN*(lval-1) */
/*           branches                                                    */
        int i;
 
        for (i=0;i<lval-2;i++) {
                a[i] = (unsigned)a[i]>>1;
                if (a[i+1]&1) a[i] ^= 1<<(BITS_IN_INT_GEN-1);
                }
        a[lval-2] = (unsigned)a[lval-2]>>1;
}


#ifdef ANSI_C
int *get_node_index_rng_int(const int *genptr)
#else
int *get_node_index_rng_int(genptr)
int *genptr;
#endif
{
      int i,*p;

      p = get_next_index_rng_int(genptr);
      if (*p == -1) return(p);
      while (!(p[0]&1)) si_halve_int(p);
      si_halve_int(p);
      return(p);
}

#ifdef ANSI_C
int *get_node_index_rng_flt(const int *genptr)
#else
int *get_node_index_rng_flt(genptr)
int *genptr;
#endif
{
        int i,*p;
 
        p = get_next_index_rng_flt(genptr);
        if (*p == -1) return(p);
        while (!(p[0]&1)) si_halve_flt(p);
        si_halve_flt(p);
        return(p);
}

#ifdef ANSI_C
int *get_node_index_rng_dbl(const int *genptr)
#else
int *get_node_index_rng_dbl(genptr)
int *genptr;
#endif
{
        int i,*p;

        p = get_next_index_rng_dbl(genptr);
        if (*p == -1) return(p);
        while (!(p[0]&1)) si_halve_dbl(p);
        si_halve_dbl(p);
        return(p);
}
 
/*************************************************************************/
/*************************************************************************/
/*                  MESSAGE PASSING ROUTINES                             */
/*************************************************************************/
/*************************************************************************/


#ifdef ANSI_C
char *pack_rng_int(int *genptr, int *size)
#else
char *pack_rng_int(genptr,size)
int *genptr,*size;
#endif
{
      int i,*p;
      struct lfgen *q;

      q = (struct lfgen *)genptr;
      *size = (2*lval+3)*sizeof(int);
      p = (int *) malloc(*size);
      if (p == NULL) {
            errprint("ERROR","pack_rng_int","insufficient memory");
            return((char *)NULL);
            }
      p[0] = lval;
      p[1] = kval;
      p[2] = gseed;
      for (i=0;i<lval;i++) p[i+3] = q->si[i];
      for (i=0;i<lval;i++) p[lval+i+3] = q->r[i];
      free(q->si);
      free(q->r);
      free(q);
      return (char *)p;
}

#ifdef ANSI_C
char *pack_rng_flt(int *genptr, int *size)
#else
char *pack_rng_flt(genptr,size)
int *genptr,*size;
#endif
{
      int i;
      char *p,*r,*s;
      struct lfgenf *q;

      q = (struct lfgenf *)genptr;
      *size = (lval+3)*sizeof(int) + lval*sizeof(float);
      p = r = (char *)malloc(*size);
      if (p == NULL) {
            errprint("ERROR","pack_rng_flt","insufficient memory");
            return(p);
            }
      s = (char *)&lval;
      for (i=0;i<sizeof(int);i++) *r++ = *s++;
      s = (char *)&kval;
      for (i=0;i<sizeof(int);i++) *r++ = *s++;
      s = (char *)&gseed;
      for (i=0;i<sizeof(int);i++) *r++ = *s++;
      s = (char *)q->si;
      for (i=0;i<(lval-1)*sizeof(int);i++) *r++ = *s++;
      s = (char *)&q->hptr;
      for (i=0;i<sizeof(int);i++) *r++ = *s++;
      s = (char *)q->r;
      for (i=0;i<lval*sizeof(float);i++) *r++ = *s++;
      free(q->si);
      free(q->r);
      free(q);
      return p;
}

#ifdef ANSI_C
char *pack_rng_dbl(int *genptr, int *size)
#else
char *pack_rng_dbl(genptr,size)
int *genptr,*size;
#endif
{
      int i;
      char *p,*r,*s;
      struct lfgend *q;
 
      q = (struct lfgend *)genptr;
      *size = (lval+3)*sizeof(int) + lval*sizeof(double);
      p = r = (char *)malloc(*size);
      if (p == NULL) {
            errprint("ERROR","pack_rng_dbl","insufficient memory");
            return(p);
            }
      s = (char *)&lval;
      for (i=0;i<sizeof(int);i++) *r++ = *s++;
      s = (char *)&kval;
      for (i=0;i<sizeof(int);i++) *r++ = *s++;
      s = (char *)&gseed;
      for (i=0;i<sizeof(int);i++) *r++ = *s++;
      s = (char *)q->si;
      for (i=0;i<(lval-1)*sizeof(int);i++) *r++ = *s++;
      s = (char *)&q->hptr;
      for (i=0;i<sizeof(int);i++) *r++ = *s++;
      s = (char *)q->r;
      for (i=0;i<lval*sizeof(double);i++) *r++ = *s++;
      free(q->si);
      free(q->r);
      free(q);
      return p;
}

#ifdef ANSI_C
int *unpack_rng_int(char *packed)
#else
int *unpack_rng_int(packed)
char *packed;
#endif
{
      int doexit=0,i,*p;
      struct lfgen *q;

/*      check values of parameters for consistency                       */
      p = (int *)packed;
      if (!lval) {
            lval = p[0];
            kval = p[1];
            gseed = p[2];
            }
      else {
            if (lval!=p[0]) doexit++;
            if (kval!=p[1]) doexit += 2;
            if (gseed!=p[2]) doexit += 4;
            }
      if (doexit) {
            if (doexit&1) errprint("ERROR","unpack_rng_int","different global L value!");
            if (doexit&2) errprint("ERROR","unpack_rng_int","different global K value!");
            if (doexit&4) errprint("ERROR","unpack_rng_int","different global seed value!");
            return((int *)NULL);
            }
      q = (struct lfgen *) malloc(sizeof(struct lfgen));
      q->si = (int *) malloc(lval*sizeof(int));
      q->r = (int *) malloc(lval*sizeof(int));
      if (q->r == NULL) {
            errprint("ERROR","unpack_rng_int","insufficient memory");
            return((int *)NULL);
            }
      for (i=0;i<lval;i++) q->si[i] = p[i+3];
      for (i=0;i<lval;i++) q->r[i] = p[i+lval+3];
      return (int *)q;
}

#ifdef ANSI_C
int *unpack_rng_flt(char *packed)
#else
int *unpack_rng_flt(packed)
char *packed;
#endif
{
      int i,temp,doexit=0;
      char *p,*s;
      struct lfgenf *q;

      p = packed;
      if (!lval) {
            s = (char *)&lval;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            s = (char *)&kval;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            s = (char *)&gseed;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            }
      else {
            s = (char *)&temp;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            if (temp!=lval) doexit++;
            s = (char *)&temp;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            if (temp!=kval) doexit += 2;
            s = (char *)&temp;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            if (temp!=gseed) doexit += 4;
            }
      if (doexit) {
            if (doexit&1) errprint("ERROR","unpack_rng_flt","different global L value!");
            if (doexit&2) errprint("ERROR","unpack_rng_flt","different global K value!");
            if (doexit&4) errprint("ERROR","unpack_rng_flt","different global seed value!");
            return((int *)NULL);
            }
      q = (struct lfgenf *) malloc(sizeof(struct lfgenf));
      q->si = (int *) malloc((lval-1)*sizeof(int));
      q->r = (float *) malloc(lval*sizeof(float));
      if (q->r == NULL) {
            errprint("ERROR","unpack_rng_flt","insufficient memory");
            return((int *)NULL);
            }
      s = (char *)q->si;
      for (i=0;i<(lval-1)*sizeof(int);i++) *s++ = *p++;
      s = (char *)&q->hptr;
      for (i=0;i<sizeof(int);i++) *s++ = *p++;
      s = (char *)q->r;
      for (i=0;i<lval*sizeof(float);i++) *s++ = *p++;
      return (int *)q;
}

#ifdef ANSI_C
int *unpack_rng_dbl(char *packed)
#else
int *unpack_rng_dbl(packed)
char *packed;
#endif
{
      int i,temp,doexit=0;
      char *p,*s;
      struct lfgend *q;

      p = packed;
      if (!lval) {
            s = (char *)&lval;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            s = (char *)&kval;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            s = (char *)&gseed;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            }
      else {
            s = (char *)&temp;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            if (temp!=lval) doexit++;
            s = (char *)&temp;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            if (temp!=kval) doexit += 2;
            s = (char *)&temp;
            for (i=0;i<sizeof(int);i++) *s++ = *p++;
            if (temp!=gseed) doexit += 4;
            }
      if (doexit) {
            if (doexit&1) errprint("ERROR","unpack_rng_dbl","different global L value!");
            if (doexit&2) errprint("ERROR","unpack_rng_dbl","different global K value!");
            if (doexit&4) errprint("ERROR","unpack_rng_dbl","different global seed value!");
            return((int *)NULL);
            }
      q = (struct lfgend *) malloc(sizeof(struct lfgend));
      q->si = (int *) malloc((lval-1)*sizeof(int));
      q->r = (double *) malloc(lval*sizeof(double));
      if (q->r == NULL) {
            errprint("ERROR","unpack_rng_dbl","insufficient memory");
            return((int *)NULL);
            }
      s = (char *)q->si;
      for (i=0;i<(lval-1)*sizeof(int);i++) *s++ = *p++;
      s = (char *)&q->hptr;
      for (i=0;i<sizeof(int);i++) *s++ = *p++;
      s = (char *)q->r;
      for (i=0;i<lval*sizeof(double);i++) *s++ = *p++;
      return (int *)q;
}


/*************************************************************************/
/*************************************************************************/
/*      FREE_RNG: remove memory for a generator                          */
/*************************************************************************/
/*************************************************************************/

#ifdef ANSI_C
void free_rng_int(int *genptr)
#else
void free_rng_int(genptr)
int *genptr;
#endif
{
      struct lfgen *q;

      q = (struct lfgen *)genptr;
      free(q->si);
      free(q->r);
      free(q);
}

#ifdef ANSI_C
void free_rng_flt(int *genptr)
#else
void free_rng_flt(genptr)
int *genptr;
#endif
{
      struct lfgenf *q;

      q = (struct lfgenf *)genptr;
      free(q->si);
      free(q->r);
      free(q);
}

#ifdef ANSI_C
void free_rng_dbl(int *genptr)
#else
void free_rng_dbl(genptr)
int *genptr;
#endif
{
      struct lfgend *q;

      q = (struct lfgend *)genptr;
      free(q->si);
      free(q->r);
      free(q);
}


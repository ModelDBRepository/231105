#ifndef SWAP
#define SWAP(a,b) {temp=(a); (a)=(b); (b)=temp;} 
#endif

static double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)

static double maxarg1,maxarg2;
#define DMAX(a,b) (maxarg1=(a),maxarg2=(b), (maxarg1) > (maxarg2) ?\
		   (maxarg1) : (maxarg2))

static double minarg1,minarg2;
#define DMIN(a,b) (minarg1=(a),minarg2=(b), (minarg1) < (minarg2) ?\
		   (minarg1) : (minarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define NR_END 1
#define FREE_ARG char*
#define EPS 1.0e-4    /* Approximate square root of the machine precision.  */ 
#define ZEPS 1.0e-10     /* A small number */
#define TINY 1.0e-20     /* A tiny number. */
#define GOLD 1.618034    /* for minimization routines */
#define GLIMIT 100.0     /* for minimization routines */
#define Rgolden 0.61803399     /* The golden ratios. */
#define Cgolden (1.0-Rgolden) 
#define SHFT(a,b,c,d)  (a)=(b);(b)=(c);(c)=(d); 
#define SHFT2(a,b,c)   (a)=(b);(b)=(c); 
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
#define ITMAX 100        /* for minimization routines */
#define CGOLD 0.3819660  /* for minimization routines */
#define MAXIT 100
#define EPSzbrent 3.0e-8
#define TINY_POWELL 1.0e-25
#define ITMAX_POWELL 200
#define TOL 2.0e-4
#define FPMIN 1.0e-30
#define EPS_GCF 3.0e-7
#define EPS_GSER 3.0e-7

#define NR_END 1
#define FREE_ARG char*

/* Function declaration */

/* nr.c */
void nrerror(char error_text[]);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
unsigned char *cvector(long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
char **cmatrix(long nrl, long nrh, long ncl, long nch);
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
     long ndl, long ndh);
int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
     long ndl, long ndh);
void correl(double data1[], double data2[], unsigned long n, double ans[]);
void convlv(double data[], unsigned long n, double respns[], unsigned long m,
     int isign, double ans[]);
void twofft(double data1[], double data2[], double fft1[], double fft2[],
     unsigned long n);
void realft(double data[], unsigned long n, int isign);
void four1(double data[], unsigned long nn, int isign);
void mrqmin(double x[], double y[], double sig[], int ndata, double a[], 
     int ia[], int ma, double **covar, double **alpha, double *chisq, 
     void (*funcs)( double, double[], double *, double[], int), 
     double *alamda);
void mrqcof(double x[], double y[], double sig[], int ndata, double a[], 
     int ia[], int ma, double **alpha, double beta[], double *chisq, 
     void (*funcs)( double, double[], double *, double[], int));
void gaussj(double **a, int n, double **b, int m);
void covsrt(double **covar, int ma, int ia[], int mfit);
void lubksb(double **a, int n, int *indx, double b[]); 
void ludcmp(double **a, int n, int *indx, double *d, int *singular); 
void mnewt(int ntrial, double x[], int n, double tolx, double tolf,
     int *singular,
     void (*usrfun)(double *x, int n, double *fvec, double **fjac, void *ptr), 
     void *ptr);
void fdjac(int n, double x[], double fvec[], double **df,
     void (*vecfunc)(int, double [], double [], void *ptr), void *ptr);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, 
     double *fc, double (*func)(double, void *), void *ptr);
double golden(double ax, double bx, double cx, double (*f)(double, void *), 
       double tol, double *xmin, void *ptr);
double brent(double ax, double bx, double cx, double (*f)(double, void *), 
       double tol, double *xmin, void *ptr);
double rtsafe(void (*funcd)(double, double *, double *, void *ptr), double x1,
        double x2, double xacc, void *ptr);
double zbrent(double (*func)(double, void *ptr), double x1, double x2, 
       double tol, void *ptr);
void balanc(double **a, int n);
void elmhes(double **a, int n);
void hqr(double **a, int n, double wr[], double wi[]);
double pythag(double a, double b);
void tred2(double **a, int n, double d[], double e[]);
void tqli(double d[], double e[], int n, double **z);
void powell(double p[], double **xi, int n, double ftol, int *iter, 
     double *fret, double (*func)(double [], void *ptr), void *ptr);
void linmin(double p[], double xi[], int n, double *fret,
     double (*func)(double [], void *ptr), void *ptr);
double f1dim(double x, void *ptr);
void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,
	double *b, double *siga, double *sigb, double *chi2, double *q);
double gammq(double a, double x);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
double gammln(double xx);

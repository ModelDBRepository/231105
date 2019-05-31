#include <stdio.h>
/* #include <stddef.h> */
/* #include <stdlib.h> */
#include <math.h>
#include "nr.h"

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate an integer matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an integer matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

char **cmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	char **m;

	/* allocate pointers to rows */
	m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
/* free a char matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate and int 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	int ***t;

	/* allocate pointers to pointers to rows */
	t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free an int i3tensor allocated by i3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

#define NRANSI

void convlv(double data[], unsigned long n, double respns[], unsigned long m,
	int isign, double ans[])
{
	unsigned long i,no2;
	double dum,mag2,*fft;

	fft=dvector(1,n<<1);
	for (i=1;i<=(m-1)/2;i++)
		respns[n+1-i]=respns[m+1-i];
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
		respns[i]=0.0;
	twofft(data,respns,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) {
		if (isign == 1) {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		} else if (isign == -1) {
			if ((mag2=SQR(ans[i-1])+SQR(ans[i])) == 0.0)
				nrerror("Deconvolving at response zero in convlv");
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		} else nrerror("No meaning for isign in convlv");
	}
	ans[2]=ans[n+1];
	realft(ans,n,-1);
	free_dvector(fft,1,n<<1);
}
#undef NRANSI

void correl(double data1[], double data2[], unsigned long n, double ans[])
{
        unsigned long no2,i;
        double dum,*fft;

        fft=dvector(1,n<<1);
        twofft(data1,data2,fft,ans,n);
        no2=n>>1;
        for (i=2;i<=n+2;i+=2) {
            ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
            ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
	}
        ans[2]=ans[n+1];
        realft(ans,n,-1);
        free_dvector(fft,1,n<<1);
}

void twofft(double data1[], double data2[], double fft1[], double fft2[],
	unsigned long n)
{
	unsigned long nn3,nn2,jj,j;
	double rep,rem,aip,aim;

	nn3=1+(nn2=2+n+n);
	for (j=1,jj=2;j<=n;j++,jj+=2) {
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}
	four1(fft1,n,1);
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

void realft(double data[], unsigned long n, int isign)
{
	void four1(double data[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}

void four1(double data[], unsigned long nn, int isign)
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

void mrqmin(double x[], double y[], double sig[], int ndata, double a[], 
int ia[], int ma, double **covar, double **alpha, double *chisq, 
void (*funcs)( double, double[], double *, double[], int), double *alamda) 
/* Levenberg-Marquardt method, attempting to reduce the value chi^2 of a fit */
/* between a set of data points x[1..ndata] , y[1..ndat a] with individual   */
/* standard deviations sig[1..ndata], and a nonlinear function dependent on  */
/* ma coeffcients a[1..ma]. The input array ia[1..ma] indicates by nonzero   */
/* entries those components of a that should be fitted for, and by zero      */
/* entries those components that should be held fixed at their input values. */
/*The program returns current best-fit values for the parameters a[1..ma] ,  */
/* and chi^2 = chisq. The arrays covar[1..ma][1. .ma], alpha[1..ma][1..ma]   */
/* are used as w rking space during most iterations. Supply a routine        */
/* funcs(x,a,yfit,dyda,ma) that evaluates the fitting function yfit, and its */
/* derivatives dyda[1..ma] with respect to the fitting parameters a at x.    */
/* On the first call p rovide an initial guess for the parameters a , and    */
/* set alamda<0 fo r initialization (which then sets alamda=.001). If a step */
/* succeeds chisq becomes smaller and alamda decreases by a facto r of 10.   */
/* If a step fails alamda grows by a factor of 10. You must call this        */
/* routine repeatedly until convergence is achieved. Then, make one final    */
/* call with alamda=0, so that covar[1..ma][1..ma] returns the covariance    */
/* matrix, and alpha the curvature matrix. (P a rameters held fixed will     */
/* return zero covariances.)                                                 */
{ 
  int j,k,l; 
  static int mfit; 
  static double ochisq,*atry ,*beta ,*da,**oneda; 

  if (*alamda < 0.0) {                                   /* Initialization.  */
    atry=dvector(1,ma); 
    beta=dvector(1,ma); 
    da=dvector(1,ma); 
      for (mfit=0,j=1;j<=ma;j++) 
        if (ia[j]) mfit++; 
      oneda=dmatrix(1,mfit,1,1); 
      *alamda=0.001; 
      mrqcof(x,y ,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs); 
      ochisq=(*chisq); 
      for (j=1;j<=ma; j++) atry[j]=a [j] ; 
    } 
    /* Alter linearized fitting matrix, by augmenting diagonal elements.    */ 
    for (j=1;j<=mfit; j++) { 
      for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k]; 
      covar[j][j]=alpha[j][j]*(1.0+(*alamda)); 
      oneda[j][1]=beta[j]; 
    } 
    gaussj(covar,mfit,oneda,1);                         /*Matrix solution.   */
    for (j=1;j<=mfit; j++) da[j]=oneda[j][1] ; 
    if (*alamda == 0.0) {    /* Once converged, evaluate covariance matrix.  */
      covsrt(covar,ma,ia,mfit); 
      free_dmatrix(oneda,1,mfit,1,1); 
      free_dvector(da,1,ma); 
      free_dvector(beta,1,ma); 
      free_dvector(atry,1,ma); 
    return; 
  } 
  for (j=0,l=1;l<=ma;l++)                         /* Did the trial succeed?  */
    if (ia[l]) atry[l]=a[l]+da[++j]; 
  mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs); 
  if (*chisq < ochisq) {               /* Success, accept the new solution.  */
    *alamda *= 0.1; 
    ochisq=(*chisq); 
    for (j=1;j<=mfit;j++) { 
      for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k] ; 
      beta[j]=da[j]; 
    } 
    for (l=1;l<=ma;l++) a[l]=atry[l]; 
  } 
  else {                            /* Failure, increase alamda and return.  */
    *alamda *= 10.0; 
    *chisq=ochisq; 
    } 
  } 

void mrqcof(double x[], double y[], double sig[], int ndata, double a[], 
int ia[], int ma, double **alpha, double beta[], double *chisq, 
void (*funcs)( double, double[], double *, double[], int)) 
/* Used by mrqmin to evaluate the linearized fitting matrix alpha, and       */
/* vector beta as in (15.5.8), and calculate chi^2.                          */
{ 
  int i,j,k,l,m,mfit=0; 
  double ymod,wt,sig2i,dy,*dyda; 

  dyda=dvector(1,ma); 
  for (j=1;j<=ma;j++) 
    if (ia[j]) mfit++; 
  for (j=1;j<=mfit; j++) {           /* Initialize (symmetric) alpha, beta.  */
    for (k=1;k<=j;k++) alpha[j][k]=0.0; 
    beta[j]=0.0; 
  } 
  *chisq=0.0; 
  for (i=1;i<=ndata ;i++) {               /* Summation loop over all data.  */
    (*funcs)(x[i],a,&ymod,dyda,ma); 
    sig2i=1.0/(sig[i]*sig[i]); 
    dy=y[i]-ymod; 
    for (j=0,l=1;l<=ma;l++) { 
      if (ia[l]) { 
        wt=dyda[l]*sig2i; 
        for (j++,k=0,m=1;m<=l;m++) 
        if (ia[m]) alpha[j][++k]+= wt*dyda[m]; 
        beta[j]+= dy*wt; 
      } 
    } 
    *chisq += dy*dy*sig2i;                               /* And find chi^2.  */
  } 
  for (j=2;j<=mfit; j++)                     /* Fill in the symmetric side.  */
  for (k=1;k<j;k++) alpha[k][j] =alpha[j][k]; 
  free_dvector(dyda,1,ma); 
}

void covsrt(double **covar, int ma, int ia[], int mfit) 
/* Expand in storage the covariance matrix covar , so as to take into    */
/* account parameters that are being held fixed. (For the latter, return */
/* zero covariances.)                                                    */
{ 
  int i,j,k; 
  double temp; 

  for (i=mfit+1; i<= ma; i++ ) 
    for (j=1;j<=i;j ++) covar[i][j]=covar[j][i]=0.0; 
  k=mfit; 
  for (j=ma;j>=1 ;j--) { 
    if (ia[j]) { 
      for (i=1;i<=ma ;i++) SWAP(covar[i][k],covar[i][j]) 
      for (i=1;i<=ma ;i++) SWAP(covar[k][i],covar[j][i]) 
      k--; 
    } 
  } 
} 

void gaussj(double **a, int n, double **b, int m) 

/* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) */
/* above. a[1..n][1..n ] is the input matrix. b[1..n][1..m ] is input     */
/* containing the m right-hand side vectors. On output, a is replaced by  */
/* its matrix inverse, and b is replaced by the corresponding set of      */
/* solution vectors.                                                      */

{ 
  int *indxc,*indxr,*ipiv; 
  int i,icol,irow,j,k,l,ll; 
  double big,dum,pivinv,temp; 
  /* The integer arrays ipiv, indxr, and indxc are used for bookkeeping   */
  /* on the pivoting.                                                     */
  indxc=ivector(1,n); 
  indxr=ivector(1,n); 
  ipiv=ivector(1,n); 
  for (j=1;j<=n; j++) ipiv[j]=0; 
  /* This is the main lo op over the columns to be reduced.               */
  for (i=1;i<=n; i++) { 
    big=0.0; 
    /* This is the outer loop of the search for a pivot element.          */
    for (j=1;j<=n;j++) 
    if (ipiv[j] != 1) 
    for (k=1;k<=n;k ++) { 
      if (ipiv[k] == 0) { 
        if (fabs(a[j][k]) >= big) { 
          big=fabs(a[j][k]); 
          irow=j; 
          icol=k; 
        } 
      } else if (ipiv[k] > 1) nrerror(" gaussj : Singular Matrix-1"); 
    }
    ++(ipiv[icol]); 
    /* We now have the pivot element, so we interchange rows, if needed,  */
    /* to put the pivot element on the diagonal. The columns are not      */
    /* physically interchanged, only relabeled: indxc[i], the column of   */
    /* the ith pivot element, is the ith column that is reduced, while    */
    /* indxr[i] is the row in which that pivot element was originally     */
    /* located. If indxr[i] != indxc[i] there is an implied column        */
    /* interchange. With this form of bookkeeping, the solution b's will  */
    /* end up in the correct order, and the inverse matrix will be        */
    /* scrambled by columns.                                              */
    if (irow != icol) { 
      for (l=1;l<=n;l++) SWAP(a[irow] [l],a [icol][l]) 
      for (l=1;l<=m;l++) SWAP(b[irow] [l],b [icol][l]) 
    } 
    /* We are now ready to divide the pivot row by the pivot element,     */
    /* located at irow and icol.                                          */
    indxr[i]=irow; 
    indxc[i]=icol; 
    if (a[icol][icol] == 0.0) nrerror(" gaussj : Singular Matrix-2"); 
    pivinv=1.0/a[icol][icol]; 
    a[icol][icol] =1.0; 
    for (l=1;l<=n;l++) a[icol][l] *= pivinv; 
    for (l=1;l<=m;l++) b[icol][l] *= pivinv; 
    /* Next, we reduce the rows...                                        */
    /*...except for the pivot one, of course.                             */
    for (ll=1;ll<=n;ll++) 
    if (ll != icol) { 
      dum=a[ll] [icol]; 
      a[ll][icol]=0.0; 
      for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum; 
      for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum; 
    } 
  } 
 /* This is the end of the main loop over columns of the reduction.       */
 /* It only remains to unscramble the solution in view of the column      */
 /* interchanges. We do this by interchanging pairs of columns in the     */
 /* reverse order that the permutation was built up.                      */
 for (l=n;l>=1; l--) { 
   if (indxr[l] != indxc[l]) 
     for (k=1;k<=n; k++) 
       SWAP(a[k] [indxr[l]],a [k][indxc[l]]); 
 } 
   /* And we are done.                                                    */
   free_ivector(ipiv ,1,n); 
   free_ivector(indxr,1,n); 
   free_ivector(indxc,1,n); 
}

void ludcmp(double **a, int n, int *indx, double *d, int *singular) 
/* Given a matrix a[1..n][ 1.. n], this routine replaces it by the LU       */
/* decomposition of a rowwise permutation of itself. a and n are input.     */
/* a is output, a rranged as in equation (2.3.14) above; indx[1..n] is an   */
/* output vector that records the row permutation effected by the partial   */
/* pivoting; d is output as +-1 depending on whether the number of row      */
/* interchanges was even or odd, respectively.                              */
/* This routine is used in combination with lubksb to solve linear          */
/* equations or invert a matrix.                                            */
{ 
  int i,imax,j,k ; 
  double big,dum,sum ,temp; 
  double *vv;               /* vv stores the implicit scaling of each row.   */

  vv=dvector(1 ,n) ; 
  *singular = 1;

  *d=1.0;                  /* No row interchanges yet.                      */
  for (i=1;i<=n; i++) {  
    /* Loop over ro ws to get the implicit scaling information.             */
    big=0.0; 
    for (j=1;j<=n;j ++) 
      if ((temp=fabs (a[ i] [j] )) > big) big=temp;
    if (big == 0.0)  /* Modification DG 27.5.2002 */
    {
      printf("Singular matrix in routine ludcmp\n");
      *singular = 0;
      return;
    }
    /* No nonzero largest element.                                          */
    vv[i]=1.0/ big ;   /* Save the scaling.                                 */
  } 
  for (j=1;j<=n; j++ ) {  
  /* This is the loop over columns of Crout's method.                       */
    for (i=1;i<j;i++) { 
    /* This is equation (2.3.12) except for i = j.                          */
      sum=a[i][j] ; 
      for (k=1;k<i;k ++) sum -= a[i][k]*a[k][j]; 
      a[i][j]=sum ; 
    } 
    big=0.0;   /* Initialize for the search for largest pivot element.      */
    for (i=j;i<=n;i ++) {   
    /* This is i = j of equation (2.3.12) and i = j + 1 ... N  of equation  */
    /* (2.3.13).                                                            */
      sum=a[i][j] ; 
      for (k=1;k<j;k ++) 
        sum -= a[i][k]*a[k][j]; 
        a[i][j]=sum ; 
        if ( (dum=vv[i]*fabs(sum) ) >= big) { 
        /* Is the figure of merit for the pivot better than the best so     */
        /* far?                                                             */ 
          big=dum; 
          imax=i; 
        } 
      } 
      if (j != imax) {   /* Do we need to interchange rows?                 */
      for (k=1;k<=n; k++) {   /* Yes, do so...                              */
        dum=a[imax][k]; 
        a[imax][k]=a [j][k] ; 
        a[j][k]=dum; 
      } 
      *d = -(*d);        /* ...and change the parity of d .                 */
      vv[imax]=vv [j] ;  /* Also interchange the scale facto r.             */
    } 
    indx[j]=imax; 
    if (a[j][j] == 0.0) a[j][j]=TINY ; 
    /* If the pivot element is zero the matrix is singula r (at least to    */
    /* the precision of the algorithm). For some applications on singular   */
    /* matrices, it is desirable to substitute TINY for zero.               */
    if (j != n) {            /* Now, finally , divide by the pivot element. */
      dum=1.0/(a[j][j]) ; 
      for (i=j+1;i<= n;i++) a[i][j] *= dum; 
    } 
  } 
  /* Go back for the next column in the reduction.                          */
  free_dvector (vv ,1, n); 
} 

void lubksb(double **a, int n, int *indx, double b[]) 
/* Solves the set of n linea r equations A * X = B . Here a[1..n][1 ..n ]   */
/* is input, not as the matrix A but rather as its LU decomposition,        */
/* determined by the routine ludcmp. indx[1..n] is input as the permutation */
/* vector returned by ludcmp. b[1..n] is input as the right-hand side       */
/* vector B, and returns with the solution vector X . a , n , and indx are  */
/* not modified by this routine and can be left in place for successive     */
/* calls with different right-hand sides b. This routine takes into account */
/* the possibility that b will begin with many zero elements, so it is      */
/* efficient for use in matrix inversion.                                   */
{ 
  int i,ii=0,ip,j; 
  double sum; 
  for (i=1;i<=n; i++ ) { 
  /* When ii is set to a positive value, it will become the index of the    */
  /* first nonvanishing element of b. We now do the forward substitution,   */
  /* equation (2.3.6). The only new wrinkle is to unscramble the            */
  /* permutation as we go.                                                  */
    ip=indx[i]; 
    sum=b[ip]; 
    b[ip]=b[i]; 
    if (ii) 
      for (j=ii;j<=i-1; j++) sum -= a[i][j]*b[j]; 
    else if (sum) ii=i; 
    /* A nonzero element w as encountered, so from now on we will have to   */
    /* do the sums in the loop above.                                       */ 
    b[i]=sum; 
  } 
  for (i=n;i>=1; i--) { 
  /* Now we do the backsubstitution, equation (2.3.7).                      */
    sum=b[i]; 
    for (j=i+1;j<=n ;j++) sum -= a[i][j]*b[j]; 
    b[i]=sum/a[i][i]; 
    /* Store a component of the solution vector X .                         */
  }   /* All done!                                                          */
}

#define FREERETURN {free_dmatrix(fjac,1,n,1,n);free_dvector(fvec,1,n);\
    free_dvector(p,1,n);free_ivector(indx,1,n);return;} 

void mnewt(int ntrial, double x[], int n, double tolx, double tolf,
     int *singular,
     void (*usrfun)(double *x, int n, double *fvec, double **fjac, void *ptr), 
     void *ptr)
/* Given an initial guess x[1..n] for a root in n dimensions, take ntrial  */
/* Newton-Raphson steps to improve the root. Stop if the root converges in */
/* either summed absolute variable increments tolx or summed absolute      */
/* function values tolf.                                                   */
{ 
  int k,i,*indx; 
  double errx,errf,d,*fvec,**fjac,*p; 
  indx=ivector(1,n) ; 
  p=dvector(1, n);
  fvec=dvector(1, n); 
  fjac=dmatrix(1, n, 1 ,n) ; 
  for (k=1;k<=ntrial;k++) { 
    usrfun(x, n ,fvec ,fjac, ptr);
/* User function supplies function values at x in fvec and Jacobian       */
/* matrix in fjac.                                                        */
    errf=0.0;
    for (i=1;i<=n;i++) errf += fabs(fvec[i]);  
/*  Check function convergence.                                           */
    if (errf <= tolf) FREERETURN 
    for (i=1;i<=n;i ++) p[i] = -fvec[i]; 
/*  Right-hand side of linear equations.                                  */
    ludcmp(fjac,n ,indx ,&d, singular); 
/*  Solve linea r equations using LU decomp osition.                      */
    lubksb(fjac,n,indx,p); 
    errx=0.0; 
/*  Check root convergence.                                               */
    for (i=1;i<=n;i ++) {           /* Update solution.                   */
      errx += fabs(p[i]); 
      x[i] += p[i]; 
    } 
    if (errx <= tolx) FREERETURN 
  } 
  FREERETURN 
}

void fdjac(int n, double x[], double fvec[], double **df,
     void (*vecfunc)(int, double [], double [], void *ptr), void *ptr)
/* Computes forward-difference approximation to Jacobian. On input, x[1..n] */ 
/* is the point at which the Jacobian is to be evaluated, fvec[1..n] is the */
/* vector of function values at the point, and vecfunc(n,x,f) is a          */
/* user-supplied routine that returns the vector of functions at x.         */
/* On output, df[1..n][ 1.. n] is the Jacobian array.                       */
{ 
  int i,j;
  double h,temp,*f; 

  f=dvector(1, n); 
  for (j=1;j<=n; j++ ) { 
    temp=x[j]; 
    h=EPS*fabs(temp); 
    if (h==0.0) h=EPS; 
    x[j]=temp+h;         /* Trick to reduce finite precision error.       */ 
    h=x[j]-temp; 
    (*vecfunc) (n,x,f,ptr); 
    x[j]=temp; 
    for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;  
/*  Forward difference formula.                                           */
  }
  free_dvector(f,1,n);
}

/* Here GOLD is the default ratio by which successive intervals are         */
/* magnified; GLIMIT is the maximum magnification allowed for a             */
/* parabolic-fit step.                                                      */

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, 
     double *fc, double (*func)(double, void *), void *ptr) 
/* Given a function func, and given distinct initial points ax and bx, this */
/* routine searches in the downhill direction (defined by the function as   */
/* evaluated at the initial points) and returns new points ax, bx, cx that  */
/* bracket a minimum of the function. Also returned are the function values */
/* at the three points, fa, fb, and fc.                                     */
{ 
  double ulim,u,r,q,fu,dum;

  *fa=(*func)(*ax, ptr); 
  *fb=(*func)(*bx, ptr); 
  /* printf("ax=%lf fa=%lf bx=%lf fb=%lf\n", *ax, *fa, *bx, *fb); */
  if (*fb > *fa) { 
  /* Switch roles of a and b so that we can go  */
  /* downhill in the direction from a to b .    */
  SHFT(dum, *ax, *bx, dum) 
  SHFT(dum, *fb, *fa, dum) 
  } 
  /* printf("ax=%lf fa=%lf bx=%lf fb=%lf\n", *ax, *fa, *bx, *fb); */
  *cx=(*bx)+GOLD*(*bx - *ax);   /* First guess for c. */
  /* printf("cx=%lf\n", *cx); */
  *fc=(*func)(*cx, ptr); 
  /* printf("ax=%lf fa=%lf bx=%lf fb=%lf cx=%lf fc=%lf\n", *ax, *fa, *bx, *fb,
   *cx, *fc);*/
  while (*fb > *fc) {  /* Keep returning here until we bracket. */
  r=(*bx - *ax)*( *fb - *fc) ;  /* Compute u by parabolic extrapolation from */
  q=(*bx - *cx)*( *fb - *fa) ;  /* a, b, c . TINY is used to prevent any     */
  u=(*bx)-((*bx - *cx)*q-(*bx - *ax)*r)/    /* possible division by zero.    */
  (2.0*SIGN(DMAX(fabs(q-r),TINY),q-r));
  ulim=(*bx)+GLIMIT*(*cx - *bx); 
  /* We won't go farther than this. Test various possibilities:  */
  if ((*bx-u)*(u - *cx) > 0.0) { /* Parabolic u is between b and c : try it. */
    fu=(*func)(u, ptr); 
    /* printf("u=%lf fu=%lf\n", u, fu); */
    if (fu < *fc) {   /* Got a minimum between b and c . */ 
      *ax=(*bx); 
      *bx=u; 
      *fa=(*fb); 
      *fb=fu; 
      return; 
    } 
    else if (fu > *fb) { /* Got a minimum between a and u . */
      *cx=u; 
      *fc=fu; 
      return; 
    } 
    u=(*cx)+GOLD*(*cx - *bx);       /* Parabolic fit was no use. Use default */
    fu=(*func)(u, ptr);             /* magnification.                        */
    /* printf("u=%lf fu=%lf\n", u, fu); */
  } 
  else if ((*cx-u)*(u-ulim) > 0.0) {  /* Parabolic fit is between c and its  */
    fu=(*func)(u, ptr);               /* allowed limit.                      */
    /* printf("u=%lf fu=%lf\n", u, fu); */
    if (fu < *fc) { 
      SHFT(*bx, *cx ,u,  *cx+GOLD*(*cx - *bx)) 
      SHFT(*fb, *fc ,fu ,(*func)(u, ptr)) 
	/* printf("u=%lf fu=%lf\n", u, fu); */
    } 
  } else if ((u-ulim)*( ulim- *cx) >= 0.0) {  /* Limit parabolic u to   */
      u=ulim;                                 /* maximum allowed value. */
      fu=(*func)(u, ptr); 
      /* printf("u=%lf fu=%lf\n", u, fu); */
    } 
    else {   /* Reject parabolic u , use default magnification. */
      u=(*cx)+GOLD*( *cx - *bx); 
      fu=(*func)(u, ptr); 
      /* printf("u=%lf fu=%lf\n", u, fu); */
    } 
    SHFT(*ax, *bx, *cx, u)  /* Eliminate oldest point and continue. */
    SHFT(*fa,* fb, *fc, fu) 
  } 
} 

double golden(double ax, double bx, double cx, double (*f)(double, void *), 
       double tol, double *xmin, void *ptr) 
/* Given a function f, and given a bracketing triplet of abscissas ax, bx,  */
/* cx (such that bx is between ax and cx, and f(bx) is less than both f(ax) */
/* and f(cx)), this routine performs a  golden section search for the       */
/* minimum, isolating it to a fractional precision of about tol. The        */
/* abscissa of the minimum is returned as xmin , and the minimum function   */
/* value is returned as golden , the returned function value.               */
{ 
  double f1,f2,x0,x1,x2,x3;  

  /* At any given time we will keep track of four points, x0,x1,x2,x3.      */ 
  x3=cx; 
  if (fabs(cx-bx) > fabs(bx-ax)) {  /* Make x0 to x1 the smaller segment,   */
    x1=bx;                          /* and fill in the new point to be      */
    x2=bx+Cgolden*(cx-bx);                /* tried.                               */
  } 
  else { 
    x2=bx; 
    x1=bx-Cgolden*(bx-ax); 
  } 
  /* The initial function evaluations. Note that we never     */
  /* need to evaluate the function at the original endpoints. */
  f1=(*f)(x1, ptr);  
  f2=(*f)(x2, ptr);  
  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) { 
    if (f2 < f1) {                      /* One possible outcome,            */ 
      SHFT3(x0,x1,x2,Rgolden*x1+Cgolden*x3)         /* its housekeeping,                */
      SHFT2(f1,f2,(*f)(x2, ptr))        /* and a new function evaluation.   */
    } else {                            /* The other outcome,               */
      SHFT3(x3,x2,x1,Rgolden*x2+Cgolden*x0) 
      SHFT2(f2,f1,(*f)(x1, ptr))        /* and its new function evaluation. */
    } 
  }                                     /* Back to see if we are done.      */
  if (f1 < f2) { /* We are done. Output the best of the two current values. */
    *xmin=x1; 
    return f1; 
  } else { 
    *xmin=x2; 
    return f2; 
  } 
} 


/* Here ITMAX is the maximum allowed number of iterations; CGOLD is the      */
/* golden ratio; ZEPS is a small number that protects against trying to      */
/* achieve fractional accuracy for a minimum that happens to be exactly zero.*/
double brent(double ax, double bx, double cx, double (*f)(double, void *), 
       double tol, double *xmin, void *ptr) 
/* Given a function f, and given a bracketing triplet of abscissas ax, bx,   */
/* cx (such that bx is between ax and cx, and f(bx) is less than both f(ax)  */
/* and f(cx)), this routine isolates the minimum to a fractional precision   */
/* of about tol using Brent's method. The abscissa of the minimum is         */
/* returned as xmin , and the minimum function value is returned as brent ,  */
/* the returned function value.                                              */
{ 
  int iter; 
  double a,b,d,etemp=0,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm; 
  double e=0.0;  /* This will be the distance moved on the step before last. */
  a=(ax < cx ? ax : cx);  /* a and b must be in ascending order,             */
  b=(ax > cx ? ax : cx);  /* but input abscissas need not be.                */
  x=w=v=bx;               /* Initializations...                              */
  fw=fv=fx=(*f)(x, ptr); 
  for (iter=1; iter<=ITMAX; iter++) {               /* Main program loop.    */
    xm=0.5*(a+b); 
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS); 
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {          /* Test for done here.    */
      *xmin=x; 
      return fx; 
    } 
    if (fabs(e) > tol1) {              /* Construct a trial parabolic fit.   */
      r=(x-w)*(fx-fv); 
      q=(x-v)*(fx-fw); 
      p=(x-v)*q-(x-w)*r; 
      q=2.0*(q-r) ; 
      if (q > 0.0) p = -p; 
      q=fabs(q); 
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) 
        d=CGOLD*(e=( x >= xm ? a-x : b-x)); 
      /* The above conditions determine the acceptability of the parabolic  */
      /* fit. Here we take the golden section step into the larger of the   */
      /* two segments.                                                      */
      else { 
        d=p/q;                                /* Take the parabolic step.   */
        u=x+d; 
        if (u-a < tol2 || b-u < tol2) 
          d=SIGN(tol1 ,xm-x ); 
        } 
    } 
    else { 
      d=CGOLD*(e= (x >= xm ? a-x : b-x)); 
    } 
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d)); 
    fu=(*f)(u, ptr); /* This is the one function evaluation per iteration.  */
    if (fu <= fx) {  /* Now decide what to do with our function evaluation. */
      if (u >= x) a=x; else b=x; 
      SHFT(v,w,x,u)                      /*Housekeepingfollows:             */ 
      SHFT(fv,fw,fx,fu) 
    } 
    else { 
      if (u < x) a=u; else b=u; 
      if (fu <= fw || w == x) { 
        v=w; 
        w=u; 
        fv=fw; 
        fw=fu; 
      } else if (fu <= fv || v == x || v == w) { 
        v=u; 
        fv=fu; 
      } 
    }                /* Done with housekeeping. Back for another iteration. */
  } 
  nrerror("Too many iterations in brent"); 

  return 0.0;
} 

double rtsafe(void (*funcd)(double, double *, double *, void *ptr), double x1,
        double x2, double xacc, void *ptr)
{
	void nrerror(char error_text[]);
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;

	(*funcd)(x1,&fl,&df, ptr);
	(*funcd)(x2,&fh,&df, ptr);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		nrerror("Root must be bracketed in rtsafe");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	(*funcd)(rts,&f,&df, ptr);
	for (j=1;j<=MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		(*funcd)(rts,&f,&df, ptr);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}

double zbrent(double (*func)(double, void *ptr), double x1, double x2, 
       double tol, void *ptr)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, ptr),fb=(*func)(b, ptr),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		nrerror("Root must be bracketed in zbrent");
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPSzbrent*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(b, ptr);
	}
	nrerror("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}

#define RADIX 2.0

void balanc(double **a, int n)
{
	int last,j,i;
	double s,r,g,f,c,sqrdx;

	sqrdx=RADIX*RADIX;
	last=0;
	while (last == 0) {
		last=1;
		for (i=1;i<=n;i++) {
			r=c=0.0;
			for (j=1;j<=n;j++)
				if (j != i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			if (c && r) {
				g=r/RADIX;
				f=1.0;
				s=c+r;
				while (c<g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g=r*RADIX;
				while (c>g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s) {
					last=0;
					g=1.0/f;
					for (j=1;j<=n;j++) a[i][j] *= g;
					for (j=1;j<=n;j++) a[j][i] *= f;
				}
			}
		}
	}
}
#undef RADIX

void elmhes(double **a, int n)
{
	int m,j,i;
	double y,x;
        double temp;

	for (m=2;m<n;m++) {
		x=0.0;
		i=m;
		for (j=m;j<=n;j++) {
			if (fabs(a[j][m-1]) > fabs(x)) {
				x=a[j][m-1];
				i=j;
			}
		}
		if (i != m) {
			for (j=m-1;j<=n;j++) SWAP(a[i][j],a[m][j])
			for (j=1;j<=n;j++) SWAP(a[j][i],a[j][m])
		}
		if (x) {
			for (i=m+1;i<=n;i++) {
				if ((y=a[i][m-1]) != 0.0) {
					y /= x;
					a[i][m-1]=y;
					for (j=m;j<=n;j++)
						a[i][j] -= y*a[m][j];
					for (j=1;j<=n;j++)
						a[j][m] += y*a[j][i];
				}
			}
		}
	}
}

#define NRANSI

void hqr(double **a, int n, double wr[], double wi[])
{
	int nn,m,l,k,j,its,i,mmin;
	double z,y,x,w,v,u,t,s,r,q,p,anorm;

	anorm=0.0;
	for (i=1;i<=n;i++)
		for (j=IMAX(i-1,1);j<=n;j++)
			anorm += fabs(a[i][j]);
	nn=n;
	t=0.0;
	while (nn >= 1) {
		its=0;
		do {
			for (l=nn;l>=2;l--) {
				s=fabs(a[l-1][l-1])+fabs(a[l][l]);
				if (s == 0.0) s=anorm;
				if ((double)(fabs(a[l][l-1]) + s) == s) break;
			}
			x=a[nn][nn];
			if (l == nn) {
				wr[nn]=x+t;
				wi[nn--]=0.0;
			} else {
				y=a[nn-1][nn-1];
				w=a[nn][nn-1]*a[nn-1][nn];
				if (l == (nn-1)) {
					p=0.5*(y-x);
					q=p*p+w;
					z=sqrt(fabs(q));
					x += t;
					if (q >= 0.0) {
						z=p+SIGN(z,p);
						wr[nn-1]=wr[nn]=x+z;
						if (z) wr[nn]=x-w/z;
						wi[nn-1]=wi[nn]=0.0;
					} else {
						wr[nn-1]=wr[nn]=x+p;
						wi[nn-1]= -(wi[nn]=z);
					}
					nn -= 2;
				} else {
					if (its == 30) nrerror("Too many iterations in hqr");
					if (its == 10 || its == 20) {
						t += x;
						for (i=1;i<=nn;i++) a[i][i] -= x;
						s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
						y=x=0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=(nn-2);m>=l;m--) {
						z=a[m][m];
						r=x-z;
						s=y-z;
						p=(r*s-w)/a[m+1][m]+a[m][m+1];
						q=a[m+1][m+1]-z-r-s;
						r=a[m+2][m+1];
						s=fabs(p)+fabs(q)+fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
						v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
						if ((double)(u+v) == v) break;
					}
					for (i=m+2;i<=nn;i++) {
						a[i][i-2]=0.0;
						if (i != (m+2)) a[i][i-3]=0.0;
					}
					for (k=m;k<=nn-1;k++) {
						if (k != m) {
							p=a[k][k-1];
							q=a[k+1][k-1];
							r=0.0;
							if (k != (nn-1)) r=a[k+2][k-1];
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							if (k == m) {
								if (l != m)
								a[k][k-1] = -a[k][k-1];
							} else
								a[k][k-1] = -s*x;
							p += s;
							x=p/s;
							y=q/s;
							z=r/s;
							q /= p;
							r /= p;
							for (j=k;j<=nn;j++) {
								p=a[k][j]+q*a[k+1][j];
								if (k != (nn-1)) {
									p += r*a[k+2][j];
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j] -= p*x;
							}
							mmin = nn<k+3 ? nn : k+3;
							for (i=l;i<=mmin;i++) {
								p=x*a[i][k]+y*a[i][k+1];
								if (k != (nn-1)) {
									p += z*a[i][k+2];
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l < nn-1);
	}
}

#undef NRANSI

void tred2(double **a, int n, double d[], double e[])
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[1]=0.0;
	e[1]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
		l=i-1;
		if (d[i]) {
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

#define NRANSI

void tqli(double d[], double e[], int n, double **z)
{
	double pythag(double a, double b);
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) nrerror("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=1;k<=n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}
#undef NRANSI


#define NRANSI

void powell(double p[], double **xi, int n, double ftol, int *iter, 
     double *fret, double (*func)(double [], void *ptr), void *ptr)
{
	int i,ibig,j;
	double del,fp,fptt,t,*pt,*ptt,*xit;

	pt=dvector(1,n);
	ptt=dvector(1,n);
	xit=dvector(1,n);
	*fret=(*func)(p,ptr);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin(p,xit,n,fret,func, ptr);
			if (fptt-(*fret) > del) {
				del=fptt-(*fret);
				ibig=i;
			}
		}
		if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY_POWELL) {
			free_dvector(xit,1,n);
			free_dvector(ptt,1,n);
			free_dvector(pt,1,n);
			return;
		}
		if (*iter == ITMAX_POWELL) nrerror("powell exceeding maximum iterations.");
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt,ptr);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func, ptr);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
#undef NRANSI

/* note #undef's at end of file */
#define NRANSI

int ncom;
double *pcom,*xicom,(*nrfunc)(double [], void *ptr);

void linmin(double p[], double xi[], int n, double *fret,
     double (*func)(double [], void *ptr), void *ptr)
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=dvector(1,n);
	xicom=dvector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim, ptr);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin, ptr);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_dvector(xicom,1,n);
	free_dvector(pcom,1,n);
}
#undef NRANSI

#define NRANSI

extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double [], void *ptr);

double f1dim(double x, void *ptr)
{
	int j;
	double f,*xt;

	xt=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt, ptr);
	free_dvector(xt,1,ncom);
	return f;
}
#undef NRANSI

#define NRANSI

void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,
	double *b, double *siga, double *sigb, double *chi2, double *q)
{
  /* double gammq(double a, double x); */
	int i;
	double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

	*b=0.0;
	if (mwt) {
		ss=0.0;
		for (i=1;i<=ndata;i++) {
			wt=1.0/SQR(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	} else {
		for (i=1;i<=ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt) {
		for (i=1;i<=ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	} else {
		for (i=1;i<=ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0;
	*q=1.0;
	if (mwt == 0) {
		for (i=1;i<=ndata;i++)
			*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	} else {
		for (i=1;i<=ndata;i++)
			*chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
		if (ndata>2) *q=gammq(0.5*(ndata-2),0.5*(*chi2));
	}
}
#undef NRANSI

double gammq(double a, double x)
{
  /* void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]); */
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

void gcf(double *gammcf, double a, double x, double *gln)
{
  /* double gammln(double xx);
     void nrerror(char error_text[]); */
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS_GCF) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

void gser(double *gamser, double a, double x, double *gln)
{
  /* double gammln(double xx);
     void nrerror(char error_text[]); */
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS_GSER) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

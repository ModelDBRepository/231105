#ifndef min
#define min(x,y) ( ((x) < (y)) ? (x) : (y) )
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

#ifndef NR_END
#define NR_END 1
#endif

#ifndef FREE_ARG
#define FREE_ARG char*
#endif

#define Mline 10000

/* Structure Declaration */

typedef struct fl_st{
  FILE *in, *dat, *cor, *out, *res;
} fl_st;

typedef struct crc_par{
  double deltat;
  int nicol2ar, *icol2ar;
  int ncol, n_line_skip, n_line_read, icol[3];
  double epsilon;
  double base_intersect;
  double Tmax, T_between_peaks;
  int kpoints, smap;
  int ntc, inter_points;
  char detcol;
} crc_par;

typedef struct peak_str{
  double amplitude, **t_isect, t_period, t_shift, phase;
  int *iarmx, *iarmn;
  int mp, np;
  int zp;
} peak_str;


/* Function Declaration */

void read_input(crc_par *crcpar, fl_st fl);
void write_input(crc_par *crcpar, fl_st fl);
void one_cor_cal_quantify(crc_par *crcpar, double **datar, double *corar,
     fl_st fl);
void read_data(crc_par *crcpar, double **datar, fl_st fl);
void calculate_correlation(crc_par *crcpar, double **datar, int idat1,
     int idat2, double *corar, fl_st fl);
void amp_ph_cal (crc_par *crcpar, int idat1, int idat2, double *corar,
     fl_st fl);
void find_local_max(crc_par *crcpar, peak_str *peakstr, double *corar,
     fl_st fl);
int find_very_local_max(int ii, crc_par *crcpar, double *corar, fl_st fl);
void find_peak_near_zero(crc_par *crcpar, peak_str *peakstr, double *corar,
     fl_st fl);
void find_local_min(crc_par *crcpar, peak_str *peakstr, double *corar,
     fl_st fl);
void find_intersections_with_zero(crc_par *crcpar, peak_str *peakstr, 
     double *corar, fl_st fl);
void find_one_intersection(int pleft, int pright, crc_par *crcpar, 
     double *corar, double *isect, fl_st fl);
int cond_intersect(int ip, crc_par *crcpar, double *corar, fl_st fl);
void linear_regression(int ipl, int ipr, crc_par *crcpar, double *corar,
     double *t_isect, fl_st fl);

void xcorrel(double data1[], double data2[],  unsigned long n, double ans[]);
void xtwofft(double data1[], double data2[], double fft1[], double fft2[],
     unsigned long n);
void xfour1(double data[], unsigned long nn, int isign);

#define Mpar 500
#define Mline 2000
#define Mword 100
#define Mdatraw 200
#define Mdatcol 1000

/* Structure Declaration */

typedef struct scan_val{
  double parmin, parmax, par_ar[Mpar];
  int npar, ipar, npt, nrepeat, irepeat;
  int seed;
  char scan_type;
  char par1[Mword], par2[Mword];
} scan_val;


/* Function Declaration */
void read_first_input_line(scan_val *sval, fl_st fl);
void update_file_old(scan_val *sval, int skip_lines, fl_st fl);
int process_line_old(scan_val *sval, char line[], FILE *ftmp);
void write_avr(par_all *parall, avr_val *av, fl_st fl);
int process_line_no_colon_old(scan_val *sval, char line[], FILE *ftmp);
int process_line_yes_colon_old(scan_val *sval, char line[], FILE *ftmp);



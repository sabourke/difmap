#ifndef vlbutil_h
#define vlbutil_h

void   radhms(double rad, int *hour, int *mins, double *secs);
char *sradhms(double rad, int precision, int colon, char *string);
void   raddms(double rad, int *sgn, int *deg,  int *mins, double *secs);
char *sraddms(double rad, int precision, int colon, char *string);
void   daydate(int year, int dayno, int *day, int *month);
char *sdaydate(int year, int dayno, char *string);
char *sutdate(int year, double vlbut, char *string);
void   dayut(double vlbut, int *dayno, int *hour, int *mins, double *secs);
int write_ut(double vlbut, int nc, char *string);
double read_ut(char *s, char **endp);
void julday(double vlbut, int year, long *jd, double *jdfrc, double *je);
char *date_str(void);
char *polname(int polcode);
void subamphs(float *amp, float *phs, float subamp, float subphs);
void addamphs(float *amp, float *phs, float addamp, float addphs);
void add_cart_to_polar(float *amp, float *phs, float re, float im);
int termstr(char *instr, int slen);
char *termcpy(char *ostr, const char *istr, int ncmax);
void fillstr(char *instr, int slen);
char *stripcpy(char *ostr, int nco, const char *istr, int nci);
char *stripstr(char *istr, int nci);
void imran(float *map, int xdim, int ydim, int xa, int xb, int ya, int yb,
	    float *mapmin, float *mapmax);
int plbeam(float bmin, float bmaj, float bpa, float xpos, float ypos,
	   float xmin, float xmax, float ymin, float ymax);

int parse_sexagesimal_string(char *string, double *value, char **endp);

typedef struct {
  enum {
    NUM_DOUBLE,   /* Number was written as an int with no exponent */
    NUM_INT       /* Number was written as floating point and/or with exponent*/
  } type;
  union {
    double dval;  /* Use if type==NUM_DOUBLE */
    int ival;     /* Use if type==NUM_INT */
  } value;
} Number;

int parse_numeric_string(char *string, char **endp, Number *number);

int write_string_arg(FILE *fp, char *fname, char *string);

int parse_date_and_time(const char *string, int tell, int nospace,
			const char **endp, int *year, int *month, int *day,
			int *hour, int *min, double *sec);
int parse_time(const char *string, int tell, const char **endp,
	       int *hour, int *min, double *sec);
int parse_date(const char *string, int tell, const char **endp,
	       int *year, int *month, int *day);
int parse_ulong(const char *string, int tell, const char **endp,
		unsigned long *ulval);
int parse_double(const char *string, int tell, const char **endp,
		 double *dval);
int parse_mjd(const char *string, int tell, const char **endp, double *mjd);

#endif

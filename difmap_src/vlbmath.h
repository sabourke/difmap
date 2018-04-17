/* Simple numeric functions */

int fnint(float fval);  /* Nearest integer for float arguments */
int dnint(double dval); /* Nearest integer for double arguments */
void frange(float *vec, int vecdim, float *vecmin, float *vecmax);
void imran(float *map, int xdim, int ydim, int xa, int xb, int ya, int yb,
	    float *mapmin, float *mapmax);
float floatmax(float a, float b);   /* Return the max float argument */
int imax(int a, int b);         /* Return the max int argument */
double dmax(double a, double b);/* Return the max double argument */
float floatmin(float a, float b);   /* Return the min float argument */
int imin(int a, int b);         /* Return the min int argument */
double dmin(double a, double b);/* Return the min double argument */
void costran(float *inparr, int ninp, float inwid, float *outarr, int nout);

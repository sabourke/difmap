#ifndef phdu_h
#define phdu_h

#ifndef fits_h
#include "fits.h"
#endif

/* Define a primary-header axis descriptor */

typedef struct {
  char *ctype;
  double crpix;
  double crval;
  double cdelt;
  double crota;
} Imaxis;

/* Define a primary-header group parameter descriptor */

typedef struct {
  char *ptype;
  double pscal;
  double pzero;
} Gpar;

/* Primary and IUE-IMAGE header descriptor. */

typedef struct PHdu {
  HDUBase
  char *origin;
  char *date_obs;
  char *telescop;
  char *instrume;
  char *observer;
  char *object;
  char *author;
  char *referenc;
  double equinox;
  double bscale;
  double bzero;
  char *bunit;
  long blank;
  Imaxis *axes;    /* Array of axis descriptors (CTYPE,CRPIX,CRVAL...) */
  Gpar *pars;      /* Array of group parameter descriptors (PTYPE,PSCAL...) */
  double datamax;
  double datamin;
  long imsize;     /* Number of elements per group image-array */
} Phdu;

Phdu *find_image(Fits *fits, char *extname, int extver, Hdu *prev);

int find_axis(Phdu *phdu, const char *ctype, int fixlen, int start);
char *axis_name(Phdu *phdu, int axis);
Imaxis *get_axis(Phdu *phdu, int axis);

int find_gpar(Phdu *phdu, const char *ptype, int fixlen, int start);
char *gpar_name(Phdu *phdu, int ipar);
Gpar *get_gpar(Phdu *phdu, int ipar);

long rgroup(Fits *fits, Phdu *phdu, long igroup, long start, long nobj,
	       Fittype type, int doscale, Fitsflag *flags, void *data);
long wgroup(Fits *fits, Phdu *phdu, long igroup, long start, long nobj,
	    Fittype type, int doscale, Fitsflag *flags, void *data);
long rimage(Fits *fits, Phdu *phdu, long igroup, long start, long nobj,
	       Fittype type, int doscale, Fitsflag *flags, void *data);
long wimage(Fits *fits, Phdu *phdu, long igroup, long start, long nobj,
	       Fittype type, int doscale, Fitsflag *flags, void *data);

#endif

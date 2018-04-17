#ifndef utils_h
#define utils_h

/* Private generic HDU implementation functions */

Hdu *new_Hdu(Hdutype type);
Hdu *ini_Hdu(Hdu *hdu, Bitpix bitpix, int *dims, int naxis,
	     int groups, int pcount, int gcount, char *extname, int extver,
	     int extlevel, int headrec, int endline);
Hdu *get_Hdu(struct Fits *fits, int headrec);

/* Declare separate functions to handle each of the derived HDU types */

#define GETFN(n) int (n)(struct Fits *fits, Hdu *hdu)
#define NEWFN(n) Hdu *(n)(Hdu *hdu)
#define DELFN(n) void (n)(Hdu *hdu)
#define SIZEFN(n) size_t (n)(void)
#define ADDFN(n) int (n)(struct Fits *fits, Hdu *hdu)
#define COPYFN(n) Hdu *(n)(Hdu *hdu)
#define ENDFN(n) int (n)(struct Fits *fits, Hdu *hdu)

typedef struct {   /* Virtual function table for derived HDU types */
  char *name;      /* Standard name to which refer to HDU type by */
  char *aips;      /* AIPS name to refer to HDU by */
  GETFN(*getfn);   /* Function for reading derived HDU parts */
  NEWFN(*newfn);   /* Function for zeroing derived HDU parts */
  DELFN(*delfn);   /* Function for deleting derived HDU parts */
  SIZEFN(*sizefn); /* Function to return size of descriptor */
  ADDFN(*addfn);   /* Function to write derived HDU header lines */
  COPYFN(*copyfn); /* Function to copy an HDU descriptor */
  ENDFN(*endfn);   /* Function to finish the data section of an HDU */
} Hdutab;

NEWFN(new_table);
DELFN(del_table);

/* Private FITS implementation functions */

size_t typesize(Fittype type);
size_t machsize(Fittype type);

/*
 * In the type conversion functions the value being converted is
 * first scaled by a given factor and then offset. Define a type to be used
 * to specify offset and scale factors.
 */
typedef struct {
  double off;    /* Offset of value */
  double mul;    /* Scale factor applied to value */
} Offscal;

int arrconv(long ndata, Fittype atype, void *adata, Offscal *os,
	    Fittype btype, void *bdata);

int typeconv(long ndata, Fittype atype, void *adata, double zero, double scal,
	     Fittype btype, void *bdata);

int matchstr(const char *sa, const char *sb, int fixlen);
char *rheadline(Fits *fits, Hdu *hdu, int lnum);
int wheadline(Fits *fits, Hdu *hdu, int lnum, char *line);

int get_data(Fits *fits, Hdu *hdu, long offset, Fittype atype, long start,
	     long nobj, Fittype btype, double zero, double scale, Offscal *os,
	     Fitsflag *flags, long blank, void *data);
int put_data(Fits *fits, Hdu *hdu, long offset, Fittype atype, long start,
	     long nobj, Fittype btype, double zero, double scale, Offscal *os,
	     Fitsflag *flags, long blank, void *data);
Fittype dat_type(Hdu *hdu);
int fits_flush(Fits *fits);
int fits_read(Fits *fits, long recnum, int doreport);
int fits_pad(Fits *fits, long recnum);

int w_extkeys(Fits *fits, Hdu *hdu);  /* Write extension name header keywords */

#endif

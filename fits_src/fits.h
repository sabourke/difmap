#ifndef fits_h
#define fits_h

#include <limits.h>

#define FITSLEN (23040/CHAR_BIT)  /* The number of chars per FITS record */

/* Enumerate header types. Each enumerator must be a unique power of two */

typedef enum {
  F_ANY=0,       /* Any type (used only as selection wildcard) */
  F_UNKNOWN=1,   /* Un-recognised HDU extension type */
  F_PRIMARY=2,   /* Primary header */
  F_IMAGE=4,     /* IUE IMAGE extension */
  F_TABLE=8,     /* ASCII table extension */
  F_BINTAB=16    /* Binary table extension */
} Hdutype;

/* Enumerated BITPIX values signaling FITS data types */

typedef enum {
  B_CHAR=8,    /* 8-bit int */
  B_16INT=16,  /* 16 bit IEEE int */
  B_32INT=32,  /* 32 bit IEEE int */
  B_FLOAT=-32, /* 32 bit IEEE floating point */
  B_DBLE=-64   /* 64 bit IEEE double precision */
} Bitpix;

/* Enumerate HDU instantiation status */
 
typedef enum {
  HDU_DESCR,    /* The HDU is currently only a descriptor */
  HDU_HEADER,   /* The HDU file header is being written */
  HDU_DATA,     /* The HDU data-segment is being written */
  HDU_INFILE    /* The HDU is completely specified in its file. */
} Hdustate;

/* Define the members of the Header-Data-Unit base class */

#define HDUBase \
  Hdutype type;     /* Type of HDU described in FITS header */\
  Bitpix bitpix;    /* BITPIX (bits per data pixel) */\
  int naxis;        /* Number of dimensions specified in header */\
  int *dims;        /* Dynamic array of 'naxis' dimensions */\
  int groups;       /* If true - random groups are present */\
  int pcount;       /* Parameter count */\
  int gcount;       /* Group count */\
  int nrec;         /* Number of records in the HDU */\
  int headrec;      /* Start record of header */\
  int datarec;      /* Start record of data */\
  int wnxtline;     /* Next header line to be written */\
  int nextline;     /* Next header line to be read */\
  int endline;      /* Line of END keyword or -1 if not yet known */\
  int pad;          /* The character to be used to pad the data-segment */\
  long grpsize;     /* Number of FITS bytes per group */\
  short complete;   /* True if the descriptor has been fully initialized */\
  Hdustate state;   /* File readiness state enumeration */\
  char *extname;    /* Extension name. */\
  int extver;       /* Extension version number */\
  int extlevel;     /* Extension level in hierachical structure */\
  struct Hdu *next; /* Next HDU in FITS file */

typedef struct Hdu { /* HDU Base-class descriptor */
  HDUBase
} Hdu;

#ifndef recio_h
struct Recio;
#endif

/* The FITS file descriptor */

typedef struct Fits {
  struct Recio *rec;/* Record I/O class object descriptor */
  char *name;       /* Name of FITS file */
  int readonly;     /* True if file has been opened without write access */
  int pedantic;     /* If true be pedantic provide extra warnings */
  int aips;         /* If true write AIPS versions of standard names */
  int modified;     /* True if the data in 'buff' has been modified */
  int complete;     /* True unless the last HDU is incompletely written */
                    /* This is only 0 between add_Hdu() and end_data() */
  int pad;          /* Padding char for records >= nullrec */
  long recnum;      /* Number of current record in buffer 'buff' */
  long nullrec;     /* Index of first un-written record in fits file */
  struct Hdu *hdu;  /* Linked list of Header-Data-Unit descriptors */
  unsigned char buff[FITSLEN]; /* FITS I/O buffer */
} Fits;

/* Enumerate the types retrievable from FITS tables */

typedef enum {
  DAT_NON,   /* void - unknown */
  DAT_SHT,   /* (short) */
  DAT_INT,   /* (int) */
  DAT_LNG,   /* (long) */
  DAT_FLT,   /* (float) */
  DAT_DBL,   /* (double) */
  DAT_CHR,   /* (char) */
  DAT_BYT,   /* (unsigned char) representation of byte. */
  DAT_BIT,   /* (unsigned char) bit array */
  DAT_LOG,   /* (char) representation of FITS logical 'T' or 'F' */
  DAT_SCMP,  /* (float)[2] representation of complex number */
  DAT_DCMP,  /* (double)[2] representation of complex number */
  DAT_COM,   /* (char *) Comment value eg HISTORY or COMMENT keyword values */
  DAT_STR    /* (char *) Terminated string */
} Fittype;

/* Define a magic value to denote that no null value has been defined */

#define NONULL 918273L

/* Define the type used for flag arrays */

typedef char Fitsflag;

/* Public functions for handling HDUs */

Hdu *find_hdu(struct Fits *fits, Hdutype type, char *extname, int extver,
	      Hdu *prev);

Hdu *new_image(Bitpix bitpix, int naxis, int *dims, char *extname, int extver,
	       int extlevel);

Hdu *new_primary(Bitpix bitpix, int naxis, int *dims, int groups, int pcount,
		 int gcount);

Hdu *new_asctab(int width, int nrow, char *extname, int extver, int extlevel,
		int tfields);

Hdu *new_bintab(int nrow, char *extname, int extver, int extlevel, int tfields,
		long heapsize);

int setaxis(Hdu *hdu, int axis, char *ctype, double crpix, double crval,
	    double cdelt, double crota);

int setgroup(Hdu *hdu, int ipar, char *ptype, double pscal, double pzero);

int setimage(Hdu *hdu, double bscale, double bzero, char *bunit, long blank,
	     double datamin, double datamax);

int setprim(Hdu *hdu, char *origin, char *date_obs, char *telescop,
	    char *instrume, char *observer, char *object, char *author,
	    char *referenc, double equinox);

int setafield(Hdu *hdu, int icol, int tbcol, double tscal, double tzero,
	      char *tform, char *tnull, char *ttype, char *tunit);

int setbfield(Hdu *hdu, int icol, double tscal, double tzero,
	      char *tform, long tnull, char *ttype, char *tunit, char *tdisp,
	      char *tdim);

Hdu *del_Hdu(Hdu *hdu);

int add_Hdu(struct Fits *fits, Hdu *hdu);

int end_header(Fits *fits, Hdu *hdu);
int end_data(Fits *fits, Hdu *hdu);

Hdu *dup_Hdu(struct Fits *afits, Hdu *ahdu, struct Fits *bfits);

Hdu *copy_Hdu(Hdu *hdu);


/* Public FITS functions */

Fits *new_Fits(const char *name, int isold, int readonly, int pedantic,
	       int aips);
Fits *del_Fits(Fits *fits);

char *typename(Fittype type);
char *fitsstr(const char *str);

#endif

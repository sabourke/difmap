#ifndef thdu_h
#define thdu_h

#ifndef fits_h
#include "fits.h"
#endif

/* Define the derived members that all table types share */

#define TABBase \
  int tfields;     /* Number of table fields per row */

/* Generic table descriptor */

typedef struct {
  HDUBase
  TABBase
} Thdu;

/* Generic table accessor functions */

#define COL_DIMFN(n) int (n)(Fits *fits, Thdu *thdu, int icol, int irow)
#define COL_TYPEFN(n) Fittype (n)(Thdu *thdu, int icol)
#define COL_VALFN(n) long (n)(Fits *fits, Thdu *thdu, int icol, int irow, \
 Fittype type, int doscale, Fitsflag *flags, int first, int ndata, void *data)
#define COL_SETFN(n) long (n)(Fits *fits, Thdu *thdu, int icol, int irow, \
 Fittype type, int doscale, Fitsflag *flags, int first, int ndata, void *data)
#define COL_FINDFN(n) int (n)(Thdu *thdu, char *ttype, int fixlen)
#define COL_NAMEFN(n) char *(n)(Thdu *thdu, int icol)

Thdu *find_table(Fits *fits, char *extname, int extver, Hdu *prev);
int setdim(Fits *fits, Thdu *thdu, int icol, int irow, int ndata);
int numrow(Thdu *thdu);
int numcol(Thdu *thdu);
int iscolvar(Thdu *thdu, int icol);

COL_FINDFN(find_column);
COL_DIMFN(col_dim);
COL_TYPEFN(col_type);
COL_VALFN(rcolumn);
COL_NAMEFN(col_name);
COL_SETFN(wcolumn);

/* Define a virtual symbol table structure to hold pointers to */
/* table-specific versions of the above functions */

typedef struct {
  COL_VALFN(*valfn);    /* A function to retrieve a column entry value */
  COL_FINDFN(*findfn);  /* A function to seek out columns by name */
  COL_TYPEFN(*typefn);  /* A function to return the data-type of a column */
  COL_DIMFN(*dimfn);    /* A function to return the dimension of a column */
  COL_NAMEFN(*namefn);  /* A function to return the name of a column */
  COL_SETFN(*setfn);    /* A function to set the value of a column entry */
} Tabfn;

#endif

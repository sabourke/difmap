#ifndef ahdu_h
#define ahdu_h

#ifndef thdu_h
#include "thdu.h"
#endif

/* ASCII-table field descriptor */

typedef struct {
  Fittype type;  /* Type of field */
  int tbcol;     /* Column at which field starts */
  double tscal;
  double tzero;
  short width;   /* Width of field in 8-bit bytes */
  char ndec;     /* Number of decimal places */
  char form;     /* Format letter */
  char *tform;
  char *tnull;
  char *ttype;
  char *tunit;
} Afield;

/* ASCII-table extension descriptor */

typedef struct {
  HDUBase
  TABBase
  Afield *fields;  /* Array of 'tfields' field descriptors */
} Ahdu;

#endif

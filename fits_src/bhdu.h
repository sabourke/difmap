#ifndef bhdu_h
#define bhdu_h

#ifndef thdu_h
#include "thdu.h"
#endif

/* Binary-table fields descriptors */

typedef struct {
  Fittype type;  /* Type of field */
  int tbcol;     /* Column at which field starts */
  double tscal;
  double tzero;
  int width;     /* Number of elements across field */
  char isvar;    /* True if the field contains a variable argument list */
  char form;     /* The FORTRAN-90 TFORM character */
  long tnull;
  char *tform;
  char *ttype;
  char *tunit;
  char *tdisp;
  char *tdim;
} Bfield;

/* Binary-table extension descriptor */

typedef struct {
  HDUBase
  TABBase
  Bfield *fields;  /* Array of 'tfields' field descriptors */
  long theap;      /* Heap offset from start of data segment (8-bit bytes) */
  long heap_nxt;   /* Next unused position in heap wrt theap */
} Bhdu;

#endif

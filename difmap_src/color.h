#ifndef color_h
#define color_h

#include "symtab.h"

/* Enumerate color table classifications. */

typedef enum {
  CM_NONE,    /* Don't install a color table */
  CM_COLOR,   /* General color map. If PGIMAG not available recolor PGGREY */
  CM_GREY     /* Grey-scale - use PGIMAG if possible - otherwise use PGGREY */
} Cmclass;

/*
 * Enumerate transfer function types.
 */
typedef enum {
  TR_LINEAR,  /* Linear transfer function */ 
  TR_LOG,     /* Logarithmic transfer function */
  TR_SQRT     /* Square root transfer function */
} Cmtran;     /* Enumeration ID of transfer function */

Cmtran get_Cmtran(char *name);
char *name_Cmtran(Cmtran tran);

/* Color map descriptor */

typedef struct {
  char *name;    /* Name of color table */
  Cmclass class; /* Type classification of color table */
  int nref;      /* Number of references to this color map. If statically */
                 /* allocate intialize to 1 to prevent deallocation */
  int nc;        /* The number of color entries in the table (can be 0) */
  float *l;      /* The normalized brightness to assign each entry to. */
  float *r;      /* Array of 'nc' red intensities */
  float *g;      /* Array of 'nc' green intensities */
  float *b;      /* Array of 'nc' blue intensities */
} Cmap;

/* Color table descriptor */

typedef struct {
  Symtab *symtab; /* Symbol table of color maps */
  Cmap *cmap;     /* Last color-map descriptor returned by get_Cmap() */
  float contra;   /* Contrast of color ramp (normally 1.0). */
  float bright;   /* Brightness at the center color index (normally 0.5). */
  Cmtran tran;    /* Color-map transfer function */
  float vmin;     /* Min data value to be displayed */
  float vmax;     /* Max data value to be displayed */
} Ctable;

Ctable *new_Ctable(void);
Ctable *del_Ctable(Ctable *ctab);
Cmap *new_Cmap(char *name, int nc);
Cmap *del_Cmap(Cmap *cmap);
int add_Cmap(Ctable *ctab, char *name, Cmap *cmap);
Cmap *get_Cmap(Ctable *ctab, char *name);
int recolor(Cmap *cmap, float contra, float bright);

#endif

#ifndef units_h
#define units_h

/* Select map and UVW plane units by the name of a recognized map-plane unit */

int skyunits(char *name);

/* Conversion functions between user and internal units */

double xytorad(double xy);
double radtoxy(double rad);
double uvtowav(double uv);
double wavtouv(double wav);

/* Unit label-type selection enumerator */

typedef enum {
  U_NAME,       /* The official name of the unit */
  U_TLAB,       /* The label to give the units in text */
  U_PLAB        /* The label to give the units in PGPLOT */
} Ultype;

/* Return the labels to refer to units by */

char *mapunits(Ultype ltype);
char *uvwunits(Ultype ltype);

/* Return the two character ordinal suffix of an integer (eg. "th" for 13) */

char *ordinal_suffix(int n);

#endif

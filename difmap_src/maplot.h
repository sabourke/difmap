#ifndef maplot_h
#define maplot_h

#include "color.h"
#include "markerlist.h"

typedef struct {
  float xc;      /* The fractional position of the beam center on the X-axis */
  float yc;      /* The fractional position of the beam center on the Y-axis */
  float minsize; /* The minimum beam size wrt the size of the plot */
  float maxsize; /* The maximum beam size wrt the size of the plot */
} MaplotBeam;

typedef struct {
  float scale;   /* The length in radians of a polarization vector of 1.0 */
                 /*  map intensity units. */
  float icut;    /* The minimum value of the normal intensity map at which */
                 /*  polarization vectors should be drawn. */
  float pcut;    /* The minimum value of the polarization intensity map at */
                 /*  which polarization vectors should be drawn. */
  unsigned dx,dy;/* Vectors will be drawn at every dx'th X-axis cell and */
                 /*  every dy'th Y-axis cell. */
} MaplotVect;

int maplot(Observation *ob, MapBeam *mb, Mapwin *mw, MaplotBeam *mpb,
	   MaplotVect *vect, int domap, Ctable *ctab, int docont, int dovect,
	   int domod, float *levs, int nlev,
	   float cmul, float *box, MarkerList *markers);

#endif

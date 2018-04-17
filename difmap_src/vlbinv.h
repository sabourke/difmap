/* Function prototypes for converting UV data into a map and beam. */
#ifndef vlbinv_h
#define vlbinv_h

#include "mapmem.h"

int uvinvert(Observation *ob, MapBeam *mb, float uvmin, float uvmax,
	     float gauval, float gaurad, int dorad, float errpow,
	     float uvbin);

void uvtrans(MapBeam *mb, int domap);

int optimal_pixel_size(Observation *ob, float uvmin, float uvmax,
		       int nx, int ny, float *xmax, float *ymax);
#endif

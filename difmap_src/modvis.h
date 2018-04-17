#ifndef modvis_h
#define modvis_h

#include "model.h"
#include "obs.h"

void cmpvis(Modcmp *cmp, Subarray *sub, int base, float freq, float uu,
	    float vv, float *amp, float *phs);

void add_cmp_to_modvis(Modcmp *cmp, Subarray *sub, int base, float freq,
		       float uu, float vv, float *re, float *im);

#endif

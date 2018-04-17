#include <stdio.h>
#include <stdlib.h>

#include "obs.h"
#include "logio.h"

/*.......................................................................
 * Change the current weight scale factor.
 * This scale factor change is applied directly to the primary stream
 * in the Observation structure.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation whose weights
 *                      are to be changed.
 *  scale      float    The new weight scale factor. This must be > 0.0f.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int wtscale(Observation *ob, float scale)
{
  Subarray *sub;  /* The descriptor of the sub-array being re-weighted */
  int isub;       /* The index of the sub-array being re-weighted */
  int ut;         /* The index of the integration being re-weighted */
  int base;       /* The index of the visibility being re-weighted */
  int cif;        /* The index of an IF */
  float mult;     /* The incremental change represented by scale */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "wtscale"))
    return 1;
  if(scale<=0.0f) {
    lprintf(stderr, "wtscale: Scale factor must be finite and positive.\n");
    return 1;
  };
/*
 * Determine the mutiplicative change in scale factor required.
 */
  mult = scale / ob->geom.wtscale;
/*
 * Scale the recorded wtscale parameter.
 */
  ob->geom.wtscale *= mult;
/*
 * Re-scale the weights in the visibilities that are currently in the
 * observation structure.
 */
  sub = ob->sub;
  for(isub=0; isub<ob->nsub; isub++,sub++) {
    Integration *integ = sub->integ;
    for(ut=0; ut<sub->ntime; ut++,integ++) {
      Visibility *vis = integ->vis;
      for(base=0; base<sub->nbase; base++,vis++) {
	if(!(vis->bad & FLAG_DEL))
	  vis->wt *= mult;
      };
    };
  };
/*
 * Scale the per-baseline weight sums by the scale factor as well.
 */
  for(sub=ob->sub,isub=0; isub<ob->nsub; isub++,sub++) {
    Baseline *bptr = sub->base;
    for(base=0; base<sub->nbase; base++,bptr++) {
      Baswt *bwt = bptr->bwt;
      for(cif=0; cif<ob->nif; cif++,bwt++)
	bwt->wtsum *= mult;
    };
  };
  return 0;
}

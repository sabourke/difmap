#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "logio.h"
#include "obs.h"
#include "vlbconst.h"
#include "visstat.h"

/*.......................................................................
 * Return the statistics of a given visibility observable, recording the
 * results in a specified container.
 *
 * Input:
 *  ob    Observation *  The parent observation.
 *  qty    VisStatQty    The visibility observable to query.
 *  uvmin       float    The minimum UV radius to take visibilities from.
 *  uvmax       float    The maximum UV radius to take visibilities from.
 * Input/Output:
 *  result    VisStat *  The results will be assigned to *result.
 * Output:
 *  return        int    0 - OK.
 *                       1 - Error.
 */
int ob_vis_stats(Observation *ob, VisStatQty qty, float uvmin, float uvmax,
		 VisStat *result)
{
  int isub;          /* The index of the sub-array being processed */
  int cif;           /* The index of the IF being processed */
  int old_if;        /* Index of current IF to be restored on exit */
  int iter;          /* The iteration number */
  int nvis;          /* The number of visibilities used */
  double sumval;     /* The sum of the values of the observable */
  double meanval;    /* The mean of the values of the observable */
  double sum_sqr_dv; /* The sum of the square difference from the mean */
  double minval;     /* The minimum value of the observable */
  double maxval;     /* The maximum value of the observable */
/*
 * Check the arguments.
 */
  if(!ob_ready(ob, OB_SELECT, "ob_vis_stats"))
    return 1;
  if(!result) {
    lprintf(stderr, "ob_vis_stats: Missing return container.\n");
    return 1;
  };
  switch(qty) {
  case VS_AMP: case VS_PHS: case VS_REAL: case VS_IMAG: case VS_UMAG:
  case VS_VMAG: case VS_UVRAD:
    break;
  default:
    lprintf(stderr, "ob_vis_stats: Unknown observable.\n");
    return 1;
  };
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Initialize the statistics.
 */
  result->nvis = 0;
  result->mean = 0.0;
  result->sigma = 0.0;
  result->scatter = 0.0;
/*
 * Work out the available UV range.
 */
  {
    UVrange *uvr = uvrange(ob, 1, 0, uvmin, uvmax);
    if(uvr==NULL)
      return 1;
    uvmin = uvr->uvrmin;
    uvmax = uvr->uvrmax;
  };
/*
 * Initialize sums.
 */
  sumval = 0.0;
  nvis = 0;
  sum_sqr_dv = 0.0;
  meanval = 0.0;
  minval = 0.0;
  maxval = 0.0;
/*
 * Two iterations will be needed to work out both the mean and the
 * scatter about the mean.
 */
  for(iter=0; iter<2; iter++) {
/*
 * Loop through all sampled IFs.
 */
    for(cif=0; (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
/*
 * Get the next IF.
 */
      if(getIF(ob, cif))
	return 1;
/*
 * Visit each subarray in turn.
 */
      for(isub=0; isub<ob->nsub; isub++) {
	Subarray *sub = &ob->sub[isub];
	int ut;
	for(ut=0; ut<sub->ntime; ut++) {
	  Visibility *vis = sub->integ[ut].vis;
	  int base;
	  for(base=0; base<sub->nbase; base++,vis++) {
/*
 * Get the square of the UV radius.
 */
	    float uu = vis->u * ob->stream.uvscale;
	    float vv = vis->v * ob->stream.uvscale;
	    float uvrad = sqrt(uu*uu+vv*vv);
/*
 * Only look at unflagged visibilities within the requested UV range.
 */
	    if(!vis->bad && uvrad >= uvmin && uvrad <= uvmax) {
	      double val = 0.0;
/*
 * Round the phases into the range -pi to pi.
 */
	      float phase = vis->phs - twopi * floor(vis->phs/twopi+0.5);
/*
 * Compute the phase for the positive U half of the conjugate symmetric
 * UV plane.
 */
	      if(uu < 0.0)
		phase = -phase;
/*
 * Get the observable that is being examined.
 */
	      switch(qty) {
	      case VS_AMP:
		val = vis->amp;
		break;
	      case VS_PHS:
		val = phase;
		break;
	      case VS_REAL:
		val = vis->amp * cos(phase);
		break;
	      case VS_IMAG:
		val = vis->amp * sin(phase);
		break;
	      case VS_UMAG:
		val = fabs(uu);
		break;
	      case VS_VMAG:
		val = fabs(vv);
		break;
	      case VS_UVRAD:
		val = uvrad;
		break;
	      };
/*
 * On the first iteration, work out the sum of the values of the
 * observable, so that we can compute the mean.
 */
	      if(iter==0) {
		sumval += val;
/*
 * Also accumulate the range of the observable.
 */
		if(nvis==0)
		  minval = maxval = val;
		else if(val < minval)
		  minval = val;
		else if(val > maxval)
		  maxval = val;
/*
 * And count the total number of visibilities used in the statistics.
 */
		nvis++;
/*
 * On the second iteration work out the sum of the square deviation
 * from the mean.
 */
	      } else {
		float dv = val - meanval;
		sum_sqr_dv += dv * dv;
	      };
	    };
	  };
	};
      };
    };
/*
 * At the end of the first iteration, compute the mean.
 */
    if(iter==0)
      meanval = nvis>0 ? sumval / nvis : 0;
  };
/*
 * Regard a lack of useable visibilities as an error.
 */
  if(nvis < 1) {
    lprintf(stderr, "ob_vis_stats: There are no useable visibilities.\n");
    return 1;
  };
/*
 * Record the results for return.
 */
  result->nvis = nvis;
  result->mean = meanval;
  result->sigma = sqrt(sum_sqr_dv) / nvis;
  result->scatter = sqrt(sum_sqr_dv / nvis);
  result->minval = minval;
  result->maxval = maxval;
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    return 1;
  return 0;
}

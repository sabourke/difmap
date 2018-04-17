#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "obs.h"

/*.......................................................................
 * Determine the range of UV radii, |U|, and |V| distances, and
 * amplitudes, of visibilities in the current IF and within the UV
 * radius range uvmin -> uvmax (inclusive).
 *
 * Input:
 *  ob     Observation *  The observation to be parameterised.
 *  doall          int    Find the range for all IFs if true. Otherwise
 *                        find the range from just the current stream IF.
 *  dores          int    If true then the amplitude range is that of
 *                        the residuals.
 *  uvmin        float    The UV radius (wavelengths) below which to ignore
 *                        data.
 *  uvmax        float    The UV radius (wavelengths) beyond which to ignore
 *                        data. If the largest of uvmin and uvmax is <= 0.0f
 *                        then the returned limits pertain to the whole data
 *                        set.
 * Output:
 *  return    UVrange *   Pointer to a static internal container, recording
 *                        the maximum UV radius, |U|,|V|, and amplitude of
 *                        Visibilities within UV radius uvmin to uvmax.
 *                        On error, NULL is returned.
 */
UVrange *uvrange(Observation *ob, int doall, int dores,
		 float uvmin, float uvmax)
{
  static UVrange uvr; /* The container of statistics to be returned */
  int ut;        /* The ut being looked at */
  int base;      /* The baseline being looked at */
  int docut;     /* If true, apply uvcut UV radius cutoff */
  int cif;       /* The index of the IF being processed */
  int bif,eif;   /* The indexes of the first and last index to be processed */
  int first = 1; /* True until the first in-range visibility is encounterred */
  int old_if;    /* The index of IF to be restored on return */
/*
 * Check the observation state.
 */
  if(!ob_ready(ob, doall ? OB_SELECT : OB_GETIF, "uvrange"))
    return NULL;
/*
 * Record the index of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Enforce positivity.
 */
  if(uvmin < 0.0f)
    uvmin = 0.0f;
  if(uvmax < 0.0f)
    uvmax = 0.0f;
/*
 * Arrange that uvmin <= uvmax.
 */
  if(uvmin > uvmax) {float ftmp = uvmin; uvmin = uvmax; uvmax = ftmp;};
/*
 * Record whether a UV cutoff range was specified.
 */
  docut = uvmax > 0.0f;
/*
 * Initialize the container.
 */
  uvr.uvrmin = uvr.uvrmax = 0.0f;
  uvr.umin = uvr.umax = 0.0f;
  uvr.vmin = uvr.vmax = 0.0f;
  uvr.ampmin = uvr.ampmax = 0.0f;
  uvr.wtmin = uvr.wtmax = 0.0f;
/*
 * Set the range of IFs to be processed.
 */
  if(doall) {
    bif = 0;
    eif = ob->nif-1;
  } else {
    bif = eif = ob->stream.cif;
  };
/*
 * Loop through all IFs, ending with the default IF.
 */
  for(cif=bif; (cif=nextIF(ob, cif, 1, 1))>=0 && cif<=eif; cif++) {
    float uvscale;  /* UV coordinate scale factor */
    Subarray *sub;  /* The descriptor of the sub-array being processed */
    int isub;       /* The index of sub in ob->sub[] */
/*
 * Get the next IF.
 */
    if(getIF(ob, cif))
      return NULL;
/*
 * Get the uvscale factor for the new IF.
 */
    uvscale = ob->stream.uvscale;
/*
 * Look at all integrations of all sub-arrays, to acertain the min and
 * max values of pertinent visibility stats.
 */
    sub = ob->sub;
    for(isub=0; isub<ob->nsub; isub++, sub++) {
      Integration *integ = sub->integ;
      for(ut=0; ut<sub->ntime; ut++,integ++) {
	Visibility *vis = integ->vis;
	for(base=0; base<sub->nbase; base++,vis++) {
/*
 * Ignore flagged visibilities.
 */
	  if(!vis->bad) {
/*
 * Get the U and V coordinates of the visbility.
 */
	    float uu = vis->u * ob->stream.uvscale;
	    float vv = vis->v * ob->stream.uvscale;
/*
 * Get the UV radius of the visibility.
 */
	    float uvrad = sqrt(uu*uu + vv*vv);
/*
 * Ignore all visibilities except those between uvmin and uvmax.
 */
	    if(!docut || (uvrad >= uvmin && uvrad <= uvmax)) {
	      float amp;   /* The amplitude of the visibility */
	      float wt;    /* The absolute weight of the visibility */
/*
 * Get the amplitude of the visibility.
 */
	      if(dores) {
		float re = vis->amp*cos(vis->phs)-vis->modamp*cos(vis->modphs);
		float im = vis->amp*sin(vis->phs)-vis->modamp*sin(vis->modphs);
		amp = sqrt(re * re + im * im);
	      } else {
		amp = vis->amp;
	      };
/*
 * Get the visibility weight.
 */
	      wt = fabs(vis->wt);
/*
 * Initialize the ranges with the first visibility that lies in range.
 */
	      if(first) {
		first = 0;
		uvr.ampmin = uvr.ampmax = amp;
		uvr.uvrmin = uvr.uvrmax = uvrad;
		uvr.umin = uvr.umax = fabs(uu);
		uvr.vmin = uvr.vmax = fabs(vv);
		uvr.wtmin = uvr.wtmax = wt;
	      } else {
/*
 * Update the amplitude range.
 */
		if(amp > uvr.ampmax)
		  uvr.ampmax = amp;
		if(amp < uvr.ampmin)
		  uvr.ampmin = amp;
/*
 * Update the error range.
 */
		if(wt > uvr.wtmax)
		  uvr.wtmax = wt;
		if(wt < uvr.wtmin)
		  uvr.wtmin = wt;
/*
 * Update the min/max UV radius.
 */
		if(uvrad < uvr.uvrmin)
		  uvr.uvrmin = uvrad;
		else if(uvrad > uvr.uvrmax)
		  uvr.uvrmax = uvrad;
/*
 * Update the min/max U distance.
 */
		if(fabs(uu) < uvr.umin)
		  uvr.umin = fabs(uu);
		else if(fabs(uu) > uvr.umax)
		  uvr.umax = fabs(uu);
/*
 * Update the min/max V distance.
 */
		if(fabs(vv) < uvr.vmin)
		  uvr.vmin = fabs(vv);
		else if(fabs(vv) > uvr.vmax)
		  uvr.vmax = fabs(vv);
	      };
	    };
	  };
	};
      };
    };
  };
/*
 * On Intel processors, performing the same floating point operation
 * with identical operands, in two different functions often return
 * values that differ by a tiny amount. Although the difference is
 * always below the advertized precision of a float, it can be enough
 * to cause a comparison between the two numbers to fail. Since we are
 * computing the limiting UV range based on the UV radii of the
 * visibilities, the UV maximum and minimum equal the radii of one
 * visibility each, and later calculations of these radii, for
 * comparison with these limits are likely to yield values that don't
 * compare equal. Thus, increment the maximum and decrement the
 * minimum by the floating point precision, to ensure the the radius
 * of the outermost visibility in the range will always be <= the
 * computed maximum radius, and the radius of the innermost visibility
 * in the range, will always be >= the minimum range.
 */
  uvr.uvrmax += uvr.uvrmax * FLT_EPSILON;
  uvr.uvrmin -= uvr.uvrmin * FLT_EPSILON;
/*
 * Restore the start IF.
 */
  if(set_cif_state(ob, old_if))
    return NULL;
/*
 * Return the container that holds the results.
 */
  return &uvr;
}


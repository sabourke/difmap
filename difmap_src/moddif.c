#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "obs.h"
#include "logio.h"

/*.......................................................................
 * Determine the goodness of fit and rms deviation between the observed
 * and model visibilities of all IFs.
 *
 * Input:
 *  ob   Observation * The data set to be examined.
 *  uvmin      float   The minimum UV radius to take visibilities from.
 *  uvmax      float   The maximum UV radius to take visibilities from.
 * Input/Output:
 *  md        Moddif * Send a pointer to the container into which to
 *                     place return values.
 * Output:
 *  return       int   0 - OK. You are guaranteed that md->ndata > 0.
 *                     1 - Error. Note that a lack of any useable points
 *                         such as when all points are flagged, is regarded
 *                         as an error.
 */
int moddif(Observation *ob, Moddif *md, float uvmin, float uvmax)
{
  int isub;        /* The index of the sub-array being processed */
  int cif;         /* The index of the IF being processed */
  int old_if;      /* State of current IF to be restored on exit */
  long nvis=0;     /* The number of visibilities used. */
  float msd=0.0f;  /* Mean square diff between model and observed data */
  float chi=0.0f;  /* Mean square number of sigma deviation */
/*
 * Sanity checks.
 */
  if(!ob_ready(ob, OB_SELECT, "moddif"))
    return 1;
  if(md==NULL) {
    lprintf(stderr, "moddif: NULL return container.\n");
    return 1;
  };
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Initialize output values.
 */
  md->ndata = 0;
  md->rms = 0.0f;
  md->chisq = 0.0f;
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
/*
 * Calculate the square modulus of the complex difference vector using the
 * cosine rule.
 */
	    float ampvis = vis->amp;
	    float phsvis = vis->phs;
	    float ampmod = vis->modamp;
	    float phsmod = vis->modphs;
	    float sqrmod = ampvis*ampvis + ampmod*ampmod -
	      2.0f * ampvis*ampmod * cos(phsvis-phsmod);
/*
 * Count the number of visibilities used.
 */
	    nvis++;
/*
 * Accumulate chi-squared.
 * vis->wt is the reciprocal of the ampitude variance.
 */
	    chi += vis->wt * sqrmod;
/*
 * Accumulate the mean-square-difference between model and data.
 */
	    msd += (sqrmod - msd) / nvis;
	  };
	};
      };
    };
  };
/*
 * Set up return values.
 */
  md->ndata = 2L * nvis;      /* Count real + imaginary as two measurements */
  md->rms = sqrt(fabs(msd));
  md->chisq = fabs(chi);
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    return 1;
  return 0;
}

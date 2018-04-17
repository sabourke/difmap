#include <stdio.h>
#include <stdlib.h>
#include "vlbconst.h"
#include "obs.h"
#include "slalib.h"

static int uvmodshift(Observation *ob, float east, float north);

/*.......................................................................
 * Shift the phase center of the observed and model data.
 * The new shifts are recorded as increments to ob->geom.east and
 * ob->geom.north.
 *
 * NB. Whenever a new un-shifted IF is paged into memory the observed
 *     visibilities of that IF must have the accumulated shift applied
 *     to them, but the model must not be shifted again. To shift just
 *     the visibilities use the uvshift() function.
 *
 * Input/Output:
 *  ob  Observation *  The observation to be modified.
 * Input:
 *  east      float    Shift eastwards (radians).
 *  north     float    Shift northwards (radians).
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int obshift(Observation *ob, float east, float north)
{
/*
 * Check the data.
 */
  if(!ob_ready(ob, OB_INDEX, "obshift"))
    return 1;
/*
 * Record the new shifts as increments to those recorded in ob->geom.
 */
  ob->geom.east += east;
  ob->geom.north += north;
/*
 * If an IF is currently in memory, shift its visibilities.
 */
  if(ob_ready(ob, OB_GETIF, NULL) && uvshift(ob, east, north))
    return 1;
/*
 * If there is a selection to maintain model visibilities for
 * shift the UV representation of the model in each IF.
 */
  if(ob_ready(ob, OB_SELECT, NULL) && uvmodshift(ob, east, north))
    return 1;
/*
 * Shift the components of the established model list.
 */
  shiftmod(ob->model, east, north);
/*
 * Shift the components of the tentative model list.
 */
  shiftmod(ob->newmod, east, north);
/*
 * Shift the components of the tentative continuum model list.
 */
  shiftmod(ob->cnewmod, east, north);
/*
 * Shift the components of the established continuum model list.
 */
  shiftmod(ob->cmodel, east, north);
  return 0;
}

/*.......................................................................
 * Undo all accumulated phase-center position shifts.
 *
 * Input:
 *  ob    Observation *  The observation to be unshifted.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int obunshift(Observation *ob)
{
/*
 * Shift back to the phase center.
 */
  if(obshift(ob, -ob->geom.east, -ob->geom.north))
    return 1;
/*
 * Ensure that the shifts are now recorded as exactly zero.
 */
  ob->geom.east = 0.0f;
  ob->geom.north = 0.0f;
  return 0;
}

/*.......................................................................
 * Please see the comments for the obshift() function before using
 * this function. This function only shifts the memory resident
 * observed visibilities.
 *
 * Apply a phase rotation to a UV data set such that the map that results
 * from it will be shifted by (dx,dy) radians.
 *
 * NB. If something is found to be incorrect in this function, note that
 *     the same algorithm is used in primdata() of uvf_write.c to remove
 *     accumulated shifts when writing to the output file. It should also
 *     be fixed.
 *
 * Input/Output:
 *  ob  Observation *  The observation to be phase rotated.
 * Input:
 *  east      float    Shift eastwards (radians).
 *  north     float    Shift northwards (radians).
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int uvshift(Observation *ob, float east, float north)
{
  Subarray *sub;    /* The sub-array being shifted */
  int ut;           /* UT number */
  int base;         /* Baseline number */
  int isub;         /* Subarray number */
/*
 * Note that this function is called by iniIF() so it needs
 * to be callable before IF visibilities are marked as being corrected.
 */
  if(!ob_ready(ob, OB_RAWIF, "uvshift"))
    return 1;
/*
 * We will need to evaluate the fourier component 2.pi.u.dx + 2.pi.v.dy
 * for every U and V in the observation. Pre-multiply 'east' and 'north' with
 * 2.pi times the UV scale factor so that they can directly multiply U and V.
 */
  east *= twopi * ob->stream.uvscale;
  north *= twopi * ob->stream.uvscale;
/*
 * Rotate by adding the pertinent amount to each visibility phase in
 * the observation.
 */
  sub = ob->sub;
  for(isub=0; isub<ob->nsub; isub++,sub++) {
    Integration *integ = sub->integ;
/*
 * Process each integration of the sub-array.
 */
    for(ut=0;ut<sub->ntime; ut++,integ++) {
      Visibility *vis = integ->vis;
      for(base=0; base<sub->nbase;base++,vis++)
	vis->phs += east * vis->u + north * vis->v;
    };
  };
  return 0;
}

/*.......................................................................
 * Apply a phase rotation to all the UV representations of the
 * established model. This entails reading each IF version of the model
 * from the uvmodel.scr paging file, shifting it by the requested
 * amount and writing it back to the paging file.
 *
 * Input/Output:
 *  ob  Observation *  The observation whose model is to be shifted.
 * Input:
 *  east      float    Shift eastwards (radians).
 *  north     float    Shift northwards (radians).
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int uvmodshift(Observation *ob, float east, float north)
{
  Subarray *sub; /* The sub-array being shifted */
  float phsinc;  /* The phase increment to implement the shift */
  float uvscale; /* The UVW coordinate scale factor for a given IF */
  int cif;       /* The index of the IF being processed */
  int ut;        /* UT number */
  int base;      /* Baseline number */
  int isub;      /* Subarray number */
  int old_if;    /* State of current IF to be restored on exit */
/*
 * Ingore this call if there is no selection to shift the model
 * components for.
 */
  if(!ob_ready(ob, OB_SELECT, NULL))
    return 0;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Is there a model to be shifted?
 */
  if(ob->hasmod) {
/*
 * We will need to evaluate the fourier component 2.pi.u.dx + 2.pi.v.dy
 * for every U and V in the observation. Pre-multiply 'east' and 'north' with
 * 2.pi .
 */
    east *= twopi;
    north *= twopi;
/*
 * Loop through all sampled IFs.
 */
    for(cif=0; (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
/*
 * Read the model of the next IF into memory.
 */
      if(getmodel(ob, cif))
	return 1;
/*
 * Get the associated UVW scaling factor.
 */
      uvscale = getuvscale(ob, cif);
/*
 * Rotate by adding the pertinent amount to each visibility phase in
 * the observation.
 */
      sub = ob->sub;
      for(isub=0; isub<ob->nsub; isub++,sub++) {
	Integration *integ = sub->integ;
/*
 * Process each integration of the sub-array.
 */
	for(ut=0;ut<sub->ntime; ut++,integ++) {
	  Visibility *vis = integ->vis;
	  for(base=0; base<sub->nbase;base++,vis++) {
	    phsinc = (east*vis->u + north*vis->v) * uvscale;
	    vis->modphs += phsinc;
	  };
	};
      };
/*
 * Write the modified model back to the uvmodel.scr file.
 */
      if(putmodel(ob, cif))
	return 1;
    };
  };
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    return 1;
  return 0;
}

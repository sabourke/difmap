#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "logio.h"
#include "obs.h"

/*.......................................................................
 * Given an array of 'nbase' visibilities and an optional U,V range,
 * return an array of flags for each of the 'nbase'
 * visibilities, in which usable visibilities are flagged as 1 and
 * un-usable visibilities as 0.
 *
 * Input:
 *  ob   Observation *  The descriptor of the parent observation.
 *  vis   Visibility *  The array of 'nbase' visibilities.
 *  nbase int           The number of baselines in 'vis'.
 *  uvmin float         The min allowed UV radius (Wavelengths).
 *  uvmax float         The max allowed UV radius (Wavelengths).
 *                      NB. the range will only be applied if the greater
 *                      of uvmin and uvmax > 0.0f.
 * Output:
 *  usable int *        Array 'nbase' flags. If usable[base] is 1 then
 *                      vis[base] is a usable visibility.
 *  return int          Number of usable baslines.
 */
int visflags(Observation *ob, Visibility *vis, int nbase,
	     float uvmin, float uvmax, int *usable)
{
  int docut;  /* If true then use UV range in flagging */
  int base;   /* Baseline number. */
  int nuse;   /* Number of baselines that are usable */
  float uvscale; /* The UVW coordinate scale factor */
/*
 * Check that the current stream is valid.
 */
  if(!ob_ready(ob, OB_GETIF, "visflags"))
    return 0;
/*
 * Enforce positivity on uvmin and uvmax.
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
 * See if a UV range for including/excluding data was specified.
 */
  docut = uvmax > 0.0f;
/*
 * Convert the given UV min and max to the units in the visibilities.
 */
  uvscale = ob->stream.uvscale;
  uvmin /= uvscale;
  uvmax /= uvscale;
/*
 * Fill the flag array (visibilities with -ve errors are flagged).
 * Also mark as usable only those visibilities that lie within the
 * the optional UV ellipse described by the user.
 */
  nuse = 0;
  for(base=0; base < nbase; base++,vis++,usable++) {
    *usable = !vis->bad;
/*
 * Exclude points outside the selected UV range.
 */
    if(docut) {
      float uu = vis->u;
      float vv = vis->v;
      float uvrad = sqrt(uu*uu+vv*vv);
      if(uvrad > uvmax || uvrad < uvmin)
	*usable = 0;
    };
/*
 * Count usable baselines.
 */
    if(*usable)
      nuse++;
  };
/*
 * Return the number of usable baselines.
 */
  return nuse;
}


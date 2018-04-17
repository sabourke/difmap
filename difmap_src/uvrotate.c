#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "obs.h"

/*.......................................................................
 * Rotate UV coordinates clockwise in the complex UV plane.
 * The rotation is recorded as an increment to ob->geom.uvangle.
 *
 * Input:
 *  ob  Observation *  The observation holding the visibilities to be
 *                     rotated.
 *  angle     float    The angle of rotation (+ve => clockwise rotation
 *                     of UV coordinates in UV plane. (radians).
 */
void uvrotate(Observation *ob, float angle)
{
/*
 * Check that we have visibilities to rotate.
 */
  if(!ob_ready(ob, OB_INDEX, "uvrotate"))
    return;
/*
 * Increment the recorded rotation.
 */
  ob->geom.uvangle += angle;
/*
 * If there is currently an IF in memory, rotate its visibilities.
 */
  if(ob_ready(ob, OB_RAWIF, NULL)) {
    Subarray *sub;   /* The descriptor of the sub-array being processed */
    int isub;        /* The index of the sub-array being processed */
    int base;        /* The baseline number being processed */
    int ut;          /* The ut number being processed */
/*
 * Precompute sin() and cos().
 */
    float sin_ang = sin(angle);
    float cos_ang = cos(angle);
/*
 * Perform the rotation.
 */
    sub = ob->sub;
    for(isub=0; isub<ob->nsub; isub++,sub++) {
      Integration *integ = sub->integ;
      for(ut=0; ut<sub->ntime; ut++,integ++) {
	Visibility *vis = integ->vis;
	for(base=0; base<sub->nbase; base++,vis++) {
	  float u = vis->u;
	  float v = vis->v;
	  vis->u = u * cos_ang + v * sin_ang;
	  vis->v = v * cos_ang - u * sin_ang;
	};
      };
    };
  };
  return;
}

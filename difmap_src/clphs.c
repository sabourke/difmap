#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "logio.h"
#include "obs.h"
#include "vlbconst.h"

#include "telspec.h"
#include "clphs.h"

/*.......................................................................
 * Construct the closure phase described in a given closure-phase-triangle
 * descriptor for the array of visibilities of a single integration.
 *
 * The returned closure phase is in the range -pi to pi.
 *
 * Input:
 *  ts     Trispec *  The descriptor of the closure-phase-triangle to
 *                    construct the phase for. 
 *  vis Visibility *  The array of visibilities of an integration
 *                    from sub-array ts->isub.
 * Output:
 *  return   Clphs *  The pointer to the static internal descriptor in
 *                    which the closure phase is described. This is
 *                    overwritten on each call.
 */
Clphs *get_clphs(Trispec *ts, Visibility *vis)
{
  static Clphs cp;   /* The return container */
  float sumvar=0.0f; /* Sum of visibility phase variances */
  int i;
/*
 * Pre-initialize the return container.
 */
  cp.wt = cp.ophs = cp.mphs = 0.0f;
  cp.bad = 0;
/*
 * Accumulate the closure phases and weight.
 */
  for(i=0; i<3; i++) {
    Visibility *v = &vis[ts->b[i].base];
    float sign = ts->b[i].sign;
    float wt = v->wt;
    float amp = v->amp;
/*
 * Accumulate the observed and model closure-phases.
 */
    cp.ophs += sign * v->phs;
    cp.mphs += sign * v->modphs;
/*
 * The visibility phase variance is equal to the amplitude variance
 * divided by the amplitude squared. So with the normal assumption that the
 * visibility weight is the reciprocal of the amplitude variance, the
 * variance of the visibility phase is the reciprocal of the
 * visibility weight times the reciprocal of the amplitude squared.
 * The variance of the closure phase is equal to the the individual
 * visibility phase variances.
 */
    if((v->bad & FLAG_DEL) || amp == 0.0f || wt == 0.0f) {
      cp.bad |= FLAG_CDEL;
    } else {
      sumvar += 1.0 / (wt * amp * amp);
      if(v->bad)
	cp.bad |= (v->bad & FLAG_BAD) ? FLAG_CBAD : FLAG_CCOR;
    };
  };
/*
 * Assign zero weight to closure phases that contain deleted visibilities.
 */
  if((cp.bad & FLAG_CDEL) || sumvar <= 0.0f)
    cp.wt = 0.0f;
  else
    cp.wt = 1.0/sumvar;
/*
 * The closure phase is known modulo 2.pi radians. Wrap the observed phase
 * into the range -pi to pi.
 */
  cp.ophs = fmod(cp.ophs, twopi);
  if(cp.ophs > pi)
    cp.ophs -= twopi;
  else if(cp.ophs < -pi)
    cp.ophs += twopi;
/*
 * Do the same for the model phase.
 */
  cp.mphs = fmod(cp.mphs, twopi);
  if(cp.mphs > pi)
    cp.mphs -= twopi;
  else if(cp.mphs < -pi)
    cp.mphs += twopi;
/*
 * Return the closure phase descriptor.
 */
  return &cp;
}


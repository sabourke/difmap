#include <math.h>
#include <stdio.h>
#include "vlbutil.h"

/*.......................................................................
 * Given two amplitude,phase pairs, subtract them in the complex plane
 * and return the result as amplitude and phase.
 *
 * Input/Output:
 *  amp   float *   The amplitude of the coordinate pair to be subtracted
 *                  from.
 *  phs   float *   The phase of the coordinate pair to be subtracted
 *                  from.
 * Input:
 *  subamp float    The amplitude of the coordinate pair to be subtracted.
 *  subphs float    The phase of the coordinate pair to be subtracted.
 */
void subamphs(float *amp, float *phs, float subamp, float subphs)
{
  static float x;  /* Real part of complex difference */
  static float y;  /* Imaginary part of complex difference */
/*
 * Subtract the complex versions of the polar numbers.
 */
  x = *amp * cos(*phs) - subamp * cos(subphs);
  y = *amp * sin(*phs) - subamp * sin(subphs);
/*
 * Convert back to amplitude and phase.
 */
  *amp = sqrt(x*x+y*y);
  *phs = (x==0.0f && y==0.0f) ? 0.0f : atan2(y,x);
  return;
}

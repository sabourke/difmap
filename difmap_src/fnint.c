#include "vlbmath.h"
/*.......................................................................
 * Round the float argument to the nearest integer.
 *
 * Input:
 *  fval  float  The float value to be rounded.
 * Output:
 *  return int   The rounded integer version of fval. If fval >= 0.0
 *               return int(fval+0.5), else if fval < 0.0 return
 *               int(fval-0.5).
 */
int fnint(float fval)
{
  return (fval >= 0.0) ? fval+0.5 : fval-0.5;
}

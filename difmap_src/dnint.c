#include "vlbmath.h"
/*.......................................................................
 * Round the double argument to the nearest integer.
 *
 * Input:
 *  dval  double The double value to be rounded.
 * Output:
 *  return int   The rounded integer version of dval. If dval >= 0.0
 *               return int(dval+0.5), else if dval < 0.0 return
 *               int(dval-0.5).
 */
int dnint(double dval)
{
  return (dval >= 0.0) ? dval+0.5 : dval-0.5;
}

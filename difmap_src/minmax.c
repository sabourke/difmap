#include "vlbmath.h"
/*.......................................................................
 * Return the max float argument.
 * Input:
 *  a   float    First of possible return values.
 *  b   float    Second of possible return values.
 * Output:
 *  return float The maximum of a and b.
 */
float floatmax(float a, float b)
{
  return (a>b) ? a:b;
}

/*.......................................................................
 * Return the max int argument.
 * Input:
 *  a   int    First of possible return values.
 *  b   int    Second of possible return values.
 * Output:
 *  return int The maximum of a and b.
 */
int imax(int a, int b)
{
  return (a>b) ? a:b;
}

/*.......................................................................
 * Return the max double argument.
 * Input:
 *  a   double    First of possible return values.
 *  b   double    Second of possible return values.
 * Output:
 *  return double The maximum of a and b.
 */
double dmax(double a, double b)
{
  return (a>b) ? a:b;
}


/*.......................................................................
 * Return the min float argument.
 * Input:
 *  a   float    First of possible return values.
 *  b   float    Second of possible return values.
 * Output:
 *  return float The minimum of a and b.
 */
float floatmin(float a, float b)
{
  return (a<b) ? a:b;
}

/*.......................................................................
 * Return the min int argument.
 * Input:
 *  a   int    First of possible return values.
 *  b   int    Second of possible return values.
 * Output:
 *  return int The minimum of a and b.
 */
int imin(int a, int b)
{
  return (a<b) ? a:b;
}

/*.......................................................................
 * Return the min double argument.
 * Input:
 *  a   double    First of possible return values.
 *  b   double    Second of possible return values.
 * Output:
 *  return double The minimum of a and b.
 */
double dmin(double a, double b)
{
  return (a<b) ? a:b;
}


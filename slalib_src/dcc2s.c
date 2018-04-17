#include "slalib.h"
#include "slamac.h"
void slaDcc2s ( double v[3], double *a, double *b )
/*
**  - - - - - - - - -
**   s l a D c c 2 s
**  - - - - - - - - -
**
**  Direction cosines to spherical coordinates.
**
**  (double precision)
**
**  Given:
**     v      double[3]   x,y,z vector
**
**  Returned:
**     *a,*b  double      spherical coordinates in radians
**
**  The spherical coordinates are longitude (+ve anticlockwise
**  looking from the +ve latitude pole) and latitude.  The
**  Cartesian coordinates are right handed, with the x axis
**  at zero longitude and latitude, and the z axis at the
**  +ve latitude pole.
**
**  If v is null, zero a and b are returned.
**  At either pole, zero a is returned.
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double x, y, z, r;

   x = v[0];
   y = v[1];
   z = v[2];
   r = sqrt ( x * x + y * y );

   *a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
   *b = ( z != 0.0 ) ? atan2 ( z, r ) : 0.0;
}

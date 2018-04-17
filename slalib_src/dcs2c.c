#include "slalib.h"
#include "slamac.h"
void slaDcs2c ( double a, double b, double v[3] )
/*
**  - - - - - - - - -
**   s l a D c s 2 c
**  - - - - - - - - -
**
**  Spherical coordinates to direction cosines.
**
**  (double precision)
**
**  Given:
**     a,b       double      spherical coordinates in radians
**                           (RA,Dec), (long,lat) etc
**
**  Returned:
**     v         double[3]   x,y,z unit vector
**
**  The spherical coordinates are longitude (+ve anticlockwise
**  looking from the +ve latitude pole) and latitude.  The
**  Cartesian coordinates are right handed, with the x axis
**  at zero longitude and latitude, and the z axis at the
**  +ve latitude pole.
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double cosb;

   cosb = cos ( b );
   v[0] = cos ( a ) * cosb;
   v[1] = sin ( a ) * cosb;
   v[2] = sin ( b );
}

#include "slalib.h"
#include "slamac.h"
void slaGeoc ( double p, double h, double *r, double *z )
/*
**  - - - - - - - -
**   s l a G e o c
**  - - - - - - - -
**
**  Convert geodetic position to geocentric.
**
**  (double precision)
**
**  Given:
**     p     double     latitude (geodetic, radians)
**     h     double     height above reference spheroid (geodetic, metres)
**
**  Returned:
**     *r    double     distance from Earth axis (AU)
**     *z    double     distance from plane of Earth equator (AU)
**
**  Notes:
**
**     1)  Geocentric latitude can be obtained by evaluating atan2(z,r).
**
**     2)  IAU 1976 constants are used.
**
**  Reference:
**     Green,R.M., Spherical Astronomy, CUP 1985, p98.
**
**  Last revision:   25 July 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double sp, cp, c, s;

/* Earth equatorial radius (metres) */
   static double a0 = 6378140.0;

/* Reference spheroid flattening factor and useful function thereof */
   static double f = 1.0 / 298.257;
   double b = ( 1.0 - f ) * ( 1.0 - f );

/* Astronomical unit in metres */
   static double au = 1.49597870e11;

/* Geodetic to geocentric conversion */
   sp = sin ( p );
   cp = cos ( p );
   c = 1.0 / sqrt ( cp * cp + b * sp * sp );
   s = b * c;
   *r = ( a0 * c + h ) * cp / au;
   *z = ( a0 * s + h ) * sp / au;
}

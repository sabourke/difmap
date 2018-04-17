#include "slalib.h"
#include "slamac.h"
double slaGmst ( double ut1 )
/*
**  - - - - - - - -
**   s l a G m s t
**  - - - - - - - -
**
**  Conversion from Universal Time to Sidereal Time.
**
**  (double precision)
**
**  Given:
**    ut1    double     Universal Time (strictly UT1) expressed as
**                      Modified Julian Date (JD-2400000.5)
**
**  The result is the Greenwich Mean Sidereal Time (double
**  precision, radians).
**
**  The IAU 1982 expression (see page S15 of the 1984 Astronomical
**  Almanac) is used, but rearranged to reduce rounding errors.
**  This expression is always described as giving the GMST at
**  0 hours UT.  In fact, it gives the difference between the
**  GMST and the UT, which happens to equal the GMST (modulo
**  24 hours) at 0 hours UT each day.  In this routine, the
**  entire UT is used directly as the argument for the
**  standard formula, and the fractional part of the UT is
**  added separately;  note that the factor 1.0027379... does
**  not appear.
**
**  See also the routine sla_GMSTA, which delivers better numerical
**  precision by accepting the UT date and time as separate arguments.
**
**  Called:  slaDranrm
**
**  Defined in slamac.h:  D2PI, DS2R, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double tu;

/* Julian centuries from fundamental epoch J2000 to this UT */
   tu = ( ut1 - 51544.5 ) / 36525.0;

/* GMST at this UT */
   return slaDranrm ( dmod ( ut1, 1.0 ) * D2PI +
                       ( 24110.54841 +
                       ( 8640184.812866 +
                       ( 0.093104 - 6.2e-6 * tu ) * tu ) * tu ) * DS2R );
}

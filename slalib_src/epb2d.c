#include "slalib.h"
#include "slamac.h"
double slaEpb2d ( double epb )
/*
**  - - - - - - - - -
**   s l a E p b 2 d
**  - - - - - - - - -
**
**  Conversion of Besselian epoch to Modified Julian Date.
**
**  (double precision)
**
**  Given:
**     epb      double       Besselian epoch
**
**  The result is the Modified Julian Date (JD - 2400000.5).
**
**  Reference:
**     Lieske,J.H., 1979. Astron. Astrophys.,73,282.
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   return 15019.81352 + ( epb - 1900.0 ) * 365.242198781;
}

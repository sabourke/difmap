#include "slalib.h"
#include "slamac.h"
double slaEpj ( double date )
/*
**  - - - - - - -
**   s l a E p j
**  - - - - - - -
**
**  Conversion of Modified Julian Date to Julian epoch.
**
**  (double precision)
**
**  Given:
**     date     double      Modified Julian Date (JD - 2400000.5)
**
**  The result is the Julian epoch.
**
**  Reference:
**     Lieske,J.H., 1979. Astron. Astrophys.,73,282.
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
  return 2000.0 + ( date - 51544.5 ) / 365.25;
}

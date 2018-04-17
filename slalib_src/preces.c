#include "slalib.h"
#include "slamac.h"
#include <ctype.h>
void slaPreces ( char sys[3], double ep0, double ep1,
                 double *ra, double *dc )
/*
**  - - - - - - - - - -
**   s l a P r e c e s
**  - - - - - - - - - -
**
**  Precession - either FK4 (Bessel-Newcomb, pre-IAU1976) or
**  FK5 (Fricke, post-IAU1976) as required.
**
**  Given:
**     sys        char[]     precession to be applied: "FK4" or "FK5"
**     ep0,ep1    double     starting and ending epoch
**     ra,dc      double     RA,Dec, mean equator & equinox of epoch ep0
**
**  Returned:
**     *ra,*dc    double     RA,Dec, mean equator & equinox of epoch ep1
**
**  Called:    slaDranrm, slaPrebn, slaPrec, slaDcs2c,
**             slaDmxv, slaDcc2s
**
**  Notes:
**
**  1)  The epochs are Besselian if sys='FK4' and Julian if 'FK5'.
**      For example, to precess coordinates in the old system from
**      equinox 1900.0 to 1950.0 the call would be:
**          slaPreces ( "FK4", 1900.0, 1950.0, &ra, &dc )
**
**  2)  This routine will not correctly convert between the old and
**      the new systems - for example conversion from B1950 to J2000.
**      For these purposes see slaFk425, slaFk524, slaFk45z and
**      slaFk54z.
**
**  3)  If an invalid sys is supplied, values of -99.0,-99.0 will
**      be returned for both ra and dc.
**
**  Last revision:   22 December 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double pm[3][3], v1[3], v2[3];

/* Validate sys */
   if ( ( toupper ( sys[0] ) != 'F' )
     || ( toupper ( sys[1] ) != 'K' )
     || ( sys[2] != '4' && sys[2] != '5' ) ) {
         *ra = -99.0;          /* Error */
         *dc = -99.0;
   } else {

   /* Generate appropriate precession matrix */
      if ( sys[2] == '4' )
         slaPrebn ( ep0, ep1, pm );
      else
         slaPrec ( ep0, ep1, pm );

   /* Convert RA,Dec to x,y,z */
      slaDcs2c ( *ra, *dc, v1 );

   /* Precess */
      slaDmxv ( pm, v1, v2 );

   /* Back to RA,Dec */
      slaDcc2s ( v2, ra, dc );
      *ra = slaDranrm ( *ra );
   }
}

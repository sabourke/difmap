#include "slalib.h"
#include "slamac.h"
void slaMapqk ( double rm, double dm, double pr, double pd,
                double px, double rv, double amprms[21],
                double *ra, double *da )
/*
**  - - - - - - - - -
**   s l a M a p q k
**  - - - - - - - - -
**
**  Quick mean to apparent place:  transform a star RA,Dec from
**  mean place to geocentric apparent place, given the
**  star-independent parameters.
**
**  Use of this routine is appropriate when efficiency is important
**  and where many star positions, all referred to the same equator
**  and equinox, are to be transformed for one epoch.  The
**  star-independent parameters can be obtained by calling the
**  slaMappa routine.
**
**  If the parallax and proper motions are zero the slaMapqkz
**  routine can be used instead.
**
**  The reference frames and timescales used are post IAU 1976.
**
**  Given:
**     rm,dm    double      mean RA,Dec (rad)
**     pr,pd    double      proper motions:  RA,Dec changes per Julian year
**     px       double      parallax (arcsec)
**     rv       double      radial velocity (km/sec, +ve if receding)
**
**     amprms   double[21]  star-independent mean-to-apparent parameters:
**
**       (0)      time interval for proper motion (Julian years)
**       (1-3)    barycentric position of the Earth (AU)
**       (4-6)    heliocentric direction of the Earth (unit vector)
**       (7)      (grav rad Sun)*2/(Sun-Earth distance)
**       (8-10)   barycentric Earth velocity in units of c
**       (11)     sqrt(1-v**2) where v=modulus(abv)
**       (12-20)  precession/nutation (3,3) matrix
**
**  Returned:
**     *ra,*da  double      apparent RA,Dec (rad)
**
**  References:
**     1984 Astronomical Almanac, pp B39-B41.
**     (also Lederle & Schwan, Astron. Astrophys. 134, 1-6, 1984)
**
**  Notes:
**
**    1)  The vectors amprms(1-3) and amprms(4-6) are referred to
**        the mean equinox and equator of epoch eq.
**
**    2)  Strictly speaking, the routine is not valid for solar-system
**        sources, though the error will usually be extremely small.
**        However, to prevent gross errors in the case where the
**        position of the Sun is specified, the gravitational
**        deflection term is restrained within about 920 arcsec of the
**        centre of the Sun's disc.  The term has a maximum value of
**        about 1.85 arcsec at this radius, and decreases to zero as
**        the centre of the disc is approached.
**
**  Called:
**     slaDcs2c       spherical to Cartesian
**     slaDvdv        dot product
**     slaDmxv        matrix x vector
**     slaDcc2s       Cartesian to spherical
**     slaDranrm      normalize angle 0-2pi
**
**  Defined in slamac.h:  DAS2R
**
**  Last revision:   21 July 1994
**
**  Copyright P.T.Wallace.  All rights reserved.
*/

#define VF 0.21094502     /* Km/s to AU/year */

{
   int i;
   double pmt, gr2e, ab1, eb[3], ehn[3], abv[3],
          q[3], pxr, w, em[3], p[3], pn[3], pde, pdep1,
          p1[3], p1dv, p1dvp1, p2[3], p3[3];

/* Unpack scalar and vector parameters */
   pmt = amprms[0];
   gr2e = amprms[7];
   ab1 = amprms[11];
   for ( i = 0; i < 3; i++ )
   {
      eb[i] = amprms[i+1];
      ehn[i] = amprms[i+4];
      abv[i] = amprms[i+8];
   }

/* Spherical to x,y,z */
   slaDcs2c ( rm, dm, q );

/* Space motion (radians per year) */
   pxr = px * DAS2R;
   w = VF * rv * pxr;
   em[0] = (-pr * q[1]) - ( pd * cos ( rm ) * sin ( dm ) ) + ( w * q[0] );
   em[1] = ( pr * q[0]) - ( pd * sin ( rm ) * sin ( dm ) ) + ( w * q[1] );
   em[2] =                ( pd * cos ( dm )              ) + ( w * q[2] );

/* Geocentric direction of star (normalized) */
   for ( i = 0; i < 3; i++ ) {
      p[i] = q[i] + ( pmt * em[i] ) - ( pxr * eb[i] );
   }
   slaDvn ( p, pn, &w );

/* Light deflection (restrained within the Sun's disc) */
   pde = slaDvdv ( pn, ehn );
   pdep1 = 1.0 + pde;
   w = gr2e / gmax ( pdep1, 1.0e-5 );
   for ( i = 0; i < 3; i++ ) {
      p1[i] = pn[i] + ( w * ( ehn[i] - pde * pn[i] ) );
   }

/* Aberration */
   p1dv = slaDvdv ( p1, abv );
   p1dvp1 = p1dv + 1.0;
   w = 1.0 + p1dv / ( ab1 + 1.0 );
   for ( i = 0; i < 3; i++ ) {
      p2[i] = ( ab1 * p1[i] + w * abv[i] ) / p1dvp1;
   }

/* Precession and nutation */
   slaDmxv ( (double(*)[3]) &amprms[12], p2, p3 );

/* Geocentric apparent RA,dec */
   slaDcc2s ( p3, ra, da );

   *ra = slaDranrm ( *ra );
}

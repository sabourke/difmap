#include "slalib.h"
#include "slamac.h"
void slaFk45z ( double r1950, double d1950, double bepoch,
                double *r2000, double *d2000 )
/*
**  - - - - - - - - -
**   s l a F k 4 5 z
**  - - - - - - - - -
**
**  Convert B1950.0 FK4 star data to J2000.0 FK5 assuming zero
**  proper motion in an inertial frame (double precision)
**
**  This routine converts stars from the old, Bessel-Newcomb, FK4
**  system to the new, IAU 1976, FK5, Fricke system, in such a
**  way that the FK5 proper motion is zero.  Because such a star
**  has, in general, a non-zero proper motion in the FK4 system,
**  the routine requires the epoch at which the position in the
**  FK4 system was determined.
**
**  The method is from Appendix 2 of Ref 1, but using the constants
**  of Ref 4.
**
**  Given:
**     r1950,d1950     double   B1950.0 FK4 RA,Dec at epoch (rad)
**     bepoch          double   Besselian epoch (e.g. 1979.3)
**
**  Returned:
**     *r2000,*d2000   double   J2000.0 FK5 RA,Dec (rad)
**
**  Notes:
**
**  1)  The epoch BEPOCH is strictly speaking Besselian, but
**      if a Julian epoch is supplied the result will be
**      affected only to a negligible extent.
**
**  2)  Conversion from Besselian epoch 1950.0 to Julian epoch
**      2000.0 only is provided for.  Conversions involving other
**      epochs will require use of the appropriate precession,
**      proper motion, and E-terms routines before and/or
**      after FK45Z is called.
**
**  3)  In the FK4 catalogue the proper motions of stars within
**      10 degrees of the poles do not embody the differential
**      E-term effect and should, strictly speaking, be handled
**      in a different manner from stars outside these regions.
**      However, given the general lack of homogeneity of the star
**      data available for routine astrometry, the difficulties of
**      handling positions that may have been determined from
**      astrometric fields spanning the polar and non-polar regions,
**      the likelihood that the differential E-terms effect was not
**      taken into account when allowing for proper motion in past
**      astrometry, and the undesirability of a discontinuity in
**      the algorithm, the decision has been made in this routine to
**      include the effect of differential E-terms on the proper
**      motions for all stars, whether polar or not.  At epoch 2000,
**      and measuring on the sky rather than in terms of dRA, the
**      errors resulting from this simplification are less than
**      1 milliarcsecond in position and 1 milliarcsecond per
**      century in proper motion.
**
**  References:
**
**     1  Aoki,S., et al, 1983.  Astron. Astrophys., 128, 263.
**
**     2  Smith, C.A. et al, 1989.  "The transformation of astrometric
**        catalog systems to the equinox J2000.0".  Astron.J. 97, 265.
**
**     3  Yallop, B.D. et al, 1989.  "Transformation of mean star places
**        from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
**        Astron.J. 97, 274.
**
**     4  Seidelmann, P.K. (ed), 1992.  "Explanatory Supplement to
**        the Astronomical Almanac", ISBN 0-935702-68-7.
**
**  Called:  slaDcs2c, slaEpj, slaEpb2d, slaDcc2s, slaDranrm
**
**  Defined in slamac.h:  D2PI
**
**  Last revision:   8 November 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double w;
   int i, j;

/* Position and position+velocity vectors */
   double r0[3], a1[3], v1[3], v2[6];

/* Radians per year to arcsec per century */
   static double pmf = 100.0 * 60.0 * 60.0 * 360.0 / D2PI;

/*
** Canonical constants  (see references)
*/

/* vectors a and adot, and matrix m (only half of which is needed here) */
   static double a[3]  = { -1.62557e-6,  -0.31919e-6, -0.13843e-6 };
   static double ad[3] = {  1.245e-3,    -1.580e-3,   -0.659e-3 };
   static double em[6][3] =
   {
     {  0.9999256782, -0.0111820611, -0.0048579477 },
     {  0.0111820610,  0.9999374784, -0.0000271765 },
     {  0.0048579479, -0.0000271474,  0.9999881997 },
     { -0.000551,     -0.238565,      0.435739     },
     {  0.238514,     -0.002667,     -0.008541     },
     { -0.435623,      0.012254,      0.002117     }
   };

/* Spherical to Cartesian */
   slaDcs2c ( r1950, d1950, r0 );

/* Adjust vector a to give zero proper motion in FK5 */
   w = ( bepoch - 1950.0 ) / pmf;
   for ( i = 0; i < 3; i++ ) {
      a1[i] = a[i] + w * ad[i];
   }

/* Remove e-terms */
   w = r0[0] * a1[0] + r0[1] * a1[1] + r0[2] * a1[2];
   for ( i = 0; i < 3; i++ ) {
      v1[i] = r0[i] - a1[i] + w * r0[i];
   }

/* Convert position vector to Fricke system */
   for ( i = 0; i < 6; i++ ) {
      w = 0.0;
      for ( j = 0; j < 3; j++ ) {
         w += em[i][j] * v1[j];
      }
      v2[i] = w;
   }

/* Allow for fictitious proper motion in FK4 */
   w = ( slaEpj ( slaEpb2d ( bepoch ) ) - 2000.0 ) / pmf;
   for ( i = 0; i < 3; i++ ) {
      v2[i] += w * v2[i+3];
   }

/* Revert to spherical coordinates */
   slaDcc2s ( v2, &w, d2000 );
   *r2000 = slaDranrm ( w );
}

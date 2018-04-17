#include "slalib.h"
#include "slamac.h"
void slaPrec ( double ep0, double ep1, double rmatp[3][3] )
/*
**  - - - - - - - -
**   s l a P r e c
**  - - - - - - - -
**
**  Form the matrix of precession between two epochs (IAU 1976, FK5).
**
**  (double precision)
**
**  Given:
**     ep0    double         beginning epoch
**     ep1    double         ending epoch
**
**  Returned:
**     rmatp  double[3][3]   precession matrix
**
**  Notes:
**
**  1)  The epochs are TDB (loosely ET) Julian epochs.
**
**  2)  The matrix is in the sense   v(ep1)  =  rmatp * v(ep0) .
**
**  3)  Though the matrix method itself is rigorous, the precession
**      angles are expressed through canonical polynomials which are
**      valid only for a limited time span.  There are also known
**      errors in the IAU precession rate.  The absolute accuracy
**      of the present formulation is better than 0.1 arcsec from
**      1960AD to 2040AD, better than 1 arcsec from 1640AD to 2360AD,
**      and remains below 3 arcsec for the whole of the period
**      500BC to 3000AD.  The errors exceed 10 arcsec outside the
**      range 1200BC to 3900AD, exceed 100 arcsec outside 4200BC to
**      5600AD and exceed 1000 arcsec outside 6800BC to 8200AD.
**      The SLALIB routine slaPrecl implements a more elaborate
**      model which is suitable for problems spanning several
**      thousand years.
**
**  References:
**     Lieske,J.H., 1979. Astron. Astrophys.,73,282.
**          equations (6) & (7), p283.
**     Kaplan,G.H., 1981. USNO circular no. 163, pa2.
**
**  Called:  slaDeuler
**
**  Defined in slamac.h:  DAS2R
**
**  Last revision:   10 July 1994
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double t0, t, tas2r, w, zeta, z, theta;

/* Interval between basic epoch J2000.0 and beginning epoch (JC) */
   t0 = ( ep0 - 2000.0 ) / 100.0;

/* Interval over which precession required (JC) */
   t =  ( ep1 - ep0 ) / 100.0;

/* Euler angles */
   tas2r = t * DAS2R;
   w = 2306.2181 + ( ( 1.39656 - ( 0.000139 * t0 ) ) * t0 );
   zeta = (w + ( ( 0.30188 - 0.000344 * t0 ) + 0.017998 * t ) * t ) * tas2r;
   z = (w + ( ( 1.09468 + 0.000066 * t0 ) + 0.018203 * t ) * t ) * tas2r;
   theta = ( ( 2004.3109 + ( - 0.85330 - 0.000217 * t0 ) * t0 )
          + ( ( -0.42665 - 0.000217 * t0 ) - 0.041833 * t ) * t ) * tas2r;

/* Rotation matrix */
   slaDeuler ( "ZYZ", -zeta, theta, -z, rmatp );
}

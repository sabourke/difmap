#include "slalib.h"
#include "slamac.h"
void slaNut ( double date, double rmatn[3][3] )
/*
**  - - - - - - -
**  s l a N u t
**  - - - - - - -
**
**  Form the matrix of nutation for a given date (IAU 1980 theory).
**
**  (double precision)
**
**  References:
**     Final report of the IAU working group on nutation,
**        chairman P.K.Seidelmann, 1980.
**     Kaplan, G.H., 1981, USNO circular no. 163, pA3-6.
**
**  Given:
**     date   double        TDB (loosely ET) as Modified Julian Date
**                                           (=JD-2400000.5)
**
**  Returned:
**     rmatn  double[3][3]  nutation matrix
**
**  The matrix is in the sense   v(true)  =  rmatn * v(mean) .
**
**  Called:   slaNutc, slaDeuler
**
**  Last revision:   14 July 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double dpsi, deps, eps0;

/* Nutation components and mean obliquity */
   slaNutc ( date, &dpsi, &deps, &eps0 );

/* Rotation matrix */
   slaDeuler ( "xzx", eps0, -dpsi, - ( eps0 + deps ), rmatn );
}

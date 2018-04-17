#include "slalib.h"
#include "slamac.h"
double slaEqeqx ( double date )
/*
**  - - - - - - - - -
**   s l a E q e q x
**  - - - - - - - - -
**
**  Equation of the equinoxes (IAU 1994, double precision).
**
**  Given:
**     date    double      TDB (loosely ET) as Modified Julian Date
**                                          (JD-2400000.5)
**
**  The result is the equation of the equinoxes (double precision)
**  in radians:
**
**  Greenwich apparent ST = Greenwich mean ST + equation of the equinoxes
**
**  References:  IAU Resolution C7, Recommendation 3 (1994)
**               Capitaine, N. & Gontier, A.-M., Astron. Astrophys.,
**               275, 645-650 (1993)
**
**  Called:  slaNutc
**
**  Last revision:   21 November 1994
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
#define T2AS 1296000.0                 /* Turns to arc seconds */
#define AS2R 0.4848136811095359949E-5  /* Arc seconds to radians */
{
   double t, om, dpsi, deps, eps0;

/* Interval between basic epoch J2000.0 and current epoch (JC) */
   t = ( date - 51544.5 ) / 36525.0;

/* Longitude of the mean ascending node of the lunar orbit on the
   ecliptic, measured from the mean equinox of date */
   om = AS2R * ( 450160.280 + ( -5.0 * T2AS - 482890.539
                               + ( 7.455 + 0.008 * t ) * t ) * t );

/* Nutation */
   slaNutc ( date, &dpsi, &deps, &eps0 );

/* Equation of the equinoxes */
   return dpsi * cos ( eps0 ) + AS2R * ( 0.00264 * sin ( om ) +
                                         0.000063 * sin ( om + om ) );
}

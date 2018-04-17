#include "slalib.h"
#include "slamac.h"
void slaPrebn ( double bep0, double bep1, double rmatp[3][3] )
/*
**  - - - - - - - - -
**   s l a P r e b n
**  - - - - - - - - -
**
**  Generate the matrix of precession between two epochs,
**  using the old, pre-IAU1976, Bessel-Newcomb model, using
**  Kinoshita's formulation (double precision)
**
**  Given:
**     BEP0    double        beginning Besselian epoch
**     BEP1    double        ending Besselian epoch
**
**  Returned:
**     RMATP   double[3][3]  precession matrix
**
**  The matrix is in the sense   v(bep1)  =  rmatp * v(bep0)
**
**  Reference:
**     Kinoshita, H. (1975) 'Formulas for precession', SAO Special
**     Report No. 364, Smithsonian Institution Astrophysical
**     Observatory, Cambridge, Massachusetts.
**
**  Called:  slaDeuler
**
**  Defined in slamac.h:  DAS2R
**
**  Last revision:   30 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double bigt, t, tas2r, w, zeta, z, theta;

/* Interval between basic epoch B1850.0 and beginning epoch in TC */
   bigt  = ( bep0 - 1850.0 ) / 100.0;

/* Interval over which precession required, in tropical centuries */
   t = ( bep1 - bep0 ) / 100.0;

/* Euler angles */
   tas2r = t * DAS2R;
   w = 2303.5548 + ( 1.39720 + 0.000059 * bigt ) * bigt;
   zeta = (w + ( 0.30242 - 0.000269 * bigt + 0.017996 * t ) * t ) * tas2r;
   z = (w + ( 1.09478 + 0.000387 * bigt + 0.018324 * t ) * t ) * tas2r;
   theta = ( 2005.1125 + ( - 0.85294 - 0.000365* bigt ) * bigt +
           ( - 0.42647 - 0.000365 * bigt - 0.041802 * t ) * t ) * tas2r;

/* Rotation matrix */
   slaDeuler ( "ZYZ", -zeta, theta, -z, rmatp );
}

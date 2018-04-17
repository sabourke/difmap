#include "slalib.h"
#include "slamac.h"
void slaDvn ( double v[3], double uv[3], double *vm )
/*
**  - - - - - - -
**   s l a D v n
**  - - - - - - -
**
**  Normalizes a 3-vector also giving the modulus.
**
**  (double precision)
**
**  Given:
**     v       double[3]      vector
**
**  Returned:
**     uv      double[3]      unit vector in direction of v
**     *vm     double         modulus of v
**
**
**  If the modulus of v is zero, uv is set to zero as well.
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int i;
   double w1, w2;

/* Modulus */
   w1 = 0.0;
   for ( i = 0; i < 3; i++ ) {
      w2 = v[i];
      w1 += w2 * w2;
   }
   w1 = sqrt ( w1 );
   *vm = w1;

/* Normalize the vector */
   w1 = ( w1 > 0.0 ) ? w1 : 1.0;

   for ( i = 0; i < 3; i++ ) {
      uv[i] = v[i] / w1;
   }
}

#include "slalib.h"
#include "slamac.h"
void slaDmxv ( double dm[3][3], double va[3], double vb[3] )
/*
**  - - - - - - - -
**   s l a D m x v
**  - - - - - - - -
**
**  Performs the 3-d forward unitary transformation:
**     vector vb = matrix dm * vector va
**
**  (double precision)
**
**  Given:
**     dm       double[3][3]    matrix
**     va       double[3]       vector
**
**  Returned:
**     vb       double[3]       result vector
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int i, j;
   double w, vw[3];

/* Matrix dm * vector va -> vector vw */
   for ( j = 0; j < 3; j++ ) {
      w = 0.0;
      for ( i = 0; i < 3; i++ ) {
         w += dm[j][i] * va[i];
      }
      vw[j] = w;
   }

/* Vector vw -> vector vb */
   for ( j = 0; j < 3; j++ ) {
      vb[j] = vw[j];
   }
}

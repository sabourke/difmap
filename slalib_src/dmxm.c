#include "slalib.h"
#include "slamac.h"
void slaDmxm ( double a[3][3], double b[3][3], double c[3][3] )
/*
**  - - - - - - - -
**   s l a D m x m
**  - - - - - - - -
**
**  Product of two 3x3 matrices:
**     matrix c  =  matrix a  x  matrix b
**
**  (double precision)
**
**  Given:
**      a      double[3][3]        matrix
**      b      double[3][3]        matrix
**
**  Returned:
**      c      double[3][3]        matrix result
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int i, j, k;
   double w, wm[3][3];

/* Multiply into scratch matrix */
   for ( i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j++ ) {
         w = 0.0;
         for ( k = 0; k < 3; k++ ) {
            w += a[i][k] * b[k][j];
         }
         wm[i][j] = w;
      }
   }

/* Return the result */
   for ( j = 0; j < 3; j++ ) {
      for ( i = 0; i < 3; i++ ) {
         c[i][j] = wm[i][j];
      }
   }
}

#include "slalib.h"
#include "slamac.h"
#include <string.h>
void slaDeuler ( char *order, double phi, double theta,
                 double psi, double rmat[3][3] )
/*
**  - - - - - - - - - -
**   s l a D e u l e r
**  - - - - - - - - - -
**
**  Form a rotation matrix from the Euler angles - three successive
**  rotations about specified Cartesian axes.
**
**  (double precision)
**
**  Given:
**    *order char     specifies about which axes the rotations occur
**    phi    double   1st rotation (radians)
**    theta  double   2nd rotation (   "   )
**    psi    double   3rd rotation (   "   )
**
**  Returned:
**    rmat   double[3][3]  rotation matrix
**
**  A rotation is positive when the reference frame rotates
**  anticlockwise as seen looking towards the origin from the
**  positive region of the specified axis.
**
**  The characters of order define which axes the three successive
**  rotations are about.  A typical value is 'zxz', indicating that
**  rmat is to become the direction cosine matrix corresponding to
**  rotations of the reference frame through phi radians about the
**  old z-axis, followed by theta radians about the resulting x-axis,
**  then psi radians about the resulting z-axis.
**
**  The axis names can be any of the following, in any order or
**  combination:  x, y, z, uppercase or lowercase, 1, 2, 3.  Normal
**  axis labelling/numbering conventions apply;  the xyz (=123)
**  triad is right-handed.  Thus, the 'zxz' example given above
**  could be written 'zxz' or '313' (or even 'zxz' or '3xz').  Order
**  is terminated by length or by the first unrecognized character.
**
**  Fewer than three rotations are acceptable, in which case the later
**  angle arguments are ignored.  Zero rotations produces a unit rmat.
**
**  Last revision:   23 November 1995
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int j, i, l, n, k;
   double result[3][3], rotn[3][3], angle, s, c , w, wm[3][3];
   char axis;

/* Initialize result matrix */
   for ( j = 0; j < 3; j++ ) {
      for ( i = 0; i < 3; i++ ) {
         result[i][j] = ( i == j ) ? 1.0 : 0.0;
      }
   }

/* Establish length of axis string */
   l = strlen ( order );

/* Look at each character of axis string until finished */
   for ( n = 0; n < 3; n++ ) {
      if ( n <= l ) {

      /* Initialize rotation matrix for the current rotation */
         for ( j = 0; j < 3; j++ ) {
            for ( i = 0; i < 3; i++ ) {
               rotn[i][j] = ( i == j ) ? 1.0 : 0.0;
            }
         }

      /* Pick up the appropriate Euler angle and take sine & cosine */
         switch ( n ) {
         case 0 :
           angle = phi;
           break;
         case 1 :
           angle = theta;
           break;
         default:
           angle = psi;
           break;
         }
         s = sin ( angle );
         c = cos ( angle );

      /* Identify the axis */
         axis =  order[n];
         if ( ( axis == 'X' ) || ( axis == 'x' ) || ( axis == '1' ) ) {

         /* Matrix for x-rotation */
            rotn[1][1] = c;
            rotn[1][2] = s;
            rotn[2][1] = -s;
            rotn[2][2] = c;
         }
         else if ( ( axis == 'Y' ) || ( axis == 'y' ) || ( axis == '2' ) ) {

         /* Matrix for y-rotation */
            rotn[0][0] = c;
            rotn[0][2] = -s;
            rotn[2][0] = s;
            rotn[2][2] = c;
         }
         else if ( ( axis == 'Z' ) || ( axis == 'z' ) || ( axis == '3' ) ) {

         /* Matrix for z-rotation */
            rotn[0][0] = c;
            rotn[0][1] = s;
            rotn[1][0] = -s;
            rotn[1][1] = c;
         } else {

         /* Unrecognized character - fake end of string */
            l = 0;
         }

      /* Apply the current rotation (matrix rotn x matrix result) */
         for ( i = 0; i < 3; i++ ) {
            for ( j = 0; j < 3; j++ ) {
               w = 0.0;
               for ( k = 0; k < 3; k++ ) {
                  w += rotn[i][k] * result[k][j];
               }
               wm[i][j] = w;
            }
         }
         for ( j = 0; j < 3; j++ ) {
            for ( i= 0; i < 3; i++ ) {
               result[i][j] = wm[i][j];
            }
         }
      }
   }

/* Copy the result */
   for ( j = 0; j < 3; j++ ) {
      for ( i = 0; i < 3; i++ ) {
         rmat[i][j] = result[i][j];
      }
   }
}

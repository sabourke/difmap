#include "slalib.h"
#include "slamac.h"
double slaDvdv ( double va[3], double vb[3] )
/*
**  - - - - - - - -
**   s l a D v d v
**  - - - - - - - -
**
**  Scalar product of two 3-vectors.
**
**  (double precision)
**
**
**  Given:
**      va      double(3)     first vector
**      vb      double(3)     second vector
**
**
**  The result is the scalar product va.vb (double precision)
**
**
**  Last revision:   31 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   return va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2];
}

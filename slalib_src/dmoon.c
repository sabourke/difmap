#include "slalib.h"
#include "slamac.h"
void slaDmoon ( double date, double pv[6] )
/*
**  - - - - - - - - -
**   s l a D m o o n
**  - - - - - - - - -
**
**  Approximate geocentric position and velocity of the Moon
**  (double precision).
**
**  Given:
**     date     double      TDB (loosely ET) as a Modified Julian Date
**                                                  (JD-2400000.5)
**
**  Returned:
**     pv       double[6]   Moon x,y,z,xdot,ydot,zdot, mean equator
**                                   and equinox of date (AU, AU/s)
**
**  Notes:
**
**  1  This routine is a full implementation of the algorithm
**     published by Meeus (see reference).
**
**  2  Meeus quotes accuracies of 10 arcsec in longitude, 3 arcsec in
**     latitude and 0.2 arcsec in HP (equivalent to about 20 km in
**     distance).  Comparison with JPL DE200 over the interval
**     1960-2025 gives RMS errors of 3.7 arcsec and 83 mas/hour in
**     longitude, 2.3 arcsec and 48 mas/hour in latitude, 11 km
**     and 81 mm/s in distance.
**
**  3  The original algorithm is expressed in terms of the obsolete
**     timescale Ephemeris Time.  Either TDB or TT can be used, but
**     not UT without incurring significant errors (30 arcsec at
**     the present time) due to the Moon's 0.5 arcsec/sec movement.
**
**  4  The algorithm is based on pre IAU 1976 standards.  However,
**     the result has been moved onto the new (FK5) equinox, an
**     adjustment which is in any case much smaller than the
**     intrinsic accuracy of the procedure.
**
**  5  Velocity is obtained by a complete analytical differentiation
**     of the Meeus model.
**
**  Reference:
**     Meeus, l'Astronomie, June 1984, p348.
**
**  Defined in slamac.h:  DD2R, DAS2R, DS2R, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/

#define CJ 3155760000.0               /* Seconds per Julian century    */
                                      /*   ( = 86400 * 36525 )         */

#define ERADAU 4.2635212653763e-5     /* Earth equatorial radius in AU */
                                      /*   ( = 6378.137 / 149597870 )  */

#define B1950 1949.9997904423         /* Julian epoch of B1950         */

{
   double t, theta, sinom, cosom, domcom, wa, dwa, wb, dwb, wom,
          dwom, sinwom, coswom, v, dv, coeff, emn, empn, dn, fn, en,
          den, dtheta, ftheta, el, del, b, db, bf, dbf, p, dp, sp, r,
          dr, x, y, z, xd, yd, zd, sel, cel, sb, cb, rcb, rbd, w, epj,
          eqcor, eps, sineps, coseps, es, ec;
   int n, i;

/*
**  Coefficients for fundamental arguments
**
**   at J1900:  T^0, T^1, T^2, T^3
**   at epoch:  T^0, T^1
**
**  Units are degrees for position and Julian centuries for time.
*/

/* Moon's mean longitude */
   static double elp0 = 270.434164,
                 elp1 = 481267.8831,
                 elp2 = -0.001133,
                 elp3 = 0.0000019;
   double elp, delp;

/* Sun's mean anomaly */
   static double em0 = 358.475833,
                 em1 = 35999.0498,
                 em2 = -0.000150,
                 em3 = -0.0000033;
   double em, dem;

/* Moon's mean anomaly */
   static double emp0 = 296.104608,
                 emp1 = 477198.8491,
                 emp2 = 0.009192,
                 emp3 = 0.0000144;
   double emp, demp;

/* Moon's mean elongation */
   static double d0 = 350.737486,
                 d1 = 445267.1142,
                 d2 = -0.001436,
                 d3 = 0.0000019;
   double d, dd;

/* Mean distance of the Moon from its ascending node */
   static double f0 = 11.250889,
                 f1 = 483202.0251,
                 f2 = -0.003211,
                 f3 = -0.0000003;
   double f, df;

/* Longitude of the Moon's ascending node */
   static double om0 = 259.183275,
                 om1 = -1934.1420,
                 om2 = 0.002078,
                 om3 = 0.0000022;
   double om, dom;

/* Coefficients for (dimensionless) E factor */
   static double e1 = -0.002495,
                 e2 = -0.00000752;
   double e, de, esq, desq;

/* Coefficients for periodic variations etc */
   static double pac = 0.000233, pa0 = 51.2,
                                 pa1 = 20.2;
   static double pbc = -0.001778;
   static double pcc = 0.000817;
   static double pdc = 0.002011;
   static double pec = 0.003964, pe0 = 346.560,
                                 pe1 = 132.870,
                                 pe2 = -0.0091731;
   static double pfc = 0.001964;
   static double pgc = 0.002541;
   static double phc = 0.001964;
   static double pic = -0.024691;
   static double pjc = -0.004328, pj0 = 275.05,
                                  pj1 = -2.30;
   static double cw1 = 0.0004664;
   static double cw2 = 0.0000754;

/*
**  Coefficients for Moon longitude, latitude, parallax series
*/
   struct term {
      double coef;     /* coefficient of L, B or P term (deg) */
      int nem;         /* multiple of M  in argument          */
      int nemp;        /*     "    "  M'  "    "              */
      int nd;          /*     "    "  D   "    "              */
      int nf;          /*     "    "  F   "    "              */
      int ne;          /* power of e to multiply term by      */
   };

/*
** Longitude                      coeff      M    M'   D    F    n
*/
   static struct term tl[] = { {  6.288750,    0,   1,   0,   0,   0 },
                               {  1.274018,    0,  -1,   2,   0,   0 },
                               {  0.658309,    0,   0,   2,   0,   0 },
                               {  0.213616,    0,   2,   0,   0,   0 },
                               { -0.185596,    1,   0,   0,   0,   1 },
                               { -0.114336,    0,   0,   0,   2,   0 },
                               {  0.058793,    0,  -2,   2,   0,   0 },
                               {  0.057212,   -1,  -1,   2,   0,   1 },
                               {  0.053320,    0,   1,   2,   0,   0 },
                               {  0.045874,   -1,   0,   2,   0,   1 },
                               {  0.041024,   -1,   1,   0,   0,   1 },
                               { -0.034718,    0,   0,   1,   0,   0 },
                               { -0.030465,    1,   1,   0,   0,   1 },
                               {  0.015326,    0,   0,   2,  -2,   0 },
                               { -0.012528,    0,   1,   0,   2,   0 },
                               { -0.010980,    0,  -1,   0,   2,   0 },
                               {  0.010674,    0,  -1,   4,   0,   0 },
                               {  0.010034,    0,   3,   0,   0,   0 },
                               {  0.008548,    0,  -2,   4,   0,   0 },
                               { -0.007910,    1,  -1,   2,   0,   1 },
                               { -0.006783,    1,   0,   2,   0,   1 },
                               {  0.005162,    0,   1,  -1,   0,   0 },
                               {  0.005000,    1,   0,   1,   0,   1 },
                               {  0.004049,   -1,   1,   2,   0,   1 },
                               {  0.003996,    0,   2,   2,   0,   0 },
                               {  0.003862,    0,   0,   4,   0,   0 },
                               {  0.003665,    0,  -3,   2,   0,   0 },
                               {  0.002695,   -1,   2,   0,   0,   1 },
                               {  0.002602,    0,   1,  -2,  -2,   0 },
                               {  0.002396,   -1,  -2,   2,   0,   1 },
                               { -0.002349,    0,   1,   1,   0,   0 },
                               {  0.002249,   -2,   0,   2,   0,   2 },
                               { -0.002125,    1,   2,   0,   0,   1 },
                               { -0.002079,    2,   0,   0,   0,   2 },
                               {  0.002059,   -2,  -1,   2,   0,   2 },
                               { -0.001773,    0,   1,   2,  -2,   0 },
                               { -0.001595,    0,   0,   2,   2,   0 },
                               {  0.001220,   -1,  -1,   4,   0,   1 },
                               { -0.001110,    0,   2,   0,   2,   0 },
                               {  0.000892,    0,   1,  -3,   0,   0 },
                               { -0.000811,    1,   1,   2,   0,   1 },
                               {  0.000761,   -1,  -2,   4,   0,   1 },
                               {  0.000717,   -2,   1,   0,   0,   2 },
                               {  0.000704,   -2,   1,  -2,   0,   2 },
                               {  0.000693,    1,  -2,   2,   0,   1 },
                               {  0.000598,   -1,   0,   2,  -2,   1 },
                               {  0.000550,    0,   1,   4,   0,   0 },
                               {  0.000538,    0,   4,   0,   0,   0 },
                               {  0.000521,   -1,   0,   4,   0,   1 },
                               {  0.000486,    0,   2,  -1,   0,   0 } };
   static int NL = ( sizeof tl / sizeof ( struct term ) );

/*
** Latitude                       coeff      M    M'   D    F    n
*/
   static struct term tb[] = { {  5.128189,    0,   0,   0,   1,   0 },
                               {  0.280606,    0,   1,   0,   1,   0 },
                               {  0.277693,    0,   1,   0,  -1,   0 },
                               {  0.173238,    0,   0,   2,  -1,   0 },
                               {  0.055413,    0,  -1,   2,   1,   0 },
                               {  0.046272,    0,  -1,   2,  -1,   0 },
                               {  0.032573,    0,   0,   2,   1,   0 },
                               {  0.017198,    0,   2,   0,   1,   0 },
                               {  0.009267,    0,   1,   2,  -1,   0 },
                               {  0.008823,    0,   2,   0,  -1,   0 },
                               {  0.008247,   -1,   0,   2,  -1,   1 },
                               {  0.004323,    0,  -2,   2,  -1,   0 },
                               {  0.004200,    0,   1,   2,   1,   0 },
                               {  0.003372,   -1,   0,  -2,   1,   1 },
                               {  0.002472,   -1,  -1,   2,   1,   1 },
                               {  0.002222,   -1,   0,   2,   1,   1 },
                               {  0.002072,   -1,  -1,   2,  -1,   1 },
                               {  0.001877,   -1,   1,   0,   1,   1 },
                               {  0.001828,    0,  -1,   4,  -1,   0 },
                               { -0.001803,    1,   0,   0,   1,   1 },
                               { -0.001750,    0,   0,   0,   3,   0 },
                               {  0.001570,   -1,   1,   0,  -1,   1 },
                               { -0.001487,    0,   0,   1,   1,   0 },
                               { -0.001481,    1,   1,   0,   1,   1 },
                               {  0.001417,   -1,  -1,   0,   1,   1 },
                               {  0.001350,   -1,   0,   0,   1,   1 },
                               {  0.001330,    0,   0,  -1,   1,   0 },
                               {  0.001106,    0,   3,   0,   1,   0 },
                               {  0.001020,    0,   0,   4,  -1,   0 },
                               {  0.000833,    0,  -1,   4,   1,   0 },
                               {  0.000781,    0,   1,   0,  -3,   0 },
                               {  0.000670,    0,  -2,   4,   1,   0 },
                               {  0.000606,    0,   0,   2,  -3,   0 },
                               {  0.000597,    0,   2,   2,  -1,   0 },
                               {  0.000492,   -1,   1,   2,  -1,   1 },
                               {  0.000450,    0,   2,  -2,  -1,   0 },
                               {  0.000439,    0,   3,   0,  -1,   0 },
                               {  0.000423,    0,   2,   2,   1,   0 },
                               {  0.000422,    0,  -3,   2,  -1,   0 },
                               { -0.000367,    1,  -1,   2,   1,   1 },
                               { -0.000353,    1,   0,   2,   1,   1 },
                               {  0.000331,    0,   0,   4,   1,   0 },
                               {  0.000317,   -1,   1,   2,   1,   1 },
                               {  0.000306,   -2,   0,   2,  -1,   2 },
                               { -0.000283,    0,   1,   0,   3,   0 } };
   static int NB = ( sizeof tb / sizeof ( struct term ) );

/*
** Parallax                       coeff      M    M'   D    F    n
*/
   static struct term tp[] = { {  0.950724,    0,   0,   0,   0,   0 },
                               {  0.051818,    0,   1,   0,   0,   0 },
                               {  0.009531,    0,  -1,   2,   0,   0 },
                               {  0.007843,    0,   0,   2,   0,   0 },
                               {  0.002824,    0,   2,   0,   0,   0 },
                               {  0.000857,    0,   1,   2,   0,   0 },
                               {  0.000533,   -1,   0,   2,   0,   1 },
                               {  0.000401,   -1,  -1,   2,   0,   1 },
                               {  0.000320,   -1,   1,   0,   0,   1 },
                               { -0.000271,    0,   0,   1,   0,   0 },
                               { -0.000264,    1,   1,   0,   0,   1 },
                               { -0.000198,    0,  -1,   0,   2,   0 },
                               {  0.000173,    0,   3,   0,   0,   0 },
                               {  0.000167,    0,  -1,   4,   0,   0 },
                               { -0.000111,    1,   0,   0,   0,   1 },
                               {  0.000103,    0,  -2,   4,   0,   0 },
                               { -0.000084,    0,   2,  -2,   0,   0 },
                               { -0.000083,    1,   0,   2,   0,   1 },
                               {  0.000079,    0,   2,   2,   0,   0 },
                               {  0.000072,    0,   0,   4,   0,   0 },
                               {  0.000064,   -1,   1,   2,   0,   1 },
                               { -0.000063,    1,  -1,   2,   0,   1 },
                               {  0.000041,    1,   0,   1,   0,   1 },
                               {  0.000035,   -1,   2,   0,   0,   1 },
                               { -0.000033,    0,   3,  -2,   0,   0 },
                               { -0.000030,    0,   1,   1,   0,   0 },
                               { -0.000029,    0,   0,  -2,   2,   0 },
                               { -0.000029,    1,   2,   0,   0,   1 },
                               {  0.000026,   -2,   0,   2,   0,   2 },
                               { -0.000023,    0,   1,  -2,   2,   0 },
                               {  0.000019,   -1,  -1,   4,   0,   1 } };
   static int NP = ( sizeof tp / sizeof ( struct term ) );



/* --------------------- */
/* Execution starts here */
/* --------------------- */

/* Centuries since J1900 */
   t = ( date - 15019.5 ) / 36525.0;

/* --------------------- */
/* Fundamental arguments */
/* --------------------- */

/* Arguments (radians) and derivatives (radians per Julian century)
   for the current epoch */

/* Moon's mean longitude */
   elp = DD2R * dmod ( elp0 + ( elp1 + ( elp2 + elp3 * t ) * t ) * t,
                                                                    360.0 );
   delp = DD2R * ( elp1 + ( 2.0 * elp2 + 3.0 *elp3 * t ) * t );

/* Sun's mean anomaly */
   em = DD2R * dmod ( em0 + ( em1 + ( em2 + em3 * t ) * t ) * t, 360.0 );
   dem = DD2R * ( em1 + ( 2.0 * em2 + 3.0 * em3 * t ) * t );

/* Moon's mean anomaly */
   emp = DD2R * dmod ( emp0 + ( emp1 + ( emp2 + emp3 * t ) * t ) * t,
                                                                    360.0 );
   demp = DD2R * ( emp1 + ( 2.0 * emp2 + 3.0 * emp3 * t ) * t );

/* Moon's mean elongation */
   d = DD2R * dmod ( d0 + ( d1 + ( d2 + d3 * t ) * t ) * t, 360.0 );
   dd = DD2R * ( d1 + ( 2.0 * d2 + 3.0 * d3 * t ) * t );

/* Mean distance of the Moon from its ascending node */
   f = DD2R * dmod ( f0 + ( f1 + ( f2 + f3 * t ) * t ) * t, 360.0 );
   df = DD2R * ( f1 + ( 2.0 * f2 + 3.0 * f3 * t ) * t );

/* Longitude of the Moon's ascending node */
   om = DD2R * dmod ( om0 + ( om1 + ( om2 + om3 * t ) * t ) * t, 360.0 );
   dom = DD2R * ( om1 + ( 2.0 * om2 + 3.0 * om3 * t ) * t );
   sinom = sin ( om );
   cosom = cos ( om );
   domcom = dom * cosom;

/* Add the periodic variations */
   theta = DD2R * ( pa0 + pa1 * t );
   wa = sin ( theta );
   dwa = DD2R * pa1 * cos ( theta );
   theta = DD2R * ( pe0 + ( pe1 + pe2 * t ) * t );
   wb = pec * sin ( theta );
   dwb = DD2R * pec*( pe1 + 2.0 * pe2 * t ) * cos ( theta );
   elp += DD2R * ( pac * wa + wb + pfc * sinom );
   delp += DD2R * ( pac * dwa + dwb + pfc * domcom );
   em += DD2R * pbc * wa;
   dem += DD2R * pbc * dwa;
   emp += DD2R * ( pcc * wa + wb + pgc * sinom );
   demp += DD2R * ( pcc * dwa + dwb + pgc * domcom );
   d += DD2R * ( pdc * wa + wb + phc * sinom );
   dd += DD2R * ( pdc * dwa + dwb + phc * domcom );
   wom = om + DD2R * ( pj0 + pj1 * t );
   dwom = dom + DD2R * pj1;
   sinwom = sin ( wom );
   coswom = cos ( wom );
   f += DD2R * ( wb + pic * sinom + pjc * sinwom );
   df += DD2R * ( dwb + pic * domcom + pjc * dwom * coswom );

/* E-factor, and square */
   e = 1.0 + ( e1 + e2 * t ) * t;
   de = e1 + 2.0 * e2 * t;
   esq = e * e;
   desq = 2.0 * e * de;

/* ----------------- */
/* Series expansions */
/* ----------------- */

/* Longitude */
   v = 0.0;
   dv = 0.0;
   for ( n = NL -1; n >= 0; n-- ) {
      coeff = tl[n].coef;
      emn = (double) tl[n].nem;
      empn = (double) tl[n].nemp;
      dn = (double) tl[n].nd;
      fn = (double) tl[n].nf;
      i = tl[n].ne;
      if ( i == 0 ) {
         en = 1.0;
         den = 0.0;
      } else if ( i == 1 ) {
         en = e;
         den = de;
      } else {
         en = esq;
         den = desq;
      }
      theta = emn * em + empn * emp + dn * d + fn * f;
      dtheta = emn * dem + empn * demp + dn * dd + fn * df;
      ftheta = sin ( theta );
      v += coeff * ftheta * en;
      dv += coeff * ( cos ( theta ) * dtheta * en + ftheta * den );
   }
   el = elp + DD2R * v;
   del = ( delp + DD2R * dv ) / CJ;

/* Latitude */
   v = 0.0;
   dv = 0.0;
   for ( n = NB - 1; n >= 0; n-- ) {
      coeff = tb[n].coef;
      emn = (double) tb[n].nem;
      empn = (double) tb[n].nemp;
      dn = (double) tb[n].nd;
      fn = (double) tb[n].nf;
      i = tb[n].ne;
      if ( i == 0 ) {
         en = 1.0;
         den = 0.0;
      } else if ( i == 1 ) {
         en = e;
         den = de;
      } else {
         en = esq;
         den = desq;
      }
      theta = emn * em + empn * emp + dn * d + fn * f;
      dtheta = emn * dem + empn * demp + dn * dd + fn * df;
      ftheta = sin ( theta );
      v += coeff * ftheta * en;
      dv += coeff * ( cos ( theta ) * dtheta * en + ftheta * den );
   }
   bf = 1.0 - cw1 * cosom - cw2 * coswom;
   dbf = cw1 * dom * sinom + cw2 * dwom * sinwom;
   b = DD2R * v * bf;
   db = DD2R * ( dv * bf + v * dbf ) / CJ;

/* Parallax */
   v = 0.0;
   dv = 0.0;
   for ( n = NP - 1; n >= 0; n-- ) {
      coeff = tp[n].coef;
      emn = (double) tp[n].nem;
      empn = (double) tp[n].nemp;
      dn = (double) tp[n].nd;
      fn = (double) tp[n].nf;
      i = tp[n].ne;
      if ( i == 0 ) {
         en = 1.0;
         den = 0.0;
      } else if ( i == 1 ) {
         en = e;
         den = de;
      } else {
         en = esq;
         den = desq;
      }
      theta = emn * em + empn * emp + dn * d + fn * f;
      dtheta = emn * dem + empn * demp + dn * dd + fn * df;
      ftheta = cos ( theta );
      v += coeff * ftheta * en;
      dv += coeff* ( - sin ( theta ) * dtheta * en + ftheta * den );
   }
   p = DD2R * v;
   dp = DD2R * dv / CJ;

/* ------------------------------ */
/* Transformation into final form */
/* ------------------------------ */

/* Parallax to distance (AU, AU/sec) */
   sp = sin ( p );
   r = ERADAU / sp;
   dr = - r * dp * cos ( p ) / sp;

/* Longitude, latitude to x, y, z (AU) */
   sel = sin ( el );
   cel = cos ( el );
   sb = sin ( b );
   cb = cos ( b );
   rcb = r * cb;
   rbd = r * db;
   w = rbd * sb - cb * dr;
   x = rcb * cel;
   y = rcb * sel;
   z = r * sb;
   xd = - y * del - w * cel;
   yd = x * del - w * sel;
   zd = rbd * cb + sb * dr;

/* Julian centuries since J2000 */
   t = ( date - 51544.5 ) / 36525.0;

/* Fricke equinox correction */
   epj = 2000.0 + t * 100.0;
   eqcor = DS2R * ( 0.035 + 0.00085 * ( epj - B1950 ) );

/* Mean obliquity (IAU 1976) */
   eps = DAS2R *
      ( 84381.448 + ( - 46.8150 + ( - 0.00059 + 0.001813 * t ) * t ) * t );

/* To the equatorial system, mean of date, FK5 system */
   sineps = sin ( eps );
   coseps = cos ( eps );
   es = eqcor * sineps;
   ec = eqcor * coseps;
   pv[0] = x - ec * y + es * z;
   pv[1] = eqcor * x + y * coseps - z * sineps;
   pv[2] = y * sineps + z * coseps;
   pv[3] = xd - ec * yd + es * zd;
   pv[4] = eqcor * xd + yd * coseps - zd * sineps;
   pv[5] = yd * sineps + zd * coseps;
}

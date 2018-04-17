#include "slalib.h"
#include "slamac.h"
void slaPlanet ( double date, int np, double pv[6], int *jstat )
/*
**  - - - - - - - - - -
**   s l a P l a n e t
**  - - - - - - - - - -
**
**  Approximate heliocentric position and velocity of a specified
**  major planet (Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus
**  or Neptune).
**
**  Given:
**     date     double      TDB (loosely ET) as a Modified Julian Date
**                                                  (JD-2400000.5)
**     np       int         planet (1=Mercury, 2=Venus, 3=EMB, ...
**                                                  ... 8=Neptune)
**
**  Returned:
**     pv       double[6]   heliocentric x,y,z,xdot,ydot,zdot, J2000
**                                           equatorial triad (AU,AU/s)
**
**     *jstat   int         status: -1 = illegal NP (outside 1-8)
**                                   0 = OK
**                                  +1 = warning: date outside 1000-3000
**                                  +2 = warning: solution didn't converge
**
**  Notes
**
**  1  The epoch, date, is in the TDB timescale and is a Modified
**     Julian Date (JD-2400000.5).
**
**  2  The reference frame is equatorial and is with respect to the
**     mean equinox and ecliptic of epoch J2000.
**
**  3  If an np value outside the range 1-8 is supplied, an error
**     status (jstat = -1) is returned and the pv vector set to zeroes.
**
**  4  The algorithm is due to J.L. Simon, P. Bretagnon, J. Chapront,
**     M. Chapront-Touze, G. Francou and J. Laskar (Bureau des
**     Longitudes, Paris, France).
**
**  5  Comparisons of the present routine with the JPL DE200 ephemeris
**     give the following RMS errors over the interval 1960-2025:
**
**                      position (km)     speed (metre/sec)
**
**        Mercury            334               0.437
**        Venus             1060               0.855
**        EMB               2010               0.815
**        Mars              7690               1.98
**        Jupiter          71700               7.70
**        Saturn          199000              19.4
**        Uranus          564000              16.4
**        Neptune         158000              14.4
**
**     From comparisons with DE102, Simon et al quote the following
**     longitude accuracies over the interval 1800-2200:
**
**        Mercury                 4"
**        Venus                   5"
**        EMB                     6"
**        Mars                   17"
**        Jupiter                71"
**        Saturn                 81"
**        Uranus                 86"
**        Neptune                11"
**
**     Over the interval 1000-3000, the accuracy is better than 1.5
**     times that over 1800-2200.  Outside the interval 1000-3000 the
**     accuracy declines.
**
**  6  The present SLALIB C implementation follows the original
**     Simon et al Fortran code closely, and delivers essentially
**     the same results.  The changes are these:
**
**       *  The date is supplied as a Modified Julian Date rather
**          than a Julian Date (MJD = JD - 2400000.5).
**
**       *  The result is returned only in equatorial Cartesian form;
**          the ecliptic longitude, latitude and radius vector are not
**          returned.
**
**       *  The result is in the J2000 equatorial frame, not ecliptic.
**
**       *  The velocity is in AU per second, not AU per day.
**
**       *  Everything is done in-line:  there are no calls to other
**          routines.
**
**       *  Different error/warning status values are used.
**
**       *  A different Kepler's-equation-solver is used (avoiding
**          use of COMPLEX*16).
**
**       *  Polynomials in T are nested to minimize rounding errors.
**
**       *  Explicit double-precision constants are used to avoid
**          mixed-mode expressions.
**
**  7  For np=3 the result is for the Earth-Moon Barycentre.  To
**     obtain the heliocentric position and velocity of the Earth,
**     either use the SLALIB routine slaEvp or use slaDmoon and
**     subtract 0.012150581 times the geocentric Moon vector from
**     the EMB vector produced by the present routine.  (The Moon
**     vector should be precessed to J2000 first, but this can
**     be omitted for modern epochs without introducing significant
**     inaccuracy.)
**
**  8  The status, jstat, indicates the most serious condition
**     encountered, where illegal np is considered the most serious,
**     followed by failure to converge, then remote epoch.
**
**  Reference:  Astron. Astrophys. 282, 663 (1994).
**
**  Defined in slamac.h:  D2PI, DAS2R, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/

#define KMAX    10                    /* Maximum number of iterations  */
                                      /* allowed when solving Kepler's */
                                      /* equation                      */

#define SINEPS  0.3977771559319137    /* Sin and cos of J2000 mean */
#define COSEPS  0.9174820620691818    /* obliquity (IAU 1976)      */

#define GKS     1.99098367477e-7      /* Gaussian gravitational constant */
                                      /* divided by seconds per day      */
                                      /* ( 0.017202098950 / 86400 )      */
{
   int i, ip, j, k;
   double t, da, dl, de, dp, di, dom, dmu, arga, argl, am, ae,
          dae, ae2, at, r, v, si2, xq, xp, tl, xsw, xcw, xm2,
          xf, ci2, xms, xmc, xpxq2, x, y, z;

/* Planetary inverse masses */
   static double amas[8] = {
      6023600.0,
       408523.5,
       328900.5,
      3098710.0,
       1047.355,
         3498.5,
        22869.0,
        19314.0
   };

/*
**    Tables giving the mean Keplerian elements, limited to T^2 terms:
**
**    a       semi-major axis (AU)
**    dlm     mean longitude (degree and arcsecond)
**    e       eccentricity
**    pi      longitude of the perihelion (degree and arcsecond)
**    dinc    inclination (degree and arcsecond)
**    omega   longitude of the ascending node (degree and arcsecond)
*/
   static double a[8][3] = {
      {  0.3870983098,           0.0,     0.0 },
      {  0.7233298200,           0.0,     0.0 },
      {  1.0000010178,           0.0,     0.0 },
      {  1.5236793419,         3e-10,     0.0 },
      {  5.2026032092,     19132e-10, -39e-10 },
      {  9.5549091915, -0.0000213896, 444e-10 },
      { 19.2184460618,     -3716e-10, 979e-10 },
      { 30.1103868694,    -16635e-10, 686e-10 }
   };
   static double dlm[8][3] = {
      { 252.25090552, 5381016286.88982,  -1.92789 },
      { 181.97980085, 2106641364.33548,   0.59381 },
      { 100.46645683, 1295977422.83429,  -2.04411 },
      { 355.43299958,  689050774.93988,   0.94264 },
      {  34.35151874,  109256603.77991, -30.60378 },
      {  50.07744430,   43996098.55732,  75.61614 },
      { 314.05500511,   15424811.93933,  -1.75083 },
      { 304.34866548,    7865503.20744,   0.21103 }
   };
   static double e[8][3] = {
      { 0.2056317526,  0.0002040653,      -28349e-10 },
      { 0.0067719164, -0.0004776521,       98127e-10 },
      { 0.0167086342, -0.0004203654,   -0.0000126734 },
      { 0.0934006477,  0.0009048438,      -80641e-10 },
      { 0.0484979255,  0.0016322542,   -0.0000471366 },
      { 0.0555481426, -0.0034664062,   -0.0000643639 },
      { 0.0463812221, -0.0002729293,    0.0000078913 },
      { 0.0094557470,  0.0000603263,             0.0 }
   };
   static double pi[8][3] = {
      {  77.45611904,  5719.11590,   -4.83016 },
      { 131.56370300,   175.48640, -498.48184 },
      { 102.93734808, 11612.35290,   53.27577 },
      { 336.06023395, 15980.45908,  -62.32800 },
      {  14.33120687,  7758.75163,  259.95938 },
      {  93.05723748, 20395.49439,  190.25952 },
      { 173.00529106,  3215.56238,  -34.09288 },
      {  48.12027554,  1050.71912,   27.39717 }
   };
   static double dinc[8][3] = {
      { 7.00498625, -214.25629,   0.28977 },
      { 3.39466189,  -30.84437, -11.67836 },
      {        0.0,  469.97289,  -3.35053 },
      { 1.84972648, -293.31722,  -8.11830 },
      { 1.30326698,  -71.55890,  11.95297 },
      { 2.48887878,   91.85195, -17.66225 },
      { 0.77319689,  -60.72723,   1.25759 },
      { 1.76995259,    8.12333,   0.08135 }
   };
   static double omega[8][3] = {
      {  48.33089304,  -4515.21727,  -31.79892 },
      {  76.67992019, -10008.48154,  -51.32614 },
      { 174.87317577,  -8679.27034,   15.34191 },
      {  49.55809321, -10620.90088, -230.57416 },
      { 100.46440702,   6362.03561,  326.52178 },
      { 113.66550252,  -9240.19942,  -66.23743 },
      {  74.00595701,   2669.15033,  145.93964 },
      { 131.78405702,   -221.94322,   -0.78728 }
   };
/*
**    Tables for trigonometric terms to be added to the mean elements
**    of the semi-major axes.
*/
   static double kp[8][9] = {
      { 69613.0, 75645.0, 88306.0, 59899.0, 15746.0, 71087.0,
                                                142173.0,  3086.0,    0.0 },
      { 21863.0, 32794.0, 26934.0, 10931.0, 26250.0, 43725.0,
                                                 53867.0, 28939.0,    0.0 },
      { 16002.0, 21863.0, 32004.0, 10931.0, 14529.0, 16368.0,
                                                 15318.0, 32794.0,    0.0 },
      {  6345.0,  7818.0, 15636.0,  7077.0,  8184.0, 14163.0,
                                                  1107.0,  4872.0,    0.0 },
      {  1760.0,  1454.0,  1167.0,   880.0,   287.0,  2640.0,
                                                    19.0,  2047.0, 1454.0 },
      {   574.0,     0.0,   880.0,   287.0,    19.0,  1760.0,
                                                  1167.0,   306.0,  574.0 },
      {   204.0,     0.0,   177.0,  1265.0,     4.0,   385.0,
                                                   200.0,   208.0,  204.0 },
      {     0.0,   102.0,   106.0,     4.0,    98.0,  1367.0,
                                                   487.0,   204.0,    0.0 }
   };
   static double ca[8][9] = {
    {       4.0,    -13.0,    11.0,    -9.0,    -9.0,    -3.0,
                                                    -1.0,     4.0,    0.0 },
    {    -156.0,     59.0,   -42.0,     6.0,    19.0,   -20.0,
                                                   -10.0,   -12.0,    0.0 },
    {      64.0,   -152.0,    62.0,    -8.0,    32.0,   -41.0,
                                                    19.0,   -11.0,    0.0 },
    {     124.0,    621.0,  -145.0,   208.0,    54.0,   -57.0,
                                                    30.0,    15.0,    0.0 },
    {  -23437.0,  -2634.0,  6601.0,  6259.0, -1507.0, -1821.0,
                                                  2620.0, -2115.0,-1489.0 },
    {   62911.0,-119919.0, 79336.0, 17814.0,-24241.0, 12068.0,
                                                  8306.0, -4893.0, 8902.0 },
    {  389061.0,-262125.0,-44088.0,  8387.0,-22976.0, -2093.0,
                                                  -615.0, -9720.0, 6633.0 },
    { -412235.0,-157046.0,-31430.0, 37817.0, -9740.0,   -13.0,
                                                 -7449.0,  9644.0,    0.0 }
   };
   static double sa[8][9] = {
      {     -29.0,    -1.0,     9.0,     6.0,    -6.0,     5.0,
                                                     4.0,     0.0,    0.0 },
      {     -48.0,  -125.0,   -26.0,   -37.0,    18.0,   -13.0,
                                                   -20.0,    -2.0,    0.0 },
      {    -150.0,   -46.0,    68.0,    54.0,    14.0,    24.0,
                                                   -28.0,    22.0,    0.0 },
      {    -621.0,   532.0,  -694.0,   -20.0,   192.0,   -94.0,
                                                    71.0,   -73.0,    0.0 },
      {  -14614.0,-19828.0, -5869.0,  1881.0, -4372.0, -2255.0,
                                                   782.0,   930.0,  913.0 },
      {  139737.0,     0.0, 24667.0, 51123.0, -5102.0,  7429.0,
                                                 -4095.0, -1976.0,-9566.0 },
      { -138081.0,     0.0, 37205.0,-49039.0,-41901.0,-33872.0,
                                                -27037.0,-12474.0,18797.0 },
      {       0.0, 28492.0,133236.0, 69654.0, 52322.0,-49577.0,
                                                -26430.0, -3593.0,    0.0 }
   };
/*
**    Tables giving the trigonometric terms to be added to the mean
**    elements of the mean longitudes.
*/
   static double kq[8][10] = {
      {  3086.0, 15746.0, 69613.0, 59899.0, 75645.0,
                                      88306.0, 12661.0, 2658.0,  0.0,   0.0 },
      { 21863.0, 32794.0, 10931.0,    73.0,  4387.0,
                                      26934.0,  1473.0, 2157.0,  0.0,   0.0 },
      {    10.0, 16002.0, 21863.0, 10931.0,  1473.0,
                                      32004.0,  4387.0,   73.0,  0.0,   0.0 },
      {    10.0,  6345.0,  7818.0,  1107.0, 15636.0,
                                       7077.0,  8184.0,  532.0, 10.0,   0.0 },
      {    19.0,  1760.0,  1454.0,   287.0,  1167.0,
                                        880.0,   574.0, 2640.0, 19.0,1454.0 },
      {    19.0,   574.0,   287.0,   306.0,  1760.0,
                                         12.0,    31.0,   38.0, 19.0, 574.0 },
      {     4.0,   204.0,   177.0,     8.0,    31.0,
                                        200.0,  1265.0,  102.0,  4.0, 204.0 },
      {     4.0,   102.0,   106.0,     8.0,    98.0,
                                       1367.0,   487.0,  204.0,  4.0, 102.0 }
   };
   static double cl[8][10] = {
    {      21.0,    -95.0,  -157.0,    41.0,    -5.0,
                                      42.0,   23.0,   30.0,     0.0,    0.0 },
    {    -160.0,   -313.0,  -235.0,    60.0,   -74.0,
                                     -76.0,  -27.0,   34.0,     0.0,    0.0 },
    {    -325.0,   -322.0,   -79.0,   232.0,   -52.0,
                                      97.0,   55.0,  -41.0,     0.0,    0.0 },
    {    2268.0,   -979.0,   802.0,   602.0,  -668.0,
                                     -33.0,  345.0,  201.0,   -55.0,    0.0 },
    {    7610.0,  -4997.0, -7689.0, -5841.0, -2617.0,
                                    1115.0, -748.0, -607.0,  6074.0,  354.0 },
    {  -18549.0,  30125.0, 20012.0,  -730.0,   824.0,
                                      23.0, 1289.0, -352.0,-14767.0,-2062.0 },
    { -135245.0, -14594.0,  4197.0, -4030.0, -5630.0,
                                   -2898.0, 2540.0, -306.0,  2939.0, 1986.0 },
    {   89948.0,   2103.0,  8963.0,  2695.0,  3682.0,
                                    1648.0,  866.0, -154.0, -1963.0, -283.0 }
   };
   static double sl[8][10] = {
    {   -342.0,    136.0,   -23.0,    62.0,    66.0,
                                 -52.0,   -33.0,    17.0,     0.0,     0.0 },
    {    524.0,   -149.0,   -35.0,   117.0,   151.0,
                                 122.0,   -71.0,   -62.0,     0.0,     0.0 },
    {   -105.0,   -137.0,   258.0,    35.0,  -116.0,
                                 -88.0,  -112.0,   -80.0,     0.0,     0.0 },
    {    854.0,   -205.0,  -936.0,  -240.0,   140.0,
                                -341.0,   -97.0,  -232.0,   536.0,     0.0 },
    { -56980.0,   8016.0,  1012.0,  1448.0, -3024.0,
                               -3710.0,   318.0,   503.0,  3767.0,   577.0 },
    { 138606.0, -13478.0, -4964.0,  1441.0, -1319.0,
                               -1482.0,   427.0,  1236.0, -9167.0, -1918.0 },
    {  71234.0, -41116.0,  5334.0, -4935.0, -1848.0,
                                  66.0,   434.0, -1748.0,  3780.0,  -701.0 },
    { -47645.0,  11647.0,  2166.0,  3194.0,   679.0,
                                   0.0,  -244.0,  -419.0, -2531.0,    48.0 }
   };


/* Validate the planet number */
   if ( np < 1 || np > 8 ) {
      *jstat = -1;
      for ( i = 0; i <= 5; i++ ) pv[i] = 0.0;
      return;
   } else {
      ip = np - 1;
   }

/* Time: Julian millennia since J2000 */
   t = ( date - 51544.5 ) / 365250.0;

/* OK status unless remote epoch */
   *jstat = ( fabs ( t ) <= 1.0 ) ? 0 : 1;

/* Compute the mean elements */
   da = a[ip][0] + ( a[ip][1] + a[ip][2] * t ) * t;
   dl = ( 3600.0 * dlm[ip][0] + ( dlm[ip][1] + dlm[ip][2] * t ) * t )
                                                                  * DAS2R;
   de = e[ip][0] + ( e[ip][1] + e[ip][2] * t ) * t;
   dp = dmod ( ( 3600.0 * pi[ip][0] + ( pi[ip][1] + pi[ip][2] * t ) * t )
                                                           * DAS2R,D2PI );
   di = ( 3600.0 * dinc[ip][0] + ( dinc[ip][1] + dinc[ip][2] * t ) * t )
                                                                  * DAS2R;
   dom = dmod( ( 3600.0 * omega[ip][0] + ( omega[ip][1]
                               + omega[ip][2] * t ) * t ) * DAS2R, D2PI );

/* Apply the trigonometric terms */
   dmu = 0.35953620 * t;
   for ( j = 0; j <= 7; j++ ) {
      arga = kp[ip][j] * dmu;
      argl = kq[ip][j] * dmu;
      da += ( ca[ip][j] * cos ( arga ) + sa[ip][j] * sin ( arga ) ) * 1e-7;
      dl += ( cl[ip][j] * cos ( argl ) + sl[ip][j] * sin ( argl ) ) * 1e-7;
   }
   arga = kp[ip][8] * dmu;
   da += t * ( ca[ip][8] * cos ( arga ) + sa[ip][8] * sin ( arga ) ) * 1e-7;
   for ( j = 8; j <= 9; j++ ) {
      argl = kq[ip][j] * dmu;
      dl += t * ( cl[ip][j] * cos ( argl ) + sl[ip][j] * sin ( argl ) )
                                                                     * 1e-7;
   }
   dl = dmod ( dl, D2PI );

/* Iterative solution of Kepler's equation to get eccentric anomaly */
   am = dl - dp;
   ae = am + de * sin ( am );
   k = 0;
   do {
      dae = ( am - ae + de * sin ( ae ) ) / ( 1.0 - de * cos ( ae ) );
      ae += dae;
      if ( k++ >= KMAX ) *jstat = 2;
      }
   while ( k < KMAX && fabs ( dae ) > 1e-12 );

/* True anomaly */
   ae2 = ae / 2.0;
   at = 2.0 * atan2 ( sqrt ( ( 1.0 + de ) / ( 1.0 - de ) )
                                           * sin ( ae2 ), cos ( ae2 ) );

/* Distance (AU) and speed (radians per second) */
   r = da * ( 1.0 - de * cos ( ae ) );
   v = GKS * sqrt ( ( 1.0 + 1.0 / amas[ip] ) / ( da * da * da ) );

   si2 = sin ( di / 2.0 );
   xq = si2 * cos ( dom );
   xp = si2 * sin ( dom );
   tl = at + dp;
   xsw = sin ( tl );
   xcw = cos ( tl );
   xm2 = 2.0 * ( xp * xcw - xq * xsw );
   xf = da / sqrt ( 1.0 - de * de );
   ci2 = cos ( di / 2 );
   xms = ( de * sin ( dp ) + xsw ) * xf;
   xmc = ( de * cos ( dp ) + xcw ) * xf;
   xpxq2 = 2.0 * xp * xq;

/* Position ( J2000 ecliptic x,y,z in AU ) */
   x = r * ( xcw - xm2 * xp );
   y = r * ( xsw + xm2 * xq );
   z = r * ( - xm2 * ci2 );

/*  Rotate to equatorial */
   pv[0] = x;
   pv[1] = y * COSEPS - z * SINEPS;
   pv[2] = y * SINEPS + z * COSEPS;

/*  Velocity ( J2000 ecliptic xdot,ydot,zdot in AU/s) */
   x = v * ( (  - 1.0 + 2.0 * xp * xp ) * xms + xpxq2 * xmc );
   y = v * ( ( 1.0 - 2.0 * xq * xq ) * xmc - xpxq2 * xms );
   z = v * ( 2.0 * ci2 * ( xp * xms + xq * xmc ) );

/*  Rotate to equatorial */
   pv[3] = x;
   pv[4] = y * COSEPS - z * SINEPS;
   pv[5] = y * SINEPS + z * COSEPS;
}

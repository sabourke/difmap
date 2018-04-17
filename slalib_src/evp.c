#include "slalib.h"
#include "slamac.h"
void slaEvp ( double date, double deqx, double dvb[3], double dpb[3],
              double dvh[3], double dph[3] )
/*
**  - - - - - - -
**   s l a E v p
**  - - - - - - -
**
**  Barycentric and heliocentric velocity and position of the Earth.
**
**  Given:
**
**     date    double     TDB (loosely ET) as a Modified Julian Date
**                                         (JD-2400000.5)
**
**     deqx    double     Julian epoch (e.g. 2000.0) of mean equator and
**                        equinox of the vectors returned.  If deqx <= 0.0,
**                        all vectors are referred to the mean equator and
**                        equinox (FK5) of epoch date.
**
**  Returned (all 3D Cartesian vectors):
**
**     dvb,dpb double[3]  barycentric velocity, position
**
**     dvh,dph double[3]  heliocentric velocity, position
**
**  (Units are AU/s for velocity and AU for position)
**
**  Called:  slaEpj, slaPrec
**
**  Accuracy:
**
**     The maximum deviations from the JPL DE96 ephemeris are as
**     follows:
**
**     barycentric velocity                  42  cm/s
**     barycentric position                6900  km
**
**     heliocentric velocity                 42  cm/s
**     heliocentric position               1600  km
**
**  This routine is adapted from the BARVEL and BARCOR Fortran
**  subroutines of P.Stumpff, which are described in
**  Astron. Astrophys. Suppl. Ser. 41, 1-8 (1980).  The present
**  routine uses double precision throughout;  most of the other
**  changes are essentially cosmetic and do not affect the
**  results.  However, some adjustments have been made so as to
**  give results that refer to the new (IAU 1976 "FK5") equinox
**  and precession, although the differences these changes make
**  relative to the results from Stumpff's original "FK4" version
**  are smaller than the inherent accuracy of the algorithm.  One
**  minor shortcoming in the original routines that has not been
**  corrected is that better numerical accuracy could be achieved
**  if the various polynomial evaluations were nested.  Note also
**  that one of Stumpff's precession constants differs by 0.001 arcsec
**  from the value given in the Explanatory Supplement to the A.E.
**
**  Defined in slamac.h:  D2PI, DS2R, dmod
**
**  Last revision:   22 September 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int ideq, i, j, k;

   double a, pertl,
          pertld, pertr, pertrd, cosa, sina, e, twoe, esq, g, twog,
          phi, f, sf, cf, phid, psid, pertp, pertpd, tl, sinlm, coslm,
          sigma, b, plon, pomg, pecc, flatm, flat;

   double dt, dlocal, dml,
          deps, dparam, dpsi, d1pdro, drd, drld, dtl, dsinls,
          dcosls, dxhd, dyhd, dzhd, dxbd, dybd, dzbd, dcosep,
          dsinep, dyahd, dzahd, dyabd, dzabd, dr,
          dxh, dyh, dzh, dxb, dyb, dzb, dyah, dzah, dyab,
          dzab, depj, deqcor;

   double sn[4], forbel[7], sorbel[17], sinlp[4], coslp[4];

   double dprema[3][3], w, vw[3];

/* Sidereal rate dcsld in longitude, rate ccsgd in mean anomaly */
   static double dcsld = 1.990987e-7;
   static double ccsgd = 1.990969e-7;

/* Some constants used in the calculation of the lunar contribution */
   static double cckm  = 3.122140e-5;
   static double ccmld = 2.661699e-6;
   static double ccfdi = 2.399485e-7;

/* Besselian epoch 1950.0 expressed as a Julian epoch */
   static double b1950 = 1949.9997904423;

/*
** ccpamv(k)=a*m*dl/dt (planets), dc1mme=1-mass(Earth+Moon)
*/
   static double ccpamv[4] = {
      8.326827e-11,
      1.843484e-11,
      1.988712e-12,
      1.881276e-12
   };
   static double dc1mme = 0.99999696;

/*
** ccpam(k)=a*m(planets)
** ccim=inclination(Moon)
*/
   static double ccpam[4] = {
      4.960906e-3,
      2.727436e-3,
      8.392311e-4,
      1.556861e-3
   };
   static double ccim = 8.978749e-2;

/*
** Constants dcfel(i,k) of fast changing elements
*/
   static double dcfel[3][8] = {
      {  1.7400353,                /* dcfel[0][0] */
         6.2565836,                /* dcfel[0][1] */
         4.7199666,                /* dcfel[0][2] */
         1.9636505e-1,             /* dcfel[0][3] */
         4.1547339,                /* dcfel[0][4] */
         4.6524223,                /* dcfel[0][5] */
         4.2620486,                /* dcfel[0][6] */
         1.4740694 },              /* dcfel[0][7] */
      {  6.2833195099091e+2,       /* dcfel[1][0] */
         6.2830194572674e+2,       /* dcfel[1][1] */
         8.3997091449254e+3,       /* dcfel[1][2] */
         8.4334662911720e+3,       /* dcfel[1][3] */
         5.2993466764997e+1,       /* dcfel[1][4] */
         2.1354275911213e+1,       /* dcfel[1][5] */
         7.5025342197656,          /* dcfel[1][6] */
         3.8377331909193 },        /* dcfel[1][7] */
      {  5.2796e-6,                /* dcfel[2][0] */
        -2.6180e-6,                /* dcfel[2][1] */
        -1.9780e-5,                /* dcfel[2][2] */
        -5.6044e-5,                /* dcfel[2][3] */
         5.8845e-6,                /* dcfel[2][4] */
         5.6797e-6,                /* dcfel[2][5] */
         5.5317e-6,                /* dcfel[2][6] */
         5.6093e-6 }               /* dcfel[2][7] */
   };

/*
** Constants dceps and ccsel(i,k) of slowly changing elements
*/
   static double dceps[3] = {
      4.093198e-1,
     -2.271110e-4,
     -2.860401e-8
   };
   static double ccsel[3][17] = {
      {  1.675104e-2,              /* ccsel[0][0]  */
         2.220221e-1,              /* ccsel[0][1]  */
         1.589963,                 /* ccsel[0][2]  */
         2.994089,                 /* ccsel[0][3]  */
         8.155457e-1,              /* ccsel[0][4]  */
         1.735614,                 /* ccsel[0][5]  */
         1.968564,                 /* ccsel[0][6]  */
         1.282417,                 /* ccsel[0][7]  */
         2.280820,                 /* ccsel[0][8]  */
         4.833473e-2,              /* ccsel[0][9]  */
         5.589232e-2,              /* ccsel[0][10] */
         4.634443e-2,              /* ccsel[0][11] */
         8.997041e-3,              /* ccsel[0][12] */
         2.284178e-2,              /* ccsel[0][13] */
         4.350267e-2,              /* ccsel[0][14] */
         1.348204e-2,              /* ccsel[0][15] */
         3.106570e-2 },            /* ccsel[0][16] */
      { -4.179579e-5,              /* ccsel[1][0]  */
         2.809917e-2,              /* ccsel[1][1]  */
         3.418075e-2,              /* ccsel[1][2]  */
         2.590824e-2,              /* ccsel[1][3]  */
         2.486352e-2,              /* ccsel[1][4]  */
         1.763719e-2,              /* ccsel[1][5]  */
         1.524020e-2,              /* ccsel[1][6]  */
         8.703393e-3,              /* ccsel[1][7]  */
         1.918010e-2,              /* ccsel[1][8]  */
         1.641773e-4,              /* ccsel[1][9]  */
        -3.455092e-4,              /* ccsel[1][10] */
        -2.658234e-5,              /* ccsel[1][11] */
         6.329728e-6,              /* ccsel[1][12] */
        -9.941590e-5,              /* ccsel[1][13] */
        -6.839749e-5,              /* ccsel[1][14] */
         1.091504e-5,              /* ccsel[1][15] */
        -1.665665e-4 },            /* ccsel[1][16] */
      { -1.260516e-7,              /* ccsel[2][0]  */
         1.852532e-5,              /* ccsel[2][1]  */
         1.430200e-5,              /* ccsel[2][2]  */
         4.155840e-6,              /* ccsel[2][3]  */
         6.836840e-6,              /* ccsel[2][4]  */
         6.370440e-6,              /* ccsel[2][5]  */
        -2.517152e-6,              /* ccsel[2][6]  */
         2.289292e-5,              /* ccsel[2][7]  */
         4.484520e-6,              /* ccsel[2][8]  */
        -4.654200e-7,              /* ccsel[2][9]  */
        -7.388560e-7,              /* ccsel[2][10] */
         7.757000e-8,              /* ccsel[2][11] */
        -1.939256e-9,              /* ccsel[2][12] */
         6.787400e-8,              /* ccsel[2][13] */
        -2.714956e-7,              /* ccsel[2][14] */
         6.903760e-7,              /* ccsel[2][15] */
        -1.590188e-7 }             /* ccsel[2][16] */
   };

/*
** Constants of the arguments of the short-period perturbations
** by the planets:   dcargs(i,k)
*/
   static double dcargs[2][15] = {
      {  5.0974222,                /* dcargs[0][0]  */
         3.9584962,                /* dcargs[0][1]  */
         1.6338070,                /* dcargs[0][2]  */
         2.5487111,                /* dcargs[0][3]  */
         4.9255514,                /* dcargs[0][4]  */
         1.3363463,                /* dcargs[0][5]  */
         1.6072053,                /* dcargs[0][6]  */
         1.3629480,                /* dcargs[0][7]  */
         5.5657014,                /* dcargs[0][8]  */
         5.0708205,                /* dcargs[0][9]  */
         3.9318944,                /* dcargs[0][10] */
         4.8989497,                /* dcargs[0][11] */
         1.3097446,                /* dcargs[0][12] */
         3.5147141,                /* dcargs[0][13] */
         3.5413158 },              /* dcargs[0][14] */
      { -7.8604195454652e+2,       /* dcargs[1][0]  */
        -5.7533848094674e+2,       /* dcargs[1][1]  */
        -1.1506769618935e+3,       /* dcargs[1][2]  */
        -3.9302097727326e+2,       /* dcargs[1][3]  */
        -5.8849265665348e+2,       /* dcargs[1][4]  */
        -5.5076098609303e+2,       /* dcargs[1][5]  */
        -5.2237501616674e+2,       /* dcargs[1][6]  */
        -1.1790629318198e+3,       /* dcargs[1][7]  */
        -1.0977134971135e+3,       /* dcargs[1][8]  */
        -1.5774000881978e+2,       /* dcargs[1][9]  */
         5.2963464780000e+1,       /* dcargs[1][10] */
         3.9809289073258e+1,       /* dcargs[1][11] */
         7.7540959633708e+1,       /* dcargs[1][12] */
         7.9618578146517e+1,       /* dcargs[1][13] */
        -5.4868336758022e+2 }      /* dcargs[1][14] */
   };

/*
** Amplitudes ccamps(n,k) of the short-period perturbations
*/
   static double ccamps[5][15] = {
      { -2.279594e-5,              /* ccamps[0][0]  */
        -3.494537e-5,              /* ccamps[0][1]  */
         6.593466e-7,              /* ccamps[0][2]  */
         1.140767e-5,              /* ccamps[0][3]  */
         9.516893e-6,              /* ccamps[0][4]  */
         7.310990e-6,              /* ccamps[0][5]  */
        -2.603449e-6,              /* ccamps[0][6]  */
        -3.228859e-6,              /* ccamps[0][7]  */
         3.442177e-7,              /* ccamps[0][8]  */
         8.702406e-6,              /* ccamps[0][9]  */
        -1.488378e-6,              /* ccamps[0][10] */
        -8.043059e-6,              /* ccamps[0][11] */
         3.699128e-6,              /* ccamps[0][12] */
         2.550120e-6,              /* ccamps[0][13] */
        -6.351059e-7 },            /* ccamps[0][14] */
      {  1.407414e-5,              /* ccamps[1][0]  */
         2.860401e-7,              /* ccamps[1][1]  */
         1.322572e-5,              /* ccamps[1][2]  */
        -2.049792e-5,              /* ccamps[1][3]  */
        -2.748894e-6,              /* ccamps[1][4]  */
        -1.924710e-6,              /* ccamps[1][5]  */
         7.359472e-6,              /* ccamps[1][6]  */
         1.308997e-7,              /* ccamps[1][7]  */
         2.671323e-6,              /* ccamps[1][8]  */
        -8.421214e-6,              /* ccamps[1][9]  */
        -1.251789e-5,              /* ccamps[1][10] */
        -2.991300e-6,              /* ccamps[1][11] */
        -3.316126e-6,              /* ccamps[1][12] */
        -1.241123e-6,              /* ccamps[1][13] */
         2.341650e-6 },            /* ccamps[1][14] */
      {  8.273188e-6,              /* ccamps[2][0]  */
         1.289448e-7,              /* ccamps[2][1]  */
         9.258695e-6,              /* ccamps[2][2]  */
        -4.747930e-6,              /* ccamps[2][3]  */
        -1.319381e-6,              /* ccamps[2][4]  */
        -8.772849e-7,              /* ccamps[2][5]  */
         3.168357e-6,              /* ccamps[2][6]  */
         1.013137e-7,              /* ccamps[2][7]  */
         1.832858e-6,              /* ccamps[2][8]  */
        -1.372341e-6,              /* ccamps[2][9]  */
         5.226868e-7,              /* ccamps[2][10] */
         1.473654e-7,              /* ccamps[2][11] */
         2.901257e-7,              /* ccamps[2][12] */
         9.901116e-8,              /* ccamps[2][13] */
         1.061492e-6 },            /* ccamps[2][14] */
      {  1.340565e-5,              /* ccamps[3][0]  */
         1.627237e-5,              /* ccamps[3][1]  */
        -4.674248e-7,              /* ccamps[3][2]  */
        -2.638763e-6,              /* ccamps[3][3]  */
        -4.549908e-6,              /* ccamps[3][4]  */
        -3.334143e-6,              /* ccamps[3][5]  */
         1.119056e-6,              /* ccamps[3][6]  */
         2.403899e-6,              /* ccamps[3][7]  */
        -2.394688e-7,              /* ccamps[3][8]  */
        -1.455234e-6,              /* ccamps[3][9]  */
        -2.049301e-7,              /* ccamps[3][10] */
        -3.154542e-7,              /* ccamps[3][11] */
         3.407826e-7,              /* ccamps[3][12] */
         2.210482e-7,              /* ccamps[3][13] */
         2.878231e-7 },            /* ccamps[3][14] */
      { -2.490817e-7,              /* ccamps[4][0]  */
        -1.823138e-7,              /* ccamps[4][1]  */
        -3.646275e-7,              /* ccamps[4][2]  */
        -1.245408e-7,              /* ccamps[4][3]  */
        -1.864821e-7,              /* ccamps[4][4]  */
        -1.745256e-7,              /* ccamps[4][5]  */
        -1.655307e-7,              /* ccamps[4][6]  */
        -3.736225e-7,              /* ccamps[4][7]  */
        -3.478444e-7,              /* ccamps[4][8]  */
        -4.998479e-8,              /* ccamps[4][9]  */
         0.0,                      /* ccamps[4][10] */
         0.0,                      /* ccamps[4][11] */
         0.0,                      /* ccamps[4][12] */
         0.0,                      /* ccamps[4][13] */
         0.0 }                     /* ccamps[4][14] */
    };

/*
** Constants of the secular perturbations in longitude
** ccsec3 and ccsec(n,k)
*/
   static double ccsec3 = -7.757020e-8;
   static double ccsec[3][4] = {
      {  1.289600e-6,              /* ccsec[0][0] */
         3.102810e-5,              /* ccsec[0][1] */
         9.124190e-6,              /* ccsec[0][2] */
         9.793240e-7 },            /* ccsec[0][3] */
      {  5.550147e-1,              /* ccsec[1][0] */
         4.035027,                 /* ccsec[1][1] */
         9.990265e-1,              /* ccsec[1][2] */
         5.508259 },               /* ccsec[1][3] */
      {  2.076942,                 /* ccsec[2][0] */
         3.525565e-1,              /* ccsec[2][1] */
         2.622706,                 /* ccsec[2][2] */
         1.559103e+1 }             /* ccsec[2][3] */
   };

/*
** Constants dcargm(i,k) of the arguments of the perturbations
** of the motion of the Moon
*/
   static double dcargm[2][3] = {
      {  5.167983,                 /* dcargm[0][0] */
         5.491315,                 /* dcargm[0][1] */
         5.959853 },               /* dcargm[0][2] */
      {  8.3286911095275e+3,       /* dcargm[1][0] */
        -7.2140632838100e+3,       /* dcargm[1][1] */
         1.5542754389685e+4 }      /* dcargm[1][2] */
   };

/*
** Amplitudes ccampm(n,k) of the perturbations of the Moon
*/
   static double ccampm[4][3] = {
      {  1.097594e-1,              /* ccampm[0][0] */
        -2.223581e-2,              /* ccampm[0][1] */
         1.148966e-2 },            /* ccampm[0][2] */
      {  2.896773e-7,              /* ccampm[1][0] */
         5.083103e-8,              /* ccampm[1][1] */
         5.658888e-8 },            /* ccampm[1][2] */
      {  5.450474e-2,              /* ccampm[2][0] */
         1.002548e-2,              /* ccampm[2][1] */
         8.249439e-3 },            /* ccampm[2][2] */
      {  1.438491e-7,              /* ccampm[3][0] */
        -2.291823e-8,              /* ccampm[3][1] */
         4.063015e-8 }             /* ccampm[3][2] */
   };

/*
**
** Execution
** ---------
**
** Control parameter ideq, and time arguments
*/
   ideq = ( deqx <= 0.0 ) ? 0 : 1;
   dt = ( date - 15019.5 ) / 36525.0;

/* Values of all elements for the instant date */
   for ( k = 0; k < 8; k++ ) {
      dlocal = dmod ( dcfel[0][k]
             + dt * ( dcfel[1][k]
               + dt * dcfel[2][k] ), D2PI );
      if ( k == 0 ) {
         dml = dlocal;
      } else {
         forbel[k-1] = dlocal;
      }
   }
   deps = dmod ( dceps[0]
        + dt * ( dceps[1]
          + dt * dceps[2] ) , D2PI );
   for ( k = 0; k < 17; k++ ) {
      sorbel[k] = dmod ( ccsel[0][k]
                + dt * ( ccsel[1][k]
                  + dt * ccsel[2][k] ), D2PI );
   }

/* Secular perturbations in longitude */
   for ( k = 0; k < 4; k++ ) {
      a = dmod ( ccsec[1][k] + dt * ccsec[2][k] , D2PI );
      sn[k] = sin ( a );
   }

/* Periodic perturbations of the EMB (Earth-Moon barycentre) */
   pertl = ccsec[0][0] * sn[0]
         + ccsec[0][1] * sn[1]
       + ( ccsec[0][2] + dt * ccsec3 ) * sn[2]
         + ccsec[0][3] * sn[3];
   pertld = 0.0;
   pertr = 0.0;
   pertrd = 0.0;
   for ( k = 0; k < 15; k++ ) {
      a = dmod ( dcargs[0][k] + dt * dcargs[1][k] , D2PI );
      cosa = cos ( a );
      sina = sin ( a );
      pertl = pertl + ccamps[0][k] * cosa + ccamps[1][k] * sina;
      pertr = pertr + ccamps[2][k] * cosa + ccamps[3][k] * sina;
      if ( k < 10 ) {
         pertld = pertld +
               ( ccamps[1][k] * cosa - ccamps[0][k] * sina ) * ccamps[4][k];
         pertrd = pertrd +
               ( ccamps[3][k] * cosa - ccamps[2][k] * sina ) * ccamps[4][k];
      }
   }

/* Elliptic part of the motion of the EMB */
   e = sorbel[0];
   twoe = e + e;
   esq = e * e;
   dparam = 1.0 - esq;
   g = forbel[0];
   twog = g + g;
   phi = twoe * ( ( 1.0 - esq / 8.0 ) * sin ( g )
                + 5.0 * e * sin ( twog ) / 8.0
                + 13.0 * esq * sin ( g + twog ) / 24.0 );
   f = forbel[0] + phi;
   sf = sin ( f );
   cf = cos ( f );
   dpsi = dparam / ( 1.0 + e * cf );
   phid = twoe * ccsgd * ( ( 1.0 + esq * 1.5 ) * cf
                         + e * ( 1.25 - sf * sf / 2.0 ) );
   psid = ccsgd * e * sf / sqrt ( dparam );

/* Perturbed heliocentric motion of the EMB */
   d1pdro = 1.0 + pertr;
   drd = d1pdro * ( psid + dpsi * pertrd );
   drld = d1pdro * dpsi * ( dcsld + phid + pertld );
   dtl = dmod ( dml + phi + pertl , D2PI );
   dsinls = sin ( dtl );
   dcosls = cos ( dtl );
   dxhd = drd * dcosls - drld * dsinls;
   dyhd = drd * dsinls + drld * dcosls;

/* Influence of eccentricity, evection and variation on the
** geocentric motion of the Moon */
   pertl = 0.0;
   pertld = 0.0;
   pertp = 0.0;
   pertpd = 0.0;
   for ( k = 0; k < 3; k++ ) {
      a = dmod ( dcargm[0][k] + dt * dcargm[1][k] , D2PI );
      sina = sin ( a );
      cosa = cos ( a );
      pertl = pertl + ccampm[0][k] * sina;
      pertld = pertld + ccampm[1][k] * cosa;
      pertp = pertp + ccampm[2][k] * cosa;
      pertpd = pertpd - ccampm[3][k] * sina;
   }

/* Heliocentric motion of the Earth */
   tl = forbel[1] + pertl;
   sinlm = sin ( tl );
   coslm = cos ( tl );
   sigma = cckm / ( 1.0 + pertp );
   a = sigma * ( ccmld + pertld );
   b = sigma * pertpd;
   dxhd  = dxhd + a * sinlm + b * coslm;
   dyhd  = dyhd - a * coslm + b * sinlm;
   dzhd  = - sigma * ccfdi * cos ( forbel[2] );

/* Barycentric motion of the Earth */
   dxbd = dxhd * dc1mme;
   dybd = dyhd * dc1mme;
   dzbd = dzhd * dc1mme;
   for ( k = 0; k < 4; k++ ) {
      plon = forbel[k+3];
      pomg = sorbel[k+1];
      pecc = sorbel[k+9];
      tl = dmod( plon + 2.0 * pecc * sin ( plon - pomg ) , D2PI );
      sinlp[k] = sin ( tl );
      coslp[k] = cos ( tl );
      dxbd = dxbd + ccpamv[k] * ( sinlp[k] + pecc * sin ( pomg ) );
      dybd = dybd - ccpamv[k] * ( coslp[k] + pecc * cos ( pomg ) );
      dzbd = dzbd - ccpamv[k] * sorbel[k+13] * cos ( plon - sorbel[k+5] );
   }

/* Transition to mean equator of date */
   dcosep = cos ( deps );
   dsinep = sin ( deps );
   dyahd  = dcosep * dyhd - dsinep * dzhd;
   dzahd  = dsinep * dyhd + dcosep * dzhd;
   dyabd  = dcosep * dybd - dsinep * dzbd;
   dzabd  = dsinep * dybd + dcosep * dzbd;

/* Heliocentric coordinates of the Earth */
   dr = dpsi * d1pdro;
   flatm = ccim * sin ( forbel[2] );
   a = sigma * cos ( flatm );
   dxh = dr * dcosls - a * coslm;
   dyh = dr * dsinls - a * sinlm;
   dzh = - sigma * sin ( flatm );

/* Barycentric coordinates of the Earth */
   dxb = dxh * dc1mme;
   dyb = dyh * dc1mme;
   dzb = dzh * dc1mme;
   for ( k = 0; k < 4; k++ ) {
      flat = sorbel[k+13] * sin ( forbel[k+3] - sorbel[k+5] );
      a = ccpam[k] * (1.0 - sorbel[k+9] * cos ( forbel[k+3] - sorbel[k+1]));
      b = a * cos(flat);
      dxb -= b * coslp[k];
      dyb -= b * sinlp[k];
      dzb -= a * sin ( flat );
   }

/* Transition to mean equator of date */
   dyah = dcosep * dyh - dsinep * dzh;
   dzah = dsinep * dyh + dcosep * dzh;
   dyab = dcosep * dyb - dsinep * dzb;
   dzab = dsinep * dyb + dcosep * dzb;

/* Copy result components into vectors, correcting for FK4 equinox */
   depj = slaEpj ( date );
   deqcor = DS2R * ( 0.035 + ( 0.00085 * ( depj - b1950 ) ) );
   dvh[0] = dxhd - deqcor * dyahd;
   dvh[1] = dyahd + deqcor * dxhd;
   dvh[2] = dzahd;
   dvb[0] = dxbd - deqcor * dyabd;
   dvb[1] = dyabd + deqcor * dxbd;
   dvb[2] = dzabd;
   dph[0] = dxh - deqcor * dyah;
   dph[1] = dyah + deqcor * dxh;
   dph[2] = dzah;
   dpb[0] = dxb - deqcor * dyab;
   dpb[1] = dyab + deqcor * dxb;
   dpb[2] = dzab;

/* Was precession to another equinox requested? */
   if ( ideq != 0 ) {

   /* Yes: compute precession matrix from MJD date to Julian Epoch deqx */
      slaPrec ( depj, deqx, dprema );

   /* Rotate dvh */
      for ( j = 0; j < 3; j++ ) {
         w = 0.0;
         for ( i = 0; i < 3; i++ ) {
            w += dprema[j][i] * dvh[i];
         }
         vw[j] = w;
      }
      for ( j = 0; j < 3; j++ ) {
         dvh[j] = vw[j];
      }

   /* Rotate dvb */
      for ( j = 0; j < 3; j++ ) {
         w = 0.0;
         for ( i = 0; i < 3; i++ ) {
            w += dprema[j][i] * dvb[i];
         }
         vw[j] = w;
      }
      for ( j = 0; j < 3; j++ ) {
         dvb[j] = vw[j];
      }

   /* Rotate dph */
      for ( j = 0; j < 3; j++ ) {
         w = 0.0;
         for ( i = 0; i < 3; i++ ) {
            w += dprema[j][i] * dph[i];
         }
         vw[j] = w;
      }
      for ( j = 0; j < 3; j++ ) {
         dph[j] = vw[j];
      }

   /* Rotate dpb */
      for ( j = 0; j < 3; j++ ) {
         w = 0.0;
         for ( i = 0; i < 3; i++ ) {
            w += dprema[j][i] * dpb[i];
         }
         vw[j] = w;
      }
      for ( j = 0; j < 3; j++ ) {
         dpb[j] = vw[j];
      }
   }
}

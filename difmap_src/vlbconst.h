/* C include file containing commonly used constants */

/* The speed of light */

#define cvel 2.99792458e8

/* pi */

#define pi 3.1415926535897932384626433832795028841971693993751

/* 2*pi */

#define twopi 6.2831853071795864769252867665590057683943387987502

/* pi/2 */

#define halfpi 1.5707963267948966192313216916397514420985846996875

/* Scale radians to degrees */

#define rtod 57.29577951308232087679815481410517033240547246656458

/* Scale arcseconds to radians */

#define astor 4.8481368110953599358991410235794797595635e-6

/* Scale milliarcseconds to radians */

#define mastor 4.8481368110953599358991410235794797595635e-9

/* Scale radians to milliarcsec */

#define rtomas 2.062648062470963551564733573307786131966597008796332528822e+8

/* Scale radians to arcseconds */

#define rtoas 2.062648062470963551564733573307786131966597008796332528822e+5

/* Scale radians to arcminutes */

#define rtoam  3.437746770784939252607889288846310219944328347993887548e+3

/* Scale degrees to radians */

#define dtor 0.0174532925199432957692369076848861271344287188854

/* Scale degrees to milliarcsec */

#define dtomas 3.6e6

/* Scale radians to hours */

#define rtoh 3.819718634205488058453210320940344688827031497771

/* Scale hours to radians */

#define htor 0.2617993877991494365385536152732919070164307832813

/* The number of seconds in a day */

#define daysec 86400.0

/* The number of minutes in a day */

#define daymin 1440.0

/* Mean sidereal seconds per UT1 second */

#define ut_to_mst 1.002737909350795

/* UT1 seconds per mean sidereal second */

#define mst_to_ut 0.997269566329084

/*
 * UT1 days per sidereal earth rotation. This differs from ut_to_mst
 * because it includes the rate of precession of right ascension.
 */
#define ut_to_rot 0.9972696632424

/*
 * Sidereal earth rotations per UT1 day. This differs from mst_to_ut
 * because it includes the rate of precession of right ascension.
 */
#define rot_to_ut 1.002737811906

/* Scale Astronomical Units to meters */

#define au_to_m 1.49597870e11

/* Scale meters to Astronomical units */

#define m_to_au (1.0/1.49597870e11)

/* Minutes per internal time unit */

#define uttomin (1.0/60.0)

/* Boltzmann's constant (Joules/Kelvin) */

#define boltzmann 1.38066e-23

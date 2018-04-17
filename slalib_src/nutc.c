#include "slalib.h"
#include "slamac.h"
void slaNutc ( double date, double *dpsi, double *deps, double *eps0 )
/*
**  - - - - - - - -
**   s l a N u t c
**  - - - - - - - -
**
**  Nutation:  longitude & obliquity components and
**             mean obliquity (IAU 1980 theory).
**
**  (double precision)
**
**  References:
**     Final report of the IAU working group on nutation,
**      chairman P.K.Seidelmann, 1980.
**     Kaplan,G.H., 1981, USNO circular No. 163, pa3-6.
**
**  Given:
**     date        double    TDB (loosely ET) as Modified Julian Date
**                                            (JD-2400000.5)
**
**  Returned:
**     *dpsi,*deps double    nutation in longitude,obliquity
**     *eps0       double    mean obliquity
**
**  Defined in slamac.h:  DAS2R, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/

#define T2AS 1296000.0                /* Turns to arc seconds */
#define U2R 0.4848136811095359949e-9  /* Units of 0.0001 arcsec to radians */

{
   double t, el, el2, el3, elp, elp2,
          f, f2, f4,
          d, d2, d4,
          om, om2,
          dp, de, a;

/* Interval between basic epoch J2000.0 and current epoch (JC) */
   t = ( date - 51544.5 ) / 36525.0;

/* Fundamental arguments in the FK5 reference system */

/* Mean longitude of the Moon minus mean longitude of the Moon's perigee */
   el = DAS2R * dmod ( 485866.733 + ( 1325.0 * T2AS + 715922.633
                   + ( 31.310 + 0.064 * t ) * t ) * t , T2AS );

/* Mean longitude of the Sun minus mean longitude of the Sun's perigee */
   elp = DAS2R * dmod ( 1287099.804 + ( 99.0 * T2AS + 1292581.224
                    + ( -0.577 - 0.012 * t ) * t ) * t, T2AS );

/* Mean longitude of the Moon minus mean longitude of the Moon's node */
   f = DAS2R * dmod ( 335778.877 + ( 1342.0 * T2AS + 295263.137
                  + ( -13.257 + 0.011 * t ) * t ) * t, T2AS );

/* Mean elongation of the Moon from the Sun */
   d = DAS2R * dmod ( 1072261.307 + ( 1236.0 * T2AS + 1105601.328
                  + ( -6.891 + 0.019 * t ) * t ) * t, T2AS );

/* Longitude of the mean ascending node of the lunar orbit on the
   ecliptic, measured from the mean equinox of date */
   om = DAS2R * dmod ( 450160.280 + ( -5.0 * T2AS - 482890.539
                   + ( 7.455 + 0.008 * t ) * t ) * t, T2AS );

/* Multiples of arguments */
   el2 = el + el;
   el3 = el2 + el;
   elp2 = elp + elp;
   f2 = f + f;
   f4 = f2 + f2;
   d2 = d + d;
   d4 = d2 + d2;
   om2 = om + om;

/* Series for the nutation */
   dp = 0.0;
   de = 0.0;

   dp += sin ( elp + d );                          /* 106  */

   dp -= sin ( f2 + d4 + om2 );                    /* 105  */

   dp += sin ( el2 + d2 );                         /* 104  */

   dp -= sin ( el - f2 + d2 );                     /* 103  */

   dp -= sin ( el + elp - d2 + om );               /* 102  */

   dp -= sin ( - elp + f2 + om );                  /* 101  */

   dp -= sin ( el - f2 - d2 );                     /* 100  */

   dp -= sin ( elp + d2 );                         /*  99  */

   dp -= sin ( f2 - d + om2 );                     /*  98  */

   dp -= sin ( - f2 + om );                        /*  97  */

   dp += sin ( - el - elp + d2 + om );             /*  96  */

   dp += sin ( elp + f2 + om );                    /*  95  */

   dp -= sin ( el + f2 - d2 );                     /*  94  */

   dp += sin ( el3 + f2 - d2 + om2 );              /*  93  */

   dp += sin ( f4 - d2 + om2 );                    /*  92  */

   dp -= sin ( el + d2 + om );                     /*  91  */

   dp -= sin ( el2 + f2 + d2 + om2 );              /*  90  */

   a = el2 + f2 - d2 + om;                         /*  89  */
   dp += sin ( a );
   de -= cos ( a );

   dp += sin ( el - elp - d2 );                    /*  88  */

   dp += sin ( - el + f4 + om2 );                  /*  87  */

   a = - el2 + f2 + d4 + om2;                      /*  86  */
   dp -= sin ( a );
   de += cos ( a );

   a  = el + f2 + d2 + om;                         /*  85  */
   dp -= sin ( a );
   de += cos ( a );

   a = el + elp + f2 - d2 + om2;                   /*  84  */
   dp += sin ( a );
   de -= cos ( a );

   dp -= sin ( el2 - d4 );                         /*  83  */

   a = - el + f2 + d4 + om2;                       /*  82  */
   dp -= 2.0 * sin ( a );
   de += cos ( a );

   a = - el2 + f2 + d2 + om2;                      /*  81  */
   dp += sin ( a );
   de = de - cos ( a );

   dp -= sin ( el - d4 );                          /*  80  */

   a = - el + om2;                                 /*  79  */
   dp += sin ( a );
   de = de - cos ( a );

   a = f2 + d + om2;                               /*  78  */
   dp += 2.0 * sin ( a );
   de = de - cos ( a );

   dp += 2.0 * sin ( el3 );                        /*  77  */

   a = el + om2;                                   /*  76  */
   dp -= 2.0 * sin ( a );
   de += cos ( a );

   a = el2 + om;                                   /*  75  */
   dp += 2.0 * sin ( a );
   de -= cos ( a );

   a = - el + f2 - d2 + om;                        /*  74  */
   dp -= 2.0 * sin ( a );
   de += cos ( a );

   a = el + elp + f2 + om2;                        /*  73  */
   dp += 2.0 * sin ( a );
   de = de - cos ( a );

   a = - elp + f2 + d2 + om2;                      /*  72  */
   dp -= 3.0 * sin ( a );
   de += cos ( a );

   a = el3 + f2 + om2;                             /*  71  */
   dp -= 3.0 * sin ( a );
   de += cos ( a );

   a = - el2 + om;                                 /*  70  */
   dp -= 2.0 * sin ( a );
   de += cos ( a );

   a = - el - elp + f2 + d2 + om2;                 /*  69  */
   dp -= 3.0 * sin ( a );
   de += cos ( a );

   a = el - elp + f2 + om2;                        /*  68  */
   dp -= 3.0 * sin ( a );
   de += cos ( a );

   dp += 3.0 * sin ( el + f2 );                    /*  67  */

   dp -= 3.0 * sin ( el + elp );                   /*  66  */

   dp -= 4.0 * sin ( d );                          /*  65  */

   dp += 4.0 * sin ( el - f2 );                    /*  64  */

   dp -= 4.0 * sin ( elp - d2 );                   /*  63  */

   a = el2 + f2 + om;                              /*  62  */
   dp -= 5.0 * sin ( a );
   de += 3.0 * cos ( a );

   dp += 5.0 * sin ( el - elp );                   /*  61  */

   a = - d2 + om;                                  /*  60  */
   dp -= 5.0 * sin ( a );
   de += 3.0 * cos ( a );

   a = el + f2 - d2 + om;                          /*  59  */
   dp += 6.0 * sin ( a );
   de -= 3.0 * cos ( a );

   a = f2 + d2 + om;                               /*  58  */
   dp -= 7.0 * sin ( a );
   de += 3.0 * cos ( a );

   a = d2 + om;                                    /*  57  */
   dp -= 6.0 * sin ( a );
   de += 3.0 * cos ( a );

   a = el2 + f2 - d2 + om2;                        /*  56  */
   dp += 6.0 * sin ( a );
   de -= 3.0 * cos ( a );

   dp += 6.0 * sin ( el + d2);                     /*  55  */

   a = el + f2 + d2 + om2;                         /*  54  */
   dp -= 8.0 * sin ( a );
   de += 3.0 * cos ( a );

   a = - elp + f2 + om2;                           /*  53  */
   dp -= 7.0 * sin ( a );
   de += 3.0 * cos ( a );

   a = elp + f2 + om2;                             /*  52  */
   dp += 7.0 * sin ( a );
   de -= 3.0 * cos ( a );

   dp -= 7.0 * sin ( el + elp - d2 );              /*  51  */

   a = - el + f2 + d2 + om;                        /*  50  */
   dp -= 10.0 * sin ( a );
   de += 5.0 * cos ( a );

   a = el - d2 + om;                               /*  49  */
   dp -= 13.0 * sin ( a );
   de += 7.0 * cos ( a );

   a = - el + d2 + om;                             /*  48  */
   dp += 16.0 * sin ( a );
   de -= 8.0 * cos ( a );

   a = - el + f2 + om;                             /*  47  */
   dp += 21.0 * sin ( a );
   de -= 10.0 * cos ( a );

   dp += 26.0 * sin ( f2 );                        /*  46  */
   de -= cos( f2 );

   a = el2 + f2 + om2;                             /*  45  */
   dp -= 31.0 * sin ( a );
   de += 13.0 * cos ( a );

   a = el + f2 - d2 + om2;                         /*  44  */
   dp += 29.0 * sin ( a );
   de -= 12.0 * cos ( a );

   dp += 29.0 * sin ( el2 );                       /*  43  */
   de -= cos( el2 );

   a = f2 + d2 + om2;                              /*  42  */
   dp -= 38.0 * sin ( a );
   de += 16.0 * cos ( a );

   a = el + f2 + om;                               /*  41  */
   dp -= 51.0 * sin ( a );
   de += 27.0 * cos ( a );

   a = - el + f2 + d2 + om2;                       /*  40  */
   dp -= 59.0 * sin ( a );
   de += 26.0 * cos ( a );

   a = - el + om;                                  /*  39  */
   dp += ( - 58.0 -  0.1 * t ) * sin ( a );
   de += 32.0 * cos ( a );

   a = el + om;                                    /*  38  */
   dp += ( 63.0 + 0.1 * t ) * sin ( a );
   de -= 33.0 * cos ( a );

   dp += 63.0 * sin ( d2 );                        /*  37  */
   de -= 2.0 * cos( d2 );

   a = - el + f2 + om2;                            /*  36  */
   dp += 123.0 * sin ( a );
   de -= 53.0 * cos ( a );

   a = el - d2;                                    /*  35  */
   dp -= 158.0 * sin ( a );
   de -= cos ( a );

   a = el + f2 + om2;                              /*  34  */
   dp -= 301.0 * sin ( a );
   de += ( 129.0 - 0.1 * t ) * cos ( a );

   a = f2 + om;                                    /*  33  */
   dp += ( - 386.0 - 0.4 * t ) * sin ( a );
   de += 200.0 * cos ( a );

   dp += ( 712.0 + 0.1 * t ) * sin ( el );         /*  32  */
   de -= 7.0 * cos( el );

   a = f2 + om2;                                   /*  31  */
   dp += ( -2274.0 - 0.2 * t ) * sin ( a );
   de += ( 977.0 - 0.5 * t ) * cos ( a );

   dp -= sin ( elp + f2 - d2 );                    /*  30  */

   dp += sin ( - el + d + om );                    /*  29  */

   dp += sin ( elp + om2 );                        /*  28  */

   dp -= sin ( elp - f2 + d2 );                    /*  27  */

   dp += sin ( - f2 + d2 + om );                   /*  26  */

   dp += sin ( el2 + elp - d2 );                   /*  25  */

   dp -= 4.0 * sin ( el - d );                     /*  24  */

   a = elp + f2 - d2 + om;                         /*  23  */
   dp += 4.0 * sin ( a );
   de -= 2.0 * cos ( a );

   a = el2 - d2 + om;                              /*  22  */
   dp += 4.0 * sin ( a );
   de -= 2.0 * cos ( a );

   a = - elp + f2 - d2 + om;                       /*  21  */
   dp -= 5.0 * sin ( a );
   de += 3.0 * cos ( a );

   a = - el2 + d2 + om;                            /*  20  */
   dp -= 6.0 * sin ( a );
   de += 3.0 * cos ( a );

   a = - elp + om;                                 /*  19  */
   dp -= 12.0 * sin ( a );
   de += 6.0 * cos ( a );

   a = elp2 + f2 - d2 + om2;                       /*  18  */
   dp += ( - 16.0 + 0.1 * t) * sin ( a );
   de += 7.0 * cos ( a );

   a = elp + om;                                   /*  17  */
   dp -= 15.0 * sin ( a );
   de += 9.0 * cos ( a );

   dp += ( 17.0 - 0.1 * t ) * sin ( elp2 );        /*  16  */

   dp -= 22.0 * sin ( f2 - d2 );                   /*  15  */

   a = el2 - d2;                                   /*  14  */
   dp += 48.0 * sin ( a );
   de += cos ( a );

   a = f2 - d2 + om;                               /*  13  */
   dp += ( 129.0 + 0.1 * t ) * sin ( a );
   de -= 70.0 * cos ( a );

   a = - elp + f2 - d2 + om2;                      /*  12  */
   dp += ( 217.0 - 0.5 * t ) * sin ( a );
   de += ( -95.0 + 0.3 * t ) * cos ( a );

   a = elp + f2 - d2 + om2;                        /*  11  */
   dp += ( - 517.0 + 1.2 * t ) * sin ( a );
   de += ( 224.0 - 0.6 * t ) * cos ( a );

   dp += ( 1426.0 - 3.4 * t ) * sin ( elp );       /*  10  */
   de += ( 54.0 - 0.1 * t) * cos ( elp );

   a = f2 - d2 + om2;                              /*   9  */
   dp += ( - 13187.0 - 1.6 * t ) * sin ( a );
   de += ( 5736.0 - 3.1 * t ) * cos ( a );

   dp += sin ( el2 - f2 + om );                    /*   8  */

   a = - elp2 + f2 - d2 + om;                      /*   7  */
   dp -= 2.0 * sin ( a );
   de +=       cos ( a );

   dp -= 3.0 * sin ( el - elp - d );               /*   6  */

   a = - el2 + f2 + om2;                           /*   5  */
   dp -= 3.0 * sin ( a );
   de +=       cos ( a );

   dp += 11.0 * sin ( el2 - f2 );                  /*   4  */

   a = - el2 + f2 + om;                            /*   3  */
   dp += 46.0 * sin ( a );
   de -= 24.0 * cos ( a );

   dp += ( 2062.0 + 0.2 * t ) * sin ( om2 );       /*   2  */
   de += ( - 895.0 + 0.5 * t ) * cos ( om2 );

   dp += ( - 171996.0 - 174.2 * t) * sin ( om );   /*   1  */
   de += ( 92025.0 + 8.9 * t ) * cos ( om );

/* Convert results to radians */
   *dpsi = dp * U2R;
   *deps = de * U2R;

/* Mean obliquity */
   *eps0 = DAS2R * ( 84381.448 +
                   ( - 46.8150 +
                   ( - 0.00059 + 0.001813 * t ) * t ) * t );
}

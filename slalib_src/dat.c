#include "slalib.h"
#include "slamac.h"
double slaDat ( double utc )
/*
**  - - - - - - -
**   s l a D a t
**  - - - - - - -
**
**  Increment to be applied to Coordinated Universal Time UTC to give
**  International Atomic Time TAI.
**
**  (double precision)
**
**  Given:
**     utc      double      UTC date as a modified JD (JD-2400000.5)
**
**  Result:  TAI-UTC in seconds
**
**  Notes:
**
**  1  The UTC is specified to be a date rather than a time to indicate
**     that care needs to be taken not to specify an instant which lies
**     within a leap second.  Though in most cases UTC can include the
**     fractional part, correct behaviour on the day of a leap second
**     can only be guaranteed up to the end of the second 23:59:59.
**
**  2  Pre 1972 January 1 a fixed value of 10 sec is returned.
**
**     :-----------------------------------------:
**     :                                         :
**     :                IMPORTANT                :
**     :                                         :
**     :  This routine must be updated on each   :
**     :     occasion that a leap second is      :
**     :                announced                :
**     :                                         :
**     :  Latest leap second:  1996 January 1    :
**     :                                         :
**     :-----------------------------------------:
**
**  Last revision:   14 November 1995
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int i;
   static int npre = 10;  /* Increment to be added to UTC to give TAI */
                          /* prior to the first UTC entry in leap_utc */
/*
** Declare a table of the UTCs at which leap seconds were announced.
** Append the UTC of each new leap second to this table as they are announced.
** One second is added for each entry in the table.
*/
   static double leap_utc[]={
    41499.0, /* 1972 July 1    TAI-UTC = 11.0 */
    41683.0, /* 1973 January 1 TAI-UTC = 12.0 */
    42048.0, /* 1974 January 1 TAI-UTC = 13.0 */
    42413.0, /* 1975 January 1 TAI-UTC = 14.0 */
    42778.0, /* 1976 January 1 TAI-UTC = 15.0 */
    43144.0, /* 1977 January 1 TAI-UTC = 16.0 */
    43509.0, /* 1978 January 1 TAI-UTC = 17.0 */
    43874.0, /* 1979 January 1 TAI-UTC = 18.0 */
    44239.0, /* 1980 January 1 TAI-UTC = 19.0 */
    44786.0, /* 1981 July 1    TAI-UTC = 20.0 */
    45151.0, /* 1982 July 1    TAI-UTC = 21.0 */
    45516.0, /* 1983 July 1    TAI-UTC = 22.0 */
    46247.0, /* 1985 July 1    TAI-UTC = 23.0 */
    47161.0, /* 1988 January 1 TAI-UTC = 24.0 */
    47892.0, /* 1990 January 1 TAI-UTC = 25.0 */
    48257.0, /* 1991 January 1 TAI-UTC = 26.0 */
    48804.0, /* 1992 July 1    TAI-UTC = 27.0 */
    49169.0, /* 1993 July 1    TAI-UTC = 28.0 */
    49534.0, /* 1994 July 1    TAI-UTC = 29.0 */
    50083.0, /* 1996 January 1 TAI-UTC = 30.0 */
    50630.0, /* 1997 July 1    TAI-UTC = 31.0 */
    51179.0, /* 1999 Jan 1     TAI-UTC = 32.0 */
    53736.0, /* 2006 Jan 1     TAI-UTC = 33.0 */
   };

/* Record the number of entries in the table */
   static int num_leap = sizeof ( leap_utc ) / sizeof ( double );

/* Find the date of the last leap second that preceded the requested UTC */
   for ( i=0;  i < num_leap && utc > leap_utc[i];  i++ );

/* Return TAI-UTC for the specified date */
   return (double) ( npre + i );
}

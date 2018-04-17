#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include "logio.h"
#include "vlbconst.h"
#include "vlbutil.h"
#include "slalib.h"

/*
 * Define calendar as days per month and name of month (0 - relative indexing)
 */

static char daytab[2][12] = {
  {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
  {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

static char *month_name[] = {
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};

/*.......................................................................
 * Convert a number in radians to hours, minutes and seconds. The
 * resulting hours will be constrained to be within 0 -> 24 hours.
 */
void radhms(double rad, int *hour, int *mins, double *secs)
{
  double hr;               /* rad converted to decimal hours */
  double hr_int,hr_frc;    /* Integral and fractional hours */
  double min_int,min_frc;  /* Integral and fractional minutes */
/*
 * Convert radians to decimal hours.
 */
  hr = rad * rtoh;
/*
 * Enforce hr to be between 0 -> 24 hours.
 */
  if(hr > 24.0) {
    hr = fmod(hr, 24.0);
  }
  else if(hr < 0.0) {
    hr = 24.0 - fmod(fabs(hr), 24.0);
    if(hr >= 24.0)
      hr = 0.0;
  };
/*
 * Split decimal hours into integral and fractional parts.
 */
  hr_frc = modf(hr, &hr_int);
/*
 * Split the fractional hours into fractional and integral minutes.
 */
  min_frc = modf(hr_frc * 60.0, &min_int);
/*
 * Assign return values.
 */
  *hour = hr_int;
  *mins = min_int;
  *secs = min_frc * 60.0;
  return;
}

/*.......................................................................
 * Given an angle in radians, compose a string of form HH MM SS.S... in
 * a provided output string, where HH is the number of degrees, MM the
 * number of minutes and SS.S... the number of seconds to a given
 * precision.
 *
 * Input:
 *  rad    double   The angle, measured in radians.
 *  precision int   The number of decimal places to represent the seconds
 *                  argument with.
 *  colon     int   If true, separate sexagesimal components with a colon.
 * Input/Output:
 *  string   char * The output string. This must be at least
 *                  10 + 'precision' characters long.
 */
char *sradhms(double rad, int precision, int colon, char *string)
{
  int width;   /* The width of the seconds field */
  int hour;    /* The integral number of hours */
  int mins;    /* The integral number of minutes */
  double secs; /* The decimal number of seconds */
/*
 * Sanity checks.
 */
  if(precision < 0)
    precision = 0;
/*
 * Hand the conversion job to radhms().
 */
  radhms(rad, &hour, &mins, &secs);
/*
 * Determine the width of the seconds field. Note that this includes
 * the decimal point character when precision != 0.
 */
  width = precision == 0 ? 2 : (3 + precision);
/*
 * Compose the string.
 */
  sprintf(string, "%2.2d%c%2.2d%c%0*.*f", hour, colon ? ':':' ',
	  mins, colon ? ':':' ', width, precision, secs);
  return string;
}

/*.......................................................................
 * Convert a number in radians to sign, degrees, minutes and seconds.
 */
void raddms(double rad, int *sgn, int *deg, int *mins, double *secs)
{
  double deg_int,deg_frc;  /* Integral and fractional degrees */
  double min_int,min_frc;  /* Integral and fractional arcminutes */
/*
 * Split decimal degrees into integral and fractional parts.
 */
  deg_frc = modf(fabs(rad * rtod), &deg_int);
/*
 * Split the fractional degrees into fractional and integral minutes.
 */
  min_frc = modf(deg_frc * 60.0, &min_int);
/*
 * Assign the return values.
 */
  *sgn = rad >= 0.0 ? 1 : -1;
  *deg = deg_int;
  *mins = min_int;
  *secs = min_frc * 60.0;
  return;
}

/*.......................................................................
 * Given an angle in radians between -pi/2 and pi/2, compose a string of
 * form DDD MM SS.S... in a provided output string, where DDD is the
 * number of degrees including the sign, MM the number of arc-minutes
 * and SS.S... the number of arc-seconds to a given precision.
 *
 * Input:
 *  rad    double   The angle, measured in radians (must be in the range
 *                  -pi/2 to pi/2).
 *  precision int   The number of decimal places to represent the seconds
 *                  argument with.
 *  colon     int   If true, separate sexagesimal components with a colon.
 * Input/Output:
 *  string   char * The output string. This must be at least
 *                  11 + 'precision' characters long.
 */
char *sraddms(double rad, int precision, int colon, char *string)
{
  int width;   /* The width of the seconds field */
  int sgn;     /* The sign of the number */
  int deg;     /* The integral number of degrees */
  int mins;    /* The integral number of minutes */
  double secs; /* The decimal number of seconds */
/*
 * Sanity checks.
 */
  if(precision < 0)
    precision = 0;
/*
 * Hand the conversion job to radhms().
 */
  raddms(rad, &sgn, &deg, &mins, &secs);
/*
 * Determine the width of the seconds field. Note that this includes
 * the decimal point character when precision != 0.
 */
  width = precision == 0 ? 2 : (3 + precision);
/*
 * Compose the string.
 */
  sprintf(string, "%c%2.2d%c%2.2d%c%0*.*f", sgn<0 ? '-' : '+', deg,
	  colon ? ':':' ', mins, colon ? ':':' ', width, precision, secs);
  return string;
}

/*.......................................................................
 * Given the Gregorian year (eg. 1991) and a day number within that year
 * return the equivalent month and day within that year.
 * 'dayno' need not be within the range 1 -> 365. If it isn't then
 * before use, year will be incremented or decremented (as
 * appropriate) and dayno will be adjusted to the corresponding day
 * within that year. 
 *
 * Input:
 *   year    int    Year.
 *   dayno   int    Day number within 'year' - starts at day 1.
 * Output:
 *   day     int  * The day within month 'month'.
 *   month   int  * The month within which 'dayno' falls.
 */
void daydate(int year, int dayno, int *day, int *month)
{
  int i, leap;
/*
 * Constrain dayno and year to legal values.
 */
  if(dayno > 365) {
    year += (dayno-1) / 365;
    dayno = (dayno-1) % 365 + 1;
  }
  else if(dayno < 1) {
    year -= abs(dayno) / 365 + 1;
    dayno = 365 - abs(dayno) % 365;
  };
/*
 * Check whether 'year' is a leap year 'leap' will be used to
 * index into the corresponding row of the daytab array.
 */
  leap = (year%4 == 0 && year%100 != 0) || year%400 == 0;
/*
 * Find the month and day corresponding to dayno.
 */
  for(i=0; dayno > daytab[leap][i]; i++)
    dayno -= daytab[leap][i];
  *month = i+1;
  *day   = dayno;
  return;
}

/*.......................................................................
 * Given the Gregorian year (eg. 1991) and a day number within that year
 * return a string representation of the corresponding year, month and
 * day. Eg. "1991 Sep 24". The string to hold the date string must be
 * provided by the user. The pointer to the date string is also returned
 * as the function's main return value.
 * 'dayno' need not be within the range 1 -> 365. If it isn't then
 * before use, year will be incremented or decremented (as
 * appropriate) and dayno will be adjusted to the corresponding day
 * within that year. 
 *
 * Input:
 *   year    int    Year.
 *   dayno   int    Day number within 'year' - starts at day 1.
 * Output:
 *   date    char * You must supply a pointer here to a char array of at
 *                  least 12 characters.
 */
char *sdaydate(int year, int dayno, char *date)
{
  int day, month;
/*
 * Hand the conversion job to daydate().
 */
  daydate(year, dayno, &day, &month);
/*
 * Encode the date in the user's string.
 */
  if(date)
    sprintf(date, "%4.4d %3.3s %2.2d", year, month_name[month-1], day);
  else
    lprintf(stderr, "ERROR: NULL string pointer given to sdaydate()\n");
  return date;
}

/*.......................................................................
 * Given a time, as the number of seconds of UTC wrt the beginning of the
 * year, and a character array of at least 12 characters - write the
 * corresponding date in the string, like "1991 Sep 24". The string
 * will be null terminated and returned.
 *
 * Input:
 *   year      int   The year from the start of which, 'ut' seconds
 *                   were counted.
 *   vlbut  double   Number of seconds since start of 'year'.
 * Output:
 *   string   char * The coresponding date will be written like
 *                   "1991 Sep 24" and null terminated. Here the caller
 *                   must provide a character array of at least 12 elements.
 */
char *sutdate(int year, double vlbut, char *string)
{
  int dayno, hour, mins;
  double secs;
/*
 * Get the day number via dayut().
 */
  dayut(vlbut, &dayno, &hour, &mins, &secs);
/*
 * Print the string via sdaydate().
 */
  return sdaydate(year, dayno, string);
}

/*.......................................................................
 * Given a UT, as the number of seconds since the beginning of
 * the year - return the current day number and the hours,mins and seconds
 * into that day. If any of dayno, hour, mins, secs are not required send NULL.
 *
 * Input:
 *   vlbut  double    UT in seconds from start of year 'year'.
 * Output:
 *   dayno     int *  Day number within year 'year'.
 *   hour      int *  Hour within day 'dayno'.
 *   mins      int *  Minute within hour 'hour'.
 *   secs    float *  Decimal seconds within minute 'mins'.
 */
void dayut(double vlbut, int *dayno, int *hour, int *mins, double *secs)
{
  double ddays, dhour, dmins, dsecs;
/*
 * How many seconds into year?
 */
  dsecs = vlbut;
/*
 * Minutes into year.
 */
  dmins = floor(dsecs/60.0);
/*
 * Seconds into minute 'dmins'.
 */
  dsecs -= dmins * 60.0;
/*
 * Hours into year.
 */
  dhour = floor(dmins/60.0);
/*
 * Minutes into hour 'dhour'.
 */
  dmins -= dhour * 60.0;
/*
 * Days into year.
 */
  ddays = floor(dhour/24.0);
/*
 * Hours into day 'ddays'.
 */
  dhour -= ddays * 24.0;
/*
 * Assign the return values.
 */
  if(dayno) *dayno = ddays + 1.0;  /* First day is day 1 */
  if(hour) *hour = dhour;
  if(mins) *mins = dmins;
  if(secs) *secs = dsecs;
  return;
}

/*.......................................................................
 * Given a user provided char array with at least 13 elements, and a
 * ut, as the number of seconds of UTC wrt the start of the year,
 * return a string representation of the ut in the user's string.
 * The format produced is:  DDD/HH:MM:SS  where DDD is the day number.
 *
 * Input:
 *  vlbut double    Number of seconds since start of year.
 *  nc       int    The number of characters available in string[]
 *                  (including space for '\0').
 * Input/Output:
 *  string  char *  User provided string. This should have at least
 *                  13 characters. If the time can not be fitted within
 *                  nc characters (including '\0') then it will not be
 *                  be written.
 * Output:
 *  return   int    The number of characters written (excluding '\0').
 *                  This will be zero if the time could not be written
 *                  without truncation.
 */
int write_ut(double vlbut, int nc, char *string)
{
  char stmp[20];         /* Temporary buffer to compose the output in */
  int slen;              /* The number of characters written */
  int dayno, hour, mins;
  double secs;
/*
 * Let dayut() do all the hard work converting vlbut to day number
 * and UT.
 */
  dayut(vlbut, &dayno, &hour, &mins, &secs);
/*
 * Initially write the string in stmp[] to find out how long it is.
 */
  sprintf(stmp, "%d/%2.2d:%2.2d:%2.2d", dayno, hour, mins, (int) secs);
  slen = strlen(stmp);
/*
 * Copy the time to the output string if there is room.
 */
  if(slen >= nc || !string) {
    if(string && nc>0)
      *string = '\0';
    return 0;
  };
  strcpy(string, stmp);
  return slen;
}

/*.......................................................................
 * Decode a string of form DDD/HH:MM:SS into a date in seconds wrt the
 * start of the year.
 *
 * Input:
 *   s        char *  The string to be decoded.
 * Input/Output:
 *   endp     char ** If endp!=NULL then *endp will point at the next
 *                    un-processed character following the date.
 *                    If endp==NULL trailing input following the date
 *                    will be treated as an error.
 * Output:
 *   return double    The date in seconds, or < 0 on error.
 *                    On error *endp will be left equal to 's'.
 */
double read_ut(char *s, char **endp)
{
  int dayno=0;     /* Day number within year */
  int hour=0;      /* Hour number */
  int min=0;       /* Minute number */
  double sec=0.0;  /* Seconds */
  char *sptr;      /* Pointer into s[] */
/*
 * Bad string?
 */
  if(!s) {
    lprintf(stderr, "read_ut: NULL date string received.\n");
    return -1.0;
  };
/*
 * Initialize the pointer to unprocessed input.
 */
  if(endp)
    *endp = s;
/*
 * Read each component in turn.
 */
  sptr = s;
/*
 * Get the day number.
 */
  while(*sptr && isspace((int)*sptr))
    sptr++;
  if(*sptr && isdigit((int)*sptr)) {
    char *eptr;
    dayno = (int) strtol(sptr, &eptr, 10);
    sptr = (*eptr=='/') ? eptr+1 : eptr;
  } else {
    lprintf(stderr, "read_ut: Missing day number in \"%s\".\n", s);
    return -1;
  };
/*
 * Get the hour.
 */
  while(*sptr && isspace((int)*sptr))
    sptr++;
  if(*sptr && isdigit((int)*sptr)) {
    hour = strtol(sptr, &sptr, 10);
    if(hour > 23) {
      lprintf(stderr, "read_ut: Hour value out of range in \"%s\".\n", s);
      return -1.0;
    };
    if(*sptr==':')
      sptr++;
  };
/*
 * Get the minute.
 */
  while(*sptr && isspace((int)*sptr))
    sptr++;
  if(*sptr && isdigit((int)*sptr)) {
    min = strtol(sptr, &sptr, 10);
    if(min > 60) {
      lprintf(stderr, "read_ut: Minute value out of range in \"%s\".\n", s);
      return -1.0;
    };
    if(*sptr==':')
      sptr++;
  };
/*
 * Get the second.
 */
  while(*sptr && isspace((int)*sptr))
    sptr++;
  if(*sptr && isdigit((int)*sptr))
    sec = strtod(sptr, &sptr);
  if(sec > 60.0) {
    lprintf(stderr, "read_ut: Seconds value out of range in \"%s\".\n", s);
    return -1.0;
  };
/*
 * Find the end of the string.
 */
  while(*sptr && isspace((int)*sptr))
    sptr++;
  if(endp) {
    *endp = sptr;
  } else if(*sptr) {
    lprintf(stderr, "Bad date string: %s\n", s);
    return -1.0;
  };
  return sec + 60.0 * (min + 60.0 * (hour + 24.0 * (dayno-1)));
}

/*.......................................................................
 * Convert a ut, measured as the number of seconds since the beginning of
 * the year 'year', into the corresponding julian day number and fractional
 * day from mean noon on that day.
 *
 * Input:
 *   vlbut  double   Number of seconds since start of year.
 *   year      int   The year from the start of which vlbut was formed.
 * Output:
 *   jd       long * The day corresponding julian day number.
 *   jdfrc  double * The fractional day since noon of julian day 'jd'.
 *   je     double * The julian day expressed as a Julian epoch.
 */
void julday(double vlbut, int year, long *jd, double *jdfrc, double *je)
{
  int dayno, hour, mins, iyear, icent;
  long jd_jan0;
  double secs;
/*
 * First decompose the UT into day number within the year and
 * hours, minutes and seconds.
 */
  dayut(vlbut, &dayno, &hour, &mins, &secs);
/*
 * Calculate the julian day number at noon on Jan 0 of 'year'.
 */
  iyear = year - 1;
  icent = iyear / 100;
  jd_jan0 = 1721425L + 365*iyear + iyear/4 - icent + icent/4;
/*
 * Work out the fraction of that day by which the ut follows
 * noon of the previous day. (Julian days start at noon).
 */
  *jdfrc = 0.5 + (hour + (mins + secs / 60.0) / 60.0 ) / 24.0;
/*
 * Turn dayno and frcday into the number of julian days since noon Jan
 * 0 and the fractional day following that julian day.
 */
  if(*jdfrc < 1.0)
    dayno -= 1;
  else
    *jdfrc -= 1.0;
/*
 * Work out the complete (integer) julian day number.
 */
  *jd = jd_jan0 + dayno;
/*
 * Turn this into the equivalent Julian Epoch.
 */
  *je = 2000.0 + (*jd-2451545L)/365.25;
  return;
}

/*.......................................................................
 * Return a string containing the current date and time.
 * This function returns a copy of the string returned by the standard
 * library function asctime(), devoid of the annoying '\n'.
 *
 * Output:
 *  return  char *  A '\0' terminated string containing the date and time.
 *                  On NULL return NULL.
 */
char *date_str(void)
{
  static char str[81]; /* The return string */
  time_t tp;
  char *nptr;  /* Pointer to newline character in return string */
/*
 * Get the time.
 */
  if(time(&tp) == -1) {
    lprintf(stderr, "date_str: No date available on this machine\n");
    return NULL;
  };
/*
 * Get a formatted string version of this and copy it into str[].
 */
  strncpy(str, ctime(&tp), sizeof(str)-1);
/*
 * Ensure that it is '\0' terminated.
 */
  str[sizeof(str)-1] = '\0';
/*
 * If there is a '\n' character (as specified in the standard),
 * terminate the string there.
 */
  nptr = strchr(str, '\n');
  if(nptr) *nptr = '\0';
  return &str[0];
}

/*.......................................................................
 * Read a sexagesimal-format number from a string.
 *
 * Input:
 *  string   char *  The string to be parsed.
 * Input/Output:
 *  value  double *  The number, in the units of the most significant
 *                   component.
 *  endp     char ** The pointer to the next unprocessed character of
 *                   string[] will be assigned to *endp. On failure
 *                   this *endp==string, indicating that no characters
 *                   were successfully processed. If endp==NULL, complain
 *                   if any characters follow the date.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
int parse_sexagesimal_string(char *string, double *value, char **endp)
{
  Number n;             /* A number read by input_number() */
  double number;        /* The number accumulated so far */
  double divisor = 1.0; /* The scale factor to divide the next component by */
  int negative = 0;     /* True if the number is negative */
  char *sptr;           /* A pointer into string[] */
/*
 * Check the arguments.
 */
  if(!string || !value) {
    lprintf(stderr, "parse_sexagesimal_string: NULL argument(s).\n");
    if(endp)
      *endp = string;
    return 1;
  };
/*
 * Arrange that if we return early due to an error, *endp will be set
 * to string.
 */
  if(endp)
    *endp = string;
/*
 * Skip leading spaces.
 */
  for(sptr=string; isspace((int) *sptr); sptr++)
    ;
/*
 * Is the number negative?
 */
  negative = *sptr == '-';
/*
 * Get the absolute value of the first number.
 */
  if(parse_numeric_string(sptr, &sptr, &n))
    return 1;
  number = fabs(n.type==NUM_DOUBLE ? n.value.dval : (double)(n.value.ival));
/*
 * Record where we have successfully read to so far.
 */
  if(endp)
    *endp = sptr;
/*
 * Read components until a non-integral number or a component that
 * is not followed by a component separator is encountered.
 */
  while(n.type==NUM_INT && *sptr == ':') {
    double d;
/*
 * Skip the separator and read the next component.
 */
    sptr++;
    if(parse_numeric_string(sptr, &sptr, &n))
      break;
/*
 * Get a double precision version of the number and domain check it.
 */
    d = n.type==NUM_DOUBLE ? n.value.dval : n.value.ival;
    if(d >= 60 || d < 0.0)
      break;
/*
 * Add the new component with the appropriate scale factor.
 */
    divisor *= 60.0;
    number += d / divisor;
/*
 * Record where we have successfully read to so far.
 */
    if(endp)
      *endp = sptr;
  };
/*
 * If requested, reject trailing input.
 */
  if(!endp && *sptr) {
    lprintf(stderr, "Unexpected characters follow a sexagesimal number.\n");
    return 1;
  };
/*
 * Assign the result for return.
 */
  *value = negative ? -number : number;
  return 0;
}

/*.......................................................................
 * Read a number from a string. Note that this function does not skip
 * leading white-space.
 *
 * Input:
 *  string        char *  The string containing the number.
 *  endp          char *  The pointer to the next unprocessed character of
 *                        string[] will be assigned to *endp. On failure
 *                        this *endp==string, indicating that no characters
 *                        were successfully processed. If endp==NULL,
 *                        complain if any characters follow the number.
 * Input/Output:
 *  number      Number *  The resulting number in a form appropriate to
 *                        way that it was written.
 * Output:
 *  return         int    0 - OK.
 *                        1 - Error.
 */
int parse_numeric_string(char *string, char **endp, Number *number)
{
  int have_mantissa = 0; /* Used to record existence of mantissa digits */
  int is_int = 1;        /* True if number is an integer with no exponent */
  char *sptr;            /* A pointer into string[] */
  char *eptr;            /* A pointer to the next character to be processed */
/*
 * Arrange that if we return early due to an error, *endp will be set
 * to string.
 */
  if(endp)
    *endp = string;
/*
 * The first stage of the process is to identify whether the number
 * is specified as an integer or as a floating point number. Having
 * done this, we can read it with the appropriate C library function.
 *
 * An initial sign character is ok.
 */
  sptr = string;
  if(*sptr == '+' || *sptr == '-')
    sptr++;
/*
 * Check for and skip the integral part of the mantissa.
 */
  if(isdigit((int)*sptr)) {
    have_mantissa = 1;
    while(isdigit((int)*sptr))
      sptr++;
  };
/*
 * Is there a decimal point?
 */
  if(*sptr == '.') {
/*
 * OK, we now know that the number isn't integral.
 */
    is_int = 0;
    sptr++;
/*
 * See if there are any digits following the decimal point.
 */
    if(isdigit((int)*sptr)) {
      have_mantissa = 1;
      while(isdigit((int)*sptr))
	sptr++;
    };
  };
/*
 * If no mantissa was found, we don't have a valid number.
 */
  if(!have_mantissa)
    return 1;
/*
 * Check for an exponent.
 */
  if(*sptr == 'e' || *sptr == 'E') {
    sptr++;
/*
 * Skip the optional sign.
 */
    if(*sptr == '+' || *sptr == '-')
      sptr++;
/*
 * If the exponent is valid, then the number is in floating point
 * format.
 */
    if(isdigit((int)*sptr))
      is_int = 0;
  };
/*
 * Now read the number, according to its format.
 */
  if(is_int) {
    number->type = NUM_INT;
    number->value.ival = strtol(string, &eptr, 10);
  } else {
    number->type = NUM_DOUBLE;
    number->value.dval = strtod(string, &eptr);
  };
/*
 * If possible return a pointer to the next character to be processed.
 * Otherwise complain if any characters remain.
 */
  if(endp) {
    *endp = eptr;
  } else if(*eptr) {
    lprintf(stderr, "Unexpected characters follow valid number.\n");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Parse a date and time from a string.
 *
 * The following are valid specifications:
 *
 *  12-DEC-1998                 ->  12-DEC-1998 00:00:00.0
 *  12-DEC-1998:13              ->  12-DEC-1998 13:00:00.0
 *  12-DEC-1998:13:34           ->  12-DEC-1998 13:34:00.0
 *  12-DEC-1998:13:34:23        ->  12-DEC-1998 13:34:23.0
 *  12-DEC-1998:13:34:23.5      ->  12-DEC-1998 13:34:23.5
 *
 * The following are also valid if nospace!=0.
 *
 *  12-DEC-1998 13              ->  12-DEC-1998 13:00:00.0
 *  12-DEC-1998 13:34           ->  12-DEC-1998 13:34:00.0
 *  12-DEC-1998 13:34:23        ->  12-DEC-1998 13:34:23.0
 *  12-DEC-1998 13:34:23.5      ->  12-DEC-1998 13:34:23.5
 *
 * Input:
 *  string    const char *  The string containing the number.
 *  tell             int    Non-zero to report errors to stderr.
 *  nospace          int    By default, the year field can be separated
 *                          from the hour field by either a colon or a space.
 *                          To disallow a space as a separator, set nospace
 *                          to non-zero.
 * Input/Output:
 *  endp      const char ** If endp!=NULL, the pointer to the next
 *                          unprocessed character of string[] will be
 *                          assigned to *endp. On failure this *endp==string,
 *                          indicating that no characters were successfully
 *                          processed. If endp==NULL, complain if any
 *                          characters follow the date and time.
 *  int             year *  The year will be assigned to *year.
 *  int            month *  The month number within the year (1,12) will
 *                          be assigned to *month.
 *  int              day *  The day number of the month (1,31) will be
 *                          assigned to *day.
 *  int             hour *  The number of hours into day 'day' (0..23).
 *  int              min *  The number of minutes into hour 'hour' (0..59).
 *  double           sec *  The number of seconds into minute 'min' (0..<60).
 * Output:
 *  return           int    0 - OK.
 *                          1 - Error.
 */
int parse_date_and_time(const char *string, int tell, int nospace,
			const char **endp, int *year, int *month, int *day,
			int *hour, int *min, double *sec)
{
  int read_time;     /* True if there is a time-specification to be read */
  const char *cptr;  /* A pointer to the next character to be processed */
/*
 * On error, endp is set to &string[0] to indicate that nothing was
 * succesfully read. Be pesimistic and set this up now, so that we
 * don't have to make this assignment for each error condition.
 */
  if(endp)
    *endp = string;
/*
 * Check the arguments.
 */
  if(!string || !year || !month || !day || !hour || !min || !sec) {
    lprintf(stderr, "input_date_and_time: NULL argument(s).\n");
    return 1;
  };
/*
 * Read the date components.
 */
  if(parse_date(string, tell, &cptr, year, month, day))
    return 1;
/*
 * See if there is a time specification.
 */
  if(*cptr == ':' || (!nospace && *cptr==' ')) {
    if(*(++cptr) == '\0') {
      if(tell)
	lprintf(stderr, "Missing time specification.\n");
      return 1;
    };
    read_time = isdigit((int)*cptr);
  } else {
    read_time = 0;
  };
/*
 * Read a time specification?
 */
  if(read_time) {
    if(parse_time(cptr, tell, &cptr, hour, min, sec))
      return 1;
  } else {
    *hour = *min = 0;
    *sec = 0.0;
  };
/*
 * If possible return a pointer to the next character to be processed.
 * Otherwise complain if any characters remain.
 */
  if(endp) {
    *endp = cptr;
  } else if(*cptr) {
    lprintf(stderr, "Unexpected characters follow valid date and time.\n");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Parse a sexagesimal time (eg. 23:34:00) from a string.
 * The recognized variations, and their meanings, are:
 *
 *  23:34:00.0  ->   23:34:00
 *  23:34:00    ->   23:34:00
 *  23:34       ->   23:34:00
 *  23          ->   23:00:00
 *
 * Input:
 *  string    const char *  The string to read from.
 *  tell             int    Non-zero to report errors to stderr.
 * Input/Output:
 *  endp      const char ** If endp!=NULL, the pointer to the next
 *                          unprocessed character of string[] will be
 *                          assigned to *endp. On failure this *endp==string,
 *                          indicating that no characters were successfully
 *                          processed. If endp==NULL, complain if any
 *                          characters follow the time.
 *  int             hour *  The number of hours (0..23).
 *  int              min *  The number of minutes into hour 'hour' (0..59).
 *  double           sec *  The number of seconds into minute 'min' (0..<60).
 * Output:
 *  return           int    0 - OK.
 *                          1 - Error.
 */
int parse_time(const char *string, int tell, const char **endp,
	       int *hour, int *min, double *sec)
{
  char *usage = "Invalid time - use hh:mm:ss.s";
  unsigned long v;    /* An hour or minute value */
  double s;           /* A seconds value */
  const char *cptr;   /* A pointer to the next character to be processed */
/*
 * On error, endp is set to &string[0] to indicate that nothing was
 * succesfully read. Be pesimistic and set this up now, so that we
 * don't have to make this assignment for each error condition.
 */
  if(endp)
    *endp = string;
/*
 * Check the arguments.
 */
  if(!string || !hour || !min || !sec) {
    lprintf(stderr, "parse_time: NULL argument(s).\n");
    return 1;
  };
/*
 * Read the hour component.
 */
  if(parse_ulong(string, 0, &cptr, &v) || v > 23) {
    if(tell)
      lputs(usage, stderr);
    return 1;
  };
  *hour = v;
/*
 * If there is a minutes component then there should be a : separator.
 */
  if(*cptr == ':') {
/*
 * Skip the separator and read the minutes component.
 */
    if(*++cptr=='\0' || parse_ulong(cptr, 0, &cptr, &v) || v > 59) {
      if(tell)
	lputs(usage, stderr);
      return 1;
    };
    *min = v;
  } else {
    *min = 0;
  };
/*
 * If there is a seconds component then there should be a : separator.
 */
  if(*cptr == ':') {
/*
 * Skip the separator and read the seconds component
 */
    if(*++cptr == '\0' || parse_double(cptr, 0, &cptr, &s) ||
       s < 0.0 || s >= 60.0) {
      if(tell)
	lputs(usage, stderr);
      return 1;
    };
    *sec = s;
  } else {
    *sec = 0.0;
  };
/*
 * If possible return a pointer to the next character to be processed.
 * Otherwise complain if any characters remain.
 */
  if(endp) {
    *endp = cptr;
  } else if(*cptr) {
    lprintf(stderr, "Unexpected characters follow valid time.\n");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Read a date of form 23JAN1997 or 23/01/1997 from a string.
 * It is assumed that the entered year includes its century.
 *
 * Input:
 *  string    const char *  The input string.
 *  tell             int    Non-zero to report errors to stderr.
 * Input/Output:
 *  endp      const char ** If endp!=NULL, the pointer to the next
 *                          unprocessed character of string[] will be
 *                          assigned to *endp. On failure this *endp==string,
 *                          indicating that no characters were successfully
 *                          processed. If endp==NULL, complain if any
 *                          characters follow the date.
 *  int             year *  The year will be assigned to *year.
 *  int            month *  The month number within the year (1,12) will
 *                          be assigned to *month.
 *  int              day *  The day number of the month (1,31) will be
 *                          assigned to *day.
 * Output:
 *  return           int    0 - OK.
 *                          1 - Error.
 */
int parse_date(const char *string, int tell, const char **endp,
	       int *year, int *month, int *day)
{
  char *usage = "Invalid date - use DD-MMM-YYYY\n";
  unsigned long dy;       /* The day of the month */
  unsigned long mn;       /* The month of the year */
  unsigned long yr;       /* The year */
  int isleap;             /* True if the parsed year is a leap year */
  enum {MONTH_LEN=3};     /* The length of a month-name abbreviation */
  char mname[MONTH_LEN+1];/* A 3-letter month-name abbreviation + '\0' */
  enum {NUM_MONTH=12};    /* The number of months in a year */
  const char *cptr;       /* A pointer to the next character to be processed */
  int i;
/*
 * Record the number of days per month, first in normal years, then in
 * leap years.
 */
  static char daytab[2][NUM_MONTH] = {
    {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
  };
/*
 * List the conventional 3-letter abbreviations for the names of months.
 */
  static char *months[NUM_MONTH] = {
    "JAN",  "FEB",  "MAR",  "APR", 
    "MAY",  "JUN",  "JUL",  "AUG",
    "SEP",  "OCT",  "NOV",  "DEC"
  };
/*
 * On error, endp is set to &string[0] to indicate that nothing was
 * succesfully read. Be pesimistic and set this up now, so that we
 * don't have to make this assignment for each error condition.
 */
  if(endp)
    *endp = string;
/*
 * Check arguments.
 */
  if(!string || !day || !month || !year) {
    lprintf(stderr, "parse_date: NULL argument(s).\n");
    return 1;
  };
/*
 * Read the day of the month.
 */
  if(parse_ulong(string, 0, &cptr, &dy)) {
    if(tell)
      lputs(usage, stderr);
    return 1;
  };
/*
 * The next character should be a '-' separator.
 */
  if(*cptr != '-') {
    if(tell)
      lputs(usage, stderr);
    return 1;
  };
/*
 * Skip the separator.
 */
  if(*++cptr == '\0') {
    if(tell)
      lputs(usage, stderr);
    return 1;
  };
/*
 * Read the 3 characters of a month name abbreviation and
 * convert them to upper case.
 */
  for(i=0; i<MONTH_LEN; i++) {
/*
 * Check the new character.
 */
    if(!isalpha((int)*cptr)) {
      if(tell)
	lputs(usage, stderr);
      return 1;
    };
/*
 * Convert the new character to lower case and add it to the
 * accumulated month name.
 */
    mname[i] = islower((int)*cptr) ? toupper((int)*cptr) : *cptr;
/*
 * Advance to the next character.
 */
    if(*++cptr == '\0') {
      if(tell)
	lputs(usage, stderr);
      return 1;
    };
  };
/*
 * Terminate the month name.
 */
  mname[MONTH_LEN] = '\0';
/*
 * Lookup the month name.
 */
  for(mn=1; mn<=NUM_MONTH && strcmp(mname, months[mn-1]) != 0; mn++)
    ;
/*
 * Month not recognized?
 */
  if(mn > NUM_MONTH) {
    if(tell)
      lprintf(stderr, "Unknown month [%s]", mname);
    return 1;
  };
/*
 * The next character should be a '-' separator.
 */
  if(*cptr != '-') {
    if(tell)
      lputs(usage, stderr);
    return 1;
  };
/*
 * Skip the separator.
 */
  if(*++cptr == '\0') {
    if(tell)
      lputs(usage, stderr);
    return 1;
  };
/*
 * Read the year.
 */
  if(parse_ulong(cptr, 0, &cptr, &yr)) {
    if(tell)
      lputs(usage, stderr);
    return 1;
  };
/*
 * Is this a leap year?
 */
  isleap = (yr%4 == 0 && yr%100 != 0) || yr%400 == 0;
/*
 * Check that the date makes sense.
 */
  if(dy < 1 || dy > daytab[isleap][mn-1]) {
    if(tell)
      lprintf(stderr, "Nonexistent date (%02lu-%.3s-%04lu)", dy,
	      months[mn-1], yr);
    return 1;
  };
/*
 * Assign the values for return.
 */
  *day = dy;
  *month = mn;
  *year = yr;
/*
 * Return a pointer to the next unprocessed character in the
 * string. If endp==NULL, treat the existence of any characters
 * following the date, as an error.
 */
  if(endp) {
    *endp = cptr;
  } else if(*cptr) {
    lprintf(stderr, "Unexpected characters follow a valid date.\n");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Read an unsigned long decimal integer from a string.
 *
 * Input:
 *  string    const char *  The string to read the integer from.
 *  tell             int    Non-zero to report errors to stderr.
 *  endp      const char ** If endp!=NULL, the pointer to the next
 *                          unprocessed character of string[] will be
 *                          assigned to *endp. On failure this *endp==string,
 *                          indicating that no characters were successfully
 *                          processed. If endp==NULL, complain if any
 *                          characters follow the number.
 * Input/Output:
 *  ulval  unsigned long *  On output *ulval will contain the unsigned
 *                          long integer read from the string.
 * Output:
 *  return           int    0 - OK.
 *                          1 - Error.
 */
int parse_ulong(const char *string, int tell, const char **endp,
		unsigned long *ulval)
{
  const char *cptr;  /* A pointer into the input string */
/*
 * On error, endp is set to &string[0] to indicate that nothing was
 * succesfully read. Be pesimistic and set this up now, so that we
 * don't have to make this assignment for each error condition.
 */
  if(endp)
    *endp = string;
/*
 * Check the arguments.
 */
  if(!string || !ulval) {
    lprintf(stderr, "parse_ulong: NULL argument(s).\n");
    return 1;
  };
/*
 * Attempt to read the number.
 */
  *ulval = strtoul(string, (char **)&cptr, 10);
/*
 * If strtoul() didn't parse any of the input string, then there
 * wasn't a valid integer there.
 */
  if(cptr == string) {
    if(tell)
      lprintf(stderr, "Missing unsigned integer.\n");
    return 1;
  };
/*
 * If possible return a pointer to the next character to be processed.
 * Otherwise complain if any characters remain.
 */
  if(endp) {
    *endp = cptr;
  } else if(*cptr) {
    lprintf(stderr, "Unexpected characters follow a valid number.\n");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Read a double precision number from a string.
 *
 * Input:
 *  string    const char *  The string to read the number from.
 *  tell             int    Non-zero to report errors to stderr.
 *  endp      const char ** If endp!=NULL, the pointer to the next
 *                          unprocessed character of string[] will be
 *                          assigned to *endp. On failure this *endp==string,
 *                          indicating that no characters were successfully
 *                          processed. If endp==NULL, complain if any
 *                          characters follow the number.
 * Input/Output:
 *  dval          double *  On output *dval will contain the floating
 *                          number read from the string.
 * Output:
 *  return           int    0 - OK.
 *                          1 - Error.
 */
int parse_double(const char *string, int tell, const char **endp,
		 double *dval)
{
  const char *cptr;  /* A pointer into the input string */
/*
 * On error, endp is set to &string[0] to indicate that nothing was
 * succesfully read. Be pesimistic and set this up now, so that we
 * don't have to make this assignment for each error condition.
 */
  if(endp)
    *endp = string;
/*
 * Check the arguments.
 */
  if(!string || !dval) {
    lprintf(stderr, "parse_double: NULL argument(s).\n");
    return 1;
  };
/*
 * Attempt to read the number.
 */
  *dval = strtod(string, (char **)&cptr);
/*
 * If strtod() didn't parse any of the input string, then there
 * wasn't a valid number there.
 */
  if(cptr == string) {
    if(tell)
      lprintf(stderr, "Missing number.\n");
    return 1;
  };
/*
 * If possible return a pointer to the next character to be processed.
 * Otherwise complain if any characters remain.
 */
  if(endp) {
    *endp = cptr;
  } else if(*cptr) {
    lprintf(stderr, "Unexpected characters follow a valid number.\n");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Given a string containing a date and time like dd-mmm-yyyy:hh:mm:ss.ss,
 * return the corresponding Modified Julian Date.
 *
 * Input:
 *  string  const char *  The string containing the Gregorian date to be parsed.
 *  tell           int    Non-zero to report errors to stderr.
 * Input/Output:
 *  endp    const char ** If endp!=NULL, the pointer to the next
 *                        unprocessed character of string[] will be
 *                        assigned to *endp. On failure this *endp==string,
 *                        indicating that no characters were successfully
 *                        processed. If endp==NULL, any characters following
 *                        the MJD will be regarded as an error, and reported
 *                        as such.
 *  mjd         double *  The Modified Julian Date will be assigned to *mjd.
 * Output:
 *  return         int    0 - OK.
 *                        1 - Error.
 */
int parse_mjd(const char *string, int tell, const char **endp, double *mjd)
{
  int year, month, day;   /* The Gregorian date */
  int hour, min;          /* The hours and minutes into the day */
  double sec;             /* The number of seconds into the minute */
  int status;             /* The return status of a slalib function */
/*
 * On error, endp is set to &string[0] to indicate that nothing was
 * succesfully read. Be pesimistic and set this up now, so that we
 * don't have to make this assignment for each error condition.
 */
  if(endp)
    *endp = string;
/*
 * Check the arguments.
 */
  if(!string || !mjd) {
    lprintf(stderr, "parse_mjd: NULL argument(s).\n");
    return 1;
  };
/*
 * Parse the string.
 */
  if(parse_date_and_time(string, tell, 1, endp, &year, &month, &day, &hour,
			 &min, &sec))
    return 1;
/*
 * Check for extra characters following the MJD.
 */
  
/*
 * Convert the Gregorian calendar date into a Modified Julian Date.
 */
  slaCldj(year, month, day, mjd, &status);
  if(status) {
    if(tell)
      lprintf(stderr, "Invalid date.\n");
    return 1;
  };
/*
 * Add the time of day to the MJD.
 */
  *mjd += (hour + (min + sec/60.0)/60.0)/24.0;
  return 0;
}


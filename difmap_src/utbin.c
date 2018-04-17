#include <math.h>
#include "utbin.h"

/*.......................................................................
 * Return a description of the time-limits of the time-bin within which
 * a given time falls.
 *
 * Input:
 *  ut        double   The time-stamp for which the date of the associated
 *                     bin is required (seconds).
 *  origin    double   The time at which the first bin in the time-grid starts
 *                     (seconds). This should normally be UT=0 on the first
 *                     day of the observation ie.
 *                       origin = ob->date.ut;
 *  binwid    double   The width of a single bin. If binwid<1 no binning
 *                     will be applied and the start, middle and end time
 *                     of the bin will all be set to the value of 'ut'.
 * Output:
 *  return     UTbin * A pointer to a static internal container that contains
 *                     the start, central and end time limits of the time
 *                     bin that 'ut' falls within.
 */
UTbin *bintime(double origin, double ut, double binwid)
{
  static UTbin utbin;  /* The return container */
  if(binwid>=1.0) {
    utbin.beg_ut = origin + binwid * floor((ut-origin)/binwid);
    utbin.mid_ut = utbin.beg_ut + binwid/2.0;
    utbin.end_ut = utbin.beg_ut + binwid;
  } else {
    utbin.end_ut = utbin.mid_ut = utbin.beg_ut = ut;
  };
  return &utbin;
}

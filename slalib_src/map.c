#include "slalib.h"
#include "slamac.h"
void slaMap ( double rm, double dm, double pr, double pd,
              double px, double rv, double eq, double date,
              double *ra, double *da )
/*
**  - - - - - - -
**   s l a M a p
**  - - - - - - -
**
**  Transform star RA,Dec from mean place to geocentric apparent.
**
**  The reference frames and timescales used are post IAU 1976.
**
**  References:
**     1984 Astronomical Almanac, pp B39-B41.
**     (also Lederle & Schwan, Astron. Astrophys. 134, 1-6, 1984)
**
**  Given:
**     rm,dm    double     mean RA,Dec (rad)
**     pr,pd    double     proper motions:  RA,Dec changes per Julian year
**     px       double     parallax (arcsec)
**     rv       double     radial velocity (km/sec, +ve if receding)
**     eq       double     epoch and equinox of star data (Julian)
**     date     double     TDB for apparent place (JD-2400000.5)
**
**  Returned:
**     *ra,*da  double     apparent RA,Dec (rad)
**
**  Called:
**     slaMappa       star-independent parameters
**     slaMapqk       quick mean to apparent
**
**  Notes:
**
**     1)  eq is the Julian epoch specifying both the reference
**         frame and the epoch of the position - usually 2000.
**         For positions where the epoch and equinox are
**         different, use the routine slaPm to apply proper
**         motion corrections before using this routine.
**
**     2)  The distinction between the required TDB and TDT is
**         always negligible.  Moreover, for all but the most
**         critical applications UTC is adequate.
**
**     3)  The proper motions in RA are dRA/dt rather than
**         cos(dec)*dra/dt.
**
**     4)  This routine may be wasteful for some applications
**         because it recomputes the Earth position/velocity and
**         the precession/nutation matrix each time, and because
**         it allows for parallax and proper motion.  Where
**         multiple transformations are to be carried out for one
**         epoch, a faster method is to call the slaMappa routine
**         once and then either the slaMapqk routine (which includes
**         parallax and proper motion) or slaMapqkz (which assumes
**         zero parallax and proper motion).
**
**  Last revision:   12 June 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double amprms[21];

/* Star-independent parameters */
   slaMappa ( eq, date, amprms );

/* Mean to apparent */
   slaMapqk ( rm, dm, pr, pd, px, rv, amprms, ra, da );
}

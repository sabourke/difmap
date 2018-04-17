#include <stdlib.h>
#include <stdio.h>

#include "logio.h"
#include "obs.h"
#include "units.h"
#include "vlbconst.h"
#include "vlbutil.h"
#include "scans.h"

/*.......................................................................
 * List useful parts of the observation header parameters from an
 * Observation class.
 *
 * Input:
 *  ob   Observation *  The observation who's header is to be listed.
 */
void vlbhead(Observation *ob)
{
  Subarray *sub;   /* The descriptor of the sub-array being described */
  If *ifp;         /* The descriptor of the IF being described. */
  Obhead *misc;    /* Pointer to ob->misc. */
  int isub;        /* The index of the sub-array being described */
  int norb;        /* The number of orbital telescopes found */
  long jd;         /* A julian day number */
  long scansum=0L; /* Sum of sub-array scan intervals */
  double jdfrc;    /* The fractional day since noon of julian day 'jd' */
  double je;       /* The julian day expressed as a Julian epoch */
  char ctmpa[20];  /* Temporary work string */
  char ctmpb[20];  /* Temporary work string */
  int i,j;
/*
 * Check that a valid observation was sent.
 */
  if(!ob_ready(ob, OB_INDEX, "vlbhead"))
    return;
/*
 * List AIPS header line keyword values.
 */
  misc = &ob->misc;
  lprintf(stdout, "\nUV FITS miscellaneous header keyword values:\n");
  lprintf(stdout, "  OBSERVER = \"%s\"\n",
	  misc->observer ? misc->observer : "(N/A)");
  lprintf(stdout, "  DATE-OBS = \"%s\"\n",
	  misc->date_obs ? misc->date_obs : "(N/A)");
  lprintf(stdout, "  ORIGIN   = \"%s\"\n",
	  misc->origin ? misc->origin  :  "(N/A)");
  lprintf(stdout, "  TELESCOP = \"%s\"\n",
	  misc->telescop ? misc->telescop : "(N/A)");
  lprintf(stdout, "  INSTRUME = \"%s\"\n",
	  misc->instrume ? misc->instrume : "(N/A)");
  lprintf(stdout, "  EQUINOX  = %.2f\n", misc->equinox);
/*
 * Describe each sub-array in the observation.
 */
  sub = ob->sub;
  for(isub=0; isub<ob->nsub; isub++,sub++) {
    lprintf(stdout, "\nSub-array %d contains:\n", isub+1);
    lprintf(stdout, " %3d baselines   %2d stations\n",sub->nbase,sub->nstat);
    lprintf(stdout, " %3d integrations   %2d scans\n", sub->ntime,
	    nscans(sub,sub->scangap));
/*
 * List the ground-based stations of the current sub-array and
 * count the number of orbital telescopes.
 */
    lprintf(stdout, "\n  Station  name               X (m)            Y (m)             Z(m)\n");
    norb = 0;
    for(i=0; i<sub->nstat; i++) {
      Station *tel = sub->tel + i;
      if(tel->type == GROUND) {
	lprintf(stdout, "    %2.2d     %-10.10s  %15e  %15e  %15e\n",
		tel->antno, tel->name, tel->geo.gnd.x, tel->geo.gnd.y,
		tel->geo.gnd.z);
      } else {
	norb++;   /* Count orbital telescopes */
      };
    };
/*
 * List any orbital telescopes that were encountered.
 */
    if(norb) {
      Binan *binan = sub->binan;
      lprintf(stdout, "\n  Station  satellite  (semi_maj eccent inclin ra_ascend arg_perig anomoly)\n");
      for(i=0; i < sub->nstat; i++) {
	Station *tel = sub->tel + i;
	if(tel->type == ORBIT) {
	  Bintel *btel = binan->bt + i;
	  lprintf(stdout, "    %2.2d     %-10.10s", tel->antno, tel->name);
/*
 * List binan->numorb orbital parameters, 3 per line.
 */
	  for(j=0; j<binan->numorb; j++) {
	    if(j>0 && j%3==0)
	      lprintf(stdout, "\n                     ");
	    lprintf(stdout, "  %15e", btel->orbparm[j]);
	  };
	  lprintf(stdout, "\n");
	};
      };
    };
/*
 * Accumulate the sum of scan durations.
 */
    scansum += timescans(sub, sub->scangap);
  };
/*
 * List the chracteristics of each IF.
 */
  lprintf(stdout, "\nThere %s %d IF%s, and a total of %d channel%s:\n",
	  ob->nif>1?"are":"is", ob->nif, ob->nif>1?"s":"",
	  ob->nctotal, ob->nctotal>1?"s":"");
  lprintf(stdout, "\n  %s\n  %s\n  %s\n",
    "IF  Channel    Frequency  Freq offset  Number of   Overall IF",
    "     origin    at origin  per channel   channels    bandwidth",
    "------------------------------------------------------------- (Hz)");
  for(i=0; i<ob->nif; i++) {
    ifp = &ob->ifs[i];
    lprintf(stdout, "  %2.2d  %7d %12g %12g    %7d %12g\n", i+1, ifp->coff+1,
	    ifp->freq, ifp->df, ob->nchan, ifp->bw);
  };
/*
 * Source parameters.
 */
  lprintf(stdout, "\nSource parameters:\n");
  lprintf(stdout, "  Source: \t %s\n", ob->source.name);
  lprintf(stdout, "  RA     = \t %s (%.1f)\t %s (apparent)\n",
	  sradhms(ob->source.ra, 3, 0, ctmpa), ob->source.epoch,
	  sradhms(ob->source.app_ra, 3, 0, ctmpb));
  lprintf(stdout, "  DEC    = \t%s        \t%s\n",
	  sraddms(ob->source.dec, 3, 0, ctmpa),
	  sraddms(ob->source.app_dec, 3, 0, ctmpb));
/*
 * If the observing center was provided in the FITS header, list it.
 */
  if(ob->source.have_obs) {
    lprintf(stdout, "\nAntenna pointing center:\n");
    lprintf(stdout, "  OBSRA  = \t %s (%.1f)\n",
	    sradhms(ob->source.obsra, 3, 0, ctmpa), ob->source.epoch);
    lprintf(stdout, "  OBSDEC = \t%s\n",
	    sraddms(ob->source.obsdec, 3, 0, ctmpa));
  };
/*
 * Describe the data.
 */
  lprintf(stdout, "\nData characteristics:\n");
  lprintf(stdout, "  Recorded units are %s.\n", misc->bunit?misc->bunit:"Jy");
  lprintf(stdout, "  Recorded polarizations:");
/*
 * List the polarizations on the same as above line.
 */
  for(i=0; i<ob->npol; i++)
    lprintf(stdout, " %s", Stokes_name(ob->pols[i]));
  lprintf(stdout, "\n");
  lprintf(stdout, "  Phases are rotated %g %s East and %g %s North.\n",
	  radtoxy(ob->geom.east), mapunits(U_NAME),
	  radtoxy(ob->geom.north), mapunits(U_NAME));
  lprintf(stdout, "  UVW coordinates are rotated by %g degrees clockwise.\n",
	  ob->geom.uvangle * rtod);
  lprintf(stdout, "  Scale factor applied to FITS data weights: %g\n",
	  ob->geom.wtscale);
  lprintf(stdout, "  Coordinate projection: %s\n", Proj_name(ob->proj));
/*
 * Decribe the overall dimensions of the data set.
 */
  lprintf(stdout, "\nSummary of overall dimensions:\n");
  lprintf(stdout,
       "  %d sub-arrays, %d IFs, %d channels, %d integrations\n",
	  ob->nsub, ob->nif, ob->nctotal, ob->nrec);
  lprintf(stdout, "  %d polarizations, and up to %d baselines per sub-array\n",
	  ob->npol, ob->nbmax);
/*
 * List time-related parameters.
 */
  lprintf(stdout, "\nTime related parameters:\n");
/*
 * Report the reference ut as yeay day-number/hours:mins:secs
 * and like (1991 Sep 24).
 */
  write_ut(ob->date.ut, sizeof(ctmpa), ctmpa);
  lprintf(stdout, "  Reference date: %d day %s  (%s)\n",
	 ob->date.year, ctmpa, sutdate(ob->date.year, ob->date.ut, ctmpb));
/*
 * Also report reference date as a Julian day number and epoch.
 */
  julday(ob->date.ut, ob->date.year, &jd, &jdfrc, &je);
  lprintf(stdout, "  Julian Date: %ld.%2.2d, Epoch J%.3f\n",
	  jd, (int) (jdfrc*100.0), je);
/*
 * Report the GAST at the above reference time.
 */
  lprintf(stdout,"  GAST at reference date: %s\n",
	  sradhms(ob->date.app_st, 3, 0, ctmpa));
/*
 * Miscellaneous information.
 */
  lprintf(stdout, "  Coherent integration time   = %.1f sec\n",
	  ob->date.cav_tim);
  lprintf(stdout, "  Incoherent integration time = %.1f sec\n",
	  ob->date.iav_tim);
  lprintf(stdout, "  Sum of scan durations = %ld sec\n", scansum);
  write_ut(ob->rec[0].integ->ut, sizeof(ctmpa), ctmpa),
  write_ut(ob->rec[ob->nrec-1].integ->ut, sizeof(ctmpb), ctmpb);
  lprintf(stdout, "  UT range: %s to %s\n", ctmpa, ctmpb);
/*
 * Calculate the ut of the middle of the observing run, convert it
 * to a julian day and epoch and report to user.
 */
  {
    double mid_ut = ob->rec[0].integ->ut +
      (ob->rec[ob->nrec-1].integ->ut - ob->rec[0].integ->ut)/2;
    julday(mid_ut, ob->date.year, &jd, &jdfrc, &je);
    lprintf(stdout, "  Mean epoch:  JD %ld.%3.3d = J%.3f\n",
	    jd, (int) (jdfrc*1000.0), je);
  };
  lprintf(stdout, "\n");
  return;
}

#include <stdlib.h>
#include <stdio.h>

#include "logio.h"
#include "obs.h"
#include "scans.h"

static const double mingap=0.1;    /* The minimum interscan gap (secs) */
static const int maxscan=20;       /* The maximum number of scans */

/*.......................................................................
 * Count the number of scans in a sub-array.
 *
 * Input:
 *  sub   Subarray *  The subarray descriptor.
 *  tsep    double    The gap that constitutes a scan separator (secs).
 *                    If tsep<=0, DEFGAP (see scans.h) is substituted.
 * Output:
 *  return     int    The number of scans in the subarray, or 0
 *                    on error.
 */
int nscans(Subarray *sub, double tsep)
{
  double lastut; /* UT of the previous integration looked at */
  double newut;  /* UT of the latest integration looked at */
  int nscan;     /* Number of scans in subarray */
  int ut;        /* The index of the latest integration */
/*
 * Subarray descriptor valid?
 */
  if(sub_bad(sub, "nscans"))
    return 0;
/*
 * Enforce a sensible minimum on the time separator.
 */
  if(tsep < mingap)
    tsep = DEFGAP;
/*
 * Count scans.
 */
  lastut = sub->integ[0].ut;
  nscan = 1;
  for(ut=0; ut<sub->ntime; ut++) {
    newut = sub->integ[ut].ut;
    if(newut-lastut > tsep)
      nscan++;
    lastut = newut;
  };
/*
 * Return the scan count.
 */
  return nscan;
}

/*.......................................................................
 * Return the index of the end UT of the scan containing a given UT
 * index. This can be used to determine the bounds of all scan by
 * repeated calls.
 *
 * Input:
 *  sub  Subarray *  The subarray descriptor.
 *  tsep   double    The gap that constitutes a scan separator (secs).
 *                   (Use the same value as sent to nscans().). 
 *                   If tsep<=0, DEFGAP (see scans.h) is substituted.
 *  uta       int    The index of an integration in the scan whose
 *                   end is to be found.
 * Output:
 *  return    int    The index of the end integration of the scan or
 *                   -1 on error.
 */
int endscan(Subarray *sub, double tsep, int uta)
{
  double lastut; /* UT of last integration looked at */
  double newut;  /* UT of latest integration looked at */
  int ut;        /* The index of the integration being checked */
/*
 * Check arguments.
 */
  if(sub_bad(sub, "endscan"))
    return -1;
  if(uta<0 || uta >= sub->ntime) {
    lprintf(stderr, "endscan: UT index out of bounds\n");
    return -1;
  };
/*
 * Enforce a sensible minimum on the time separator.
 */
  if(tsep < mingap)
    tsep = DEFGAP;
/*
 * Find the end of the current scan.
 */
  lastut = sub->integ[uta].ut;
  for(ut=uta; ut<sub->ntime; ut++) {
    newut = sub->integ[ut].ut;
    if(newut-lastut > tsep)
      break;
    lastut=newut;
  };
  return ut-1;
}

/*.......................................................................
 * Determine the total duration of all scans.
 *
 * Input:
 *  sub    Subarray *  The subarray descriptor.
 *  tsep     double    The gap that constitutes a scan separator (secs).
 *                     If tsep<=0, DEFGAP (see scans.h) is substituted.
 * Output:
 *  return      int    The subarray duration in seconds - or 0 on
 *                      error.
 */
int timescans(Subarray *sub, double tsep)
{
  double utdiff;  /* The difference between two UTs (seconds) */
  int uta,utb;    /* Indexes of start and end integrations in a scan */
  int scan;       /* The index of the latest scan */
  int nscan;      /* The number of scans in the subarray */
  int duration=0; /* Sum of scan durations */
/*
 * Subarray descriptor valid?
 */
  if(sub_bad(sub, "timescans"))
    return 0;
/*
 * Count the number of scans in the subarray.
 */
  nscan = nscans(sub, tsep);
/*
 * Determine the sum of scan durations.
 */
  for(uta=utb=scan=0; uta<sub->ntime && scan<nscan; scan++,uta=utb+1) {
    utb = endscan(sub, tsep, uta);
    utdiff = sub->integ[utb].ut - sub->integ[uta].ut;
    if(utdiff>0.0)
      duration += utdiff;
  };
  return duration;
}

/*.......................................................................
 * Check and record a new interscan gap in a single sub-array, or all
 * sub-arrays of an observation.
 *
 * If the selected gap would produce more than 'maxscan' scans in any one
 * sub-array, then it will be rejected and an error condition returned.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 *  gap      double    The interscan gap (seconds). A value of gap<=1.0
 *                     selects the default interscan gap.
 *  isub        int    The index of the sub-array to assign the gap
 *                     to, or -1 to select all sub-arrays.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int scangap(Observation *ob, double gap, int isub)
{
  int isa,isb;   /* The range of sub-arrays to be checked */
  int is;        /* The index of the sub-array being checked */
/* 
 * Check the validity of the observation.
 */
  if(!ob_ready(ob, OB_INDEX, "scangap"))
    return 1;
/*
 * Check the sub-array index.
 */
  if(isub >= ob->nsub || isub < -1) {
    lprintf(stderr, "scangap: Sub-array index out of range.\n");
    return 1;
  };
/*
 * Determine the range of sub-arrays to be checked.
 */
  if(isub<0) {
    isa = 0;
    isb = ob->nsub-1;
  } else {
    isa = isb = isub;
  };
/*
 * Set the default interscan gap if requested.
 */
  if(gap < mingap)
    gap = DEFGAP;
/*
 * Check each of the selected sub-arrays.
 */
  for(is=isa; is<=isb; is++) {
    if(nscans(&ob->sub[is], gap) > maxscan) {
      lprintf(stderr, "Scangap: Interscan gap too short.\n");
      return 1;
    };
  };
/*
 * Assign the new scan interval.
 */
  for(is=isa; is<=isb; is++)
    ob->sub[is].scangap = gap;
/*
 * Report the change.
 */
  if(isub == -1) {
    lprintf(stdout,
    "Delimiting interscan gap changed to %g seconds in all sub-arrays.\n", gap);
  } else {
    lprintf(stdout,
	    "Delimiting interscan gap changed to %g seconds in sub-array %d.\n",
	    gap, isub+1);
  };
  return 0;
}

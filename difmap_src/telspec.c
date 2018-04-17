#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "obs.h"
#include "telspec.h"
#include "logio.h"

enum {MAX_TS_LEN=80};  /* Max length of user input line */
enum {MAX_NTEL=5};     /* Max number of telescope name components in a */
		       /* string */

/* Declare a container type for the return value of read_tspec() */

typedef struct {
  int nfix;         /* The number of component items read from spec string */
  int isub;          /* The index of the given, or default sub-array */
  int tel[MAX_NTEL]; /* The 'nfix-1' indexes of specified telescopes */
} Tspec;

static Tspec *read_tspec(Observation *ob, char *s, char **endp, int d_sub,
			 int maxtel, char *name);
static int write_tspec(Observation *ob, Tspec *tspec, int nref, int fixref,
		       int n, char *s);
static void clr_tspec(Tspec *tspec);

/*.......................................................................
 * Read a user sub-array specification. Note that you should call
 * next_sub(ob, FIND_FIRST, ...) before using the returned specification.
 *
 * Input:
 *  ob    Observation *   The descriptor of the observation.
 *  s            char *   The string to read the specification from, or
 *                        NULL if the specification should be read from
 *                        stdin.
 *  endp         char **  Pointer to the next unprocessed char in 's',
 *                        or NULL if s==NULL. If endp==NULL then trailing
 *                        characters following a valid specification will
 *                        be interpretted as an error and reported thus.
 *  d_sub         int     The sub-array specification to be substituted
 *                        if the user neglects to specify one.
 * Output:
 *  return    Subspec *   A pointer to a static internal result container,
 *                        or NULL on error.
 */
Subspec *read_Subspec(Observation *ob, char *s, char **endp, int d_sub)
{
  static Subspec ss;      /* Sub-array description container to be returned */
/*
 * Read the sub-array spec.
 */
  Tspec *tspec = read_tspec(ob, s, endp, d_sub, 0, "sub-array");
  if(tspec) {
    ss.nfix = tspec->nfix;
    ss.isub = tspec->isub;
    return &ss;
  };
/*
 * Invalid sub-array spec.
 */
  return NULL;
}

/*.......................................................................
 * Write a sub-array specification string.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation to which the
 *                     specification refers.
 *  ss      Subspec *  The sub-array specification to be written.
 *  nref        int    If nref > ss->nfix and fixref is true, then nref
 *                     initial indexes will be displayed instead of ss->nfix.
 *  fixref      int    See 'nref'.
 *  n           int    The max number of characters (including '\0') to
 *                     place in s[].
 * Input/Output:
 *  s          char *  An array of at least 'n' characters for the return
 *                     string.
 */
int write_Subspec(Observation *ob, Subspec *ss, int nref, int fixref,
		  int n, char *s)
{
  Tspec tspec;
  clr_tspec(&tspec);
  tspec.nfix = ss->nfix;
  tspec.isub  = ss->isub;
  return write_tspec(ob, &tspec, nref, fixref, n, s);
}

/*.......................................................................
 * Search for the first/last sub-array that matches the given
 * specification. This is a convenient interface to the FIND_FIRST
 * operation of next_sub(). See its header comments for details.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 *  nfix        int    The number of fixed indexes, in the order isub.
 *                     Note that nfix must lie in the range 0..1.
 *  isub        int    The index of the sub-array to locate the sub-array in.
 *                     (fixed if nfix>=1).
 *  forward     int    The direction in which to search for sub-arrays.
 *                      0 - search backward.
 *                      1 - search forward.
 *  nref        int    This is the number of indexes in the initial
 *                     reference prefix wrt which sub-arrays are
 *                     sought and should normally equal nfix.
 *  fixref      int    If true and nref>nfix, then treat the
 *                     first nref indexes as fixed (without
 *                     changing the recorded nfix).
 *  report      int    If report is true and no sub-array is found report
 *                     the failure via stderr.
 * Output:
 *  return Subspec *  The pointer to an internal static object containing
 *                     the telescope and sub-array indexes of the latest
 *                     sub-array, or NULL if no sub-array was found.
 */
Subspec *find_sub(Observation *ob, int nfix, int isub, int forward,
		  int nref, int fixref, int report)
{
  static Subspec ss;  /* The new sub-array spec */
/* 
 * Initialize the sub-array specification descriptor.
 */
  ss.nfix = nfix;
  ss.isub = isub;
/*
 * Handle the specified operation.
 */
  if(next_sub(ob, FIND_FIRST, forward, nref, fixref, report, &ss)==0)
    return &ss;
  return NULL;
}

/*.......................................................................
 * Return the indexes of a sub-array consistent with given search limits. 
 *
 * The sub-array index will be held fixed if ss->nfix>0.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 *  oper     Findop    Specify the type of sub-array search wanted.
 *                      FIND_FIRST -  Locate the first/last sub-array.
 *                      SKIP_SUB   -  Skip to the previous/next sub-array.
 *                      FIND_NEXT  -  Find the prev/next sub-array matching
 *                                    the sub-array specification.
 *  forward     int    The direction to search for sub-arrays.
 *                      0 - search backward.
 *                      1 - search forward.
 *  nref        int    This is the number of indexes in the initial
 *                     reference prefix wrt which sub-arrays are
 *                     sought and should normally equal ss->nfix.
 *  fixref      int    If true and nref>ss->nfix, then treat the
 *                     first nref indexes as fixed (without
 *                     changing the recorded ss->nfix).
 *  report      int    If report is true and no sub-array is found report
 *                     the failure via stderr.
 * Input/Output:
 *  ss      Subspec *  The sub-array specification to be processed.
 *                     This will not be changed if an error is
 *                     detected, or no new sub-array is found.
 * Output:
 *  return             0 - Sub-Array found (recorded in *ss).
 *                     1 - Hit limits of fixed indexes before finding
 *                         sub-array.
 *                    -1 - Error.
 */
int next_sub(Observation *ob, Findop oper, int forward, int nref, int fixref,
	     int report, Subspec *ss)
{
  Subspec ct; /* Internal work container */
  int nfixed; /* The number of fixed indexes (isub) */
/*
 * Sanity check the inputs.
 */
  if(!ob_ready(ob, OB_INDEX, "next_sub"))
    return -1;
  if(ss==NULL) {
    lprintf(stderr, "next_sub: NULL Specification descriptor.\n");
    return -1;
  };
/*
 * Check the number of fixed indexes.
 */
  if(ss->nfix < 0 || ss->nfix > 1) {
    lprintf(stderr, "next_sub: Can\'t handle nfix=%d.\n", ss->nfix);
    return -1;
  };
/*
 * Copy the input specification.
 */
  ct = *ss;
/*
 * How many indexes are fixed by default.
 */
  nfixed = (fixref && nref > ss->nfix) ? nref : ss->nfix;
/*
 * Pre-position?
 */
  switch(oper) {
  default:
  case FIND_FIRST:
    break;
  case SKIP_SUB:
  case FIND_NEXT:
    if(forward)
      ct.isub++;
    else
      ct.isub--;
    break;
  };
/*
 * If all indexes are fixed, only the specified sub-array will do.
 */
  if(!(oper==FIND_NEXT && nfixed>=1) &&
     (ct.isub >= 0 && ct.isub < ob->nsub)) {
    *ss = ct;
    return 0;
  };
/*
 * Report the failed search.
 */
  if(report) {
    switch(oper) {
    case FIND_FIRST:
      lprintf(stderr, "No sub-arrays match");
      if(ss->nfix > 0)
	lprintf(stderr, " %d:", ss->isub+1);
      lprintf(stderr, ".\n");
      break;
    case SKIP_SUB:
      lprintf(stderr, "No sub-arrays found %s sub-array %d.\n",
	      forward ? "beyond":"prior to", ss->isub+1);
      break;
    case FIND_NEXT:
      lprintf(stderr, "All sub-arrays processed.\n");
      break;
    default:
      lprintf(stderr, "next_sub: Unknown operation.\n");
      break;
    };
  };
  return 1;  /* Hit fixed index limits */
}

/*.......................................................................
 * Read a user telescope specification. Note that you should call
 * next_tel(ob, FIND_FIRST, ...) before using the returned specification.
 *
 * Input:
 *  ob    Observation *   The descriptor of the observation.
 *  s            char *   The string to read the specification from, or
 *                        NULL if the specification should be read from
 *                        stdin.
 *  endp         char **  Pointer to the next unprocessed char in 's',
 *                        or NULL if s==NULL. If endp==NULL then trailing
 *                        characters following a valid specification will
 *                        be interpretted as an error and reported thus.
 *  d_sub         int     The sub-array specification to be substituted
 *                        if the user neglects to specify one.
 * Output:
 *  return    Telspec *   A pointer to a static internal result container,
 *                        or NULL on error.
 */
Telspec *read_Telspec(Observation *ob, char *s, char **endp, int d_sub)
{
  static Telspec ts;      /* Telescope description container to be returned */
/*
 * Read the telescope spec.
 */
  Tspec *tspec = read_tspec(ob, s, endp, d_sub, 1, "telescope");
  if(tspec) {
    ts.nfix = tspec->nfix;
    ts.isub = tspec->isub;
    ts.ta = tspec->tel[0];
    return &ts;
  };
/*
 * Not a valid telescope spec.
 */
  return NULL;
}

/*.......................................................................
 * Write a telescope specification string.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation to which the
 *                     specification refers.
 *  ts      Telspec *  The telescope specification to be written.
 *  nref        int    If nref > ts->nfix and fixref is true, then nref
 *                     initial indexes will be displayed instead of ts->nfix.
 *  fixref      int    See 'nref'.
 *  n           int    The max number of characters (including '\0') to
 *                     place in s[].
 * Input/Output:
 *  s          char *  An array of at least 'n' characters for the return
 *                     string.
 */
int write_Telspec(Observation *ob, Telspec *ts, int nref, int fixref,
		  int n, char *s)
{
  Tspec tspec;
  clr_tspec(&tspec);
  tspec.nfix  = ts->nfix;
  tspec.isub   = ts->isub;
  tspec.tel[0] = ts->ta;
  return write_tspec(ob, &tspec, nref, fixref, n, s);
}

/*.......................................................................
 * Search for the first/last telescope that matches the given
 * specification. This is a convenient interface to the FIND_FIRST
 * operation of next_tel(). See its header comments for details.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 *  nfix        int    The number of fixed indexes, in the order isub,ta.
 *                     Note that nfix must lie in the range 0..2.
 *  isub        int    The index of the sub-array to locate the telescope in.
 *                     (fixed if nfix>=1).
 *  ta          int    The 1st telescope index (fixed if nfix>=2).
 *  forward     int    The direction to search for telescopes.
 *                      0 - search backward.
 *                      1 - search forward.
 *  nref        int    This is the number of indexes in the initial
 *                     reference prefix wrt which telescopes are
 *                     sought and should normally equal nfix.
 *  fixref      int    If true and nref>nfix, then treat the
 *                     first nref indexes as fixed (without
 *                     changing the recorded nfix).
 *  report      int    If report is true and no telescope is found report
 *                     the failure via stderr.
 * Output:
 *  return  Telspec *  The pointer to an internal static object containing
 *                     the sub-array and telescope indexes of the latest
 *                     telescope, or NULL if no telescope was found.
 */
Telspec *find_tel(Observation *ob, int nfix, int isub, int ta, int forward,
		  int nref, int fixref, int report)
{
  static Telspec ts;  /* The new telescope spec */
/* 
 * Initialize the telescope specification descriptor.
 */
  ts.nfix = nfix;
  ts.isub = isub;
  ts.ta = ta;
/*
 * Handle the specified operation.
 */
  if(next_tel(ob, FIND_FIRST, forward, nref, fixref, report, &ts)==0)
    return &ts;
  return NULL;
}

/*.......................................................................
 * Return the sub-array and telescope indexes of a telescope
 * consistent with given search limits.
 *
 * The first nfix sub-array and telescope indexes (isub,ta) will
 * be held fixed, while the remaining (2 - nfix) indexes will optionally
 * be incremented or decremented as necessary until a valid telescope is
 * located.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 *  oper     Findop    Specify the type of telescope search wanted.
 *                      FIND_FIRST -  Locate the first/last telescope that
 *                                    matches the given telescope
 *                                    specification.
 *                      SKIP_SUB   -  Skip to the previous/next sub-array.
 *                      SKIP_TA    -  Skip to the previous/next ta index.
 *                      FIND_NEXT  -  Find the prev/next matching telescope.
 *                                    (excluding the current telescope).
 *  forward     int    The direction to search for telescopes.
 *                      0 - search backward.
 *                      1 - search forward.
 *  nref        int    This is the number of indexes in the initial
 *                     reference prefix wrt which telescopes are
 *                     sought and should normally equal ts->nfix.
 *  fixref      int    If true and nref>ts->nfix, then treat the
 *                     first nref indexes as fixed (without
 *                     changing the recorded ts->nfix).
 *  report      int    If report is true and no telescope is found report
 *                     the failure via stderr.
 * Input/Output:
 *  ts     Telspec *  The telescope specification to be processed.
 *                     This will not be changed if an error is
 *                     detected, or no new telescope is found.
 * Output:
 *  return             0 - Telescope found (recorded in *ts).
 *                     1 - Hit limits of fixed indexes before finding
 *                         telescope.
 *                    -1 - Error.
 */
int next_tel(Observation *ob, Findop oper, int forward, int nref, int fixref,
	     int report, Telspec *ts)
{
  Telspec ct; /* Internal work container */
  int nfixed; /* The number of fixed indexes (isub,ta) */
/*
 * Sanity check the inputs.
 */
  if(!ob_ready(ob, OB_INDEX, "next_tel"))
    return -1;
  if(ts==NULL) {
    lprintf(stderr, "next_tel: NULL Specification descriptor.\n");
    return -1;
  };
/*
 * Check the number of fixed indexes.
 */
  if(ts->nfix < 0 || ts->nfix > 2) {
    lprintf(stderr, "next_tel: Can\'t handle nfix=%d.\n", ts->nfix);
    return -1;
  };
/*
 * Check the validity of the sub-array index.
 */
  if(ts->isub <0 || ts->isub >= ob->nsub) {
    lprintf(stderr, "next_tel: Out of range sub-array index.\n");
    return -1;
  };
/*
 * Copy the input specification.
 */
  ct = *ts;
/*
 * How many indexes are fixed by default.
 */
  nfixed = (fixref && nref > ts->nfix) ? nref : ts->nfix;
/*
 * Pre-position?
 */
  switch(oper) {
  default:
  case FIND_FIRST:
/*
 * Fill in the start telescope index if necessary.
 */
    if(nfixed <= 1)
      ct.ta = forward ? 0 : ob->sub[ts->isub].nstat - 1;
    break;
  case SKIP_SUB:
    ct.nfix = nfixed = 0;
    ct.ta = forward ? ob->sub[ts->isub].nstat : -1;
    break;
  case SKIP_TA:
    if(nfixed > 1)
      nfixed = 1;
    if(ct.nfix > 1)
      ct.nfix = 1;
    if(forward)
      ct.ta++;
    else
      ct.ta--;
    break;
  case FIND_NEXT:
    if(forward)
      ct.ta++;
    else
      ct.ta--;
    break;
  };
/*
 * If all telescopes are fixed, only the specified telescope will do.
 */
  if(!(oper==FIND_NEXT && nfixed>=2)) {
/*
 * Search forward through the remaining sub-arrays.
 */
    if(forward) {
      for( ; ct.isub < ob->nsub; ct.isub++) {
	Subarray *sub = ob->sub + ct.isub;
/*
 * Start a new sub-array?
 */
	if(ct.isub != ts->isub)
	  ct.ta = 0;
/*
 * Search forward for the next telescope.
 */
	if(ct.ta >= 0 && ct.ta < sub->nstat) {
	  *ts = ct;
	  return 0;
	};
	if(nfixed >= 1)
	  break;
      };
/*
 * Search backward through the preceding sub-arrays.
 */
    } else {
      for( ; ct.isub >= 0; ct.isub--) {
	Subarray *sub = ob->sub + ct.isub;
/*
 * Start a new sub-array?
 */
	if(ct.isub != ts->isub)
	  ct.ta = sub->nstat-1;
/*
 * Search backwards for the next telescope.
 */
	if(ct.ta >= 0 && ct.ta < sub->nstat) {
	  *ts = ct;
	  return 0;
	};
	if(nfixed >= 1)
	  break;
      };
    };
  };
/*
 * Report the failed search.
 */
  if(report) {
    switch(oper) {
    case FIND_FIRST:
      lprintf(stderr, "No telescopes match");
      if(ts->nfix > 0)
	lprintf(stderr, " %d:", ts->isub+1);
      if(ts->nfix > 1)
	lprintf(stderr, "%s", ob->sub[ts->isub].tel[ts->ta].name);
      lprintf(stderr, ".\n");
      break;
    case SKIP_SUB:
      lprintf(stderr, "No telescopes found in sub-arrays %s sub-array %d.\n",
	      forward ? "beyond":"prior to", ts->isub+1);
      break;
    case SKIP_TA:
    case FIND_NEXT:
      lprintf(stderr, "No telescopes found %s telescope %d:%s.\n",
	    forward?"beyond":"prior to", ts->isub+1,
	      ob->sub[ts->isub].tel[ts->ta].name);
      break;
    default:
      lprintf(stderr, "next_tel: Unknown operation.\n");
      break;
    };
  };
  return 1;  /* Hit fixed index limits */
}

/*.......................................................................
 * Read a user baseline specification. Note that you should call
 * next_base(ob, FIND_FIRST, ...) before using the returned specification.
 *
 * Input:
 *  ob    Observation *   The descriptor of the observation.
 *  s            char *   The string to read the specification from, or
 *                        NULL if the specification should be read from
 *                        stdin.
 *  endp         char **  Pointer to the next unprocessed char in 's',
 *                        or NULL if s==NULL. If endp==NULL then trailing
 *                        characters following a valid specification will
 *                        be interpretted as an error and reported thus.
 *  d_sub         int     The sub-array specification to be substituted
 *                        if the user neglects to specify one.
 * Output:
 *  return   Basespec *   A pointer to a static internal result container,
 *                        or NULL on error.
 */
Basespec *read_Basespec(Observation *ob, char *s, char **endp, int d_sub)
{
  static Basespec bs;     /* Baseline description container to be returned */
/*
 * Read the telescope spec.
 */
  Tspec *tspec = read_tspec(ob, s, endp, d_sub, 2, "baseline");
  if(tspec) {
    bs.nfix = tspec->nfix;
    bs.isub = tspec->isub;
    bs.ta = tspec->tel[0];
    bs.tb = tspec->tel[1];
    bs.base = 0;
    return &bs;
  };
/*
 * Invalid specification string.
 */
  return NULL;
}

/*.......................................................................
 * Write a baseline specification string.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation to which the
 *                     specification refers.
 *  bs     Basespec *  The baseline specification to be written.
 *  nref        int    If nref > bs->nfix and fixref is true, then nref
 *                     initial indexes will be displayed instead of bs->nfix.
 *  fixref      int    See 'nref'.
 *  n           int    The max number of characters (including '\0') to
 *                     place in s[].
 * Input/Output:
 *  s          char *  An array of at least 'n' characters for the return
 *                     string.
 */
int write_Basespec(Observation *ob, Basespec *bs, int nref, int fixref,
		   int n, char *s)
{
  Tspec tspec;
  clr_tspec(&tspec);
  tspec.nfix  = bs->nfix;
  tspec.isub   = bs->isub;
  tspec.tel[0] = bs->ta;
  tspec.tel[1] = bs->tb;
  return write_tspec(ob, &tspec, nref, fixref, n, s);
}

/*.......................................................................
 * Search for the first/last baseline that matches the given
 * specification. This is a convenient interface to the FIND_FIRST
 * operation of next_base(). See its header comments for details.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 *  nfix        int    The number of fixed indexes, in the order isub,ta,tb.
 *                     Note that nfix must lie in the range 0..3.
 *  isub        int    The index of the sub-array to locate the baseline in.
 *                     (fixed if nfix>=1).
 *  ta          int    The 1st telescope index (fixed if nfix>=2).
 *  tb          int    The 2nd telescope index (fixed if nfix>=3).
 *  forward     int    The direction to search for baselines.
 *                      0 - search backward.
 *                      1 - search forward.
 *  nref        int    This is the number of indexes in the initial
 *                     reference prefix wrt which baselines are
 *                     sought and should normally equal nfix.
 *                     See the discussion of 'allref'.
 *  allref      int    If true (1), all unique baselines that match
 *                     in the first 'nref' indexes are acceptable
 *                     (ie. if nref <= nfix, and index nref+1 is a
 *                     telescope index, then the latter will be 
 *                     iterated throughout its whole range).
 *                     If false (0), only baselines in which all
 *                     free telescope indexes are in ascending
 *                     order with respect to the 'nref'th index,
 *                     are acceptable.
 *  fixref      int    If true and nref>nfix, then treat the
 *                     first nref indexes as fixed (without
 *                     changing the recorded nfix).
 *  report      int    If report is true and no baseline is found report
 *                     the failure via stderr.
 * Output:
 *  return Basespec *  The pointer to an internal static object containing
 *                     the telescope and baseline indexes of the latest
 *                     baseline, or NULL if no baseline was found.
 */
Basespec *find_base(Observation *ob, int nfix, int isub, int ta, int tb,
		    int forward, int nref, int allref, int fixref, int report)
{
  static Basespec bs;  /* The new baseline spec */
/* 
 * Initialize the baseline specification descriptor.
 */
  bs.nfix = nfix;
  bs.isub = isub;
  bs.ta = ta;
  bs.tb = tb;
/*
 * Handle the specified operation.
 */
  if(next_base(ob, FIND_FIRST, forward, nref, allref, fixref, report, &bs)==0)
    return &bs;
  return NULL;
}

/*.......................................................................
 * Return the indexes of three valid telescopes that form a baseline
 * consistent with given search limits. 
 *
 * The first ts->nfix indexes (or optionally 'nref' indexes), in the
 * order [isub,ta,tb], are treated as immutable unless explicitly
 * overriden via one of the SKIP_ operations. Searches for new
 * baselines will fail if they require one of these indexes to be
 * incremented or decremented. Telescope indexes are maintained in
 * increasing order, except for the special case of the nref+1'th
 * index, which if allref is true, is allowed to range over all
 * stations in the particular sub-array.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 *  oper     Findop    Specify the type of baseline search wanted.
 *                      FIND_FIRST -  Locate the first/last baseline that
 *                                    matches the given baseline
 *                                    specification.
 *                      SKIP_SUB   -  Skip to the previous/next sub-array.
 *                                    If bs->nfix>0, it will be assigned 0.
 *                      SKIP_TA    -  Skip to the previous/next ta index.
 *                                    If bs->nfix>1, it will be assigned 1.
 *                      SKIP_TB    -  Skip to the previous/next tb index.
 *                                    If bs->nfix>2, it will be assigned 2.
 *                      FIND_NEXT  -  Find the prev/next matching baseline.
 *                                    (excluding the current baseline).
 *  forward     int    The direction to search for baselines.
 *                      0 - search backward.
 *                      1 - search forward.
 *  nref        int    This is the number of indexes in the initial
 *                     reference prefix wrt which baselines are
 *                     sought and should normally equal bs->nfix.
 *                     See above on 'allref'.
 *  allref      int    If true (1), all unique baselines that match
 *                     in the first 'nref' indexes are acceptable
 *                     (ie. if nref <= nfix, and index nref+1 is a
 *                     telescope index, then the latter will be 
 *                     iterated throughout its whole range).
 *                     If false (0), only baselines in which all
 *                     free telescope indexes are in ascending
 *                     order with respect to the 'nref'th index,
 *                     are acceptable.
 *  fixref      int    If true and nref>bs->nfix, then treat the
 *                     first nref indexes as fixed (without
 *                     changing bs->nfix).
 *  report      int    If report is true and no baseline is found report
 *                     the failure via stderr.
 * Input/Output:
 *  bs     Basespec *  The baseline specification to be processed.
 *                     This will not be changed if an error is
 *                     detected, or no new baseline is found.
 * Output:
 *  return             0 - Baseline found (recorded in *bs).
 *                     1 - Hit limits of fixed indexes before finding
 *                         baseline.
 *                    -1 - Error.
 */
int next_base(Observation *ob, Findop oper, int forward, int nref, int allref,
	      int fixref, int report, Basespec *bs)
{
  Basespec ct; /* Internal work container */
  int align_b; /* True if tb should equal at least ta+1 */
  int nfixed;  /* The number of fixed indexes (isub,ta,tb) */
/*
 * Sanity check the inputs.
 */
  if(!ob_ready(ob, OB_INDEX, "next_base"))
    return -1;
  if(bs==NULL) {
    lprintf(stderr, "next_base: NULL Specification descriptor.\n");
    return -1;
  };
/*
 * Check the number of fixed indexes.
 */
  if(bs->nfix < 0 || bs->nfix > 3) {
    lprintf(stderr, "next_base: Can\'t handle nfix=%d.\n", bs->nfix);
    return -1;
  };
/*
 * Check the validity of the sub-array index.
 */
  if(bs->isub <0 || bs->isub >= ob->nsub) {
    lprintf(stderr, "next_base: Out of range sub-array index.\n");
    return -1;
  };
/*
 * Copy the input specification.
 */
  ct = *bs;
/*
 * Determine which telescopes must be aligned.
 */
  align_b = 1;
  switch(nref) {
  case 2:
    align_b = !allref;
    break;
  };
/*
 * How many indexes are fixed by default.
 */
  nfixed = (fixref && nref > bs->nfix) ? nref : bs->nfix;
/*
 * Pre-position?
 */
  switch(oper) {
  default:
  case FIND_FIRST:
/*
 * Fill in starting values for any free telescope indexes.
 */
    if(forward) {
      switch(nfixed) { /* Note that each case falls through to the next */
      case 0:
      case 1:
	ct.ta = 0;
      case 2:
	ct.tb = align_b ? ct.ta+1 : 0;
	break;
      };
    } else {
      int maxtel = ob->sub[bs->isub].nstat - 1;
      switch(nfixed) { /* Note that each case falls through to the next */
      case 2:
	ct.tb = maxtel;
      case 1:
	ct.ta = align_b ? ct.tb-1 : maxtel;
	break;
      };
    };
    break;
  case SKIP_SUB:
    ct.nfix = nfixed = 0;
    ct.ta = ct.tb = forward ? ob->sub[bs->isub].nstat : -1;
    break;
  case SKIP_TA:
    if(nfixed > 1)
      nfixed = 1;
    if(ct.nfix > 1)
      ct.nfix = 1;
    ct.tb = forward ? ob->sub[bs->isub].nstat : -1;
    break;
  case SKIP_TB:
    if(nfixed > 2)
      nfixed = 2;
    if(ct.nfix > 2)
      ct.nfix = 2;
    if(forward) {
      ct.tb = (align_b && ct.tb<ct.ta ? ct.ta:ct.tb) + 1;
    } else {
      ct.tb--;
    };
    break;
  case FIND_NEXT:
    if(forward) {
      ct.tb = (align_b && ct.tb<ct.ta ? ct.ta:ct.tb) + 1;
    } else {
      ct.tb--;
    };
    break;
  };
/*
 * If all telescopes are fixed, only the specified baseline will do.
 */
  if(!(oper==FIND_NEXT && nfixed>=3)) {
/*
 * Search forward through the remaining sub-arrays.
 */
    if(forward) {
      for( ; ct.isub < ob->nsub; ct.isub++) {
	Subarray *sub = ob->sub + ct.isub;
/*
 * Start a new sub-array?
 */
	if(ct.isub != bs->isub) {
	  ct.ta = 0;
	  ct.tb = 1;
	};
/*
 * Search forward for the next baseline.
 */
	for( ; ct.ta<sub->nstat; ct.ta++, ct.tb=(align_b?ct.ta+1:0)) {
	  for( ; ct.tb<sub->nstat; ct.tb++) {
	    if((ct.base=loc_base(sub, ct.ta, ct.tb)) >= 0) {
	      *bs = ct;
	      return 0;
	    };
	    if(nfixed >= 3)
	      break;
	  };
	  if(nfixed >= 2)
	    break;
	};
	if(nfixed >= 1)
	  break;
      };
/*
 * Search backward through the preceding sub-arrays.
 */
    } else {
      for( ; ct.isub >= 0; ct.isub--) {
	Subarray *sub = ob->sub + ct.isub;
/*
 * Start a new sub-array?
 */
	if(ct.isub != bs->isub) {
	  ct.tb = sub->nstat-1;
	  ct.ta = align_b ? ct.tb-1 : sub->nstat-1;
	};
/*
 * Search backwards for the next baseline.
 */
	for( ; ct.ta>=0; ct.ta--, ct.tb=sub->nstat-1) {
	  for( ; ct.tb>(align_b?ct.ta:-1); ct.tb--) {
	    if((ct.base=loc_base(sub, ct.ta, ct.tb)) >= 0) {
	      *bs = ct;
	      return 0;
	    };
	    if(nfixed >= 3)
	      break;
	  };
	  if(nfixed >= 2)
	    break;
	};
	if(nfixed >= 1)
	  break;
      };
    };
  };
/*
 * Report the failed search.
 */
  if(report) {
    switch(oper) {
    case FIND_FIRST:
      lprintf(stderr, "No baselines match");
      if(bs->nfix > 0)
	lprintf(stderr, " %d:", bs->isub+1);
      if(bs->nfix > 1)
	lprintf(stderr, "%s", ob->sub[bs->isub].tel[bs->ta].name);
      if(bs->nfix > 2)
	lprintf(stderr, " %s", ob->sub[bs->isub].tel[bs->tb].name);
      lprintf(stderr, ".\n");
      break;
    case SKIP_SUB:
      lprintf(stderr, "No baselines found in sub-arrays %s sub-array %d.\n",
	      forward ? "beyond":"prior to", bs->isub+1);
      break;
    case SKIP_TA:
      lprintf(stderr, "No baselines found for telescopes %s %d:%s.\n",
	    forward?"beyond":"prior to", bs->isub+1,
	      ob->sub[bs->isub].tel[bs->ta].name);
      break;
    case SKIP_TB:
    case FIND_NEXT:
      lprintf(stderr, "No baselines found %s baseline %d:%s %s.\n",
	    forward?"beyond":"prior to", bs->isub+1,
	      ob->sub[bs->isub].tel[bs->ta].name,
	      ob->sub[bs->isub].tel[bs->tb].name);
      break;
    default:
      lprintf(stderr, "next_base: Unknown operation.\n");
      break;
    };
  };
  return 1;  /* Hit fixed index limits */
}

/*.......................................................................
 * Find the 0-relative baseline number in a given sub-array corresponding
 * to the two 0-relative sub-array telescope numbers in tel_a and tel_b.
 *
 * Input:
 *   sub  Subarray *  The Subarray in which tel_a and tel_b are defined.
 *   tel_a     int    The zero-relative number of either telescope
 *                    on the required baseline.
 *   tel_b     int    The other telescope to tel_a on the sought
 *                    baseline.
 * Output:
 *   return    int    The zero-relative baseline number made up of
 *                    tel_a and tel_b. If not found -1 is returned.
 */
int loc_base(Subarray *sub, int tel_a, int tel_b)
{
  int base;   /* The current baseline number being checked */
  Baseline *bptr; /* Pointer to the current baseline index element */
/*
 * Search through the baseline index array for tel_a and tel_b.
 */
  bptr = sub->base;
  for(base=0; base < sub->nbase; base++,bptr++) {
    if((bptr->tel_a == tel_a && bptr->tel_b == tel_b) ||
       (bptr->tel_a == tel_b && bptr->tel_b == tel_a))
      return base;
  };
/*
 * No match.
 */
  return -1;
}

/*.......................................................................
 * Read a user closure triangle specification. Note that you should call
 * next_tri(ob, FIND_FIRST, ...) before using the returned specification.
 *
 * Input:
 *  ob  Observation *   The descriptor of the observation.
 *  s          char *   The string to read the specification from, or
 *                      NULL if the specification should be read from
 *                      stdin.
 *  endp       char **  Pointer to the next unprocessed char in 's',
 *                      or NULL if s==NULL. If endp==NULL then trailing
 *                      characters following a valid specification will
 *                      be interpretted as an error and reported thus.
 *  d_sub       int     The sub-array specification to be substituted
 *                      if the user neglects to specify one.
 * Output:
 *  return  Trispec *   A pointer to a static internal result container,
 *                      or NULL on error.
 */
Trispec *read_Trispec(Observation *ob, char *s, char **endp, int d_sub)
{
  static Trispec ts;     /* Triangle description container to be returned */
  int i;
/*
 * Read the telescope spec.
 */
  Tspec *tspec = read_tspec(ob, s, endp, d_sub, 3, "closure triangle");
  if(tspec) {
    ts.nfix = tspec->nfix;
    ts.isub = tspec->isub;
    ts.ta = tspec->tel[0];
    ts.tb = tspec->tel[1];
    ts.tc = tspec->tel[2];
    for(i=0; i<3; i++) {
      ts.b[i].base = 0;
      ts.b[i].sign = 1;
    };
    return &ts;
  };
/*
 * Invalid triangle specification.
 */
  return NULL;
}

/*.......................................................................
 * Write a closure triangle specification string.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation to which the
 *                     specification refers.
 *  ts      Trispec *  The closure triangle specification to be written.
 *  nref        int    If nref > ts->nfix and fixref is true, then nref
 *                     initial indexes will be displayed instead of ts->nfix.
 *  fixref      int    See 'nref'.
 *  n           int    The max number of characters (including '\0') to
 *                     place in s[].
 * Input/Output:
 *  s          char *  An array of at least 'n' characters for the return
 *                     string.
 */
int write_Trispec(Observation *ob, Trispec *ts, int nref, int fixref,
		  int n, char *s)
{
  Tspec tspec;
  clr_tspec(&tspec);
  tspec.nfix  = ts->nfix;
  tspec.isub   = ts->isub;
  tspec.tel[0] = ts->ta;
  tspec.tel[1] = ts->tb;
  tspec.tel[2] = ts->tc;
  return write_tspec(ob, &tspec, nref, fixref, n, s);
}

static int endtri(Subarray *sub, Trispec *ct, Trispec *ts);

/*.......................................................................
 * Return the indexes of three valid baselines that form a closure
 * triangle consistent with given search limits.
 *
 * The first ts->nfix indexes (or optionally 'nref' indexes), in the
 * order [isub,ta,tb,tc], are treated as immutable unless explicitly
 * overriden via one of the SKIP_ operations. Searches for new
 * triangles will fail if they require one of these indexes to be
 * incremented or decremented. Telescope indexes are maintained in
 * increasing order, except for the special case of the nref+1'th
 * index, which if allref is true, is allowed to range over all
 * stations in the particular sub-array.
 *
 * Input:
 *  ob Observation *  The descriptor of the observation.
 *  oper    Findop    Specify the type of triangle search wanted.
 *                     FIND_FIRST -  Locate the first/last triangle that
 *                                   matches the given triangle specification.
 *                     SKIP_SUB   -  Skip to the previous/next sub-array.
 *                                   If bs->nfix>0, it will be assigned 0.
 *                     SKIP_TA    -  Skip to the previous/next ta index.
 *                                   If bs->nfix>1, it will be assigned 1.
 *                     SKIP_TB    -  Skip to the previous/next tb index.
 *                                   If bs->nfix>2, it will be assigned 2.
 *                     SKIP_TC    -  Skip to the previous/next tc index.
 *                                   If bs->nfix>3, it will be assigned 3.
 *                     FIND_NEXT  -  Find the prev/next matching triangle.
 *                                   (excluding the current triangle).
 *  forward    int    The direction to search for triangles.
 *                     0 - search backward.
 *                     1 - search forward.
 *  nref       int    This is the number of indexes in the initial
 *                    reference prefix wrt which triangles are
 *                    sought and should normally equal bs->nfix.
 *                    See above on 'allref'.
 *  allref     int    If true (1), all unique triangles that match
 *                    in the first 'nref' indexes are acceptable
 *                    (ie. if nref <= nfix, and index nref+1 is a
 *                    telescope index, then the latter will be 
 *                    iterated throughout its whole range).
 *                    If false (0), only triangles in which all
 *                    free telescope indexes are in ascending
 *                    order with respect to the 'nref'th index,
 *                    are acceptable.
 *  fixref     int    If true and nref>bs->nfix, then treat the
 *                    first nref indexes as fixed (without
 *                    changing bs->nfix).
 *  report     int    If report is true and no triangle is found report
 *                    the failure via stderr.
 * Input/Output:
 *  ts     Trispec *  The triangle specification to be processed.
 *                    This will not be changed if an error is
 *                    detected, or no new triangle is found.
 * Output:
 *  return            0 - Triangle found (recorded in *ts).
 *                    1 - Hit limits of fixed indexes before finding
 *                        triangle.
 *                   -1 - Error.
 */
int next_tri(Observation *ob, Findop oper, int forward, int nref, int allref,
	     int fixref, int report, Trispec *ts)
{
  Trispec ct;   /* Internal work container */
  int align_b;  /* True if tb should equal at least ta+1 */
  int align_c;  /* True if tc should equal at least tb+1 */
  int nfixed;   /* The number of fixed indexes (isub,ta,tb,tc) */
/*
 * Sanity check the inputs.
 */
  if(!ob_ready(ob, OB_INDEX, "next_tri"))
    return -1;
  if(ts==NULL) {
    lprintf(stderr, "next_tri: NULL Specification descriptor.\n");
    return -1;
  };
/*
 * Check the number of fixed indexes.
 */
  if(ts->nfix < 0 || ts->nfix > 4) {
    lprintf(stderr, "next_tri: Can\'t handle nfix=%d.\n", ts->nfix);
    return -1;
  };
/*
 * Check the validity of the sub-array index.
 */
  if(ts->isub <0 || ts->isub >= ob->nsub) {
    lprintf(stderr, "next_tri: Out of range sub-array index.\n");
    return -1;
  };
/*
 * Copy the input specification.
 */
  ct = *ts;
/*
 * Determine which telescopes must be aligned.
 */
  align_b = align_c = 1;
  switch(nref) {
  case 2:
    align_b = !allref;
    break;
  case 3:
    align_c = !allref;
    break;
  };
/*
 * How many indexes are fixed by default.
 */
  nfixed = (fixref && nref > ts->nfix) ? nref : ts->nfix;
/*
 * Pre-position?
 */
  switch(oper) {
  default:
  case FIND_FIRST:
/*
 * Fill in starting values for any free telescope indexes.
 */
    if(forward) {
      switch(nfixed) { /* Note that each case falls through to the next */
      case 0:
      case 1:
	ct.ta = 0;
      case 2:
	ct.tb = align_b ? ct.ta+1 : 0;
      case 3:
	ct.tc = align_c ? ct.tb+1 : 0;
	break;
      };
    } else {
      int maxtel = ob->sub[ts->isub].nstat - 1;
      switch(nfixed) { /* Note that each case falls through to the next */
      case 3:
	ct.tc = maxtel;
      case 2:
	ct.tb = align_c ? ct.tc-1 : maxtel;
      case 1:
	ct.ta = align_b ? ct.tb-1 : maxtel;
	break;
      };
    };
    break;
  case SKIP_SUB:
    ct.nfix = nfixed = 0;
    ct.ta = ct.tb = ct.tc = forward ? ob->sub[ts->isub].nstat : -1;
    break;
  case SKIP_TA:
    if(nfixed > 1)
      nfixed = 1;
    if(ct.nfix > 1)
      ct.nfix = 1;
    ct.tb = ct.tc = forward ? ob->sub[ts->isub].nstat : -1;
    break;
  case SKIP_TB:
    if(nfixed > 2)
      nfixed = 2;
    if(ct.nfix > 2)
      ct.nfix = 2;
    ct.tc = forward ? ob->sub[ts->isub].nstat : -1;
    break;
  case SKIP_TC:
    if(nfixed > 3)
      nfixed = 3;
    if(ct.nfix > 3)
      ct.nfix = 3;
    if(forward) {
      ct.tc = (align_c && ct.tc<ct.tb ? ct.tb:ct.tc) + 1;
    } else {
      ct.tc--;
    };
    break;
  case FIND_NEXT:
    if(forward) {
      ct.tc = (align_c && ct.tc<ct.tb ? ct.tb:ct.tc) + 1;
    } else {
      ct.tc--;
    };
    break;
  };
/*
 * If all telescopes are fixed, only the specified triangle will do.
 */
  if(!(oper==FIND_NEXT && nfixed>=4)) {
/*
 * Search forward through the remaining sub-arrays.
 */
    if(forward) {
      for( ; ct.isub < ob->nsub; ct.isub++) {
	Subarray *sub = ob->sub + ct.isub;
/*
 * Start a new sub-array?
 */
	if(ct.isub != ts->isub) {
	  ct.ta = 0;
	  ct.tb = 1;
	  ct.tc = 2;
	};
/*
 * Search forward for the next triangle.
 */
	for( ; ct.ta<sub->nstat; ct.ta++, ct.tb=(align_b?ct.ta+1:0),
	    ct.tc=(align_c?ct.tb+1:0)) {
	  for( ; ct.tb<sub->nstat; ct.tb++, ct.tc=(align_c?ct.tb+1:0)) {
	    if((ct.b[0].base=loc_base(sub, ct.ta, ct.tb)) >= 0 ) {
	      for( ; ct.tc<sub->nstat; ct.tc++) {
		if((ct.b[1].base=loc_base(sub, ct.tb, ct.tc)) >= 0 &&
		   (ct.b[2].base=loc_base(sub, ct.ta, ct.tc)) >= 0 ) {
		  return endtri(sub, &ct, ts);
		};
		if(nfixed >= 4)
		  break;
	      };
	    };
	    if(nfixed >= 3)
	      break;
	  };
	  if(nfixed >= 2)
	    break;
	};
	if(nfixed >= 1)
	  break;
      };
/*
 * Search backward through the preceding sub-arrays.
 */
    } else {
      for( ; ct.isub >= 0; ct.isub--) {
	Subarray *sub = ob->sub + ct.isub;
/*
 * Start a new sub-array?
 */
	if(ct.isub != ts->isub) {
	  ct.tc = sub->nstat-1;
	  ct.tb = align_c ? ct.tc-1 : sub->nstat-1;
	  ct.ta = align_b ? ct.tb-1 : sub->nstat-1;
	};
/*
 * Search backwards for the next triangle.
 */
	for( ; ct.ta>=0; ct.ta--, ct.tc=sub->nstat-1,
	    ct.tb=(align_c ? ct.tc-1:sub->nstat-1)) {
	  for( ; ct.tb>(align_b?ct.ta:-1); ct.tc=sub->nstat-1, ct.tb--) {
	    if((ct.b[0].base=loc_base(sub, ct.ta, ct.tb)) >= 0 ) {
	      for( ; ct.tc>(align_c?ct.tb:-1); ct.tc--) {
		if((ct.b[1].base=loc_base(sub, ct.tb, ct.tc)) >= 0 &&
		   (ct.b[2].base=loc_base(sub, ct.ta, ct.tc)) >= 0 ) {
		  return endtri(sub, &ct, ts);
		};
		if(nfixed >= 4)
		  break;
	      };
	    };
	    if(nfixed >= 3)
	      break;
	  };
	  if(nfixed >= 2)
	    break;
	};
	if(nfixed >= 1)
	  break;
      };
    };
  };
/*
 * Report the failed search.
 */
  if(report) {
    switch(oper) {
    case FIND_FIRST:
      if(ts->nfix < 1) {
	lprintf(stderr, "No triangles found.\n");
      } else {
	lprintf(stderr, "No triangles match");
	if(ts->nfix > 0)
	  lprintf(stderr, " %d:", ts->isub+1);
	if(ts->nfix > 1)
	  lprintf(stderr, "%s", ob->sub[ts->isub].tel[ts->ta].name);
	if(ts->nfix > 2)
	  lprintf(stderr, " %s", ob->sub[ts->isub].tel[ts->tb].name);
	if(ts->nfix > 3)
	  lprintf(stderr, " %s", ob->sub[ts->isub].tel[ts->tc].name);
	lprintf(stderr, ".\n");
      };
      break;
    case SKIP_SUB:
      lprintf(stderr, "No triangles found in sub-arrays %s sub-array %d.\n",
	      forward ? "beyond":"prior to", ts->isub+1);
      break;
    case SKIP_TA:
      lprintf(stderr, "No triangles found for telescopes %s %d:%s.\n",
	    forward?"beyond":"prior to", ts->isub+1,
	      ob->sub[ts->isub].tel[ts->ta].name);
      break;
    case SKIP_TB:
      lprintf(stderr, "No triangles found for baselines %s %d:%s %s.\n",
	    forward?"beyond":"prior to", ts->isub+1,
	      ob->sub[ts->isub].tel[ts->ta].name,
	      ob->sub[ts->isub].tel[ts->tb].name);
      break;
    case SKIP_TC:
    case FIND_NEXT:
      lprintf(stderr, "No triangles found %s triangle %d:%s %s %s.\n",
	    forward?"beyond":"prior to", ts->isub+1,
	      ob->sub[ts->isub].tel[ts->ta].name,
	      ob->sub[ts->isub].tel[ts->tb].name,
	      ob->sub[ts->isub].tel[ts->tc].name);
      break;
    default:
      lprintf(stderr, "next_tri: Unknown operation.\n");
      break;
    };
  };
  return 1;  /* Hit fixed index limits */
}

/*.......................................................................
 * Private return function of next_tri(). If the passed descriptor is not
 * NULL, its baseline phase signs are filled in.
 */
static int endtri(Subarray *sub, Trispec *ct, Trispec *ts)
{
/*
 * If any of the baselines are directed in the opposite direction to
 * the requested directions, arrange for their signs to be toggled.
 */
  ct->b[0].sign = sub->base[ct->b[0].base].tel_a == ct->ta ? 1 : -1;
  ct->b[1].sign = sub->base[ct->b[1].base].tel_a == ct->tb ? 1 : -1;
  ct->b[2].sign = sub->base[ct->b[2].base].tel_a == ct->tc ? 1 : -1;
/*
 * Copy the located triangle specification into the return container.
 */
  *ts = *ct;
  return 0;
}

/*.......................................................................
 * Search for the first/last triangle that matches the given
 * specification. This is a convenient interface to the FIND_FIRST
 * operation of next_tri(). See its header comments for details.
 *
 * Input:
 *  ob Observation *  The descriptor of the observation.
 *  nfix      int     The number of fixed indexes, in the order isub,ta,tb,tc.
 *                    Note that nfix must lie in the range 0..4.
 *  isub       int    The index of the sub-array to locate the triangle in.
 *                    (fixed if nfix>=1).
 *  ta         int    The 1st telescope index (fixed if nfix>=2).
 *  tb         int    The 2nd telescope index (fixed if nfix>=3).
 *  tc         int    The 3rd telescope index (fixed if nfix>=4).
 *  forward    int    The direction to search for triangles.
 *                     0 - search backward.
 *                     1 - search forward.
 *  nref       int    This is the number of indexes in the initial
 *                    reference prefix wrt which triangles are
 *                    sought and should normally equal nfix.
 *                    See the discussion on 'allref'.
 *  allref     int    If true (1), all unique triangles that match
 *                    in the first 'nref' indexes are acceptable
 *                    (ie. if nref <= nfix, and index nref+1 is a
 *                    telescope index, then the latter will be 
 *                    iterated throughout its whole range).
 *                    If false (0), only triangles in which all
 *                    free telescope indexes are in ascending
 *                    order with respect to the 'nref'th index,
 *                    are acceptable.
 *  fixref     int    If true and nref>nfix, then treat the
 *                    first nref indexes as fixed (without
 *                    changing the recorded nfix).
 *  report     int    If report is true and no triangle is found report
 *                    the failure via stderr.
 * Output:
 *  return Trispec *  The pointer to an internal static object containing
 *                    the telescope and baseline indexes of the latest
 *                    triangle, or NULL if no triangle was found.
 */
Trispec *find_tri(Observation *ob, int nfix, int isub, int ta, int tb, int tc,
		  int forward, int nref, int allref, int fixref, int report)
{
  static Trispec ts;  /* The new triangle spec */
/* 
 * Initialize the triangle specification descriptor.
 */
  ts.nfix = nfix;
  ts.isub = isub;
  ts.ta = ta;
  ts.tb = tb;
  ts.tc = tc;
/*
 * Handle the specified operation.
 */
  if(next_tri(ob, FIND_FIRST, forward, nref, allref, fixref, report, &ts)==0)
    return &ts;
  return NULL;
}

/*.......................................................................
 * Read a telescope specfication string of the form:
 *   "[sub-array:]station1[-station2][-station3]....".
 *
 * Input:
 *  ob  Observation *  The observation to which the specification refers.
 *  s          char *  The string to read the specification from, or
 *                     NULL if the specification should be read from stdin.
 *  endp       char ** Pointer to the next unprocessed char in 's',
 *                     or NULL if s==NULL. If endp==NULL then trailing
 *                     characters following a valid specification will
 *                     be interpretted as an error and reported thus.
 *  d_sub       int    The index of the sub-array to be substituted
 *                     if the user neglects to specify one.
 *  maxtel      int    The maximum number of telescope name elements
 *                     expected. Note that this must not exceed MAX_NTEL.
 *  name       char *  The name of the required telescope aggregration.
 *                     This is used in error reporting and prompting only.
 *  endp
 * Output:
 *  return    Tspec *  A pointer to an internal static structure
 *                     containing the decoded specification, or NULL
 *                     on error.
 */
static Tspec *read_tspec(Observation *ob, char *s, char **endp, int d_sub,
			 int maxtel, char *name) 
{
  static Tspec tspec;     /* Telescope specification container to be returned */
  char buf[MAX_TS_LEN+1]; /* Buffer used to read user input */
  char *sptr;             /* Pointer into the specification string */
  Subarray *sub;          /* The descriptor of the chosen sub-array */
  int finished=0;         /* True after reaching end of valid specification */
  int i;
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "read_tspec"))
    return NULL;
  if (d_sub<0 || d_sub>=ob->nsub) {
    lprintf(stderr, "read_tspec: The default sub-array 'd_sub' is invalid.\n");
    return NULL;
  };
/*
 * Initialize return values.
 */
  clr_tspec(&tspec);
  tspec.isub = d_sub;
  if(endp) *endp = s;
  if(s==NULL) endp=NULL;
/*
 * Prompt for telescope specification?
 */
  if(s) {
    sptr = s;
  } else {
    sptr = buf;
    printf("Enter %s ([sub-array:]", name);
    for(i=0; i<maxtel; i++)
      printf("[%sstation%d]", i==0 ? "":"-", i+1);
    printf("): ");
    if(fgets(buf, MAX_TS_LEN, stdin)==NULL) {
      lprintf(stderr, "read_tspec: Error reading input.\n");
      return NULL;
    };
/*
 * Remove the terminal newline character, if present.
 */
    {
      char *cptr = strchr(buf, '\n');
      if(cptr) *cptr = '\0';
    };
  };
/*
 * Skip leading white-space in the input string.
 */
  while(*sptr && isspace((int)*sptr))
    sptr++;
/*
 * Read a number from the string.
 */
  if(isdigit((int)*sptr)) {
    tspec.nfix = 1;
    tspec.isub = strtol(sptr, &sptr, 10) - 1;
    if(endp) *endp = sptr;
/*
 * Skip white-space.
 */
    while(*sptr && isspace((int)*sptr))
      sptr++;
/*
 * If the next character is not a ':' or the end of the string, then
 * we have reached the end of the specification.
 */
    if(*sptr==':') {
      sptr++;
      if(endp) *endp = sptr;
/*
 * Skip white-space.
 */
      while(*sptr && isspace((int)*sptr))
	sptr++;
    } else {
      finished = 1;
    };
/*
 * If the first character in the string is a colon, then this means
 * substitute the default sub-array as a fixed item.
 */
  } else {
    if(*sptr==':') {
      tspec.nfix = 1;
      sptr++;
      if(endp) *endp = sptr;
/*
 * Skip white-space.
 */
      while(*sptr && isspace((int)*sptr))
	sptr++;
    };    
  };
/*
 * Check the range of the sub-array prefix.
 */
  if(tspec.isub < 0 || tspec.isub >= ob->nsub) {
    lprintf(stderr, "read_tspec: Sub-array pre-fix (%d) out of range 1-%d.\n",
	    tspec.isub+1, ob->nsub);
    return NULL;
  };
/*
 * Get the descriptor of the selected sub-array.
 */
  sub = &ob->sub[tspec.isub];
/*
 * Read up to maxtel telescope name components from the specification string.
 */
  for(i=0; *sptr && !finished && i<maxtel; i++) {
    int nmatch = 0;  /* Number of (un)ambiguous matches with telescope name */
    int name_len=0;  /* The length of the new telescope name specification */
    int tel=0;       /* The index of the telescope being compared */
/*
 * Locate the end of the next telescope name component.
 */
    char *eptr = sptr;
    while(*eptr && *eptr!='-' && *eptr!='+' && *eptr!='!' &&
	  isgraph((int)*eptr)) {
      if(*eptr == '\\')   /* Skip the escape character */
	eptr++;
      if(*eptr)           /* Include the current character in the name */
	eptr++;
    };
    name_len = (eptr - sptr);
/*
 * A telescope name of '*' terminates the list and is equivalent to
 * to omitting trailing telescope names.
 */
    if(name_len==1 && *sptr=='*') {
      while(*eptr && isspace((int)*eptr))
	eptr++;
      if(endp) *endp = eptr;
      sptr = eptr;
      finished = 1;
      break;
    };
/*
 * Make a case-insensitive comparison of the new telescope name with the
 * name of each telescope in the sub->tel array, and report ambiguous matches.
 */
    for(tel=0; tel<sub->nstat; tel++) {
      char *cptr = sptr;
      char *tptr = &sub->tel[tel].name[0];
      if(*cptr == '\\')  /* Remove any initial escape character */
	  cptr++;
      while(cptr<eptr && *tptr && (*cptr == *tptr ||
           (islower((int)*cptr) && toupper((int)*cptr) == *tptr) ||
           (isupper((int)*cptr) && tolower((int)*cptr) == *tptr))) {
	cptr++;
	tptr++;
	if(*cptr == '\\')  /* Don't compare the escape character */
	  cptr++;
      };
/*
 * Count successful matches.
 */
      if(cptr==eptr) {
/*
 * Exact match?
 */
	if(*tptr == '\0') {
	  nmatch = 1;
	  tspec.tel[i] = tel;
	  break;
/*
 * First match?
 */
	} else if(++nmatch == 1) {
	  tspec.tel[i] = tel;
	} else {
/*
 * Don't start the list of ambiguous matches until a second match is seen.
 */
	  if(nmatch == 2 ) {
	    lprintf(stderr, "\'%.*s\' is ambiguous with telescopes:\n",
		    name_len, sptr);
	    lprintf(stderr, "  %s\n", sub->tel[tspec.tel[i]].name);
	  };
	  lprintf(stderr, "  %s\n", sub->tel[tel].name);
	};
      };
    };
/*
 * Abort on ambiguous match.
 */
    if(nmatch > 1) {
      return NULL;
    }
/*
 * Report if no telescope matched.
 */
    else if(nmatch == 0) {
      lprintf(stderr, "No telescope name matches \'%.*s\'.\n", name_len, sptr);
      return NULL;
    } else {
/*
 * Even if the sub-array was not explicitly specified, it becomes implicitly
 * fixed if at least one station from the default sub-array is cited.
 */
      if(!tspec.nfix)
	tspec.nfix = 1;
/*
 * Count successfully acquired specification items.
 */
      tspec.nfix++;
      if(endp) *endp = eptr;
      sptr = eptr;
/*
 * Telescope names may be separated by white-space and/or one hyphen.
 * Skip the separator if valid.
 */
      if(!isspace((int)*eptr) && *eptr != '-') {
	finished = 1;
      } else {
	while(*sptr && isspace((int)*sptr))
	  sptr++;
	if(*sptr == '-')
	  sptr++;
	while(*sptr && isspace((int)*sptr))
	  sptr++;
/*
 * Check for baseline selection list separators.
 */
	if(*sptr=='+' || *sptr=='!')
	  finished = 1;
      };
    };
  };
/*
 * If endp==NULL, then regard trailing non-white-space characters
 * as erroneous. Otherwise, assume that the caller wants to deal with
 * trailing input.
 */
  if(endp==NULL) {
    while(*sptr && isspace((int)*sptr))
      sptr++;
    if(*sptr) {
      lprintf(stderr,
	     "read_tspec: Garbage follows %s specification (\"%s\").\n",
	      name, sptr);
      return NULL;
    };
  };
  return &tspec;
}

/*.......................................................................
 * Write a telescope-aggregate specification string.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation to which the
 *                     specification refers.
 *  tspec     Tspec *  The telescope-aggregate specification to be written.
 *  nref        int    If nref > tspec->nfix and fixref is true, then nref
 *                     initial indexes will be displayed instead of tspec->nfix.
 *  fixref      int    See 'nref'.
 *  n           int    The max number of characters (including '\0') to
 *                     place in s[].
 * Input/Output:
 *  s          char *  An array of at least 'n' characters for the return
 *                     string.
 * Output:
 *  return      int    The number of characters (excluding '\0') written
 *                     to s[], or -1 on error, and -2 on truncation.
 */
static int write_tspec(Observation *ob, Tspec *tspec, int nref, int fixref,
		       int n, char *s)
{
  char buf[MAX_TS_LEN+1]; /* Work buffer */
  int ndone=0;  /* The total number of characters written */
  int nnew;     /* The latest additional number of characters written */
  int nfix;     /* The number of fixed indexes to display */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "write_tspec"))
    return -1;
  if(tspec==NULL || s==NULL)
    return -1;
/*
 * How many indexes are to be displayed?
 */
  nfix = fixref && nref > tspec->nfix ? nref : tspec->nfix;
/*
 * Is the specification empty?
 */
  if(nfix <= 0) {
    if(n < 2) {
      lprintf(stderr, "write_tspec: String too short.\n");
      return -2;    /* Truncated string */
    };
    strcpy(s, "*");
    ndone = 1;
  }
/*
 * Does the specification contain a sub-array index?
 */
  else if(nfix >= 1) {
    int isub = tspec->isub;
/*
 * Check the sub-array index spec.
 */
    if(isub < 0 || isub >= ob->nsub) {
      lprintf(stderr, "write_tspec: Sub-array index out of range.\n");
      return -1;
    };
/*
 * Write the sub-array spec if there is room.
 */
    sprintf(buf, "%d:", isub+1);
    nnew = strlen(buf);
    if(ndone+nnew >= n) {
      return -2;           /* Truncated string */
    } else {
      Subarray *sub = ob->sub + isub;
      int tel;
      strcpy(s+ndone, buf);
      ndone += nnew;
/*
 * Write nfix-1 telescope components.
 */
      for(tel=0; ndone>=0 && tel<nfix-1; tel++) {
	int itel = tspec->tel[tel];
/*
 * Check the telescope index spec.
 */
	if(itel < 0 || itel >= sub->nstat) {
	  lprintf(stderr, "write_tspec: Telescope index out of range.\n");
	  return -1;
	};
/*
 * Write the latest telescope name spec if there is room.
 */
	sprintf(buf, "%s%s", tel==0?"":"-", sub->tel[itel].name);
	nnew = strlen(buf);
	if(ndone+nnew >= n)
	  return -2;           /* Truncated string */
	strcpy(s+ndone, buf);
	ndone += nnew;
      };
    };
  };
  return ndone;
}

/*.......................................................................
 * Clear a Tspec descriptor.
 *
 * Input/Output:
 *  tspec    Tspec *  The descriptor to be cleared.
 */
static void clr_tspec(Tspec *tspec)
{
  int i;
  tspec->nfix = 0;
  tspec->isub = 0;
  for(i=0; i<MAX_NTEL; i++)
    tspec->tel[i] = 0;
  return;
}


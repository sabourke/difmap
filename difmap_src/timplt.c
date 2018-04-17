#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "obs.h"
#include "vlbconst.h"
#include "telspec.h"
#include "visplot.h"
#include "vlbutil.h"
#include "vlbmath.h"
#include "scans.h"
#include "cpgplot.h"
#include "logio.h"
#include "freelist.h"

/*
 * Define a structure to use to store the details of the scans
 * in a given sub-array.
 */
typedef struct Scan Scan;
struct Scan {
  Scan *next;        /* The next scan in the list */
  float vxa,vxb;     /* Min/max NDC X-coords of scan sub-plot */
  float stmin,stmax; /* The time range in the scan */
  float tmin,tmax;   /* The visible part of the time range */
  int view;          /* True if any of the scan is visible */
};

#define SCAN_BLK_SIZE 100  /* The number of scans per freelist block */

/*
 * Define a structure to be used to flag the status of a given
 * station in a given integration.
 */
typedef struct {
  int used;    /* Count of sampling by good visibilities */
  int flagged; /* Count of sampling by flagged visibilities */
  int badcor;  /* Count of sampling by correction-flagged visibilities */
} Telstat;

/* Declare a type used to store details about the last cursor input */

typedef struct {
  int key;        /* The upper case version of the key that the user pressed */
  int waslow;     /* True if the key the user pressed was lower case */
  float tval;     /* The time value selected with the cursor */
  float yval;     /* The y-axis value selected with the cursor */
  Scan *scan;     /* The scan containing tval */
} Tcurs;

/*
 * Elements of the following type are used to associate sorted times,
 * in the currently selected time system, with integrations.
 */
typedef struct {
  Integration *integ; /* The integration of this time sample */
  float t;            /* The time in the form used in plotting the X-axis */
} TimeSample;

static int t_cmp_time_samples(const void *v1, const void *v2);

/*
 * Define a structure used to encapsulate state-info regading the
 * current plot.
 */
typedef struct {
  Subspec ss;        /* The specification of the sub-array being plotted */
  double utref;      /* Reference UT of observation (seconds) */
  double stref;      /* The apparent sidereal time at utref (seconds) */
  Observation *ob;   /* The descriptor of the observation being plotted */
  Subarray *sub;     /* The descriptor of sub-array ss->isub in ob->sub[] */
  Telstat *ts;       /* Used to record the status of each telescope */
  TimeSample *times; /* The time samples in X-axis plot order */
  int times_stale;   /* True if the times in times[] or the scans in 'scans' */
                     /*  are out of date. */
  FreeList *scan_mem;/* Memory for allocating Scan elements */
  Scan *scans;       /* List of scans */
  float wxa,wxb;     /* The plotted range of X-axis world coordinates */
  float wya,wyb;     /* The plotted range of Y-axis world coordinates */
  float vxa,vxb;     /* The X-axis NDC coords of the surrounding viewport */
  float vya,vyb;     /* The Y-axis NDC coords of the surrounding viewport */
  int ta,tb;         /* The indexes of the first and last plotted */
                     /*  integrations in times[]. */
  int docurs;        /* True when cursor control is in effect */
  int dobig;         /* If true, plot with larger dot size */
  int doscan;        /* True if the plot is to be separated into scans */
  int docross;       /* True to enable cross-hair mode */
  int doutc;         /* True to plot UTC along the X-axis, false for sidereal */
                     /*  time. */
  Tcurs cursor;      /* Cursor entry descriptor */
  int modified;      /* 0 - until first edit operation */
  int if_ed;         /* If true, edits are restricted to the current IF */
  int ch_ed;         /* If true, edits are restricted to current freq */
                     /*  channels. */
  int npage;         /* The sequential number of the page being plotted */
} Tpar;

/* Set keys for interactive display-editing */

enum {
  KEY_NONE='\0',  /* Null key press */
  KEY_DOT ='.',   /* Key to toggle marker symbol size */
  KEY_CUR ='A',   /* Cursor position select. */
  KEY_BRK ='B',   /* Toggle separating plot into scans */
  KEY_CUT ='C',   /* Flag cursor selected region */
  KEY_CAN ='D',   /* Key to cancel time range */
  KEY_HELP='H',   /* Key to request help */
  KEY_IF  ='I',   /* Toggle IF editing mode */
  KEY_DIS ='L',   /* Redisplay the plot */
  KEY_NEXT='N',   /* Next sub-array */
  KEY_PREV='P',   /* Prev sub-array */
  KEY_REST='R',   /* Unflag cursor selected region */
  KEY_SUB ='T',   /* Select a sub-array by number */
  KEY_UT  ='U',   /* Select time range via the cursor */
  KEY_CH  ='W',   /* Toggle channel editing mode */
  KEY_QUIT='X',   /* Quit */
  KEY_PRVIF='[',  /* Show the previous IF */
  KEY_NXTIF=']',  /* Show the next IF */
  KEY_CROSS='+',  /* Toggle cross-hair cursor mode */
  KEY_GST  ='G'   /* Toggle between UTC and GST */
};

static const float xmarg=0.05; /* The fraction of the X range for margin */
static const int datcol=10;    /* The color of good data points */
static const int badcol=2;     /* The color of flagged data points */
static const int badccol=11;   /* The color of correction flagged data points */
static const int parcol=7;     /* The color for partially flagged data */
static const int cutcol=2;     /* PGPLOT color index for cut edit window */
static const int rescol=10;    /* PGPLOT color index for restore edit window */
static const int zoomcol=5;    /* PGPLOT color index for zoom cursor window */
static const int dotsym= -1;   /* Default - small plot marker symbol */
static const int bigsym=1;     /* Alternate - big plot marker symbol */

/* Internal plotting functions */

static Tpar *new_Tpar(Observation *ob, Subspec *ss, int cif, int docurs, int doscan, int dobig);
static Tpar *del_Tpar(Tpar *tp);
static int t_pldata(Tpar *tp, int ta, int tb, int erase);
static int t_label(Tpar *tp);
static float t_time(Tpar *tp, double ut);
static int t_redisp(Tpar *tp);
static int t_yrange(Tpar *tp);
static int t_time_range(Tpar *tp);
static int t_vpwin(Tpar *tp);
static int t_plaxes(Tpar *tp, int erase);
static int t_new_time_range(Tpar *tp);
static int t_flags(Tpar *tp, char key, int waslow);
static int t_edbox(Tpar *tp, int flag);
static int t_newmode(Tpar *tp, int if_ed, int ch_ed);
static int t_mlab(Tpar *tp, int erase);
static int t_sampling(Tpar *tp, Subarray *sub, Integration *integ);
static int t_toggle_timesys(Tpar *tp);

typedef enum {T_ALLNEW, T_SKIP_SUB, T_NXT_SUB} Subop;

static int t_newsub(Tpar *tp, Subop oper, int forward, int report,
		    Subspec *init);

static Scan *new_Scans(Tpar *tp);
static Scan *del_Scans(Tpar *tp);
static int t_get_times(Tpar *tp);
static int t_update_times(Tpar *tp, int force);

/* Describe possible cursor band types */

typedef enum {
  B_NORM=0,   /* Normal cursor - no banding required */
  B_LINE=1,   /* Line drawn between reference position and cursor */
  B_RECT=2,   /* Rectangle drawn between reference position and cursor */
  B_YRNG=3,   /* Two horizontal lines bracketing a Y-axis range */
  B_XRNG=4,   /* Two vertical lines bracketing an X-axis range */
  B_YVAL=5,   /* Vertical line through the cursor */
  B_XVAL=6,   /* Horizontal line through the cursor */
  B_CROSS=7   /* Cross-hair */
} Bandmode;

static int t_cursor(Tpar *tp, int noout, Bandmode mode, float xref, float yref,
		    int ci);

/*.......................................................................
 * Create a new Tpar descriptor and an array of scans, then assign
 * given defaults.
 *
 * Input:
 *  ob Observation *  The data descriptor.
 *  ss     Subspec *  The specification of the first sub-array, or
 *                    NULL for the default.
 *  cif        int    The index of the first IF to plot, or -1 for the
 *                    first unsampled IF.
 *  docurs     int    True if cursor control is required. If the current
 *                    device has no cursor, this will be ignored.
 *  dobig      int    True to plot with a larger dot size.
 * Output:
 *  return    Tpar *  The new Tpar descriptor, or NULL on error.
 */
static Tpar *new_Tpar(Observation *ob, Subspec *ss, int cif, int docurs,
		      int doscan, int dobig)
{
  Tpar *tp;      /* Pointer to the new Tpar descriptor */
  char answer[10]; /* Answer from pgqinf() */
  int slen;        /* Length of string */
/*
 * If no sub-array specification was provided, substitute a default.
 */
  if((ss && next_sub(ob, FIND_FIRST, 1, ss->nfix, 0, 1, ss)) ||
     (!ss && !(ss=find_sub(ob, 0, 0, 1, 0, 0, 1))))
    return NULL;
/*
 * An IF index of -1 (0 on the command line) requests the default IF,
 * substitute the first unsampled IF.
 */
  if(cif == -1) {
    if((cif = nextIF(ob, 0, 1, 1)) < 0) {
      lprintf(stderr, "tplot: There are no selected IFs available.\n");
      return NULL;
    };
  } else if(cif < 0 || cif >= ob->nif) {
    lprintf(stderr, "tplot: IF %d does not exist.\n", cif+1);
    return NULL;
  };
/*
 * Attempt to read the start IF.
 */
  if(getIF(ob, cif))
    return NULL;
/*
 * Allocate the Tpar descriptor.
 */
  tp = (Tpar *) malloc(sizeof(Tpar));
  if(tp == NULL) {
    lprintf(stderr, "new_Tpar: Insufficient memory for plot descriptor.\n");
    return tp;
  };
/*
 * NULLify pointer members so that del_Tpar can be called before the
 * struct has been fully initialized.
 */
  tp->sub = NULL;
  tp->ts = NULL;
  tp->times = NULL;
  tp->scan_mem = NULL;
  tp->scans = NULL;
/*
 * Record the descriptor of the observation.
 */
  tp->ob = ob;
/*
 * Record the initial telescope specification.
 */
  tp->ss = *ss;
/*
 * Get the descriptor of the initial sub-array.
 */
  tp->sub = ob->sub + ss->isub;
/*
 * We need a work array with a size equal to the maxmimum number of stations
 * in any sub-array. Find the max number of stations.
 */
  {
    int maxstat = 0;
    int isub;
    for(isub=0; isub<ob->nsub; isub++) {
      Subarray *sub = ob->sub + isub;
      if(sub->nstat > maxstat)
	maxstat = sub->nstat;
    };
/*
 * Allocate an integer work array of this size.
 */
    tp->ts = malloc(sizeof(Telstat) * maxstat);
    if(tp->ts==NULL) {
      lprintf(stderr, "new_Tpar: Insufficient memory.\n");
      return del_Tpar(tp);
    };
  };
/*
 * We need a time-sample array with a size equal to the maximum number
 * of times in any sub-array.
 */
  {
    int maxtime = 0;
    int isub;
/*
 * Find out the max number of times needed in any subarray.
 */
    for(isub=0; isub<ob->nsub; isub++) {
      Subarray *sub = ob->sub + isub;
      if(sub->ntime > maxtime)
	maxtime = sub->ntime;
    };
/*
 * Allocate a time-sample array.
 */
    tp->times = malloc(sizeof(TimeSample) * maxtime);
    if(!tp->times) {
      lprintf(stderr, "new_Tpar: Insufficient memory.\n");
      return del_Tpar(tp);
    };
/*
 * Mark the times as needing updating.
 */
    tp->times_stale = 1;
  };
/*
 * Allocate a freelist for Scan elements.
 */
  tp->scan_mem = new_FreeList("new_Tpar", sizeof(Scan), SCAN_BLK_SIZE);
  if(!tp->scan_mem)
    return del_Tpar(tp);
/*
 * Assign starting values to the rest of the members.
 */
  tp->docross = 0;
  tp->doutc = 1;
  tp->dobig = 0;
  tp->wxa = tp->wxb = 0.0f;
  tp->ta = 0;
  tp->tb = tp->sub->ntime - 1;
  tp->utref = ob->date.ut;
  tp->stref = ob->date.app_st * rtoh * 3600;
/*
 * If cursor interaction is required, check if the device has a cursor.
 */
  if(docurs) {
    slen = sizeof(answer)-1;
    cpgqinf("CURSOR", answer, &slen);
    docurs = strncmp(answer,"YES",3) == 0;
  };
  tp->docurs = docurs;
  tp->cursor.key = KEY_NONE;
  tp->dobig = dobig;
  tp->doscan = doscan;
/*
 * Data not modified yet.
 */
  tp->modified = 0;
/*
 * Set default editing scope.
 */
  tp->if_ed = 0; /* Global IF editing */
  tp->ch_ed = 0; /* Global channel editing */
/*
 * No pages plotted yet.
 */
  tp->npage = 0;
/*
 * Return the new descriptor.
 */
  return tp;
}

/*.......................................................................
 * Tpar (visibility plot descriptor) destructor function.
 *
 * Input:
 *  tp     Tpar *  Tpar pointer returned by new_Tpar().
 * Output:
 *  return Tpar *  Always NULL, so that you can write tp=del_Tpar(tp);
 */
static Tpar *del_Tpar(Tpar *tp)
{
  if(tp != NULL) {
    if(tp->ts)
      free(tp->ts);
    if(tp->times)
      free(tp->times);
    tp->scan_mem = del_FreeList("del_Tpar", tp->scan_mem, 1);
    tp->scans = NULL;   /* NB. Already deleted by deleting scan_mem */
    free(tp);
  };
  return NULL;
}

/*.......................................................................
 * Return a list of scans in tp->scans.
 *
 * Input:
 *  tp       Tpar *  The plot parameter container.
 *                     tp->doscan: If 0 the whole observation is
 *                                 treated as one scan.
 *                     tp->scan_mem:  The freelist of scans.
 *                     tp->scans: This is deleted first.
 *
 * Input/Output:
 *  return   Scan *  The the value of tp->scans, which will be NULL
 *                   on error.
 */
static Scan *new_Scans(Tpar *tp)
{
  Subarray *sub;    /* The descriptor of the current sub-array */
  int ta,tb;        /* Integration indexes of start and end of a scan */
/*
 * Get the current sub-array.
 */
  sub = tp->sub;
/*
 * Start by deleting any existing scan list.
 */
  (void) del_Scans(tp);
/*
 * If the times[] array needs updating do it first.
 */
  if(tp->times_stale && t_get_times(tp))
    return NULL;
/*
 * If not in scan mode, only one scan is needed covering the whole
 * range.
 */
  if(!tp->doscan) {
    Scan *scan = new_FreeListNode("new_Scans", tp->scan_mem);
    if(!scan)
      return NULL;
    scan->next = NULL;
    scan->vxa = tp->vxa;
    scan->vxb = tp->vxb;
    scan->stmin = tp->times[0].t;
    scan->stmax = tp->times[sub->ntime-1].t;
    scan->tmin = scan->tmax = 0.0;
    tp->scans = scan;
  } else {
/*
 * Get the minimum time interval that is taken to separate scans.
 */
    double tsep = sub->scangap <= 0.0 ? DEFGAP : sub->scangap;
/*
 * The tail of the list of scans.
 */
    Scan *tail = NULL;
/*
 * Find the start and end of each scan, and record them in the list.
 */
    for(ta=0; ta < sub->ntime; ta = tb) {
      float prev_time;  /* The previous time seen */
/*
 * Allocate a new scan container.
 */
      Scan *scan = new_FreeListNode("new_Scans", tp->scan_mem);
      if(!scan)
	return del_Scans(tp);
      scan->next = NULL;
      scan->vxa = 0.0;
      scan->vxb = 0.0;
      scan->tmin = scan->tmax = 0.0;
/*
 * Get the start and end of the scan.
 */
      prev_time = scan->stmin = tp->times[ta].t;
      for(tb=ta;
	  tb < sub->ntime && tp->times[tb].t - prev_time < tsep;
	  prev_time=tp->times[tb].t, tb++)
	;
/*
 * Record the end time.
 */
      scan->stmax = prev_time;
/*
 * Append the scan to the list.
 */
      if(!tp->scans)
	tp->scans = scan;
      else
	tail->next = scan;
      tail = scan;
    };
  };
/*
 * Return the new list, to show that we succeeded.
 */
  return tp->scans;
}

/*.......................................................................
 * Delete the list of scans in tp->scans.
 *
 * Input:
 *  tp     Tpar *   The resource object of this module.
 * Output:
 *  return Scan *   The deleted list (always NULL).
 */
static Scan *del_Scans(Tpar *tp)
{
  Scan *next = tp->scans;  /* The next scan in the list */
/*
 * Delete the nodes of the list.
 */
  while(next) {
    Scan *scan = next;
    next = scan->next;
    scan = del_FreeListNode("del_Scans", tp->scan_mem, scan);
  };
  tp->scans = NULL;
  return NULL;
}

/*.......................................................................
 * Create a list of integrations sorted into the order of the x-axis
 * timescale.
 *
 * Input:
 *  tp       Tpar *  The resource object of this module.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int t_get_times(Tpar *tp)
{
  int t;   /* An index into tp->times[] */
/*
 * Get the current subarray.
 */
  Subarray *sub = tp->sub;
/*
 * Copy the pertinent times into the tp->times[] array.
 */
  for(t=0; t<sub->ntime; t++) {
    TimeSample *sample = tp->times + t;
    sample->integ = sub->integ + t;
    sample->t = t_time(tp, sample->integ->ut);
  };
/*
 * For UTC, the times are already sorted. Sidereal times however, need
 * to be sorted.
 */
  if(!tp->doutc)
    qsort(tp->times, sub->ntime, sizeof(tp->times[0]), t_cmp_time_samples);
/*
 * Mark the time array as up to date.
 */
  tp->times_stale = 0;
  return 0;
}

/*.......................................................................
 * This is a qsort() comparison function for sorting arrays of TimeSample
 * elements.
 */
static int t_cmp_time_samples(const void *v1, const void *v2)
{
  const TimeSample *s1 = v1;
  const TimeSample *s2 = v2;
  if(s1->t < s2->t)
    return -1;
  else if(s1->t == s2->t)
    return 0;
  else
    return 1;
}

/*.......................................................................
 * Determine the Y-axis plot range in the current sub-array. Assign the
 * new range in the Tpar descriptor.
 *
 * Input:
 *  tp        Tpar * This contains existing plotting attributes.
 * Output:
 *  return    int   0 - OK.
 *                  On error -1 is returned and no changes are made
 *                  to tp.
 */
static int t_yrange(Tpar *tp)
{
/*
 * Check inputs.
 */
  if(tp==NULL) {
    lprintf(stderr, "t_yrange: NULL plot intercepted\n");
    return -1;
  };
/*
 * The Y-axis plot range for a given sub-array is simply the number
 * of stations in that sub-array, with room above and below.
 */
  tp->wya = -1.0f;
  tp->wyb = (float) tp->sub->nstat;
  return 0;
}

/*.......................................................................
 * Return the time plot range for the time range and plot options in the
 * passed Tpar object.
 *
 * Input/Output:
 *  tp        Tpar * On entry this contains existing plotting attributes.
 *                   Currently only ta and tb need be initialized.
 *                   On output tp->wxa and tp->wxb will contain
 *                   the min and max times of the X-axis.
 * Output:
 *  return    int    0 - OK.
 *                   On error -1 is returned and no changes are made
 *                   to *tmin or *tmax.
 */
static int t_time_range(Tpar *tp)
{
  float xa;    /* Start time of range */
  float xb;    /* End time of range */
  Scan *scan;  /* The scan being processed */
/*
 * Check inputs.
 */
  if(tp==NULL) {
    lprintf(stderr, "t_time_range: NULL plot intercepted\n");
    return -1;
  };
/*
 * Valid ta and tb?
 */
  if(tp->ta < 0 || tp->ta > tp->tb || tp->tb >= tp->sub->ntime) {
    lprintf(stderr, "t_time_range: ta and tb are invalid\n");
    return -1;
  };
/*
 * Make sure that the times[] array and scan list are up to date.
 */
  if(t_update_times(tp, 0))
    return -1;
/*
 * Determine the times corresponding to integrations ta and tb.
 */
  tp->wxa = tp->times[tp->ta].t;
  tp->wxb = tp->times[tp->tb].t;
/*
 * Determine the displayed time ranges within the scans.
 * sc->view flags whether any of the scan is visible.
 */
  for(scan=tp->scans; scan; scan=scan->next) {
    scan->view = tp->wxb >= scan->stmin && tp->wxa <= scan->stmax;
    if(scan->view) {
      xa = (tp->wxa < scan->stmin) ? scan->stmin : tp->wxa;
      xb = (tp->wxb > scan->stmax) ? scan->stmax : tp->wxb;
/*
 * Leave a fractional margin around time range. (Also ensure that the min
 * range is 30 seconds to avoid precision problems).
 */
      if(fabs(xb - xa) > 30.0f) {
	scan->tmin = xa - (xb-xa)*xmarg;
	scan->tmax = xb + (xb-xa)*xmarg;
      } else {
	scan->tmin = xa - 15.0f;
	scan->tmax = xb + 15.0f;
      };
    } else {
      scan->tmin = scan->tmax = 0.0f;
    };
  };
  return 0;
}

/*.......................................................................
 * Set up the viewport limits for the plots leaving 4 char heights on
 * each side of plot for labelling.
 *
 * Input/Output:
 *  tp     Tpar *  The correction parameter struct.
 *                 On output the vxa,vxb,vya,vyb fields will be
 *                 initialized. All other fields are ignored.
 * Output:
 *  return  int     0 - OK.
 */
static int t_vpwin(Tpar *tp)
{
  float vxa,vxb;  /* X viewport limits enclosing whole plot */
  float vya,vyb;  /* Y viewport limits enclosing whole plot */
  float tsum;     /* Sum of scan times within current time range */
  Scan *scan;     /* The scan being processed */
/*
 * Get the standard viewport.
 */
  cpgsch(1.0f);
  cpgvstd();
  cpgqvp(0, &vxa, &vxb, &vya, &vyb);
/*
 * Store this in the plot descriptor, for subsequent use during labelling.
 */
  tp->vxa = vxa;
  tp->vxb = vxb;
  tp->vya = vya;
  tp->vyb = vyb;
/*
 * Make sure that the times[] array and scan list are up to date.
 */
  if(t_update_times(tp, 0))
    return 1;
/*
 * Apportion viewports horizontally for different scans.
 * First find the sum of time ranges covered by all scans within the
 * current time range.
 */
  tsum = 0.0f;
  for(scan=tp->scans; scan; scan=scan->next)
    tsum += scan->tmax - scan->tmin;
/*
 * Use the fraction of the sum of time ranges taken up by each scan
 * to determine the fraction of the horizontal viewport range taken
 * up by that scan.
 */
  vxa = tp->vxa;
  for(scan=tp->scans; scan; scan=scan->next) {
    scan->vxa = vxa;
    if(scan->view)
      scan->vxb = vxa + (tp->vxb - tp->vxa) * (scan->tmax - scan->tmin) / tsum;
    else
      scan->vxb = scan->vxa;    /* Scan not visible */
    vxa = scan->vxb;
  };
  return 0;
}

/*.......................................................................
 * Draw plot axes.
 *
 * Input:
 *  tp      Tpar *  The plot descriptor.
 *  erase    int    If true erase current axes instead of plotting.
 * Output:
 *  return   int    0 - OK.
 */
static int t_plaxes(Tpar *tp, int erase)
{
  Subarray *sub; /* The descriptor of the sub-array to be plotted */
  char label[81];/* Labelling string */
  Scan *scan;    /* The scan being labelled */
  float tmin;    /* Start time */
  float tmax;    /* End time */
  float ch;      /* Character height to use */
  int oldcol;    /* Color index on entry to function */
  int itel;      /* Index of the telescope being labeled */
/*
 * Check arguments.
 */
  if(tp==NULL) {
    lprintf(stderr, "t_plaxes: NULL plot descriptor intercepted\n");
    return -1;
  };
/*
 * Make sure that the times[] array and scan list are up to date.
 */
  if(t_update_times(tp, 0))
    return 1;
/*
 * Store current color index, to be restored on return.
 */
  cpgqci(&oldcol);
/*
 * Buffer while plotting the axes.
 */
  cpgbbuf();
/*
 * Set color for drawing or erasing.
 */
  cpgsci(erase ? 0:1);
/*
 * Set the character height.
 */
  ch = 0.8f;
/*
 * Get the descriptor of the sub-array being anottated.
 */
  sub = tp->sub;
/*
 * Plot the two Y-axes at each end of the frame enclosing the scans.
 */
  cpgsch(ch);
/*
 * Y-axes.
 */
  cpgsvp(tp->vxa, tp->vxb, tp->vya, tp->vyb);
  cpgswin(0.0f, 1.0f, tp->wya, tp->wyb);
  cpgbox(" ", 0.0f, 0, "BCT", 1.0f, 0);
/*
 * Label each telescope tick mark on the Y axis.
 */
  cpgsch(0.6f);
  for(itel=0; itel<sub->nstat; itel++) {
    sprintf(label, "%s\\(0699)", sub->tel[itel].name);
    cpgptxt(0.0f, (float) itel, 0.0f, 1.0f, label);
  };
  cpgsch(ch);
/*
 * Do internal and external X-axes for each visible scan.
 */
  for(scan=tp->scans; scan; scan=scan->next) {
/*
 * If the scan isn't visible, skip to the next one.
 */
    if(!scan->view)
      continue;
/*
 * Calculate the start and end time in seconds. For UTC,
 * add one day such that days in the year start from 1 rather than 0.
 */
    if(tp->doutc) {
      tmin = tp->utref + scan->tmin + daysec;
      tmax = tp->utref + scan->tmax + daysec;
    } else {
      tmin = scan->tmin;
      tmax = scan->tmax;
    };
/*
 * Draw internal Y-axes as unadorned vertical lines.
 */
    cpgsvp(tp->vxa, tp->vxb, tp->vya, tp->vyb);
    cpgswin(tp->vxa, tp->vxb, tp->vya, tp->vyb);    
    if(scan->next && scan->next->view) {
      cpgmove(scan->vxb, tp->vya);
      cpgdraw(scan->vxb, tp->vyb);
    };
/*
 * Draw the X-axes.
 */
    cpgsvp(scan->vxa, scan->vxb, tp->vya, tp->vyb);
    cpgswin(tmin, tmax, 0.0f, 1.0f);
    cpgtbox("ZHBCNST", 0.0f, 0, " ", 0.0f, 0);
  };
/*
 * Restore entry color.
 */
  cpgsci(oldcol);
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Plot or erase sampling data between given integrations.
 *
 * Input:
 *  tp       Tpar * The plot descriptor.
 *  ta        int   Index of first time to be plotted from vp->times[].
 *  tb        int   Index of last time to be plotted from vp->times[].
 *  erase     int   If true erase points instead of plotting them.
 * Output:
 *  return    int   0 - OK.
 */
static int t_pldata(Tpar *tp, int ta, int tb, int erase)
{
  Subarray *sub;      /* Local pointer to the sub-array being displayed */
  Scan *scan;         /* The scan being plotted */
  Telstat *ts;        /* Pointer into tp->ts[] */
  int itel;           /* The index of a telescope */
  int oldcol;         /* Color index on entry to function */
  int t;              /* Index within tp->times[] of the integration that is */
                      /*  being plotted. */
  int first;          /* True until the first point has been plotted */
  int marker;         /* Marker symbol */
/*
 * Start pgplot buffering.
 */
  cpgbbuf();
/*
 * Store current color index, to be restored on return.
 */
  cpgqci(&oldcol);
  cpgsch(1.0f);
/*
 * Get the desriptor of the current sub-array.
 */
  sub = tp->sub;
/*
 * Determine which marker symbol to use.
 */
  marker = tp->dobig ? bigsym : dotsym;
/*
 * Make sure that the times[] array and scan list are up to date.
 */
  if(t_update_times(tp, 0))
    return 1;
/*
 * Draw each point with the appropriate symbol and color for its
 * correction status.
 */
  scan = tp->scans;
  first=1;
  for(t=ta; t <= tb; t++) {
    TimeSample *sample = tp->times + t;
    float tval = sample->t;
/*
 * Skip to the right scan for this integration.
 */
    if(first || tval > scan->stmax) {
      first = 0;
      while(tval > scan->stmax)
	scan = scan->next;
      cpgsvp(scan->vxa, scan->vxb, tp->vya, tp->vyb);
      cpgswin(scan->tmin, scan->tmax, tp->wya, tp->wyb);
    };
/*
 * Deposit telescope sampling statistics in tp->ts[0..sub->nstat-1].
 */
    t_sampling(tp, sub, sample->integ);
/*
 * Plot a red point for fully flagged telescope, a green point for
 * a totally unflagged but used telescope, and orange for a partially
 * flagged telescope.
 */
    ts = tp->ts;
    for(itel=0; itel<sub->nstat; itel++,ts++) {
      if(ts->used || ts->flagged || ts->badcor) {
	float ypos = (float) itel;
	int icol;
	if(erase)
	  icol = 0;
	else if(ts->used)
	  icol = ts->flagged || ts->badcor ? parcol : datcol;
	else if(ts->flagged)
	  icol = badcol;
	else
	  icol = badccol;
	cpgsci(icol);
	cpgpt(1, &tval, &ypos, marker);
      };
    };
  };    /* End of  loop over integrations */
/*
 * Restore entry color and terminate pgplot buffering.
 */
  cpgsci(oldcol);
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Read the cursor, and record its position etc.. in tp->cursor.
 *
 * Input:
 *  tp       Tpar *  The plot descriptor.
 *  noout     int    If true then don't return until the cursor is
 *                   pressed inside a sub-plot.
 *  mode Bandmode    The desired type of cursor, from:
 *                    B_NORM  -  A single point is required - no banding.
 *                    B_LINE  -  Line band between xref,yref and the cursor.
 *                    B_RECT  -  Rectangular band between xref,yref and the
 *                               cursor.
 *                    B_YRNG - Two horizontal lines bracketing a Y-axis range.
 *                    B_XRNG - Two vertical lines bracketing an X-axis range.
 *                    B_YVAL - Vertical line through the cursor.
 *                    B_XVAL - Horizontal line through the cursor.
 *                    B_CROSS- Cross hair centered on cursor.
 *  xref    float    xref and yref denote the reference position for
 *  yref    float    band-type cursors. They are ignored if mode==B_NORM.
 *  ci        int    The color index to draw band cursors with.
 * Output:
 *  return    int    0 - OK.
 */
static int t_cursor(Tpar *tp, int noout, Bandmode mode, float xref, float yref,
		    int ci)
{
  static float xpos; /* The X world coordinate of the cursor */
  static float ypos; /* The Y world coordinate of the cursor */
  Tcurs *tc;         /* Pointer to tp->cursor cursor descriptor */
  Scan *scan;        /* The scan being processed */
/*
 * Get the cursor descriptor.
 */
  tc = &tp->cursor;  
/*
 * Set the viewport around the whole viewsurface and make the world
 * coords the same as NDC so that the returned cursor position
 * is measured in NDC.
 */
  cpgsvp(0.0f,1.0f,0.0f,1.0f);
  cpgswin(0.0f,1.0f,0.0f,1.0f);
/*
 * Make sure that the times[] array and scan list are up to date.
 */
  if(t_update_times(tp, 0))
    return 1;
/*
 * If this is the first call this plot session initialize the position
 * at which to bring up the cursor. Otherwise use the values retained from
 * the previous call in xpos and ypos.
 */
  if(tc->key == KEY_NONE) {
    xpos = 0.5f;
    ypos = 0.5f;
  };
/*
 * Substitute the cross-hair cursor for the normal cursor is requested.
 */
  if(mode == B_NORM && tp->docross)
    mode = B_CROSS;
/*
 * Initialize the return value.
 */
  tc->key = KEY_NONE;
  tc->waslow = 0;
  tc->tval = 0.0f;
  tc->yval = 0.0f;
  tc->scan = NULL;
/*
 * Convert the cursor reference positions into NDC.
 */
  switch(mode) {
  case B_RECT: case B_XRNG: case B_YRNG:
    {
      Scan *tail = tp->scans;  /* The tail of the list of scans */
/*
 * Locate the scan that contains the reference time.
 */
      for(scan=tp->scans; scan; scan=scan->next) {
	if(xref >= scan->tmin && xref <= scan->tmax)
	  break;
	tail = scan;
      };
      if(!scan)
	scan = xref < tp->scans->tmin ? tp->scans : tail;
/*
 * Convert the reference UT and station index to the equivalent NDC position.
 */
      xref = scan->vxa + (xref - scan->tmin) * (scan->vxb - scan->vxa) /
	(scan->tmax - scan->tmin);
      yref = tp->vya + (yref - tp->wya) * (tp->vyb - tp->vya) /
	(tp->wyb - tp->wya);
    };
    break;
  default:
    xref = yref = 0.0f;
    break;
  };
/*
 * Read the cursor.
 */
  do {
    char key;
    cpgsci(ci);
    if(!cpgband((int) mode, 0, xref, yref, &xpos, &ypos, &key))
      return 1;
/*
 * Convert key to upper case and record it for return.
 */
    tc->waslow = islower((int)key);
    tc->key = tc->waslow ? toupper((int)key) : key;
/*
 * See if the point is in the plot.
 */
    if((xpos >= tp->vxa && xpos <= tp->vxb) && 
       (ypos >= tp->vya && ypos <= tp->vyb)) {
/*
 * Convert from NDC to the respective Y-value selected.
 */
      tc->yval = tp->wya + (ypos - tp->vya)/(tp->vyb - tp->vya) *
	                   (tp->wyb - tp->wya);
/*
 * Identify the scan that the cursor was in and use this to
 * determine the selected time value.
 */
      for(scan=tp->scans; tc->scan==NULL && scan; scan=scan->next) {
	if(xpos >= scan->vxa && xpos <= scan->vxb) {
	  tc->tval = scan->tmin + (xpos - scan->vxa) /
	    (scan->vxb - scan->vxa) * (scan->tmax - scan->tmin);
	  tc->scan = scan;
	};
      };
    };
    if(!tc->scan && noout)
      printf("The cursor must be in the plots.\n");
  } while(!tc->scan && noout); /* While not in a plot and noout is true */
  return 0;
}

/*.......................................................................
 * Write labels around the frame enclosing the plot.
 *
 * Input:
 *  tp       Tpar *  The plot descriptor.
 * Output:
 *  return      int    0 - OK.
 */
static int t_label(Tpar *tp)
{
  Observation *ob;      /* The descriptor of the observation being plotted */
  Subarray *sub;        /* The descriptor of the sub-array being plotted */
  char awrk[81];        /* Work string for labelling */
  char bwrk[81];        /* Work string for labelling */
/*
 * Check arguments.
 */
  if(tp==NULL) {
    lprintf(stderr, "t_label: NULL plot intercepted\n");
    return -1;
  };
/*
 * Get the descriptors of the observation and sub-array being plotted.
 */
  ob = tp->ob;
  sub = tp->sub;
/*
 * Set the viewport around the plot grid.
 */
  cpgsvp(tp->vxa, tp->vxb, tp->vya, tp->vyb);
/*
 * Compose and write main titles.
 */
  cpgsci(1);
  cpgsch(1.0f);
  sprintf(awrk, "%s  %s", ob->source.name,
	  sutdate(ob->date.year, ob->date.ut, bwrk));
  cpgmtxt("T", 1.7f, 0.0f, 0.0f, awrk);
  sprintf(awrk, "Time sampling in IF %d  Sub-array %d  Pol %s",
	  ob->stream.cif+1, tp->ss.isub+1, Stokes_name(ob->stream.pol.type));
  cpgmtxt("T", 0.5f, 0.0f, 0.0f, awrk);
/*
 * In non-interative mode, tell the user what is being plotted.
 */
  if(!tp->docurs)
    lprintf(stdout, "Page %02.2d: Subarray %d\n", tp->npage, tp->ss.isub+1);
/*
 * Write the X-axis label.
 */
  if(tp->doutc)
    strcpy(awrk, "Universal Time");
  else
    strcpy(awrk, "Greenwhich Mean Sidereal Time");
  cpgmtxt("B", 3.0f, 0.5f, 0.5f, awrk);
  return 0;
}

/*.......................................................................
 * Replot the current scans.
 *
 * Input:
 *  tp        Tpar *  The plot descriptor.
 * Output:
 *  return       int    0 - OK.
 */
static int t_redisp(Tpar *tp)
{
  int ierr=0; /* True if an error occurs */
/*
 * Cursory check of arguments.
 */
  if(tp==NULL) {
    lprintf(stderr, "t_redisp: NULL plot intercepted\n");
    return -1;
  };
/*
 * Clear page.
 */
  cpgpage();
/*
 * Count pages.
 */
  tp->npage++;
/*
 * If necessary, update the time array and scan list.
 */
  ierr = ierr || t_update_times(tp, 0);
/*
 * Determine the UT plot range for all plots.
 */
  ierr = ierr || t_time_range(tp);
/*
 * Set up viewport.
 */
  ierr = ierr || t_vpwin(tp);
/*
 * Draw scans.
 */
  cpgbbuf();
  ierr = ierr || t_yrange(tp);
  ierr = ierr || t_plaxes(tp, 0);
  ierr = ierr || t_pldata(tp, tp->ta, tp->tb, 0);
  ierr = ierr || (tp->docurs && t_mlab(tp, 0));
  ierr = ierr || t_label(tp);
  cpgebuf();
  return ierr;
}

/*.......................................................................
 * Provide an interactive means of displaying telescope sampling in
 * successive sub-arrays of the current observation.
 *
 * Input:
 *  ob  Observation *   The observation who's corrections are to be
 *                      displayed.
 *  ss      Subspec *   The specification of the first sub-array, or
 *                      NULL for the default.
 *  cif         int     The index of the first IF to plot.
 *  docurs      int     If true enter interactive cursor mode.
 *  opts        char *  An optional string of flag toggling keys that
 *                      toggle after the values below have been applied.
 *                      Send NULL or empty string if not required.
 *  modified    int *   If data are editted, *modified will be assigned true.
 *                      If no data are editted, *modified will be assigned 0.
 * Output:
 *  return      int     0 - ok.
 *                      1 - error.
 */
int timplt(Observation *ob, Subspec *ss, int cif, int docurs,
	   char *opts, int *modified)
{
  int ierr=0;      /* Error flag */
  int old_if;      /* State of current IF to be restored on exit */
  Tpar *tp;        /* The plot descriptor */
  int i;
/*
 * No data have been edited yet.
 */
  if(modified!=NULL)
    *modified = 0;
/*
 * Check the state of the observation.
 */
  if(!ob_ready(ob, OB_SELECT, "timplt"))
    return 1;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Allocate a plot descriptor.
 */
  tp = new_Tpar(ob, ss, cif, docurs, 0, 0);
  if(tp==NULL)
    return 1;
/*
 * If a string of flag options was given, interpret them here.
 */
  if(opts != NULL) {
    int slen = strlen(opts);
    for(i=0; i<slen; i++) {
      int key = opts[i];
      int waslow = islower(key);
      if(waslow) key = toupper(key);
      t_flags(tp, key, waslow);
/*
 * Handle flags that are normally treated as commands.
 */
      switch(key) {
      case KEY_DOT:
	tp->dobig = !tp->dobig;
	break;
      case KEY_CROSS:
	tp->docross = !tp->docross;
	break;
      case KEY_GST:
	tp->doutc = !tp->doutc;
	tp->times_stale = 1;
	break;
      };
    };
  };
/*
 * Interactive plotting?
 */
  if(tp->docurs) {
/*
 * Inform user of the help option.
 */
    lprintf(stdout,
	    "Move the cursor into the plot window and press \'%c\' for help\n",
	    KEY_HELP);
/*
 * Display the sampling for the specified sub-array.
 */
    ierr = ierr || t_redisp(tp);
/*
 * Start interactive editing loop.
 */
    while(!ierr && tp->cursor.key!=KEY_QUIT) {
      int nflag=0;     /* Number of flag toggling operations done */
/*
 * Read the cursor and obey toggle flags until the first non-toggle
 * key is intercepted.
 */
      while(!(ierr=t_cursor(tp, 0, B_NORM, 0.0f, 0.0f, 1)) &&
	    t_flags(tp, tp->cursor.key, tp->cursor.waslow)==0)
	nflag++;
/*
 * Take action appropriate to the key that the user pressed.
 */
      if(!ierr && nflag > 0) {  /* Update display after a flag toggling */
	ierr = t_redisp(tp);
      } else if(!ierr) {
/*
 * Obey the cursor request.
 */
	switch(tp->cursor.key) {
	case KEY_NEXT:                    /* Plot next page */
	  ierr = t_newsub(tp, T_SKIP_SUB, 1, 1, NULL) < 0;
	  break;
	case KEY_PREV:                    /* Plot previous page */
	  ierr = t_newsub(tp, T_SKIP_SUB, 0, 1, NULL) < 0;
	  break;
	case KEY_CUT:
	  ierr = t_edbox(tp, 1);
	  break;
	case KEY_REST:
	  ierr = t_edbox(tp, 0);
	  break;
	case KEY_SUB:
	  {
	    Subspec *ss = read_Subspec(tp->ob, NULL, NULL, tp->ss.isub);
	    ierr = ss && t_newsub(tp, T_ALLNEW, 1, 1, ss) < 0;
	  };
	  break;
	case KEY_DIS:
	  ierr = ierr || t_redisp(tp);
	  break;
	case KEY_DOT:
	  tp->dobig = !tp->dobig;
	  ierr = ierr || t_redisp(tp);
	  break;
	case KEY_UT:
	  ierr = ierr || t_new_time_range(tp);
	  break;
	case KEY_IF:       /* Toggle IF editing mode */
	  ierr = t_newmode(tp, !tp->if_ed, tp->ch_ed);
	  break;
	case KEY_CH:       /* Toggle channel editing mode */
	  ierr = t_newmode(tp, tp->if_ed, !tp->ch_ed);
	  break;
	case KEY_PRVIF:
	case KEY_NXTIF:
	  {
	    int step = tp->cursor.key==KEY_NXTIF ? 1 : -1;
	    int cif = nextIF(ob, ob->stream.cif + step, 1, step);
	    ierr = cif >= 0 && (getIF(ob, cif) || t_redisp(tp));
	  };
	  break;
	case KEY_CROSS:   /* Toggle cross-hair cursor mode */
	  tp->docross = !tp->docross;
	  break;
	case KEY_GST:
	  ierr = ierr || t_toggle_timesys(tp) || t_redisp(tp);
	  break;
	case KEY_HELP:
	  printf("List of keys to enter via cursor.\n");
	  printf(" %c - Quit this session.\n", KEY_QUIT);
	  printf(" %c - Redisplay current plot.\n", KEY_DIS);
	  printf(" %c - Display the Next sub-array.\n", KEY_NEXT);
	  printf(" %c - Display the Previous sub-array.\n", KEY_PREV);
	  printf(" %c - Display the Next IF.\n", KEY_NXTIF);
	  printf(" %c - Display the Previous IF.\n", KEY_PRVIF);
	  printf(" %c - Select a sub-array from the keyboard.\n", KEY_SUB);
	  printf(" %c - Select new UT range with cursor key %c.\n", KEY_UT,
		 KEY_CUR);
          printf(" %c - Initiate selection of an area to flag.\n", KEY_CUT);
          printf(" %c - Initiate selection of an area to un-flag.\n", KEY_REST);
	  printf(" %c - Toggle breaking of display into scans.\n", KEY_BRK);
	  printf(" %c - Toggle IF editing scope.\n", KEY_IF);
	  printf(" %c - Toggle spectral-line channel editing scope.\n", KEY_CH);
	  printf(" %c - Toggle whether to use a cross-hair cursor if available.\n", KEY_CROSS);
	  printf(" %c - Toggle between UTC and Greenwhich sidereal time.\n",
		 KEY_GST);
	  break;
	};
      };
    };
  }
/*
 * Non-interactive plotting?
 */
  else if(!ierr) {
/*
 * Plot the first page.
 */
    ierr = t_redisp(tp);
    if(!ierr) {
      int iret=0;
/*
 * If requested, plot the rest of the available pages.
 */
      while((iret = t_newsub(tp, T_NXT_SUB, 1, 0, NULL)) == 0)
	;
      ierr = iret < 0;
    };
  };
/*
 * If the data were edited, return this information.
 */
  if(tp->modified && modified!=NULL)
    *modified = 1;
/*
 * Flush pending edits.
 */
  ed_flush(tp->ob);
/*
 * Delete the plot descriptor before returning.
 */
  tp = del_Tpar(tp);
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    ierr = 1;
  return ierr;
}

/*.......................................................................
 * Receive input of new time range via the cursor and redisplay the plot
 * within that range. If the user presses the KEY_UT key then display
 * the full time range.
 *
 * Input:
 *  tp      Tpar *  The plot descriptor.
 * Output:
 *  return   int    0 - OK. Anything else means fatal error.
 */
static int t_new_time_range(Tpar *tp)
{
  int dofull=0;   /* True if full time range requested */
  int iter;       /* Iterate over getting two valid keypresses */
  float tval[2]={0.0f,0.0f}; /* The two time end points */
/*
 * Make sure that the times[] array and scan list are up to date.
 */
  if(t_update_times(tp, 0))
    return 1;
/*
 * Inform user.
 */
  printf("For help selecting a new time display range press 'H'.\n");
/*
 * Get the first cursor position for the new time range.
 */
  for(iter=0; iter<2 && !dofull; iter++) {
    int accepted = 0;
    while(!accepted) {
      t_cursor(tp, 1, iter==0 ? B_XVAL : B_XRNG, tval[0], 0.0f, zoomcol);
      switch(tp->cursor.key) {
      case KEY_UT:        /* Display full range and end time selection */
	accepted = dofull = 1;
	break;
      case KEY_QUIT: case KEY_CAN:     /* Abort time selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start time */
	tval[iter] = tp->cursor.tval;
	accepted=1;
	break;
      default:                  /* Unexpected cursor input key - show usage */
	printf("To select a new time display range use keys:\n");
	printf(" %c - Select the %s time.\n", KEY_CUR, iter==0 ? "start":"end");
	printf(" %c - Cancel time display range selection.\n", KEY_CAN);
	printf(" %c - Display the full time display range available.\n", KEY_UT);
	break;
      };
    };
  };
/*
 * Get the time indexes.
 */
  if(dofull) {
    tp->ta = 0;
    tp->tb = tp->sub->ntime - 1;
  } else {
    int t;                   /* tp->times[] index loop variable */
    double tmin = tval[0];
    double tmax = tval[1];
/*
 * Swap the range limits such that tmin < tmax.
 */
    if(tmin > tmax) {double dtmp = tmin; tmin = tmax; tmax = dtmp;};
/*
 * Locate the equivalent time indexes to the selected time values.
 */
    for(t=tp->ta; t<tp->tb && tp->times[t].t < tmin; t++);
    tp->ta = t;
    for(t=tp->ta; t<=tp->tb && tp->times[t].t <= tmax; t++);
    tp->tb = (tp->ta<t) ? (t-1):tp->ta;
  };
/*
 * Display the result.
 */
  return t_redisp(tp);
}

/*.......................................................................
 * Toggle plotting flags given an option key.
 *
 * Input:
 *  tp      Tpar *  The plot descriptor.
 *  key     char    The command key (upper case).
 *  waslow   int    If true the key was lower case.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Key unrecognized.
 */
static int t_flags(Tpar *tp, char key, int waslow)
{
  switch (key) {
  case KEY_BRK:
    tp->doscan = !tp->doscan;
    tp->times_stale = 1;
    break;
  default:
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Display the sampling for a new sub-array.
 *
 * Input:
 *  tp           Tpar *   Plot-parameter block.
 *  oper        Subop     T_ALLNEW  - Display the sub-array specified in
 *                                    *init.
 *                        T_SKIP_SUB - Start plotting from the next sub-array.
 *                        T_NXT_SUB - Plot the next sub-array if allowed
 *                                    by the current sub-array specification.
 *  forward       int     0 - Locate the next sub-array in order of
 *                            decreasing sub-array index.
 *                        1 - Locate the next sub-array in order of
 *                            increasing sub-array index.
 *  init      Subspec *  The new sub-array specification to be used when
 *                       oper == T_ALLNEW.
 * Output:
 *  return        int     0 - OK.
 *                        1 - No more sub-arrays in specified direction.
 *                       -1 - Error.
 */
static int t_newsub(Tpar *tp, Subop oper, int forward, int report,
		    Subspec *init)
{
  Subspec ss;      /* The new sub-array specification */
/*
 * Handle the specified change in sub-array.
 */
  switch(oper) {
  case T_ALLNEW:
    ss = *init;
    if(next_sub(tp->ob, FIND_FIRST, forward, 0, 0, report, &ss))
      return 0;
    break;
  case T_SKIP_SUB:
    ss = tp->ss;
    if(next_sub(tp->ob, SKIP_SUB, forward, 0, 0, report, &ss))
      return 1;
    break;
  case T_NXT_SUB:
    ss = tp->ss;
    if(next_sub(tp->ob, FIND_NEXT, forward, 0, 0, report, &ss))
      return 1;
    break;
  default:
    lprintf(stderr, "t_newsub: Unrecognised opcode.\n");
    return -1;
  };
/*
 * Record the new sub-array.
 */
  tp->ss = ss;
  tp->sub = tp->ob->sub + ss.isub;
/*
 * Mark the times as needing to be updated.
 */
  tp->times_stale = 1;
/*
 * Assign starting values to the rest of the members.
 */
  tp->ta = 0;
  tp->tb = tp->sub->ntime-1;  /* Default to show all data */
/*
 * Display the new data.
 */
  return t_redisp(tp)==0 ? 0 : -1;
}

/*.......................................................................
 * Edit data within a rectangular box specified by the user.
 *
 * Input:
 *  tp    Tpar *  The plot descriptor.
 *  flag   int    If true, flag data within the box.
 *                If false, un-flag data within the box.
 * Output:
 *  return int    0 - OK.
 */
static int t_edbox(Tpar *tp, int flag)
{
  int iter;                   /* Iteration count to get two valid keypresses */
  float yval[2]={0.0f,0.0f};  /* The two Y-axis limits */
  float tval[2]={0.0f,0.0f};  /* The two UT end points wrt tp->refut. */
  int ta,tb;       /* The indexes of the first and last enclosed integrations */
  int tel_a,tel_b; /* The indexes of the first and last enclosed stations */
/*
 * Get the current sub-array.
 */
  Subarray *sub = tp->sub;
/*
 * Make sure that the times[] array and scan list are up to date.
 */
  if(t_update_times(tp, 0))
    return 1;
/*
 * Get the first cursor position for the new UT range.
 */
  for(iter=0; iter<2; iter++) {
    int accepted = 0;
    while(!accepted) {
      t_cursor(tp, 1, iter==0 ? B_NORM : B_RECT, tval[0], yval[0],
	       flag ? cutcol:rescol);
      switch(tp->cursor.key) {
      case KEY_QUIT: case KEY_CAN:     /* Abort UT selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start UT */
	tval[iter] = tp->cursor.tval;
	yval[iter] = tp->cursor.yval;
	accepted=1;
	break;
      default:                  /* Unexpected cursor input key - show usage */
	printf("To edit a rectangular area use keys:\n");
	printf(" %c - Select the %s corner of the area.\n", KEY_CUR,
	       iter==0 ? "start":"end");
	printf(" %c - Abort the edit.\n", KEY_CAN);
	break;
      };
    };
  };
/*
 * Get the indexes of the first and last integrations selected and
 * record them in ta and tb.
 */
  {
    int t;                /* Time index loop variable */
    double tmin = tval[0];
    double tmax = tval[1];
/*
 * Swap the range limits such that tmin < tmax.
 */
    if(tmin>tmax) {double dtmp=tmin; tmin=tmax; tmax=dtmp;};
/*
 * Locate the equivalent ut indexes to the selected ut values.
 */
    for(t=tp->ta; t<tp->tb && tp->times[t].t < tmin; t++);
    ta = t;
    for(t=ta; t<=tp->tb && tp->times[t].t <= tmax; t++);
    tb = ta<t ? (t-1) : tp->ta;
  };
/*
 * Determine the enclosed stations.
 */
  if(yval[0] > yval[1]) {float ytmp=yval[0]; yval[0]=yval[1]; yval[1]=ytmp;};
  tel_a = ceil(yval[0]);
  tel_b = floor(yval[1]);
/*
 * Limit the telescope indexes to legal values.
 */
  if(tel_a < 0) tel_a = 0;
  if(tel_a >= sub->nstat) tel_a = sub->nstat - 1;
  if(tel_b < 0) tel_b = 0;
  if(tel_b >= sub->nstat) tel_b = sub->nstat - 1;
/*
 * tel_a must be <= tel_b.
 */
  if(tel_b < tel_a) {int ttmp=tel_a; tel_a = tel_b; tel_b=ttmp;};
/*
 * Mark the data as modified.
 */
  tp->modified = 1;
/*
 * Edit data associated with the station(s) over the selected time interval.
 */
  {
    int t;
    for(t=ta; t<=tb; t++) {
      Integration *integ = tp->times[t].integ;
      int tel;
/*
 * Get current sampling to see what needs editing.
 */
      t_sampling(tp, sub, integ);
/*
 * Edit stations that need it.
 */
      for(tel=tel_a; tel<=tel_b; tel++) {
	Telstat *ts = &tp->ts[tel];
	if(flag ? (ts->used || !ts->flagged) : (!ts->used || ts->flagged)) {
	  if(ed_integ(tp->ob, sub, integ - sub->integ, tp->ob->stream.cif,
		      flag, 0, 1, tp->ch_ed, tp->if_ed, tel))
	    return 1;
	};
      };
    };
  };
/*
 * Re-display the edited region.
 */
  return t_pldata(tp, ta, tb, 0);
}

/*.......................................................................
 * Plot an extra mode label for editting sessions.
 *
 * Input:
 *  tp     Tpar *  The plot descriptor.
 *  erase   int    If true, erase existing mode label.
 * Output:
 *  return  int    0 - OK.
 */
static int t_mlab(Tpar *tp, int erase)
{
  Observation *ob;  /* The descriptor of the observation being plotted */
  int oldcol;       /* Temporary storage for entry color index */
  char label[81];   /* Temporary work string to compose mode label in */
/*
 * Get the descriptor of the observation.
 */
  ob = tp->ob;
/*
 * Store the existing plot color.
 */
  cpgqci(&oldcol);
/*
 * Set the viewport around the plot.
 */
  cpgsvp(tp->vxa, tp->vxb, tp->vya, tp->vyb);
/*
 * Set color to erase or draw.
 */
  cpgsci(erase ? 0:1);
/*
 * Compose the mode label.
 */
  sprintf(label, "Edit %s channels of %s.",
	  tp->ch_ed ? "selected" : "all",
	  tp->if_ed ? "the displayed IF" : "all IFs");
/*
 * Plot mode line.
 */
  cpgsch(1.0f);
  cpgmtxt("T", 2.9f, 0.0f, 0.0f, label);
/*
 * Restore entry color.
 */
  cpgsci(oldcol);
  return 0;
}

/*.......................................................................
 * Re-plot the mode line to reflect changes in edit mode.
 *
 * Input:
 *  tp       Tpar *  The plot descriptor.
 *  if_ed     int    Select IF based editing if true.
 *  ch_ed     int    Select channel based editing if true.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int t_newmode(Tpar *tp, int if_ed, int ch_ed)
{
/*
 * Buffer until the new text has been plotted.
 */
  cpgbbuf();
/*
 * Erase the existing mode line.
 */
  t_mlab(tp, 1);
/*
 * Install the new editing modes.
 */
  tp->if_ed = if_ed;
  tp->ch_ed = ch_ed;
/*
 * Draw the new mode line.
 */
  t_mlab(tp, 0); /* Plot new mode line */
/*
 * reveal the changes.
 */
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Determine the telescope sampling of a given integration.
 * The sampling will be returned in tp->ts[0..sub->nstat-1].
 *
 * Input:
 *  tp           Tpar *  The plot descriptor.
 *  sub      Subarray *  The descriptor of the sub-array to look at.
 *  integ Integration *  The integration to be processed.
 * Output:
 *  return     int    0 - OK.
 */
static int t_sampling(Tpar *tp, Subarray *sub, Integration *integ)
{
  Visibility *vis;    /* Pointer into integ->vis[] */
  Baseline *bptr;     /* Pointer into sub->base[] */
  int base;           /* Index of the baseline in *bptr */
  int itel;           /* The index of a telescope */
  Telstat *ts;        /* Pointer into tp->ts[] */
/*
 * Clear the flag array used to record which stations are sampled.
 */
  ts = tp->ts;
  for(itel=0; itel<sub->nstat; itel++,ts++)
    ts->badcor = ts->flagged = ts->used = 0;
/*
 * Look at the visibility on each baseline of the current integration
 * and add its status to the status of the constituent telescopes.
 */
  ts = tp->ts;
  bptr = sub->base;
  vis = integ->vis;
  for(base=0; base<sub->nbase; base++, vis++, bptr++) {
    if(!vis->bad) {
      ts[bptr->tel_a].used++;
      ts[bptr->tel_b].used++;
    } else if(!(vis->bad & FLAG_DEL)) {
      if(vis->bad & FLAG_BAD) {
	ts[bptr->tel_a].flagged++;
	ts[bptr->tel_b].flagged++;
      };
      if(vis->bad & FLAG_TA)
	ts[bptr->tel_a].badcor++;
      if(vis->bad & FLAG_TB)
	ts[bptr->tel_b].badcor++;
    };
  };
  return 0;
}

/*.......................................................................
 * Given the UT timestamp of an integration, return the time in the
 * timescale specified by tp->doutc.
 *
 * Input:
 *  tp         Tpar *  The resource object of this module.
 *  ut       double    An integration UTC.
 * Output:
 *  return    float    The time ready for plotting.
 */
static float t_time(Tpar *tp, double ut)
{
  if(tp->doutc)
    return ut - tp->utref;
  else
    return fmod(tp->stref + (ut - tp->utref) * ut_to_mst, daysec);
}

/*.......................................................................
 * If needed, update the array of times, and the list of scans.
 *
 * Input:
 *  tp     Tpar *    The resourced object of the program.
 *  force   int      True to force an update regardless of the
 *                   value of tp->times_stale.
 * Output:
 *  return  int      0 - OK.
 *                   1 - Error.
 */
static int t_update_times(Tpar *tp, int force)
{
/*
 * Is an update needed?
 */
  if(force || tp->times_stale) {
    if(t_get_times(tp) || new_Scans(tp)==NULL)
      return 1;
    tp->times_stale = 0;
  };
  return 0;
}

/*.......................................................................
 * Toggle between displaying UTC and sidereal time along the X-axis.
 *
 * Input:
 *  tp    Tpar *   The resource object of the program.
 * Output:
 *  return int     0 - OK.
 *                 1 - Error.
 */
static int t_toggle_timesys(Tpar *tp)
{
  tp->doutc = !tp->doutc;
  tp->times_stale = 1;
  tp->ta = 0;
  tp->tb = tp->sub->ntime - 1;
  return 0;
}


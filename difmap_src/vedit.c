#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "obs.h"
#include "vlbconst.h"
#include "vlbutil.h"
#include "telspec.h"
#include "visplot.h"
#include "vplot.h"
#include "cpgplot.h"
#include "logio.h"

/*
 * Define the selection keys.
 */
enum {
  KEY_NONE='\0',  /* Null key press */
  KEY_HELP='H',   /* Key to list usage information */
  KEY_UT  ='U',   /* Key to select time input mode */
  KEY_INT =' ',   /* Enter integration based flagging mode */
  KEY_REST='R',   /* Key to introduce boxed-restore */
  KEY_CUT ='C',   /* Key to introduce cutting */
  KEY_FLG ='F',   /* Toggle display of flagged data */
  KEY_ERR ='E',   /* Toggle display of error bars */
  KEY_GST ='G',   /* Toggle between UTC and GST */
  KEY_IF  ='I',   /* Toggle IF editing mode */
  KEY_CH  ='W',   /* Toggle channel editing mode */
  KEY_FUL ='V',   /* Key to select full display range */
  KEY_DIS ='L',   /* Key to redisplay the current plot */
  KEY_NXT ='N',   /* Key to display the next baseline */
  KEY_ORDER='O',  /* Toggle baseline ordering */
  KEY_PRV ='P',   /* Key to display the previous baseline */
  KEY_TEL ='T',   /* Key to request keyboard input reference-telescope */
  KEY_CUR ='A',   /* Key for cursor position input */
  KEY_CAN ='D',   /* Key to cancel incomplete select range */
  KEY_MOD ='M',   /* Toggle model plotting */
  KEY_QUIT='X',   /* Key to quit from this function */
  KEY_NUMB='S',   /* Split into keyboard specified number of sub-plots */
  KEY_BRK ='B',   /* Toggle breaking into scans */
  KEY_ZAP ='K',   /* Zap points within a scan (only one baseline at a time) */
  KEY_ZOOM='Z',   /* Zoom in or out wrt amplitude or phase */
  KEY_AMP ='1',   /* Display only amplitudes */
  KEY_PHS ='2',   /* Display only phases */
  KEY_BOTH='3',   /* Display both amplitudes and phases */
  KEY_PRVIF='[',  /* Show the previous IF */
  KEY_NXTIF=']',  /* Show the next IF */
  KEY_CROSS='+',  /* Toggle cross-hair cursor mode */
  KEY_DIFF='-'    /* Toggle to and from plotting the difference between */
                  /*  the data and the model. */
};

static const int cutcol=2;   /* PGPLOT color index for cut edit window */
static const int rescol=10;  /* PGPLOT color index for restore edit window */
static const int zoomcol=5;  /* PGPLOT color index for zoom cursor window */

static int v_box(Vedpar *vp, int doflag);
static int v_find(Vedpar *vp, Vissub *vs, float tval, float value, int isamp);
static int v_edit(Vedpar *vp, Vissub *vs, int flag, int t);
static int v_toggle(Vedpar *vp, Vissub *vs, float tval, float value,
		    int wasamp);
static int v_new_time_range(Vedpar *vp);
static int v_zoom(Vedpar *vp);
static int v_newnum(Vedpar *vp);
static int v_flags(Vedpar *vp, int key, int waslow);

/*.......................................................................
 * Receive input of new time range via the cursor and redisplay the plot
 * within that range. If the user presses the KEY_UT key then display
 * the full time range.
 *
 * Input:
 *  vp      Vedpar *  The visibility plot attributes to supply to visplt().
 * Output:
 *  return     int    0 - OK. Anything else denotes a fatal error.
 */
static int v_new_time_range(Vedpar *vp)
{
  int accepted;              /* True when cursor entry accepted */
  int dofull=0;              /* True if full time range requested */
  int iter;                  /* Iterate over getting two valid keypresses */
  float tval[2]={0.0f,0.0f}; /* The two time end points */
/*
 * Get the first cursor position for the new time range.
 */
  for(iter=0; iter<2 && !dofull; iter++) {
    do {
      if(v_cursor(vp, 1, iter==0?B_XVAL:B_XRNG, 0, NULL, tval[0], 0.0f,
		  zoomcol))
	return 1;
      accepted = 0;
      switch(vp->cursor.key) {
      case KEY_UT:              /* Revert to the full time range */
	accepted = dofull = 1;
	break;
      case KEY_QUIT: case KEY_CAN:     /* Abort time selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected time */
	accepted=1;
	tval[iter] = vp->cursor.tval;
	break;
      default:                  /* Unexpected cursor input key - show usage */
	printf("\nTime range selection:\n");
	printf(" %c - Select the %s time.\n", KEY_CUR, iter==0 ? "start" : "end");
	printf(" %c - Abort selection.\n", KEY_CAN);
	printf(" %c - Revert to the full time range.\n", KEY_UT);
	break;
      };
    } while(!accepted);
  };
/*
 * Get the time indexes.
 */
  if(dofull) {
    vp->ta = 0;
    vp->tb = vp->sub->ntime - 1;
  } else {
    int t;                /* Time index loop variable */
    double tmin = tval[0];
    double tmax = tval[1];
/*
 * Swap the range limits such that tmin < tmax.
 */
    if(tmin>tmax) {double dtmp=tmin; tmin=tmax; tmax=dtmp;};
/*
 * Locate the equivalent t indexes to the selected time values.
 */
    for(t=vp->ta; t<vp->tb && vp->times[t].t < tmin; t++);
    vp->ta = t;
    for(t=vp->ta; t<=vp->tb && vp->times[t].t <= tmax; t++);
    vp->tb = (vp->ta<t) ? (t-1):vp->ta;
  };
/*
 * Display the result.
 */
  return v_redisp(vp);
}

/*.......................................................................
 * Allow the user to zoom in or out in amplitude or phase.
 *
 * Input:
 *  vp      Vedpar *  The visibility plot attributes to supply to visplt().
 * Output:
 *  return     int    0 - OK. Anything else means fatal error.
 */
static int v_zoom(Vedpar *vp)
{
  int dofull=0;   /* True if full phase range requested */
  int iter;       /* Iterate over getting two valid keypresses */
  float value[2]={0.0f,0.0f}; /* The two selected values */
  float wasamp=0; /* True if the first value was an amplitude */
  Vissub *vs=NULL;/* The descriptor of the first sub-plot selected */
/*
 * Get the two cursor selections.
 */
  for(iter=0; iter<2 && !dofull; iter++) {
    int accepted = 0;
    while(!accepted) {
      if(v_cursor(vp, 1, iter==0?B_YVAL:B_YRNG, wasamp, vs, vp->wxa, value[0],
		  zoomcol))
	return 1;
      switch(vp->cursor.key) {
      case KEY_ZOOM:    /* Display full range and end phase range selection */
	accepted = dofull = 1;
	break;
      case KEY_QUIT: case KEY_CAN:     /* Abort selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected value */
	if(iter==1 && vp->cursor.value == value[0]) {
	  printf("Second value identical to first. Selection aborted.\n");
	  return 0;
	} else if(iter==1 && (!wasamp != !vp->cursor.wasamp ||
			      vp->cursor.vs != vs)) {
	  printf("Second selection in a different sub-plot. Selection aborted.\n");
	  return 0;
	} else {
	  vs = vp->cursor.vs;
	  wasamp = vp->cursor.wasamp;
	  value[iter] = vp->cursor.value;
	  accepted=1;
	};
	break;
      default:                  /* Unexpected cursor input key - show usage */
	printf("\nAmplitude or phase range selection:\n");
	printf(" %c - Select the %s value of the range.\n", KEY_CUR,
	       iter==0 ? "start":"end");
	printf(" %c - Abort selection.\n", KEY_CAN);
	printf(" %c - Revert to the full range.\n", KEY_ZOOM);
	break;
      };
    };
  };
/*
 * Assign the selected range.
 */
  if(dofull) {   /* Full data range? */
    vp->ampmin = vp->ampmax = 0.0f; /* This enables autoscaling */
    vp->phsmin = -pi;
    vp->phsmax = pi;
  } else {
/*
 * Sort the values into ascending order.
 */
    if(value[0] > value[1]) {
      float newval = value[0];
      value[0] = value[1];
      value[1] = newval;
    };
/*
 * Install the new limits.
 */
    if(wasamp) {
      vp->ampmin = value[0];
      vp->ampmax = value[1];
    } else {
      vp->phsmin = value[0];
      vp->phsmax = value[1];
    };
  };
/*
 * Display the result.
 */
  return v_redisp(vp);
}

/*.......................................................................
 * Allow a range box to be selected and optionally either flag all points
 * above and below the box, or restore all points inside the box.
 *
 * Input:
 *  vp      Vedpar *  The visibility plot attributes to supply to visplt().
 *  doflag     int    If true flag data within box. If false unflag
 *                    everything within the box.
 * Output:
 *  return     int    0 - OK. Anything else means fatal error.
 */
static int v_box(Vedpar *vp, int doflag)
{
  int t;              /* Time index loop variable */
  Vissub *vs=NULL;    /* Plot in which cursor was pressed */
  int accepted;       /* True when cursor entry accepted */
  double tmin=0.0;    /* Start of the selected time range */
  double tmax=0.0;    /* End of the selected time range */
  float minval=0.0f;  /* Lower of value[0] and value[1] */
  float maxval=0.0f;  /* higher of value[0] and value[1] */
  float xref = 0.0f;  /* The x-axis reference position of the cursor */
  float yref = 0.0f;  /* The y-axis reference position of the cursor */
  int wasamp=0;       /* True when first selection was an amplitude */
  int iter;           /* Iterate over getting two valid keypresses */
/*
 * Get the first cursor position.
 */
  for(iter=0; iter<2; iter++) {
    do {
      if(v_cursor(vp, 1, iter==0 ? B_NORM:B_RECT, wasamp, vs,
			   xref, yref, doflag ? cutcol:rescol))
	return 1;
      accepted = 0;
      switch(vp->cursor.key) {
      case KEY_QUIT: case KEY_CAN:     /* Abort box selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected vertex */
	if(iter==0) {
	  xref = tmin = tmax = vp->cursor.tval;
	  yref = minval = maxval = vp->cursor.value;
	  vs = vp->cursor.vs;
	  wasamp = vp->cursor.wasamp;
	} else {
	  if(vp->cursor.vs != vs || !wasamp != !vp->cursor.wasamp) {
	    fprintf(stderr, "Select box spans more than one plot.\n");
	    return 0;
	  };
	  if(vp->cursor.tval < tmin)
	    tmin = vp->cursor.tval;
	  else
	    tmax = vp->cursor.tval;
	  if(vp->cursor.value < minval)
	    minval = vp->cursor.value;
	  else
	    maxval = vp->cursor.value;
	};
	accepted=1;
	break;
      default:                  /* Unexpected cursor input key - show usage */
	printf("\nSelect %s box.\n", doflag ? "flagging":"restore");
	printf(" %c - Select %s corner.\n", KEY_CUR,
	       iter==0 ? "first":"opposite");
	printf(" %c - Abort selection.\n", KEY_CAN);
	break;
      };
    } while(!accepted);
  };
/*
 * Buffer PGPLOT instructions while the points are redrawn.
 */
  cpgbbuf();
/*
 * Only flag points in the time range selected.
 */
  for(t=vp->ta;  t<=vp->tb;  t++) {
    TimeSample *sample = vp->times + t;
    float tval = sample->t;
    if(tval >= tmin && tval <= tmax) {  /* Only flag in given time range */
      Visibility *vis = sample->integ->vis + vs->base;
/*
 * Ignore deleted points.
 */
      if(!(vis->bad & FLAG_DEL)) {
	int inside;     /* True if the visibility is inside select-box */
/*
 * Get the appropriate amplitude and phase of the visibility.
 */
	float amp,phs;
	v_data_point(vp, vis, &amp, &phs);
/*
 * Determine whether the current point is inside or outside the
 * select-box.
 */
	if(wasamp)
	  inside = amp >= minval && amp <= maxval;
	else {
	  phs -= twopi * floor(phs/twopi + 0.5);
	  inside = phs >= minval && phs <= maxval;
	};
/*
 * Edit and redisplay the point in all relevant sub-plots.
 */
	if(inside && v_edit(vp, vs, doflag, t))
	  return 1;
      };
    };
  };
  cpgebuf(); /* Plotting complete - end plot buffering */
  return 0;
}

/*.......................................................................
 * Flag all points within the scan and baseline selected by the cursor.
 *
 * Input:
 *  vp      Vedpar *  The visibility plot attributes to supply to visplt().
 *  vs      Vissub *  The sub-plot descriptor returned by v_cursor().
 *                    (If vs==NULL, -1 is returned).
 *  scan      Scan *  The scan returned by v_cursor().
 *                    (If sc==NULL, -1 is returned).
 *  doflag     int    If true flag data within box. If false unflag
 *                    everything within the box.
 * Output:
 *  return     int    0 - OK. Anything else means fatal error.
 */
static int v_zap(Vedpar *vp, Vissub *vs, Scan *scan, int doflag)
{
  int t;          /* Time index loop variable */
  int save_mode;  /* Used to record station/baseline editing mode */
  double tmin;    /* Start of the selected time range */
  double tmax;    /* End of the selected time range */
  int ierr=0;     /* Error status flag */
/*
 * Ignore call if the cursor was not pressed within a plot.
 */
  if(vs==NULL || scan==NULL)
    return 0;
/*
 * Get the time values of the start and end of the scan.
 */
  tmin = scan->tmin;
  tmax = scan->tmax;
/*
 * Buffer PGPLOT instructions while the points are redrawn.
 */
  cpgbbuf();
/*
 * Turn off station editing.
 */
  save_mode = vp->stat_ed;
  vp->stat_ed = 0;
/*
 * Only flag points in the time range selected.
 */
  for(t=vp->ta; !ierr && t<=vp->tb; t++) {
    TimeSample *sample = vp->times + t;
    float tval = sample->t;
    if(tval >= tmin && tval <= tmax) {  /* Only flag in given time range */
      Visibility *vis = sample->integ->vis + vs->base;
/*
 * Ignore deleted points.
 */
      if(vis->bad & FLAG_DEL)
	continue;
/*
 * Edit and redisplay the point.
 */
      ierr = v_edit(vp, vs, doflag, t);
    };
  };
  vp->stat_ed = save_mode; /* Restore station/baseline editing mode. */
  cpgebuf(); /* Plotting complete - end plot buffering */
  return ierr;
}

/*.......................................................................
 * Take a cursor position returned by v_cursor() and locate the index
 * of the closest plotted point.
 *
 * Input:
 *  vp       Vedpar *   The plot descriptor.
 *  vs       Vissub *   The sub-plot descriptor returned by v_cursor().
 *                      (If vs==NULL, -1 is returned).
 *  tval      float     The time returned by v_cursor().
 *  value     float     The amp or phase returned by v_cursor().
 *  isamp       int     The point type returned by v_cursor().
 * Output:
 *  return      int     The vp->time[] index of the nearest point or
 *                      -1 if no displayed data in zone where cursor
 *                      was pressed.
 */
static int v_find(Vedpar *vp, Vissub *vs, float tval, float value, int isamp)
{
  int t;             /* The index of the time sample being checked */
  float xtomm,ytomm; /* Conversion factor between world coords and mm */
  int best_time=0;   /* The time index of the point closest to the cursor */
  float dist;        /* The squared dist (mm) from data point to cursor */
  float mindist=0.0f;/* The min value of 'dist' */
  int first=1;       /* True until end of first iteration of search loop */
/*
 * Cursor pressed outside of any plot?
 */
  if(vs==NULL)
    return -1;
/*
 * Determine conversion factors from world coords to mm.
 */
  if(v_scale(vp, vs, isamp, &xtomm, &ytomm))
    return -1;
/*
 * Locate the nearest point.
 */
  for(t=vp->ta; t<=vp->tb; t++) {
    TimeSample *sample = vp->times + t;
    Visibility *vis = sample->integ->vis + vs->base;
    float xdif = xtomm * (sample->t - tval);
    float ydif;
    float amp,phs;
/*
 * Ignore deleted points. Also ignore flagged data if not displayed..
 */
    if((vis->bad & FLAG_DEL) || (vis->bad && !vp->doflag))
      continue;
/*
 * Get the amplitude and phase of the visibility.
 */
    v_data_point(vp, vis, &amp, &phs);
/*
 * Determine the amplitude or phase difference dependant upon which
 * partition of the plot recieved the cursor press.
 */
    if(isamp)
      ydif = ytomm * (value - amp);
    else {
      phs -= twopi * floor(phs/twopi + 0.5);
      ydif = ytomm * (value - phs);
    };
/*
 * Compare the squared distance from this point with the that of
 * the previous closest point.
 */
    dist = xdif*xdif + ydif*ydif;
    if(first || dist < mindist) {
      first = 0;
      best_time = t;
      mindist = dist;
    };
  };
/*
 * No points in plot!
 */
  if(first)
    return -1;
  return best_time;
}

/*.......................................................................
 * Edit one point (by baseline or station) and redisplay on each
 * respective displayed sub-plot.
 *
 * Input:
 *  vp       Vedpar *  The plot descriptor.
 *  vs       Vissub *  The sole baseline to be editted if vp->stat_ed==0.
 *  flag        int    If true, flag unflagged data. If false, restore
 *                     flagged data.
 *  t           int    The index of the time sample to be editted.
 * Output:
 *  return      int    0 - OK.
 */
static int v_edit(Vedpar *vp, Vissub *vs, int flag, int t)
{
  Observation *ob;  /* The descriptor of the observation being edited */
  Subarray *sub;    /* The descriptor of the displayed sub-array */
  int iplot;        /* Sub-plot index */
  int ierr=0;       /* Error status */
/*
 * Check descriptors.
 */
  if(vp==NULL || (!vp->stat_ed && vs==NULL)) {
    lprintf(stderr, "v_edit: NULL %s descriptor intercepted\n",
	    vp==NULL?"plot":"sub-plot");
    return 1;
  };
/*
 * Get the descriptor of the observation and sub-array being edited.
 */
  ob = vp->ob;
  sub = vp->sub;
/*
 * This function modifies the data.
 */
  vp->modified = 1;
/*
 * Start by erasing the given integration from all relevant sub-plots.
 */
  cpgbbuf();  /* Buffer changes */
  if(vp->stat_ed) {
    for(iplot=0; !ierr && iplot<vp->nplot; iplot++)
      ierr = v_pldata(vp, &vp->vplots[iplot], t, t, 1);
  } else {
    ierr = v_pldata(vp, vs, t, t, 1);
  };
/*
 * Now edit the point on all relevant baselines.
 */
  if(!ierr) {
    ierr = ed_integ(ob, vp->sub, vp->times[t].integ - vp->sub->integ,
		    ob->stream.cif, flag, !vp->stat_ed, vp->stat_ed,
		    vp->ch_ed, vp->if_ed, vp->stat_ed ? vp->bs_beg.ta:vs->base);
  };
/*
 * Re-plot the given integration on all relevant sub-plots.
 */
  if(!ierr) {
    if(vp->stat_ed) {
      for(iplot=0; !ierr && iplot<vp->nplot; iplot++)
	ierr = v_pldata(vp, &vp->vplots[iplot], t, t, 0);
    } else {
      ierr = v_pldata(vp, vs, t, t, 0);
    };
  };  /* End of  If no error occurred while editing points */
/*
 * End pgplot buffering.
 */
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Given the values returned by v_cursor(), toggle the flagged status of
 * the integration closest to the cursor.
 *
 * Input:
 *  vp       Vedpar *  The plot descriptor.
 *  vs       Vissub *  The sub-plot descriptor returned by v_cursor().
 *  tval      float    The time value returned by v_cursor().
 *  value     float    The amp or phase returned by v_cursor().
 *  wasamp      int    Type of 'value', returned by v_cursor().
 * Output:
 *  return      int    0 - OK.
 */
static int v_toggle(Vedpar *vp, Vissub *vs, float tval,
	     float value, int wasamp)
{
  Visibility *vis;  /* The nearest visibility */
  int t;            /* The index of the nearest time sample */
/*
 * Find the nearest integration on the specified baseline.
 */
  t = v_find(vp, vs, tval, value, wasamp);
  if(t < 0)
    return 0; /* Ignore failed locations */
/*
 * Toggle the status of the visibility and propogate effects
 * into other effected baselines.
 */
  vis = vp->times[t].integ->vis + vs->base;
  return v_edit(vp, vs, !(vis->bad & FLAG_BAD), t);
}

/*.......................................................................
 * Request a new value for the number of sub-plots per page.
 *
 * Input:
 *  vp       Vedpar *  The plot descriptor.
 * Output:
 *  return      int    0 - OK.
 */
static int v_newnum(Vedpar *vp)
{
  char awrk[81];  /* Temporary work array */
  char *cptr;     /* Pointer into awrk[] */
  long nrow;       /* New number of plots. */
/*
 * Get the new number from the user.
 */
  printf("Enter the required number of plots per page: ");
  if(fgets(awrk, sizeof(awrk)-1, stdin) == NULL) {
    fprintf(stderr, "Error reading input.\n");
    return 0;
  };
/*
 * Skip leading white-space.
 */
  cptr = awrk;
  while(*cptr && isspace((int)*cptr))
    cptr++;
/*
 * Read the number.
 */
  nrow = *cptr ? strtol(cptr, &cptr, 10) : 0;
/*
 * Skip trailing white-space.
 */
  while(*cptr && isspace((int)*cptr))
    cptr++;
  if(*cptr != '\0' || nrow < 0) {
    fprintf(stderr, "Unexpected input (not a positive integer).\n");
    return 0;
  };
/*
 * Register the requested change in numbers of plot slots.
 */
  v_setnrow(vp, (int) nrow);
/*
 * Plot the new number of plots.
 */
  return v_plot(vp, V_REPLOT, 1, NULL) < 0;
}

/*.......................................................................
 * Cursor control interface to plotting and editing functions.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation to be plotted.
 *  bs     Basespec *  The descriptor of the first baseline to be plotted.
 *                     If NULL, the default is the first baseline in the
 *                     observation.
 *  cif         int    The index of the start IF, or -1 for the default.
 *  nrow        int    Number of sub-plots in Y.
 *  npage       int    The max number of pages to plot if non-interactive.
 *                     A value of 0 means no page limit.
 *  docurs      int    If true allow cursor control if device has a cursor.
 *  opts       char *  An optional string of flag toggling keys that
 *                     toggle after the values below have been applied.
 *                     Send NULL or empty string if not required.
 *  doscan      int    If true break up plot into scans if present.
 *  doamp       int    If true plot amplitude part of sub-plots.
 *  dophs       int    If true plot phase part of sub-plots.
 *  doflag      int    If true plot flagged data in addition to unflagged.
 *  domod       int    If true plot model data in addition to unflagged.
 *  dobars      int    If true plot error bars.
 *  showall     int    If true, take account of flagged data in autoscaling.
 * Input/Output:
 *  modified    int *  If modified!=NULL then *modified will be assigned
 *                     with 0 if no data were edited and 1 if data were
 *                     edited.
 * Output:
 *  return      int    0 - OK.
 */
int vedit(Observation *ob, Basespec *bs, int cif, int nrow, int npage,
	  int docurs, char *opts, int doscan, int doamp, int dophs,
	  int doflag, int domod, int dobars, int showall, int *modified)
{
  Vedpar *vp;     /* Plot descriptor */
  int ierr=0;     /* Error status */
  int finished=0; /* True when user quits */
  int old_if;     /* State of current IF to be restored on exit */
  int i;
/*
 * Data not modified yet.
 */
  if(modified!=NULL)
    *modified = 0;
/*
 * Is the observation ready to be plotted?
 */
  if(!ob_ready(ob, OB_SELECT, "vedit"))
    return 1;
/*
 * If no baseline specification has been provided, override with the
 * default.
 */
  if((bs && next_base(ob, FIND_FIRST, 1, 2, 1, 0, 1, bs)) ||
     (!bs && !(bs=find_base(ob, 0, 0, 0, 0, 1, 2, 1, 0, 1))))
    return 1;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Allocate and initialize the plot descriptor.
 */
  vp = new_Vedpar(ob, cif, docurs, doscan, doamp, dophs, doflag, domod, dobars,
		  showall, nrow); 
  if(vp==NULL)
    return 1;
/*
 * If a string of flag options was given, interpret them here.
 */
  if(opts != NULL) {
    size_t slen = strlen(opts);
    for(i=0; i<slen; i++) {
      int key = opts[i];
      int waslow = islower(key);
      if(waslow)
	key = toupper(key);
      v_flags(vp, key, waslow);
/*
 * Handle flags that are normally treated as commands.
 */
      switch(key) {
      case KEY_INT:
	vp->stat_ed = !vp->stat_ed;   /* Toggle station/baseline flagging */
	break;
      case KEY_IF:
	vp->if_ed = !vp->if_ed;       /* Toggle IF flagging */
	break;
      case KEY_CH:
	vp->ch_ed = !vp->ch_ed;       /* Toggle channel flagging */
	break;
      case KEY_ORDER:
	vp->doall = !vp->doall;
	break;
      case KEY_CROSS:
	vp->docross = !vp->docross;
	break;
      case KEY_GST:
	v_toggle_timesys(vp);
	break;
      };
    };
  };
/*
 * Plot the first page.
 */
  ierr = v_plot(vp, V_ALLNEW, 1, bs) <= 0;
/*
 * Interactive plotting?
 */
  if(vp->docurs) {
/*
 * Interactive mode - Inform user of the way to receive usage information.
 */
    lprintf(stdout,
	  "For help move the cursor into the plot window and press \'%c\'.\n",
	  KEY_HELP);
/*
 * Start the interactive display/editing loop.
 */
    while(!finished && ierr==0) {
      int wasflag=0;  /* True if last key stroke toggled a flag */
      int nflag=0;    /* Number of flag toggling operations done */
/*
 * Read the cursor.
 */
      do {
	ierr = v_cursor(vp, 0, B_NORM, 0, NULL, 0.0f, 0.0f, 1);
/*
 * Toggle flags where appropriate.
 */
	wasflag = v_flags(vp, vp->cursor.key, vp->cursor.waslow) == 0;
	nflag += wasflag;
      } while(wasflag);   /* Don't do anything more if a flag was toggled */
/*
 * Take action appropriate to the key that the user pressed.
 */
      if(nflag > 0) {  /* Update display after a sequence of flag toggling */
	nflag = 0;
	ierr = v_plot(vp, V_REPLOT, 1, NULL) < 0;
      } else {
	switch (vp->cursor.key) {
	case KEY_INT:                   /* Toggle flagging mode */
	  ierr = v_newmode(vp, !vp->stat_ed, vp->if_ed, vp->ch_ed);
	  break;
	case KEY_IF:                    /* Toggle IF editing mode */
	  ierr = v_newmode(vp, vp->stat_ed, !vp->if_ed, vp->ch_ed);
	  break;
	case KEY_CH:                    /* Toggle channel editing mode */
	  ierr = v_newmode(vp, vp->stat_ed, vp->if_ed, !vp->ch_ed);
	  break;
	case KEY_DIS:      /* Re-display current baseline */
	  ierr = v_plot(vp, V_REPLOT, 1, NULL) < 0;
	  break;
	case KEY_NXT:                      /* Plot next page */
	  ierr = v_plot(vp, vp->cursor.waslow ? V_NXT_TB:V_NXTSUB, 1, NULL) < 0;
	  break;
	case KEY_PRV:                     /* Plot previous page */
	  ierr = v_plot(vp, vp->cursor.waslow ? V_NXT_TB:V_NXTSUB, 0, NULL) < 0;
	  break;
	case KEY_PRVIF:
	case KEY_NXTIF:
	  {
	    int step = vp->cursor.key==KEY_NXTIF ? 1 : -1;
	    int cif = nextIF(ob, ob->stream.cif + step, 1, step);
	    ierr = cif >= 0 && (getIF(ob, cif) || v_redisp(vp));
	  };
	  break;
	case KEY_TEL:                     /* New reference telescope */
/*
 * Get the user telescope/baseline request and plot the result.
 */
	  bs = read_Basespec(ob, NULL, NULL, vp->bs_beg.isub);
	  ierr = bs && v_plot(vp, V_ALLNEW, 1, bs) < 0;
	  break;
	case KEY_UT:      /* Following cursor input selects time range */
	  ierr = v_new_time_range(vp);
	  break;
	case KEY_ZOOM:    /* Following cursor input selects amp/phase range */
	  ierr = v_zoom(vp);
	  break;
	case KEY_CUT:     /* Flag points within a select-box */
	  ierr = v_box(vp, 1);
	  break;
	case KEY_REST:    /* Restore points inside a select-box */
	  ierr = v_box(vp, 0);
	  break;
	case KEY_NUMB:    /* Change number of plots per page */
	  ierr = v_newnum(vp);
	  break;
	case KEY_CUR:     /* Toggle flag status of nearest point */
	  ierr = v_toggle(vp, vp->cursor.vs, vp->cursor.tval,
			  vp->cursor.value, vp->cursor.wasamp);
	  break;
	case KEY_ZAP:     /* Zap all points of a given scan/baseline */
	  ierr = v_zap(vp, vp->cursor.vs, vp->cursor.scan, 1);
	  break;
	case KEY_CROSS:
	  vp->docross = !vp->docross;
	  break;
	case KEY_ORDER:
	  vp->doall = !vp->doall;
	  if(v_plot(vp, V_RESET, 1, NULL) <= 0)
	    vp->doall = !vp->doall;
	  break;
	case KEY_GST:
	  ierr = ierr || v_toggle_timesys(vp) || v_redisp(vp);
	  break;
	case KEY_HELP:    /* Print usage info */
	  printf("Vplot key bindings:\n");
	  printf(" %c - List the following key bindings.\n", KEY_HELP);
	  printf(" %c - Exit vplot (right-mouse-button).\n", KEY_QUIT);
	  printf(" %c - Flag or un-flag the visibility nearest the cursor (left-mouse-button).\n", KEY_CUR);
	  printf(" %c - Select a new time range (hit %c again for the full range).\n", KEY_UT, KEY_UT);
	  printf(" %c - Select a new amplitude or phase range (hit %c twice for full range).\n", KEY_ZOOM, KEY_ZOOM);
	  printf(" %c - Flag all data inside a specified rectangular box.\n",
		 KEY_CUT);
	  printf(" %c - Restore data inside a specified rectangular box.\n",
		 KEY_REST);
	  printf(" %c - Flag all visibilities of a selected baseline and scan.\n", KEY_ZAP);
	  printf(" %c - Redisplay the current plot.\n", KEY_DIS);
	  printf(" %c - Display the next set of baselines.\n",
		 tolower(KEY_NXT));
	  printf(" %c - Display the preceding set of baselines.\n",
		 tolower(KEY_PRV));
	  printf(" %c - Display the next sub-array.\n", KEY_NXT);
	  printf(" %c - Display the preceding sub-array.\n", KEY_PRV);
	  printf(" %c - Plot from the next IF.\n", KEY_NXTIF);
	  printf(" %c - Plot from the preceding IF.\n", KEY_PRVIF);
	  printf(" %c - Toggle whether to display model visibilities.\n",
		 KEY_MOD);
	  printf(" %c - Toggle whether to display flagged visibilities.\n",
		 KEY_FLG);
	  printf(" %c - Toggle whether to display error bars.\n", KEY_ERR);
	  printf(" %c - Toggle between GST and UTC times along the X-axis.\n",
		 KEY_GST);
	  printf(" %c - Select the number of sub-plots per page.\n", KEY_NUMB);
	  printf(" %c - Toggle between seeing all or just upper baselines.\n",
		 KEY_ORDER);
	  printf(" %c - Plot only amplitudes.\n", KEY_AMP);
	  printf(" %c - Plot only phases.\n", KEY_PHS);
	  printf(" %c - Plot both amplitudes and phases.\n", KEY_BOTH);
	  printf(" %c - Toggle whether to display residuals.\n", KEY_DIFF);
	  printf(" %c - Toggle whether to break the plot into scans (where present).\n", KEY_BRK);
	  printf(" %c - Toggle whether to use flagged data in autoscaling.\n",
		 KEY_FUL);
	  printf(" %c - Toggle whether to use a cross-hair cursor if available.\n", KEY_CROSS);
	  printf(" %c - Request a new reference telescope/baseline.\n",KEY_TEL);
	  printf(" %c - (SPACE BAR) Toggle station based vs. baseline based editing.\n", KEY_INT);
	  printf(" %c - Toggle IF editing scope.\n", KEY_IF);
	  printf(" %c - Toggle spectral-line channel editing scope.\n", KEY_CH);
	  printf("\n");
	  break;
	case KEY_QUIT:
	  finished = 1;
	  break;
	default:
	  break;
	};
      };
    };
  }
/*
 * Non-interactive plotting?
 */
  else if(!ierr) {
    int page;
/*
 * Plot as many pages as required, and keep the user informed.
 * Note that one page has already been plotted.
 */
    for(page=1; !ierr && (npage<=0 || page<npage); page++) {
      int nplotted = v_plot(vp, V_NEXT, 1, NULL);
      if(nplotted < 0)
	ierr = 1;
      else if(nplotted == 0)
	break;
    };
  };
/*
 * Flush any pending edits.
 */
  ed_flush(ob);
/*
 * Have the data been modified.
 */
  if(modified!=NULL)
    *modified = vp->modified;
/*
 * Clean up.
 */
  vp = del_Vedpar(vp);
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    ierr = 1;
  return ierr;
}

/*.......................................................................
 * Toggle plotting flags given a command key.
 *
 * Input:
 *  vp    Vedpar *  The plot descriptor.
 *  key      int    The command key (upper case).
 *  waslow   int    If true the key was lower case.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Key unrecognized.
 */
static int v_flags(Vedpar *vp, int key, int waslow)
{
  switch (key) {
  case KEY_MOD:      /* Toggle the display-model flag */
    vp->domod = !vp->domod;
    break;
  case KEY_FLG:      /* Toggle the display-flagged-data flag */
    vp->doflag = !vp->doflag;
    break;
  case KEY_ERR:
    vp->dobars = !vp->dobars;
    break;
  case KEY_FUL:      /* Toggle full - restricted plot */
    vp->showall = !vp->showall;
    break;
  case KEY_AMP:
    vp->doamp = 1;
    vp->dophs = 0;
    break;
  case KEY_PHS:
    vp->doamp = 0;
    vp->dophs = 1;
    break;
  case KEY_BOTH:
    vp->doamp = vp->dophs = 1;
    break;
  case KEY_DIFF:
    vp->dodiff = !vp->dodiff;
    break;
  case KEY_BRK:
    vp->doscan = !vp->doscan;
    if(vp->sub)
      v_update_scans(vp);
    break;
  default:
    return 1;
  };
  return 0;
}

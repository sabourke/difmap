#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "obs.h"
#include "vlbconst.h"
#include "vlbutil.h"
#include "telspec.h"
#include "visplot.h"
#include "vlbmath.h"
#include "vplot.h"
#include "scans.h"
#include "cpgplot.h"
#include "logio.h"

static const float ymarg=0.1;  /* The fraction of the Y range for margin */
static const float xmarg=0.05; /* The fraction of the X range for margin */
static const int datcol=10;    /* The color of unflagged data points */
static const int badcol=2;     /* The color of flagged data points */
static const int badccol=11;   /* The color of correction flagged data points */
static const int modcol=5;     /* The color to plot the model with */
static const int datsym=1;     /* The marker of good points */
static const int badsym=2;     /* The marker of flagged points */
static const int badcsym=5;    /* The marker of correction flagged points */

static int v_arange(Vedpar *vp, Vissub *vs);
static int v_time_range(Vedpar *vp);
static int v_vpwin(Vedpar *vp, int nrow, int nplot);
static int v_plaxes(Vedpar *vp, Vissub *vs, int dobot, int dotop, int erase);
static int v_mlab(Vedpar *vp, int erase);

static int v_cmp_time_samples(const void *v1, const void *v2);

static Scan *new_Scans(Vedpar *vp);
static Scan *del_Scans(Vedpar *vp);
static int v_get_times(Vedpar *vp);
static float v_time(Vedpar *vp, double ut);

/*.......................................................................
 * Create a new Vedpar descriptor and contained Vissub array, and
 * assign given defaults.
 *
 * Input:
 *  ob Observation *  The descriptor of the observation.
 *  cif        int    The index of the start IF, or -1 for the default.
 *  docurs     int    True if cursor control is required. If the current
 *                    device has no cursor, this will be ignored.
 *  doamp      int    Default for request for amplitudes to be plotted.
 *  dophs      int    Default for request for phases to be plotted.
 *  doflag     int    Default for request for flagged data to be plotted.
 *  domod      int    Default for request for model to be plotted.
 *  dobars     int    Default for request for error bars to be plotted.
 *  showall    int    Default for request for full range of data+model.
 *  nrow       int    The initial number of sub-plots. If nrow<=0
 *                    the default sub->nstat-1 will be used.
 * Output:
 *  return  Vedpar *  The new Vedpar descriptor, or NULL on error.
 */
Vedpar *new_Vedpar(Observation *ob, int cif, int docurs, int doscan, int doamp,
		   int dophs, int doflag, int domod, int dobars, int showall,
		   int nrow)
{
  Vedpar *vp;     /* Pointer to the new Vedpar descriptor */
  char awrk[81];  /* General work string */
  int slen;       /* Length of string */
/*
 * An IF index of -1 (0 on the command line) requests the default IF,
 * substitute the first unsampled IF.
 */
  if(cif == -1) {
    if((cif = nextIF(ob, 0, 1, 1)) < 0) {
      lprintf(stderr, "vplot: There are no selected IFs available.\n");
      return NULL;
    };
  } else if(cif < 0 || cif >= ob->nif) {
    lprintf(stderr, "vplot: IF %d does not exist.\n", cif+1);
    return NULL;
  };
/*
 * Attempt to read the start IF.
 */
  if(getIF(ob, cif))
    return NULL;
/*
 * Allocate the Vedpar descriptor.
 */
  vp = (Vedpar *) malloc(sizeof(Vedpar));
  if(vp == NULL) {
    lprintf(stderr, "new_Vedpar: Insufficient memory for plot descriptor.\n");
    return vp;
  };
/*
 * NULLify pointer members so that del_Vedpar can be called before the
 * struct has been fully initialized.
 */
  vp->vplots = NULL;
  vp->times = NULL;
  vp->scan_mem = NULL;
  vp->scans = NULL;
/*
 * Record the descriptor of the observation.
 */
  vp->ob = ob;
  vp->sub = NULL;
/*
 * Determine the max number of antennas per sub-array. This will be
 * used to set the max number of plots per page.
 */
  vp->maxplot = 0;
  {
    Subarray *sub;
    for(sub=ob->sub; sub < ob->sub + ob->nsub; sub++) {
      if(sub->nstat > vp->maxplot)
	vp->maxplot = sub->nstat;
    };
  };
/*
 * Allocate the subplot descriptor array.
 */
  vp->vplots = (Vissub *) malloc(vp->maxplot * sizeof(Vissub));
  if(vp->vplots == NULL) {
    lprintf(stderr,"new_Vedpar: Insufficient memory for subplot descriptors\n");
    return del_Vedpar(vp);
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
    vp->times = malloc(sizeof(TimeSample) * maxtime);
    if(!vp->times) {
      lprintf(stderr, "new_Vedpar: Insufficient memory.\n");
      return del_Vedpar(vp);
    };
  };
/*
 * Allocate a freelist for Scan elements.
 */
  vp->scan_mem = new_FreeList("new_Vedpar", sizeof(Scan), SCAN_BLK_SIZE);
  if(!vp->scan_mem)
    return del_Vedpar(vp);
/*
 * Assign defaults to the rest of the members.
 */
  vp->utref = ob->date.ut;
  vp->stref = ob->date.app_st * rtoh * 3600;
  vp->phsmin = -pi;
  vp->phsmax = pi;
  vp->ampmin = vp->ampmax = 0.0f; /* Default to auto-scaled amplitudes */
  vp->vxa = vp->vya = 0.0f;
  vp->vxb = vp->vyb = 1.0f;
  vp->modified = 0;
  vp->stat_ed = 1;                /* Default to statio based editing */
  vp->if_ed = 0;                  /* Default to global IF editing */
  vp->ch_ed = 0;                  /* Default to global channel editing */
  vp->tb = vp->ta = 0;
  vp->doamp = doamp;
  vp->dophs = dophs;
  vp->doflag = doflag;
  vp->domod = domod;
  vp->dobars = dobars;
  vp->docross = 0;
  vp->doutc = 1;
  vp->doall = 1;
  vp->showall = showall;
  vp->dodiff = 0;
  vp->nrow = 0;
  vp->nplot = 0;
  vp->nreq = nrow;
  vp->doscan = doscan;
  vp->npage = 0;
/*
 * If cursor interaction is required, check if the device has a cursor.
 */
  if(docurs) {
    slen = sizeof(awrk)-1;
    cpgqinf("CURSOR", awrk, &slen);
    docurs = strncmp(awrk,"YES",3) == 0;
  };
  vp->docurs = docurs;
  vp->cursor.key = '\0';
/*
 * Return the new descriptor.
 */
  return vp;
}

/*.......................................................................
 * Vedpar (visibility plot descriptor) destructor function.
 *
 * Input:
 *  vp     Vedpar *  Vedpar pointer returned by new_Vedpar().
 * Output:
 *  return Vedpar *  Always NULL, so that you can write vp=del_Vedpar(vp);
 */
Vedpar *del_Vedpar(Vedpar *vp)
{
  if(vp != NULL) {
    if(vp->vplots != NULL)
      free(vp->vplots);
    if(vp->times)
      free(vp->times);
    vp->scan_mem = del_FreeList("del_Vedpar", vp->scan_mem, 1);
    vp->scans = NULL;   /* NB. Already deleted by deleting scan_mem */
    free(vp);
  };
  return NULL;
}

/*.......................................................................
 * Return a list of scans in vp->scans.
 *
 * Input:
 *  vp       Vedpar *  The plot parameter container.
 *                     vp->doscan: If 0 the whole observation is
 *                                 treated as one scan.
 *                     vp->scan_mem:  The freelist of scans.
 *                     vp->scans: This is deleted first.
 *
 * Input/Output:
 *  return   Scan *  The the value of vp->scans, which will be NULL
 *                   on error.
 */
static Scan *new_Scans(Vedpar *vp)
{
  Subarray *sub;    /* The descriptor of the current sub-array */
  int ta,tb;        /* Integration indexes of start and end of a scan */
/*
 * Get the current sub-array.
 */
  sub = vp->sub;
/*
 * Start by deleting any existing scan list.
 */
  (void) del_Scans(vp);
/*
 * If not in scan mode, only one scan is needed covering the whole
 * range.
 */
  if(!vp->doscan) {
    Scan *scan = new_FreeListNode("new_Scans", vp->scan_mem);
    if(!scan)
      return NULL;
    scan->next = NULL;
    scan->vxa = vp->vxa;
    scan->vxb = vp->vxb;
    scan->stmin = vp->times[0].t;
    scan->stmax = vp->times[sub->ntime-1].t;
    scan->tmin = scan->tmax = 0.0;
    vp->scans = scan;
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
      Scan *scan = new_FreeListNode("new_Scans", vp->scan_mem);
      if(!scan)
	return del_Scans(vp);
      scan->next = NULL;
      scan->vxa = 0.0;
      scan->vxb = 0.0;
      scan->tmin = scan->tmax = 0.0;
/*
 * Get the start and end of the scan.
 */
      prev_time = scan->stmin = vp->times[ta].t;
      for(tb=ta;
	  tb < sub->ntime && vp->times[tb].t - prev_time < tsep;
	  prev_time=vp->times[tb].t, tb++)
	;
/*
 * Record the end time.
 */
      scan->stmax = prev_time;
/*
 * Append the scan to the list.
 */
      if(!vp->scans)
	vp->scans = scan;
      else
	tail->next = scan;
      tail = scan;
    };
  };
/*
 * Return the new list, to show that we succeeded.
 */
  return vp->scans;
}

/*.......................................................................
 * Delete the list of scans in vp->scans.
 *
 * Input:
 *  vp     Vedpar *   The resource object of this module.
 * Output:
 *  return Scan *   The deleted list (always NULL).
 */
static Scan *del_Scans(Vedpar *vp)
{
  Scan *next = vp->scans;  /* The next scan in the list */
/*
 * Delete the nodes of the list.
 */
  while(next) {
    Scan *scan = next;
    next = scan->next;
    scan = del_FreeListNode("del_Scans", vp->scan_mem, scan);
  };
  vp->scans = NULL;
  return NULL;
}

/*.......................................................................
 * Create a list of integrations sorted into the order of the x-axis
 * timescale.
 *
 * Input:
 *  vp       Vedpar *  The resource object of this module.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int v_get_times(Vedpar *vp)
{
  int t;   /* An index into vp->times[] */
/*
 * Get the current subarray.
 */
  Subarray *sub = vp->sub;
/*
 * Copy the pertinent times into the vp->times[] array.
 */
  for(t=0; t<sub->ntime; t++) {
    TimeSample *sample = vp->times + t;
    sample->integ = sub->integ + t;
    sample->t = v_time(vp, sample->integ->ut);
  };
/*
 * For UTC, the times are already sorted. Sidereal times however, need
 * to be sorted.
 */
  if(!vp->doutc)
    qsort(vp->times, sub->ntime, sizeof(vp->times[0]), v_cmp_time_samples);
  return 0;
}

/*.......................................................................
 * This is a qsort() comparison function for sorting arrays of TimeSample
 * elements.
 */
static int v_cmp_time_samples(const void *v1, const void *v2)
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
 * Determine the amplitude plot range for the time range and plot options
 * in a passed Vedpar descriptor. Return the new range in a passed
 * Vissub descriptor.
 *
 * Input:
 *  vp        Vedpar * This contains existing plotting attributes.
 *                     All the doxxx flags, the showall flag and uta and
 *                     utb must be initialized before calling this function.
 * Input/Output:
 *  vs        Vissub * On entry, only vs->base need be initialized.
 *                     On return vs->ampmin and vs->ampmax will contain
 *                     sugested min and max amplitudes for the plot,
 *                     including margins.
 * Output:
 *  return       int   0 - OK.
 *                     On error -1 is returned and no changes are made
 *                     to vs.
 */
static int v_arange(Vedpar *vp, Vissub *vs)
{
  int t;          /* The index of a time sample in vp->times[] */
  int flagged;    /* True if a visibility is flagged */
  int nflag=0;    /* The number of flagged data points used */
  int ngood=0;    /* The number of unflagged data points used */
  float amin;     /* Minimum amplitude in plot */
  float amax;     /* Maximum amplitude in plot */
  float adif;     /* Amplitude range before adding margins */
  float uamax;    /* Max unflagged amplitude found */
  float famax;    /* Max flagged amplitude found */
  float uamin;    /* Min unflagged amplitude found */
  float famin;    /* Min flagged amplitude found */
  float amp;      /* An amplitude */
  float amplo;    /* The bottom of the plotted amplitude point */
  float amphi;    /* The top of the plotted amplitude point */
/*
 * Check inputs.
 */
  if(vp==0) {
    lprintf(stderr, "v_arange: NULL Vedpar descriptor intercepted\n");
    return -1;
  };
/*
 * Has a fixed range been requested?
 */
  if(vp->ampmax > vp->ampmin) {
    amin = vp->ampmin;
    amax = vp->ampmax;
  } else {
/*
 * Find the max and min amplitudes of the unflagged, flagged and model data in
 * the plot range. Also count each type of datum point within the
 * given time range.
 */
    amin = 0.0f; /* Always plot zero amplitude */
    famax = uamax = 0.0f;
    famin = uamin = 0.0f;
    for(t=vp->ta; t<=vp->tb; t++) {
      Integration *integ = vp->times[t].integ;
      Visibility *vis  = integ->vis + vs->base;
/*
 * Only consider points with physical errors.
 */
      if(!(vis->bad & FLAG_DEL)) {
/*
 * Get the amplitude to be plotted.
 */
	v_data_point(vp, vis, &amp, NULL);
/*
 * Flag status of the new point?
 */
	flagged = vis->bad;
/*
 * We want the max amplitude to include the error bar when displaying
 * error bars.
 */
	if(vp->dobars) {
	  float bar = 1.0f/sqrt(fabs(vis->wt));
	  amphi = amp + bar;
	  amplo = amp - bar;
	} else {
	  amphi = amp;
	  amplo = amp;
	};
/*
 * When displaying the model, expand the range to include the model.
 */
	if(!vp->dodiff && vp->domod) {
	  float modamp = vis->modamp;
	  if(modamp > amphi)
	    amphi = modamp;
	  else if(modamp < amplo)
	    amplo = modamp;
	};
	
/*
 * Update the unflagged data range.
 */
	if(!flagged) {
	  ngood++;
	  if(amphi > uamax)
	    uamax = amphi;
	  if(amplo < uamin)
	    uamin = amplo;
/*
 * Update the flagged data range if flagged visibilities are being plotted.
 */
	} else if(vp->doflag) {
	  nflag++;
	  if(amphi > famax)
	    famax = amphi;
	  if(amplo < famin)
	    famin = amplo;
	};
      };    /* end of if(errors are physical) */
    };      /* end loop over t */
/*
 * Use the appropriate max amplitude.
 */
    if(vp->showall)                   /* Max amplitude of all data */
      amax = floatmax(famax,uamax);
    else if(ngood)                    /* Max amplitude of unflagged data */
      amax = uamax;
    else                              /* No good data so use max of flagged */
      amax = famax;
/*
 * Use the appropriate min amplitude.
 */
    if(vp->dodiff) {
      if(vp->showall)
	amin = floatmin(famin,uamin);
      else if(ngood)                    /* Min amplitude of unflagged data */
	amin = uamin;
      else                              /* No good data so use min of flagged */
	amin = famin;
    } else {
      amin = 0.0;
    };
/*
 * Add margins to the amplitude range.
 */
    adif = amax-amin;
    amin -= adif * ymarg;
    amax += adif * ymarg;
/*
 * If the derived amplitude range is zero set a default range.
 */
    if(adif==0.0f)
      amax=1.0f;
  };
/*
 * Set up return values.
 */
  vs->ampmin = amin;
  vs->ampmax = amax;
/*
 * Return the count of the total number of points in the plot range.
 */
  return 0;
}

/*.......................................................................
 * Return the X-axis plot range for the time range and plot options
 * in a passed Vedpar descriptor.
 *
 * Input/Output:
 *  vp        Vedpar * On entry this contains existing plotting attributes.
 *                     Currently only uta and utb need be initialized.
 *                     On output vp->wxa and vp->wxb will contain
 *                     the min and max times of the plot.
 * Output:
 *  return      int    0 - OK.
 *                    -1 - Error.
 */
static int v_time_range(Vedpar *vp)
{
  float xa;    /* Start time of range */
  float xb;    /* End time of range */
  Scan *scan;  /* The scan being processed */
/*
 * Check inputs.
 */
  if(vp==0) {
    lprintf(stderr, "v_time_range: NULL Vedpar descriptor intercepted\n");
    return -1;
  };
/*
 * Valid uta and utb?
 */
  if(vp->ta < 0 || vp->ta>vp->tb || vp->tb >= vp->sub->ntime) {
    lprintf(stderr, "v_time_range: uta and utb are invalid\n");
    return -1;
  };
/*
 * Determine the times corresponding to integrations uta and utb
 * with respect to the reference time vlb->ut (seconds).
 */
  vp->wxa = vp->times[vp->ta].t;
  vp->wxb = vp->times[vp->tb].t;
/*
 * Determine the displayed time ranges within the scans.
 * sc->view flags whether any of the scan is visible.
 */
  for(scan=vp->scans; scan; scan=scan->next) {
    scan->view = vp->wxb >= scan->stmin && vp->wxa <= scan->stmax;
    if(scan->view) {
      xa = (vp->wxa < scan->stmin) ? scan->stmin : vp->wxa;
      xb = (vp->wxb > scan->stmax) ? scan->stmax : vp->wxb;
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
 * Set up the viewport limits for the stack of plots leaving 4 char
 * heights on each side of plot for labelling.
 *
 * Input/Output:
 *  vp    Vedpar *  The Vis-edit parameter struct. On input vp->nplot
 *                  must be set with the number of amp+phase plots
 *                  required. vp->vplots[] must have been pre-allocated
 *                  and on output the vxa,vxb,vya,vyb fields will be
 *                  initialized. All other fields are ignored.
 *  nrow     int    The number of sub-plot slots on the display
 *  nplot    int    The actual number of sub-plots to be plotted.
 *                  nplot <= nrow.
 * Output:
 *  return   int     0 - OK.
 *                  -1 - Neither vp->doamp nor vp->dophs is true.
 */
static int v_vpwin(Vedpar *vp, int nrow, int nplot)
{
  const float phsfrc=0.3f; /* The fraction of sub-plot to assign to phase */
  Vissub *vs;     /* The sub-plot being assigned a viewport. */
  float vxa,vxb;  /* X viewport limits enclosing whole stack */
  float vya,vyb;  /* Y viewport limits enclosing whole stack */
  float tsum;     /* Sum of scan times within current time range */
  Scan *scan;     /* The scan being processed */
  int i;
/*
 * Check arguments.
 */
  if(nplot > nrow || nplot > vp->maxplot) {
    lprintf(stderr, "v_vpwin: Too many plots requested.\n");
    return -1;
  };
  if(nplot <= 0) {
    lprintf(stderr, "v_vpwin: %d plots requested\?", nplot);
    return -1;
  };
/*
 * Get the standard viewport.
 */
  cpgsch(1.0f);
  cpgvstd();
  cpgqvp(0,&vxa, &vxb, &vya, &vyb);
/*
 * Store this in the plot descriptor, for subsequent use during labelling.
 */
  vp->vxa = vxa;
  vp->vxb = vxb;
  vp->vya = vya;
  vp->vyb = vyb;
  vp->nplot = nplot;
  vp->nrow = nrow;
/*
 * Divide it into vp->nplot vertically adjacent viewports.
 */
  for(i=0; i<nplot; i++) {
    vs = &vp->vplots[i];
    vp->vxa = vxa;
    vp->vxb = vxb;
    vs->vyb = vyb - i*(vyb-vya)/nrow;
    vs->vya = vs->vyb - (vyb-vya)/nrow;
/*
 * vs->vymid specifies the bottom of the amplitude plot and the top of
 * the phase plot.
 */
    if(vp->doamp && vp->dophs)
      vs->vymid = vs->vya + phsfrc * (vs->vyb - vs->vya); /* amp+phase plot */
    else if(vp->dophs)
      vs->vymid = vs->vyb;  /* Only phase plot */
    else if(vp->doamp)
      vs->vymid = vs->vya;  /* Only amplitude plot */
    else {
      lprintf(stderr, "v_vpwin: Neither amplitude nor phase plot requested\n");
      return -1;
    };
  };
/*
 * Apportion viewports horizontally for different scans.
 * First find the sum of time ranges covered by all scans within the
 * current time range.
 */
  tsum = 0.0f;
  for(scan=vp->scans; scan; scan=scan->next)
    tsum += scan->tmax - scan->tmin;
/*
 * Use the fraction of the sum of time ranges taken up by each scan
 * to determine the fraction of the horizontal viewport range taken
 * up by that scan.
 */
  vxa = vp->vxa;
  for(scan=vp->scans; scan; scan=scan->next) {
    scan->vxa = vxa;
    if(scan->view)
      scan->vxb = vxa + (vp->vxb - vp->vxa) * (scan->tmax - scan->tmin) / tsum;
    else
      scan->vxb = scan->vxa;    /* Scan not visible */
    vxa = scan->vxb;
  };
/*
 * Scale the character height with the number of plots.
 */
  cpgsch(3.0f/vp->nplot);
  return 0;
}

/*.......................................................................
 * Draw axes for a given sub-plot.
 *
 * Input:
 *  vp      Vedpar *  The plot descriptor.
 *  vs      Vissub *  The sub-plot descriptor.
 *  dotop      int    If true draw ticked axis along top axis of viewport.
 *  dobot      int    If true draw ticked and labelled axis along bottom
 *                    of the sub-plot viewport.
 *  erase      int    If true erase current axes instead of plotting.
 * Output:
 *  return     int    0 - OK.
 *                    Anything else if an error occured.
 */
static int v_plaxes(Vedpar *vp, Vissub *vs, int dotop, int dobot, int erase)
{
  Subarray *sub; /* Local pointer to sub-array descriptor */
  Scan *scan;    /* The scan being labelled */
  float tmin;    /* Start time */
  float tmax;    /* End time */
  float ch;      /* Character height to use */
  int oldcol;    /* Color index on entry to function */
  char label[25];/* Temporary string to compose baseline label in */
/*
 * Check arguments.
 */
  if(vp==NULL || vs==NULL) {
    lprintf(stderr, "v_plaxes: NULL %s descriptor intercepted\n",
	    (vp==NULL)?"plot":"sub-plot");
    return -1;
  };
/*
 * Get a local pointer to the current sub-array descriptor.
 */
  sub = vp->sub;
/*
 * Store current color index, to be restored on return.
 */
  cpgqci(&oldcol);
/*
 * Set color for drawing or erasing.
 */
  cpgsci(erase ? 0:1);
/*
 * Work out an appropriate character height, given the plot density.
 */
  ch = 1.0f/sqrt((double) vp->nplot);
/*
 * Plot the two Y-axes at each end of the frame enclosing the scans.
 */
  cpgsch(ch);
  if(vp->doamp) {
    cpgsvp(vp->vxa, vp->vxb, vs->vymid, vs->vyb);
    cpgswin(0.0f, 1.0f, vs->ampmin, vs->ampmax);
    cpgbox(" ", 0.0f, 0, "BCVNST", 0.0f, 0);
  };
  if(vp->dophs) {
    cpgsvp(vp->vxa, vp->vxb, vs->vya, vs->vymid);
    cpgswin(0.0f, 1.0f, vp->phsmin * rtod, vp->phsmax * rtod);
    cpgbox(" ", 0.0f, 0, "BCVNST", 0.0f, 0);
  };
/*
 * Do internal and X-axes for each visible scan.
 */
  for(scan=vp->scans; scan; scan=scan->next) {
/*
 * If the scan isn't visible, skip to the next one.
 */
    if(!scan->view)
      continue;
/*
 * Calculate the start and end time in seconds. For UTC,
 * add one day such that days in the year start from 1 rather than 0.
 */
    if(vp->doutc) {
      tmin = vp->utref + scan->tmin + daysec;
      tmax = vp->utref + scan->tmax + daysec;
    } else {
      tmin = scan->tmin;
      tmax = scan->tmax;
    };
/*
 * Draw internal Y-axes as unadorned vertical lines.
 */
    cpgsvp(vp->vxa, vp->vxb, vp->vya, vp->vyb);
    cpgswin(vp->vxa, vp->vxb, vp->vya, vp->vyb);    
    if(scan->next && scan->next->view) {
      cpgmove(scan->vxb, vs->vya);
      cpgdraw(scan->vxb, vs->vyb);
    };
/*
 * Draw X-axes.
 */
    if(vp->doamp && vp->dophs) {
      cpgmove(scan->vxa, vs->vymid);
      cpgdraw(scan->vxb, vs->vymid);
    };
/*
 * Write numeric labels under the last plot.
 */
    cpgsvp(scan->vxa, scan->vxb, vs->vya, vs->vyb);
    cpgswin(tmin, tmax, 0.0f, 1.0f);
    cpgsch(dotop?0.7f:ch);
    cpgtbox("ZHCST", 0.0f, 0, " ", 0.0f, 0);
    cpgsch(dobot ? 0.7f : ch);
    cpgtbox(dobot ? "ZHBNST" : "ZHBST", 0.0f, 0, " ", 0.0f, 0);
  };
/*
 * Set viewport around whole sub-plot and write a baseline label
 * inside the top right hand corner.
 */
  cpgsvp(vp->vxa, vp->vxb, vs->vya, vs->vyb);
  sprintf(label, "%.10s-%.10s", sub->tel[sub->base[vs->base].tel_a].name,
	  sub->tel[sub->base[vs->base].tel_b].name);
  cpgsch(0.5f);
  cpgmtxt("T", -1.5f, 0.99f, 1.0f, label);
/*
 * Restore entry color.
 */
  cpgsci(oldcol);
  return 0;
}

/*.......................................................................
 * Plot or erase amplitude and phase points.
 *
 * Input:
 *  vp       Vedpar * The plot descriptor.
 *  vs       Vissub * The subplot descriptor.
 *  ta          int   Index of first sample from vp->times[] to be plotted.
 *  tb          int   Index of second sample from vp->times[] to be plotted.
 *  erase       int   If true erase points instead of plotting them.
 * Output:
 *  return      int   0 - OK.
 */
int v_pldata(Vedpar *vp, Vissub *vs, int ta, int tb, int erase)
{
  Subarray *sub;      /* Local pointer to the sub-array being displayed */
  int oldcol;         /* Color index on entry to function */
  float amp,phs;      /* Amp and phase in normalised viewport coordinates */
  float phserr;       /* The phase error in normalised viewport coords */
  float amperr;       /* The amplitude error in normalised viewport coords */
  int isflagd;        /* If 1 then the visibility is flagged */
  int isym;           /* PGPLOT marker symbol to be plotted */
  int icol;           /* PGPLOT color index of plotted points */
  int base;           /* The index of the baseline being plotted */
  int t;              /* Index of integration being plotted */
/*
 * Start pgplot buffering.
 */
  cpgbbuf();
/*
 * Store current color index, to be restored on return.
 */
  cpgqci(&oldcol);
/*
 * Set color for drawing or erasing.
 */
  if(erase)
    cpgsci(0);
  cpgsch(1.0f);
/*
 * Get a local copies of the baseline index and the pointer to the sub-array
 * descriptor.
 */
  base = vs->base;
  sub = vp->sub;
/*
 * Draw each point with the appropriate symbol and color for its
 * flag status.
 *
 * Start with the amplitude points.
 */
  if(vp->doamp) {
    Scan *scan = vp->scans;
    int first=1;
    for(t=ta; t<= tb; t++) {
      TimeSample *sample = vp->times + t;
      Visibility *vis = sample->integ->vis + base;
      float tval = sample->t;
/*
 * Skip to the right scan for this point.
 */
      if(first || tval > scan->stmax) {
	first = 0;
	while(tval > scan->stmax)
	  scan = scan->next;
	cpgsvp(scan->vxa, scan->vxb, vs->vymid, vs->vyb);
	cpgswin(scan->tmin, scan->tmax, vs->ampmin, vs->ampmax);
      };
/*
 * Ignore data with unphysical zero errors.
 */
      if(vis->bad & FLAG_DEL)
	continue;
/*
 * Is the input visibility flagged?
 */
      isflagd = vis->bad;
/*
 * Should we plot this point and if so, which color and symbol should
 * be used?
 */
      if(!isflagd || vp->doflag) {
	if(isflagd) {               /* Flagged correction */
	  if(vis->bad & FLAG_BAD) { /* Visibility flag */
	    isym = badsym;
	    icol = badcol;
	  } else {                  /* Selfcal correction flag */
	    isym = badcsym;
	    icol = badccol;
	  };
	} else {                    /* Not flagged */
	  isym = datsym;
	  icol = datcol;
	};
/*
 * Install the new color.
 */
	cpgsci(erase ? 0 : icol);
/*
 * Get the world coordinates to be plotted.
 */
	v_data_point(vp, vis, &amp, NULL);
	amperr = 1.0f/sqrt(fabs(vis->wt));
/*
 * Plot the point.
 */
	cpgpt(1, &tval, &amp, isym);
/*
 * Plot error bars if requested.
 */
	if(vp->dobars) {
	  cpgmove(tval, amp - amperr);
	  cpgdraw(tval, amp + amperr);
	};
      };  /* End of  if not flagged or allowed to plot flagged data */
    };    /* End of  loop over integrations */
  };      /* End of  if plot amplitudes */
/*
 * Now plot the phases.
 */
  if(vp->dophs) {
    Scan *scan = vp->scans;
    int first = 1;
    for(t=ta; t<= tb; t++) {
      TimeSample *sample = vp->times + t;
      Visibility *vis = sample->integ->vis + base;
      float tval = sample->t;
/*
 * Skip to the right scan for this point.
 */
      if(first || tval > scan->stmax) {
	first = 0;
	while(tval > scan->stmax)
	  scan = scan->next;
	cpgsvp(scan->vxa, scan->vxb, vs->vya, vs->vymid);
	cpgswin(scan->tmin, scan->tmax, vp->phsmin, vp->phsmax);
      };
/*
 * Ignore data with unphysical zero errors.
 */
      if(vis->bad & FLAG_DEL)
	continue;
/*
 * Is the input visibility flagged?
 */
      isflagd = vis->bad;
/*
 * Should we plot this point and if so, which color and symbol should
 * be used?
 */
      if(!isflagd || vp->doflag) {
	if(isflagd) {               /* Flagged correction */
	  if(vis->bad & FLAG_BAD) { /* Visibility flag */
	    isym = badsym;
	    icol = badcol;
	  } else {                  /* Selfcal correction flag */
	    isym = badcsym;
	    icol = badccol;
	  };
	} else {                    /* Not flagged */
	  isym = datsym;
	  icol = datcol;
	};
/*
 * Install the new color.
 */
	cpgsci(erase ? 0 : icol);
/*
 * Get the world coordinates to be plotted.
 */
	v_data_point(vp, vis, &amp, &phs);
	phs -= twopi*floor(phs/twopi+0.5);  /* Wrap into -pi to pi range */
	if(amp > 1.0e-20f)
	  phserr = 1.0f/sqrt(fabs(vis->wt)) / amp;
	else
	  phserr = 0.0f;
/*
 * Plot the point.
 */
	cpgpt(1, &tval, &phs, isym);
/*
 * Plot error bars if requested.
 */
	if(vp->dobars) {
	  cpgmove(tval, phs - phserr);
	  cpgdraw(tval, phs + phserr);
	};
      };  /* End of  if not flagged or allowed to plot flagged data */
    };    /* End of  loop over integrations */
  };      /* End of  if plot phases */
/*
 * Restore entry color and terminate pgplot buffering.
 */
  cpgsci(oldcol);
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Plot or erase amplitude and phase model lines.
 *
 * Input:
 *  vp       Vedpar * The plot descriptor.
 *  vs       Vissub * The sub-plot descriptor.
 *  erase       int   If true erase points instead of plotting them.
 * Output:
 *  return      int   0 - OK.
 */
int v_plmodel(Vedpar *vp, Vissub *vs, int erase)
{
  Subarray *sub;  /* The descriptor of the current sub-array */
  int oldcol;     /* Color index on entry to function */
  int t;          /* The vp->times[] index of the integration being plotted */
  int base;       /* The baseline being plotted */
  float amp,phs;  /* The amplitude and phase to be plotted */
/*
 * Do nothing if no model exists or vp->domod==0, or vp->dodiff!=0.
 */
  if(!vp->ob->hasmod || !vp->domod || vp->dodiff)
    return 0;
/*
 * Store current color index, to be restored on return.
 */
  cpgqci(&oldcol);
/*
 * Set color for drawing or erasing.
 */
  if(erase)
    cpgsci(0);
  else
    cpgsci(modcol);
/*
 * Get the sub-array descriptor.
 */
  sub = vp->sub;
/*
 * Get a local copy of the baseline index to aid optimization.
 */
  base = vs->base;
/*
 * Turn on pgplot buffering.
 */
  cpgbbuf();
/*
 * Start with the amplitude model.
 */
  if(vp->doamp) {
    Scan *scan = vp->scans;
    int first=1;
    float prevut =0.0f;   /* Time of previous plotted amplitude */
    for(t=vp->ta; t<=vp->tb; t++) {
      TimeSample *sample = vp->times + t;
      Visibility *vis = sample->integ->vis + base;
      float tval = sample->t;
/*
 * Ignore deleted data. Also ignore the model of flagged data when
 * flagged points aren't being displayed.
 */
      if(!(vis->bad & FLAG_DEL) && (vp->doflag || !vis->bad)) {
	amp = vis->modamp;
/*
 * If not in the current scan - skip to the right scan and position
 * for the start of the new model line.
 */
	if(first || tval > scan->stmax || tval-prevut > sub->scangap) {
	  while(tval > scan->stmax)
	    scan = scan->next;
	  cpgsvp(scan->vxa, scan->vxb, vs->vymid, vs->vyb);
	  cpgswin(scan->tmin, scan->tmax, vs->ampmin, vs->ampmax);
	  first = 0;
	  cpgmove(tval, amp);
	} else {
	  cpgdraw(tval, amp); /* Continue the model line */
	};
	prevut = tval;
      };
    };
  };
/*
 * Now do the phase model.
 */
  if(vp->dophs) {
    Scan *scan = vp->scans;
    int first=1;
    float prevphs=0.0f;   /* Previous plotted phase */
    float prevut =0.0f;   /* Time of previous plotted phase */
    for(t=vp->ta; t<=vp->tb; t++) {
      TimeSample *sample = vp->times + t;
      Visibility *vis = sample->integ->vis + base;
      float tval = sample->t;
/*
 * Ignore deleted data. Also ignore the model of flagged data when
 * flagged points aren't being displayed.
 */
      if(!(vis->bad & FLAG_DEL) && (vp->doflag || !vis->bad)) {
	phs = vis->modphs;
/*
 * Wrap the phase into the range -pi to pi.
 */
	phs -= twopi * floor(phs/twopi + 0.5);
/*
 * If not in the current scan - skip to the right scan and position
 * for the start of the new model line.
 */
	if(first || tval > scan->stmax || tval-prevut > sub->scangap) {
	  while(tval > scan->stmax)
	    scan = scan->next;
	  cpgsvp(scan->vxa, scan->vxb, vs->vya, vs->vymid);
	  cpgswin(scan->tmin, scan->tmax, vp->phsmin, vp->phsmax);
	  first = 0;
	  cpgmove(tval, phs);
	} else {
	  float phsdif = phs - prevphs; /* Phase excursion wrt previous phase */
/*
 * Because the phase is only known modulo 2.pi radians, there are three
 * possible paths between the previous and current phase point within the
 * chosen -pi to pi degree plot range. Take the shortest of these.
 */
	  if(phsdif > pi) {
	    cpgdraw(tval, phs-twopi);
	    cpgmove(prevut, prevphs+twopi);
	    cpgdraw(tval, phs);
	  } else if(phsdif < -pi) {
	    cpgdraw(tval, phs+twopi);
	    cpgmove(prevut, prevphs-twopi);
	    cpgdraw(tval, phs);
	  } else {
	    cpgdraw(tval, phs); /* Continue the model line */
	  };
	};
	prevut = tval;
	prevphs = phs;
      };
    };
  };
/*
 * Restore entry color and terminate pgplot buffering.
 */
  cpgsci(oldcol);
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Read the cursor position and return the cursor selection details in
 * vp->cursor.
 *
 * Input:
 *  vp     Vedpar *  The plot descriptor.
 *  noout     int    If true then don't return until the cursor is
 *                   pressed inside a sub-plot.
 *  mode Bandmode    The desired type of cursor, from:
 *                    B_NORM - A single point is required - no banding.
 *                    B_LINE - Line band between vcref and the cursor.
 *                    B_RECT - Rectangular band between xref,yref and the
 *                             cursor.
 *                    B_YRNG - Two horizontal lines bracketing a Y-axis range.
 *                    B_XRNG - Two vertical lines bracketing an X-axis range.
 *                    B_YVAL - Vertical line through the cursor.
 *                    B_XVAL - Horizontal line through the cursor.
 *                    B_CROSS- Cross hair centered on cursor.
 *  isamp     int    0 - yref denotes a phase.
 *                   1 - yref denotes an amplitude.
 *  vsref  Vissub *  The descriptor of the sub-plot to which xref,yref refer.
 *                   This can be NULL if yref is not relevant.
 *  xref    float    xref and yref denote the reference position for
 *  yref    float    band-type cursors. They are ignored if mode==B_NORM.
 *  ci        int    The color index to draw band cursors with.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
int v_cursor(Vedpar *vp, int noout, Bandmode mode, int isamp, Vissub *vsref,
		    float xref, float yref, int ci)
{
  static float xpos=0.5f; /* The X NDC position of the cursor */
  static float ypos=0.5f; /* The Y NDC position of the cursor */
  char key;               /* The cursor selection key */
  Vissub *vs=NULL;        /* The subplot being checked */
  Vcurs *vc;              /* Pointer to vp->cursor */
  Scan *scan;             /* The time scan being checked */
  int i;
/*
 * Get the cursor descriptor.
 */
  vc = &vp->cursor;  
/*
 * Set the viewport around the whole viewsurface and make the world
 * coords the same as NDC so that the returned cursor position
 * is measured in NDC.
 */
  cpgsvp(0.0f,1.0f,0.0f,1.0f);
  cpgswin(0.0f,1.0f,0.0f,1.0f);
/*
 * If this is the first call of the plot session initialize the position
 * at which to bring up the cursor. Otherwise use the values retained from
 * the previous call in xpos and ypos.
 */
  if(vc->key == '\0') {
    xpos = 0.5f;
    ypos = 0.5f;
  };
/*
 * Initialize the return value.
 */
  vc->key = '\0';
  vc->waslow = 0;
  vc->wasamp = 0;
  vc->vs = NULL;
  vc->iplot = 0;
  vc->scan = NULL;
  vc->tval = 0.0f;
  vc->value = 0.0f;
/*
 * Substitute the cross-hair cursor for the normal cursor is requested.
 */
  if(mode == B_NORM && vp->docross)
    mode = B_CROSS;
/*
 * Convert the cursor reference positions into NDC.
 */
  switch(mode) {
  case B_RECT: case B_XRNG: case B_YRNG:
    {
      Scan *tail = vp->scans;  /* The tail of the list of scans */
/*
 * Locate the scan that contains the reference time.
 */
      for(scan=vp->scans; scan; scan=scan->next) {
	if(xref >= scan->tmin && xref <= scan->tmax)
	  break;
	tail = scan;
      };
      if(!scan)
	scan = xref < vp->scans->tmin ? vp->scans : tail;
/*
 * Convert the reference time and amp/phase to the equivalent NDC position.
 */
      xref = scan->vxa + (xref - scan->tmin) * (scan->vxb - scan->vxa) /
	(scan->tmax - scan->tmin);
/*
 * Get the Y-axis reference value.
 */
      if(vsref==NULL) {
	yref = 0.0;
      } else {
	yref = isamp ?
	  vsref->vymid + (yref - vsref->ampmin) * (vsref->vyb - vsref->vymid) /
	    (vsref->ampmax - vsref->ampmin):
          vsref->vya   + (yref - vp->phsmin) * (vsref->vymid - vsref->vya) /
	    (vp->phsmax - vp->phsmin);
      };
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
    cpgsci(ci);
    if(!cpgband((int) mode, 0, xref, yref, &xpos, &ypos, &key))
      return 1;
/*
 * Convert key to upper case.
 */
    vc->waslow = islower((int)key);
    vc->key = (vc->waslow) ? toupper((int)key) : key;
/*
 * See if the point is in any sub-plot.
 */
    for(i=0; i<vp->nplot; i++) {
      vs = &vp->vplots[i];
      if((xpos >= vp->vxa && xpos <= vp->vxb) && 
	 (ypos >= vs->vya && ypos <= vs->vyb) )
	break;
    };
/*
 * Was the cursor in a subplot?
 */
    if(i < vp->nplot) {
      vc->vs = vs;
      vc->wasamp = vp->doamp && ypos > vs->vymid; /* In amp part of plot? */
/*
 * Convert from NDC to the respective amplitude or phase selected.
 */
      if(vc->wasamp) {
	vc->value = vs->ampmin + 
	  (ypos - vs->vymid)/(vs->vyb - vs->vymid) * (vs->ampmax - vs->ampmin);
      } else {
	vc->value = vp->phsmin +
	  (ypos - vs->vya)/(vs->vymid - vs->vya) * (vp->phsmax - vp->phsmin);
      };
/*
 * Identify the scan that the cursor was in and use this to
 * determine the selected time value.
 */
      for(scan=vp->scans; vc->scan==NULL && scan; scan=scan->next) {
	if(xpos >= scan->vxa && xpos <= scan->vxb) {
	  vc->tval = scan->tmin + (xpos - scan->vxa) /
	    (scan->vxb - scan->vxa) * (scan->tmax - scan->tmin);
	  vc->scan = scan;
	};
      };
/*
 * Cursor not in any sub-plot.
 */
    } else {
      vc->vs = vs = NULL;
      vc->tval = vc->value = 0.0f;
      vc->wasamp = 0;
      vc->scan = NULL;
    };
    if(!vs && noout)
      printf("The cursor must be in one of the plots.\n");
  } while(vs==NULL && noout); /* Repeat if outside plots and noout is true */
  return 0;
}

/*.......................................................................
 * Write labels around the frame enclosing all sub-plots.
 *
 * Input:
 *  vp       Vedpar *  The plot descriptor.
 * Output:
 *  return      int    0 - OK.
 */
int v_label(Vedpar *vp)
{
  Observation *ob;      /* The descriptor of the observation being plotted */
  Subarray *sub;        /* The descriptor of the sub-array being plotted */
  char awrk[81];        /* Work string for labelling */
  char bwrk[81];        /* Work string for labelling */
/*
 * Check arguments.
 */
  if(vp==NULL) {
    lprintf(stderr, "v_label: NULL plot descriptor intercepted\n");
    return -1;
  };
/*
 * Get the descriptors of the observation and sub-array being plotted.
 */
  ob = vp->ob;
  sub = vp->sub;
/*
 * Set the viewport around the plot grid.
 */
  cpgsvp(vp->vxa, vp->vxb, vp->vplots[vp->nplot-1].vya, vp->vyb);
/*
 * Compose and write main title.
 */
  cpgsci(1);
  cpgsch(1.0f);
  sprintf(awrk, "%s  %s", ob->source.name,
	  sutdate(ob->date.year, ob->date.ut, bwrk));
  cpgmtxt("T", 1.7f, 0.0f, 0.0f, awrk);
  sprintf(awrk,"%s of %d:%s in IF %d, Pol %s",
	  vp->doall ? "Baselines":"Upper baselines",
          vp->bs_beg.isub+1, sub->tel[vp->bs_beg.ta].name, ob->stream.cif+1,
	  Stokes_name(ob->stream.pol.type));
  cpgmtxt("T",0.5f,0.0f,0.0f,awrk);
/*
 * In non-interactive mode tell the user what is being plotted.
 */
  if(!vp->docurs) {
    lprintf(stdout, "Page %02.2d: %s of %d:%s\n", vp->npage,
	    vp->doall ? "Baselines":"Upper baselines",
	    vp->bs_beg.isub+1, sub->tel[vp->bs_beg.ta].name);
  };
/*
 * Write Y labels.
 */
  sprintf(awrk, "%s%s%s%s", vp->dophs?"Phase":"",
	  vp->dophs&&vp->doamp?" and ":"",
	  vp->doamp?"Amplitude":"",
	  vp->dodiff ? " residuals":"");
  cpgmtxt("L",3.0f,0.5f,0.5f,awrk);
/*
 * Write the X-axis label.
 */
  if(vp->doutc)
    strcpy(awrk, "Universal Time");
  else
    strcpy(awrk, "Greenwhich Mean Sidereal Time");
  cpgmtxt("B", 2.5f, 0.5f, 0.5f, awrk);
  return 0;
}

/*.......................................................................
 * Plot an extra mode label for editting sessions.
 *
 * Input:
 *  vp       Vedpar *  The plot descriptor.
 *  erase       int    If true, erase existing mode label.
 * Output:
 *  return      int    0 - OK.
 */
static int v_mlab(Vedpar *vp, int erase)
{
  Observation *ob;  /* The descriptor of the observation being plotted */
  int oldcol;       /* Temporary storage for entry color index */
  char label[81];   /* Temporary work string to compose mode label in */
/*
 * Get the descriptor of the observation.
 */
  ob = vp->ob;
/*
 * Store the existing plot color.
 */
  cpgqci(&oldcol);
/*
 * Set color to erase or draw.
 */
  if(erase)
    cpgsci(0);
  else
    cpgsci(1);
/*
 * Set viewport around sub-plots.
 */
  cpgsvp(vp->vxa, vp->vxb, vp->vya, vp->vyb);
/*
 * Compose the mode label.
 */
  sprintf(label, "%s editing of %s channels of %s.",
	  vp->stat_ed ? "Station" : "Baseline",
	  vp->ch_ed ? "selected":"all",
	  vp->if_ed ? "the displayed IF" : "all IFs");
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
 *  vp       Vedpar *  The plot descriptor.
 *  stat_ed     int    Select station based editing if true.
 *  if_ed       int    Select IF based editing if true.
 *  ch_ed       int    Select channel based editing if true.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int v_newmode(Vedpar *vp, int stat_ed, int if_ed, int ch_ed)
{
/*
 * Buffer until the new text has been plotted.
 */
  cpgbbuf();
/*
 * Erase the existing mode line.
 */
  v_mlab(vp, 1);
/*
 * Install the new editing modes.
 */
  vp->stat_ed = stat_ed;
  vp->if_ed = if_ed;
  vp->ch_ed = ch_ed;
/*
 * Draw the new mode line.
 */
  v_mlab(vp, 0); /* Plot new mode line */
/*
 * reveal the changes.
 */
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Replot the current plots to reflect new attribute selections such as
 * a new time range. This function should not be called until the first
 * succesful call to v_plot() with oper=V_ALLNEW has been made.
 *
 * Input:
 *  vp        Vedpar *  The plot descriptor.
 * Output:
 *  return       int    0 - OK.
 */
int v_redisp(Vedpar *vp)
{
  Vissub *vs; /* Pointer to current sub-plot descriptor */
  int iplot;  /* Number of sub-plot being drawn */
  int ierr=0; /* True if an error occurs */
/*
 * Cursory check of arguments.
 */
  if(vp==NULL) {
    lprintf(stderr, "v_redisp: NULL plot descriptor intercepted\n");
    return -1;
  };
/*
 * Nothing to plot?
 */
  if(vp->nplot<=0) {
    lprintf(stderr, "v_redisp: No plot rows have been initialized.\n");
    return -1;
  };
/*
 * Clear page.
 */
  cpgpage();
/*
 * Count pages.
 */
  vp->npage++;
/*
 * Determine the time plot range for all plots.
 */
  ierr = ierr || v_time_range(vp);
/*
 * Set up viewport slots for each sub-plot.
 */
  ierr = ierr || v_vpwin(vp, vp->nrow, vp->nplot);
/*
 * Plot each sub-plot.
 */
  for(iplot=0; iplot<vp->nplot && !ierr; iplot++) {
    vs = &vp->vplots[iplot];
    cpgbbuf();
    ierr = ierr || v_arange(vp, vs);
    ierr = ierr || v_plaxes(vp, vs, iplot==0, iplot==vp->nplot-1, 0);
    ierr = ierr || v_pldata(vp, vs, vp->ta, vp->tb, 0);
    ierr = ierr || v_plmodel(vp, vs, 0);
    ierr = ierr || (iplot==0 && v_label(vp));
    if(vp->docurs)
      ierr = ierr || v_mlab(vp, 0);
    cpgebuf();
  };
  return ierr;
}

/*.......................................................................
 * Display baselines of a given reference telescope in a given sub-array
 * and in a specified order. Optionally move to the next reference
 * telescope and optionally the next sub-array in the specified direction
 * when no baselines remain to be plotted to the given specification.
 *
 * Input:
 *  vp       Vedpar *  The plot descriptor.
 *  oper      Vedop    The action to take to display the next page.
 *                      V_ALLNEW - Start new plot using the given *init
 *                                 baseline specification descriptor.
 *                      V_REPLOT - Re-initialize the current plot to account
 *                                 for attribute changes.
 *                      V_RESET  - Replot from the start/end of the current
 *                                 telescope.
 *                      V_NXTSUB - Start plotting from the next available
 *                                 sub-array in the given direction.
 *                      V_NXT_TA - Start plotting for the next available
 *                                 reference telescope in the given direction.
 *                      V_NXT_TB - Start plotting the next page of
 *                                 baselines in direction dir. 
 *                      V_NEXT   - Start plotting the next page of
 *                                 baselines in direction dir if there are
 *                                 any remaining baselines matching the last
 *                                 baseline specification.
 *  forward     int    The order in which to search for plottable baselines.
 *                       0 - Search in the direction of decreasing
 *                           baseline order.
 *                       1 - Search in the direction of increasing
 *                           baseline order.
 *  init   Basespec *  If oper==V_ALLNEW, then *init will be used to
 *                     replace the current start baseline spec. Otherwise
 *                     'init' will be ignored, and should be sent as NULL.
 * Output:
 *  return      int    The number of sub-plots plotted. Or -1 on
 *                     error.
 */
int v_plot(Vedpar *vp, Vedop oper, int forward, Basespec *init)
{
  Observation *ob; /* The descriptor of the parent observation */
  Basespec bs;     /* The new trial baseline specification */
/*
 * Check arguments.
 */
  if(vp==NULL) {
    lprintf(stderr, "v_plot: NULL plot descriptor intercepted\n");
    return -1;
  };
/*
 * Get the descriptor of the parent observation.
 */
  ob = vp->ob;
/*
 * First call must use V_ALLNEW.
 */
  if(vp->nplot<1 && oper!=V_ALLNEW) {
    lprintf(stderr, "v_plot: First call must use V_ALLNEW.\n");
    return -1;
  };
/*
 * Perform any positioning operations required for the particular
 * opcode.
 */
  switch(oper) {
  case V_ALLNEW:
    if(!init) {
      lprintf(stderr, "v_plot: Invalid baseline specification.\n");
      return -1;
    };
    bs = *init;
    if(next_base(ob, FIND_FIRST, forward, 2, vp->doall, 0, 1, &bs))
      return 0;
    break;
  case V_REPLOT:
    bs = vp->bs_beg;
    break;
  case V_RESET:
    bs = vp->bs_beg;
    if(next_base(ob, FIND_FIRST, forward, 2, vp->doall, 1, 1, &bs))
      return 0;
    break;
  case V_NXTSUB:
    bs = vp->bs_beg;
    if(next_base(ob, SKIP_SUB, forward, 2, vp->doall, 0, 1, &bs))
      return 0;
    break;
  case V_NXT_TA:
    bs = vp->bs_beg;
    if(next_base(ob, SKIP_TA, forward, 2, vp->doall, 0, 0, &bs) &&
       next_base(ob, SKIP_SUB,forward, 2, vp->doall, 0, 1, &bs))
      return 0;
    break;
  case V_NXT_TB:
    bs = forward ? vp->bs_end : vp->bs_beg;
    if(next_base(ob, SKIP_TB, forward, 2, vp->doall, 0, 0, &bs) &&
       next_base(ob, SKIP_TA, forward, 2, vp->doall, 0, 0, &bs) &&
       next_base(ob, SKIP_SUB,forward, 2, vp->doall, 0, 1, &bs)) 
      return 0;
    break;
  case V_NEXT:
    bs = forward ? vp->bs_end : vp->bs_beg;
    if(next_base(ob, FIND_NEXT, forward, 2, vp->doall, 0, 0, &bs))
      return 0;
    break;
  default:
    lprintf(stderr, "v_plot: Unknown opcode intercepted.\n");
    return -1;
    break;
  };
/*
 * Record the start baseline spec.
 */
  vp->bs_beg = bs;
  vp->bs_end = bs;
/*
 * If the sub-array has changed, prepare for the new sub-array.
 */
  if(vp->nplot<1 || vp->sub != ob->sub + bs.isub) {
/*
 * Record the index and descriptor of the new sub-array.
 */
    vp->sub = ob->sub + bs.isub;
/*
 * Set up the number of plots per page.
 */
    v_setnrow(vp, vp->nreq);
/*
 * Set up for the full time range of the new sub-array.
 */
    vp->ta = 0;
    vp->tb = vp->sub->ntime-1;
/*
 * Update the array of times and the list of scans.
 */
    if(v_update_times(vp))
      return -1;
  };
/*
 * Locate the rest of the baselines of the current sub-array and reference
 * telescope that can be plotted on the new page.
 */
  vp->nplot = 0;
  do {
    vp->vplots[vp->nplot].base = bs.base; /* Record the sub-plot baseline */
    vp->nplot++;                          /* Record addition to sub-plot list */
/*
 * Record the spec of the latest baseline.
 */
    if(forward)
      vp->bs_end = bs;
    else
      vp->bs_beg = bs;
  } while(vp->nplot < vp->nrow && next_base(ob, FIND_NEXT, forward, 2,
					    vp->doall, 1, 0, &bs)==0);
/*
 * If we were searching in reverse, the baselines in vp->vplots[] will now
 * be in reverse plotting order. Rearrange them back into the forward order.
 */
  if(!forward) {
    Vissub *vsa = &vp->vplots[0];
    Vissub *vsb = &vp->vplots[vp->nplot-1];
    for( ; vsa<vsb; vsa++,vsb--) {
      int base = vsa->base;
      vsa->base = vsb->base;
      vsb->base = base;
    };
  };
/*
 * Display the new baselines.
 */
  if(v_redisp(vp))
    return -1;
/*
 * Return the number of baselines plotted.
 */
  return vp->nplot;
}

/*.......................................................................
 * Determine scaling factors required to convert from world coordinates
 * to mm in the gievn partition of a given sub-plot.
 *
 * Input:
 *  vp   Vedpar *  The plot descriptor.
 *  vs   Vissub *  The sub-plot descriptor.
 *  doamp   int    If true make ytomm the conversion from amplitude
 *                 to mm. Otherwise yield the phase conversion.
 * Output:
 *  xtomm float *  t_shift * *xtomm yields the physical size
 *                 of the t shift in mm.
 *  ytomm float *  Amplitude_or_phase_shift * *xtomm yields the physical
 *                 size of the amplitude or phase shift in mm.
 *  return  int    0 - OK.
 */
int v_scale(Vedpar *vp, Vissub *vs, int doamp, float *xtomm, float *ytomm)
{
  float xa,xb,ya,yb; /* Physical coordinates of viewport */
  Scan *scan;        /* The first visible scan */
/*
 * Find the first displayed scan.
 */
  for(scan=vp->scans; scan && !scan->view; scan=scan->next)
    ;
  if(!scan) {
    lprintf(stderr, "v_scale: No scans visible\n");
    return -1;
  };
/*
 * Determine the size of the viewport in physical device coordinates
 * (millimeters).
 */
  if(doamp)
    cpgsvp(scan->vxa, scan->vxb, vs->vymid, vs->vyb);
  else
    cpgsvp(scan->vxa, scan->vxb, vs->vya, vs->vymid);
  cpgqvp(2,&xa,&xb,&ya,&yb);
/*
 * Calculate factors to convert world coords into mm.
 */
  *xtomm = fabs((xb - xa)/(scan->tmax - scan->tmin));
  if(doamp)
    *ytomm = fabs((yb - ya)/(vs->ampmax - vs->ampmin));
  else
    *ytomm = fabs((yb - ya)/(vp->phsmax - vp->phsmin));
  return 0;
}

/*.......................................................................
 * Handle a user request for a new number of plot slots.
 *
 * Input:
 *  vp      Vedpar *  The plot descriptor.
 *  nreq       int    The requested number of plots per page.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
int v_setnrow(Vedpar *vp, int nreq)
{
  int nmax;  /* The max number of plots per page in the current sub-array */
/*
 * Check arguments.
 */
  if(vp==NULL) {
    lprintf(stderr, "v_setnrow: NULL plot descriptor intercepted.\n");
    return 1;
  };
/*
 * Record the requested number of plots for future sub-array selections.
 */
  vp->nreq = nreq;
/*
 * How many plots per page are possible in the current sub-array?
 */
  nmax = vp->sub->nstat - 1;
/*
 * Now apply bounds for the current sub-array, recording the actual number
 * of slots to be used, in vp->nrow.
 */
  vp->nrow = (nreq>0 && nreq<=nmax) ? nreq : nmax;
  return 0;
}

/*.......................................................................
 * Return the appropriate amplitude and/or phase of a visibility for
 * plotting.
 *
 * Input:
 *  vp        Vedpar *   The plot descriptor.
 *  vis   Visibility *   The visibility who's amplitude is required.
 *  amp        float *   The value to be plotted in the amplitude plot.
 *                       Pass amp=NULL if not required.
 *  phs        float *   The value to be plotted in the amplitude plot.
 *                       Pass phs=NULL if not required.
 */
void v_data_point(Vedpar *vp, Visibility *vis, float *amp, float *phs)
{
/*
 * Are we plotting difference between the model and the data?
 */
  if(vp->dodiff) {
    float re = vis->amp * cos(vis->phs) - vis->modamp * cos(vis->modphs);
    float im = vis->amp * sin(vis->phs) - vis->modamp * sin(vis->modphs);
    if(amp)
      *amp = sqrt(re * re + im * im);
    if(phs)
      *phs = re==0.0 && im==0.0 ? 0.0: atan2(im, re);
  } else {
    if(amp)
      *amp = vis->amp;
    if(phs)
      *phs = vis->phs;
  };
}

/*.......................................................................
 * Given the UT timestamp of an integration, return the time in the
 * timescale specified by vp->doutc.
 *
 * Input:
 *  vp       Vedpar *  The resource object of this module.
 *  ut       double    An integration UTC.
 * Output:
 *  return    float    The time ready for plotting.
 */
static float v_time(Vedpar *vp, double ut)
{
  if(vp->doutc)
    return ut - vp->utref;
  else
    return fmod(vp->stref + (ut - vp->utref) * ut_to_mst, daysec);
}

/*.......................................................................
 * Update the array of times, and the associated list of scans.
 *
 * Input:
 *  vp     Vedpar *    The resourced object of the program.
 * Output:
 *  return  int      0 - OK.
 *                   1 - Error.
 */
int v_update_times(Vedpar *vp)
{
  return v_get_times(vp) || new_Scans(vp)==NULL;
}

/*.......................................................................
 * Update the list of scans.
 *
 * Input:
 *  vp     Vedpar *    The resourced object of the program.
 * Output:
 *  return  int      0 - OK.
 *                   1 - Error.
 */
int v_update_scans(Vedpar *vp)
{
  return new_Scans(vp)==NULL;
}

/*.......................................................................
 * Toggle between displaying UTC and sidereal time along the X-axis.
 *
 * Input:
 *  vp    Vedpar *   The resource object of the program.
 * Output:
 *  return int     0 - OK.
 *                 1 - Error.
 */
int v_toggle_timesys(Vedpar *vp)
{
/*
 * Switch time systems.
 */
  vp->doutc = !vp->doutc;
/*
 * Has a subarray been selected yet?
 */
  if(vp->sub) {
/*
 * Arrange to start displaying the full time range.
 */
    vp->ta = 0;
    vp->tb = vp->sub->ntime - 1;
/*
 * Update the array of times and the list of scans.
 */
    return v_update_times(vp);
  };
  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "obs.h"
#include "vlbconst.h"
#include "vlbutil.h"
#include "vlbfft.h"
#include "telspec.h"
#include "baselist.h"
#include "pollist.h"
#include "spectra.h"
#include "enumpar.h"
#include "units.h"
#include "specplot.h"
#include "logio.h"
#include "cpgplot.h"

static const float ymarg=0.1;   /* The fraction of the Y range for margin */
static const float phsfrc=0.3f; /* Fraction of amp+phase plot devoted to phase*/
static const int datcol=10;     /* The color of unflagged data points */
static const int datsym=1;      /* The marker of good points */
static const int zoomcol=5;     /* PGPLOT color index for zoom cursor window */
static const float labsep=1.3f; /* Separation between label lines (characters)*/
static const float lmarg=3.0f;  /* Left axis margin (characters) */
static const float rmarg=0.5f;  /* Right axis margin (characters) */
static const float bmarg=3.0f;  /* Bottom axis margin (characters) */
static const float tmarg=0.5f;  /* Top axis margin (characters) */
static const float nsigma=4.5f; /* The smoothing function width wrt HWHM */ 

/* Define macros used to translate between X-axis coordinates and channels */

#define CHAN_TO_X(spx, chan) ((spx)->xoff + (spx)->xmul * (chan))
#define X_TO_CHAN(spx, x) (((x) - (spx)->xoff)/(spx)->xmul)

/* List distinguishing attributes of a sub-plot */

typedef struct {
  int uta,utb;    /*  Time: Start/end ob->rec[] integration indexes */
  int isub,base;  /* Split: Sub-array and baseline being plotted */
  int isel;       /*  Base: Index of the baseline selection in sp->sa->bgl */
  int ipol;       /*   Pol: Index of the polarization in sp->sa->pols */
  int iuv;        /*   UVR: UV range iteration index */
} SpAttr;

/* Spectrum sub-plot descriptor */

typedef struct {
  Spectrum *spec;     /* The spectrum being plotted */
  float vya, vyb;     /* Y axis limits of viewport (NDC) */
  float vymid;        /* The Y axis boundary between amp and phase plots */
  float amin, amax;   /* The range of amplitudes in the plot */
  float pmin, pmax;   /* The range of phases in the plot (degrees) */
  SpAttr spa;         /* Distinguishing characteristics of the plot */
} SpSubplot;

/* IF-specific x-axis descriptors - shared between sub-plots */

typedef struct {
  int doplot;       /* True if this IF is to be plotted */
  int slot;         /* The sequential position of this IF on the display */
  int cmin,cmax;    /* Index range of available channels wrt other IFs */
  float vxa, vxb;   /* X axis limits of viewport (NDC) */
  float xoff,xmul;  /* X-axis world-coord = xoff + chan * xmul */
  float xmin,xmax;  /* X-axis world-coordinate limits */
  int ca, cb;       /* Range of channels to be plotted */
} SpXdim;

/* Declare a container for cursor selection details */

typedef struct {
  int key;       /* The upper-case value of the key used to return the cursor */
  int waslow;    /* True if 'key' was lower case */
  int wasamp;    /* True if the cursor was pressed in an amplitude plot */
  int iplot;     /* The index of the selected sub-plot in sp->splots[] */
  int cif;       /* The index of the selected IF */
  float x;       /* The X-axis value selected (User X coordinates) */
  float y;       /* The amplitude or phase (degrees) selected */
} SpCurs;

/* Declare a container describing plot attributes */

typedef struct {
  Observation *ob;    /* The observation being plotted */
  SpCurs cursor;      /* Cursor selection descriptor */
  Specattr *sa;       /* Specplot attributes */
  float vxa, vxb;     /* X axis limits of the area enclosing all plots (NDC) */
  float vya, vyb;     /* Y axis limits of the area enclosing all plots (NDC) */
  float xch;          /* X-axis character height for sub-plot labelling */
  float ych;          /* Y-axis character height for sub-plot labelling */
  int docurs;         /* If true, allow cursor interaction */
  int nif;            /* The number of IFs in the observation. */
  int nifplot;        /* The number of IFs to be plotted */
  int nslot;          /* Requested number of plots per page */
  int nplot;          /* The actual number of plots on the page */
  int npage;          /* The ordinal number of the current page */
  Spectra *spectra;   /* The list of allocated spectra */
  Basegrp *scr_bgrp;  /* Scratch baseline selection list container */
  SpXdim *spx;        /* Array of IF-specific shared X-axis plot dimensions */
  char *xlabel;       /* The X-axis label string */
  int nsplot;         /* The allocated number of elements in splots[] */
  SpSubplot *splots;  /* Sub-plot descriptors for the current page */
  float *ywork;       /* A work array for plotted re,im or amp,phase pairs */
  float *wwork;       /* A work array for plotted weights */
  int nwork;          /* The number of elements in ywork[],wwork[] and work[] */
  float *smfn;        /* Half-even array of nsmth smoothing function values */
  int nchan;          /* The max number of channels in any spectrum */
  int nsmth;          /* The number of channels corresponding to nsigma*/
} Specplot;

/* 
 * List key bindings. If you change these, check whether s_flags(),
 * s_interract() and sp_set_flags() need changing.
 */
typedef enum {
  KEY_NONE = '\0', /* No cursor selection has been made */
  KEY_CROSS= '+',  /* Toggle the cross-hair cursor */
  KEY_CUR  = 'A',  /* Key for cursor position input */
  KEY_CAN  = 'D',  /* Key to cancel incomplete select range */
  KEY_ERR  = 'E',  /* Toggle plotting error bars */
  KEY_HELP = 'H',  /* List key bindings */
  KEY_DISP = 'L',  /* Re-display the current plot */
  KEY_NEXT = 'N',  /* Display the next page of spectra */
  KEY_ORDER= 'O',  /* Change the selection sort-order */
  KEY_PREV = 'P',  /* Display the previous page of spectra */
  KEY_NUMB = 'S',  /* Change the number of plots per page */
  KEY_SPEC = 'T',  /* Change indexing keys */
  KEY_QUIT = 'X',  /* End the plot session */
  KEY_JOIN = 'J',  /* Toggle the drawing of lines between data-points */
  KEY_XAXIS= 'U',  /* Delimit a new channel range with the cursor */
  KEY_VECT = 'V',  /* Toggle vector vs. scalar averaging */
  KEY_ZOOM = 'Z',  /* Change the plotted amplitude and/or phase range */
  KEY_AMP  = '1',  /* Request just amplitude plots */
  KEY_PHS  = '2',  /* Request just phase plots */
  KEY_BOTH = '3'   /* Request both amplitude and phase plots */
} Spkey;

static int s_flags(Specattr *sa, int key, int waslow);

/* New page plotting options */

typedef enum {
  S_ALLNEW,   /* Initialize for a new set of plots. */
  S_RESET,    /* Re-construct and plot the current spectra. */
  S_NEXT      /* Plot the next page in sequence. */
} SpNext;

static int s_page(Specplot *sp, SpNext oper, int forward);

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

static int s_cursor(Specplot *sp, Bandmode mode, int isamp, int iplot,
		    int cif, float xref, float yref, int ci);

/* Iterators */

typedef struct {
  SpAttr spa;
  Stokes pol;
  Basegrp *bgrp;
  int uta, utb;
  float uvmin, uvmax;
} Specposn;

static int s_iterate(Specplot *sp, SpKey key, int reset, int advance,
		     int forward, Specposn *spp);
static Basegrp *s_base_iter(Specplot *sp, int reset, int advance, int forward,
			    SpAttr *spa);
static int s_time_iter(Specplot *sp, int reset, int advance, int forward,
		       SpAttr *spa);
static Stokes s_pol_iter(Specplot *sp, int reset, int advance, int forward,
			 SpAttr *spa);
static int s_uvr_iter(Specplot *sp, int reset, int advance, int forward,
			 SpAttr *spa, float *uvmin, float *uvmax);

static Specplot *new_Specplot(Observation *ob, Specattr *sa, int docurs);
static Specplot *del_Specplot(Specplot *sp);
static int s_redisp(Specplot *sp);
static int s_getspec(Specplot *sp, int iplot, int cif);
static int s_yrange(Specplot *sp, int iplot);
static int s_vpwin(Specplot *sp);
static int s_setwin(Specplot *sp, int iplot, int cif, int doamp);
static int s_xrange(Specplot *sp);
static int s_plaxes(Specplot *sp, int iplot, int cif, int erase);
static int s_label(Specplot *sp);
static int s_plamp(Specplot *sp, int iplot, int cif, int erase);
static int s_plphs(Specplot *sp, int iplot, int cif, int erase);
static int s_label(Specplot *sp);
static int s_interact(Specplot *sp);
static int s_auto(Specplot *sp, int npage);
static int s_get_smfn(Specplot *sp, int cif);
static int bad_Specattr(Specattr *sa, char *fname);
static int s_get_xrange(Specplot *sp);
static int s_get_yrange(Specplot *sp);
static int s_get_xaxis(Specplot *sp);
static int s_get_smooth(Specplot *sp);
static int s_get_order(Specplot *sp);
static char *s_get_arg(char *string, char **endp);
static char *s_getline(char *prompt);
static int s_get_sel(Specplot *sp);
static int s_get_bgl(Specplot *sp, char *args);
static int s_get_pol(Specplot *sp, char *args);
static int s_get_times(Specplot *sp, char *args);
static int s_get_uvr(Specplot *sp, char *args);
static int s_newnum(Specplot *sp);
static SpSubplot *new_SpSubplot(Specplot *sp, int nnew);

static int s_title(Specplot *sp, SpKey key, SpAttr *spa, int nc, char *s);
static int s_bgrp_title(Specplot *sp, SpAttr *spa, int nc, char *s);
static int s_pol_title(Specplot *sp, SpAttr *spa, int nc, char *s);
static int s_time_title(Specplot *sp, SpAttr *spa, int nc, char *s);
static int s_uvr_title(Specplot *sp, SpAttr *spa, int nc, char *s);

static int s_coords(Specplot *sp, int cif, SpXunit xunit,
		    float *xoff, float *xmul);

/*.......................................................................
 * Plot spectra for a given spectral-line observation.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation.
 *  sa      Specattr *  A container of desired specplot attributes,
 *                      allocated and default initialized via
 *                      new_Specattr() and modified via sp_*() functions.
 *  docurs       int    If false disallow interactive plotting.
 *  npage        int    The maximum number of pages to plot when plotting
 *                      to a non-interactive device. A value of 0 is
 *                      means no limit.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int specplot(Observation *ob, Specattr *sa, int docurs, int npage)
{
  Specplot *sp;  /* The spectrum plot descriptor */
  int iret;      /* Return code */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "specplot") ||
     bad_Specattr(sa, "specplot"))
    return 1;
/*
 * Interpret the argument and return a descriptor of the spectrum plot.
 */
  sp = new_Specplot(ob, sa, docurs);
  if(!sp)
    return 1;
/*
 * Handle interactive and hardcopy plotting?
 */
  iret = sp->docurs ? s_interact(sp) : s_auto(sp, npage);
/*
 * Delete the descriptor.
 */
  sp = del_Specplot(sp);
  return iret;
}

/*.......................................................................
 * Create and initialize a specplot state descriptor.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation.
 *  sa      Specattr *  A container of desired specplot attributes,
 *                      allocated and default initialized via
 *                      new_Specattr() and modified via sp_*() functions.
 *  docurs       int    If false disallow interactive plotting.
 * Output:
 *  return       int    The new descriptor, or NULL on error.
 */
static Specplot *new_Specplot(Observation *ob, Specattr *sa, int docurs)
{
  Specplot *sp;  /* The container to be returned */
  char answer[10]; /* Buffer for PGQINF() replies */
  int slen;        /* The input/output length of a PGQINF() answer string */
  int i;
/*
 * Is there a plot device open?
 */
  slen = sizeof(answer)-1;
  cpgqinf("STATE", answer, &slen);
  if(strncmp(answer, "OPEN", 4) != 0) {
    lprintf(stderr, "specplot: No plot device open.\n");
    return NULL;
  };
/*
 * Make sure that we have at least one baseline selection.
 */
  if((!sa->bgl || sa->bgl->nsel<1) && sp_set_bgl(ob, sa, SP_SPLIT, NULL))
    return NULL;
/*
 * Allocate the descriptor container.
 */
  sp = (Specplot *) malloc(sizeof(Specplot));
  if(!sp) {
    lprintf(stderr, "specplot: Insufficient memory.\n");
    return NULL;
  };
/*
 * Before doing anything that might fail, initialize the descriptor at least
 * to the point at which it can safely be passed to del_Specplot() for
 * deletion.
 */
  sp->ob = ob;
  sp->cursor.key = KEY_NONE;
  sp->sa = sa;
  sp->vxa = 0.0f;
  sp->vxb = 1.0f;
  sp->vya = 0.0f;
  sp->vyb = 1.0f;
  sp->xch = 1.0f;
  sp->ych = 1.0f;
  sp->docurs = docurs;
  sp->nif = ob->nif;
  sp->nifplot = 0;
  sp->nslot = sp->sa->nplot <= 0 ? 3 : sp->sa->nplot;
  sp->nplot = 0;
  sp->npage = 0;
  sp->spectra = NULL;
  sp->scr_bgrp = NULL;
  sp->spx = NULL;
  sp->xlabel = "";
  sp->nsplot = 0;
  sp->splots = NULL;
  sp->ywork = NULL;
  sp->wwork = NULL;
  sp->nwork = 0;
  sp->smfn = NULL;
  sp->nchan = 0;
  sp->nsmth = 0;
/*
 * If cursor interaction is required, check if the device has a cursor.
 */
  slen = sizeof(answer)-1;
  cpgqinf("CURSOR", answer, &slen);
  sp->docurs = strncmp(answer, "YES", 3) == 0;
/*
 * Allocate a spectrum list container.
 */
  sp->spectra = new_Spectra(ob);
  if(!sp->spectra)
    return del_Specplot(sp);
/*
 * Allocate a scratch baseline group.
 */
  sp->scr_bgrp = new_Basegrp();
  if(!sp->scr_bgrp)
    return del_Specplot(sp);
/*
 * Allocate an array of sp->nif X-axis descriptors.
 */
  sp->spx = (SpXdim *) malloc(sizeof(SpXdim) * sp->nif);
  if(!sp->spx) {
    lprintf(stderr, "specplot: Insufficient memory.\n");
    return del_Specplot(sp);
  };
/*
 * Initialize each descriptor.
 */
  for(i=0; i<sp->nif; i++) {
    SpXdim *spx = sp->spx + i;
    spx->doplot = 1;
    spx->slot = i;
    spx->cmin = ob->ifs[i].coff;
    spx->cmax = spx->cmin + ob->nchan - 1;
    spx->vxa = spx->vxb = 0.0f;
    spx->xoff = 0.0f;
    spx->xmul = 1.0f;
    spx->xmin = spx->xmax = 0.0f;
    spx->ca = 0;
    spx->cb = ob->nchan-1;
  };
/*
 * Allocate sp->nslot sub-plot descriptors.
 */
  if(new_SpSubplot(sp, sp->nslot) == NULL)
    return del_Specplot(sp);
/*
 * Allocate two scratch arrays to use when constructing plot spectra from
 * raw spectra. For each element we need two float's for complex spectra.
 * Round the number of elements up to the next positive 
 * finite power of two at or above ob->nchan, so that we can use FFTs to
 * convert to cross-correlation spectra.
 */
  {
    unsigned int n;
    unsigned int ipow = 0;
    for(n = ob->nchan-1; n>0; n >>= 1U)
      ipow++;
    sp->nwork = 2 * (1U << ipow);
  };
  sp->ywork = (float *) malloc(sizeof(float) * sp->nwork);
  sp->wwork = (float *) malloc(sizeof(float) * sp->nwork);
  if(!sp->ywork || !sp->wwork) {
    lprintf(stderr, "specplot: Insufficient memory for scratch spectrum.\n");
    return del_Specplot(sp);
  };
/*
 * Allocate an array in which to cache the upper half of the current
 * smoothing function.
 */
  sp->nchan = ob->nchan;
  sp->smfn = (float *) malloc(sizeof(float) * sp->nchan);
  if(!sp->smfn) {
    lprintf(stderr, "specplot: Insufficient memory for smoothing function.\n");
    return del_Specplot(sp);
  };
  return sp;
}

/*.......................................................................
 * Delete a Specplot descriptor and its contents.
 *
 * Input:
 *  sp     Specplot *  The descriptor to be deleted.
 * Output:
 *  return Specplot *  Allways NULL.
 */
static Specplot *del_Specplot(Specplot *sp)
{
  if(sp) {
    sp->spectra = del_Spectra(sp->spectra);
    sp->scr_bgrp = del_Basegrp(sp->scr_bgrp, NULL);
    if(sp->spx)
      free(sp->spx);
    if(sp->splots)
      free(sp->splots);
    if(sp->ywork)
      free(sp->ywork);
    if(sp->wwork)
      free(sp->wwork);
    if(sp->smfn)
      free(sp->smfn);
    free(sp);
  };
  return sp;
}

/*.......................................................................
 * Toggle plotting flags given a command key.
 *
 * Input:
 *  sa  Specattr *  The specplot attributes descriptor.
 *  key      int    The command key (upper case).
 *  waslow   int    If true the key was lower case.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Key unrecognized.
 */
static int s_flags(Specattr *sa, int key, int waslow)
{
  switch (key) {
  case KEY_AMP:
    sa->doamp = 1;
    sa->dophs = 0;
    break;
  case KEY_PHS:
    sa->doamp = 0;
    sa->dophs = 1;
    break;
  case KEY_BOTH:
    sa->doamp = sa->dophs = 1;
    break;
  case KEY_JOIN:
    if(waslow)
      sa->dojoin = !sa->dojoin;
    else
      sa->dohist = !sa->dohist;
    break;
  case KEY_ERR:
    sa->dobars = !sa->dobars;
    break;
  default:
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Start the next page of plots.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  oper    SpNext     The action to take to display the next page.
 *                      S_ALLNEW - Initialize for a new set of plots.
 *                      S_RESET  - Re-construct and plot the current
 *                                 spectra using revised attributes.
 *                      S_NEXT   - Plot the next page in sequence.
 *  forward    int     The ordinal direction to plot baselines from.
 * Output:
 *  return     int     The number of sub-plots plotted, or -1 on error.
 */
static int s_page(Specplot *sp, SpNext oper, int forward)
{
  Specattr *sa;   /* Specplot attributes container */
  Specposn spp;   /* Specplot attribute descriptors */
  SpKey *key;     /* Array of iterator key types */
  int i;
/*
 * Get the specplot attributes descriptor.
 */
  sa = sp->sa;
  key = sa->key;
/*
 * The first page must be plotted with S_ALLNEW.
 */
  if(sp->npage < 1)
    oper = S_ALLNEW;
/*
 * Initialize the position descriptor.
 */
  spp.spa.uta = spp.spa.utb = 0;
  spp.spa.isub = spp.spa.base = 0;
  spp.spa.isel = 0;
  spp.spa.ipol = 0;
  spp.spa.iuv = 0;
  spp.pol = NO_POL;
  spp.bgrp = NULL;
  spp.uta = spp.utb = 0;
  spp.uvmin = spp.uvmax = 0.0f;
/*
 * Locate the first of the requested plots.
 */
  switch(oper) {
  case S_ALLNEW:
    for(i=0; i<SP_NKEY; i++) {
      if(s_iterate(sp, key[i], 1, 0, 1, &spp)) {
	printf("No spectra selected.\n");
	return -1;
      };
    };
    break;
  case S_RESET:
    spp.spa = sp->splots[0].spa;
    forward = 1;
    for(i=0; i<SP_NKEY; i++) {
      if(s_iterate(sp, key[i], 0, 0, 1, &spp)) {
	printf("No spectra selected.\n");
	return -1;
      };
    };
    break;
  case S_NEXT:
    {
      int found = 0;  /* True after a succesful increment */
/*
 * Locate the plot following the last existing plot in the given direction.
 */
      spp.spa = sp->splots[forward ? sp->nplot-1 : 0].spa;
/*
 * Attempt to increment each of the variable keys until we run out of
 * keys, or an increment is succesful. Reset unincrementable keys
 * to the appropriate end of their ranges before trying to increment
 * the next key.
 */
      for(i=0; !found && i<sa->nkey; i++) {
	if(s_iterate(sp, key[i], 0, 1, forward, &spp)==0)
	  found = 1;
	else if(s_iterate(sp, key[i], 1, 0, forward, &spp)) {
	  lprintf(stderr, "s_page: Unable to reset search key.\n");
	  return -1;
	};
      };
/*
 * All increments failed?
 */
      if(!found) {
	printf("No spectra remain to be plotted.\n");
	return 0;
      };
/*
 * Acquire the details of the current positions of the rest of the keys
 * without.
 */
      for( ; i<SP_NKEY; i++) {
	if(s_iterate(sp, key[i], 0, 0, forward, &spp)) {
	  lprintf(stderr, "s_page: Error locating search key.\n");
	  return -1;
	};
      };
    };
    break;
  default:
    lprintf(stderr, "s_page: Unrecognised SpNext code.\n");
    return -1;
    break;
  };
/*
 * Install up to sp->nslot spectra.
 */
  sp->nplot = 0;
  do {
    SpSubplot *sps = sp->splots + sp->nplot;
/*
 * Set up the next spectrum container.
 */
    if(sps->spec) {
      if(spc_set_bgrp(sp->spectra, sps->spec, spp.bgrp) ||
	 spc_set_pol(sp->spectra, sps->spec, spp.pol) ||
	 spc_set_ut(sp->spectra, sps->spec, spp.uta, spp.utb) ||
	 spc_set_uvrange(sp->spectra, sps->spec, spp.uvmin, spp.uvmax) ||
	 spc_set_avmode(sp->spectra, sps->spec, sa->avmode==SP_VECTOR))
	return -1;
    } else {
      sps->spec = add_Spectrum(sp->spectra, sa->avmode==SP_VECTOR,
			       spp.pol, spp.uta, spp.utb,
			       spp.uvmin, spp.uvmax, spp.bgrp);
      if(!sps->spec)
	return -1;
    };
/*
 * Record details of the new spectrum.
 */
    sps->spa = spp.spa;
/*
 * Is there a new sub-plot to be plotted.
 */
  } while(++sp->nplot < sp->nslot &&
	  s_iterate(sp, key[0], 0, 1, forward, &spp)==0);
/*
 * When plotting a previous page, the spectra sub-plots end up in reverse
 * order. Put them back into forward order.
 */
  if(!forward) {
    SpSubplot *spa = sp->splots;
    SpSubplot *spb = sp->splots + sp->nplot-1;
    for( ; spa<spb; spa++,spb--) {
      SpSubplot sptmp = *spa;
      *spa = *spb;
      *spb = sptmp;
    };
  };
/*
 * Read the spectra.
 */
  if(get_Spectra(sp->spectra))
    return -1;
/*
 * Display the new spectra.
 */
  if(s_redisp(sp))
    return -1;
/*
 * Return the number of plots plotted.
 */
  return sp->nplot;
}

/*.......................................................................
 * Display the sp->nplot spectra currently recorded in sp->splots[].
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 * Output:
 *  return     int     0 - OK.
 */
static int s_redisp(Specplot *sp)
{
  int iplot;      /* The index of the sub-plot being processed */
  int ierr = 0;   /* Non-zero if an error occurs */
/*
 * Nothing to plot?
 */
  if(sp->nplot<=0) {
    lprintf(stderr, "s_redisp: No plot rows have been initialized.\n");
    return 1;
  };
/*
 * Clear page.
 */
  cpgpage();
/*
 * Count pages.
 */
  sp->npage++;
/*
 * Set the X-axis world-coordinate ranges.
 */
  ierr = ierr || s_xrange(sp);
/*
 * Determine viewports for each sub-plot.
 */
  ierr = ierr || s_vpwin(sp);
/*
 * Plot each sub-plot.
 */
  for(iplot=0; iplot<sp->nplot && !ierr; iplot++) {
    int cif;
    cpgbbuf();
    ierr = ierr || s_yrange(sp, iplot);        /* Get the sub-plot data range */
    for(cif=0; cif<sp->nif; cif++) {
      if(sp->spx[cif].doplot) {
	ierr = ierr || s_getspec(sp, iplot, cif);  /* Get sp->ywork,sp->wwork */
	ierr = ierr || s_plaxes(sp, iplot, cif, 0);/* Plot axes */
	ierr = ierr || s_plamp(sp, iplot, cif, 0); /* Plot amplitude */
	ierr = ierr || s_plphs(sp, iplot, cif, 0); /* Plot phase */
      };
    };
    ierr = ierr || (iplot==0 && s_label(sp));     /* Label around the plot */
    cpgebuf();
  };
  return ierr;
}

/*.......................................................................
 * Convert the complex spectrum of a given IF into the form required
 * for plotting. This includes application of optional smoothing, and
 * conversion from complex to polar coordinates. The final spectrum
 * is placed in sp->ywork[] as pairs of amplitude and phase values. The
 * final weights are placed in sp->wwork[] in the first element of each
 * pair of elements parallel to sp->ywork[].
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  iplot      int     The index of the sub-plot in sp->splots[] whose
 *                     spectrum is to be prepared.
 *  cif        int     The index of the IF within the sub-plot.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int s_getspec(Specplot *sp, int iplot, int cif)
{
  SpXdim *spx = sp->spx + cif;
  SpSubplot *sps = sp->splots + iplot;
  int nchan;       /* The number of channels in the spectrum */
  int i;
/*
 * Process this IF?
 */
  if(!spx->doplot)
    return 0;
/*
 * How many channels are there in the given IF?
 */
  nchan = sps->spec->ifs[cif].nchan;
/*
 * Copy the complex spectrum without smoothing.
 */
  if(sp->sa->smooth.type == SM_NONE) {
    Cvis *chan = sps->spec->ifs[cif].chan;
    for(i=0; i<nchan; i++) {
      sp->ywork[2*i] = chan[i].re;
      sp->ywork[2*i+1] = chan[i].im;
      sp->wwork[2*i] = chan[i].wt;
      sp->wwork[2*i+1] = 0.0f;
    };
/*
 * Construct a smoothed version of the spectrum.
 */
  } else {
    Cvis *chan = sps->spec->ifs[cif].chan;
    float *smfn = sp->smfn;
/*
 * Initialize the smoothing function.
 */
    if(s_get_smfn(sp, cif))
      return 1;
/*
 * Calculate smoothed real and imaginary parts and the resulting weights
 * for each spectrum channel.
 */
    for(i=0; i<nchan; i++,chan++) {
      float re_w_s = 0.0f; /* Sum of real.weight.smfn */
      float im_w_s = 0.0f; /* Sum of imag.weight.smfn */
      float w_s = 0.0f;    /* Sum of weight.smfn */
      float w_ss = 0.0f;   /* Sum of weight.smfn.smfn */
      int jlim;            /* The max index to use in smfn[] */
      Cvis *cvis;          /* Pointer to the next channel to be included */
      int j;               /* Index into smfn[] */
/*
 * Get the contribution to the convolutions from the current
 * output channel and the channels below it.
 */
      jlim = sp->nsmth-1;
      if(i-jlim < 0)
	jlim = i;
      for(j=jlim,cvis=chan-j; j>=0; j--,cvis++) {
	float s = smfn[j];
	float ws = cvis->wt * s;
	re_w_s += cvis->re * ws;
	im_w_s += cvis->im * ws;
	w_s += ws;
	w_ss += ws * s;
      };
/*
 * Get the contribution to the convolutions from channels above the
 * current output channel.
 */
      jlim = sp->nsmth-1;
      if(i+jlim >= nchan)
	jlim = nchan-i-1;
      for(j=1,cvis=chan+j; j<=jlim; j++,cvis++) {
	float s = smfn[j];
	float ws = cvis->wt * s;
	re_w_s += cvis->re * ws;
	im_w_s += cvis->im * ws;
	w_s += ws;
	w_ss += ws * s;
      };
/*
 * Get the output weight and data values.
 */
      {
	float wt = w_ss != 0.0f ? w_s * w_s / w_ss : 0.0f;
	if(wt < 0.0f) wt = 0.0f;
	sp->wwork[2*i]   = wt;
	sp->ywork[2*i]   = wt > 0.0f ? re_w_s / w_s : 0.0f;
	sp->ywork[2*i+1] = wt > 0.0f ? im_w_s / w_s : 0.0f;
      };
    };
  };
/*
 * Zero-pad any remaining elements of the work array. This makes it
 * it suitable for processing with an FFT.
 */
  for(i=2*nchan; i<sp->nwork; i++) {
    sp->ywork[i] = 0.0f;
    sp->wwork[i] = 0.0f;
  };
/*
 * Convert complex values to amplitude and phase.
 */
  for(i=0; i<2*nchan; i += 2) {
    float re = sp->ywork[i];
    float im = sp->ywork[i+1];
    if(sp->sa->avmode==SP_VECTOR) {
      sp->ywork[i] = sqrt(re*re + im*im);
      sp->ywork[i+1] = (re!=0.0f || im!=0.0f) ? atan2(im,re) : 0.0f;
    } else {
      sp->ywork[i] = re;
      sp->ywork[i+1] = im;
    };
  };
  return 0;
}

/*.......................................................................
 * Ascertain the amplitude and phase range of a given sub-plot.
 *
 * Input:
 *  sp   Specplot *   The spectrum plot descriptor.
 *  iplot     int     The index of the sub-plot in sp->splots[] whose
 *                    spectrum is to be processed.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
static int s_yrange(Specplot *sp, int iplot)
{
  SpSubplot *sps =  sp->splots + iplot; /* The descriptor of the sub-plot */
  float amin,amax; /* The amplitude range of the spectra */
  float pmin,pmax; /* The phase range of the spectra */
  int first=1;     /* True until amin and amax have been first defined */
  int cif;         /* The index of an IF */
  int i;
/*
 * Amplitude.
 * Get the optional fixed limits.
 */
  amin = sp->sa->amin;
  amax = sp->sa->amax;
/*
 * If no valid limits have been specified, autoscale.
 */
  if(amin >= amax) {
/*
 * Determine the overall amplitude range from the spectra of all IFs.
 */
    for(cif=0; cif<sps->spec->nif; cif++) {
      SpXdim *spx = sp->spx + cif;
      if(spx->doplot) {
	if(s_getspec(sp, iplot, cif))
	  return 1;
/*
 * Determine the range of spectrum amplitudes in sp->ywork[].
 */
	for(i=spx->ca; i<=spx->cb; i++) {
	  float amp = sp->ywork[2*i];
	  if(first) {
	    first = 0;
	    amax = amin = amp;
	  } else if(amp > amax) {
	    amax = amp;
	  } else if(amp < amin) {
	    amin = amp;
	  };
	};
      };
    };
/*
 * Make sure that 0 amplitude appears in the plot.
 */
    if(amin > 0.0f)
      amin = 0.0f;
    if(amax < 0.0f)
      amax = 0.0f;
/*
 * Cater for the case where spectra contain no data.
 */
    if(amin==0.0f && amax==0.0f)
      amax = 1.0f;
/*
 * Leave margins around the plots.
 */
    {
      float margin = ymarg * (amin >= amax ? fabs(amax) : fabs(amax-amin));
      amin -= margin;
      amax += margin;
    };
  };
/*
 * Phase.
 * Get the optional fixed limits.
 */
  pmin = sp->sa->pmin;
  pmax = sp->sa->pmax;
/*
 * If no valid limits have been specified, autoscale.
 */
  if(pmin >= pmax) {
    pmin = -180.0f;
    pmax = 180.0f;
  };
/*
 * Record the new amplitude and phase ranges.
 */
  sps->amin = amin;
  sps->amax = amax;
  sps->pmin = pmin;
  sps->pmax = pmax;
  return 0;
}

/*.......................................................................
 * Determine the viewport coordinates of all sub-plots.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int s_vpwin(Specplot *sp)
{
  float hch;  /* The NDC height of unity character height horizontal text */
  float vch;  /* The NDC height of unity character height vertical text */
  float hgap; /* The NDC horizontal gap between adjacent IF plots */
  float vgap; /* The NDC vertical gap above each sub-plot */
  float xwidth; /* The sum of the widths of each plot in world coordinates */
  int i;
/*
 * Determine the NDC height of text characters.
 */
  cpgsch(1.0);
  cpgqcs(0, &hch, &vch);
/*
 * Delimit the part of the viewport that encloses all plots, exluding
 * all labelling.
 */
  sp->vxa = hch * (labsep + lmarg);
  sp->vxb = 1.0f - hch * rmarg;
  sp->vya = vch * (labsep + bmarg);
  sp->vyb = 1.0f - vch * (SP_NKEY * labsep + tmarg);
/*
 * Count the number of IFs that are to be plotted and the total
 * world-coordinate width involved.
 */
  sp->nifplot = 0;
  xwidth = 0.0f;
  for(i=0; i<sp->nif; i++) {
    SpXdim *spx = sp->spx + i;
    if(spx->doplot) {
      spx->slot = sp->nifplot++;
      xwidth += fabs(spx->xmax - spx->xmin);
    };
  };
/*
 * Ensure that we have something to plot.
 */
  if(sp->nifplot < 1 || xwidth <= 0.0f) {
    lprintf(stderr, "s_vpwin: No IFs selected for plotting.\n");
    return 1;
  };
/*
 * Set a horizontal gap between horizontally adjacent IF sub-plots.
 */
  hgap = sp->nifplot < 2 ? 0.0 : 0.05 * (sp->vxb - sp->vxa) / (sp->nifplot - 1);
/*
 * Set the vertical gap above each sub-plot.
 */
  {
    float vfrac = 0.05 + 0.03 * (sp->nslot-1); /* Fractional size of sub-plot */
    if(vfrac > 0.25)
      vfrac = 0.25;
    vgap = (sp->vyb - sp->vya) * vfrac / sp->nslot;
  };
/*
 * Divide the X-axis between IFs.
 */
  {
    float hoff = sp->vxa;     /* Horizontal offset accumulated */
    float hwid = (sp->vxb - sp->vxa - (sp->nifplot-1) * hgap); /* Total width */
    for(i=0; i<sp->nif; i++) {
      SpXdim *spx = sp->spx + i;
      if(spx->doplot) {
	spx->vxa = hoff;
	hoff += hwid * fabs(spx->xmax - spx->xmin) / xwidth;
	spx->vxb = hoff;
	hoff += hgap;
      } else {
	spx->vxa = spx->vxb = 0.0f;
      };
    };
  };
/*
 * Divide the Y-axis between sub-plots.
 */
  {
    float vtop = sp->vyb;        /* The top edge of the used plot area */
/*
 * Determine the height of each sub-plot.
 */
    float ysize = (sp->vyb - sp->vya) / sp->nslot - vgap;
/*
 * Work out the character height needed to comfortably fit text within
 * the gap left above each sub-plot.
 */
    sp->xch = 0.6f * vgap/vch;
    sp->ych = 0.8f * sp->xch;
/*
 * Assign Y-axis ranges to each sub-plot.
 */
    for(i=0; i<sp->nplot; i++) {
      SpSubplot *sps = sp->splots + i;
      vtop -= vgap;
      sps->vyb = vtop;
      vtop -= ysize;
      sps->vya = vtop;
      if(sp->sa->doamp && sp->sa->dophs)
	sps->vymid = sps->vya + phsfrc * (sps->vyb - sps->vya);
      else if(sp->sa->doamp)
	sps->vymid = sps->vya;
      else
	sps->vymid = sps->vyb;
    };
/*
 * Shrink the plot area to fit the actual number of plots to be plotted.
 */
    sp->vya = vtop;
  };
  return 0;
}

/*.......................................................................
 * Set the viewport around a given sub-plot and set up its world
 * coordinates.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  iplot      int     The index of the sub-plot in sp->splots[].
 *  cif        int     The index of the IF within the sub-plot.
 *  doamp      int     If true set up the ampitude plot, otherwise for
 *                     the phase plot.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int s_setwin(Specplot *sp, int iplot, int cif, int doamp)
{
  SpSubplot *sps = sp->splots + iplot;
  SpXdim *spx = sp->spx + cif;
/*
 * Process this IF?
 */
  if(!spx->doplot)
    return 0;
/*
 * Check that we are plotting the type of plot requested.
 */
  if((doamp && !sp->sa->doamp) || (!doamp && !sp->sa->dophs)) {
    lprintf(stderr, "s_setwin: Can't set %s viewport.\n",
	    doamp ? "amplitude":"phase");
    return 1;
  };
/*
 * Set up the viewport and world-coordinates.
 */
  if(doamp) {
    cpgsvp(spx->vxa, spx->vxb, sps->vymid, sps->vyb);
    cpgswin(spx->xmin, spx->xmax, sps->amin, sps->amax);
  } else {
    cpgsvp(spx->vxa, spx->vxb, sps->vya, sps->vymid);
    cpgswin(spx->xmin, spx->xmax, sps->pmin, sps->pmax);
  };
  return 0;
}

/*.......................................................................
 * Plot the axes of a given sub-plot.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  iplot      int     The index of the sub-plot in sp->splots[].
 *  cif        int     The index of the IF within the sub-plot.
 *  erase      int     If true, erase the axes, otherwise draw them anew.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int s_plaxes(Specplot *sp, int iplot, int cif, int erase)
{
  char awrk[81];        /* Work string for labelling */
  SpSubplot *sps = sp->splots + iplot;
  SpXdim *spx = sp->spx + cif;
  enum {PG_MAXOPT=15};
  char xopts[PG_MAXOPT];   /* Array of PGBOX X-axis options */
  char yopts[PG_MAXOPT];   /* Array of PGBOX Y-axis options */
  char *yopt;              /* Pointer into yopts[] */
  char *xopt;              /* Pointer into xopts[] */
/*
 * Process this IF?
 */
  if(!spx->doplot)
    return 0;
/*
 * Set the appropriate plot color.
 */
  cpgsci(erase ? 0 : 1);
/*
 * Set the PGBOX options for the Y axes.
 */
  yopt = yopts;
  *yopt++ = 'S';      /* Draw minor tick marks */
  *yopt++ = 'T';      /* Draw major tick marks */
  *yopt++ = 'B';      /* Draw left axis */
  *yopt++ = 'C';      /* Draw right axis */
  *yopt++ = 'V';      /* Draw labels perpendicular to axis */
  if(spx->slot==0)
    *yopt++ = 'N';    /* Draw numeric labels on the left-most plot */
  *yopt='\0';
/*
 * Set the PGBOX options for the X axes.
 */
  xopt = xopts;
  *xopt++ = 'S';      /* Draw minor tick marks */
  *xopt++ = 'T';      /* Draw major tick marks */
  *xopt++ = 'B';      /* Draw lower axis */
  *xopt++ = 'C';      /* Draw upper axis */
  if(iplot==sp->nplot-1)
    *xopt++ = 'N';    /* Draw numeric labels under the last plot */
  *xopt = '\0';
/*
 * Plot amplitude-plot Y-axes.
 */
  if(sp->sa->doamp) {
    if(s_setwin(sp, iplot, cif, 1))
      return 1;
    cpgsch(sp->ych);
    cpgbox("", 0.0f, 0, yopts, 0.0f, 0);
  };
/*
 * Plot phase-plot Y-axes.
 */
  if(sp->sa->dophs) {
    if(s_setwin(sp, iplot, cif, 0))
      return 1;
    cpgsch(sp->ych);
    cpgbox("", 0.0f, 0, yopts, 0.0f, 0);
  };
/*
 * Plot X-axes.
 */
  cpgsvp(spx->vxa, spx->vxb, sps->vya, sps->vyb);
  cpgswin(spx->xmin, spx->xmax, 0.0f, 1.0f);
  cpgsch(sp->xch);
  if(spx->ca == spx->cb) {
    cpgbox(xopts, 0.5f * (spx->xmin+spx->xmax), 1, "", 0.0f, 0);
  } else {
    cpgbox(xopts, 0.0f, 0, "", 0.0f, 0);
  };
/*
 * Draw a line separating the phase and amplitude plots.
 */
  if(sp->sa->doamp && sp->sa->dophs) {
    cpgmove(spx->xmin, phsfrc);
    cpgdraw(spx->xmax, phsfrc);
  };
/*
 * Label the sub-plot.
 */
  if(spx->slot==sp->nifplot-1) {
    if(s_title(sp, sp->sa->key[0], &sps->spa, sizeof(awrk), awrk) > 0)
      cpgmtxt("T", 0.4f, 1.0f, 1.0f, awrk);
  };
  if(iplot==0) {
    sprintf(awrk, "IF %d", cif+1);
    cpgmtxt("T", -1.5f, sp->nifplot * 0.01, 0.0f, awrk);
  };
/*
 * Restore the default plot color.
 */
  cpgsci(1);
  cpgsch(1.0f);
  return 0;
}

/*.......................................................................
 * Set up the X-axis display range of each IF and compute the offset
 * and scale factors required to convert from channel coordinates to
 * the user specified X-axis coordinates.
 *
 * Input:
 *  sp     Specplot *   The spectrum plot descriptor.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
static int s_xrange(Specplot *sp)
{
  SpXdim *spx = sp->spx; /* The array of sp->nif X-axis descriptors */
  int cif;               /* The index of the IF being processed */
/*
 * Set up each the X-axis range for each IF.
 */
  for(cif=0; cif<sp->nif; cif++,spx++) {
    spx->doplot = (sp->sa->ca <= spx->cmax && sp->sa->cb >= spx->cmin);
    if(spx->doplot) {
/*
 * Initialize spx->xoff and spx->xmul.
 */
      if(s_coords(sp, cif, sp->sa->xunit, &spx->xoff, &spx->xmul))
	return 1;
/*
 * Ascertain the plotted channel range of the current IF.
 */
      spx->ca = (sp->sa->ca < spx->cmin ? spx->cmin : sp->sa->ca) - spx->cmin;
      spx->cb = (sp->sa->cb > spx->cmax ? spx->cmax : sp->sa->cb) - spx->cmin;
/*
 * Get the world coordinate plot range of the current IF.
 */
      spx->xmin = CHAN_TO_X(spx, spx->ca);
      spx->xmax = CHAN_TO_X(spx, spx->cb);
/*
 * Add a half-channel margin to each end of the axis range.
 */
      {
	float margin = 0.5f * spx->xmul;
	spx->xmin -= margin;
	spx->xmax += margin;
      };
    };
  };
/*
 * Get a label for the X-axis coordinates.
 */      
  switch(sp->sa->xunit) {
  case SP_CHAN:
    sp->xlabel = "Channels";
    break;
  case SP_FREQ:
    sp->xlabel = "Frequency (GHz)";
    break;
  default:
    sp->xlabel = "?";
    lprintf(stderr, "s_xrange: X-axis type not recognized.\n");
    return 1;
    break;
  };
  return 0;
}

/*.......................................................................
 * Plot the amplitude spectrum in a given IF sub-plot.
 * The spectrum is assumed to have already been placed in sp->ywork[],
 * so it is imperitive that s_getspec() have been called prior to this
 * function so that the spectrum data.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  iplot      int     The index of the sub-plot in sp->splots[].
 *  cif        int     The index of the IF within the sub-plot.
 *  erase      int     If true, erase the data, otherwise draw them anew.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int s_plamp(Specplot *sp, int iplot, int cif, int erase)
{
  SpXdim *spx = sp->spx + cif;
  int started=0;     /* False to start a new vector */
  float old_x = 0.0f;/* The x-coordinate of the last point plotted */
  float old_y = 0.0f;/* The y-coordinate of the last point plotted */
  int dojoin;        /* True to draw vectors between adjacent lines */
  int dohist;        /* True if connecting lines are histogrammed */
  int i;
/*
 * Process this IF?
 */
  if(!spx->doplot)
    return 0;
/*
 * No amplitude plot?
 */
  if(!sp->sa->doamp)
    return 0;
/*
 * Set the viewport and world coordinates of the selected IF sub-plot.
 */
  if(s_setwin(sp, iplot, cif, 1))
    return 1;
/*
 * Set the data plot color.
 */
  cpgsci(erase ? 0 : datcol);
/*
 * Only attempt to plot vectors between points if there are at least two points.
 */
  dojoin = sp->sa->dojoin && spx->ca != spx->cb;
  dohist = sp->sa->dohist;
/*
 * Plot the data.
 */
  for(i=spx->ca; i<=spx->cb; i++) {
    float xa = CHAN_TO_X(spx, i-0.5);
    float xb = CHAN_TO_X(spx, i+0.5);
    float x =  (xa+xb)/2.0f;
    float y = sp->ywork[2*i];
    float wt = sp->wwork[2*i];
    if(wt > 0.0f) {
      if(dojoin) {
	if(started) {
	  if(dohist) {
	    cpgmove(xa, old_y);
	    cpgdraw(xa, y);
	    cpgdraw(xb, y);
	  } else {
	    cpgmove(old_x, old_y);
	    cpgdraw(x,y);
	  };
	} else {
	  started = 1;
	  if(dohist) {
	    cpgmove(xa, y);
	    cpgdraw(xb, y);
	  } else {
	    cpgpt(1, &x, &y, datsym);
	  };
	};
      } else {
	cpgpt(1, &x, &y, datsym);
      };
/*
 * Draw error bars?
 */
      if(sp->sa->dobars) {
	float amperr = 1.0f/sqrt(wt);
	cpgmove(x, y-amperr);
	cpgdraw(x, y+amperr);
      };
      old_x = x;
      old_y = y;
    } else {
      started = 0;
    };
  };
/*
 * Restore the default plot color.
 */
  cpgsci(1);
  return 0;
}

/*.......................................................................
 * Plot the phase spectrum in a given IF sub-plot.
 * The spectrum is assumed to have already been placed in sp->ywork[],
 * so it is imperitive that s_getspec() have been called prior to this
 * function so that the spectrum data.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  iplot      int     The index of the sub-plot in sp->splots[].
 *  cif        int     The index of the IF within the sub-plot.
 *  erase      int     If true, erase the data, otherwise draw them anew.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int s_plphs(Specplot *sp, int iplot, int cif, int erase)
{
  SpXdim *spx = sp->spx + cif;
  int dojoin;           /* True to draw vectors between adjacent lines */
  int dohist;           /* True if connecting lines are to be histogrammed */
  float old_x = 0.0f;   /* Previous X-coordinate associated with old_y */
  float old_y = 0.0f;   /* Previous plotted phase */
  int started = 0;      /* False to start a new vector */
  int i;
/*
 * Process this IF?
 */
  if(!spx->doplot)
    return 0;
/*
 * No phase plot?
 */
  if(!sp->sa->dophs)
    return 0;
/*
 * Set the viewport and world coordinates of the selected IF sub-plot.
 */
  if(s_setwin(sp, iplot, cif, 0))
    return 1;
/*
 * Set the data plot color.
 */
  cpgsci(erase ? 0 : datcol);
/*
 * Only attempt to plot vectors between points if there are at least two points.
 */
  dojoin = sp->sa->dojoin && spx->ca != spx->cb;
  dohist = sp->sa->dohist;
/*
 * Plot the phases.
 */
  for(i=spx->ca; i<=spx->cb; i++) {
    float xa = CHAN_TO_X(spx, i-0.5);
    float xb = CHAN_TO_X(spx, i+0.5);
    float x = (xa + xb)/2.0f;
    float y = rtod * sp->ywork[2*i+1];
    float wt = sp->wwork[2*i];
/*
 * Only plot good data.
 */
    if(wt > 0.0f) {
/*
 * Plot the data joined by lines. Given that phases at -180 and +180 are
 * equivalent, join phases that straddle this boundary by a distance shorter
 * than the distance within the plot, by a line that crosses both the
 * top and bottom axes.
 */
      if(dojoin) {
/*
 * Start the vector.
 */
	if(!started) {
	  started = 1;
	  if(dohist) {
	    cpgmove(xa, y);
	    cpgdraw(xb, y);
	  } else {
	    cpgpt(1, &x, &y, datsym);
	  };
/*
 * Continue a phase line.
 */
	} else if(dohist) {
	  float ydif = (y - old_y);  /* Phase offset between adjacent phases */
	  cpgmove(xa, old_y);
/*
 * Is it shorter to join the phases by a line that exits the bottom of the
 * phase plot and re-emerges from the top?
 */
	  if(ydif > 180.0f) {
	    cpgdraw(xa, y-360.0f);       /* Line exits at -180 */
	    cpgmove(xa, old_y+360.0f);   /* Line re-enters at +180 */
/*
 * Is it shorter to join the phases by a line that exits the top of the
 * phase plot and re-emerges from the bottom?
 */
	  } else if(ydif < -180.0f) {
	    cpgdraw(xa, y+360.0f);        /* Line exits at +180 */
	    cpgmove(xa, old_y-360.0f);    /* Line re-enters at -180 */
	  };
	  cpgdraw(xa, y);
	  cpgdraw(xb, y);
	} else {
	  float ydif = (y - old_y);  /* Phase offset between adjacent phases */
	  cpgmove(old_x, old_y);
/*
 * Is it shorter to join the phases by a line that exits the bottom of the
 * phase plot and re-emerges from the top?
 */
	  if(ydif > 180.0f) {
	    cpgdraw(x, y-360.0f);           /* Line exits at -180 */
	    cpgmove(old_x, old_y+360.0f);   /* Line re-enters at +180 */
/*
 * Is it shorter to join the phases by a line that exits the top of the
 * phase plot and re-emerges from the bottom?
 */
	  } else if(ydif < -180.0f) {
	    cpgdraw(x, y+360.0f);            /* Line exits at +180 */
	    cpgmove(old_x, old_y-360.0f);    /* Line re-enters at -180 */
	  };
	  cpgdraw(x, y);
	};
/*
 * Record the current X and Y coordinates so that they can be
 * connected to the next phase.
 */
	old_y = y;
	old_x = x;
/*
 * Plot each phase as an isolated point.
 */
      } else {
	cpgpt(1, &x, &y, datsym);
      };
/*
 * Draw error bars?
 */
      if(sp->sa->dobars) {
	float amp = sp->ywork[2*i];
	if(amp != 0.0f) {
	  float phserr = rtod/sqrt(wt)/amp;
	  cpgmove(x, y-phserr);
	  cpgdraw(x, y+phserr);
	};
      };
    } else {
      started = 0;
    };
  };
/*
 * Restore the default plot color.
 */
  cpgsci(1);
  return 0;
}

/*.......................................................................
 * Draw labels around the viewport that encloses all of the plots.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int s_label(Specplot *sp)
{
  Observation *ob;      /* The descriptor of the observation being plotted */
  Specattr *sa;         /* Specplot attributes descriptor */
  char awrk[81];        /* Work string for labelling */
  char bwrk[81];        /* Work string for labelling */
  int ititle = SP_NKEY; /* The number of next title line to be written */
  int i;
/*
 * Get the observation descriptor.
 */
  ob = sp->ob;
  sa = sp->sa;
/*
 * Set the viewport around the area that encloses the plot axes.
 */
  cpgsvp(sp->vxa, sp->vxb, sp->vya, sp->vyb);
/*
 * Compose and write the main title.
 */
  cpgsci(1);
  sprintf(awrk,"%s  %s  %s averaged spectra.", ob->source.name,
	  sutdate(ob->date.year, ob->date.ut, bwrk),
	  sp->sa->avmode==SP_VECTOR ? "Vector":"Scalar");
  cpgsch(1.0f);
  cpgmtxt("T", tmarg + labsep * --ititle, 0.0f, 0.0f, awrk);
/*
 * In non-interactive mode tell the user what is being plotted.
 */
  if(!sp->docurs)
    lprintf(stdout, "Page %02.2d.\n", sp->npage);
/*
 * Write extra title lines describing the features that are common
 * to all plots.
 */
  for(i=1; i<SP_NKEY; i++) {
    if(s_title(sp, sp->sa->key[i], &sp->splots[0].spa, sizeof(awrk), awrk)>0)
      cpgmtxt("T", tmarg + labsep * --ititle, 0.0f, 0.0f, awrk);
  };
/*
 * Write the Y-axis label(s).
 */
  sprintf(awrk, "%s%s%s", sp->sa->dophs ? "Phase" : "",
	                  sp->sa->dophs && sp->sa->doamp ? " and " : "",
	                  sp->sa->doamp ? "Amplitude" : "");
  cpgmtxt("L", lmarg, 0.5f, 0.5f, awrk);
/*
 * Write the X-axis label.
 */
  cpgmtxt("B", bmarg, 0.5f, 0.5f, sp->xlabel);
  return 0;
}

/*.......................................................................
 * Read the cursor position and return the cursor selection details in
 * sp->cursor.
 *
 * Input:
 *  sp     Vedpar *  The plot descriptor.
 *  mode Bandmode    If cpgband() is supported, mode specifies the following
 *                   cursor types:
 *                    B_NORM - A single point is required - no banding.
 *                    B_LINE - Line band between xref,yref and the cursor.
 *                    B_RECT - Rectangular band between xref,yref and the
 *                             cursor.
 *                    B_YRNG - Two horizontal lines bracketing a Y-axis range.
 *                    B_XRNG - Two vertical lines bracketing an X-axis range.
 *                    B_YVAL - Vertical line through the cursor.
 *                    B_XVAL - Horizontal line through the cursor.
 *                    B_CROSS- Cross hair centered on cursor.
 *  isamp     int    0 - yref denotes a phase.
 *                   1 - yref denotes an amplitude.
 *  iplot     int    The index of the sub-plot to which xref,yref refer.
 *  cif       int    The index of the IF to which xref,yref refer.
 *  xref    float    xref and yref denote the reference position for
 *  yref    float    band-type cursors. They are ignored if mode==B_NORM.
 *  ci        int    The color index to draw band cursors with.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int s_cursor(Specplot *sp, Bandmode mode, int isamp, int iplot,
		    int cif, float xref, float yref, int ci)
{
  static float xpos=0.5f; /* The X NDC position of the cursor */
  static float ypos=0.5f; /* The Y NDC position of the cursor */
  char key;               /* The user cursor key selection */
  SpCurs *sc;             /* Pointer to sp->cursor */
/*
 * Get the cursor descriptor.
 */
  sc = &sp->cursor;  
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
  if(sc->key == '\0') {
    xpos = 0.5f;
    ypos = 0.5f;
  };
/*
 * Initialize the return value.
 */
  sc->key = KEY_NONE;
  sc->waslow = 0;
  sc->wasamp = 0;
  sc->iplot = 0;
  sc->cif = 0;
  sc->x = 0.0f;
  sc->y = 0.0f;
/*
 * Substitute the cross-hair cursor for the normal cursor is requested.
 */
  if(mode == B_NORM && sp->sa->docross)
    mode = B_CROSS;
/*
 * Convert the cursor reference position into NDC.
 */
  switch(mode) {
  case B_LINE: 
  case B_RECT:
  case B_XRNG:
  case B_YRNG:
    {
      SpSubplot *sps = sp->splots + iplot;
      SpXdim *spx = sp->spx + cif;
      xref = spx->vxa + (spx->vxb - spx->vxa) *
	(xref - spx->xmin) / (spx->xmax - spx->xmin);
      if(isamp) {
	yref = sps->vymid + (sps->vyb - sps->vymid) *
	  (yref - sps->amin) / (sps->amax - sps->amin);
      } else {
	yref = sps->vya + (sps->vymid - sps->vya) *
	  (yref - sps->pmin) / (sps->pmax - sps->pmin);
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
  cpgsci(ci);
  if(!cpgband((int) mode, 0, xref, yref, &xpos, &ypos, &key))
    return 1;
/*
 * Convert key to upper case.
 */
  sc->waslow = islower((int)key);
  sc->key = (sc->waslow) ? toupper(key) : key;
/*
 * Locate the sub-plot selected by the cursor. This is either the
 * plot in which the cursor falls, or the sub-plot with the nearest X-axis
 * if the cursor is not within any sub-plot.
 */
  {
    int first=1;         /* True until mindiff has been initialized */
    float mindiff=0.0f;  /* Minimum distance between an edge and the cursor */
    for(iplot=0; iplot<sp->nplot; iplot++) {
      SpSubplot *sps = sp->splots + iplot;
/*
 * Found inside a sub-plot?
 */
      if(ypos >= sps->vya && ypos <= sps->vyb) {
	sc->iplot = iplot;
	break;
      } else {
/*
 * See if the current sub-plot has a closer X-axis to the cursor than
 * previous sub-plots.
 */
	float adiff = fabs(ypos - sps->vya);
	float bdiff = fabs(ypos - sps->vyb);
	if(first || adiff < mindiff) {
	  first=0;
	  sc->iplot = iplot;
	  mindiff = adiff;
	};
	if(bdiff < mindiff) {
	  sc->iplot = iplot;
	  mindiff = bdiff;
	};
      };
    };
  };
/*
 * Locate the IF sub-plot selected by the cursor. This is either the
 * plot in which the cursor falls, or the IF plot with the nearest Y-axis
 * if the cursor is not within any IF sub-plot.
 */
  {
    int first=1;         /* True until mindiff has been initialized */
    float mindiff=0.0f;  /* Minimum distance between an edge and the cursor */
    for(cif=0; cif<sp->nif; cif++) {
      SpXdim *spx = sp->spx + cif;
      if(spx->doplot) {
/*
 * Found inside a plot?
 */
	if(xpos >= spx->vxa && xpos <= spx->vxb) {
	  sc->cif = cif;
	  break;
	} else {
/*
 * See if the current IF sub-plot has a closer Y-axis to the cursor than
 * previous sub-plots.
 */
	  float adiff = fabs(xpos - spx->vxa);
	  float bdiff = fabs(xpos - spx->vxb);
	  if(first || adiff < mindiff) {
	    first=0;
	    sc->cif = cif;
	    mindiff = adiff;
	  };
	  if(bdiff < mindiff) {
	    sc->cif = cif;
	    mindiff = bdiff;
	  };
	};
      };
    };
  };
/*
 * Convert xpos from NDC to world coordinates.
 */
  {
    SpXdim *spx = sp->spx + sc->cif;
    sc->x = spx->xmin + (xpos - spx->vxa) * (spx->xmax - spx->xmin) /
            (spx->vxb - spx->vxa);
/*
 * Apply window limits. Do this in channels to prevent channels of other
 * IF plots being selected when the cursor is in the plot margin.
 */
    if(X_TO_CHAN(spx,sc->x) < spx->ca)
      sc->x = CHAN_TO_X(spx,spx->ca);
    if(X_TO_CHAN(spx,sc->x) > spx->cb)
      sc->x = CHAN_TO_X(spx,spx->cb);
  };
/*
 * Convert ypos from NDC to world coordinates.
 */
  {
    SpSubplot *sps = sp->splots + sc->iplot;
    sc->wasamp = (!sp->sa->dophs || (sp->sa->doamp && ypos >= sps->vymid));
    if(sc->wasamp) {
      sc->y = sps->amin + (ypos - sps->vymid) * (sps->amax - sps->amin) /
	      (sps->vyb - sps->vymid);
      if(sc->y < sps->amin)
	sc->y = sps->amin;
      if(sc->y > sps->amax)
	sc->y = sps->amax;
    } else {
      sc->y = sps->pmin + (ypos - sps->vya) * (sps->pmax - sps->pmin) /
	      (sps->vymid - sps->vya);
      if(sc->y < sps->pmin)
	sc->y = sps->pmin;
      if(sc->y > sps->pmax)
	sc->y = sps->pmax;
    };
  };
  return 0;
}

/*.......................................................................
 * Display plots interactively, according to cursor input from the user.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int s_interact(Specplot *sp)
{
  int finished=0; /* True when the user quits */
/*
 * Plot the first page.
 */
  if(s_page(sp, S_ALLNEW, 1) <= 0)
    return 1;
/*
 * Inform user of the way to receive usage information.
 */
  lprintf(stdout,
	  "For help move the cursor into the plot window and press \'%c\'.\n",
	  KEY_HELP);
/*
 * Enter the interactive cursor entry loop.
 */
  while(!finished) {
    int wasflag=0;  /* True if last key stroke toggled a flag */
    int nflag=0;    /* The number of flags enterred */
/*
 * Read the cursor.
 */
    do {
      if(s_cursor(sp, B_NORM, 0, 0, 0, 0.0f, 0.0f, 1))
	return 1;
      wasflag = s_flags(sp->sa, sp->cursor.key, sp->cursor.waslow)==0;
      if(wasflag)
	nflag++;
    } while(wasflag);
/*
 * Take action appropriate to the key that the user pressed.
 */
    if(nflag > 0) {     /* Update display after a sequence of flag toggling */
      if(s_redisp(sp))
	return 1;
    } else {
      switch (sp->cursor.key) {
      case KEY_QUIT:
	finished = 1;
	break;
      case KEY_CROSS:
	sp->sa->docross = !sp->sa->docross;
	break;
      case KEY_DISP:
	if(s_redisp(sp))
	  return 1;
	break;
      case KEY_XAXIS:
	if(sp->cursor.waslow ? s_get_xrange(sp) : s_get_xaxis(sp))
	  return 1;
	break;
      case KEY_VECT:
	sp->sa->avmode = sp->sa->avmode==SP_VECTOR ? SP_SCALAR : SP_VECTOR;
	if(s_page(sp, S_RESET, 1) < 0)
	  return 1;
	break;
      case KEY_ZOOM:
	if(s_get_yrange(sp))
	  return 1;
	break;
      case KEY_NEXT:
	if(s_page(sp, S_NEXT, 1) < 0)
	  return 1;
	break;
      case KEY_PREV:
	if(s_page(sp, S_NEXT, 0) < 0)
	  return 1;
	break;
      case KEY_NUMB:
	if(sp->cursor.waslow ? s_newnum(sp) : s_get_smooth(sp))
	  return 1;
	break;
      case KEY_SPEC:
	if(s_get_sel(sp))
	  return 1;
	break;
      case KEY_ORDER:
	if(s_get_order(sp))
	  return 1;
	break;
      case KEY_HELP:
	printf("Specplot key bindings:\n");
	printf(" %c - List the following key bindings.\n", KEY_HELP);
	printf(" %c - Exit specplot (right-mouse-button).\n", KEY_QUIT);
	printf(" %c - Redisplay the current plot.\n", KEY_DISP);
	printf(" %c - Display the next page of spectra.\n", KEY_NEXT);
	printf(" %c - Display the preceding page of spectra.\n", KEY_PREV);
	printf(" %c - Change the number of plots per page.\n",
	       tolower(KEY_NUMB));
	printf(" %c - Change the smoothing parameters.\n", KEY_NUMB);
	printf(" %c - Plot only amplitudes.\n", KEY_AMP);
	printf(" %c - Plot only phases.\n", KEY_PHS);
	printf(" %c - Plot both amplitudes and phases.\n", KEY_BOTH);
	printf(" %c - Toggle error bars on/off.\n", KEY_ERR);
	printf(" %c - Toggle between vector and scalar averaging.\n", KEY_VECT);
	printf(" %c - Change the baselines, polarization or times used.\n",
	       KEY_SPEC);
	printf(" %c - Change the sort-order of selections.\n", KEY_ORDER);
	printf(" %c - Delimit a new channel range with the cursor. (hit %c twice for full range).\n",
	       tolower(KEY_XAXIS), tolower(KEY_XAXIS));
	printf(" %c - Change the X-axis type, smoothing function and FWHM.\n",
	       KEY_XAXIS);
	printf(" %c - Select a new amplitude or phase range (hit %c twice for full range).\n", KEY_ZOOM, KEY_ZOOM);
	printf(" %c - Toggle whether to join adjacent points with lines.\n",
	       tolower(KEY_JOIN));
	printf(" %c - Toggle whether to draw lines as bins or vectors.\n",
	       KEY_JOIN);
	printf(" %c - Toggle whether to use a crosshair cursor if available.\n",
	       KEY_CROSS);
	printf("\n");
	break;
      default:
	break;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Handle non-interactive plotting.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  npage      int     The maximum number of pages to plot when plotting
 *                     to a non-interactive device. A value of 0 is
 *                     means no limit.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int s_auto(Specplot *sp, int npage)
{
  int ierr = 0;  /* Error return status */
/*
 * Plot the first page.
 */
  if(s_page(sp, S_ALLNEW, 1) <= 0)
    return 1;
/*
 * Plot the rest of the pages.
 */
  while(!ierr && (npage <= 0 || sp->npage < npage)) {
    int nplotted = s_page(sp, S_NEXT, 1);
    if(nplotted < 0)
      ierr = 1;
    else if(nplotted == 0)
      break;
  };
  return (ierr < 0);
}

/*.......................................................................
 * Create a new specplot attributes container and initialize it to
 * defaults.
 *
 * Input:
 *  ob    Observation *  The observation to be plotted.
 * Output:
 *  return   Specattr *  The new descriptor, or NULL on error.
 */
Specattr *new_Specattr(Observation *ob)
{
  Specattr *sa;         /* The container to be returned */
  int i;
/*
 * List X-axis type SpXunit enumerations versus their names.
 */
  static Enumpar xttab[]   = {{"channels", SP_CHAN}, {"frequency", SP_FREQ}};
/*
 * List smoothing function SmType enumerations versus their names.
 */
  static Enumpar smtab[]   = {{"none", SM_NONE}, {"hanning", SM_HANNING},
			      {"gaussian", SM_GAUSSIAN}, {"boxcar", SM_BOXCAR},
			      {"sinc", SM_SINC}};
/*
 * List major mode SpMode enumerations versus their names.
 */
  static Enumpar keytab[] = {{"baseline", SP_BASE}, {"polarization", SP_POL},
			     {"time", SP_TIME},  {"uvrange", SP_UVR}};
/*
 * List available averaging modes.
 */
  static Enumpar avtab[] = {{"vector", SP_VECTOR}, {"scalar", SP_SCALAR}};
/*
 * List available baseline selection modes.
 */
  static Enumpar bmtab[] = {{"split", SP_SPLIT}, {"group", SP_GROUP}};
/*
 * Allocate the container.
 */
  sa = (Specattr *) malloc(sizeof(Specattr));
  if(!sa) {
    lprintf(stderr, "new_Specattr: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize all members of the container, at least to the point at which
 * it can safely be sent to del_Specattr().
 */
  sa->stime = ob->rec[0].integ->ut;
  sa->etime = ob->rec[ob->nrec-1].integ->ut;
  sa->scan = sa->etime - sa->stime;
  sa->amin = sa->amax = 0.0f;
  sa->pmin = -180.0f;
  sa->pmax = 180.0f;
  sa->ca = 0;
  sa->cb = ob->nchan * ob->nif - 1;
  sa->nplot = 3;
  sa->avmode = SP_VECTOR;
  sa->doamp = 1;
  sa->dophs = 1;
  sa->docross = 0;
  sa->dojoin = 1;
  sa->dohist = 1;
  sa->dobars = 1;
  sa->pl = NULL;
  sa->bgl = NULL;
  sa->uvr.uvrlim = -1.0f;  /* -1 signifies an uninitialized state */
  sa->uvr.uvmin = sa->uvr.uvmax = sa->uvr.uvstep = 0.0f;
  for(i=0; i<SP_NKEY; i++)
    sa->key[i] = (SpKey) i;
  sa->nkey = 1;
  sa->xunit = SP_CHAN;
  sa->smooth.type = SM_NONE;
  sa->smooth.fwhm = 0.0;
  sa->bmode = SP_SPLIT;
  sa->xtsym = NULL;
  sa->keysym = NULL;
  sa->smsym = NULL;
/*
 * Set the default list of baseline selection groups.
 */
  if(sp_set_bgl(ob, sa, SP_SPLIT, NULL))
    return del_Specattr(sa);
/*
 * Set the default list of polarizations.
 */
  if(sp_set_pol(ob, sa, NULL))
    return del_Specattr(sa);
/*
 * Set the default UV range.
 */
  if(sp_set_uvrange(ob, sa, 0.0f, 0.0f, 0.0f))
    return del_Specattr(sa);
/*
 * Initialize type-name symbol tables.
 */
  sa->xtsym = new_Enumtab(xttab, sizeof(xttab)/sizeof(Enumpar),
			  "Specplot x-axis type");
  if(!sa->xtsym)
    return del_Specattr(sa);
  sa->keysym = new_Enumtab(keytab, sizeof(keytab)/sizeof(Enumpar),
			   "Specplot selection");
  if(!sa->keysym)
    return del_Specattr(sa);
  sa->smsym = new_Enumtab(smtab, sizeof(smtab)/sizeof(Enumpar),
			  "Specplot smoothing function");
  if(!sa->smsym)
    return del_Specattr(sa);
  sa->avsym = new_Enumtab(avtab, sizeof(avtab)/sizeof(Enumpar),
			  "Specplot averaging mode");
  if(!sa->avsym)
    return del_Specattr(sa);
  sa->bmsym = new_Enumtab(bmtab, sizeof(bmtab)/sizeof(Enumpar),
			  "Specplot baseline selection mode");
  if(!sa->bmsym)
    return del_Specattr(sa);
/*
 * Return the initialized container.
 */
  return sa;
}

/*.......................................................................
 * Delete a specplot attributes container and its contents.
 *
 * Input:
 *  sa     Specattr *  The container to be deleted.
 * Output:
 *  return Specattr *  Allways NULL.
 */
Specattr *del_Specattr(Specattr *sa)
{
  if(sa) {
    sa->pl = del_Pollist(sa->pl);
    sa->bgl = del_Bgrplist(sa->bgl);
    sa->xtsym = del_Enumtab(sa->xtsym);
    sa->keysym = del_Enumtab(sa->keysym);
    sa->smsym = del_Enumtab(sa->smsym);
    sa->avsym = del_Enumtab(sa->avsym);
    sa->smsym = del_Enumtab(sa->smsym);
  };
  return NULL;
}

/*.......................................................................
 * Replace the current list of polarizations.
 *
 * Input:
 *  sa     Specattr *   The specplot attributes container.
 *  pl      Pollist *   The new list of polarizations, or NULL for the default.
 *                      Note that ownership of the list is transfered
 *                      at this point. It should not hereafter be deleted or
 *                      changed by the caller.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
int sp_set_pol(Observation *ob, Specattr *sa, Pollist *pl)
{
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "sp_set_pol"))
    return 1;
  if(bad_Specattr(sa, "sp_set_pol"))
    return 1;
/*
 * If the list contains no polarizations the container is superfluous.
 */
  if(pl && pl->npol < 1)
    pl = del_Pollist(pl);
/*
 * Discard the previous list, if any.
 */
  sa->pl = del_Pollist(sa->pl);
/*
 * Install the new list.
 */
  sa->pl = pl;
  return 0;
}

/*.......................................................................
 * Replace the current list of baseline selection groups.
 *
 * Input:
 *  sa     Specattr *   The specplot attributes container.
 *  bmode   SpBMode     A baseline selection mode, from:
 *                       SP_SPLIT   -  Show each baseline of the first group.
 *                       SP_GROUP   -  Show each baseline group.
 *  bgl    Bgrplist *   The new list of baseline groups, or NULL
 *                      to select all baselines.
 *                      Note that ownership of the list is transfered
 *                      at this point. It should not hereafter be deleted or
 *                      changed by the caller.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
int sp_set_bgl(Observation *ob, Specattr *sa, SpBMode bmode, Bgrplist *bgl)
{
  Basegrp *bgrp;  /* A baseline group */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "sp_set_bgl"))
    return 1;
  if(bad_Specattr(sa, "sp_set_bgl"))
    return 1;
/*
 * Ensure that there is at least one baseline group.
 */
  if(!bgl && (bgl = new_Bgrplist()) == NULL)
    return 1;
  if(bgl->nsel<1 && add_Basegrp(ob, bgl, NULL, "")==NULL) {
    bgl = del_Bgrplist(bgl);
    return 1;
  };
/*
 * Check each selection to ensure that at least one baseline is selected.
 */
  for(bgrp=bgl->bgrp; bgrp; bgrp = bgrp->next) {
    if(size_Basegrp(ob, bgrp, -1) < 1) {
      lprintf(stderr, "sp_set_bgl: Empty baseline selection.\n");
      bgl = del_Bgrplist(bgl);
      return 1;
    };
  };
/*
 * Discard the previous list, if any.
 */
  sa->bgl = del_Bgrplist(sa->bgl);
/*
 * Install the new list and selection mode.
 */
  sa->bmode = bmode;
  sa->bgl = bgl;
  return 0;
}

/*.......................................................................
 * Record the start time, end time and scan delimiter to be used
 * in specplot.
 *
 * Input:
 *  ob  Observation *   The observation to be plotted with specplot.
 *  sa     Specattr *   The specplot attributes container.
 *  stime    double     The first time to sample (seconds into year).
 *  etime    double     The last time to sample (seconds into year).
 *  scan     double     If scan >= 0 then scans are groups of
 *                      integrations that cover up to 'scan' minutes
 *                      each.
 *                      If scan < 0 then scans are groups of
 *                      integrations that are separated from their
 *                      neigbouring groups by at least -scan seconds.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
int sp_set_times(Observation *ob, Specattr *sa, double stime, double
		 etime, double scan) 
{
/*
 * Check the descriptor.
 */
  if(bad_Specattr(sa, "sp_set_times"))
    return 1;
/*
 * Ensure that stime <= etime.
 */
  if(stime > etime) {
    double tmp_time = stime;
    stime = etime;
    etime = tmp_time;
  };
/*
 * Search for the closest times that match the specification.
 */
  {
    int ut_a = ob_find_ut(ob, stime, UT_GE);
    int ut_b = ob_find_ut(ob, etime, UT_LE);
    if(ut_a < 0 || ut_b < 0 || ut_a > ut_b) {
      lprintf(stderr, "sp_set_times: Time range unsampled.\n");
      return 1;
    };
    stime = ob->rec[ut_a].integ->ut;
    etime = ob->rec[ut_b].integ->ut;
  };
/*
 * Record the adjusted time parameters.
 */
  sa->stime = stime;
  sa->etime = etime;
  sa->scan = scan;
  return 0;
}

/*.......................................................................
 * Record the overall range and iterator increment of UV radii to be
 * sampled by specplot spectra. All radius arguments are measured in
 * wavelengths.
 *
 * Input:
 *  ob  Observation *   The observation to be plotted with specplot.
 *  sa     Specattr *   The specplot attributes container.
 *  uvmin     float     The smallest UV radius to sample.
 *  uvmax     float     The largest UV radius to sample, or 0.0 to
 *                      request the max UV radius available.
 *  uvstep    float     The annular width to break the range uvmin -> uvmax
 *                      into.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
int sp_set_uvrange(Observation *ob, Specattr *sa, float uvmin, float uvmax,
		   float uvstep)
{
  SpUV *uvr;
/*
 * Check the descriptor.
 */
  if(bad_Specattr(sa, "sp_set_uvrange"))
    return 1;
/*
 * Get the UV range container.
 */
  uvr = &sa->uvr;
/*
 * If the overall UV radius range hasn't been determined yet, do so now?
 */
  if(uvr->uvrlim < 0.0f) {
    double maxfreq;         /* Max frequency observed */
    double uvrmax=0.0;      /* Max UV radius (light seconds) */
    Subarray *sub;
    for(sub=ob->sub; sub<ob->sub + ob->nsub; sub++) {
      Integration *integ;
      for(integ=sub->integ; integ<sub->integ + sub->ntime; integ++) {
	Visibility *vis;
	for(vis=integ->vis; vis<integ->vis + sub->nbase; vis++) {
	  float uvrad = sqrt(vis->u * vis->u + vis->v * vis->v);
	  if(uvrad > uvrmax)
	    uvrmax = uvrad;
	};
      };
    };
/*
 * Find the maximum frequency sampled so that we can convert the
 * max uv radius from light seconds to the max numbers of wavelengths.
 */
    {
      If *ifp;
      maxfreq = 0.0;
      for(ifp=ob->ifs; ifp<ob->ifs + ob->nif; ifp++) {
	float fmax = ifp->freq + (ifp->df > 0.0 ? ifp->bw : 0.0);
	if(fmax > maxfreq)
	  maxfreq = fmax;
      };
    };
/*
 * Record the new max UV radius in wavelengths.
 */
    uvr->uvrlim = maxfreq * uvrmax;
  };
/*
 * If all UV values are zero, set a UV range of zero.
 */
  if(uvr->uvrlim == 0.0f) {
    uvmin = uvmax = uvstep = 0.0f;
  } else {
/*
 * Limit the entered UV radii.
 */
    if(uvmin < 0.0f)
      uvmin = 0.0f;
/*
 * If uvmax is 0 or greater than the available range, substitute the
 * max available radius.
 */
    if(uvmax <= 0.0f || uvmax > uvr->uvrlim)
      uvmax = uvr->uvrlim;
/*
 * Interpret the step size.
 */
    if(uvstep <= 0.0f || uvstep > uvmax-uvmin)
      uvstep = uvmax - uvmin;
/*
 * Check the new range.
 */
    if(uvstep <= 0.0f) {
      lprintf(stderr,
	      "sp_set_uvrange: No data lie within the requested UV range.\n");
      return 1;
    };
  };
/*
 * Record the new UV range parameters.
 */
  uvr->uvmin = uvmin;
  uvr->uvmax = uvmax;
  uvr->uvstep = uvstep;
  return 0;
}

/*.......................................................................
 * Record specplot smoothing parameters.
 *
 * Input:
 *  sa     Specattr *   The specplot attributes container.
 *  xunit   SpXunit     The units of fwhm.
 *                       SP_CHAN   - Channel indexes.
 *                       SP_FREQ   - Frequency.
 *  smtype     char *   The type of smoothing window to be applied.
 *                      Let v=channel, hwhm= HWHM of smoothing fn, and x=v/hwhm.
 *                       SM_NONE     - 1.0 (No smoothing).
 *                       SM_HANNING  - 0.5*(1 + cos(2*pi*x)) for x=-0.5..0.5,
 *                                     and 0 beyond.
 *                       SM_GAUSSIAN - exp(-ln(2)*x^2).
 *                       SM_BOXCAR   - 1.0 for x=-1..1, and 0 elsewhere.
 *                       SM_SINC     - sin(c*x)/(c*x), with c=1.8954942670340.
 *  fwhm      float     The full-width at half maximum to give the smoothing
 *                      function.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
int sp_set_smooth(Specattr *sa, SpXunit xunit, SmType smtype, float fwhm)
{
/*
 * Check arguments.
 */
  if(bad_Specattr(sa, "sp_set_smooth"))
    return 1;
  sa->smooth.xunit = xunit;
  sa->smooth.type = smtype;
  sa->smooth.fwhm = fwhm >= 0.0f ? fwhm : 0.0f;
  return 0;
}

/*.......................................................................
 * Input:
 *  sa     Specattr *   The specplot attributes container.
 *  nplot        int    The number of plots per page (0 -> 3).
 *  xunit   SpXunit     The units along the X axis.
 *                       SP_CHAN   - Channel indexes.
 *                       SP_FREQ   - Frequency.
 *  avmode  SpAvMode    Visibility time-averaging mode.
 *                       SP_VECTOR - Vector average.
 *                       SP_SCALAR - Scalar average.
 */
int sp_set_options(Specattr *sa, int nplot, SpXunit xunit, SpAvMode avmode)
{
/*
 * Check arguments.
 */
  if(bad_Specattr(sa, "sp_set_options"))
    return 1;
/*
 * Record the modes.
 */
  sa->nplot = nplot <= 0 ? 3 : nplot;
  sa->xunit = xunit;
  sa->avmode = avmode;
  return 0;
}

/*.......................................................................
 * Set the order of the specplot selection keys.
 *
 * Input:
 *  sa     Specattr *   The specplot attributes container.
 *  keys      SpKey     An array of unique selection keys in the
 *                      required sort order, from:
 *                       SP_BASE   - Display each baseline selection.
 *                       SP_POL    - Plot each specified polarization.
 *                       SP_TIME   - Plot each time scan.
 *  nkey        int     The number of keys in keys[].
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
int sp_set_order(Specattr *sa, SpKey *keys, int nkey)
{
  SpKey s_keys[SP_NKEY];  /* An array in which to construct the list of keys */
  int i;
/*
 * Check arguments.
 */
  if(bad_Specattr(sa, "sp_set_order"))
    return 1;
  if(keys==NULL || nkey < 1) {
    lprintf(stderr, "sp_set_order: No selection keys received.\n");
    return 1;
  };
/*
 * Make a temporary copy the current array of keys.
 */
  for(i=0; i<SP_NKEY; i++)
    s_keys[i] = sa->key[i];
/*
 * For each new key displace the previous key that occupied the specified
 * ordinal position to the position previously occupied by the new key.
 */
  for(i=0; i<nkey; i++) {
    SpKey oldkey = s_keys[i];
    SpKey newkey = keys[i];
    int j;
/*
 * Search for the position previously occupied by the new key.
 */
    for(j=i; j<SP_NKEY && newkey != s_keys[j]; j++)
      ;
    if(j >= SP_NKEY) {
      lprintf(stderr, "sp_set_order: \"%s\" is cited more than once.\n",
	      name_enum(sa->keysym, newkey, "(unknown)"));
      return 1;
    };
/*
 * Displace the old key with the new one, and preserve the old one
 * at the position previously held by the new one.
 */
    s_keys[j] = oldkey;
    s_keys[i] = newkey;
  };
/*
 * Install the new keys.
 */
  sa->nkey = nkey;
  for(i=0; i<SP_NKEY; i++)
    sa->key[i] = s_keys[i];
  return 0;
}

/*.......................................................................
 * Record specplot axis ranges.
 *
 * Input:
 *  ob  Observation *   The observation who's spectra will be plotted.
 *  sa     Specattr *   The specplot attributes container.
 *  ca,cb       int     The range of channel indexes to plot.
 *                      Channel 0, is the first channel of the first IF.
 *                      Channels increase across IF boundaries.
 *                      -1,-1 selects the full range.
 *  amin      float     The minimum amplitude to plot.
 *  amax      float     The maximum amplitude to plot.
 *                      For autoscaling send amin == amax = 0.
 *  pmin      float     Minimum phase to plot (degrees).
 *  pmax      float     Maximum phase to plot (degrees).
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
int sp_set_axes(Observation *ob, Specattr *sa, int ca, int cb,
		float amin, float amax, float pmin, float pmax)
{
  int max_chan;
/*
 * Check arguments.
 */
  if(bad_Specattr(sa, "sp_set_axes"))
    return 1;
/*
 * Get the maximum channel index.
 */
  max_chan = ob->nchan * ob->nif - 1;
/*
 * Limit the channel index range to 0 -> max_chan.
 */
  if(ca < 0)
    ca = 0;
  if(cb < 0)
    cb = max_chan;
  if(ca > max_chan)
    ca = max_chan;
  if(cb > max_chan)
    cb = max_chan;
/*
 * Keep the indexes in increasing order.
 */
  if(ca > cb) {
    int ctmp = ca;
    ca = cb;
    cb = ctmp;
  };
/*
 * Record the new channel range.
 */
  sa->ca = ca;
  sa->cb = cb;
/*
 * Install the new amplitude range.
 */
  sa->amin = amin;
  sa->amax = amax;
/*
 * Limit the phase range to -180 -> +180 degrees.
 */
  if(pmin < -180.0f) pmin = -180.0f;
  if(pmax < -180.0f) pmax = -180.0f;
  if(pmin > 180.0f) pmin = 180.0f;
  if(pmax > 180.0f) pmax = 180.0f;
/*
 * Install the new phase range.
 */
  if(pmin > pmax) {
    float ptmp = pmin;
    pmin = pmax;
    pmax = ptmp;
  };
  if(pmin == pmax) {
    pmin = -180.0f;
    pmax = 180.0f;
  };
  sa->pmin = pmin;
  sa->pmax = pmax;
  return 0;
}

/*.......................................................................
 * Report an error if the given Specattr descriptor is invalid.
 * Currently this only checks for sa being non-NULL.
 *
 * Input:
 *  sa    Specattr *  The descriptor to be checked.
 *  fname     char *  The name of the calling function.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int bad_Specattr(Specattr *sa, char *fname)
{
  if(!sa) {
    lprintf(stderr, "%s: NULL Specattr descriptor.\n",
	    fname ? fname : "bad_Specattr");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Allow interactive selection of the channel range to be plotted.
 *
 * Input:
 *  sp   Specplot *   The spectrum plot descriptor.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
static int s_get_xrange(Specplot *sp)
{
  int dofull=0;   /* True if full channel range requested */
  int iter;       /* Iterate over getting two valid keypresses */
  int chval[2];   /* The selected channel indexes */
  SpCurs ref;     /* Reference cursor selection */
/*
 * Get the two channel cursor selections.
 */
  ref = sp->cursor;
  for(iter=0; iter<2 && !dofull; iter++) {
    int accepted = 0;
    while(!accepted) {
      if(s_cursor(sp, iter==0 ? B_XVAL:B_XRNG, ref.wasamp, ref.iplot, ref.cif,
		  ref.x, ref.y, zoomcol))
	return 1;
      switch(sp->cursor.key) {
      case KEY_XAXIS:        /* Display full range and end range selection */
	accepted = dofull = 1;
	break;
      case KEY_QUIT:
      case KEY_CAN:         /* Abort selection */
	return 0;
	break;
      case KEY_CUR:         /* Accept the selected limit */
	ref = sp->cursor;
	{
	  SpXdim *spx = sp->spx + ref.cif;
	  chval[iter] = X_TO_CHAN(spx, ref.x) + spx->cmin;
	};
	accepted=1;
	break;
      default:              /* Unexpected cursor input key - show usage */
	printf("%c - Select the position of the %s X-axis limit.\n", KEY_CUR,
	       iter==0 ? "start":"end");
	printf("%c - Abort selection.\n", KEY_CAN);
	printf("%c - Revert to the full range.\n", KEY_XAXIS);
	break;
      };
    };
  };
/*
 * Set up the new channel index range.
 */
  {
    Specattr *sa = sp->sa;
    int max_chan = sp->ob->nchan * sp->ob->nif - 1;
    int ca, cb;
    if(dofull) {
      ca = 0;
      cb = max_chan;
    } else if(chval[0] <= chval[1]) {
      ca = chval[0];
      cb = chval[1];
    } else {
      ca = chval[1];
      cb = chval[0];
    };
/*
 * Install the new channel index range.
 */
    if(sp_set_axes(sp->ob, sa, ca, cb, sa->amin, sa->amax, sa->pmin, sa->pmax))
      return 0;
  };
/*
 * Display the result.
 */
  return s_redisp(sp);
}

/*.......................................................................
 * Allow interactive selection of the plotted amplitude or phase range.
 *
 * Input:
 *  sp   Specplot *   The spectrum plot descriptor.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
static int s_get_yrange(Specplot *sp)
{
  int dofull=0;   /* True to unzoom the amplitude and phase ranges */
  int iter;       /* Iterate over getting two valid keypresses */
  float yval[2];  /* The selected Y values */
  SpCurs ref;     /* Reference cursor selection */
/*
 * Get the two channel cursor selections.
 */
  ref = sp->cursor;
  for(iter=0; iter<2 && !dofull; iter++) {
    int accepted = 0;
    while(!accepted) {
      if(s_cursor(sp, iter==0 ? B_YVAL:B_YRNG, ref.wasamp, ref.iplot, ref.cif,
		  ref.x, ref.y, zoomcol))
	return 1;
      switch(sp->cursor.key) {
      case KEY_ZOOM:        /* Display full range and end selection */
	accepted = dofull = 1;
	break;
      case KEY_QUIT:
      case KEY_CAN:         /* Abort channel selection */
	return 0;
	break;
      case KEY_CUR:         /* Accept the selected y value */
	if(iter==1 && (!sp->cursor.wasamp != !ref.wasamp ||
		       sp->cursor.iplot != ref.iplot)) {
	  printf("Second selection in a different sub-plot. Selection aborted.\n");
	  return 0;
	} else {
	  yval[iter] = sp->cursor.y;
	  ref = sp->cursor;
	};
	accepted=1;
	break;
      default:              /* Unexpected cursor input key - show usage */
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
 * Set up the new channel amplitude and/or phase ranges.
 */
  {
    Specattr *sa = sp->sa;
/*
 * Start with the current ranges.
 */
    float amin = sa->amin;
    float amax = sa->amax;
    float pmin = sa->pmin;
    float pmax = sa->pmax;
/*
 * Override the changed ranges.
 */
    if(dofull) {
      amin = amax = 0.0f;
      pmin = pmax = 0.0f;
    } else if(ref.wasamp) {
      amin = yval[0];
      amax = yval[1];
    } else {
      pmin = yval[0];
      pmax = yval[1];
    };
/*
 * Install the new ranges.
 */
    sp_set_axes(sp->ob, sa, sa->ca,sa->cb, amin,amax, pmin,pmax);
  };
/*
 * Display the result.
 */
  return s_redisp(sp);
}

/*.......................................................................
 * Allow the user to select a new X-axis type and smoothing function
 * from their terminal.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Fatal error.
 */
static int s_get_xaxis(Specplot *sp)
{
  char *cptr;     /* Pointer into the returned line */
  char *arg;      /* Pointer to a copy of the next argument */
  Enumpar *epar;  /* A symbol table entry */
  Specattr *sa = sp->sa;
  SpXunit xunit=SP_CHAN;
/*
 * Get the new number from the user.
 */
  printf("Enter the required X-axis type.\n");
  if((cptr=s_getline("[Default = channels]: ")) == NULL)
    return 0;
/*
 * Get the X-axis type.
 */
  if(*cptr) {
    if(!(arg=s_get_arg(cptr, &cptr)))
      return 0;
    if((epar=find_enum(sa->xtsym, arg)) == NULL)
      return 0;
    xunit = (SpXunit) epar->id;
  };
  if(*cptr) {
    fprintf(stderr, "Unexpected input at end of string: %s\n", cptr);
    return 0;
  };
  sp_set_options(sa, sa->nplot, xunit, sa->avmode);
  return s_redisp(sp);
}

/*.......................................................................
 * Allow the user to select new smoothing parameters from their terminal.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Fatal error.
 */
static int s_get_smooth(Specplot *sp)
{
  char *cptr;     /* Pointer into the returned line */
  char *arg;      /* Pointer to a copy of the next argument */
  Enumpar *epar;  /* A symbol table entry */
  SpXunit xunit=SP_CHAN;
  SmType smtype=SM_NONE;
  float fwhm = 0.0f;
/*
 * Get the new number from the user.
 */
  printf("Enter zero or more of the fwhm units, smoothing function and fwhm.\n");
  if((cptr=s_getline("[Default = channels, none, 0.0]: ")) == NULL)
    return 0;
/*
 * Get the X-axis type.
 */
  if(*cptr) {
    if(!(arg=s_get_arg(cptr, &cptr)))
      return 0;
    if((epar=find_enum(sp->sa->xtsym, arg)) == NULL)
      return 0;
    xunit = (SpXunit) epar->id;
  };
/*
 * Get the smoothing function type.
 */
  if(*cptr) {
    if(*cptr == ',')
      cptr++;
    if(!(arg=s_get_arg(cptr, &cptr)))
      return 0;
    if((epar=find_enum(sp->sa->smsym, arg)) == NULL)
      return 0;
    smtype = (SmType) epar->id;
  };
/*
 * Get the smoothing function fwhm.
 */
  if(*cptr) {
    if(*cptr == ',')
      cptr++;
    if(!(arg=s_get_arg(cptr, &cptr)))
      return 0;
    {
      char *eptr = arg;
      fwhm = strtod(arg, &eptr);
      if(arg==eptr || fwhm < 0.0) {
	fprintf(stderr, "Unacceptable fwhm: %s\n", arg);
	return 0;
      };
    };
  };
  if(*cptr) {
    fprintf(stderr, "Unexpected input at end of string: %s\n", cptr);
    return 0;
  };
  return sp_set_smooth(sp->sa, xunit, smtype, fwhm) || s_redisp(sp);
}

/*.......................................................................
 * Allow the user to select a new selection key order.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Fatal error.
 */
static int s_get_order(Specplot *sp)
{
  char *cptr;           /* Pointer into the returned line */
  char *arg;            /* Pointer to a copy of the next argument */
  Enumpar *epar;        /* A symbol table entry */
  SpKey keys[SP_NKEY];  /* The array of narg keys */
  int nkey;             /* The number of keys selected */
/*
 * Get the new number from the user.
 */
  printf("Enter one or more selection types in the desired sort-order.\n");
  if((cptr=s_getline("[eg. baseline]: ")) == NULL)
    return 0;
/*
 * Get the X-axis type.
 */
  nkey = 0;
  do {
    if(!(arg=s_get_arg(cptr, &cptr)))
      return 0;
    if((epar=find_enum(sp->sa->keysym, arg)) == NULL)
      return 0;
    keys[nkey++] = (SpKey) epar->id;
    if(*cptr==',')
      cptr++;
  } while(*cptr && nkey<SP_NKEY);
/*
 * Too many arguments?
 */
  if(*cptr) {
    fprintf(stderr, "Too many arguments.\n");
    return 0;
  };
  if(sp_set_order(sp->sa, keys, nkey))
    return 0;
  return s_page(sp, S_ALLNEW, 1) < 0;
}

/*.......................................................................
 * Utility routine to return a copy of the next comma separated argument
 * from a given string. Leading white-space is ignored. Then all input
 * up to the next comma, or to the end of the string if there is no
 * comma, will be copied into the returned string. Trailing
 * white-space will then be removed from the copy.
 *
 * Input:
 *  string    char *   The string to read from.
 *  endp      char **  If endp==NULL then a comma in the string is
 *                     treated as an error and reported accordingly.
 *                     If endp!=NULL then the pointer to the next
 *                     unprocessed character in string[] will be
 *                     assigned to *endp. This will either point at a
 *                     comma or the end of the string.
 * Output:
 *  return    char *   A pointer to an internal static character array
 *                     containing the copy of the new argument, or
 *                     NULL on error.
 */
static char *s_get_arg(char *string, char **endp)
{
  enum {S_MAX_ARG=80};
  static char outstr[S_MAX_ARG+1];  /* The output buffer */
  char *sptr;    /* Pointer to the first non-white-space character */
  char *eptr;    /* Pointer to the end of the argument */
  int nc;        /* The number of characters to be copied */
  int i;
/*
 * Initialize the returned record of the next unprocessed character.
 */
  if(endp) *endp = string;
  sptr = string;
/*
 * Skip leading white-space.
 */
  while(*sptr && isspace((int)*sptr))
    sptr++;
/*
 * Locate the end of the argument and the last white-space character.
 */
  nc = 0;
  for(eptr=sptr; *eptr && *eptr != ','; eptr++) {
    if(!isspace((int)*eptr))
      nc = eptr - sptr + 1;
  };
/*
 * Missing argument?
 */
  if(nc==0) {
    fprintf(stderr, "Missing argument.\n");
    return NULL;
  };
/*
 * If terminated by a comma and no more input is expected, report an
 * error and quit.
 */
  if(!endp && *eptr == ',') {
    lprintf(stderr, "Unexpected input at end of string: \"..%s\".\n",
	    eptr);
    return NULL;
  };
/*
 * Argument too long?
 */
  if(nc > S_MAX_ARG) {
    fprintf(stderr, "Argument too long.\n");
    return NULL;
  };
/*
 * Copy up to the last non-white-space character.
 */
  for(i=0; i<nc; i++)
    outstr[i] = sptr[i];
  outstr[i] = '\0';
  if(endp)
    *endp = eptr;
  return outstr;
}

/*.......................................................................
 * Prompt for and read a line of input from the user. The string is
 * terminated at the newline character.
 *
 * Input:
 *  prompt   char *  A mandatory prompt to present to the user.
 * Output:
 *  return   char *  A pointer to an internal static copy of the input
 *                   line, or NULL on error.
 */
static char *s_getline(char *prompt)
{
  enum {S_MAX_LINE=132};
  static char outstr[S_MAX_LINE+1]; /* The output buffer */
/*
 * Prompt the user and read their response.
 */
  printf("%s", prompt);
  if(fgets(outstr, (int) S_MAX_LINE, stdin)==NULL) {
    lprintf(stderr, "Error reading input.\n");
    return NULL;
  };
/*
 * Locate the newline character and replace terminate the string there.
 */
  {
    char *cptr = strchr(outstr, '\n');
    if(!cptr) {
      int c;
      lprintf(stderr, "Input line too long.\n");
      do c=getc(stdin); while(c != '\n' && c != EOF);
      return NULL;
    };
    *cptr = '\0';
  };
  return outstr;
}

/*.......................................................................
 * Allow a new selection list to be entered.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Fatal error.
 */
static int s_get_sel(Specplot *sp)
{
  Enumpar *epar;  /* A mode symbol table entry */
  SpKey key;      /* The selection type */
  char *cptr;     /* A pointer into a copy of the input line */
  char *arg;      /* An argument from the input line */
/*
 * Read the user's specifications.
 */
  printf("Enter a new selection specification.\n");
  do {
    cptr = s_getline("[? for a list of options]: ");
    if(!cptr)
      return 0;
/*
 * Get the first argument.
 */
    if(!(arg=s_get_arg(cptr, &cptr)))
      return 0;
    if(*arg == '?') {
      printf("Selection options include:\n");
      printf(" time, <start-time>, <end-time>, <scan_time>\n");
      printf(" polarization, <polarization-name>, ...\n");
      printf(" baseline, group, <baseline-group1>, ...\n");
      printf(" baseline, split, <baseline-group>\n");
      printf(" uvrange, <uvmin>, <uvmax>, <uvstep>\n");
    };
  } while(*arg == '?');
/*
 * Decode the mode argument.
 */
  if((epar=find_enum(sp->sa->keysym, arg)) == NULL)
    return 0;
  key = (SpKey) epar->id;
/*
 * Decode mode-specific arguments.
 */
  switch(key) {
  case SP_BASE:
    if(*cptr && s_get_bgl(sp, cptr+1))
      return 0;
    break;
  case SP_POL:
    if(*cptr && s_get_pol(sp, cptr+1))
      return 0;
    break;
  case SP_TIME:
    if(*cptr && s_get_times(sp, cptr+1))
      return 0;
    break;
  case SP_UVR:
    if(*cptr && s_get_uvr(sp, cptr+1))
      return 0;
    break;
  default:
    lprintf(stderr, "s_get_sel: Unknown selection type.\n");
    return 0;
    break;
  };
  return s_page(sp, S_ALLNEW, 1) < 0;
}

/*.......................................................................
 * Get a new list of comma separated baseline selections from a string,
 * and install them in the specplot attributes descriptor.
 *
 * Input:
 *  sp   Specplot *  The spectrum plot descriptor.
 *  args     char *  A string of arguments.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int s_get_bgl(Specplot *sp, char *args)
{
  Enumpar *epar;  /* A mode symbol table entry */
  char *arg;      /* A selection argument */
/*
 * Are there any arguments?
 */
  if(*args) {
    SpBMode bmode;         /* Baseline selection mode */
    Bgrplist *bgl = NULL;  /* List of baseline groups */
/*
 * Read the baseline selection mode.
 */
    if(!(arg=s_get_arg(args, &args)))
      return 0;
    if((epar=find_enum(sp->sa->bmsym, arg)) == NULL)
      return 0;
    bmode = (SpBMode) epar->id;
/*
 * If there are more arguments, compile them into a list of baseline
 * groups.
 */
    if(*args) {
      bgl = new_Bgrplist();
      if(!bgl)
	return 1;
/*
 * Parse out each baseline group argument from args[] and
 * add each to the list of groups.
 */
      while(*args==',' && (arg=s_get_arg(args+1, &args)) != NULL) {
	if(!add_Basegrp(sp->ob, bgl, NULL, arg)) {
	  bgl = del_Bgrplist(bgl);
	  return 1;
	};
      };
    };
/*
 * Install the new selection.
 */
    if(sp_set_bgl(sp->ob, sp->sa, bmode, bgl))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Get a new list of comma separated polarizations from a string,
 * and install them in the specplot attributes descriptor.
 *
 * Input:
 *  sp   Specplot *  The spectrum plot descriptor.
 *  args     char *  A string of polarizations.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int s_get_pol(Specplot *sp, char *args)
{
  char *arg;  /* A baseline selection list argument */
/*
 * Are there any arguments?
 */
  if(*args) {
/*
 * Get an initially empty list of polarizations.
 */
    Pollist *pl = new_Pollist();
    if(!pl)
      return 1;
/*
 * Parse out each polarization argument from args[] and
 * add each to the list of polarizations.
 */
    for( ; *args && (arg=s_get_arg(args, &args)) != NULL; args++) {
      Stokes pol = Stokes_id(arg);
      if(pol==NO_POL || !add_Polnode(sp->ob, pl, pol)) {
	pl = del_Pollist(pl);
	return 1;
      };
    };
/*
 * Install the new list.
 */
    if(sp_set_pol(sp->ob, sp->sa, pl))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Get a new time range from a string, install it in the specplot
 * attributes descriptor. 
 *
 * Input:
 *  sp   Specplot *  The spectrum plot descriptor.
 *  args     char *  A string containing the new time range.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int s_get_times(Specplot *sp, char *args)
{
  char *arg;       /* An argument string */
  double stime;    /* Start time of range */
  double etime;    /* End time of range */
  double scan;     /* Scan duration or interval */
/*
 * Get the start time.
 */
  if(*args) {
    if((arg=s_get_arg(args, &args)) == NULL ||
       (stime = read_ut(arg, NULL)) < 0)
      return 1;
  } else {
    stime = sp->ob->rec[0].integ->ut;
  };
/*
 * Get the end time.
 */
  if(*args) {
    if((arg=s_get_arg(args+1, &args)) == NULL ||
       (etime = read_ut(arg, NULL)) < 0)
      return 1;
  } else {
    etime = sp->ob->rec[sp->ob->nrec-1].integ->ut;
  };
/*
 * Read the scan description.
 */
  if(*args) {
    char *endp;
    if((arg=s_get_arg(args+1, &args)) == NULL)
      return 1;
    scan = strtod(arg, &endp);
    if(arg == endp) {
      fprintf(stderr, "Unacceptable scan time: %s\n", arg);
      return 1;
    };
    scan *= 60.0;  /* Convert from minutes to seconds */
  } else {
    scan = fabs(etime - stime);
  };
/*
 * There should not be any more arguments.
 */
  if(*args) {
    fprintf(stderr, "Too many arguments.\n");
    return 1;
  };
/*
 * Install the new time range.
 */
  if(sp_set_times(sp->ob, sp->sa, stime, etime, scan))
    return 1;
  return 0;
}

/*.......................................................................
 * Get a new UV range from a string, and install them in the specplot
 * attributes descriptor. 
 *
 * Input:
 *  sp   Specplot *  The spectrum plot descriptor.
 *  args     char *  A string of arguments.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int s_get_uvr(Specplot *sp, char *args)
{
  char *arg;   /* A baseline selection list argument */
  int i;
/*
 * Are there any arguments?
 */
  if(*args) {
    float uv_values[3] = {0.0f, 0.0f, 0.0f};  /* uvmin, uvmax, uvstep */
/*
 * Parse out each UV range argument from args[].
 */
    for(i=0; i<3 && *args && (arg=s_get_arg(args, &args)) != NULL; i++,args++) {
      char *eptr = arg;
      uv_values[i] = uvtowav(strtod(arg, &eptr));
      if(*eptr != '\0') {
	fprintf(stderr, "Error reading argument: %s\n", arg);
	return 1;
      };
    };
/*
 * Install the new list.
 */
    if(sp_set_uvrange(sp->ob, sp->sa, uv_values[0], uv_values[1], uv_values[2]))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Allow the user to specify a new number of plots per page.
 *
 * Input:
 *  sp   Specplot *  The spectrum plot descriptor.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int s_newnum(Specplot *sp)
{
  int nslot;  /* The new number for sp->nslot */
  char *endp; /* A pointer to any trailing input in arg */
  char *arg;
/*
 * Ask the user to input the required new number of plot slots.
 */
  arg = s_getline("Enter the new number of plots per page: ");
  if(!arg)
    return 0;
/*
 * Read the user's reply.
 */
  while(isspace((int)*arg))
    arg++;
  nslot = strtol(arg, &endp, 10);
  while(isspace((int)*endp))
    endp++;
  if(endp == arg) {
    fprintf(stderr, "Not an integer: %s\n", arg);
    return 0;
  };
/*
 * Set limits.
 */
  if(nslot < 0 || nslot > 100) {
    fprintf(stderr, "The number must be between 0 and 100.\n");
    return 0;
  };
/*
 * If zero then leave things as they are.
 */
  if(nslot==0)
    return 0;
/*
 * Allocate more sub-plot descriptors if necessary?
 */
  if(nslot > sp->nsplot && new_SpSubplot(sp, nslot)==NULL)
    return 0;
/*
 * Record the new number of plots per page and re-display the current
 * plot using the new number.
 */
  sp->sa->nplot = sp->nslot = nslot;
  return s_page(sp, S_RESET, 1) <= 0;
}

/*.......................................................................
 * Allocate or resize the array of sub-plot descriptors sp->splot,
 * and default initialize each new descriptor.
 *
 * Input:
 *  sp      Specplot *  The plot descriptor whos sub-plot array is to be
 *                      allocated or resized.
 *  nnew         int    The required number of elements. The existing number
 *                      will be taken from sp->nsplot, unless sp->splots==NULL.
 * Output:
 *  return SpSubplot *  The pointer to newly allocated memory. On error
 *                      NULL will be returned and sp->splots will not be
 *                      changed.
 */
static SpSubplot *new_SpSubplot(Specplot *sp, int nnew)
{
  int nold; /* The existing number of sub-plot elements */
  int i;
/*
 * Check arguments.
 */
  if(nnew <= 0 ) {
    lprintf(stderr, "new_SpSubplot: %d sub-plots requested.\n");
    return NULL;
  };
/*
 * Allocate from scratch?
 */
  if(sp->splots==NULL) {
    nold = 0;
    sp->splots = (SpSubplot *) malloc(sizeof(SpSubplot) * nnew);
    if(!sp->splots) {
      lprintf(stderr, "new_SpSubplot: Insufficient memory.\n");
      return NULL;
    };
    sp->nsplot = nnew;
/*
 * Resize an existing array?
 */
  } else if(nnew > sp->nsplot) {
    SpSubplot *sps = (SpSubplot *) realloc(sp->splots, sizeof(SpSubplot)*nnew);
    if(!sps) {
      lprintf(stderr, "Insufficient memory for %d plots.\n", nnew);
      return NULL;
    };
    nold = sp->nsplot;
    sp->splots = sps;
    sp->nsplot = nnew;
/*
 * There are already enough sub-plot elements.
 */
  } else {  /* nnew <= sp->nslot */
    nold = sp->nslot;
  };
/*
 * Initialize each of the new sub-plot descriptors.
 */
  for(i=nold; i<nnew; i++) {
    SpSubplot *sps = sp->splots + i;
    sps->spec = NULL;
    sps->vya = sps->vyb = 0.0f;
    sps->vymid = 0.0f;
    sps->amin = sps->amax = 0.0f;
    sps->pmin = sps->pmax = 0.0f;
    sps->spa.uta = sps->spa.utb = 0;
    sps->spa.isub = sps->spa.base = 0;
    sps->spa.isel = 0;
    sps->spa.ipol = 0;
    sps->spa.iuv = 0;
  };
  return sp->splots;
}

/*.......................................................................
 * Write a baseline selection title into a string of given length.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 *  spa     SpAttr *  The attributes of the sub-plot being labelled.
 *  nc         int    The max number of characters (including '\0') to
 *                    place in s[].
 *  s         char *  The string to write the title into.
 * Output:
 *  return     int    The number of characters written (excluding '\0'.
 *                    If the title didn't fit it won't have been written
 *                    and 0 will be returned.
 */
static int s_bgrp_title(Specplot *sp, SpAttr *spa, int nc, char *s)
{
  int nused = 0;   /* The number of characters written */
  int nnew;        /* The number of characters to be written next */
  char *prefix;    /* Title prefix */
  Basegrp *bgrp;
/*
 * Get the baseline group to be described.
 */
  if((bgrp = s_base_iter(sp, 0, 0, 1, spa))) {
/*
 * Work out the title prefix.
 */
    prefix = sp->sa->bmode==SP_GROUP ? "Baseline group " : "Baseline ";
/*
 * Write the prefix.
 */
    nnew = strlen(prefix);
    if(nnew < nc-nused) {
      sprintf(s, prefix);
      nused += nnew;
/*
 * Write the baseline selection.
 */
      nnew = write_Basegrp(sp->ob, bgrp, nc-nused, s + nused);
      if(nnew >= 0) {
	nused += nnew;
	return nused;
      };
    };
  };
/*
 * String too short.
 */
  if(nc>0)
    s[0] = '\0';
  return 0;
}

/*.......................................................................
 * Write a polarization selection title into a string of given length.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 *  prefix    char *  An optional prefix to write before the specification.
 *  spa     SpAttr *  The attributes of the sub-plot being labelled.
 *  nc         int    The max number of characters (including '\0') to
 *                    place in s[].
 *  s         char *  The string to write the title into.
 * Output:
 *  return     int    The number of characters written (excluding '\0'.
 *                    If the title didn't fit it won't have been written
 *                    and 0 will be returned.
 */
static int s_pol_title(Specplot *sp, SpAttr *spa, int nc, char *s)
{
  char *pname;     /* The name of the polarization */
  int nused = 0;   /* The number of characters written */
  int nnew;        /* The number of characters to be written next */
  Stokes pol;      /* The stokes parameter enumeration id */
  char *prefix = "Polarization ";
/*
 * Get the stokes parameter id.
 */
  if((pol = s_pol_iter(sp, 0, 0, 1, spa)) != NO_POL) {
/*
 * Write a prefix?
 */
    nnew = strlen(prefix);
    if(nnew < nc-nused) {
      sprintf(s, prefix);
      nused += nnew;
/*
 * Write the polarization selection.
 */
      pname = Stokes_name(pol);
      nnew = strlen(pname);
      if(nnew+8 < nc-nused) {
	sprintf(s+nused, "\\fr%s\\fn", pname);
	return strlen(s);
      };
    };
  };
/*
 * String too short.
 */
  if(nc>0)
    s[0] = '\0';
  return 0;
}

/*.......................................................................
 * Write a UV radius range selection title into a string of given length.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 *  prefix    char *  An optional prefix to write before the specification.
 *  spa     SpAttr *  The attributes of the sub-plot being labelled.
 *  nc         int    The max number of characters (including '\0') to
 *                    place in s[].
 *  s         char *  The string to write the title into.
 * Output:
 *  return     int    The number of characters written (excluding '\0'.
 *                    If the title didn't fit it won't have been written
 *                    and 0 will be returned.
 */
static int s_uvr_title(Specplot *sp, SpAttr *spa, int nc, char *s)
{
  char awrk[80];     /* Work string */
  int nused = 0;     /* The number of characters written */
  int nnew;          /* The number of characters to be written next */
  float uvmin,uvmax; /* The UV range limits (wavelengths) */
  char *prefix = "UV range ";
/*
 * Get the current UV radius range.
 */
  if(s_uvr_iter(sp, 0, 0, 1, spa, &uvmin, &uvmax)==0) {
/*
 * Write a prefix?
 */
    nnew = strlen(prefix);
    if(nnew < nc-nused) {
      sprintf(s, prefix);
      nused += nnew;
/*
 * Write the range selection.
 */
      sprintf(awrk, "%g -> %g (%s)", wavtouv(uvmin), wavtouv(uvmax),
	      uvwunits(U_PLAB));
      nnew = strlen(awrk);
      if(nnew+8 < nc-nused) {
	sprintf(s+nused, "%s", awrk);
	return strlen(s);
      };
    };
  };
/*
 * String too short.
 */
  if(nc>0)
    s[0] = '\0';
  return 0;
}

/*.......................................................................
 * Write a time-range selection title into a string of given length.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 *  spa     SpAttr *  The attributes of the sub-plot being labelled.
 *  nc         int    The max number of characters (including '\0') to
 *                    place in s[].
 *  s         char *  The string to write the title into.
 * Output:
 *  return     int    The number of characters written (excluding '\0'.
 *                    If the title didn't fit it won't have been written
 *                    and 0 will be returned.
 */
static int s_time_title(Specplot *sp, SpAttr *spa, int nc, char *s)
{
  int nused = 0;     /* The number of characters written */
  int nnew;          /* The number of characters to be written next */
  int uta,utb;       /* The indexes of the two times */
  char *prefix = "Time ";
/*
 * Get the integration record indexes.
 */
  uta = spa->uta;
  utb = spa->utb;
/*
 * Write a prefix?
 */
  nnew = strlen(prefix);
  if(nnew < nc-nused) {
    sprintf(s, prefix);
    nused += nnew;
/*
 * Write the time-range.
 */
    nnew = write_ut(sp->ob->rec[uta].integ->ut, nc-nused, s + nused);
    if(nnew > 0 && nnew < nc-nused) {
      nused += nnew;
/*
 * One time or two?
 */
      if(uta != utb && 3 < nc-nused) {   /* Time range */
	strcat(s + nused, " - ");
	nused += 3;
	if(nc-nused > 0) {
	  nnew = write_ut(sp->ob->rec[utb].integ->ut, nc-nused, s + nused);
	  if(nnew > 0) {
	    return strlen(s);
	  };
	};
      } else {                           /* Single time */
	return strlen(s);
      };
    };
  };
/*
 * String too short.
 */
  if(nc>0)
    s[0] = '\0';
  return 0;
}

/*.......................................................................
 * Calculate regularly gridded values of the current smoothing function
 * and record them in sp->smfn[] and their number in sp->nsmth.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor that records the
 *                    details of the required smoothing function.
 *  cif        int    The IF to which the smoothing function is to be
 *                    calculated for.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int s_get_smfn(Specplot *sp, int cif)
{
  SmType smooth_type = sp->sa->smooth.type;
  float xmul;   /* The scale factor to convert the fwhm to channels */
  float fwhm;   /* The full-width at half maxmimum in channels */
  float *smfn = sp->smfn;
  int i;
/*
 * Get the full-width at half maximum of the smoothing function.
 */
  if(s_coords(sp, cif, sp->sa->smooth.xunit, NULL, &xmul))
    return 1;
  fwhm = xmul != 0.0f ? fabs(sp->sa->smooth.fwhm/xmul) : 0.0f;
/*
 * No smoothing selected?
 */
  if(smooth_type==SM_NONE || fwhm<=0.0f) {
    sp->nsmth = 1;
    sp->smfn[0] = 1.0f;
  } else {
/*
 * Smoothing a data array involves convolving it with a chosen
 * smoothing function s(x). This could be acheived using FFTs but this
 * is complicated by propogation of data weights and due to potential
 * aliasing such that it would require 5 large complex arrays, for the FFTs
 * of weight, data*weight, smoothing_function, smoothing_function^2,
 * and an extra intermediate array. Note that analytic Fourier transforms of
 * the smoothing function are not appropriate (I have tried this).
 *
 * So instead, the convolution will be done directly. Each smoothing
 * function s(x) will be described in terms of its full-width
 * at half maximum value, denoted as fwhm. For convenience the
 * dependent coordinate of each function will be defined such that it is
 * unity where the function falls to half its maximum value. In terms of
 * the spectrum channel--based coordinates, i, this requires the
 * following change of coordinates.
 *
 *  dx = 2.di/fwhm
 */
    float dx = 2.0f/fwhm;
/*
 * How many elements should be initialized with the smoothing function.
 */
    sp->nsmth = fwhm * nsigma / dx;
    if(sp->nsmth < 1)
      sp->nsmth = 1;
    else if(sp->nsmth > sp->nchan)
      sp->nsmth = sp->nchan;
/*
 * Calculate the requested smoothing function values and record them
 * in sp->smfn[0..sp->nsmth-1].
 */
    switch(smooth_type) {
/*
 * Sinc smoothing function: s[x] = B/pi * sin(B.x)/(B.x)
 * Where B = 1.89549426703340.
 */
    case SM_SINC:
      {
	const double b = 1.89549426703340;
	for(i=0; i<sp->nsmth; i++) {
	  float bx = b * i * dx;
	  smfn[i] = b/pi * (bx > 1e-18f ? sin(bx)/bx : 1.0f);
	};
      };
      break;
/*
 * Hanning window:  s[x] = sin(pi.x) / (2.pi.x.(1-x^2))
 */
    case SM_HANNING:
      for(i=0; i<sp->nsmth; i++) {
	float x = i * dx;
	if(x < 1e-18f)
	  smfn[i] = 0.5f;
	else if(fabs(x - 1.0f) < 1e-18f)
	  smfn[i] = 0.25f;
	else
	  smfn[i] = sin(pi * x) / (twopi * x * (1.0f - x*x));
      };
      break;
/*
 * Gaussian window function: sqrt(ln(2)/pi).exp(-ln(2).x^2).
 */
    case SM_GAUSSIAN:
      {
	float ln2 = log(2.0);
	float gscale = sqrt(ln2/pi);
	for(i=0; i<sp->nsmth; i++) {
	  float x = i * dx;
	  smfn[i] = x < 5.0f ? gscale * exp(-ln2 * x * x) : 0.0f;
	};
      };
      break;
/*
 * Boxcar window function: s[x] = / 0.5    |x| <= 1
 *                                \ 0.0    |x|  > 1
 */
    case SM_BOXCAR:
      {
	for(i=0; i<sp->nsmth; i++) {
	  float x = i * dx;
	  smfn[i] = (x <= 1.0f) ? 0.5f : 0.0f;
	};
      };
      break;
    default:
      lprintf(stderr, "s_get_smfn: Unrecognised smoothing function type.\n");
      return 1;
    };
  };
  return 0;
}

/*.......................................................................
 * Return the linear factors relating a given X-axis coordinate and
 * channels.
 *
 * Displayed coordinate = xoff + channel_index * xmul.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 *  cif        int    The IF to characterise.
 *  xunit  SpXunit    The type of X-axis coordinate.
 * Input/Output:
 *  xoff     float *  If xoff!=NULL, then *xoff will be assigned the
 *                    coordinate of the first channel.
 *  xmul     float *  If xmul!=NULL, then *xmul will be assigned the
 *                    coordinate increment between adjacent channels.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error *xoff and *xmul will be assigned 0.0f
 *                        and 1.0f respectively.
 */
static int s_coords(Specplot *sp, int cif, SpXunit xunit,
		    float *xoff, float *xmul)
{
  SpXdim *spx = sp->spx + cif;
  If *ifs = sp->ob->ifs + cif;
  float s_xoff;
  float s_xmul;
  int waserr = 0;
/*
 * Get the required transformation factors.
 */
  switch(xunit) {
  case SP_CHAN:
    s_xoff = spx->cmin + 1.0f;
    s_xmul = 1.0f;
    break;
  case SP_FREQ:
    s_xoff = ifs->freq * 1.0e-9f;
    s_xmul = ifs->df * 1.0e-9f;
    break;
  default:
    lprintf(stderr, "s_coords: X-axis type not recognized.\n");
    waserr = 1;
    s_xoff = 0.0f;
    s_xmul = 1.0f;
    break;
  };
/*
 * Assign the return values.
 */
  if(xoff)
    *xoff = s_xoff;
  if(xmul)
    *xmul = s_xmul;
  return waserr;
}

/*.......................................................................
 * Set and/or return the current specplot attribute flags.
 *
 * Input:
 *  sa   Specattr *   The specplot attributes descriptor.
 *  flags    char *   The string of flags to install, or NULL if this
 *                    call is being made solely to query the current
 *                    set of flags.
 * Output:
 *  return   char *   A string containing the current set of flags, or
 *                    NULL on error.
 */
char *sp_set_flags(Specattr *sa, char *flags)
{
  enum {MAX_FLAG = 20};         /* Max number of flags */
  static char fstr[MAX_FLAG+1]; /* Return string of flags */
  int i;
/*
 * Check the attributes descriptor.
 */
  if(bad_Specattr(sa, "sp_flags"))
    return NULL;
/*
 * Install a new set of flags?
 */
  if(flags) {
    char *flag;
/*
 * Start with everything turned off.
 */
    sa->doamp = 0;
    sa->dophs = 0;
    sa->docross = 0;
    sa->dojoin = 0;
    sa->dohist = 0;
    sa->dobars = 0;
/*
 * Now obey the flag string.
 */
    for(flag=flags; *flag; flag++) {
      int waslow = islower((int)*flag);
      int key = waslow ? toupper(*flag) : *flag;
      s_flags(sa, key, waslow);
/*
 * Handle flags that are normally treated as commands.
 */
      switch(key) {
      case KEY_CROSS:
	sa->docross = !sa->docross;
	break;
      default:
	break;
      };
    };
/*
 * Ensure that at least one of doamp and dophs is true.
 */
    if(!sa->doamp && !sa->dophs)
      sa->doamp = sa->dophs = 1;
  };
/*
 * Compile a string from the current set of flags.
 */
  i = 0;
  if(sa->doamp || sa->dophs)
    fstr[i++] = sa->doamp ? (sa->dophs ? KEY_BOTH : KEY_AMP) : KEY_PHS;
  if(sa->docross)
    fstr[i++] = '+';
  if(sa->dojoin)
    fstr[i++] = 'j';
  if(sa->dohist)
    fstr[i++] = 'J';
  if(sa->dobars)
    fstr[i++] = 'e';
  fstr[i] = '\0';
  return fstr;
}

/*.......................................................................
 * Iterate over baselines.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  reset      int     If true start searching from the start or end
 *                     of the range of baselines, according to the direction
 *                     implied by 'forward'.
 *  advance    int     If true, the current position is to be skipped.
 *  forward    int     The search direction.
 * Input/Output:
 *  spa     SpAttr *   The input/output iterator ordinal state.
 * Output:
 *  return Basegrp *   A temporary pointer to a baseline group descriptor
 *                     describing the new state, or NULL if not found.
 */
static Basegrp *s_base_iter(Specplot *sp, int reset, int advance, int forward,
			     SpAttr *spa)
{
  Bgrplist *bgl;  /* The list of baseline groups */
  int i;
/*
 * Get the current list of baseline groups.
 */
  bgl = sp->sa->bgl;
/*
 * Get the baseline selection group appropriate to the plot mode.
 */
  switch(sp->sa->bmode) {
  case SP_GROUP:
    {
      int isel;
      if(reset) {
	isel = forward ? 0 : (bgl->nsel-1);
      } else {
	isel = spa->isel + (advance ? (forward ? 1 : -1) : 0);
      };
      if(isel < 0 || isel >= bgl->nsel) {
	return NULL;
      } else {
	Basegrp *bgrp;
/*
 * Find the isel'th group.
 */
	for(i=0,bgrp=bgl->bgrp; bgrp; bgrp = bgrp->next,i++) {
	  if(i == isel) {
	    spa->isel = isel;
	    return bgrp;
	  };
	};
      };
    };
    break;
  case SP_SPLIT:
    {
      Basegrp *bgrp = bgl->bgrp;   /* The baseline selection to use */
      int found = 0;               /* True if a valid baseline is located */
      int isub, base;              /* Sub-array and baseline indexes */
/*
 * Get the sub-array and baseline to search from.
 */
      if(reset) {
	if(forward) {
	  isub = base = 0;
	} else {
	  isub = sp->ob->nsub - 1;
	  base = sp->ob->sub[isub].nbase - 1;
	};
	found = in_Basegrp(sp->ob, isub, base, bgrp) ||
	        srch_Basegrp(sp->ob, bgrp, forward, &isub, &base);
      } else {
	isub = spa->isub;
	base = spa->base;
/*
 * Locate a suitable baseline.
 */
	found = advance ? srch_Basegrp(sp->ob, bgrp, forward, &isub, &base) :
	                  in_Basegrp(sp->ob, isub, base, bgrp);
      };
/*
 * If a valid baseline was found, initialize sp->scr_bgrp to select it.
 */
      if(found) {
	Baseline *bptr = sp->ob->sub[isub].base + base;
	Basespec *bs = find_base(sp->ob, 3, isub, bptr->tel_a, bptr->tel_b,
				 1, 3, 1, 1, 1);
	if(!bs || !clr_Basegrp(sp->scr_bgrp) ||
	   !add_Basesel(sp->ob, sp->scr_bgrp, bs, 1))
	  return NULL;
	spa->isub = isub;
	spa->base = base;
	return sp->scr_bgrp;
      };
    };
    break;
  default:     /* Take the first baseline group from the list */
    lprintf(stderr, "base_iter: Bad baseline mode.\n");
    return NULL;
    break;
  };
  return NULL;
}

/*.......................................................................
 * Iterate over time ranges.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  reset      int     If true start searching from the start or end
 *                     of the ranges of times, according to the direction
 *                     implied by 'forward'.
 *  advance    int     If true, the current position is to be skipped.
 *  forward    int     The search direction.
 * Input/Output:
 *  spa     SpAttr *   The input/output iterator ordinal state.
 * Output:
 *  return     int     0 - ok.
 *                     1 - Not found.
 */
static int s_time_iter(Specplot *sp, int reset, int advance, int forward,
		       SpAttr *spa)
{
  Observation *ob; /* The descriptor of the observation being plotted */
  Specattr *sa;    /* Specplot plot attributes */
  double stime;    /* Start time of slot */
  double etime;    /* End time of slot */
  double scan;     /* Scan interval or duration */
  int ut_a;        /* The index of the lowest integration >= stime */
  int ut_b;        /* The index of the highest integration <= etime */
/*
 * Get the specplot delimiting attributes and the observation descriptor.
 */
  sa = sp->sa;
  ob = sp->ob;
/*
 * Get the scan interval/duration.
 */
  scan = fabs(sp->sa->scan);
/*
 * Use the existing time range?
 */
  if(!reset && !advance) {
    ut_a = spa->uta;
    ut_b = spa->utb;
  } else {
/*
 * Find the time range to start the search from.
 */
    if(reset) {
      ut_a = ut_b = forward ? 0 : ob->nrec - 1;
    } else {
      ut_a = ut_b = forward ? (spa->utb + 1) : (spa->uta - 1);
    };
/*
 * Search for the next time slot in the forward direction?
 */
    if(forward) {
/*
 * Find the start integration of the new time slot.
 */
      if(ut_a >= ob->nrec)
	return 1;
      if(ut_a < 0 || ob->rec[ut_a].integ->ut < sa->stime)
	ut_a = ob_find_ut(ob, sa->stime, UT_GE);
      if(ut_a < 0 || ut_a >= ob->nrec)
	return 1;
      stime = ob->rec[ut_a].integ->ut;
/*
 * Find the end of the time slot.
 */
      if(sp->sa->scan < 0.0f) {   /* Treat scan as a scan separator */
	for(ut_b=ut_a;
	    (ut_b<ob->nrec-1 &&
	     ob->rec[ut_b+1].integ->ut - ob->rec[ut_b].integ->ut < scan);
	    ut_b++)
	  ;
      } else {                    /* Treat scan as a scan duration */
	etime = stime + scan;
	if(etime > sa->etime)
	  etime = sa->etime;
	if(stime > etime)
	  return 1;
	ut_b = ob_find_ut(ob, etime, UT_LE);
      };
/*
 * Search for the next time slot in the backward direction?
 */
    } else {
/*
 * Find the end integration of the new time slot.
 */
      if(ut_b < 0)
	return 1;
      if(ut_b >= ob->nrec || ob->rec[ut_b].integ->ut > sa->etime)
	ut_b = ob_find_ut(ob, sa->etime, UT_LE);
      if(ut_b < 0 || ut_b >= ob->nrec)
	return 1;
      etime = ob->rec[ut_b].integ->ut;
/*
 * Find the start of the time slot.
 */
      if(sp->sa->scan < 0.0f) {   /* Treat scan as a scan separator */
	for(ut_a=ut_b;
	    (ut_a>0 &&
	     ob->rec[ut_a].integ->ut - ob->rec[ut_a-1].integ->ut < scan);
	    ut_a--)
	  ;
      } else {                   /* Treat scan as a scan period */
	stime = etime - scan;
	if(stime < sa->stime)
	  stime = sa->stime;
	if(stime > etime)
	  return 1;
	ut_a = ob_find_ut(ob, stime, UT_GE);
      };
    };
  };
/*
 * Did we acquire a valid time slot?
 */
  if(ut_a > ut_b ||
     (ut_a < 0 || ut_a >= ob->nrec || ob->rec[ut_a].integ->ut > sa->etime) ||
     (ut_b < 0 || ut_b >= ob->nrec || ob->rec[ut_b].integ->ut < sa->stime))
    return 1;
/*
 * Record the new time slot limits.
 */
  spa->uta = ut_a;
  spa->utb = ut_b;
  return 0;
}

/*.......................................................................
 * Iterate over polarizations.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  reset      int     If true start searching from the start or end
 *                     of the list of polarizations, according to the
 *                     direction implied by 'forward'.
 *  advance    int     If true, the current position is to be skipped.
 *  forward    int     The search direction.
 * Input/Output:
 *  spa     SpAttr *   The input/output iterator ordinal state.
 * Output:
 *  return  Stokes     The new polarization id, or NO_POL if not found.
 */
static Stokes s_pol_iter(Specplot *sp, int reset, int advance, int forward,
			 SpAttr *spa)
{
  Observation *ob = sp->ob; /* The observation being plotted */
  Pollist *pl=sp->sa->pl;   /* The list of polarizations */
  Stokes pol = NO_POL;      /* The id of the selected polarization */
  Obpol obpol;              /* The descriptor of the selected polarization */
  int ipol;                 /* The index of the new polarization */
  int npol;                 /* The number of polarizations in the list */
  int i;
/*
 * Determine the number of available polarizations.
 */
  npol = pl ? pl->npol : 1;
/*
 * Get the ordinal polsition of the next polarization in the list.
 */
  if(reset) {
    ipol = forward ? 0 : (npol - 1);
  } else {
    ipol = spa->ipol + (advance ? (forward ? 1 : -1) : 0);
  };
/*
 * Has a valid polarization been selected?
 */
  if(ipol < 0 || ipol >= npol)
    return NO_POL;
/*
 * Locate the selected polarization slot.
 */
  if(pl) {
    Polnode *pn;
    for(i=0, pn=pl->head; pol==NO_POL && pn; i++, pn=pn->next) {
      if(i==ipol)
	pol = pn->pol;
    };
    if(pol==NO_POL)
      return NO_POL;
  };
/*
 * Check the polarization selected.
 */
  if(get_Obpol(ob, pol, 1, &obpol)==0) {
    spa->ipol = ipol;
    return obpol.type;
  };
  return NO_POL;
}

/*.......................................................................
 * Iterate over UV radius annuli.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  reset      int     If true start searching from the start or end
 *                     of the ranges of UV radii, according to the
 *                     direction implied by 'forward'.
 *  advance    int     If true, the current position is to be skipped.
 *  forward    int     The search direction.
 * Input/Output:
 *  spa     SpAttr *   The input/output iterator ordinal state.
 *  uvmin    float *   If uvmin!=NULL, then *uvmin will contain the
 *                     minimum UV radius (wavelengths).
 *  uvmax    float *   If uvmax!=NULL, then *uvmax will contain the
 *                     maximum UV radius (wavelengths).
 * Output:
 *  return     int     0 - OK.
 *                     1 - Not found.
 */
static int s_uvr_iter(Specplot *sp, int reset, int advance, int forward,
			 SpAttr *spa, float *uvmin, float *uvmax)
{
  SpUV *uvr;  /* UV range container */
  int iuv;    /* Iterator ordinal position index */
  int nuv;    /* The number of annulus ranges */
/*
 * Get the UV range descriptor.
 */
  uvr = &sp->sa->uvr;
/*
 * If the UV range has not been initialized do so now.
 */
  if(uvr->uvrlim < 0.0f && sp_set_uvrange(sp->ob, sp->sa, 0.0f, 0.0f, 0.0f))
    return 1;
/*
 * Determine the number of available radius ranges.
 */
  if(uvr->uvstep <= 0.0f)
    nuv = 1;
  else 
    nuv = abs(ceil((uvr->uvmax - uvr->uvmin) / uvr->uvstep));
/*
 * Get the ordinal polsition of the next UV radius annulus.
 */
  if(reset) {
    iuv = forward ? 0 : (nuv - 1);
  } else {
    iuv = spa->iuv + (advance ? (forward ? 1 : -1) : 0);
  };
/*
 * Has a valid range been selected?
 */
  if(iuv < 0 || iuv >= nuv)
    return 1;
/*
 * Install the new oridnal position and return the new range.
 */
  spa->iuv = iuv;
  if(uvmin)
    *uvmin = uvr->uvmin + uvr->uvstep * spa->iuv;
  if(uvmax)
    *uvmax = uvr->uvmin + uvr->uvstep * (1 + spa->iuv);
  return 0;
}

/*.......................................................................
 * Plot an index key title.
 *
 * Input:
 *  sp    Specplot *  The spectrum plot descriptor.
 *  key      SpKey    The key to plot a title for.
 *  spa     SpAttr *  The attributes of the sub-plot being labelled.
 *  nc         int    The max number of characters (including '\0') to
 *                    place in s[].
 *  s         char *  The string to write the title into.
 * Output:
 *  return     int    The number of characters written (excluding '\0'.
 *                    If the title didn't fit it won't have been written
 *                    and 0 will be returned.
 */
static int s_title(Specplot *sp, SpKey key, SpAttr *spa, int nc, char *s)
{
  int nused = 0;
/*
 * Dispatch to the appropriate function.
 */
  switch(key) {
  case SP_BASE:
    nused = s_bgrp_title(sp, spa, nc, s);
    break;
  case SP_TIME:
    nused = s_time_title(sp, spa, nc, s);
    break;
  case SP_POL:
    nused = s_pol_title(sp, spa, nc, s);
    break;
  case SP_UVR:
    nused = s_uvr_title(sp, spa, nc, s);
    break;
  default:
    s[0] = '\0';
    nused = 0;
    break;
  };
  return nused;
}

/*.......................................................................
 * Iterate a given key.
 *
 * Input:
 *  sp    Specplot *   The spectrum plot descriptor.
 *  key      SpKey     The key to iterate.
 *  reset      int     If true start searching from the start or end
 *                     of the range of items, according to the direction
 *                     implied by 'forward'.
 *  advance    int     If true, the current position is to be skipped.
 *  forward    int     The search direction.
 * Input/Output:
 *  spp   Specposn *   The descriptor of the located item if found.
 * Output:
 *  return     int     0 - ok.
 *                     1 - Not found.
 */
static int s_iterate(Specplot *sp, SpKey key, int reset, int advance,
		     int forward, Specposn *spp)
{
/*
 * Dispatch to the appropriate iterator function.
 */
  switch(key) {
  case SP_BASE:
    spp->bgrp = s_base_iter(sp, reset, advance, forward, &spp->spa);
    if(!spp->bgrp)
      return 1;
    break;
  case SP_TIME:
    if(s_time_iter(sp, reset, advance, forward, &spp->spa))
      return 1;
    spp->uta = spp->spa.uta;
    spp->utb = spp->spa.utb;
    break;
  case SP_POL:
    spp->pol = s_pol_iter(sp, reset, advance, forward, &spp->spa);
    if(spp->pol==NO_POL)
      return 1;
    break;
  case SP_UVR:
    if(s_uvr_iter(sp, reset, advance, forward, &spp->spa,
		  &spp->uvmin, &spp->uvmax))
      return 1;
    break;
  default:
    lprintf(stderr, "s_iterate: Unknown iterator code.\n");
    return 1;
    break;
  };
  return 0;
}


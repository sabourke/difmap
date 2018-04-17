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

typedef struct {
  float vxa,vxb;  /* Min/max NDC X-coords of scan sub-plot */
  float sutmin;   /* UT range in scan is utmin -> utmax */
  float sutmax;
  float utmin;    /* The UT range to be displayed from this scan */
  float utmax;
  int view;       /* Flags whether any of scan is visible */
} Scans;

/* Declare a type used to store details about the last cursor input */

typedef struct {
  int key;         /* The upper-case version of the cursor key selected */
  int waslow;      /* True if the selected cursor key was lower case */
  int wasamp;      /* True if 'value' denotes an amplitude */
  float utval;     /* The time coordinate selected */
  float value;     /* The amplitude [or phase if wasamp==0] selected */
  Scans *sc;       /* The descriptor of the scan selected */
} Ccurs;

typedef struct {
  Telspec ts;      /* The specification of the telescope being plotted */
  double utref;    /* Reference UT - seconds */
  Observation *ob; /* The descriptor of the observation being plotted */
  Subarray *sub;   /* The descriptor of sub-array ts.isub in ob->sub[] */
  float utmin;     /* World min X coordinate (seconds wrt utref) */
  float utmax;     /* World max X coordinate (seconds wrt utref) */
  float utsum;     /* Sum of scan UT ranges currently visible */
  float vxa,vxb;   /* NDC coords of viewport surrounding grid of sub-plots */
  float vya,vyb;
  float vymid;     /* NDC coord of line separating amplitude and phase plots */
  float ampmin;    /* Minimum amplitude plotted */
  float ampmax;    /* Maximum amplitude plotted */
  int modified;    /* This remains 0 unless the data are edited */
  int uta,utb;     /* The indexes of the first and last plotted integrations */
  int cif;         /* The raw index of the current IF */
  int docurs;      /* True when cursor control is in effect */
  int docross;     /* True to enable cross-hair mode */  
  int doscan;      /* True if the plot is to be separated into scans */
  int nscan;       /* Number of scans in 'scans' array */
  Scans *scans;    /* Array of scan descriptors */
  Ccurs cursor;    /* Cursor entry descriptor */
  int npage;       /* The sequential number of the page being plotted */
} Corpar;

/* Set keys for interactive display-editing */

enum {
  KEY_NONE='\0',  /* Null key press */
  KEY_HELP='H',   /* Key to request help */
  KEY_QUIT='X',   /* Quit */
  KEY_CUR ='A',   /* Edit a point */
  KEY_NEXT='N',   /* Next telescope */
  KEY_PREV='P',   /* Prev telescope */
  KEY_TEL ='T',   /* Telescope entry via the keyboard */
  KEY_UT  ='U',   /* Select UT range via the cursor */
  KEY_CAN ='D',   /* Key to cancel UT range */
  KEY_DIS ='L',   /* Redisplay the plot */
  KEY_BRK ='B',   /* Toggle separating plot into scans */
  KEY_PRVIF='[',  /* Show the previous IF */
  KEY_NXTIF=']',  /* Show the next IF */
  KEY_CROSS='+'   /* Toggle cross-hair cursor mode */
};

static const float ymarg=0.1;  /* The fraction of the Y range for margin */
static const float xmarg=0.05; /* The fraction of the X range for margin */
static const int corcol=10;    /* The color of unflagged corrections */
static const int badcol=11;    /* The color of flagged corrections */
static const int zoomcol=5;    /* PGPLOT color index for zoom cursor window */
static const int corsym=1;     /* The PGPLOT marker for unflagged corrections */
static const int badsym=1;     /* The PGPLOT marker for flagged corrections */

/* Plotting functions */

static Corpar *new_Corpar(Observation *ob, Telspec *ts, int cif, int docurs,
			  int doscan);
static Corpar *del_Corpar(Corpar *cp);
static int set_Scans(Corpar *cp);
static int c_pldata(Corpar *cp, int uta, int utb, int erase);
static int c_label(Corpar *cp);
static int c_redisp(Corpar *cp);
static int c_scale(Corpar *cp, int doamp, float *xtomm, float *ytomm);
static float c_labinc(float tspan, int ntry);
static int c_arange(Corpar *cp);
static int c_utrange(Corpar *cp);
static int c_vpwin(Corpar *cp);
static int c_plaxes(Corpar *cp, int erase);

typedef enum {C_ALLNEW, C_NXT_ISUB, C_NXT_TA, C_NXT_TEL} Telop;

static int c_newtel(Corpar *cp, Telop oper, int forward, int report,
		    Telspec *init);

/* Editing functions */

static int c_newut(Corpar *cp);
static int c_find(Corpar *cp, float utval, float yval,
		  int wasamp);
static int c_flags(Corpar *cp, char key, int waslow);

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

static int c_cursor(Corpar *cp, int noout, Bandmode mode, int isamp,
		    float xref, float yref, int ci);

/*
 * Enumerate correction editing modes.
 */
typedef enum {ED_RESET, ED_FLAG} Edmode;
static int c_edit_cor(Corpar *cp, int ut, Edmode mode);

/*.......................................................................
 * Create a new Corpar descriptor and an array of scans, then assign
 * given defaults.
 *
 * Input:
 *  ob Observation *  The data descriptor.
 *  ts     Telspec *  The specification of the first telescope to plot,
 *                    or NULL (or an empty specification) for the default.
 *  cif        int    The index of the start IF, or -1 for the default.
 *  docurs     int    True if cursor control is required. If the current
 *                    device has no cursor, this will be ignored.
 * Output:
 *  return  Corpar *  The new Corpar descriptor, or NULL on error.
 */
static Corpar *new_Corpar(Observation *ob, Telspec *ts, int cif, int docurs,
			  int doscan)
{
  Corpar *cp;      /* Pointer to the new Corpar descriptor */
  char answer[10]; /* Answer from pgqinf() */
  int slen;        /* Length of string */
/*
 * If no telescope specification was provided, substitute a default.
 */
  if((ts && next_tel(ob, FIND_FIRST, 1, 0, 0, 1, ts) != 0) ||
     (ts==NULL && (ts=find_tel(ob, 0, 0, 0, 1, 0, 0, 1))==NULL))
    return NULL;
/*
 * An IF index of -1 (0 on the command line) requests the default IF,
 * substitute the first IF.
 */
  if(cif == -1)
    cif = 0;
  if(cif < 0 || cif >= ob->nif) {
    lprintf(stderr, "corplot: IF %d does not exist.\n", cif+1);
    return NULL;
  };
/*
 * Allocate the Corpar descriptor.
 */
  cp = (Corpar *) malloc(sizeof(Corpar));
  if(cp == NULL) {
    lprintf(stderr, "new_Corpar: Insufficient memory for plot descriptor.\n");
    return NULL;
  };
/*
 * NULLify pointer members so that del_Corpar can be called before the
 * struct has been fully initialized.
 */
  cp->scans = NULL;
  cp->sub = NULL;
  cp->ob = ob;
/*
 * Record the initial telescope specification.
 */
  cp->ts = *ts;
/*
 * The data have not been modified yet.
 */
  cp->modified = 0;
/*
 * Assign starting values to the rest of the members.
 */
  cp->nscan = 0;
  cp->doscan = doscan;
  cp->utmin = cp->utmax = 0.0f;
  cp->uta = cp->utb = 0;
  cp->utref = ob->date.ut;
  cp->cif = cif;
/*
 * If cursor interaction is required, check if the device has a cursor.
 */
  if(docurs) {
    slen = sizeof(answer)-1;
    cpgqinf("CURSOR", answer, &slen);
    docurs = strncmp(answer,"YES",3) == 0;
  };
  cp->docurs = docurs;
  cp->cursor.key = KEY_NONE;
  cp->docross = 0;
  cp->npage = 0;
/*
 * Initialize for the first sub-array.
 */
  cp->sub = ob->sub + ts->isub;
/*
 * Set up scan info for the new sub-array.
 */
  if(set_Scans(cp)==0)
    return del_Corpar(cp);
/*
 * Assign starting values to the rest of the members.
 */
  cp->uta = 0;
  cp->utb = cp->sub->ntime-1;  /* Default to show all data */
/*
 * Return the new descriptor.
 */
  return cp;
}

/*.......................................................................
 * Corpar (visibility plot descriptor) destructor function.
 *
 * Input:
 *  cp     Corpar *  Corpar pointer returned by new_Corpar().
 * Output:
 *  return Corpar *  Always NULL, so that you can write cp=del_Corpar(cp);
 */
static Corpar *del_Corpar(Corpar *cp)
{
  if(cp) {
    if(cp->scans)
      free(cp->scans);
    free(cp);
  };
  return NULL;
}

/*.......................................................................
 * Determine a new set of scans from a new time separator.
 *
 * Input:
 *  cp         Corpar *  The plot parameter container.
 *                        cp->doscan: If 0 the whole observation is
 *                                    treated as one scan.
 *                        cp->scans:  Must be NULL on first call.
 *                        cp->nscan:  Must be 0 on first call.
 * Input/Output:
 *  cp         Corpar *  The plot descriptor.
 *                       cp->scans will be allocated via malloc if
 *                       (cp->scans==NULL or cp->nscan==0).
 *                       Otherwise cp->scans will be realloc'd to
 *                       the new array size. cp->nscans will be set with
 *                       the number of scans initialized.
 *                       On error cp will be left unchanged.
 *  return        int    The number of scans initialized - this may be
 *                       0 if a memory-allocation failure occurs.
 */
static int set_Scans(Corpar *cp)
{
  Subarray *sub;/* The descriptor of the current sub-array */
  Scans *scans; /* Array of 'nscan' Scans */
  int nscan;    /* Number of scans required. */
  int scan;     /* Scan number */
  int uta,utb;  /* Integration indexes of start and end of a scan */
/*
 * Check arguments.
 */
  if(cp==NULL) {
    lprintf(stderr, "set_Scans: NULL plot descriptor intercepted.\n");
    return cp->nscan;
  };
/*
 * Get a convenient pointer to the current sub-array.
 */
  sub = cp->sub;
/*
 * Determine the number of scans required.
 */
  if(!cp->doscan)
    nscan=1;
  else
    nscan=nscans(sub, sub->scangap);
  if(nscan==0)
    return cp->nscan;
/*
 * Allocate memory or realloc?
 */
  if(cp->nscan==0 || cp->scans==NULL)
    scans = (Scans *) malloc(nscan * sizeof(Scans));
  else if(cp->nscan != nscan)
    scans = (Scans *) realloc(cp->scans, nscan * sizeof(Scans));
  else
    scans = cp->scans;
  if(scans==NULL) {
    lprintf(stderr, "set_Scans: Insufficient memory for new Scans\n");
    return cp->nscan;
  };
/*
 * Copy into plot descriptor.
 */
  cp->scans = scans;
  cp->nscan = nscan;
/*
 * Assign scan UT limits in the scan descriptor array.
 */
  if(cp->doscan) {
    for(uta=scan=0; scan<nscan; scan++) {
      scans[scan].sutmin = (sub->integ[uta].ut - cp->utref);
      utb = endscan(sub, sub->scangap, uta);
      scans[scan].sutmax = (sub->integ[utb].ut - cp->utref);
      uta = utb+1;
    };
  } else {
    scans[0].sutmin = sub->integ[0].ut - cp->utref;
    scans[0].sutmax = sub->integ[sub->ntime-1].ut - cp->utref;
  };
  return cp->nscan;
}

/*.......................................................................
 * Determine the amplitude plot range for the ut range and plot options
 * in a passed Corpar descriptor. Assign the new range in the Corpar
 * descriptor.
 *
 * Input:
 *  cp        Corpar * This contains existing plotting attributes.
 *                     All the doxxx flags, the showall flag and uta and
 *                     utb must be initialized before calling this function.
 * Output:
 *  return       int   0 - OK.
 *                     On error -1 is returned and no changes are made
 *                     to cp.
 */
static int c_arange(Corpar *cp)
{
  Integration *integ; /* The descriptor of integration 'ut' */
  int ut;         /* The index of an integration */
  float amin;     /* Minimum amplitude-gain in plot */
  float amax;     /* Maximum amplitude-gain in plot */
  float adif;     /* Amplitude range before adding margins */
  float amp;      /* An amplitude-gain */
/*
 * Check inputs.
 */
  if(cp==NULL) {
    lprintf(stderr, "c_arange: NULL plot intercepted\n");
    return -1;
  };
/*
 * Find the max amplitude gain in the plot range.
 */
  amin = amax = 0.0f;                  /* Always plot zero amplitude-gain */
  for(ut=cp->uta, integ = &cp->sub->integ[ut]; ut<=cp->utb; ut++,integ++) {
    amp = fabs(integ->icor[cp->cif].tcor[cp->ts.ta].amp_cor);
    if(amax < amp) amax = amp;
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
/*
 * Set up return values.
 */
  cp->ampmin = amin;
  cp->ampmax = amax;
/*
 * Return the count of the total number of points in the plot range.
 */
  return 0;
}

/*.......................................................................
 * Return the UT plot range for the ut range and plot options
 * in a passed Corpar descriptor.
 *
 * Input/Output:
 *  cp        Corpar * On entry this contains existing plotting attributes.
 *                     Currently only uta and utb need be initialized.
 *                     On output cp->utmin and cp->utmax will contain
 *                     the min and max UTs of the plot in seconds
 *                     wrt cp->utref, including margins.
 * Output:
 *  return      int    0 - OK.
 *                     On error -1 is returned and no changes are made
 *                     to *utmin or *utmax.
 */
static int c_utrange(Corpar *cp)
{
  float xa;   /* Start UT of range */
  float xb;   /* End UT of range */
  int scan;   /* The index of the scan being processed */
  Scans *sc;  /* The scan descriptor being processed */
/*
 * Check inputs.
 */
  if(cp==NULL) {
    lprintf(stderr, "c_utrange: NULL plot intercepted\n");
    return -1;
  };
/*
 * Valid uta and utb?
 */
  if(cp->uta < 0 || cp->uta>cp->utb || cp->utb >= cp->sub->ntime) {
    lprintf(stderr, "c_utrange: uta and utb are invalid\n");
    return -1;
  };
/*
 * Determine the times corresponding to integrations uta and utb
 * with respect to the reference time vlb->ut, in second units.
 */
  cp->utmin = cp->sub->integ[cp->uta].ut - cp->utref;
  cp->utmax = cp->sub->integ[cp->utb].ut - cp->utref;
/*
 * Determine the displayed UT ranges within the scans.
 * sc->view flags whether any of the scan is visible.
 */
  sc = &cp->scans[0];
  for(scan=0; scan<cp->nscan; scan++,sc++) {
    sc->view = cp->utmax >= sc->sutmin && cp->utmin <= sc->sutmax;
    if(sc->view) {
      xa = (cp->utmin<sc->sutmin)?sc->sutmin:cp->utmin;
      xb = (cp->utmax>sc->sutmax)?sc->sutmax:cp->utmax;
/*
 * Leave a fractional margin around UT range. (Also ensure that the min
 * range is 30 seconds to avoid precision problems).
 */
      if(fabs(xb-xa) > 30.0f) {
	sc->utmin = xa - (xb-xa)*xmarg;
	sc->utmax = xb + (xb-xa)*xmarg;
      } else {
	sc->utmin = xa - 15.0f;
	sc->utmax = xb + 15.0f;
      };
    } else {
      sc->utmin = sc->utmax = 0.0f;
    };
  };
  return 0;
}

/*.......................................................................
 * Set up the viewport limits for the plots leaving 4 char heights on
 * each side of plot for labelling.
 *
 * Input/Output:
 *  cp    Corpar *  The correction parameter struct.
 *                  On output the vxa,vxb,vya,vyb fields will be
 *                  initialized. All other fields are ignored.
 * Output:
 *  return   int     0 - OK.
 */
static int c_vpwin(Corpar *cp)
{
  const float phsfrc=0.5f; /* The fraction of sub-plot to assign to phase */
  float vxa,vxb;  /* X viewport limits enclosing whole plot */
  float vya,vyb;  /* Y viewport limits enclosing whole plot */
  float utsum;    /* Sum of scan UTs within current UT range */
  Scans *sc;      /* A scan descriptor from cp->scans[] */
  int scan;       /* Number of scan */
/*
 * Get the standard viewport.
 */
  cpgsch(1.0f);
  cpgvstd();
  cpgqvp(0, &vxa, &vxb, &vya, &vyb);
/*
 * Store this in the plot descriptor, for subsequent use during labelling.
 */
  cp->vxa = vxa;
  cp->vxb = vxb;
  cp->vya = vya;
  cp->vyb = vyb;
/*
 * cp->vymid specifies the bottom of the amplitude plot and the top of
 * the phase plot.
 */
  cp->vymid = cp->vya + phsfrc * (cp->vyb - cp->vya);
/*
 * Apportion viewports horizontally for different scans.
 * First find the sum of UT ranges covered by all scans within the
 * current UT range.
 */
  sc = &cp->scans[0];
  utsum = 0.0f;
  for(scan=0; scan<cp->nscan; scan++,sc++)
    utsum += sc->utmax - sc->utmin;
  cp->utsum=utsum;
/*
 * Use the fraction of the sum of ut ranges taken up by each scan
 * to determine the fraction of the horizontal viewport range taken
 * up by that scan.
 */
  sc = &cp->scans[0];
  vxa=cp->vxa;
  for(scan=0; scan<cp->nscan; scan++,sc++) {
    sc->vxa=vxa;
    if(sc->view)
      sc->vxb = vxa + (cp->vxb-cp->vxa) * (sc->utmax-sc->utmin)/utsum;
    else
      sc->vxb = sc->vxa; /* Scan not visible */
    vxa = sc->vxb;
  };
  return 0;
}

/*.......................................................................
 * Draw plot axes.
 *
 * Input:
 *  cp      Corpar *  The plot descriptor.
 *  erase      int    If true erase current axes instead of plotting.
 * Output:
 *  return     int    0 - OK.
 *                    Anything else if an error occured.
 */
static int c_plaxes(Corpar *cp, int erase)
{
  const float labsep=2.0f; /* Character separation of labels from axes */
  Scans *sc;     /* The scan being labelled */
  float utmin;   /* Start UT + 1 day, in seconds since start of year */
  float utmax;   /* End UT + 1 day, in seconds since start of year */
  float timinc;  /* Time increment in seconds for labelling */
  float ch;      /* Character height to use */
  int oldcol;    /* Color index on entry to function */
  int scan;      /* The scan index being processed */
  int first;     /* Index of first visible scan */
  int last;      /* Index of final visible scan */
/*
 * Check arguments.
 */
  if(cp==NULL) {
    lprintf(stderr, "c_plaxes: NULL plot descriptor intercepted\n");
    return -1;
  };
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
 * Determine the tick increment to use on the time axis.
 */
  timinc = c_labinc(cp->utsum, 12);
/*
 * Find first and last visible scans.
 */
  sc = &cp->scans[0];
  last = first = -1;
  for(scan=0; scan<cp->nscan; scan++,sc++) {
    if(sc->view) {
      if(first < 0)
	first = scan;
      last = scan;
    };
  };
  if(last < 0) {
    lprintf(stderr, "c_plaxes: No scans visible - can't plot axes\n");
    return -1;
  };
/*
 * Plot the two Y-axes at each end of the frame enclosing the scans.
 */
  cpgsch(ch);
/*
 * Amplitude axis.
 */
  cpgsvp(cp->vxa, cp->vxb, cp->vymid, cp->vyb);
  cpgswin(0.0f, 1.0f, cp->ampmin, cp->ampmax);
  cpgbox(" ", 0.0f, 0, "BCNST", 0.0f, 0);
  cpgmtxt("L", labsep, 0.5f, 0.5f, "Gain");
/*
 * Phase axis.
 */
  cpgsvp(cp->vxa, cp->vxb, cp->vya, cp->vymid);
  cpgswin(0.0f, 1.0f, -180.0f, 180.0f);
  cpgbox(" ", 0.0f, 0, "BCNST", 0.0f, 0);
  cpgmtxt("L", labsep, 0.5f, 0.5f, "Phase  (Degrees)");
/*
 * Do internal and X-axes for each visible scan.
 */
  sc = &cp->scans[first];
  for(scan=first; scan<=last; scan++,sc++) {
/*
 * Calculate the start and end UT in seconds. Add one day such that
 * days in the year start from 1 rather than 0.
 */
    utmin = daysec + cp->utref + sc->utmin;
    utmax = daysec + cp->utref + sc->utmax;
/*
 * Draw internal Y-axes as unadorned vertical lines.
 */
    cpgsvp(cp->vxa, cp->vxb, cp->vya, cp->vyb);
    cpgswin(cp->vxa, cp->vxb, cp->vya, cp->vyb);    
    if(scan != first) {
      cpgmove(sc->vxa, cp->vya);
      cpgdraw(sc->vxa, cp->vyb);
    };
    if(scan != last) {
      cpgmove(sc->vxb, cp->vya);
      cpgdraw(sc->vxb, cp->vyb);
    };
/*
 * Draw the X-axes of the amplitude plot.
 */
    cpgsvp(sc->vxa, sc->vxb, cp->vymid, cp->vyb);
    cpgswin(utmin, utmax, 0.0f, 1.0f);
    cpgtbox("ZHBCST",  /*timinc*/0.0f, 0, " ", 0.0f, 0);
/*
 * Draw the X-axes of the phase plot and enumerate.
 */
    cpgsvp(sc->vxa, sc->vxb, cp->vya, cp->vymid);
    cpgswin(utmin, utmax, 0.0f, 1.0f);
    cpgtbox("ZHBCNST", /*timinc*/0.0f, 0, " ", 0.0f, 0);
  };
/*
 * Restore entry color.
 */
  cpgsci(oldcol);
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Plot or erase amplitude-gain and phase-gain points.
 *
 * Input:
 *  cp       Corpar * The plot descriptor.
 *  uta         int   Index of first integration to be plotted.
 *  utb         int   Index of second integration to be plotted.
 *  erase       int   If true erase points instead of plotting them.
 * Output:
 *  return      int   0 - OK.
 */
static int c_pldata(Corpar *cp, int uta, int utb, int erase)
{
  Subarray *sub;        /* Local pointer to the sub-array being displayed */
  Integration *integ;   /* The descriptor of the integration being displayed */
  Scans *sc;            /* The scan being plotted */
  static int oldcol;    /* Color index on entry to function */
  static int isym;      /* PGPLOT marker symbol to be plotted */
  static int ut;        /* Index of integration being plotted */
  static int first;     /* True until the first point has been plotted */
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
  sub = cp->sub;
/*
 * Draw each point with the appropriate symbol and color for its
 * correction status.
 *
 * Start with the amplitude points.
 */
  sc = &cp->scans[0];
  first=1;
  for(ut=uta,integ = &sub->integ[ut]; ut<= utb; ut++,integ++) {
    Telcor *tcor = &integ->icor[cp->cif].tcor[cp->ts.ta];
    float utval = integ->ut - cp->utref;
    float amp = tcor->amp_cor;
/*
 * Skip to the right scan for this point.
 */
    if(first || utval > sc->sutmax) {
      first = 0;
      while(utval > sc->sutmax)
	sc++;
      cpgsvp(sc->vxa, sc->vxb, cp->vymid, cp->vyb);
      cpgswin(sc->utmin, sc->utmax, cp->ampmin, cp->ampmax);
    };
/*
 * Select the color and marker symbol appropriate for the correction status
 * of the point.
 */
    if(tcor->bad) {
      isym = badsym;
      cpgsci(erase ? 0:badcol);
    } else {
      isym = corsym;
      cpgsci(erase ? 0:corcol);
    };
/*
 * Plot the point.
 */
    cpgpt(1, &utval, &amp, isym);
  };    /* End of  loop over integrations */
/*
 * Now plot the phases.
 */
  sc = &cp->scans[0];
  first = 1;
  for(ut=uta, integ = &sub->integ[ut]; ut<= utb; ut++,integ++) {
    Telcor *tcor = &integ->icor[cp->cif].tcor[cp->ts.ta];
    float utval = integ->ut - cp->utref;
    float phs = tcor->phs_cor;
/*
 * Skip to the right scan for this point.
 */
    if(first || utval > sc->sutmax) {
      first = 0;
      while(utval > sc->sutmax)
	sc++;
      cpgsvp(sc->vxa, sc->vxb, cp->vya, cp->vymid);
      cpgswin(sc->utmin, sc->utmax, -pi, pi);
    };
/*
 * Select the color and marker symbol appropriate for the correction
 * status of the point.
 */
    if(tcor->bad) {
      isym = badsym;
      cpgsci(erase ? 0:badcol);
    } else {
      isym = corsym;
      cpgsci(erase ? 0:corcol);
    };
/*
 * Wrap the phase correction into the range -pi to pi.
 */
    phs -= twopi*floor(phs/twopi+0.5);  /* Wrap into -pi to pi range */
/*
 * Plot the point.
 */
    cpgpt(1, &utval, &phs, isym);
  };                             /* End of loop over integrations */
/*
 * Restore entry color and terminate pgplot buffering.
 */
  cpgsci(oldcol);
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Read the cursor, and record its position etc.. in cp->cursor.
 *
 * Input:
 *  cp     Corpar *  The plot descriptor.
 *  noout     int    If true then don't return until the cursor is
 *                   pressed inside a sub-plot.
 *  mode Bandmode    The desired type of cursor, from:
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
 *  xref    float    xref and yref denote the reference position for
 *  yref    float    band-type cursors. They are ignored if mode==B_NORM.
 *  ci        int    The color index to draw band cursors with.
 * Output:
 *  return    int    0 - OK.
 */
static int c_cursor(Corpar *cp, int noout, Bandmode mode, int isamp,
		    float xref, float yref, int ci)
{
  static float xpos; /* The X world coordinate of the cursor */
  static float ypos; /* The Y world coordinate of the cursor */
  Ccurs *cc;         /* Pointer to cp->cursor cursor descriptor */
  int scan;          /* The scan number being processed */
/*
 * Get the cursor descriptor.
 */
  cc = &cp->cursor;  
/*
 * Set the viewport around the whole viewsurface and make the world
 * coords the same as NDC so that the returned cursor position
 * is measured in NDC.
 */
  cpgsvp(0.0f,1.0f,0.0f,1.0f);
  cpgswin(0.0f,1.0f,0.0f,1.0f);
/*
 * If this is the first call this plot session initialize the position
 * at which to bring up the cursor. Otherwise use the values retained from
 * the previous call in xpos and ypos.
 */
  if(cc->key == KEY_NONE) {
    xpos = 0.5f;
    ypos = 0.5f;
  };
/*
 * Substitute the cross-hair cursor for the normal cursor is requested.
 */
  if(mode == B_NORM && cp->docross)
    mode = B_CROSS;
/*
 * Initialize the return value.
 */
  cc->key = KEY_NONE;
  cc->waslow = 0;
  cc->wasamp = 0;
  cc->utval = 0.0f;
  cc->value = 0.0f;
  cc->sc = NULL;
/*
 * Convert the cursor reference positions into NDC.
 */
  switch(mode) {
  case B_RECT: case B_XRNG: case B_YRNG:
    {
/*
 * Locate the scan that contains the reference UT.
 */
      Scans *sc = cp->scans;
      for(scan=0; scan<cp->nscan; scan++,sc++) {
	if(xref >= sc->utmin && xref <= sc->utmax)
	  break;
      };
      if(scan >= cp->nscan)
	sc = xref<cp->scans[0].utmin ? &cp->scans[0] : &cp->scans[cp->nscan-1];
/*
 * Convert the reference UT and amp/phase to the equivalent NDC position.
 */
      xref = sc->vxa + (xref - sc->utmin) * (sc->vxb - sc->vxa) /
	(sc->utmax - sc->utmin);
/*
 * Get the Y-axis reference value.
 */
      yref = isamp ?
	cp->vymid+(yref-cp->ampmin)*(cp->vyb-cp->vymid)/(cp->ampmax-cp->ampmin):
        cp->vya   + (yref - -pi) * (cp->vymid-cp->vya)/twopi;
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
    cc->waslow = islower((int)key);
    cc->key = cc->waslow ? toupper((int)key) : key;
/*
 * See if the point is in the plot.
 */
    if((xpos >= cp->vxa && xpos <= cp->vxb) && 
       (ypos >= cp->vya && ypos <= cp->vyb)) {
/*
 * In amp part of plot?
 */
      cc->wasamp = ypos > cp->vymid;
/*
 * Convert from NDC to the respective amplitude or phase selected.
 */
      if(cc->wasamp) {
	cc->value = cp->ampmin + 
	  (ypos - cp->vymid)/(cp->vyb - cp->vymid) * (cp->ampmax - cp->ampmin);
      } else {
	cc->value = -pi + (ypos - cp->vya)/(cp->vymid - cp->vya) * twopi;
      };
/*
 * Identify the scan that the cursor was in and use this to
 * determine the selected UT value.
 */
      for(scan=0; scan<cp->nscan; scan++) {
	Scans *sc = &cp->scans[scan];
	if(xpos >= sc->vxa && xpos <= sc->vxb) {
	  cc->utval = sc->utmin + (xpos-sc->vxa)/(sc->vxb-sc->vxa) *
	    (sc->utmax - sc->utmin);
	  cc->sc = sc;
	  break;
	};
      };
/*
 * Cursor not in plot.
 */
    } else {
      cc->utval = cc->value = 0.0f;
      cc->wasamp = 0;
    };
    if(!cc->sc && noout)
      printf("The cursor must be in one of the plots.\n");
  } while(!cc->sc && noout); /* While not in a plot and noout is true */
  return 0;
}

/*.......................................................................
 * Write labels around the frame enclosing all scan sub-plots.
 *
 * Input:
 *  cp       Corpar *  The plot descriptor.
 * Output:
 *  return      int    0 - OK.
 */
static int c_label(Corpar *cp)
{
  Observation *ob;      /* The descriptor of the observation being plotted */
  Subarray *sub;        /* The descriptor of the sub-array being plotted */
  char awrk[81];        /* Work string for labelling */
  char bwrk[81];        /* Work string for labelling */
/*
 * Check arguments.
 */
  if(cp==NULL) {
    lprintf(stderr, "c_label: NULL plot intercepted\n");
    return -1;
  };
/*
 * Get the descriptors of the observation and sub-array being plotted.
 */
  ob = cp->ob;
  sub = cp->sub;
/*
 * Set the viewport around the plot grid.
 */
  cpgsvp(cp->vxa, cp->vxb, cp->vya, cp->vyb);
/*
 * Compose and write main titles.
 */
  cpgsci(1);
  cpgsch(1.0f);
  sprintf(awrk, "%s  %s", ob->source.name,
          sutdate(ob->date.year, ob->date.ut, bwrk));
  cpgmtxt("T", 1.7f, 0.0f, 0.0f, awrk);
  sprintf(awrk, "Corrections for IF %d  Pol %s  Station %d:%.20s",
           cp->cif+1, Stokes_name(ob->stream.pol.type), cp->ts.isub+1,
	  sub->tel[cp->ts.ta].name);
  cpgmtxt("T", 0.5f, 0.0f, 0.0f, awrk);
  sprintf(awrk, "%d of %d", cp->ts.ta+1, sub->nstat);
  cpgmtxt("T", 0.5f, 1.0f, 1.0f, awrk);
/*
 * In non-interative mode, tell the user what is being plotted.
 */
  if(!cp->docurs)
    lprintf(stdout, "Page %02.2d: Station %d:%s\n", cp->npage,
	    cp->ts.isub+1, sub->tel[cp->ts.ta].name);
/*
 * Write X-axis label.
 */
  cpgmtxt("B",3.0f,0.5f,0.5f,"Correction UT");
  return 0;
}

/*.......................................................................
 * Replot the current scans.
 *
 * Input:
 *  cp        Corpar *  The plot descriptor.
 * Output:
 *  return       int    0 - OK.
 */
static int c_redisp(Corpar *cp)
{
  int ierr=0; /* True if an error occurs */
/*
 * Cursory check of arguments.
 */
  if(cp==NULL) {
    lprintf(stderr, "c_redisp: NULL plot intercepted\n");
    return -1;
  };
/*
 * Clear page.
 */
  cpgpage();
/*
 * Count pages.
 */
  cp->npage++;
/*
 * Find scan limits.
 */
  set_Scans(cp);
/*
 * Determine the UT plot range for all plots.
 */
  ierr = ierr || c_utrange(cp);
/*
 * Set up viewport.
 */
  ierr = ierr || c_vpwin(cp);
/*
 * Draw scans.
 */
  cpgbbuf();
  ierr = ierr || c_arange(cp);
  ierr = ierr || c_plaxes(cp, 0);
  ierr = ierr || c_pldata(cp, cp->uta, cp->utb, 0);
  ierr = ierr || c_label(cp);
  cpgebuf();
  return ierr;
}

/*.......................................................................
 * Determine scaling factors required to convert from world coordinates
 * to mm.
 *
 * Input:
 *  cp   Corpar *  The plot descriptor.
 *  doamp   int    If true make ytomm the conversion from amplitude
 *                 to mm. Otherwise yield the phase conversion.
 * Output:
 *  xtomm float *  ut_shift * *xtomm yields the physical size
 *                 of the ut shift in mm.
 *  ytomm float *  Amplitude_or_phase_shift * *xtomm yields the physical
 *                 size of the amplitude or phase shift in mm.
 *  return  int    0 - OK.
 */
static int c_scale(Corpar *cp, int doamp, float *xtomm, float *ytomm)
{
  float xa,xb,ya,yb; /* Physical coordinates of viewport */
  int scan;   /* The number of a scan */
  Scans *sc;  /* The first displayed scan */
/*
 * Find the first displayed scan.
 */
  sc = &cp->scans[0];
  for(scan=0; scan<cp->nscan && !sc->view; scan++, sc++);
  if(scan >= cp->nscan) {
    lprintf(stderr, "c_scale: No scans visible\n");
    return -1;
  };
/*
 * Determine the size of the viewport in physical device coordinates
 * (millimeters).
 */
  if(doamp)
    cpgsvp(sc->vxa, sc->vxb, cp->vymid, cp->vyb);
  else
    cpgsvp(sc->vxa, sc->vxb, cp->vya, cp->vymid);
  cpgqvp(2,&xa,&xb,&ya,&yb);
/*
 * Calculate factors to convert world coords into mm.
 */
  *xtomm = fabs((xb-xa)/(sc->utmax-sc->utmin));
  if(doamp)
    *ytomm = fabs((yb-ya)/(cp->ampmax-cp->ampmin));
  else
    *ytomm = fabs((yb-ya)/twopi);
  return 0;
}

/*.......................................................................
 * Split a time range into suitable increments for labelling with pgtbox().
 *
 * Input:
 *  tspan  float   The time span to cover, in seconds.
 *  ntry     int   The number of tick intervals to try to split tspan
 *                 into.
 * Output:
 *  return float   The tick interval chosen (seconds). On error this is
 *                 0.0f.
 */
static float c_labinc(float tspan, int ntry)
{
  /* Rounded fractions of seconds allowed (after scaling by 10) */
  static const int fracinc[]={1,2,5,10};
  static const int nfrac=sizeof(fracinc)/sizeof(int);

  /* Numbers divisible into 60 */
  static const int sixtyinc[]={ 1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60};
  static const int nsixty=sizeof(sixtyinc)/sizeof(int);

  /* Numbers divisible into 24 */
  static const int hourinc[]={ 1, 2, 3, 4, 6, 8, 12, 24};
  static const int nhour=sizeof(hourinc)/sizeof(int);

  static const float hrsec=3600.0f;
  static const float minsec=60.0f;
  static const float frcsec=10.0f;

  static const float minspan=0.1f;  /* Min time increment */
  float newinc; /* The final estimate of the tic increment */
  float tinc;   /* The initial estimate of the increment from the input */
  int itinc;    /* Input increment scaled into other units */
  int i;
/*
 * Check arguments.
 */
  if(ntry < 2) {
    lprintf(stderr, "labinc: Too few (%d) ticks requested\n", ntry);
    return 0.0f;
  };
/*
 * Get an initial estimate of the tick increment.
 */
  tinc = tspan/ntry;
/*
 * Find the rounded tick increment that gives the same or a higher
 * number of ticks than found in the input approximation.
 */
  if(tinc < minspan) {
    lprintf(stderr, "labinc: Time span too small for labelling\n");
    newinc=0.0f;
  }
/*
 * Less than a second?
 */
  else if(tinc < 1.0f) {
    itinc = fnint(tinc*frcsec);
    for(i=1; i<nfrac && fracinc[i]<=itinc; i++);
    newinc=fnint(fracinc[i]/frcsec);
  }
/*
 * Less than a minute?
 */
  else if(tinc < minsec) {
    itinc = fnint(tinc);
    for(i=1; i<nsixty && sixtyinc[i]<=itinc; i++);
    newinc=sixtyinc[i];
  }
/*
 * Less than an hour?
 */
  else if(tinc < hrsec) {
    itinc = fnint(tinc/minsec);
    for(i=1; i<nsixty && sixtyinc[i]<=itinc; i++);
    newinc=fnint(sixtyinc[i]*minsec);
  }
/*
 * Less than a day?
 */
  else if(tinc < daysec) {
    itinc = fnint(tinc/hrsec);
    for(i=1; i<nhour && hourinc[i]<=itinc; i++);
    newinc=fnint(hourinc[i]*hrsec);
  }
/*
 * Less than a year?
 */
  else if(tinc < 3.1536e+07) {
    newinc = fnint(daysec*fnint(tinc/daysec));
  }
  else {
    lprintf(stderr, "labinc: Time range too large for algorithm\n");
    newinc=0.0f;
  };
/*
 * Return the increment.
 */
  return newinc;
}

/*.......................................................................
 * Provide an interactive means of perusing and editing telescope
 * corrections using the cursor.
 *
 * Input:
 *  ob  Observation *  The observation who's corrections are to be
 *                     displayed.
 *  ts      Telspec *  The specification of the first telescope to plot,
 *                     or NULL (or an empty specification) for the default.
 *  cif         int    The index of the first IF to plot.
 *  docurs      int    If true enter interactive cursor mode.
 * Input/Output:
 *  modified    int *  If modified!=NULL then *modified will be assigned
 *                     the value 1 if the data were edited, or 0 otherwise.
 * Output:
 *  return      int    0 - ok.
 *                     1 - error.
 */
int corplot(Observation *ob, Telspec *ts, int cif, int docurs, int *modified)
{
  int ierr=0;      /* Error flag */
  int wasflag=0;   /* True if last key stroke toggled a flag */
  int nflag=0;     /* Number of flag toggling operations done */
  int ut;          /* A integration index */
  Corpar *cp;      /* Plot descriptor */
/*
 * Data not edited yet.
 */
  if(modified)
    *modified = 0;
/*
 * Check the state of the observation.
 */
  if(!ob_ready(ob, OB_INDEX, "corplot"))
    return 1;
/*
 * Allocate a plot descriptor.
 */
  cp = new_Corpar(ob, ts, cif, docurs, 1);
  if(cp==NULL)
    return 1;
/*
 * Interactive plotting?
 */
  if(cp->docurs) {
/*
 * Inform user of the help option.
 */
    lprintf(stdout,
	    "Move the cursor into the plot window and press \'%c\' for help\n",
	    KEY_HELP);
/*
 * Display the corrections for the telescope specified.
 */
    ierr = ierr || c_redisp(cp);
/*
 * Start interactive editing loop.
 */
    cp->cursor.key = KEY_NONE;
    while(!ierr && cp->cursor.key!=KEY_QUIT) {
/*
 * Read the cursor.
 */
      do {
	ierr = c_cursor(cp, 0, B_NORM, 0, 0.0f, 0.0f, 1);
/*
 * Toggle flags where appropriate.
 */
	wasflag = c_flags(cp, cp->cursor.key, cp->cursor.waslow) == 0;
	nflag += wasflag;
      } while(wasflag);   /* Don't do anything more if a flag was toggled */
/*
 * Take action appropriate to the key that the user pressed.
 */
      if(nflag > 0) {  /* Update display after a sequence of flag toggling */
	nflag = 0;
	ierr = c_redisp(cp);
      } else {
/*
 * Obey the cursor request.
 */
	switch(cp->cursor.key) {
	case KEY_NEXT:                    /* Plot next page */
	  ierr = c_newtel(cp, cp->cursor.waslow ? C_NXT_TA:C_NXT_ISUB, 1, 1,
			  NULL) < 0;
	  break;
	case KEY_PREV:                    /* Plot previous page */
	  ierr = c_newtel(cp, cp->cursor.waslow ? C_NXT_TA:C_NXT_ISUB, 0, 1,
			  NULL) < 0;
	  break;
	case KEY_DIS:
	  ierr = ierr || c_redisp(cp);
	  break;
	case KEY_UT:
	  ierr = ierr || c_newut(cp);
	  break;
	case KEY_PRVIF:
	case KEY_NXTIF:
	  {
	    int step = cp->cursor.key==KEY_NXTIF ? 1 : -1;
	    int cif = nextIF(ob, cp->cif + step, 0, step);
	    if(cif >= 0) {
	      cp->cif = cif;
	      ierr = c_redisp(cp);
	    };
	  };
	  break;
	case KEY_TEL:
/*
 * Identify the telescope and its sub-array and plot the selection.
 */
	  {
	    Telspec *ts = read_Telspec(cp->ob, NULL, NULL, cp->ts.isub);
	    c_newtel(cp, C_ALLNEW, 1, 1, ts);
	  };
	  break;
	case KEY_CUR:
	  if(cp->cursor.sc) {
	    ut = c_find(cp, cp->cursor.utval, cp->cursor.value,
			cp->cursor.wasamp);
	    ierr = ierr || c_edit_cor(cp, ut, ED_FLAG);
	  };
	  break;
	case KEY_CAN:
	  if(cp->cursor.sc) {
	    ut = c_find(cp, cp->cursor.utval, cp->cursor.value,
			cp->cursor.wasamp);
	    ierr = ierr || c_edit_cor(cp, ut, ED_RESET);
	  };
	  break;
	case KEY_CROSS:   /* Toggle cross-hair cursor mode */
	  cp->docross = !cp->docross;
	  break;
	case KEY_HELP:
	  printf("List of keys to enter via cursor.\n");
	  printf(" %c - Quit this session.\n", KEY_QUIT);
	  printf(" %c - Display corrections of the next telescope.\n",
		 tolower(KEY_NEXT));
	  printf(" %c - Display corrections of the previous telescope.\n",
		 tolower(KEY_PREV));
	  printf(" %c - Display corrections of the Next sub-array.\n",
		 KEY_NEXT);
	  printf(" %c - Display corrections of the Previous sub-array.\n",
		 KEY_PREV);
	  printf(" %c - Display corrections of the Next IF.\n", KEY_NXTIF);
	  printf(" %c - Display corrections of the Previous IF.\n", KEY_PRVIF);
	  printf(" %c - Select new UT range with cursor key %c.\n", KEY_UT,
		 KEY_CUR);
	  printf(" %c - Redisplay current plot.\n", KEY_DIS);
	  printf(" %c - Select displayed telescope from keyboard.\n", KEY_TEL);
	  printf(" %c - Toggle the correction flag of the nearest point.\n",
		 KEY_CUR);
	  printf(" %c - Uncorrect the telescope correction of the nearest point.\n", KEY_CAN);
	  printf(" %c - Toggle breaking of display into scans.\n", KEY_BRK);
	  printf(" %c - Toggle whether to use a cross-hair cursor if available.\n", KEY_CROSS);
	  break;
	};
      };
    };
  }
/*
 * Non-interactive plotting?
 */
  else if(!ierr) {
    int iret=0;
/*
 * Plot the first page.
 */
    ierr = c_redisp(cp);
/*
 * If requested, plot the rest of the available pages.
 */
    while((iret = c_newtel(cp, C_NXT_TEL, 1, 0, NULL)) == 0)
      ;
    ierr = iret < 0;
  };
/*
 * Flush any pending edits.
 */
  ed_flush(ob);
/*
 * Data modified?
 */
  if(modified!=NULL)
    *modified = cp->modified;
/*
 * Delete the plot descriptor.
 */
  cp = del_Corpar(cp);
  return ierr;
}

/*.......................................................................
 * Receive input of new UT range via the cursor and redisplay the plot
 * within that range. If the user presses the KEY_UT key then display
 * the full UT range.
 *
 * Input:
 *  cp      Corpar *  The plot descriptor.
 * Output:
 *  return     int    0 - OK. Anything else means fatal error.
 */
static int c_newut(Corpar *cp)
{
  int accepted;   /* True when cursor entry accepted */
  int dofull=0;   /* True if full UT range requested */
  int iter;       /* Iterate over getting two valid keypresses */
  float utval[2]={0.0f,0.0f}; /* The two UT end points wrt cp->refut. */
/*
 * Get the first cursor position for the new UT range.
 */
  for(iter=0; iter<2 && !dofull; iter++) {
    do {
      accepted = 0;
      if(c_cursor(cp, 0,iter==0?B_XVAL:B_XRNG, 0, utval[0], 0.0f, zoomcol))
	return 1;
      switch(cp->cursor.key) {
      case KEY_UT:        /* Display full range and end UT selection */
	accepted = dofull = 1;
	break;
      case KEY_QUIT: case KEY_CAN:     /* Abort UT selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected start UT */
	utval[iter] = cp->cursor.utval;
	accepted=1;
	break;
      default:                  /* Unexpected cursor input key - show usage */
	printf("To select a new UT display range use keys:\n");
	printf(" %c - Select the %s UT.\n", KEY_CUR, iter==0 ? "start":"end");
	printf(" %c - Cancel UT display range selection.\n", KEY_CAN);
	printf(" %c - Display the full UT display range available.\n", KEY_UT);
	break;
      };
    } while(!accepted);
  };
/*
 * Get the UT indexes.
 */
  if(dofull) {
    cp->uta = 0;
    cp->utb = cp->sub->ntime - 1;
  } else {
    int ut;               /* UT index loop variable */
    Integration *integ;   /* Descriptor of integration = sub->integ[ut] */
    double utmin = utval[0] + cp->utref;
    double utmax = utval[1] + cp->utref;
/*
 * Swap the range limits such that utmin < utmax.
 */
    if(utmin>utmax) {double dtmp=utmin; utmin=utmax; utmax=dtmp;};
/*
 * Locate the equivalent ut indexes to the selected ut values.
 */
    for(ut=cp->uta,integ = &cp->sub->integ[ut]; ut<cp->utb && integ->ut < utmin; ut++,integ++);
    cp->uta = ut;
    for(ut=cp->uta; ut<=cp->utb && integ->ut <= utmax; ut++,integ++);
    cp->utb = (cp->uta<ut) ? (ut-1):cp->uta;
  };
/*
 * Display the result.
 */
  return c_redisp(cp);
}

/*.......................................................................
 * Take a cursor position returned by c_cursor() and locate the index
 * of the closest plotted point.
 *
 * Input:
 *  cp       Corpar *   The plot descriptor.
 *  utval     float     The relative UT returned by c_cursor().
 *  yval      float     The amp or phase returned by c_cursor().
 *  isamp       int     The point type returned by c_cursor().
 * Output:
 *  return      int     The integration index of the nearest point.
 */
static int c_find(Corpar *cp, float utval, float yval, int isamp)
{
  double vlbut;       /* The ut selected, in same units as integration UTs */
  Integration *integ; /* Descriptor of an integration */
  Telcor *tcor;       /* Pointer to correction being considered */
  int ut;             /* The integration index being checked */
  float xtomm,ytomm;  /* Conversion factor between world coords and mm */
  int bestut=0;       /* The ut index of the point closest to the cursor */
  float xdif;         /* X Distance in mm between data point and cursor */
  float ydif;         /* Y Distance in mm between data point and cursor */
  float dist;         /* The squared distance (mm) from data point to cursor */
  float mindist=0.0f; /* The min value of 'dist' */
  float phs;          /* A visibility phase wrapped into the range -pi to pi */
  int first=1;        /* True until end of first iteration of search loop */
/*
 * Determine conversion factors from world coords to mm.
 */
  if(c_scale(cp, isamp, &xtomm, &ytomm))
    return -1;
/*
 * Calculate the UT wrt the start of the year.
 */
  vlbut = utval + cp->utref;
/*
 * Locate the nearest point.
 */
  for(ut=cp->uta, integ = &cp->sub->integ[ut]; ut<=cp->utb; ut++,integ++) {
    tcor = &integ->icor[cp->cif].tcor[cp->ts.ta];
    xdif = xtomm * (integ->ut - vlbut);
    if(isamp)
      ydif = ytomm * (yval - fabs(tcor->amp_cor));
    else {
      phs = tcor->phs_cor - twopi * floor(tcor->phs_cor/twopi + 0.5);
      ydif = ytomm * (yval - phs);
    };
/*
 * Compare the squared distance from this point with the that of
 * the previous closest point.
 */
    dist = xdif*xdif + ydif*ydif;
    if(first || dist < mindist) {
      first = 0;
      bestut = ut;
      mindist = dist;
    };
  };
/*
 * No points in plot!
 */
  if(first)
    return -1;
  return bestut;
}

/*.......................................................................
 * Flag or reset the corrections of the current telescope at a given
 * integration.
 *
 * Input:
 *  cp       Corpar *   The plot descriptor.
 *  ut          int     The integration index of the point to be zapped.
 *  mode     Edmode     ED_RESET - Reset the correction.
 *                      ED_FLAG  - Toggle the correction flag.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
static int c_edit_cor(Corpar *cp, int ut, Edmode mode)
{
  Telcor *tcor;  /* The correction being modified */
/*
 * Flag the data as modified.
 */
  cp->modified = 1;
/*
 * Get the correction to be modified.
 */
  tcor = &cp->sub->integ[ut].icor[cp->cif].tcor[cp->ts.ta];
/*
 * Erase the displayed corrections.
 */
  if(c_pldata(cp, ut, ut, 1))
    return 1;
/*
 * Flag or reset the correction.
 */
  switch(mode) {
  case ED_RESET:
    if(clr_Telcor(cp->ob, cp->sub, cp->cif, ut, cp->ts.ta))
      return 1;
    break;
  case ED_FLAG:
    if(ed_Telcor(cp->ob, cp->sub, cp->cif, ut, cp->ts.ta, !tcor->bad))
      return 1;
    break;
  };
/*
 * Plot the reset corrections.
 */
  if(c_pldata(cp, ut, ut, 0))
    return 1;
  return 0;
}

/*.......................................................................
 * Toggle plotting flags given a command key.
 *
 * Input:
 *  cp    Corpar *  The plot descriptor.
 *  key     char    The command key (upper case).
 *  waslow   int    If true the key was lower case.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Key unrecognized.
 */
static int c_flags(Corpar *cp, char key, int waslow)
{
  switch (key) {
  case KEY_BRK:
    cp->doscan = !cp->doscan;
    break;
  default:
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Plot the corrections for a new telescope.
 *
 * Input:
 *  cp         Corpar *   Plot-parameter block.
 *  oper        Telop     C_ALLNEW  - Display the telescope specified in
 *                                    *init.
 *                        C_NXT_ISUB - Start plotting from the next sub-array.
 *                        C_NXT_TA  - Plot the next telescope.
 *                        C_NXT_TEL - Plot the next telescope if allowed
 *                                    by the current telescope specification.
 *  forward       int     0 - Locate the next telescope in order of
 *                            decreasing telescope index.
 *                        1 - Locate the next telescope in order of
 *                            increasing telescope index.
 *  init      Telspec *  The new telescope specification to be used when
 *                       oper == C_ALLNEW.
 * Output:
 *  return        int     0 - OK.
 *                        1 - No more telescopes in specified direction.
 *                       -1 - Error.
 */
static int c_newtel(Corpar *cp, Telop oper, int forward, int report,
		    Telspec *init)
{
  Telspec ts;      /* The new telescope specification */
/*
 * Handle the specified change in reference telescope.
 */
  switch(oper) {
  case C_ALLNEW:
    ts = *init;
    if(next_tel(cp->ob, FIND_FIRST, forward, 0, 0, report, &ts))
      return 1;
    break;
  case C_NXT_ISUB:
    ts = cp->ts;
    if(next_tel(cp->ob, SKIP_SUB, forward, 0, 0, report, &ts))
      return 1;
    break;
  case C_NXT_TA:
    ts = cp->ts;
    if(next_tel(cp->ob, SKIP_TA, forward, 0, 0, 0, &ts) &&
       next_tel(cp->ob, SKIP_SUB, forward, 0, 0, report, &ts))
      return 1;
    break;
  case C_NXT_TEL:
    ts = cp->ts;
    if(next_tel(cp->ob, FIND_NEXT, forward, 0, 0, report, &ts))
      return 1;
    break;
  default:
    lprintf(stderr, "c_newtel: Unrecognised opcode.\n");
    return -1;
  };
/*
 * If the sub-array changed, initialize for the next sub-array.
 */
  if(cp->sub==NULL || cp->ts.isub != ts.isub) {
    cp->sub = cp->ob->sub + ts.isub;
/*
 * Set up scan info for the new sub-array.
 */
    if(set_Scans(cp)==0)
      return -1;
/*
 * Assign starting values to the rest of the members.
 */
    cp->uta = 0;
    cp->utb = cp->sub->ntime-1;  /* Default to show all data */
  };
/*
 * Record the new telescope specification.
 */
  cp->ts = ts;
/*
 * Display the the new telscope's corrections.
 */
  return c_redisp(cp) ? -1 : 0;
}


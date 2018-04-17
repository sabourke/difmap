#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "obs.h"
#include "units.h"
#include "vlbutil.h"
#include "telspec.h"
#include "visplot.h"
#include "cpgplot.h"
#include "logio.h"

static const int datcol=10;  /* PGPLOT color index for observed data */
static const int altcol=1;   /* Alternative to datcol for highlighting data */
static const int axcol=1;    /* PGPLOT color index for axes */
static const int zoomcol=5;  /* PGPLOT color index for zoom cursor window */
static const int cutcol=2;   /* PGPLOT color index for edit cursor window */

/* Define selection keys */

enum {
  KEY_NONE='\0',  /* Null key */
  KEY_DOT ='.',   /* Key to toggle marker symbol size */
  KEY_CUR ='A',   /* Key to enter positions */
  KEY_CUT ='C',   /* Key to initiate cut area selection */
  KEY_CAN ='D',   /* Cancel selection */
  KEY_HELP='H',   /* Key to request help */
  KEY_DIS ='L',   /* Key to request re-display of plot */
  KEY_NXT ='N',   /* Key to highlight next station */
  KEY_PRV ='P',   /* Key to highlight previous station */
  KEY_SHOW='S',   /* Report on the nearest point to the cursor */
  KEY_TEL ='T',   /* Key to select highlighted telescope by name */
  KEY_CH  ='W',   /* Toggle channel editing mode */
  KEY_QUIT='X',   /* Key to quit interactive session */
  KEY_ZOOM='Z',   /* Key to zoom in(out) on a selected area */
  KEY_CROSS='+',  /* Toggle cross-hair cursor mode */
  KEY_CONJ='%'    /* Toggle whether to display the conjugate visibility */
};

/* Type used to return cursor selection */

typedef struct {
  float uu;      /* Cursor selected U coordinate (wavelengths) */
  float vv;      /* Cursor selected V coordinate (wavelengths) */
  int key;       /* The upper-case character corresponding to the key pressed */
  int waslow;    /* The original case of the selection key */
} Keypos;

typedef struct {
  Observation *ob; /* The descriptor of the observation being plotted */
  Keypos kp;       /* The descriptor of the last cursor selection */
  Telspec init;    /* The specification of the first available telescope */
  Telspec ts;      /* Specification of the reference telescope */
  int highlight;   /* True when telescope highlighting is enabled */
  int fixu;        /* If true don't autoscale the U axis range */
  int fixv;        /* If true don't autoscale the V axis range */
  float umin;      /* The plot encloses umin to umax on the U axis */
  float umax;      /* The plot encloses umin to umax on the U axis */
  float vmin;      /* The plot encloses vmin to vmax on the V axis */
  float vmax;      /* The plot encloses vmin to vmax on the V axis */
  int docurs;      /* If true allow cursor interaction */
  int dobig;       /* If true, plot with larger dot size */
  int docross;     /* True to enable cross-hair mode */  
  int doconj;      /* If true, display conjugated visibilities */
  int modified;    /* Remains 0 unless the data are edited */
  int ch_ed;       /* If true, edits are restricted to current freq channels */
} U_par;

static int u_redisp(U_par *up);
static int u_axes(U_par *up, int axcol);
static void u_namplt(U_par *up, int erase);
static int u_uvplot(U_par *up, int othcol, int refcol);
static void u_basepl(U_par *up, Subarray *sub, int base, int color);
static int u_zoom(U_par *up);
static int u_mlab(U_par *up, int erase);
static int u_newmode(U_par *up, int ch_ed);
static int u_edbox(U_par *up);
static int u_editpt(U_par *up, int cif, int isub, int base, int ut, int flag);
static int u_setrange(U_par *up, int fixu, float umin, float umax,
		                 int fixv, float vmin, float vmax);
static int u_getrange(U_par *up);

typedef enum {U_ALLNEW, U_NXTSUB, U_NXTTEL} Telop;

static int u_newtel(U_par *up, Telop oper, int forward, Telspec *init);

/* Define a structure for the return type of u_findpt() */

typedef struct {
  int found;     /* True only if a point was selected */
  int ut;        /* The index of the integration of the closest point */
  int base;      /* The index of the baseline of the closest point */
  int isub;      /* The index of the sub-array of the closest point */
  int cif;       /* The index of the IF of the closest point */
} Bestvis;

static Bestvis u_findpt(U_par *up, float xpos, float ypos);

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

static int u_cursor(U_par *up, Bandmode mode, float xref, float yref, int ci);

/*.......................................................................
 * Plot observed visibility points in a plot of the UV plane.
 * NB. This function uses the currently set viewport to set its window
 * within, so the caller must ensure that this is valid.
 *
 * Input:
 *  ob  Observation *  The observation to be plotted.
 *  ts      Telspec *  The specification of the first telescope to highlight,
 *                     or NULL (or an empty specification) for no highlighting.
 *  docurs      int    If true, provide interaction via cursor.
 *  opts       char *  An optional string of flag toggling keys that
 *                     toggle after the values below have been applied.
 *                     Send NULL or empty string if not required.
 *  umax      float    The maximum U coordinate to be displayed
 *                     (wavelengths).
 *  vmax      float    The maximum U coordinate to be displayed
 *                     (wavelengths).
 * Input/Output:
 *  modified   int *   If modified!=NULL then *modified will be assigned
 *                     the value 1 if the data were edited, or 0 otherwise.
 * Output:
 *  return      int    0 - OK.
 */
int uvplot(Observation *ob, Telspec *ts, int docurs, char *opts,
	   float umax, float vmax, int *modified)
{
  U_par up;            /* Internal parameter container */
  Bestvis best;        /* Container for the results of a visibility search */
  char answer[10];     /* String to hold answers to PGPLOT inquiries */
  int slen;            /* Length of answer in 'answer' */
  int oldcol;          /* Color on entry to this function */
  int old_if;          /* State of current IF to be restored on exit */
  int ierr=0;          /* Error status flag */
  int i;
/*
 * Data not yet modified.
 */
  if(modified)
    *modified = 0;
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_SELECT, "uvplot"))
    return 1;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Record the Observation descriptor.
 */
  up.ob = ob;
/*
 * Set the UV range to be displayed.
 */
  if(u_setrange(&up, 1, -umax, umax, 1, -vmax, vmax))
    return 1;
/*
 * Get the specification of the first available telescope.
 */
  {
    Telspec *init = find_tel(ob, 0, 0, 0, 1, 0, 0, 1);
    if(!init)
      return 1;
    up.init = *init;
    up.init.nfix = 2;
  };
/*
 * Record the initial highlight telescope specification.
 */
  if(ts && ts->nfix>0) {
    if(next_tel(ob, FIND_FIRST, 1, 0, 0, 1, ts))
      return 1;
    up.ts = *ts;
    up.highlight = 1;
  } else {
    up.ts = up.init;
    up.highlight = 0;
  };
/*
 * Flag the data as un-modified.
 */
  up.modified = 0;
  up.ch_ed = 0;
/*
 * Default to small dot markers.
 */
  up.dobig = 0;
  up.docross = 0;
  up.doconj=1;
/*
 * If a string of flag options was given, interpret them here.
 */
  if(opts != NULL) {
    slen = strlen(opts);
    for(i=0; i<slen; i++) {
      int key = opts[i];
/*
 * Handle flags that are normally treated as commands.
 */
      switch(key) {
      case KEY_DOT:
	up.dobig = !up.dobig;
	break;
      case KEY_CROSS:
	up.docross = !up.docross;
	break;
      case KEY_CONJ:
	up.doconj = !up.doconj;
	break;
      };
    };
  };
/*
 * If cursor interaction has been requested, check whether the current
 * device has a cursor.
 */
  if(docurs) {
    slen = sizeof(answer)-1;
    cpgqinf("CURSOR", answer, &slen);
    docurs = strncmp(answer,"YES",3) == 0;
  };
  up.docurs = docurs;
/*
 * Tell the user how to list keys.
 */
  if(docurs) {
    lprintf(stdout,
	    "Move the cursor into the plot window and press \'%c\' for help\n",
	    KEY_HELP);
  };
/*
 * Store the entry color.
 */
  cpgqci(&oldcol);
/*
 * Initial plot.
 */
  ierr = u_redisp(&up);
/*
 * Non-interactive mode?
 */
  if(!docurs) {
    cpgsci(oldcol);
    return ierr;
  };
/*
 * Start the interactive loop.
 */
  up.kp.key = KEY_NONE;
  while(!ierr && up.kp.key != KEY_QUIT) {
/*
 * Read a cursor position and the key that was pushed to return the
 * position.
 */
    ierr = u_cursor(&up, B_NORM, 0.0f, 0.0f, 1);
    if(ierr)
      break;
/*
 * Interpret the key pressed.
 */
    switch(up.kp.key) {
    case KEY_DIS:
      ierr = u_redisp(&up);
      break;
    case KEY_DOT:
      up.dobig = !up.dobig;
      ierr = u_redisp(&up);
      break;      
    case KEY_NXT:
      ierr = u_newtel(&up, up.kp.waslow ? U_NXTTEL:U_NXTSUB, 1, NULL);
      break;
    case KEY_PRV:
      ierr = u_newtel(&up, up.kp.waslow ? U_NXTTEL:U_NXTSUB, 0, NULL);
      break;
    case KEY_TEL:
      {
	Telspec *ts = read_Telspec(up.ob, NULL, NULL, up.ts.isub);
	ierr = ts && u_newtel(&up, U_ALLNEW, 1, ts);
      };
      break;
    case KEY_SHOW:
      best = u_findpt(&up, up.kp.uu, up.kp.vv);
      if(best.found) {
	Subarray *subptr = &up.ob->sub[best.isub];
	Baseline *bptr = &subptr->base[best.base];
	char date_str[20];
	write_ut(subptr->integ[best.ut].ut, sizeof(date_str), date_str);
	printf("Visibility on baseline %d:%s-%s (IF %d) at UT %s\n",
	       best.isub+1, subptr->tel[bptr->tel_a].name,
	       subptr->tel[bptr->tel_b].name, best.cif+1, date_str);
      };
      break;
    case KEY_CUT:  /* Initiate cut area selection */
      ierr = u_edbox(&up);
      break;
    case KEY_CH:       /* Toggle channel editing mode */
      ierr = u_newmode(&up, !up.ch_ed);
      break;
    case KEY_CROSS:   /* Toggle cross-hair cursor mode */
      up.docross = !up.docross;
      break;
    case KEY_CONJ:   /* Toggle conjugate reflection mode */
      up.doconj = !up.doconj;
      ierr = u_redisp(&up);
      break;
    case KEY_HELP:
      printf("You requested help by pressing \'%c\'.\n", KEY_HELP);
      printf("The following keys are defined when pressed inside the plot:\n");
      printf(" %c - Quit uvplot\n", KEY_QUIT);
      printf(" %c - Re-display plot.\n", KEY_DIS);
      printf(" %c - Zoom in on a rectangular sub-plot.\n", KEY_ZOOM);
      printf(" %c - Re-display plot with alternate marker symbol.\n", KEY_DOT);
      printf(" %c - Highlight next telescope\n", tolower(KEY_NXT));
      printf(" %c - Highlight previous telescope\n", tolower(KEY_PRV));
      printf(" %c - Step to the next sub-array to highlight.\n", KEY_NXT);
      printf(" %c - Step to the preceding sub-array to highlight.\n", KEY_PRV);
      printf(" %c - Specify highlighted telescope from keyboard\n", KEY_TEL);
      printf(" %c - Show the baseline and time of the nearest point to the cursor\n", KEY_SHOW);
      printf(" %c - Initiate selection of an area to flag.\n", KEY_CUT);
      printf(" %c - Toggle spectral-line channel based editing.\n",KEY_CH);
      printf(" %c - Toggle whether to use a cross-hair cursor if available.\n",
	     KEY_CROSS);
      printf(" %c - Toggle whether to display conjugate symmetric visibilities.\n",
	     KEY_CONJ);
      break;
    case KEY_ZOOM:
      ierr = u_zoom(&up);
      break;
    };
  };
/*
 * Reinstate entry color.
 */
  cpgsci(oldcol);
/*
 * Data modified.
 */
  if(modified!=NULL)
    *modified = up.modified;
/*
 * Flush any pending edits.
 */
  ed_flush(ob);
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    ierr = 1;
  return ierr;
}

/*.......................................................................
 * Private function of uvradplt() to (re-)plot the axes,data and model
 * points.
 *
 * Input:
 *  up        U_par *  Plot-parameter block.
 * Output:
 *  return      int    0 - OK.
 */
static int u_redisp(U_par *up)
{
  int ierr=0;      /* Error status flag */
/*
 * Start a new page.
 */
  cpgpage();
/*
 * Buffer PGPLOT commands until finished.
 */
  cpgbbuf();
/*
 * Assign the standard viewport.
 */
  cpgsch(1.0f);
  cpgvstd();
/*
 * Plot axes.
 */
  ierr = ierr || u_axes(up, axcol);
/*
 * Plot the mode line if interactive.
 */
  ierr = ierr || (up->docurs && u_mlab(up, 0));
/*
 * Display the visibility UV sampling.
 */
  ierr = ierr || u_uvplot(up, datcol, altcol);
/*
 * Reveal the re-displayed plot.
 */
  cpgebuf();
  return ierr;
}

/*.......................................................................
 * Private function of uvradplt to determine axis limits and draw and
 * label the plot axes.
 *
 * Input:
 *  up        U_par *  Plot-parameter block.
 *  axcol       int    The PGPLOT color index to use.
 * Output:
 *  return      int    0 - OK.
 */
static int u_axes(U_par *up, int axcol)
{
  const float margin=0.03; /* Fractional margin to leave around data points */
  const float topsep=0.7f; /* Separation of title from frame */
  Observation *ob;         /* The descriptor of the observation being plotted */
  float wmargin;   /* World-coordinate margin */
  float wxa,wxb;   /* Min,max X world coordinates */
  float wya,wyb;   /* Min,max Y world coordinates */
  char awrk[80];   /* Work string */
  char bwrk[80];   /* Work string */
/*
 * Get the descriptor of the Observation being plotted.
 */
  ob = up->ob;
/*
 * Get the data ranges to be plotted.
 */
  if(u_getrange(up))
    return 1;
/*
 * Check limits.
 */
  if(up->umax <= up->umin || up->vmax <= up->vmin) {
    fprintf(stderr, "uvplot: No data within proscribed ranges.\n");
    return 1;
  };
/*
 * Determine the window coordinates, leaving a margin around the
 * data. Note that by convention the (RA) X-axis is displayed going
 * from the its most positive to most negative coordinates.
 */
  wmargin = (up->umax - up->umin) * margin;
  wxa = up->umax + wmargin;
  wxb = up->umin - wmargin;
  wmargin = (up->vmax - up->vmin) * margin;
  wyb = up->vmax + wmargin;
  wya = up->vmin - wmargin;
/*
 * Set the plot window in the units required for the labels.
 */
  cpgsci(axcol);
  cpgsch(1.0f);
  cpgwnad(wavtouv(wxa), wavtouv(wxb), wavtouv(wya), wavtouv(wyb));
  cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
/*
 * Compose a title.
 */
  sprintf(awrk, "%.16s\\fr at \\fn%.3f GHz in %s  %s", ob->source.name,
	  getfreq(ob,-1)/1.0e9, Stokes_name(ob->stream.pol.type),
	  sutdate(ob->date.year, ob->date.ut, bwrk));
/*
 * Label the plot.
 */
  cpgmtxt("T", topsep, 0.0f, 0.0f, awrk);
/*
 * Compose axis labels.
 */
  sprintf(awrk, "U (%s)", uvwunits(U_PLAB));
  sprintf(bwrk, "V (%s)", uvwunits(U_PLAB));
  cpglab(awrk, bwrk, "");
/*
 * Plot the reference telescope name and its sub-array as a seconday title.
 */
  u_namplt(up, 0);
/*
 * Set the window with the units of the data.
 */
  cpgswin(wxa, wxb, wya, wyb);
  return 0;
}

/*.......................................................................
 * Given the indexes of a previously highlighted telescope and a new
 * telescope to highlight, re-plot the previous telescope in the normal
 * color and plot the baselines of the new telescope in highlighted color.
 *
 * Input:
 *  up          U_par *   Plot-parameter block.
 *  oper        Telop     U_ALLNEW  - Instate the new telescope spec from *init.
 *                        U_NXTSUB  - Highlight telescope of next sub-array.
 *                        U_NXTTEL  - Highlight the next telescope.
 *  forward       int     0 - Highlight the next telescope in order of
 *                            decreasing telescope index.
 *                        1 - Highlight the next telescope in order of
 *                            increasing telescope index.
 *  init      Telspec *  The new telescope specification to be used when
 *                       oper == U_ALLNEW.
 * Output:
 *  return        int    0 - OK.
 *                       1 - Error.
 */
static int u_newtel(U_par *up, Telop oper, int forward, Telspec *init)
{
  Telspec ts;      /* The new telescope specification */
  Observation *ob; /* The descriptor of the observation */
  int oldtel;      /* The index of the previous reference telescope */
  int oldsub;      /* The index of the previous reference telescope sub-array */
  int newtel;      /* The index of the new reference telescope */
  int newsub;      /* The index of the new reference telescope sub-array */
  int cif;         /* The index of the IF being processed */
  int ierr=0;      /* Error status */
/*
 * Get the descriptor of the observation.
 */
  ob = up->ob;
/*
 * Handle the specified change in reference telescope.
 */
  switch(oper) {
  case U_ALLNEW:
    ts = *init;
    if(next_tel(ob, FIND_FIRST, 1, 0, 0, 1, &ts))
      return 0;    
    break;
  case U_NXTSUB:
  case U_NXTTEL:
/*
 * If highlighting is currently disabled, re-enable it if the search
 * direction is forward. Otherwise ignore the request.
 */
    if(!up->highlight) {
      if(forward)
	ts = up->init;
      else
	return 0;
/*
 * Supplant the currently highlighted telescope.
 */
    } else {
      int iret=1;
      ts = up->ts;
/*
 * Locate the new telescope. Note the fallthrough between cases.
 */
      switch(oper) {
      case U_NXTTEL:
	iret = next_tel(ob, SKIP_TA, forward, 0, 0, 0, &ts);
      case U_NXTSUB:
	if(iret==1)
	  iret = next_tel(ob, SKIP_SUB, forward, 0, 0, 1, &ts);
      default:
	break;
      };
/*
 * Successful search?
 */
      if(iret == 0) {
	ts.nfix = 2;
      } else if(iret==1) {
	ts.nfix = forward ? 2:0;  /* Turn off highlighting? */
      } else {
	return 1;                 /* Error during search */
      };
    };
    break;
  default:
    lprintf(stderr, "u_newtel: Unrecognised opcode.\n");
    return 1;
  };
/*
 * Get the details of the currently hightlighted station and the
 * new station.
 * Note that an empty spec (nfix==0) denotes no highlighting.
 */
  oldtel = up->highlight ? up->ts.ta : -1;
  oldsub = up->highlight ? up->ts.isub : -1;
  newtel = ts.nfix != 0 ? ts.ta : -1;
  newsub = ts.nfix != 0 ? ts.isub : -1;
/*
 * Buffer PGPLOT commands.
 */
  cpgbbuf();
/*
 * Loop through all sampled IFs.
 */
  for(cif=0; (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
    Subarray *sub;  /* The descriptor of the sub-array being processed */
    int isub;       /* The index of sub in ob->sub[] */
/*
 * Get the next IF.
 */
    if(getIF(ob, cif)) {
      cpgebuf();         /* Release PGPLOT resources */
      return 1;          /* Error return */
    };
/*
 * Loop through sub-arrays of the new IF.
 */
    sub = ob->sub;
    for(isub=0; isub<ob->nsub; isub++,sub++) {
/*
 * Is either sub-array one of those to be changed?
 */
      int isoldsub = (isub == oldsub);
      int isnewsub = (isub == newsub);
      if(isoldsub || isnewsub) {
/*
 * Plot one baseline at a time.
 */
	Baseline *bptr = sub->base;
	int base;
	for(base=0; base<sub->nbase; base++,bptr++) {
/*
 * Is the new baseline associated with the previous reference telescope?
 */
	  int isold = isoldsub && (oldtel==bptr->tel_a || oldtel==bptr->tel_b);
/*
 * Is the new baseline associated with the new reference telescope?
 */
	  int isnew = isnewsub && (newtel==bptr->tel_a || newtel==bptr->tel_b);
/*
 * If the status of the baseline has changed, replot it in its new color.
 */
	  if(isold != isnew)
	    u_basepl(up, sub, base, isold ? datcol : altcol);
	};
      };
    };
  };
/*
 * Erase the old station title.
 */
  if(up->highlight)
    u_namplt(up, 1);
/*
 * Record the new highlighted telescope.
 */
  up->highlight = ts.nfix != 0;
  up->ts = ts;
/*
 * Plot the new station title.
 */
  if(up->highlight)
    u_namplt(up, 0);
/*
 * Terminate buffering.
 */
  cpgebuf();
  return ierr;
}

/*.......................................................................
 * Erase/plot reference telescope name.
 *
 * Input:
 *  up      U_par *  The plot descriptor.
 *  erase     int    If true, erase the telescope name.
 */
static void u_namplt(U_par *up, int erase)
{
/*
 * Is there a highlighted telescope?
 */
  if(up->highlight) {
    char title[81];  /* Temporary buffer to compose label in */
    int reftel;      /* The index of the telescope to be named */
    int refsub;      /* The index of the sub-array to be named */
/*
 * Get local copies of the indexes of the reference telescope and sub-array.
 */
    reftel = up->ts.ta;
    refsub = up->ts.isub;
/*
 * Telescope to name?
 */
    cpgsci(erase?0:1);
    sprintf(title, "%d:%s", refsub+1, up->ob->sub[refsub].tel[reftel].name);
    cpgmtxt("T", 1.0f, 1.0f, 1.0f, title);
    cpgsci(1);
  };
  return;
}

/*.......................................................................
 * Private function of uvradplt to plot the (U,V) locii of selected
 * visibilities.
 *
 * Input:
 *  up        U_par *  Plot-parameter block.
 *  othcol      int    The color to use for unhighlighted visibilities.
 *  refcol      int    The color to highlight visibilities of the reference
 *                     telescope.
 * Output:
 *  return      int    0 - OK.
 */
static int u_uvplot(U_par *up, int othcol, int refcol)
{
  Observation *ob;    /* The descriptor of the observation */
  int reftel;     /* The index of the reference telescope */
  int refsub;     /* The index of the reference sub-array */
  int cif;        /* The index of the IF being processed */
/*
 * Get the descriptor of the observation.
 */
  ob = up->ob;
/*
 * Get the highlighted reference telescope.
 */
  reftel = up->highlight ? up->ts.ta : -1;
  refsub = up->highlight ? up->ts.isub : -1;
/*
 * Buffer plotting operations.
 */
  cpgbbuf();
/*
 * Loop through all sampled IFs.
 */
  for(cif=0; (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
    Subarray *sub;  /* The descriptor of the sub-array being processed */
    int isub;       /* The index of sub in ob->sub[] */
/*
 * Get the next IF.
 */
    if(getIF(ob, cif)) {
      cpgebuf();         /* Release PGPLOT resources */
      return 1;          /* Error return */
    };
/*
 * Loop through sub-arrays of the new IF.
 */
    sub = ob->sub;
    for(isub=0; isub<ob->nsub; isub++,sub++) {
/*
 * Is the new sub-array the reference sub-array?
 */
      int isrefsub = (isub == refsub);
/*
 * Plot one baseline at a time.
 */
      Baseline *bptr = sub->base;
      int base;
      for(base=0; base<sub->nbase; base++,bptr++) {
/*
 * Is the new baseline associated with the reference telescope?
 */
	int isref = isrefsub && (reftel==bptr->tel_a || reftel==bptr->tel_b);
/*
 * Plot this baseline in the appropriate color.
 */
	u_basepl(up, sub, base, isref ? refcol : othcol);
      };
    };
  };
  cpgebuf();
  return 0;
}


/*.......................................................................
 * Private function of u_uvplot() to plot the UV positions of all
 * visibilities in the current IF and given baseline, in a specified
 * color.
 *
 * Input:
 *  up        U_par *  Plot-parameter block.
 *  sub    Subarray *  The descriptor of the sub-array to plot from.
 *  base        int    The baseline to plot.
 *  color       int    The PGPLOT color index to use.
 */
static void u_basepl(U_par *up, Subarray *sub, int base, int color)
{
  Integration *integ;   /* The descriptor of the integration being plotted */
  float uvscale;        /* Factor to scale UV corrds by to get wavelengths */
  float umin,umax;      /* The U axis range to plot */
  float vmin,vmax;      /* The V axis range to plot */
  int datsym;           /* Marker symbol */
  int ut;               /* The index of the integration being plotted */
/*
 * Get the conversion factor from UVW in light-seconds to wavelengths.
 */
  uvscale = up->ob->stream.uvscale;
/*
 * Get the range to be plotted.
 */
  umin = up->umin;
  umax = up->umax;
  vmin = up->vmin;
  vmax = up->vmax;
/*
 * Determine which marker symbol to use.
 */
  datsym = up->dobig ? 1 : -1;
/*
 * Plot the data.
 */
  cpgbbuf();
  cpgsci(color);
/*
 * Plot one integration's worth at a time.
 */
  integ = sub->integ;
  for(ut=0; ut<sub->ntime; ut++,integ++) {
    Visibility *vis = &integ->vis[base];
/*
 * Only plot good data.
 */
    if(!vis->bad) {
/*
 * Only plot data that lies inwards of the margins.
 */
      float uu = vis->u * uvscale;
      float vv = vis->v * uvscale;
/*
 * Plot the new point.
 */
      if(uu<umax && uu>umin && vv<vmax && vv>vmin)
	cpgpt(1, &uu, &vv, datsym);  /* Plot point */
/*
 * Plot the conjugate point.
 */
      if(up->doconj) {
	uu = -uu;
	vv = -vv;
	if(uu<umax && uu>umin && vv<vmax && vv>vmin)
	  cpgpt(1, &uu, &vv, datsym);
      };
    };
  };
  cpgebuf();
  return;
}


/*.......................................................................
 * Read the cursor and return the key that was entered. Also return the
 * position at which the cursor was pressed, bounded by the axis edges.
 *
 * Input:
 *  up      U_par *  The plot descriptor.
 *                   up->kp descriptor will be initialized with the
 *                   position and key selected by the user. On the first
 *                   call, set up->kp.key=KEY_NONE so that the cursor can
 *                   be positioned in the center of the plot.
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
 *  xref    float    xref and yref denote the reference position for
 *  yref    float    band-type cursors. They are ignored if mode==B_NORM.
 *  ci        int    The color index to draw band cursors with.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int u_cursor(U_par *up, Bandmode mode, float xref, float yref, int ci)
{
  static float xpos=0.0f; /* The X-world-coordinate of the cursor */
  static float ypos=0.0f; /* The Y-world-coordinate of the cursor */
  char key;               /* The key that caused cpgcurs() to return */
/*
 * Position the cursor in the center of the plot?
 */
  if(up->kp.key==KEY_NONE) {
    xpos = (up->umin + up->umax) / 2.0f;
    ypos = (up->vmin + up->vmax) / 2.0f;
  } else {
    xpos = up->kp.uu;
    ypos = up->kp.vv;
  };
/*
 * Substitute the cross-hair cursor for the normal cursor is requested.
 */
  if(mode == B_NORM && up->docross)
    mode = B_CROSS;
/*
 * Limit the given X and Y coordinates such that the cursor will appear within
 * the plot.
 */
  if(xpos < up->umin)
    xpos = up->umin;
  else if(xpos > up->umax)
    xpos = up->umax;
  if(ypos < up->vmin)
    ypos = up->vmin;
  else if(ypos > up->vmax)
    ypos = up->vmax;  
/*
 * Read a cursor position and the key that was pushed to return the
 * position.
 */
  cpgsci(ci);
  if(!cpgband((int) mode, 0, xref, yref, &xpos, &ypos, &key))
    return 1;
/*
 * Record the case of the enterred key.
 */
  up->kp.waslow = islower((int)key);
/*
 * Convert the selected key to upper case to simplify pending comparisons.
 */
  if(up->kp.waslow)
    key = toupper((int)key);
/*
 * Enforce plot bounds on the returned cursor position.
 */
  if(xpos < up->umin)
    xpos = up->umin;
  else if(xpos > up->umax)
    xpos = up->umax;
  if(ypos < up->vmin)
    ypos = up->vmin;
  else if(ypos > up->vmax)
    ypos = up->vmax;
/*
 * Record the results in the Keypos structure.
 */
  up->kp.key = key;
  up->kp.uu = xpos;
  up->kp.vv = ypos;
  return 0;
}

/*.......................................................................
 * Allow user selection of a sub-plot to zoom in to.
 *
 * Input:
 *  up        U_par *  Plot-parameter block.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int u_zoom(U_par *up)
{
  int npts=0;       /* The number of vertexes so far selected */
  float umin=0.0f;  /* The selected lower bound on U */
  float umax=0.0f;  /* The selected lower bound on U */
  float vmin=0.0f;  /* The selected lower bound on V */
  float vmax=0.0f;  /* The selected lower bound on V */
  float uref=0.0f;  /* The U coordinate of the first vertex selected */
  float vref=0.0f;  /* The V coordinate of the first vertex selected */
/*
 * Tell the user what to do.
 */
  printf(
    "Set the two opposite corners of the sub-plot. Press '%c' for help.\n",
    KEY_HELP);
/*
 * Have the user selected two points within the plot.
 */
  while(npts < 2) {
/*
 * Read the cursor.
 */
    if(u_cursor(up, npts==0 ? B_NORM:B_RECT, uref, vref, zoomcol))
      return 1;
/*
 * Act on the key used to return the cursor.
 */
    switch(up->kp.key) {
    case KEY_CUR:         /* Select position in plot */
/*
 * First position?
 */
      if(npts==0) {
	uref = umin = umax = up->kp.uu;
	vref = vmin = vmax = up->kp.vv;
      } else {
/*
 * Sort the selected U limits into ascending order.
 */
	if(up->kp.uu > umin)
	  umax = up->kp.uu;
	else
	  umin = up->kp.uu;
/*
 * Sort the selected V limits into ascending order.
 */
	if(up->kp.vv > vmin)
	  vmax = up->kp.vv;
	else
	  vmin = up->kp.vv;
      };
/*
 * Record the successful acquisition of a new corner.
 */
      npts++;
      break;
    case KEY_CAN:  /* Cancel selection */
      printf("Sub-plot selection cancelled.\n");
      return 0;
      break;
    case KEY_ZOOM: /* Select the full UV range */
      return u_setrange(up, 0, 0.0f, 0.0f, 0, 0.0f, 0.0f) || u_redisp(up);
      break;
    default:
      printf("You are currently selecting a sub-plot to display - use keys:\n");
      printf(" %c - Select %s corner of the area with this key.\n",
	     KEY_CUR, npts==0 ? "a" : "the second (opposite)");
      printf(" %c - Revert to the full plot range.\n", KEY_ZOOM);
      printf(" %c - Abort the selection with this key.\n", KEY_CAN);
      break;
    };
  };
/*
 * Check the selected range.
 */
  if(umin >= umax || vmin >= vmax) {
    printf("The sub-plot is too small to plot. Selection aborted.\n");
    return 0;
  };
/*
 * Install the new limits and re-display the plot.
 */
  return u_setrange(up, 1, umin, umax, 1, vmin, vmax) || u_redisp(up);
}

/*.......................................................................
 * Search for the data-point nearest the cursor.
 *
 * Input:
 *  up        U_par *  Plot-parameter block.
 *  xpos      float    The X-world coordinate of the cursor.
 *  ypos      float    The Y-world coordinate of the cursor.
 *  isamp       int    True if ypos is an amplitude, false if a phase.
 * Output:
 *  return   Bestvis   Container, containing the details of the
 *                     nearest point. Only use the results if best->found
 *                     is true.
 */
static Bestvis u_findpt(U_par *up, float xpos, float ypos)
{
  Observation *ob; /* The descriptor of the plotted observation */
  Bestvis best;    /* Container for the position of the closest point */
  int cif;         /* The index of the IF being searched */
  float uvscale;   /* Factor to scale UV coords by to get wavelengths */
  float minrsq=0.0f;/* Min distance squared between cursor and data point */
  float wxa,wxb;   /* World X-coordinate viewport limits */
  float wya,wyb;   /* World Y-coordinate viewport limits */
  float vxa,vxb;   /* Physical X viewport limits */
  float vya,vyb;   /* Physical X viewport limits */
  float xtomm;     /* Conversion factor between X world-coords and mm */
  float ytomm;     /* Conversion factor between Y world-coords and mm */
/*
 * No visibility has been found yet.
 */
  best.found = 0;
/*
 * Find the conversion factor between physical coordinates on each
 * axis and physical device coordinates in mm.
 */
  cpgqwin(&wxa, &wxb, &wya, &wyb);
  cpgqvp(2, &vxa, &vxb, &vya, &vyb);
  xtomm = (vxb-vxa) / (wxb-wxa);
  ytomm = (vyb-vya) / (wyb-wya);
/*
 * Get the descriptor of the plotted observation.
 */
  ob = up->ob;
/*
 * Loop through all sampled IFs.
 */
  for(cif=0; (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
    Subarray *sub;  /* The descriptor of the sub-array being processed */
    int isub;       /* The index of sub in ob->sub[] */
/*
 * Get the next IF.
 */
    if(getIF(ob, cif)) {
      cpgebuf();         /* Release PGPLOT resources */
      best.found = 0;
      return best;
    };
/*
 * Get the conversion factor from UVW in light-seconds to wavelengths.
 */
    uvscale = ob->stream.uvscale;
/*
 * Search each integration of each sub-array in the latest IF
 * for the closest point to the cursor.
 */
    sub = ob->sub;
    for(isub=0; isub<ob->nsub; isub++,sub++) {
      Integration *integ = sub->integ;
      int ut;
/*
 * Search integraions.
 */
      for(ut=0; ut<sub->ntime; ut++,integ++) {
	Visibility *vis = integ->vis;
	int base;
	for(base=0; base<sub->nbase; base++,vis++) {
/*
 * Only search good data.
 */
	  if(!vis->bad) {
/*
 * First the non-conjugate point.
 */
	    {
/*
 * Get the new UV coordinates.
 */
	      float uu = vis->u * uvscale;
	      float vv = vis->v * uvscale;
/*
 * Get the physical X offset between the cursor and the data point.
 */
	      float xoff = xtomm * (uu - xpos);
/*
 * Get the physical Y offset between the cursor and the data point.
 */
	      float yoff = ytomm * (vv - ypos);
/*
 * Compute the physical radius (squared) of the cursor from the point.
 */
	      float newrsq = xoff*xoff + yoff*yoff;
/*
 * Is the point within the plot bounds?
 */
	      if(uu >= up->umin && uu <= up->umax &&
		 vv >= up->vmin && vv <= up->vmax) {
/*
 * Is this closer to the cursor than any previous point.
 */
		if(!best.found || newrsq < minrsq) {
		  best.found = 1;
		  minrsq = newrsq;
		  best.ut = ut;
		  best.base = base;
		  best.isub = isub;
		  best.cif = cif;
		};
	      };
	    };
/*
 * Now check the conjugate point.
 */
	    if(up->doconj) {
/*
 * Get the new UV coordinates.
 */
	      float uu = -vis->u * uvscale;
	      float vv = -vis->v * uvscale;
/*
 * Get the physical X offset between the cursor and the data point.
 */
	      float xoff = xtomm * (uu - xpos);
/*
 * Get the physical Y offset between the cursor and the data point.
 */
	      float yoff = ytomm * (vv - ypos);
/*
 * Compute the physical radius (squared) of the cursor from the point.
 */
	      float newrsq = xoff*xoff + yoff*yoff;
/*
 * Is the point within the plot bounds?
 */
	      if(uu >= up->umin && uu <= up->umax &&
		 vv >= up->vmin && vv <= up->vmax) {
/*
 * Is this closer to the cursor than any previous point.
 */
		if(!best.found || newrsq < minrsq) {
		  best.found = 1;
		  minrsq = newrsq;
		  best.ut = ut;
		  best.base = base;
		  best.isub = isub;
		  best.cif = cif;
		};
	      };
	    };
	  };
	};
      };
    };
  };
/*
 * Nothing found?
 */
  if(!best.found)
    lprintf(stderr, "u_findpt: No data in range.\n");
/*
 * Return details of the nearest visibility.
 */
  return best;
}

/*.......................................................................
 * Flag a given visibility and erase it from the display.
 *
 * Input:
 *  up        U_par *  Plot-parameter block.
 *  cif         int    The index of the IF to be edited.
 *  isub        int    The index of the sub-array in which the visibility
 *                     resides.
 *  base        int    The index of the baseline of the visibility to be
 *                     flagged.
 *  ut          int    The index of the integration of the visibility
 *                     to be flagged.
 *  flag        int    If true, flag the point. If false, unflag the point.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int u_editpt(U_par *up, int cif, int isub, int base, int ut, int flag)
{
  Observation *ob; /* The descriptor of the observation being edited */
  Subarray *sub;   /* The descriptor of the sub-array being edited */
  Visibility *vis; /* The visibility being editted. */
  int flagged;     /* True if the visibility was originally flagged */
/*
 * Sanity check the arguments.
 */
  if(up==NULL) {
    lprintf(stderr, "u_editpt: NULL plot descriptor intercepted\n");
    return 1;
  };
/*
 * Get the descriptor of the observation.
 */
  ob = up->ob;
/*
 * Check the given sub-array index.
 */
  if(isub<0 || isub>=ob->nsub) {
    lprintf(stderr, "u_editpt: Out of range sub-array index.\n");
    return 1;
  };
/*
 * Check the given IF index.
 */
  if(cif<0 || cif>=ob->nif) {
    lprintf(stderr, "u_editpt: Out of range IF index.\n");
    return 1;
  };
/*
 * Get the desriptor of the associated sub-array.
 */
  sub = &ob->sub[isub];
/*
 * Check the baseline and integration indexes.
 */
  if(base<0 || base>=sub->nbase) {
    lprintf(stderr, "u_editpt: Out of range baseline index.\n");
    return 1;
  };
  if(ut<0 || ut>=sub->ntime) {
    lprintf(stderr, "u_editpt: Out of range integration index.\n");
    return 1;
  };
/*
 * Make sure that we have the required IF.
 */
  if(getIF(ob, cif))
    return 1;
/*
 * Mark the data as modified.
 */
  up->modified = 1;
/*
 * Get the visibility refered to.
 */
  vis = &sub->integ[ut].vis[base];
/*
 * Check its current status.
 */
  flagged = vis->bad;
/*
 * Flag or unflag the visibility if required.
 */
  if((flag!=0) == !flagged) {
    int oldcol;                         /* Used to preserve entry plot color */
    float uvscale = ob->stream.uvscale; /* UV dist Conversion to wavelengths */
    float uu = vis->u * uvscale;        /* The U coord of the visibility */
    float vv = vis->v * uvscale;        /* The V coord of the visibility */
    int datsym = up->dobig ? 1 : -1;    /* Symbol to plot point with */
/*
 * Edit the data.
 */
    ed_integ(ob, sub, ut, cif, flag, 1, 0, up->ch_ed, 1, base);
/*
 * Store the entry color.
 */
    cpgqci(&oldcol);
/*
 * Choose the color to erase flagged, or draw unflagged point.
 */
    if(flag)
      cpgsci(0);
    else if(up->highlight && isub==up->ts.isub &&
       (sub->base[base].tel_a==up->ts.ta || sub->base[base].tel_b==up->ts.ta))
      cpgsci(altcol);
    else
      cpgsci(datcol);
/*
 * Plot the point.
 */
    if(uu<up->umax && uu>up->umin && vv<up->vmax && vv>up->vmin)
      cpgpt(1, &uu, &vv, datsym);
/*
 * Plot its conjugate point.
 */
    uu = -uu;
    vv = -vv;
    if(uu<up->umax && uu>up->umin && vv<up->vmax && vv>up->vmin)
      cpgpt(1, &uu, &vv, datsym);
/*
 * Restore entry color.
 */
    cpgsci(oldcol);
  };
  return 0;
}

/*.......................................................................
 * Re-plot the mode line to reflect changes in edit mode.
 *
 * Input:
 *  up        U_par *  The plot descriptor.
 *  ch_ed       int    Select channel based editing if true.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int u_newmode(U_par *up, int ch_ed)
{
/*
 * Buffer until the new text has been plotted.
 */
  cpgbbuf();
/*
 * Erase the existing mode line.
 */
  u_mlab(up, 1);
/*
 * Install the new editing mode.
 */
  up->ch_ed = ch_ed;
/*
 * Draw the new mode line.
 */
  u_mlab(up, 0); /* Plot new mode line */
/*
 * reveal the changes.
 */
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Plot an extra mode label for editting sessions.
 *
 * Input:
 *  up       U_par *  The plot descriptor.
 *  erase       int    If true, erase existing mode label.
 * Output:
 *  return      int    0 - OK.
 */
static int u_mlab(U_par *up, int erase)
{
  Observation *ob;  /* The descriptor of the observation being plotted */
  int oldcol;       /* Temporary storage for entry color index */
  char label[81];   /* Temporary work string to compose mode label in */
/*
 * Get the descriptor of the observation.
 */
  ob = up->ob;
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
 * Compose the mode label.
 */
  sprintf(label, "Edit %s channels.", up->ch_ed ? "selected" : "all");
/*
 * Plot mode line.
 */
  cpgsch(1.0f);
  cpgmtxt("T", 2.5f, 0.0f, 0.0f, label);
/*
 * Restore entry color.
 */
  cpgsci(oldcol);
  return 0;
}

/*.......................................................................
 * Flag points within a user selected rectangular box.
 *
 * Input:
 *  up        U_par *  Plot-parameter block.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int u_edbox(U_par *up)
{
  Observation *ob;  /* The descriptor of the observation */
  int cif;          /* The index of the IF being processed */
  int npts=0;       /* The number of points so far selected */
  float umin=0.0f;  /* The U coordinate of the base of the selected area */
  float umax=0.0f;  /* The U coordinate of the top of the selected area */
  float vmin=0.0f;  /* The V coordinate of the base of the selected area */
  float vmax=0.0f;  /* The V coordinate of the top of the selected area */
  float uref=0.0f;  /* The U coordinate of the first vertex selected */
  float vref=0.0f;  /* The V coordinate of the first vertex selected */
  float uvscale;    /* Factor to scale UV coords by to get wavelengths */
/*
 * Tell the user what to do.
 */
  printf(
    "Set the two opposite corners of the area to flag. Press '%c' for help.\n",
    KEY_HELP);
/*
 * Have the user selected two points within the plot.
 */
  while(npts < 2) {
/*
 * Read the cursor.
 */
    if(u_cursor(up, npts==0 ? B_NORM:B_RECT, uref, vref, cutcol))
      return 1;
/*
 * Act on the key used to return the cursor.
 */
    switch(up->kp.key) {
    case KEY_CUR:         /* Select position in plot */
/*
 * First position?
 */
      if(npts==0) {
	uref = umin = umax = up->kp.uu;
	vref = vmin = vmax = up->kp.vv;
      } else {
/*
 * Record the new U and V coordiantes as the minimum or maximum of the
 * two chosen.
 */
	if(up->kp.uu > umin)
	  umax = up->kp.uu;
	else
	  umin = up->kp.uu;
/*
 * Record the new amplitude as the minimum or maximum of the two chosen.
 */
	if(up->kp.vv > vmin)
	  vmax = up->kp.vv;
	else
	  vmin = up->kp.vv;
      };
/*
 * Record the successful acquisition of a new corner.
 */
      npts++;
      break;
    case KEY_CAN:  /* Cancel selection */
      printf("Cut area selection cancelled.\n");
      return 0;
      break;
    default:
      printf("You are currently selecting a rectangular area to flag - use keys:\n");
      printf(" %c - Select %s corner of the area with this key.\n",
	     KEY_CUR, npts==0 ? "a" : "the second (opposite)");
      printf(" %c - Abort the selection with this key.\n", KEY_CAN);
      break;
    };
  };
/*
 * Get the descriptor of the observation.
 */
  ob = up->ob;
/*
 * Buffer PGPLOT display operations for speed.
 */
  cpgbbuf();
/*
 * Loop through all sampled IFs.
 */
  for(cif=0; (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
    Subarray *sub;  /* The descriptor of the sub-array being processed */
    int isub;       /* The index of sub in ob->sub[] */
/*
 * Get the next IF.
 */
    if(getIF(ob, cif)) {
      cpgebuf();         /* Release PGPLOT resources */
      return 1;          /* Error return */
    };
/*
 * Get the conversion factor from UVW in light-seconds to wavelengths.
 */
    uvscale = ob->stream.uvscale;
/*
 * Flag data within the selected area.
 */
    sub = ob->sub;
    for(isub=0; isub<up->ob->nsub; isub++,sub++) {
      Integration *integ = sub->integ;
      int ut;
/*
 * Search all integrations of subarray isub.
 */
      for(ut=0; ut<sub->ntime; ut++,integ++) {
	Visibility *vis = integ->vis;
	int base;
/*
 * Search visibilities of all baselines of integration ut of sub-array isub.
 */
	for(base=0; base<sub->nbase; base++,vis++) {
/*
 * Only search good data.
 */
	  if(!vis->bad) {
/*
 * Get the new UV coordinates.
 */
	    float uu = vis->u * uvscale;
	    float vv = vis->v * uvscale;
/*
 * Is the point within the selected rectangular area?
 * This requires checking both the normal and conjugate points.
 */
	    if(( uu >= umin &&  uu <= umax &&  vv >= vmin &&  vv <= vmax) ||
	       (-uu >= umin && -uu <= umax && -vv >= vmin && -vv <= vmax)) {
/*
 * Edit the point and erase it from the plot.
 */
	      if(u_editpt(up, cif, isub, base, ut, 1)) {
		cpgebuf();
		return 1;
	      };
	    };
	  };
	};
      };
    };
  };
/*
 * Reveal the edited plot.
 */
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Change plot ranges. Note that u_redisp() must be called after this
 * function. 
 *
 * Input:
 *  up     U_par *  The plot descriptor.
 *  fixu   int      If fixu != 0 and umin and umax describe a valid U
 *                  distance range then the U distance range will be
 *                  set by umin and umax. Otherwise it will be
 *                  determined from the data.
 *  umin   float    Minimum U distance (wavelengths).
 *  umax   float    Maximum U distance (wavelengths).
 *  fixv   int      If fixv != 0 and vmin and vmax describe a valid V
 *                  distance range then the V distance range will be
 *                  set by vmin and vmax. Otherwise it will be
 *                  determined from the data.
 *  vmin   float    Minimum V distance (wavelengths).
 *  vmax   float    Maximum V distance (wavelengths).
 * Output:
 *  return    int   0 - OK.
 *                  1 - Error.
 */
static int u_setrange(U_par *up, int fixu, float umin, float umax,
		                 int fixv, float vmin, float vmax)
{
/*
 * Check any given fixed U distance range.
 */
  if(fixu) {
/*
 * Ensure that umin < umax.
 */
    if(umin > umax) {float newmin=umax; umax=umin; umin=newmin;};
/*
 * If the range is too small, revert to autoscaling.
 */
    if(umin==umax)
      fixu = 0;
  };
/*
 * Check any given fixed V distance range.
 */
  if(fixv) {
/*
 * Ensure that vmin < vmax.
 */
    if(vmin > vmax) {float newmin=vmax; vmax=vmin; vmin=newmin;};
/*
 * If the range is too small, revert to autoscaling.
 */
    if(vmin==vmax)
      fixv = 0;
  };
/*
 * Record the requested ranges in the plot descriptor.
 */
  up->umin = fixu ? umin : 0.0f;
  up->umax = fixu ? umax : 0.0f;
  up->vmin = fixv ? vmin : 0.0f;
  up->vmax = fixv ? vmax : 0.0f;
/*
 * Record the autoscale status of the ranges.
 */
  up->fixu = fixu;
  up->fixv = fixv;
/*
 * Re-display the plot.
 */
  return 0;
}

/*.......................................................................
 * Private function of u_axes(), used to update the amplitude, phase and UV
 * radius plotting ranges from the values previously set by
 * u_setrange().
 *
 * Input:
 *  up    U_par *   The plot descriptor.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Error.
 */
static int u_getrange(U_par *up)
{
/*
 * The U and V ranges can be autoscaled, so if either quantity is
 * not fixed call uvrange() to determine the data range.
 */
  if(!up->fixu || !up->fixv) {
/*
 * Get the UV range of visibilities in all IFs.
 */
    UVrange *uvr = uvrange(up->ob, 1, 0, 0.0f ,0.0f);
    if(uvr==NULL)
      return 1;
/*
 * Use the determined ranges to set the un-specified parameters.
 */
    if(!up->fixu) {
      up->umin = -uvr->uvrmax;
      up->umax =  uvr->uvrmax;
    };
    if(!up->fixu) {
      up->vmin = -uvr->uvrmax;
      up->vmax =  uvr->uvrmax;
    };
  };
  return 0;
}

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "obs.h"
#include "units.h"
#include "vlbconst.h"
#include "telspec.h"
#include "visplot.h"
#include "vlbutil.h"
#include "cpgplot.h"
#include "logio.h"

static const double inc_pa=5.0*pi/180.0;/* Projection angle increment radians*/
static const int datcol=10;  /* PGPLOT color index for observed data */
static const int altcol=1;   /* Alternative to datcol for highlighting data */
static const int modcol=2;   /* PGPLOT color index for model data */
static const int axcol=1;    /* PGPLOT color index for axes */
static const int zoomcol=5;  /* PGPLOT color index for zoom cursor window */
static const int cutcol=2;   /* PGPLOT color index for edit cursor window */
static const int statcol=3;  /* PGPLOT color index for statistics box */
static const int dotsym= -1; /* Default - small plot marker symbol */
static const int bigsym=1;   /* Alternate - big plot marker symbol */

/* Define selection keys */

enum {
  KEY_NONE='\0',  /* Null key */
  KEY_DOT ='.',   /* Key to toggle marker symbol size */
  KEY_AMP ='1',   /* Key to select just an amplitude plot */
  KEY_PHS ='2',   /* Key to select just a phase plot */
  KEY_BOTH='3',   /* Key to select both amplitude and phase plots */
  KEY_INC ='>',   /* Increment projection position angle */
  KEY_DEC ='<',   /* Decrement projection position angle */
  KEY_ANG ='?',   /* Key to request keyboard entry of projection angle */
  KEY_CUR ='A',   /* Key to enter positions */
  KEY_CUT ='C',   /* Key to initiate cut area selection */
  KEY_CAN ='D',   /* Key to cancel UV display range selection */
  KEY_ERR ='E',   /* Toggle whether to show an error sub-plot */
  KEY_HELP='H',   /* Key to request help */
  KEY_DIS ='L',   /* Key to request re-display of plot */
  KEY_MOD ='M',   /* Key to toggle model display */
  KEY_NXT ='N',   /* Key to highlight next station */
  KEY_PRV ='P',   /* Key to highlight previous station */
  KEY_SHOW='S',   /* Report on the nearest point to the cursor */
  KEY_TEL ='T',   /* Key to select highlighted telescope by name */
  KEY_UVR ='U',   /* Select a new UV range */
  KEY_VEC ='V',   /* Show the vector average stats of a given region */
  KEY_CH  ='W',   /* Toggle channel editing mode */
  KEY_QUIT='X',   /* Key to quit interactive session */
  KEY_ZOOM='Z',   /* Key to select a new subplot display range */
  KEY_CROSS='+',  /* Toggle cross-hair cursor mode */
  KEY_DIFF='-'    /* Toggle to and from plotting the difference between */
                  /*  the data and the model. */
};

typedef enum {   /* Enumerate the supported types of sub-plots */
  AMP_PLOT,      /* True when refering to an amplitude plot */
  PHS_PLOT,      /* True when refering to a phase plot */
  ERR_PLOT       /* True when refering to a error plot */
} RpType;

/* Type used to return cursor selection */

typedef struct {
  float xpos;    /* NDC X coord of last entered cursor position */
  float ypos;    /* NDC Y coord of last entered cursor position */
  float uvdist;  /* Cursor selected UV radius (wavelengths) */
  float value;   /* Cursor selected amplitude or phase (raw data units) */
  int key;       /* The upper-case character corresponding to the key pressed */
  int waslow;    /* The original case of the selection key */
  RpType plot;   /* The plot of the selected position */
} Keypos;

typedef struct {
  Observation *ob; /* The descriptor of the observation being plotted */
  Keypos kp;       /* Descriptor of the last enterred cursor selection  */
  Telspec init;    /* The specification of the first available telescope */
  Telspec ts;      /* Specification of the reference telescope */
  int highlight;   /* True when telescope highlighting is enabled */
  float vxa,vxb;   /* NDC Y-coordinates of the viewport edges */
  float vya,vyb;   /* NDC X-coordinates of the viewport edges */
  float vatop;     /* NDC Y-coord of the top of the amplitude viewport */
  float vabot;     /* NDC Y-coord of the bottom of the amplitude viewport */
  float vptop;     /* NDC Y-coord of the top of the phase viewport */
  float vpbot;     /* NDC Y-coord of the bottom of the phase viewport */
  float vetop;     /* NDC Y-coord of the top of the error viewport */
  float vebot;     /* NDC Y-coord of the bottom of the error viewport */
  float wxa,wxb;   /* World X-coordinates corresponding to vxa,vxb */
  float wyaa,wyab; /* World Y-coordinates of ampitudes at vabot,vatop */
  float wypa,wypb; /* World Y-coordinates of phases at vpbot,vptop */
  float wyea,wyeb; /* World Y-coordinates of errors at vebot,vetop */
  float uvmin;     /* Minimum displayed UV radius (excludes margins) */
  float uvmax;     /* Maximum displayed UV radius (excludes margins) */
  float ampmin;    /* Minimum displayed amplitude (excludes margins) */
  float ampmax;    /* Maximum displayed amplitude (excludes margins) */
  float phsmin;    /* Minimum displayed phase (excludes margins) */
  float phsmax;    /* Maximum displayed phase (excludes margins) */
  float errmin;    /* Minimum displayed error (excludes margins) */
  float errmax;    /* Maximum displayed error (excludes margins) */
  struct {         /* UV-radius projection attributes filled by r_newphi() */
    double phi;    /* Projection angle (radians) */
    float sinphi;  /* sin(phi) */
    float cosphi;  /* sin(phi) */
  } proj;
  int doproj;      /* If true, project UV radii onto position angle 'phi' */
  int doamp;       /* True if amplitudes are plotted */
  int dophs;       /* True if phases are plotted */
  int doerr;       /* True if error are plotted */
  int fixuvr;      /* If true don't autoscale the UV radius range */
  int fixamp;      /* If true don't autoscale the amplitude range */
  int fixphs;      /* If true don't autoscale the phase range */
  int fixerr;      /* If true don't autoscale the error range */
  int docurs;      /* If true allow cursor interaction */
  int domod;       /* If true, plot the model */
  int dobig;       /* If true, plot with larger dot size */
  int docross;     /* True to enable cross-hair mode */  
  int dodiff;      /* If true, plot the difference between data and model */
  int modified;    /* Remains 0 unless the data are edited */
  int ch_ed;       /* If true, edits are restricted to current freq channels */
} R_par;

static int r_redisp(R_par *rp, int newpage);
static int r_plvis(R_par *rp);
static int r_axes(R_par *rp, int axcol);
static int r_basepl(R_par *rp, Subarray *sub, int base, int color);
static int r_modpl(R_par *rp, Subarray *sub, int base, int modcol);
static int r_zoom(R_par *rp);
static int r_newuvr(R_par *rp);
static void r_namplt(R_par *rp, int erase);
static int r_editpt(R_par *rp, int cif, int isub, int base, int ut, int flag);
static int r_edbox(R_par *rp);
static int r_scalar_stats(R_par *rp);
static int r_vector_stats(R_par *rp);
static int r_setrange(R_par *rp, int fixuvr, float uvmin, float uvmax,
		      int fixamp, float ampmin, float ampmax,
		      int fixphs, float phsmin, float phsmax,
		      int fixerr, float errmin, float errmax);
static int r_newphi(R_par *rp, double phi, int update);
static int r_getrange(R_par *rp);
static int r_newmode(R_par *rp, int ch_ed);
static int r_mlab(R_par *rp, int erase);
static double wrapphi(double phi);
static int r_ampwin(R_par *rp);
static int r_phswin(R_par *rp);
static int r_errwin(R_par *rp);
static int r_vpwin(R_par *rp);
static float r_uvdist(R_par *rp, float u, float v);
static float r_vis_amp(R_par *rp, Visibility *vis);
static float r_vis_phs(R_par *rp, Visibility *vis, float u, float v);
static float r_vis_err(R_par *rp, Visibility *vis);
static float r_mod_phs(R_par *rp, Visibility *vis, float u, float v);
static int r_flags(R_par *rp, int key, int waslow);
static int r_getphi(R_par *rp);

typedef enum {R_ALLNEW, R_NXTSUB, R_NXTTEL} Telop;

static int r_newtel(R_par *rp, Telop oper, int forward, Telspec *init);

/* Define a structure for the return type of r_findpt() */

typedef struct {
  int found;     /* True only if a point was selected */
  int ut;        /* The index of the integration of the closest point */
  int base;      /* The index of the baseline of the closest point */
  int isub;      /* The index of the sub-array of the closest point */
  int cif;       /* The index of the IF of the closest point */
} Bestvis;

static Bestvis r_findpt(R_par *rp, float xpos, float ypos, RpType plot);

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

static int r_cursor(R_par *rp, Bandmode mode, RpType plot,
		    float xref, float yref, int ci);

/*.......................................................................
 * Plot observed and model visibilities versus UV radius.
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
 *  doproj      int    If true, plot projected radii.
 *  phi       float    The position angle (radians) to project UV radii onto
 *                     if doproj is true, North through East wrt the V axis.
 *  uvmin     float    The minimum UV radius to be displayed (wavelengths)
 *  uvmax     float    The maximum UV radius to be displayed (wavelengths)
 *                     If uvmin>=uvmax, auto-ranging will be substituted.
 *  ampmin    float    The minimum amplitude to be displayed.
 *  ampmax    float    The maximum amplitude to be displayed.
 *                     If ampmin>=ampmax, auto-ranging will be substituted.
 *  phsmin    float    The minimum phase to be plotted (radians).
 *  phsmax    float    The maximum phase to be plotted (radians).
 *                     If phsmin>=phsmin, -pi -> pi will be subsituted.
 * Input/Output:
 *  modified    int *   If modified!=NULL then *modified will be assigned
 *                      the value 1 if the data were edited, or 0 otherwise.
 * Output:
 *  return      int    0 - OK.
 */
int uvradplt(Observation *ob, Telspec *ts, int docurs, char *opts,
	     int doproj, float phi, float uvmin, float uvmax,
	     float ampmin, float ampmax, float phsmin, float phsmax,
	     int *modified)
{
  R_par rp;         /* Internal parameter container */
  Bestvis best;     /* Container for the results of a visibility search */
  char answer[10];  /* String to hold answers to PGPLOT inquiries */
  int slen;         /* Length of answer in 'answer' */
  int oldcol;       /* Color on entry to this function */
  int old_if;       /* State of current IF to be restored on exit */
  int ierr=0;       /* Error status flag */
  int i;
/*
 * Data not yet modified.
 */
  if(modified)
    *modified = 0;
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_SELECT, "uvradplt"))
    return 1;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Record the Observation descriptor.
 */
  rp.ob = ob;
/*
 * Parse input paramaters.
 */
  rp.domod = 0;
  rp.doproj = doproj;
  rp.doamp = 1;
  rp.dophs = 0;
  rp.doerr = 0;
/*
 * Get the specification of the first available telescope.
 */
  {
    Telspec *init = find_tel(ob, 0, 0, 0, 1, 0, 0, 1);
    if(!init)
      return 1;
    rp.init = *init;
    rp.init.nfix = 2;
  };
/*
 * Record the initial highlight telescope specification.
 */
  if(ts && ts->nfix>0) {
    if(next_tel(ob, FIND_FIRST, 1, 0, 0, 1, ts))
      return 1;
    rp.ts = *ts;
    rp.highlight = 1;
  } else {
    rp.ts = rp.init;
    rp.highlight = 0;
  };
/*
 * Flag the data as un-modified.
 */
  rp.modified = 0;
/*
 * Default to global channel editing.
 */
  rp.ch_ed = 0;
/*
 * Default to small dot markers.
 */
  rp.dobig = 0;
  rp.docross = 0;
/*
 * Default to plotting the data and the model separately.
 */
  rp.dodiff = 0;
/*
 * If cursor interaction has been requested, check whether the current
 * device has a cursor.
 */
  if(docurs) {
    slen = sizeof(answer)-1;
    cpgqinf("CURSOR", answer, &slen);
    docurs = strncmp(answer,"YES",3) == 0;
  };
  rp.docurs = docurs;
/*
 * Store the entry color.
 */
  cpgqci(&oldcol);
/*
 * Assign the start projection angle.
 */
  if(r_newphi(&rp, phi, 0))
    return 1;
/*
 * If a string of flag options was given, interpret them here.
 */
  if(opts != NULL) {
    slen = strlen(opts);
    for(i=0; i<slen; i++) {
      int key = opts[i];
      int waslow = islower(key);
      if(waslow) key = toupper(key);
      r_flags(&rp, key, waslow);
/*
 * Handle flags that are normally treated as commands.
 */
      switch(key) {
      case KEY_DOT:
	rp.dobig = !rp.dobig;
	break;
      case KEY_CROSS:
	rp.docross = !rp.docross;
	break;
      };
    };
  };
/*
 * Set the initial amplitude and UV radius plot ranges and display the
 * initial plot.
 */
  if(r_setrange(&rp, 1, uvmin, uvmax, 1, ampmin, ampmax, 1, phsmin, phsmax,
		0, 0.0f, 0.0f) ||
     r_redisp(&rp, 1))
    return 1;
/*
 * Interactive mode?
 */
  if(rp.docurs) {
/*
 * Tell user how to obtain key listing.
 */
    lprintf(stdout,
	    "Move the cursor into the plot window and press \'%c\' for help\n",
	    KEY_HELP);
/*
 * Start the interactive loop.
 */
    rp.kp.key = KEY_NONE;
    while(!ierr && rp.kp.key != KEY_QUIT) {
      int wasflag=0;  /* True if last key stroke toggled a flag */
      int nflag=0;    /* Number of flag toggling operations done */
/*
 * Read the cursor.
 */
      do {
/*
 * Read a cursor position and the key that was pushed to return the
 * position.
 */
	if(r_cursor(&rp, B_NORM, 0, 0.0f, 0.0f, 1)) {
	  ierr = 1;
	  wasflag = 0;
	  nflag = 0;
	} else {
/*
 * Toggle flags where appropriate.
 */
	  wasflag = r_flags(&rp, rp.kp.key, rp.kp.waslow) == 0;
	  nflag += wasflag;
	};
      } while(wasflag);   /* Don't do anything more if a flag was toggled */
/*
 * Take action appropriate to the key that the user pressed.
 */
      if(nflag > 0) {  /* Update display after a sequence of flag toggling */
	nflag = 0;
	ierr = r_redisp(&rp, 1);
      } else if(!ierr) {
/*
 * Interpret the key pressed.
 */
	switch(rp.kp.key) {
	case KEY_DIS:
	  ierr = r_redisp(&rp, 1);
	  break;
	case KEY_DOT:
	  rp.dobig = !rp.dobig;
	  ierr = r_redisp(&rp, 1);
	  break;
	case KEY_INC:
	  ierr = rp.doproj && r_newphi(&rp, rp.proj.phi + inc_pa, 1);
	  break;
	case KEY_DEC:
	  ierr = rp.doproj && r_newphi(&rp, rp.proj.phi - inc_pa, 1);
	  break;
	case KEY_ANG:
	  ierr = rp.doproj && r_getphi(&rp);
	  break;
	case KEY_NXT:
	  ierr = r_newtel(&rp, rp.kp.waslow ? R_NXTTEL:R_NXTSUB, 1, NULL);
	  break;
	case KEY_PRV:
	  ierr = r_newtel(&rp, rp.kp.waslow ? R_NXTTEL:R_NXTSUB, 0, NULL);
	  break;
	case KEY_TEL:
	  {
	    Telspec *ts = read_Telspec(rp.ob, NULL, NULL, rp.ts.isub);
	    ierr = ts && r_newtel(&rp, R_ALLNEW, 1, ts);
	  };
	  break;
	case KEY_SHOW:
	  if(rp.kp.waslow) {
	    best = r_findpt(&rp, rp.kp.uvdist, rp.kp.value, rp.kp.plot);
	    if(best.found) {
	      Subarray *subptr = &rp.ob->sub[best.isub];
	      Baseline *bptr = &subptr->base[best.base];
	      char date_str[20];
	      write_ut(subptr->integ[best.ut].ut, sizeof(date_str), date_str);
	      printf("Visibility on baseline %d:%s-%s (IF %d) at UT %s\n",
		     best.isub+1, subptr->tel[bptr->tel_a].name,
		     subptr->tel[bptr->tel_b].name, best.cif+1, date_str);
	    };
	  } else {
	    r_scalar_stats(&rp);
	  };
	  break;
	case KEY_VEC:
	  r_vector_stats(&rp);
	  break;
	case KEY_CUR:  /* Flag the nearest point to the cursor */
	  best = r_findpt(&rp, rp.kp.uvdist, rp.kp.value, rp.kp.plot);
	  ierr = best.found && r_editpt(&rp, best.cif, best.isub, best.base, best.ut, 1);
	  break;
	case KEY_CUT:  /* Initiate cut area selection */
	  ierr = r_edbox(&rp);
	  break;
	case KEY_CH:       /* Toggle channel editing mode */
	  ierr = r_newmode(&rp, !rp.ch_ed);
	  break;
	case KEY_CROSS:   /* Toggle cross-hair cursor mode */
	  rp.docross = !rp.docross;
	  break;
	case KEY_HELP:
	  printf("You requested help by pressing \'%c\'.\n", KEY_HELP);
	  printf("The following keys are defined when pressed inside the plot:\n");
	  printf(" %c - Quit radplt\n", KEY_QUIT);
	  printf(" %c - Re-display whole plot\n", KEY_DIS);
	  printf(" %c - Re-display plot with alternate marker symbol.\n", KEY_DOT);
	  printf(" %c - Highlight next telescope\n", tolower(KEY_NXT));
	  printf(" %c - Highlight previous telescope\n", tolower(KEY_PRV));
	  printf(" %c - Step to the next sub-array to highlight.\n", KEY_NXT);
	  printf(" %c - Step to the preceding sub-array to highlight.\n", KEY_PRV);
	  printf(" %c - Specify highlighted telescope from keyboard\n", KEY_TEL);
	  printf(" %c - Show the baseline and time of the nearest point to the cursor\n", tolower(KEY_SHOW));
	  printf(" %c - Show the amp/phase statistics of the data within a selected area.\n", KEY_SHOW);
	  printf(" %c - Show the real/imag statistics of the data within a selected area.\n", KEY_VEC);
	  printf(" %c - (Left-mouse-button) Flag the point closest to the cursor\n", KEY_CUR);
	  printf(" %c - Initiate selection of an area to flag.\n", KEY_CUT);
	  printf(" %c - Toggle spectral-line channel based editing.\n",KEY_CH);
	  printf(" %c - Select a new amplitude or phase display range.\n", KEY_ZOOM);
	  printf(" %c - Select a new UV-radius display range.\n", KEY_UVR);
	  if(rp.doproj) {
	    printf("Projection angle selection:\n");
	    printf(" %c - Enter a new projection angle from the keyboard.\n",
		   KEY_ANG);
	    printf(" %c - Decrease the projection angle by 5 degrees.\n",
		   KEY_DEC);
	    printf(" %c - Increase the projection angle by 5 degrees.\n",
		   KEY_INC);
	  };
	  printf("Display mode options:\n");
	  printf(" %c - Toggle model plotting.\n", KEY_MOD);
	  printf(" %c - Display amplitude only.\n", KEY_AMP);
	  printf(" %c - Display phase only.\n", KEY_PHS);
	  printf(" %c - Display amplitude and phase.\n", KEY_BOTH);
	  printf(" %c - Toggle whether to display an error plot.\n", KEY_ERR);
	  printf(" %c - Toggle whether to display residuals.\n", KEY_DIFF);
	  printf(" %c - Toggle whether to use a cross-hair cursor if available.\n", KEY_CROSS);
	  break;
	case KEY_UVR:  /* Allow user selection of a new UV display range */
	  ierr = r_newuvr(&rp);
	  break;
	case KEY_ZOOM:
	  ierr = r_zoom(&rp);
	  break;
	};
      };
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
    *modified = rp.modified;
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
 *  rp        R_par *  Plot-parameter block.
 *  newpage     int    If true, start a new page before plotting.
 *                     Otherwise erase what is there before plotting.
 * Output:
 *  return      int    0 - OK.
 */
static int r_redisp(R_par *rp, int newpage)
{
  int ierr=0;      /* Error status flag */
/*
 * Buffer PGPLOT commands until finished.
 */
  cpgbbuf();
/*
 * Start a new page?
 */
  if(newpage) {
    cpgpage();
/*
 * Erase the plot area by plotting an erase rectangle over it.
 */
  } else {
    cpgsvp(0.0,1.0,0.0,1.0);
    cpgsfs(1);
    cpgsci(0);
    cpgswin(0.0,1.0,0.0,1.0);
    cpgrect(0.0,1.0,0.0,1.0);
    cpgsci(1);
  };
/*
 * Assign the standard viewport.
 */
  cpgsch(1.0f);
  cpgvstd();
/*
 * Plot axes.
 */
  ierr = ierr || r_axes(rp, axcol);
/*
 * Plot the mode line if interactive.
 */
  ierr = ierr || (rp->docurs && r_mlab(rp, 0));
/*
 * Plot the model and observed data.
 */
  ierr = ierr || r_plvis(rp);
/*
 * Reveal the re-displayed plot.
 */
  cpgebuf();
  return ierr;
}

/*.......................................................................
 * Plot the model and data within the current plot.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 * Output:
 *  return      int    0 - OK.
 */
static int r_plvis(R_par *rp)
{
  Observation *ob; /* The descriptor of the observation */
  int refsub;      /* The existing reference sub-array index */
  int reftel;      /* The existing reference telescope index */
  int cif;         /* The index of the IF being processed */
  int ierr=0;      /* Error status flag */
/*
 * Get the descriptor of the observation.
 */
  ob = rp->ob;
/*
 * Get the highlighted reference telescope.
 */
  reftel = rp->highlight ? rp->ts.ta : -1;
  refsub = rp->highlight ? rp->ts.isub : -1;
/*
 * Buffer PGPLOT commands until finished.
 */
  cpgbbuf();
/*
 * Loop through all sampled IFs.
 */
  for(cif=0; !ierr && (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
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
    for(isub=0; !ierr && isub<ob->nsub; isub++,sub++) {
/*
 * Is the new sub-array the reference sub-array?
 */
      int isrefsub = (isub == refsub);
/*
 * Plot one baseline at a time.
 */
      Baseline *bptr = sub->base;
      int base;
      for(base=0; !ierr && base<sub->nbase; base++,bptr++) {
/*
 * Is the new baseline associated with the reference telescope?
 */
	int isref = isrefsub && (reftel==bptr->tel_a || reftel==bptr->tel_b);
/*
 * Plot the model.
 */
	ierr = r_modpl(rp, sub, base, modcol);
/*
 * Plot this baseline in the appropriate color.
 */
	ierr = ierr || r_basepl(rp, sub, base, isref ? altcol : datcol);
      };
    };
  };
/*
 * Reveal the plot.
 */
  cpgebuf();
  return ierr;
}

/*.......................................................................
 * Private function of uvradplt to determine axis limits and draw and
 * label the plot axes. Calling this function initializes the following
 * members of *rp: vxa,vxb,vya,vyb,vabot,vatop,vpbot,vptop,vebot,vetop,
 * wxa,wxb,wyaa,wyab,wypa,wypb,uvmin,uvmax,ampmin,ampmax.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 *  axcol       int    The PGPLOT color index to use.
 * Output:
 *  return      int    0 - OK.
 */
static int r_axes(R_par *rp, int axcol)
{
  const float margin=0.1f; /* Fractional margin to leave around data points */
  const float topsep=0.7f; /* Separation of title from frame */
  const float botsep=2.7f; /* Separation of X-axis label from frame */
  const float lhssep=2.7f; /* Separation of Y-axis label from frame */
  Observation *ob;         /* The descriptor of the observation being plotted */
  char awrk[80];   /* Work string */
  char bwrk[80];   /* Work string */
/*
 * Get the descriptor of the Observation being plotted.
 */
  ob = rp->ob;
/*
 * Get the viewport limits.
 */
  if(r_vpwin(rp))
    return 1;
/*
 * Get the data limits to be plotted.
 */
  if(r_getrange(rp))
    return 1;
/*
 * Check limits.
 */
  if(rp->uvmax <= 0.0f || rp->ampmax <= 0.0f) {
    fprintf(stderr, "uvradplt: No data within selected UV range.\n");
    return 1;
  };
/*
 * Determine the raw UV coordinate limits including margins.
 */
  rp->wxa = rp->uvmin - (rp->uvmax-rp->uvmin)*margin;
  rp->wxb = rp->uvmax + (rp->uvmax-rp->uvmin)*margin;
/*
 * Ensure that the UV axis range is finite.
 */
  if(rp->wxb - rp->wxa < rp->wxb/1000.0f) {
    float extra = rp->wxb/500.0f;
    if(extra <= 0.0f) extra = 1.0f;
    rp->wxa -= extra;
    rp->wxb += extra;
  };
/*
 * Get the amplitude axis limits.
 */
  rp->wyaa = rp->ampmin - (rp->ampmax - rp->ampmin) * margin;
  rp->wyab = rp->ampmax + (rp->ampmax - rp->ampmin) * margin;
/*
 * Ensure that the amplitude axis range is finite.
 */
  if(rp->wyab - rp->wyaa < rp->wyab/1000.0f) {
    float extra = rp->wyab/500.0f;
    if(extra <= 0.0f) extra = 1.0f;
    rp->wyaa -= extra;
    rp->wyab += extra;
  };
/*
 * Get the phase axis limits.
 */
  rp->wypa = rp->phsmin - (rp->phsmax - rp->phsmin) * margin;
  rp->wypb = rp->phsmax + (rp->phsmax - rp->phsmin) * margin;
/*
 * Ensure that the phase axis range is finite.
 */
  if(rp->wypb - rp->wypa < 0.1f * dtor) {
    float extra = 1.0f * dtor;
    rp->wypa -= extra;
    rp->wypb += extra;
  };
/*
 * Get the error axis limits.
 */
  rp->wyea = rp->errmin - (rp->errmax - rp->errmin) * margin;
  rp->wyeb = rp->errmax + (rp->errmax - rp->errmin) * margin;
/*
 * Ensure that the error axis range is finite.
 */
  if(rp->wyeb - rp->wyea < rp->wyeb/1000.0f) {
    float extra = rp->wyeb/500.0f;
    if(extra <= 0.0f) extra = 1.0f;
    rp->wyea -= extra;
    rp->wyeb += extra;
  };
/*
 * Set a viewport enclosing all of the plots.
 */
  cpgsci(axcol);
  cpgsch(1.0f);
  cpgsvp(rp->vxa, rp->vxb, rp->vya, rp->vyb);
/*
 * Compose and plot a title.
 */
  sprintf(awrk, "%.16s\\fr at \\fn%.3f GHz in %s  %s",
	  ob->source.name, getfreq(ob,-1)/1.0e9,
	  Stokes_name(ob->stream.pol.type),
	  sutdate(ob->date.year, ob->date.ut, bwrk));
  cpgmtxt("T", topsep, 0.0f, 0.0f, awrk);
/*
 * Plot the reference telescope name as a secondary title.
 */
  r_namplt(rp, 0);
/*
 * Label the X axis.
 */
  if(rp->doproj) {
    sprintf(awrk, "Radial UV distance along P.A. %g^  (%s)",
	    rp->proj.phi*rtod, uvwunits(U_PLAB));
  } else {
    sprintf(awrk, "UV radius  (%s)", uvwunits(U_PLAB));
  };
  cpgmtxt("B", botsep, 0.5f, 0.5f, awrk);
/*
 * Draw the error plot axes if required.
 */
  if(rp->doerr) {
    r_errwin(rp);
    cpgswin(wavtouv(rp->wxa), wavtouv(rp->wxb), rp->wyea, rp->wyeb);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgmtxt("L", lhssep, 0.5f, 0.5f, "Error");
  };
/*
 * Draw the phase plot axes if required.
 */
  if(rp->dophs) {
    r_phswin(rp);
    cpgswin(wavtouv(rp->wxa), wavtouv(rp->wxb), rp->wypa*rtod, rp->wypb*rtod);
    cpgbox(rp->doerr ? "BCST":"BCNST",0.0,0,"BCNST",0.0,0);
    cpgmtxt("L", lhssep, 0.5f, 0.5f, rp->dodiff ? "Residual phase" : "Phase");
  };
/*
 * Draw the amplitude plot axes if required.
 */
  if(rp->doamp) {
    r_ampwin(rp);
    cpgswin(wavtouv(rp->wxa), wavtouv(rp->wxb), rp->wyaa, rp->wyab);
    cpgbox(rp->dophs || rp->doerr ? "BCST":"BCNST", 0.0, 0, "BCNST", 0.0, 0);
    cpgmtxt("L", lhssep, 0.5f, 0.5f, rp->dodiff ? "Residual amplitude" :
	    "Amplitude");
  };
  return 0;
}

/*.......................................................................
 * Private function of uvradplt to plot the observed amplitude of
 * a specified baseline versus UV radius, using a specified color.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 *  sub    Subarray *  The descriptor of the sub-array to plot from.
 *  base        int    The baseline to plot.
 *  color       int    The PGPLOT color index to use.
 * Output:
 *  return      int    0 - OK.
 */
static int r_basepl(R_par *rp, Subarray *sub, int base, int color)
{
  Integration *integ; /* The descriptor of the integration being plotted */
  float uvscale;      /* Factor to scale UV corrds by to get wavelengths */
  int datsym;         /* Marker symbol */
  int ut;             /* The ut being looked at */
/*
 * Get the conversion factor from UVW in light-seconds to wavelengths.
 */
  uvscale = rp->ob->stream.uvscale;
/*
 * Determine which marker symbol to use.
 */
  datsym = rp->dobig ? bigsym : dotsym;
/*
 * Plot the data.
 */
  cpgbbuf();
  cpgsci(color);
/*
 * Start with the amplitude plot.
 */
  if(rp->doamp && r_ampwin(rp)==0) {
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
 * Compare the UV radius and visibility amplitude with the cut-off limits.
 */
	float uu = vis->u * uvscale;
	float vv = vis->v * uvscale;
	float uvdist = r_uvdist(rp, uu, vv);
	float amp = r_vis_amp(rp, vis);
	if(uvdist >= rp->uvmin && uvdist <= rp->uvmax &&
	   amp >= rp->ampmin && amp <= rp->ampmax) {
	  cpgpt(1, &uvdist, &amp, datsym);  /* Plot point */
	};
      };
    };
  };
/*
 * Now do the phase plot.
 */
  if(rp->dophs && r_phswin(rp)==0) {
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
 * Compare the UV radius and visibility phase with the cut-off limits.
 */
	float uu = vis->u * uvscale;
	float vv = vis->v * uvscale;
	float uvdist = r_uvdist(rp, uu, vv);
	float phs = r_vis_phs(rp, vis, uu, vv);
	if(uvdist >= rp->uvmin && uvdist <= rp->uvmax &&
	   phs >= rp->phsmin && phs <= rp->phsmax) {
	  cpgpt(1, &uvdist, &phs, datsym);  /* Plot point */
	};
      };
    };
  };
/*
 * Now do the error plot.
 */
  if(rp->doerr && r_errwin(rp)==0) {
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
 * Compare the UV radius and visibility error with the cut-off limits.
 */
	float uu = vis->u * uvscale;
	float vv = vis->v * uvscale;
	float uvdist = r_uvdist(rp, uu, vv);
	float err = r_vis_err(rp, vis);
	if(uvdist >= rp->uvmin && uvdist <= rp->uvmax &&
	   err >= rp->errmin && err <= rp->errmax) {
	  cpgpt(1, &uvdist, &err, datsym);  /* Plot point */
	};
      };
    };
  };
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Private function of uvradplt to plot the model visbility amplitudes of
 * a specified baseline versus UV radius, using a specified color.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 *  sub    Subarray *  The subarray to plot the model from.
 *  base        int    The baseline to plot.
 *  modcol      int    The PGPLOT color index to use.
 *  uvmin     float    The minimum UV radius to be displayed (wavelengths)
 *  uvmax     float    The maximum UV radius to be displayed (wavelengths)
 * Output:
 *  return      int    0 - OK.
 */
static int r_modpl(R_par *rp, Subarray *sub, int base, int modcol)
{
  Integration *integ; /* The descriptor of the integration being plotted */
  float uvscale;      /* Factor to scale UV corrds by to get wavelengths */
  int modsym;         /* Plot marker to be used */
  int ut;             /* The ut being looked at */
/*
 * If there is no model to be plotted, or model plotting has not been
 * requested, do nothing.
 */
  if(!rp->ob->hasmod || !rp->domod || rp->dodiff)
    return 0;
/*
 * Get the conversion factor from UVW in light-seconds to wavelengths.
 */
  uvscale = rp->ob->stream.uvscale;
/*
 * Determine which marker symbol to use.
 */
  modsym = rp->dobig ? bigsym : dotsym;
/*
 * Plot the model.
 */
  cpgbbuf();
  cpgsci(modcol);
/*
 * Start with the amplitude plot.
 */
  if(rp->doamp && r_ampwin(rp)==0) {
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
 * Compare the UV radius and visibility amplitude with the cut-off limits.
 */
	float uu = vis->u * uvscale;
	float vv = vis->v * uvscale;
	float uvdist = r_uvdist(rp, uu, vv);
	float amp = vis->modamp;
	if(uvdist >= rp->uvmin && uvdist <= rp->uvmax &&
	   amp >= rp->ampmin && amp <= rp->ampmax) {
	  cpgpt(1, &uvdist, &amp, modsym);  /* Plot point */
	};
      };
    };
  };
/*
 * Now do the phase plot.
 */
  if(rp->dophs && r_phswin(rp)==0) {
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
 * Compare the UV radius and visibility phase with the cut-off limits.
 */
	float uu = vis->u * uvscale;
	float vv = vis->v * uvscale;
	float uvdist = r_uvdist(rp, uu, vv);
	float phs = r_mod_phs(rp, vis, uu, vv);
	if(uvdist >= rp->uvmin && uvdist <= rp->uvmax &&
	   phs >= rp->phsmin && phs <= rp->phsmax) {
	  cpgpt(1, &uvdist, &phs, modsym);  /* Plot point */
	};
      };
    };
  };
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Given the indexes of a previously highlighted telescope and a new
 * telescope to highlight, re-plot the previous telescope in the normal
 * color and plot the baselines of the new telescope in highlighted color.
 *
 * Input:
 *  rp           R_par *  Plot-parameter block.
 *  oper        Telop     R_ALLNEW  - Instate the new telescope spec from *init.
 *                        R_NXTSUB  - Highlight telescope of next sub-array.
 *                        R_NXTTEL  - Highlight the next telescope.
 *  forward       int     0 - Highlight the next telescope in order of
 *                            decreasing telescope index.
 *                        1 - Highlight the next telescope in order of
 *                            increasing telescope index.
 *  init      Telspec *  The new telescope specification to be used when
 *                       oper == R_ALLNEW.
 * Output:
 *  return         int    0 - OK.
 */
static int r_newtel(R_par *rp, Telop oper, int forward, Telspec *init)
{
  Telspec ts;      /* The new telescope specification */
  Observation *ob; /* The descriptor of the observation */
  int oldtel;      /* THe index of the previous reference telescope */
  int oldsub;      /* THe index of the previous reference telescope sub-array */
  int newtel;      /* The index of the new reference telescope */
  int newsub;      /* The index of the new reference telescope sub-array */
  int cif;         /* The index of the IF being processed */
/*
 * Get the descriptor of the observation.
 */
  ob = rp->ob;
/*
 * Handle the specified change in reference telescope.
 */
  switch(oper) {
  case R_ALLNEW:
    ts = *init;
    if(next_tel(ob, FIND_FIRST, 1, 0, 0, 1, &ts))
      return 0;
    break;
  case R_NXTSUB:
  case R_NXTTEL:
/*
 * If highlighting is currently disabled, re-enable it if the search
 * direction is forward. Otherwise ignore the request.
 */
    if(!rp->highlight) {
      if(forward)
	ts = rp->init;
      else
	return 0;
/*
 * Supplant the currently highlighted telescope.
 */
    } else {
      int iret=1;
      ts = rp->ts;
/*
 * Locate the new telescope. Note the fallthrough between cases.
 */
      switch(oper) {
      case R_NXTTEL:
	iret = next_tel(ob, SKIP_TA, forward, 0, 0, 0, &ts);
      case R_NXTSUB:
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
    lprintf(stderr, "r_newtel: Unrecognised opcode.\n");
    return 1;
  };
/*
 * Get the details of the currently hightlighted station and the
 * new station.
 * Note that an empty spec (nfix==0) denotes no highlighting.
 */
  oldtel = rp->highlight ? rp->ts.ta : -1;
  oldsub = rp->highlight ? rp->ts.isub : -1;
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
	    r_basepl(rp, sub, base, isold ? datcol : altcol);
	};
      };
    };
  };
/*
 * Erase the old station title.
 */
  if(rp->highlight)
    r_namplt(rp, 1);
/*
 * Record the new highlighted telescope.
 */
  rp->highlight = ts.nfix != 0;
  rp->ts = ts;
/*
 * Plot the new station title.
 */
  if(rp->highlight)
    r_namplt(rp, 0);
/*
 * Terminate buffering.
 */
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Erase obsolete telescope title and plot new one.
 *
 * Input:
 *  rp      R_par *  The plot descriptor.
 *  erase     int    If true, erase the telescope name.
 */
static void r_namplt(R_par *rp, int erase)
{
  char title[81];  /* Temporary buffer to compose label in */
  int reftel;      /* The index of the telescope to be named */
  int refsub;      /* The index of the sub-array to be named */
/*
 * Get local copies of the indexes of the reference telescope and sub-array.
 */
  reftel = rp->ts.ta;
  refsub = rp->ts.isub;
/*
 * Set a viewport enclosing all plots.
 */
  cpgsvp(rp->vxa, rp->vxb, rp->vya, rp->vyb);
/*
 * Telescope to name?
 */
  if(rp->highlight) {
    cpgsci(erase?0:1);
    sprintf(title, "%d:%s", refsub+1, rp->ob->sub[refsub].tel[reftel].name);
    cpgmtxt("T", 1.0f, 1.0f, 1.0f, title);
    cpgsci(1);
  };
  return;
}

/*.......................................................................
 * Search for the data-point nearest the cursor.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 *  xpos      float    The X-world coordinate of the cursor.
 *  ypos      float    The Y-world coordinate of the cursor.
 *  plot     RpType    The sub-plot to which ypos refers.
 * Ourput:
 *  return   Bestvis   Container, containing the details of the
 *                     nearest point. Only use the results if best->found
 *                     is true.
 */
static Bestvis r_findpt(R_par *rp, float xpos, float ypos, RpType plot)
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
 * Set a window around the sub-plot in which the cursor was pressed.
 */
  switch(plot) {
  case AMP_PLOT:
    r_ampwin(rp);
    break;
  case PHS_PLOT:
    r_phswin(rp);
    break;
  case ERR_PLOT:
    r_errwin(rp);
    break;
  };
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
  ob = rp->ob;
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
 * Get the new UV radius.
 */
	    float uu = vis->u * uvscale;
	    float vv = vis->v * uvscale;
	    float uvrad = r_uvdist(rp, uu, vv);
	    float value = plot==AMP_PLOT ? r_vis_amp(rp,vis) :
	                 (plot==PHS_PLOT ? r_vis_phs(rp, vis, uu, vv) :
			 (plot==ERR_PLOT ? r_vis_err(rp,vis) : 0.0f));
/*
 * Get the physical X offset between the cursor and the data point.
 */
	    float xoff = xtomm * (uvrad - xpos);
/*
 * Get the physical Y offset between the cursor and the data point.
 */
	    float yoff = ytomm * (value - ypos);
/*
 * Compute the physical radius (squared) of the cursor from the point.
 */
	    float newrsq = xoff*xoff + yoff*yoff;
/*
 * Is the point within the plot bounds?
 */
	    if(uvrad >= rp->uvmin && uvrad <= rp->uvmax &&
	       (plot==AMP_PLOT ? (value >= rp->ampmin && value <= rp->ampmax) :
	       (plot==PHS_PLOT ? (value >= rp->phsmin && value <= rp->phsmax) :
	       (plot==ERR_PLOT ? (value >= rp->errmin && value <= rp->errmax) :
		0.0f)))) {
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
/*
 * Nothing found?
 */
  if(!best.found)
    lprintf(stderr, "r_findpt: No data in range.\n");
/*
 * Return details of the nearest visibility.
 */
  return best;
}

/*.......................................................................
 * Flag a given visibility and erase it from the display.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
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
static int r_editpt(R_par *rp, int cif, int isub, int base, int ut, int flag)
{
  Observation *ob; /* The descriptor of the observation being edited */
  Subarray *sub;   /* The descriptor of the sub-array being edited */
  Visibility *vis; /* The visibility being editted. */
  int flagged;     /* True if the visibility was originally flagged */
/*
 * Sanity check the arguments.
 */
  if(rp==NULL) {
    lprintf(stderr, "r_editpt: NULL plot descriptor intercepted\n");
    return 1;
  };
/*
 * Get the descriptor of the observation.
 */
  ob = rp->ob;
/*
 * Check the given sub-array index.
 */
  if(isub<0 || isub>=ob->nsub) {
    lprintf(stderr, "r_editpt: Out of range sub-array index.\n");
    return 1;
  };
/*
 * Check the given IF index.
 */
  if(cif<0 || cif>=ob->nif) {
    lprintf(stderr, "r_editpt: Out of range IF index.\n");
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
    lprintf(stderr, "r_editpt: Out of range baseline index.\n");
    return 1;
  };
  if(ut<0 || ut>=sub->ntime) {
    lprintf(stderr, "r_editpt: Out of range integration index.\n");
    return 1;
  };
/*
 * Get the required IF.
 */
  if(getIF(ob, cif))
    return 1;
/*
 * Mark the data as modified.
 */
  rp->modified = 1;
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
    float uu = vis->u * uvscale;        /* The U coord of the vsibility */
    float vv = vis->v * uvscale;        /* The V coord of the vsibility */
    float uvrad = r_uvdist(rp, uu, vv);
/*
 * Edit the data.
 */
    ed_integ(ob, sub, ut, cif, flag, 1, 0, rp->ch_ed, 1, base);
/*
 * Store the entry color.
 */
    cpgqci(&oldcol);
/*
 * Erase flagged or draw unflagged point.
 */
    cpgsci(flag ? 0 : (isub==rp->ts.isub && (sub->base[base].tel_a==rp->ts.ta || sub->base[base].tel_b==rp->ts.ta) ? altcol : datcol));
    if(rp->doamp) {
      float amp = r_vis_amp(rp, vis);
      r_ampwin(rp);
      cpgpt(1, &uvrad, &amp, rp->dobig ? bigsym : dotsym);
    };
    if(rp->dophs) {
      float phs = r_vis_phs(rp, vis, uu, vv);
      r_phswin(rp);
      cpgpt(1, &uvrad, &phs, rp->dobig ? bigsym : dotsym);
    };
    if(rp->doerr) {
      float err = r_vis_err(rp, vis);
      r_errwin(rp);
      cpgpt(1, &uvrad, &err, rp->dobig ? bigsym : dotsym);
    };
/*
 * Restore entry color.
 */
    cpgsci(oldcol);
  };
  return 0;
}

/*.......................................................................
 * Allow user selection of a new amplitude or phase display range.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int r_zoom(R_par *rp)
{
  RpType plot=AMP_PLOT; /* The plot in which the first selection was made */
  int npts=0;           /* The number of vertices selected with the cursor */
  float valmin=0.0f;    /* The selected min amplitude/phase */
  float valmax=0.0f;    /* The selected max amplitude/phase */
/*
 * Tell the user what to do.
 */
  printf("Select two y-axis limits with the cursor (Press '%c' for help).\n",
	 KEY_HELP);
/*
 * Get the two amplitude or phase selections.
 */
  while(npts < 2) {
/*
 * Read the cursor.
 */
    if(r_cursor(rp, npts==0 ? B_YVAL:B_YRNG, plot, rp->uvmin, valmin, zoomcol))
      return 1;
/*
 * Act on the key selected.
 */
    switch(rp->kp.key) {
    case KEY_ZOOM:  /* Select full available amplitude range */
      return r_setrange(rp, rp->fixuvr,rp->uvmin,rp->uvmax, 0,0.0f,0.0f,
			0,0.0f,0.0f, 0,0.0f,0.0f) || r_redisp(rp, 1);
      break;
    case KEY_CUR:  /* Select amplitude */
      if(npts==0) {
	valmin = valmax = rp->kp.value;
	plot = rp->kp.plot;
      } else {
	if(plot != rp->kp.plot) {
	  lprintf(stderr,"zoom: Area spans two windows - selection aborted.\n");
	  return 0;
	};
	if(rp->kp.value > valmin)
	  valmax = rp->kp.value;
	else
	  valmin = rp->kp.value;
      };
      npts++;
      break;
    case KEY_CAN:  /* Cancel selection */
      printf("Display range selection cancelled.\n");
      return 0;
      break;
    default:
      printf("You are currently selecting a new %s display range.\n",
	     npts==0 ? "y-axis" :
	     (plot==AMP_PLOT ? "amplitude" :
	     (plot==PHS_PLOT ? "phase" :
	     (plot==ERR_PLOT ? "error" : "unknown"))));
      printf(" %c - Select the %s limit of the range.\n",
	     KEY_CUR, npts==0 ? "lower":"upper");
      printf(" %c - Select the full available ranges.\n", KEY_ZOOM);
      printf(" %c - Abort the display range selection.\n", KEY_CAN);
      break;
    };
  };
/*
 * Instate the new range.
 */
  switch(plot) {
  case AMP_PLOT:
    if(r_setrange(rp, rp->fixuvr, rp->uvmin, rp->uvmax,
		  1, valmin, valmax,
		  rp->fixphs, rp->phsmin, rp->phsmax,
		  rp->fixerr, rp->errmin, rp->errmax))
      return 1;
    break;
  case PHS_PLOT:
    if(r_setrange(rp, rp->fixuvr, rp->uvmin, rp->uvmax,
		  rp->fixamp, rp->ampmin, rp->ampmax,
		  1, valmin, valmax,
		  rp->fixerr, rp->errmin, rp->errmax))
      return 1;
    break;
  case ERR_PLOT:
    if(r_setrange(rp, rp->fixuvr, rp->uvmin, rp->uvmax,
		  rp->fixamp, rp->ampmin, rp->ampmax,
		  rp->fixphs, rp->phsmin, rp->phsmax,
		  1, valmin, valmax))
      return 1;
    break;
  };
/*
 * Replot the range.
 */
  return r_redisp(rp, 1);
}

/*.......................................................................
 * Handle interactive selection of a new UV distance display range.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int r_newuvr(R_par *rp)
{
  int npts=0;       /* The number of vertices selected with the cursor */
  float uvmin=0.0f; /* The selected UV min radius */
  float uvmax=0.0f; /* The selected UV max radius */
/*
 * Tell the user what to do.
 */
  printf("Select two UV distances with the cursor (Press '%c' for help).\n",
	 KEY_HELP);
/*
 * Get two UV distance selections.
 */
  while(npts < 2) {
/*
 * Read the cursor.
 */
    if(r_cursor(rp, npts==0 ? B_XVAL:B_XRNG, 0, uvmin, rp->ampmin, zoomcol))
      return 1;
/*
 * Act on the key selected.
 */
    switch(rp->kp.key) {
    case KEY_UVR:  /* Select full available UV range */
      return r_setrange(rp, 0,0.0,0.0, rp->fixamp, rp->ampmin, rp->ampmax,
			rp->fixphs, rp->phsmin, rp->phsmax,
			rp->fixerr, rp->errmin, rp->errmax) || r_redisp(rp, 1);
      break;
    case KEY_CUR:  /* Select UV radius */
      if(npts==0) {
	uvmin = uvmax = rp->kp.uvdist;
      } else {
	if(rp->kp.uvdist > uvmin)
	  uvmax = rp->kp.uvdist;
	else
	  uvmin = rp->kp.uvdist;
      };
      npts++;
      break;
    case KEY_CAN:  /* Cancel selection */
      printf("UV range selection cancelled.\n");
      return 0;
      break;
    default:
      printf("You are currently selecting a new UV radius display range.\n");
      printf(" %c - Select the %s UV radius of the range.\n",
	     KEY_CUR, npts==0 ? "start":"end");
      printf(" %c - Select the full available UV range.\n",
	     KEY_UVR);
      printf(" %c - Abort the UV range selection.\n", KEY_CAN);
      break;
    };
  };
/*
 * Instate the new range.
 */
  if(r_setrange(rp, 1, uvmin, uvmax, rp->fixamp, rp->ampmin, rp->ampmax,
		rp->fixphs, rp->phsmin, rp->phsmax,
		rp->fixerr, rp->errmin, rp->errmax))
    return 1;
/*
 * Redisplay with the new range.
 */
  return r_redisp(rp, 1);
}

/*.......................................................................
 * Read the cursor and return the key that was entered. Also return the
 * position at which the cursor was pressed, bounded by the axis edges.
 *
 * Input:
 *  rp      R_par *  The plot descriptor.
 *                   rp->kp descriptor will be initialized with the
 *                   position and key selected by the user. On the first
 *                   call, set rp->kp.key=KEY_NONE so that the cursor can
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
 *  plot   RpType    The plot that yref refers to.
 *  xref    float    xref and yref denote the reference position for
 *  yref    float    band-type cursors. They are ignored if mode==B_NORM.
 *  ci        int    The color index to draw band cursors with.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int r_cursor(R_par *rp, Bandmode mode, RpType plot,
		    float xref, float yref, int ci)
{
  static float xpos=0.0f; /* The X-world-coordinate of the cursor */
  static float ypos=0.0f; /* The Y-world-coordinate of the cursor */
  char key;               /* The key that caused cpgcurs() to return */
/*
 * Check arguments.
 */
  if(rp==NULL) {
    lprintf(stderr, "get_curs: NULL plot descriptor intercepted\n");
    return 1;
  };
/*
 * Set a window around the whole plot with NDC world coordinates.
 */
  cpgsvp(rp->vxa, rp->vxb, rp->vya, rp->vyb);
  cpgswin(rp->vxa, rp->vxb, rp->vya, rp->vyb);
/*
 * Position the cursor in the center of the plot?
 */
  if(rp->kp.key==KEY_NONE) {
    xpos = (rp->vxa + rp->vxb) / 2.0f;
    ypos = (rp->vya + rp->vyb) / 2.0f;
  } else {
    xpos = rp->kp.xpos;
    ypos = rp->kp.ypos;
  };
/*
 * Substitute the cross-hair cursor for the normal cursor is requested.
 */
  if(mode == B_NORM && rp->docross)
    mode = B_CROSS;
/*
 * Limit the given X and Y coordinates such that the cursor will appear within
 * the plot.
 */
  if(xpos < rp->vxa)
    xpos = rp->vxa;
  else if(xpos > rp->vxb)
    xpos = rp->vxb;
  if(ypos < rp->vya)
    ypos = rp->vya;
  else if(ypos > rp->vyb)
    ypos = rp->vyb;
/*
 * Convert the cursor reference positions into NDC.
 */
  xref = rp->vxa + (xref - rp->wxa) * (rp->vxb-rp->vxa)/(rp->wxb-rp->wxa);
  switch(plot) {
  case AMP_PLOT:
    yref = rp->vabot + (yref - rp->wyaa) *
      (rp->vatop - rp->vabot) / (rp->wyab - rp->wyaa);
    break;
  case PHS_PLOT:
    yref = rp->vpbot + (yref - rp->wypa) *
      (rp->vptop - rp->vpbot) / (rp->wypb - rp->wypa);
    break;
  case ERR_PLOT:
    yref = rp->vebot + (yref - rp->wyea) *
      (rp->vetop - rp->vebot) / (rp->wyeb - rp->wyea);
    break;
  };
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
  rp->kp.waslow = islower((int)key);
/*
 * Convert the selected key to upper case to simplify pending comparisons.
 */
  if(rp->kp.waslow)
    key = toupper(key);
/*
 * Enforce plot bounds on the returned cursor position.
 */
  if(xpos < rp->vxa)
    xpos = rp->vxa;
  else if(xpos > rp->vxb)
    xpos = rp->vxb;
  if(ypos < rp->vya)
    ypos = rp->vya;
  else if(ypos > rp->vyb)
    ypos = rp->vyb;  
/*
 * Record the results in the Keypos structure.
 */
  rp->kp.xpos = xpos;
  rp->kp.ypos = ypos;
  rp->kp.key = key;
  rp->kp.uvdist = rp->wxa+(rp->wxb-rp->wxa)*(xpos-rp->vxa)/(rp->vxb-rp->vxa);
/*
 * Work out which sub-plot the cursor was in, and what the corresponding
 * world coordinate selected was.
 */
  if(rp->doamp && ypos >= rp->vabot && ypos <= rp->vatop) {
    rp->kp.plot = AMP_PLOT;
    rp->kp.value = rp->wyaa + (rp->wyab - rp->wyaa) * 
      (ypos - rp->vabot) / (rp->vatop - rp->vabot);
  } else if(rp->dophs && ypos >= rp->vpbot && ypos <= rp->vptop) {
    rp->kp.plot = PHS_PLOT;
    rp->kp.value = rp->wypa + (rp->wypb - rp->wypa) * 
      (ypos - rp->vpbot) / (rp->vptop - rp->vpbot);
  } else if(rp->doerr && ypos >= rp->vebot && ypos <= rp->vetop) {
    rp->kp.plot = ERR_PLOT;
    rp->kp.value = rp->wyea + (rp->wyeb - rp->wyea) * 
      (ypos - rp->vebot) / (rp->vetop - rp->vebot);
  } else {
    rp->kp.plot = AMP_PLOT;
    rp->kp.value = 0.0f;
  };
  return 0;
}


/*.......................................................................
 * Flag points within a user selected rectangular box.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int r_edbox(R_par *rp)
{
  Observation *ob;      /* The descriptor of the observation */
  int cif;              /* The index of the IF being processed */
  int npts=0;           /* The number of points so far selected */
  RpType plot=AMP_PLOT; /* The plot that the user first selected */
  float valmin=0.0f;    /* The value of the base of the selected area */
  float valmax=0.0f;    /* The value of the top of the selected area */
  float uvmin=0.0f;     /* The UV dist of the left edge of the selected area */
  float uvmax=0.0f;     /* The UV dist of the right edge of the selected area */
  float uvscale;        /* Factor to scale UV coords by to get wavelengths */
  float uvref=0.0f;     /* The UV distance of the first vertex selected */
  float valref=0.0f;    /* The amp/phase of the first vertex selected */
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
    if(r_cursor(rp, npts==0?B_NORM:B_RECT, plot, uvref, valref, cutcol))
      return 1;
/*
 * Act on the key used to return the cursor.
 */
    switch(rp->kp.key) {
    case KEY_CUR:         /* Select position in plot */
/*
 * First position?
 */
      if(npts==0) {
	plot = rp->kp.plot;
	uvref = uvmin = uvmax = rp->kp.uvdist;
	valref = valmin = valmax = rp->kp.value;
      } else {
	if(plot != rp->kp.plot) {
	  lprintf(stderr, "Select box spans sub-plots - edit cancelled.\n");
	  return 1;
	};
/*
 * Record the new UV distance as the minimum or maximum of the two chosen.
 */
	if(rp->kp.uvdist > uvmin)
	  uvmax = rp->kp.uvdist;
	else
	  uvmin = rp->kp.uvdist;
/*
 * Record the new amplitude as the minimum or maximum of the two chosen.
 */
	if(rp->kp.value > valmin)
	  valmax = rp->kp.value;
	else
	  valmin = rp->kp.value;
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
  ob = rp->ob;
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
    for(isub=0; isub<rp->ob->nsub; isub++,sub++) {
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
 * Get the new UV radius.
 */
	    float uu = vis->u * uvscale;
	    float vv = vis->v * uvscale;
	    float uvrad = r_uvdist(rp, uu, vv);
	    float value = plot==AMP_PLOT ? r_vis_amp(rp,vis) :
	                 (plot==PHS_PLOT ? r_vis_phs(rp, vis, uu, vv) :
			 (plot==ERR_PLOT ? r_vis_err(rp,vis) : 0.0f));
/*
 * Is the point within the selected rectangular area?
 */
	    if(uvrad >= uvmin && uvrad <= uvmax &&
	       value >= valmin && value <= valmax) {
/*
 * Edit the point and erase it from the plot.
 */
	      if(r_editpt(rp, cif, isub, base, ut, 1)) {
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
 * Display the statistics of the visibilities that lie within a box that
 * the user draws.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int r_scalar_stats(R_par *rp)
{
  Observation *ob;      /* The descriptor of the observation */
  int cif;              /* The index of the IF being processed */
  int npts=0;           /* The number of points so far selected */
  int nvis=0;           /* The number of visibilities within the selected */
                        /*  area */
  RpType plot=AMP_PLOT; /* The plot that the user first selected */
  int iter;             /* An iteration count */
  float valmin=0.0f;    /* The value of the base of the selected area */
  float valmax=0.0f;    /* The value of the top of the selected area */
  float uvmin=0.0f;     /* The UV dist of the left edge of the selected area. */
  float uvmax=0.0f;     /* The UV dist of the right edge of the selected area.*/
  float uvscale;        /* Factor to scale UV coords by to get wavelengths */
  float uvref=0.0f;     /* The UV distance of the first vertex selected */
  float valref=0.0f;    /* The amp/phase of the first vertex selected */
  double ampsum=0.0f;   /* The sum of the amp's of the selected visibilities */
  double phssum=0.0f;   /* The sum of the phases of the selected visibilities */
  double mean_amp=0.0f; /* The mean ampitude of the selected visibilities */
  double mean_phs=0.0f; /* The mean phase of the selected visibilities */
  double ampsumsq=0.0;  /* The sum of the squares of the (amp - mean_amp) */
  double phssumsq=0.0;  /* The sum of the squares of the (phase - mean_phase) */
/*
 * Tell the user what to do.
 */
  printf("Set the two opposite corners of the area of interest. Press '%c' for help.\n",
	 KEY_HELP);
/*
 * Have the user select two points within the plot.
 */
  while(npts < 2) {
/*
 * Read the cursor.
 */
    if(r_cursor(rp, npts==0?B_NORM:B_RECT, plot, uvref, valref, statcol))
      return 1;
/*
 * Act on the key used to return the cursor.
 */
    switch(rp->kp.key) {
    case KEY_CUR:         /* Select position in plot */
/*
 * First position?
 */
      if(npts==0) {
	plot = rp->kp.plot;
	uvref = uvmin = uvmax = rp->kp.uvdist;
	valref = valmin = valmax = rp->kp.value;
      } else {
	if(plot != rp->kp.plot) {
	  lprintf(stderr, "Select box spans amp and phase boxes - edit cancelled.\n");
	  return 1;
	};
/*
 * Record the new UV distance as the minimum or maximum of the two chosen.
 */
	if(rp->kp.uvdist > uvmin)
	  uvmax = rp->kp.uvdist;
	else
	  uvmin = rp->kp.uvdist;
/*
 * Record the new amplitude as the minimum or maximum of the two chosen.
 */
	if(rp->kp.value > valmin)
	  valmax = rp->kp.value;
	else
	  valmin = rp->kp.value;
      };
/*
 * Record the successful acquisition of a new corner.
 */
      npts++;
      break;
    case KEY_CAN:  /* Cancel selection */
      printf("Selection cancelled.\n");
      return 0;
      break;
    default:
      printf("You are currently selecting a rectangular area to be described - use keys:\n");
      printf(" %c - Select %s corner of the area with this key.\n",
	     KEY_CUR, npts==0 ? "a" : "the second (opposite)");
      printf(" %c - Abort the selection with this key.\n", KEY_CAN);
      break;
    };
  };
/*
 * Get the descriptor of the observation.
 */
  ob = rp->ob;
/*
 * Loop through the data twice, once to find the mean of the visibilities
 * within the selected area, then again to find the scatter.
 */
  for(iter=0; iter<2; iter++) {
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
      sub = ob->sub;
      for(isub=0; isub<rp->ob->nsub; isub++,sub++) {
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
 * Get the new UV radius.
 */
	      float uu = vis->u * uvscale;
	      float vv = vis->v * uvscale;
	      float uvrad = r_uvdist(rp, uu, vv);
	      float amp = r_vis_amp(rp, vis);
	      float phs = r_vis_phs(rp, vis, uu, vv);
	      float value = plot==AMP_PLOT ? amp :
	                   (plot==PHS_PLOT ? phs :
		   	   (plot==ERR_PLOT ? r_vis_err(rp,vis) : 0.0f));
/*
 * Is the point within the selected rectangular area?
 */
	      if(uvrad >= uvmin && uvrad <= uvmax &&
		 value >= valmin && value <= valmax) {
		if(iter==0) {
		  ampsum += amp;
		  phssum += phs;
		  nvis++;
		} else {
		  float ampdiff = amp - mean_amp;
		  float phsdiff = phs - mean_phs;
		  ampsumsq += ampdiff * ampdiff;
		  phssumsq += phsdiff * phsdiff;
		};
	      };
	    };
	  };
	};
      };
    };
/*
 * Stop now if no visibilities were found.
 */
    if(nvis < 1) {
      lprintf(stdout, "No visibilities lie within the selected area.\n");
      return 0;
    };
/*
 * At the end of the first iteration, compute the mean.
 */
    if(iter==0) {
      mean_amp = ampsum / nvis;
      mean_phs = phssum / nvis;
    };
  };
/*
 * Report the results.
 */
  lprintf(stdout,
	  "\nThe statistics of the %d visibilities within the box are:\n",
	  nvis);
  lprintf(stdout, " Amp mean=%f +/- %f  RMS scatter=%f\n",
	  mean_amp, sqrt(ampsumsq)/nvis, sqrt(ampsumsq/nvis));
  lprintf(stdout, " Phs mean=%f +/- %f  RMS scatter=%f  (degrees)\n",
	  mean_phs * rtod, sqrt(phssumsq)/nvis * rtod,
	  sqrt(phssumsq/nvis) * rtod);
  return 0;
}

/*.......................................................................
 * Display the real/imaginary statistics of the visibilities that lie
 * within a box that the user draws.
 *
 * Input:
 *  rp        R_par *  Plot-parameter block.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int r_vector_stats(R_par *rp)
{
  Observation *ob;      /* The descriptor of the observation */
  int cif;              /* The index of the IF being processed */
  int npts=0;           /* The number of points so far selected */
  int nvis=0;           /* The number of visibilities within the selected */
                        /*  area */
  RpType plot=AMP_PLOT; /* The plot that the user first selected */
  int iter;             /* An iteration count */
  float valmin=0.0f;    /* The value of the base of the selected area */
  float valmax=0.0f;    /* The value of the top of the selected area */
  float uvmin=0.0f;     /* The UV dist of the left edge of the selected area. */
  float uvmax=0.0f;     /* The UV dist of the right edge of the selected area.*/
  float uvscale;        /* Factor to scale UV coords by to get wavelengths */
  float uvref=0.0f;     /* The UV distance of the first vertex selected */
  float valref=0.0f;    /* The amp/phase of the first vertex selected */
  float resum=0.0f;     /* The sum of the real parts of the selected */
                        /*  visibilities */
  float imsum=0.0f;     /* The sum of the imaginary parts of the selected */
                        /*  visibilities */
  float mean_re=0.0f;   /* The mean real part of the selected visibilities */
  float mean_im=0.0f;   /* The mean imag part of the selected visibilities */
  float resumsq=0.0;    /* The sum of the squares of the (re - mean_re) */
  float imsumsq=0.0;    /* The sum of the squares of the (im - mean_im) */
/*
 * Tell the user what to do.
 */
  printf("Set the two opposite corners of the area of interest. Press '%c' for help.\n",
	 KEY_HELP);
/*
 * Have the user select two points within the plot.
 */
  while(npts < 2) {
/*
 * Read the cursor.
 */
    if(r_cursor(rp, npts==0?B_NORM:B_RECT, plot, uvref, valref, statcol))
      return 1;
/*
 * Act on the key used to return the cursor.
 */
    switch(rp->kp.key) {
    case KEY_CUR:         /* Select position in plot */
/*
 * First position?
 */
      if(npts==0) {
	plot = rp->kp.plot;
	uvref = uvmin = uvmax = rp->kp.uvdist;
	valref = valmin = valmax = rp->kp.value;
      } else {
	if(plot != rp->kp.plot) {
	  lprintf(stderr, "Select box spans amp and phase boxes - edit cancelled.\n");
	  return 1;
	};
/*
 * Record the new UV distance as the minimum or maximum of the two chosen.
 */
	if(rp->kp.uvdist > uvmin)
	  uvmax = rp->kp.uvdist;
	else
	  uvmin = rp->kp.uvdist;
/*
 * Record the new real part as the minimum or maximum of the two chosen.
 */
	if(rp->kp.value > valmin)
	  valmax = rp->kp.value;
	else
	  valmin = rp->kp.value;
      };
/*
 * Record the successful acquisition of a new corner.
 */
      npts++;
      break;
    case KEY_CAN:  /* Cancel selection */
      printf("Selection cancelled.\n");
      return 0;
      break;
    default:
      printf("You are currently selecting a rectangular area to be described - use keys:\n");
      printf(" %c - Select %s corner of the area with this key.\n",
	     KEY_CUR, npts==0 ? "a" : "the second (opposite)");
      printf(" %c - Abort the selection with this key.\n", KEY_CAN);
      break;
    };
  };
/*
 * Get the descriptor of the observation.
 */
  ob = rp->ob;
/*
 * Loop through the data twice, once to find the mean of the visibilities
 * within the selected area, then again to find the scatter.
 */
  for(iter=0; iter<2; iter++) {
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
      sub = ob->sub;
      for(isub=0; isub<rp->ob->nsub; isub++,sub++) {
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
 * Get the new UV radius.
 */
	      float uu = vis->u * uvscale;
	      float vv = vis->v * uvscale;
	      float uvrad = r_uvdist(rp, uu, vv);
	      float amp = r_vis_amp(rp, vis);
	      float phs = r_vis_phs(rp, vis, uu, vv);
	      float value = plot==AMP_PLOT ? amp :
	                   (plot==PHS_PLOT ? phs :
		   	   (plot==ERR_PLOT ? r_vis_err(rp,vis) : 0.0f));
/*
 * Is the point within the selected rectangular area?
 */
	      if(uvrad >= uvmin && uvrad <= uvmax &&
		 value >= valmin && value <= valmax) {
		float re = amp * cos(phs);
		float im = amp * sin(phs);
		if(iter==0) {
		  resum += re;
		  imsum += im;
		  nvis++;
		} else {
		  float rediff = re - mean_re;
		  float imdiff = im - mean_im;
		  resumsq += rediff * rediff;
		  imsumsq += imdiff * imdiff;
		};
	      };
	    };
	  };
	};
      };
    };
/*
 * Stop now if no visibilities were found.
 */
    if(nvis < 1) {
      lprintf(stdout, "No visibilities lie within the selected area.\n");
      return 0;
    };
/*
 * At the end of the first iteration, compute the mean.
 */
    if(iter==0) {
      mean_re = resum / nvis;
      mean_im = imsum / nvis;
    };
  };
/*
 * Report the results.
 */
  lprintf(stdout,
	  "\nThe statistics of the %d visibilities within the box are:\n",
	  nvis);
  lprintf(stdout, " Real mean=%f +/- %f  RMS scatter=%f\n",
	  mean_re, sqrt(resumsq)/nvis, sqrt(resumsq/nvis));
  lprintf(stdout, " Imag mean=%f +/- %f  RMS scatter=%f\n",
	  mean_im, sqrt(imsumsq)/nvis, sqrt(imsumsq/nvis));
  return 0;
}

/*.......................................................................
 * Change the plotted amplitude phase and UV radius ranges.
 * Note that r_redisp() must be called after this function.
 *
 * Input:
 *  rp     R_par *  The plot descriptor.
 *  fixuvr   int    If fixuvr != 0 and uvmin and uvmax describe a
 *                  valid UV radius range (after positivity has been
 *                  enforced and the values sorted such that
 *                  uvmin<=uvmax), then the UV radius range will be
 *                  set by uvmin and uvmax. Otherwise it will be
 *                  determined from the data.
 *  uvmin  float    Minimum UV radius (wavelengths).
 *  uvmax  float    Maximum UV radius (wavelengths).
 *  fixamp   int    If fixamp != 0 and ampmin and ampmax describe a
 *                  valid amplitude range (after positivity has been
 *                  enforced and the values sorted such that
 *                  ampmin<=ampmax), then the amplitude range will be
 *                  set by ampmin and ampmax. Otherwise it will be
 *                  determined from the data.
 *  ampmin  float   Minimum amplitude.
 *  ampmax  float   Maximum amplitude.
 *  fixphs   int    If fixphs != 0 and phsmin and phsmax describe a
 *                  valid phase range (after positivity has been
 *                  enforced and the values sorted such that
 *                  phsmin<=phsmax), then the phase range will be
 *                  set by phsmin and phsmax. Otherwise it will be
 *                  determined from the data.
 *  phsmin  float   Minimum phase.
 *  phsmax  float   Maximum phase.
 *  fixerr   int    If fixerr != 0 and errmin and errmax describe a
 *                  valid error range (after positivity has been
 *                  enforced and the values sorted such that
 *                  errmin<=errmax), then the error range will be
 *                  set by errmin and errmax. Otherwise it will be
 *                  determined from the data.
 *  errmin  float   Minimum error.
 *  errmax  float   Maximum error.
 * Output:
 *  return    int   0 - OK.
 *                  1 - Error.
 */
static int r_setrange(R_par *rp, int fixuvr, float uvmin, float uvmax,
		      int fixamp, float ampmin, float ampmax,
		      int fixphs, float phsmin, float phsmax,
		      int fixerr, float errmin, float errmax)
{
/*
 * Check fixed amplitude ranges where given.
 */
  if(fixamp) {
/*
 * The amplitude range should not encompass unphysical negative amplitudes.
 */
    if(ampmin < 0.0f)
      ampmin = 0.0f;
    if(ampmax < 0.0f)
      ampmin = 0.0f;
/*
 * Ensure that ampmin < ampmax.
 */
    if(ampmin > ampmax) {float newmin=ampmax; ampmax=ampmin; ampmin=newmin;};
/*
 * If the range is too small, revert to autoscaling.
 */
    if(ampmin==ampmax)
      fixamp = 0;
  };
/*
 * Check fixed phase ranges where given.
 */
  if(fixphs) {
/*
 * Ensure that phsmin < phsmax.
 */
    if(phsmin > phsmax) {float newmin=phsmax; phsmax=phsmin; phsmin=newmin;};
/*
 * Limit the phase range to +/-pi.
 */
    if(phsmin < -pi)
      phsmin = -pi;
    if(phsmin > pi)
      phsmin = pi;
    if(phsmax < -pi)
      phsmax = -pi;
    if(phsmax > pi)
      phsmax = pi;
/*
 * If the range is too small, revert to autoscaling.
 */
    if(phsmax - phsmin < 1.0e-5f * dtor)
      fixphs = 0;
  };
/*
 * Check fixed error ranges where given.
 */
  if(fixerr) {
/*
 * The error range should not encompass unphysical negative errors.
 */
    if(errmin < 0.0f)
      errmin = 0.0f;
    if(errmax < 0.0f)
      errmin = 0.0f;
/*
 * Ensure that errmin < errmax.
 */
    if(errmin > errmax) {float newmin=errmax; errmax=errmin; errmin=newmin;};
/*
 * If the range is too small, revert to autoscaling.
 */
    if(errmax - errmin < 1.0e-5f * dtor)
      fixerr = 0;
  };
/*
 * Check fixed UV radius ranges where given.
 */
  if(fixuvr) {
/*
 * The UV radius range should not encompass unphysical negative UV radii.
 */
    if(uvmin < 0.0f)
      uvmin = 0.0f;
    if(uvmax < 0.0f)
      uvmin = 0.0f;
/*
 * Ensure that uvmin < uvmax.
 */
    if(uvmin > uvmax) {float newmin=uvmax; uvmax=uvmin; uvmin=newmin;};
/*
 * If the range is too small, revert to autoscaling.
 */
    if(uvmin==uvmax)
      fixuvr = 0;
  };
/*
 * Record the requested ranges in the plot descriptor.
 */
  rp->uvmin = fixuvr ? uvmin : 0.0f;
  rp->uvmax = fixuvr ? uvmax : 0.0f;
  rp->ampmin = fixamp ? ampmin : 0.0f;
  rp->ampmax = fixamp ? ampmax : 0.0f;
  rp->phsmin = fixphs ? phsmin : 0.0f;
  rp->phsmax = fixphs ? phsmax : 0.0f;
  rp->errmin = fixerr ? errmin : 0.0f;
  rp->errmax = fixerr ? errmax : 0.0f;
/*
 * Record the autoscale status of the ranges.
 */
  rp->fixuvr = fixuvr;
  rp->fixamp = fixamp;
  rp->fixphs = fixphs;
  rp->fixerr = fixerr;
/*
 * Re-display the plot.
 */
  return 0;
}

/*.......................................................................
 * Private function of r_axes(), used to update the amplitude, phase and UV
 * radius plotting ranges from the values previously set by
 * r_setrange().
 *
 * Input:
 *  rp    R_par *   The plot descriptor.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Error.
 */
static int r_getrange(R_par *rp)
{
/*
 * The amplitude, UV and error ranges are autoscaled, so if any of these
 * quantities is not fixed call uvrange() to determine the data range.
 */
  if(!rp->fixuvr || !rp->fixamp || !rp->fixerr) {
/*
 * Get the UV range of visibilities in all IFs.
 */
    UVrange *uvr = uvrange(rp->ob, 1, rp->dodiff, 0.0f ,0.0f);
    if(uvr==NULL)
      return 1;
/*
 * Use the determined ranges to set the un-specified parameters.
 */
    if(!rp->fixuvr) {
      rp->uvmin = rp->doproj ? 0.0f : uvr->uvrmin;
      rp->uvmax = uvr->uvrmax;
    };
    if(!rp->fixamp) {
      rp->ampmin = 0.0f;
      rp->ampmax = uvr->ampmax;
    };
    if(!rp->fixerr) {
      rp->errmin = 0.0f;
      rp->errmax = uvr->wtmin!=0.0f ? 1.0f/sqrt(fabs(uvr->wtmin)) : 0.0f;
    };
  };
/*
 * The phase range if not fixed is not autoranged, but set to -180 to
 * 180 degrees.
 */
  if(!rp->fixphs) {
    rp->phsmin = -pi;
    rp->phsmax = +pi;
  };
  return 0;
}

/*.......................................................................
 * Re-plot the mode line to reflect changes in edit mode.
 *
 * Input:
 *  rp        R_par *  The plot descriptor.
 *  ch_ed       int    Select channel based editing if true.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int r_newmode(R_par *rp, int ch_ed)
{
/*
 * Buffer until the new text has been plotted.
 */
  cpgbbuf();
/*
 * Erase the existing mode line.
 */
  r_mlab(rp, 1);
/*
 * Install the new editing modes.
 */
  rp->ch_ed = ch_ed;
/*
 * Draw the new mode line.
 */
  r_mlab(rp, 0); /* Plot new mode line */
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
 *  rp       R_par *  The plot descriptor.
 *  erase       int    If true, erase existing mode label.
 * Output:
 *  return      int    0 - OK.
 */
static int r_mlab(R_par *rp, int erase)
{
  Observation *ob;  /* The descriptor of the observation being plotted */
  int oldcol;       /* Temporary storage for entry color index */
  char label[81];   /* Temporary work string to compose mode label in */
/*
 * Get the descriptor of the observation.
 */
  ob = rp->ob;
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
  sprintf(label, "Edit %s channels.", rp->ch_ed ? "selected" : "all");
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
 * Wrap the given projection angle into the range -pi to pi.
 *
 * Input:
 *  phi    double  The projection angle (radians).
 * Output:
 *  return double  'phi' wrapped into the range -pi to pi.
 */
static double wrapphi(double phi)
{
  phi = fmod(phi, twopi);
  if(phi > pi)
    phi -= twopi;
  else if(phi < -pi)
    phi += twopi;
  return phi;
}

/*.......................................................................
 * Change the angle of projection if in projection mode. This initializes
 * the contents of the rp->proj structure and optionally redisplays the
 * plot.
 *
 * Input:
 *  rp       R_par *  The plot descriptor.
 *  phi     double    The new projection angle.
 *  update     int    If true, redisplay the plot with the new projection
 *                    angle.
 * Output:
 *  return     int    0 - OK.
 */
static int r_newphi(R_par *rp, double phi, int update)
{
/*
 * Install the new value of phi.
 */
  rp->proj.phi = wrapphi(phi);
  rp->proj.sinphi = sin(rp->proj.phi);
  rp->proj.cosphi = cos(rp->proj.phi);
/*
 * Re-display the plot if requested.
 */
  return update && r_redisp(rp, 0);
}

/*.......................................................................
 * Determine the Normalized-Device-Coordinate (NDC) dimensions of the
 * viewports enclosing the plots. The result depends on the values of
 * rp->doamp and rp->dophs. On output rp->vxa,rp->vxb,rp->vya,vy->vya
 * will be assigned the limits of the rectangular area enclosing both the
 * amplitude and phase plots, and rp->vymid will be assigned the coordinate
 * of the line separating the amplitude and phase plots; under the
 * convention that the phase plot goes from rp->vya to rp->vymid and the
 * amplitude plot from rp->vymid to rp->vyb.
 *
 * Input:
 *  rp       R_par *  The plot descriptor.
 * Output:
 *  return      int    0 - OK.
 */
static int r_vpwin(R_par *rp)
{
  float vxa,vxb;  /* X viewport limits enclosing both plots. */
  float vya,vyb;  /* Y viewport limits enclosing both plots. */
  float dy;       /* The height of each sub-plot */
  float vtop;     /* The top of the latest sub-plot */
  int nplot;      /* The number of sub-plots */
/*
 * Get the standard viewport.
 */
  cpgsch(1.0f);
  cpgvstd();
  cpgqvp(0,&vxa, &vxb, &vya, &vyb);
/*
 * Store this in the plot descriptor.
 */
  rp->vxa = vxa;
  rp->vxb = vxb;
  rp->vya = vya;
  rp->vyb = vyb;
/*
 * Count the number of plots.
 */
  nplot = (rp->doamp!=0) + (rp->dophs!=0) + (rp->doerr!=0);
  if(nplot==0) {
    lprintf(stderr, "r_vpwin: No sub-plots selected.\n");
    return 1;
  };
/*
 * Compute the height of each sub-plot.
 */
  dy = (rp->vyb - rp->vya) / nplot;
/*
 * Work out the viewport Y coordinates of the top and bottom of
 * each sub-plot.
 */
  vtop = rp->vya;
  if(rp->doerr) {
    rp->vebot = vtop;
    rp->vetop = vtop += dy;
  } else {
    rp->vebot = rp->vetop = 0.0f;
  };
  if(rp->dophs) {
    rp->vpbot = vtop;
    rp->vptop = vtop += dy;
  } else {
    rp->vpbot = rp->vptop = 0.0f;
  };
  if(rp->doamp) {
    rp->vabot = vtop;
    rp->vatop = vtop += dy;
  } else {
    rp->vabot = rp->vatop = 0.0f;
  };
  return 0;
}

/*.......................................................................
 * Set the viewport and world coordinates appropriately for the amplitude
 * sub-plot.
 *
 * Input:
 *  rp       R_par *  The plot descriptor.
 * Output:
 *  return     int    0 - OK.
 */
static int r_ampwin(R_par *rp)
{
  cpgsvp(rp->vxa, rp->vxb, rp->vabot, rp->vatop);
  cpgswin(rp->wxa, rp->wxb, rp->wyaa, rp->wyab);
  return 0;
}

/*.......................................................................
 * Set the viewport and world coordinates appropriately for the phase
 * sub-plot.
 *
 * Input:
 *  rp       R_par *  The plot descriptor.
 * Output:
 *  return     int    0 - OK.
 */
static int r_phswin(R_par *rp)
{
  cpgsvp(rp->vxa, rp->vxb, rp->vpbot, rp->vptop);
  cpgswin(rp->wxa, rp->wxb, rp->wypa, rp->wypb);
  return 0;
}

/*.......................................................................
 * Set the viewport and world coordinates appropriately for the error
 * sub-plot.
 *
 * Input:
 *  rp       R_par *  The plot descriptor.
 * Output:
 *  return     int    0 - OK.
 */
static int r_errwin(R_par *rp)
{
  cpgsvp(rp->vxa, rp->vxb, rp->vebot, rp->vetop);
  cpgswin(rp->wxa, rp->wxb, rp->wyea, rp->wyeb);
  return 0;
}

/*.......................................................................
 * Return the required projected or unprojected UV radial distance
 * corresponding to give U and V coordinates.
 *
 * Input:
 *  rp     R_par *  The plot descriptor.
 *  u,v    float    The U and V coordinates.
 * Output:
 *  return float    The UV distance in the same units as 'u' and 'v'.
 *                  If rp->doproj is true, then the absolute distance
 *                  projected onto angle rp->proj.phi will be returned.
 *                  Otherwise the radial distance of u,v from the origin
 *                  will be returned.
 */
static float r_uvdist(R_par *rp, float u, float v)
{
  return rp->doproj ? fabs(u*rp->proj.sinphi + v*rp->proj.cosphi) :
                      sqrt(u*u+v*v);
}

/*.......................................................................
 * Return the appropriate amplitude of a visibility for plotting.
 *
 * Input:
 *  rp         R_par *   The plot descriptor.
 *  vis   Visibility *   The visibility being plotted.
 * Output:
 *  return     float     The requested amplitude.
 */
static float r_vis_amp(R_par *rp, Visibility *vis)
{
/*
 * Are we plotting difference between the model and the data?
 */
  if(rp->dodiff) {
    float re = vis->amp * cos(vis->phs) - vis->modamp * cos(vis->modphs);
    float im = vis->amp * sin(vis->phs) - vis->modamp * sin(vis->modphs);
    return sqrt(re * re + im * im);
/*
 * Or are we plotting the amplitude directly?
 */
  } else {
    return vis->amp;
  };
}

/*.......................................................................
 * Return the appropriate phase of a visibility for plotting.
 *
 * Input:
 *  rp         R_par *   The plot descriptor.
 *  vis   Visibility *   The visibility being plotted.
 *  u,v        float     The U and V coordinates of the visibility.
 * Output:
 *  return     float     The requested phase.
 */
static float r_vis_phs(R_par *rp, Visibility *vis, float u, float v)
{
  float phs;  /* The phase being computed. */
/*
 * Are we plotting difference between the model and the data?
 */
  if(rp->dodiff) {
    float re = vis->amp * cos(vis->phs) - vis->modamp * cos(vis->modphs);
    float im = vis->amp * sin(vis->phs) - vis->modamp * sin(vis->modphs);
    phs = re==0.0 && im==0.0 ? 0.0: atan2(im, re);
  } else {
    phs = vis->phs;
  };
/*
 * Wrap the phase into the range -pi to pi.
 */
  phs -= twopi * floor(phs/twopi+0.5);
/*
 * Correct for phase conjugation. For projected distances flip the sign
 * if the distance is negative. For unprojected radii, flip the sign
 * of negative U data, as done in invert.
 */
  if(rp->doproj) {
    if(u*rp->proj.sinphi + v*rp->proj.cosphi < 0.0f)
      phs = -phs;
  } else {
    if(u < 0.0f)
      phs = -phs;
  };
  return phs;
}

/*.......................................................................
 * Return the model phase of a visibility for plotting.
 *
 * Input:
 *  rp         R_par *   The plot descriptor.
 *  vis   Visibility *   The visibility being plotted.
 *  u,v        float     The U and V coordinates of the visibility.
 * Output:
 *  return     float     The requested phase.
 */
static float r_mod_phs(R_par *rp, Visibility *vis, float u, float v)
{
  float phs = vis->modphs;  /* The phase being computed. */
/*
 * Wrap the phase into the range -pi to pi.
 */
  phs -= twopi * floor(phs/twopi+0.5);
/*
 * Correct for phase conjugation. For projected distances flip the sign
 * if the distance is negative. For unprojected radii, flip the sign
 * of negative U data, as done in invert.
 */
  if(rp->doproj) {
    if(u*rp->proj.sinphi + v*rp->proj.cosphi < 0.0f)
      phs = -phs;
  } else {
    if(u < 0.0f)
      phs = -phs;
  };
  return phs;
}

/*.......................................................................
 * Return the uncertainty of a visibility.
 *
 * Input:
 *  rp         R_par *   The plot descriptor.
 *  vis   Visibility *   The visibility being plotted.
 * Output:
 *  return     float     The requested phase.
 */
static float r_vis_err(R_par *rp, Visibility *vis)
{
  return vis->wt!=0.0f ? 1.0f/sqrt(fabs(vis->wt)) : 0.0f;
}

/*.......................................................................
 * Toggle plotting flags given a command key.
 *
 * Input:
 *  rp     R_par *  The plot descriptor.
 *  key      int    The command key (upper case).
 *  waslow   int    If true the key was lower case.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Key unrecognized.
 */
static int r_flags(R_par *rp, int key, int waslow)
{
  switch (key) {
  case KEY_MOD:      /* Toggle the display-model flag */
    rp->domod = !rp->domod;
    break;
  case KEY_AMP:
    rp->doamp = 1;
    rp->dophs = 0;
    break;
  case KEY_PHS:
    rp->doamp = 0;
    rp->dophs = 1;
    break;
  case KEY_BOTH:
    rp->doamp = rp->dophs = 1;
    break;
  case KEY_DIFF:
    rp->dodiff = !rp->dodiff;
    break;
  case KEY_ERR:
    rp->doerr = !rp->doerr;
    break;
  default:
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Allow the user to enter a new projection angle.
 *
 * Input:
 *  rp    R_par *   The plot descriptor.
 * Output:
 *  return  int     0 - OK.
 */
static int r_getphi(R_par *rp)
{
  enum {BLEN=40};  /* The length of the input buffer */
  char buf[BLEN];  /* The buffer to read the user's reponse into */
  double newphi;   /* The new projection angle (degrees) */
  char *endp;      /* Pointer to first unconverted character in input string */
/*
 * Prompt the user to enter an angle.
 */
  printf("Enter a new projection angle (degrees): ");
/*
 * Read user specification.
 */
  if(fgets(buf, BLEN, stdin)==NULL) {
    lprintf(stderr, "Error reading projection angle.\n");
    return 0;
  };
/*
 * Input too long?
 */
  if(!strchr(buf, '\n')) {
    lprintf(stderr, "Projection angle input too long.\n");
/*
 * Discard input up to the next newline.
 */
    {
      int c;
      do c=getc(stdin); while(c != EOF && c!='\n');
    };
    return 0;
  };
/*
 * Skip leading white-space.
 */
  endp = &buf[0];
  while(isspace((int)*endp))
    endp++;
/*
 * Nothing in string?
 */
  if(*endp == '\0') {
    lprintf(stdout, "Projection angle unchanged.\n");
    return 0;
  };
/*
 * Attempt to read the angle from the buffer.
 */
  newphi = strtod(endp, &endp);
/*
 * Skip trailing white-space.
 */
  while(isspace((int)*endp))
    endp++;
/*
 * Anything left in the string counts as an error.
 */
  if(*endp != '\0') {
    lprintf(stderr, "Bad projection angle input: %s\n", buf);
    return 0;
  };
/*
 * Assign the new projection angle.
 */
  return r_newphi(rp, newphi * dtor, 1);
}


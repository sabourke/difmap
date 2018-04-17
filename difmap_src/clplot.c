#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "obs.h"
#include "vlbconst.h"
#include "vlbutil.h"
#include "vlbmath.h"
#include "telspec.h"
#include "visplot.h"
#include "clphs.h"
#include "scans.h"
#include "cpgplot.h"
#include "logio.h"

static const float xmarg=0.05; /* The fraction of the X range for margin */
static const int datcol=10;    /* The color of unflagged data points */
static const int badcol=2;     /* The color of flagged data points */
static const int badccol=11;   /* The color of correction flagged data points */
static const int modcol=5;     /* The color to plot the model with */
static const int datsym=1;     /* The PGPLOT marker for good points */
static const int badsym=2;     /* The PGPLOT marker for flagged points */
static const int badcsym=5;    /* The marker of correction flagged points */
static const int cutcol=2;     /* PGPLOT color index for cut edit window */
static const int rescol=10;    /* PGPLOT color index for restore edit window */
static const int zoomcol=5;    /* PGPLOT color index for zoom cursor window */

/*
 * Define the selection keys.
 */
enum {
  KEY_NONE='\0',  /* Null key press */
  KEY_MODE=' ',   /* Toggle between triangle and baseline editing modes */
  KEY_CUR ='A',   /* Key for cursor position input */
  KEY_BRK ='B',   /* Toggle breaking into scans */
  KEY_CUT ='C',   /* Key to introduce cutting */
  KEY_CAN ='D',   /* Key to cancel incomplete select range */
  KEY_HELP='H',   /* Key to list usage information */
  KEY_ERR ='E',   /* Toggle display of error bars */
  KEY_FLG ='F',   /* Toggle display of flagged data */
  KEY_IF  ='I',   /* Toggle IF editing mode */
  KEY_DIS ='L',   /* Key to redisplay the current plot */
  KEY_MOD ='M',   /* Toggle model plotting */
  KEY_NXT ='N',   /* Key to display the next set of closure triangles */
  KEY_ORDER='O',  /* Key to toggle baseline order */
  KEY_PRV ='P',   /* Key to display the previous set of closure triangles */
  KEY_REST='R',   /* Key to introduce boxed-restore */
  KEY_NUMB='S',   /* Split into keyboard specified number of sub-plots */
  KEY_TEL ='T',   /* Key to request keyboard input of new triangle set */
  KEY_UT  ='U',   /* Key to select UT input mode */
  KEY_CH  ='W',   /* Toggle channel editing mode */
  KEY_QUIT='X',   /* Key to quit from this function */
  KEY_ZOOM='Z',   /* Zoom in on a user selected phase range or unzoom */
  KEY_PRVIF='[',  /* Show the previous IF */
  KEY_NXTIF=']',  /* Show the next IF */
  KEY_CROSS='+'   /* Toggle cross-hair cursor mode */
};

/* Scan dependent information container */

typedef struct {
  float vxa,vxb;  /* Min/max NDC X-coords of scan sub-plot */
  float sutmin;   /* UT range in scan is utmin -> utmax */
  float sutmax;
  float utmin;    /* The UT range to be displayed from this scan */
  float utmax;
  int view;       /* Flags whether any of scan is visible */
} Scans;

/* Sub-plot dependent information container */

typedef struct {
  float vya,vyb;  /* Min,max NDC Y-coords of viewport */
  Trispec ts;      /* The descriptor of the referenced closure triangle */
} Clssub;

/*
 * Enumerate the possible display selection modes in terms of the number
 * of reference indexes that they imply, ordered [isub,ta,tb,tc].
 */

typedef enum {REF_TEL=2, REF_BAS=3, REF_TRI=4} Nref;

/* Cursor selection container */

typedef struct {
  int key;       /* The upper-case value of the key used to return the cursor */
  int waslow;    /* True if 'key' was lower case */
  Clssub *cs;    /* Descriptor of the plot in which the cursor was pressed */
  int iplot;     /* The index of 'cs' in cp->cplots[] */
  Scans *sc;     /* The descriptor of the scan under the cursor */
  float utval;   /* The UT value under the cursor */
  float clphs;   /* The closure-phase value under the cursor */
} Clscurs;

/* The plot descriptor */

typedef struct {
  double utref;    /* Reference UT of UT axis */
  Observation *ob; /* The descriptor of the observation being plotted */
  Subarray *sub;   /* The descriptor of sub-array refsub */
  float utmin;     /* World min X coordinate (seconds wrt utref) */
  float utmax;     /* World max X coordinate (seconds wrt utref) */
  float utsum;     /* Sum of scan UT ranges currently visible */
  float phsmin;    /* World min phase coordinate (radians) */
  float phsmax;    /* World max phase coordinate (radians) */
  float vxa,vxb;   /* Viewport surrounding grid of sub-plots */
  float vya,vyb;
  Nref nref;       /* The number of indexes in the reference aggregate */
  int modified;    /* This remains 0 unless the data are edited */
  int if_ed;       /* If true, edits are restricted to the current IF */
  int ch_ed;       /* If true, edits are restricted to current freq channels */
  int tri_ed;      /* If true restrict edits to the 3 baselines of a triangle */
  int uta,utb;     /* The indexes of the first and last plotted integrations */
  int docurs;      /* True when cursor control is in effect */
  int doflag;      /* True if flagged points are or will be plotted */
  int domod;       /* True if model lines are or will be plotted */
  int dobars;      /* True if error bars are or will be plotted */
  int docross;     /* True to enable cross-hair mode */
  int doall;       /* Display all unique triangles to the fixed stations */
  int nreq;        /* The last requested number of plots (>= nrow) */
  int nrow;        /* The number of sub-plot slots on the display */
  int nplot;       /* The actual number of sub-plots plotted */
  Clssub *cplots;  /* Dynamic array of nrow subplot descriptors */
  int doscan;      /* Break into scans if true */
  int nscan;       /* Number of scans in 'scans' array */
  Scans *scans;    /* Array of scan descriptors */
  Clscurs cursor;  /* The descriptor of the last cursor selection */
  int npage;       /* The sequential number of the page being plotted */
} Clspar;

static int c_newut(Clspar *cp);
static int c_newphs(Clspar *cp);
static int c_newnum(Clspar *cp);
static int c_flags(Clspar *cp, char key, int waslow);

static Clspar *new_Clspar(Observation *ob, int docurs, int doscan, int doflag,
		   int domod, int dobars, int nrow);
static Clspar *del_Clspar(Clspar *cp);

static int get_Scans(Clspar *cp);

static int c_pldata(Clspar *cp, Clssub *cs, int uta, int utb, int erase);
static int c_plmodel(Clspar *cp, Clssub *cs, int erase);

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

static int c_cursor(Clspar *cp, int noout, Bandmode mode,
		    Clssub *csref, float xref, float yref, int ci);

static int c_label(Clspar *cp);
static int c_redisp(Clspar *cp);
static int c_mlab(Clspar *cp, int erase);
static int c_newmode(Clspar *cp, int tri_ed, int if_ed, int ch_ed);
static int c_edit(Clspar *cp, Clssub *cs, int flag, int ut);
static int c_toggle(Clspar *cp, Clscurs *cc);
static int c_box(Clspar *cp, int doflag);

/* Define an enum used to specify what c_plot should plot next */

typedef enum {
  C_ALLNEW,     /* Start new plot */
  C_REPLOT,     /* Replot the current page */
  C_RESET,      /* Reset to the first/end triangle of the current selection */
  C_NXT_SUB,    /* Skip to the start of the next sub-array */
  C_NXT_TC,     /* Plot next page of triangles */
  C_NXT_TRI     /* As T_NXT_TC, but stop when all selected triangles */
                /* have been seen. */
} Clsop;

static int c_plot(Clspar *cp, Clsop oper, int forward, Trispec *refts);
static int c_setnrow(Clspar *cp, int nreq);
static int c_utrange(Clspar *cp);
static int c_cpwin(Clspar *cp, int nrow, int nplot);
static int c_plaxes(Clspar *cp, Clssub *cs, int dobot, int dotop, int erase);
static int c_scale(Clspar *cp, Clssub *cs, float *xtomm, float *ytomm);
static int c_find(Clspar *cp, Clscurs *cc);

/*.......................................................................
 * Receive input of a new UT range via the cursor and redisplay the plot
 * within that range. If the user presses the KEY_UT key then display
 * the full UT range.
 *
 * Input:
 *  cp      Clspar *  The plot descriptor.
 * Output:
 *  return     int    0 - OK. Anything else means fatal error.
 */
static int c_newut(Clspar *cp)
{
  int dofull=0;   /* True if full UT range requested */
  int iter;       /* Iterate over getting two valid keypresses */
  float utval[2]={0.0f,0.0f}; /* The two UT end points wrt cp->utref. */
/*
 * Get the two UT cursor selections.
 */
  for(iter=0; iter<2 && !dofull; iter++) {
    int accepted = 0;
    while(!accepted) {
      if(c_cursor(cp, 1, iter==0?B_XVAL:B_XRNG, NULL, utval[0], 0.0f, zoomcol))
	return 1;
      switch(cp->cursor.key) {
      case KEY_UT:        /* Display full range and end UT selection */
	accepted = dofull = 1;
	break;
      case KEY_QUIT: case KEY_CAN:     /* Abort UT selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected UT */
	utval[iter] = cp->cursor.utval;
	accepted=1;
	break;
      default:                  /* Unexpected cursor input key - show usage */
	printf("%c - Select the position of the %s UT.\n", KEY_CUR,
	       iter==0 ? "start":"end");
	printf("%c - Abort selection.\n", KEY_CAN);
	printf("%c - Revert to the full UT range.\n", KEY_UT);
	break;
      };
    };
  };
/*
 * Get the UT indexes corresponding to the selected UTs.
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
 * Receive input of a new closure-phase range via the cursor and
 * redisplay the plot within that range. If the user presses the KEY_ZOOM
 * key then display the full -180 to 180 degree range.
 *
 * Input:
 *  cp      Clspar *  The plot descriptor.
 * Output:
 *  return     int    0 - OK. Anything else means fatal error.
 */
static int c_newphs(Clspar *cp)
{
  int dofull=0;   /* True if full phase range requested */
  int iter;       /* Iterate over getting two valid keypresses */
  float phase[2]={0.0f,0.0f}; /* The two phases selected. */
  int iplot[2];   /* The indexes of the sub-plots of the cursor selection */
  Clssub *cs[2]={NULL,NULL}; /* Subarray descriptors of the cursor selections */
/*
 * Get the two phase cursor selections.
 */
  for(iter=0; iter<2 && !dofull; iter++) {
    int accepted = 0;
    while(!accepted) {
      if(c_cursor(cp, 1, iter==0?B_YVAL:B_YRNG, cs[0], 0.0f, phase[0], zoomcol))
	return 1;
      switch(cp->cursor.key) {
      case KEY_ZOOM:    /* Display full range and end phase range selection */
	accepted = dofull = 1;
	break;
      case KEY_QUIT: case KEY_CAN:     /* Abort selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected phase */
	if(iter==1 && cp->cursor.clphs == phase[0]) {
	  printf("Second phase identical to first - please redo the second.\n");
	} else {
	  phase[iter] = cp->cursor.clphs;
	  iplot[iter] = cp->cursor.iplot;
	  cs[iter] = cp->cursor.cs;
	  accepted=1;
	};
	break;
      default:                  /* Unexpected cursor input key - show usage */
	printf("%c - Select the %s phase of the phase range.\n", KEY_CUR,
	       iter==0 ? "start":"end");
	printf("%c - Abort selection.\n", KEY_CAN);
	printf("%c - Revert to the full -180 to 180 degree range.\n", KEY_ZOOM);
	break;
      };
    };
  };
/*
 * Assign the selected range.
 */
  if(dofull) {
    cp->phsmin = -pi;
    cp->phsmax = pi;
  } else {
/*
 * If the second cursor selection was outside the sub-plot of the first
 * then substitute the phase of the straddled edge.
 */
    if(iplot[1] < iplot[0])
      phase[1] = cp->phsmax;
    else if(iplot[1] > iplot[0])
      phase[1] = cp->phsmin;
/*
 * Sort phase[] into ascending order.
 */
    if(phase[1] < phase[0]) {
      float phs = phase[0];
      phase[0] = phase[1];
      phase[1] = phs;
    };
/*
 * Set limits.
 */
    if(phase[0] < cp->phsmin)
      phase[0] = cp->phsmin;
    if(phase[1] > cp->phsmax)
      phase[1] = cp->phsmax;
/*
 * The above limits may have forced the two phases to be equal.
 */
    if(phase[0] == phase[1]) {
      printf("The two phases are indentical - selection aborted.\n");
      return 0;
    };
/*
 * Install the new limits.
 */
    cp->phsmin = phase[0];
    cp->phsmax = phase[1];
  };
/*
 * Display the result.
 */
  return c_redisp(cp);
}

/*.......................................................................
 * Request a new value for the number of sub-plots per page.
 *
 * Input:
 *  cp       Clspar *  The plot descriptor.
 * Output:
 *  return      int    0 - OK.
 */
static int c_newnum(Clspar *cp)
{
  char awrk[20]; /* Temporary work array */
  int nreq;      /* New number of plots. */
/*
 * Get the new number from the user.
 */
  printf("Enter the required number of plots per page: ");
  if(fgets(awrk, sizeof(awrk)-1, stdin) == NULL) {
    fprintf(stderr, "Error reading input.\n");
    return 0;
  };
/*
 * Read the number.
 */
  nreq = atoi(awrk);
/*
 * Register the requested change in numbers of plot slots.
 */
  c_setnrow(cp, nreq);
/*
 * Plot the new number of plots.
 */
  return c_plot(cp, C_REPLOT, 1, NULL) < 0;
}

/*.......................................................................
 * Cursor control interface to plotting and editing functions.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation to be plotted.
 *  ts      Trispec *  The descriptor of the initial triangle to plot, or
 *                     NULL to select the first triangle.
 *                     In non-interactive mode - (docurs=0) - all
 *                     triangles corresponding to unique combinations of the
 *                     variable telescope indexes will be plotted.
 *  cif         int    The index of the start IF.
 *  nrow        int    Number of sub-plots in Y.
 *  npage       int    The max number of pages to plot if non-interactive.
 *                     A value of 0 means no page limit.
 *  docurs      int    If 0 disallow cursor control.
 *  opts       char *  An optional string of flag toggling keys that
 *                     toggle after the values below have been applied.
 *                     Send NULL or empty string if not required.
 *  doscan      int    If true break up plot into scans if present.
 *  doflag      int    If true plot flagged data in addition to unflagged.
 *  domod       int    If true plot model data in addition to unflagged.
 *  dobars      int    If true plot error bars.
 * Input/Output:
 *  modified    int *  If modified!=NULL then *modified will be assigned
 *                     with 0 if no data were edited and 1 if data were
 *                     edited.
 * Output:
 *  return      int    0 - OK.
 */
int clsplot(Observation *ob, Trispec *ts, int cif, int nrow, int npage,
	    int docurs, char *opts, int doscan, int doflag, int domod,
	    int dobars, int *modified)
{
  Clspar *cp;     /* Plot descriptor */
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
 * Check inputs.
 */
  if(!ob_ready(ob, OB_SELECT, "cpplot"))
    return 1;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Substitute default starting triangle?
 */
  if((ts && next_tri(ob, FIND_FIRST, 1, ts->nfix<2?2:ts->nfix, 1, 0, 1, ts)) ||
     (!ts && !(ts=find_tri(ob, 0, 0, 0, 0, 0, 1, 2, 1, 0, 1))))
    return 1;
/*
 * An IF index of -1 (0 on the command line) requests the default IF,
 * substitute the first unsampled IF.
 */
  if(cif == -1) {
    if((cif = nextIF(ob, 0, 1, 1)) < 0) {
      lprintf(stderr, "cpplot: There are no selected IFs available.\n");
      return 1;
    };
  } else if(cif < 0 || cif >= ob->nif) {
    lprintf(stderr, "cpplot: IF %d does not exist.\n", cif+1);
    return 1;
  };
/*
 * Attempt to read the start IF.
 */
  if(getIF(ob, cif))
    return 1;
/*
 * Allocate and initialize the plot descriptor.
 */
  cp = new_Clspar(ob, docurs, doscan, doflag, domod, dobars, nrow);
  if(cp==NULL)
    return 1;
/*
 * If a string of flag options was given, interpret them here.
 */
  if(opts != NULL) {
    int slen = strlen(opts);
    for(i=0; i<slen; i++) {
      char key = opts[i];
      int waslow = islower((int)key);
      if(waslow)
	key = toupper((int)key);
      c_flags(cp, key, waslow);
/*
 * Handle flags that are normally treated as commands.
 */
      switch(key) {
      case KEY_ORDER:
	cp->doall = !cp->doall;
	break;
      case KEY_MODE:
	cp->tri_ed = !cp->tri_ed;     /* Toggle triangle flagging mode */
	break;
      case KEY_IF:
	cp->if_ed = !cp->if_ed;       /* Toggle IF flagging */
	break;
      case KEY_CH:
	cp->ch_ed = !cp->ch_ed;       /* Toggle channel flagging */
	break;
      case KEY_CROSS:
	cp->docross = !cp->docross;
	break;
      };
    };
  };
/*
 * Plot the first page.
 */
  ierr = c_plot(cp, C_ALLNEW, 1, ts) <= 0;
/*
 * Interactive plotting?
 */
  if(cp->docurs) {
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
      int nflag=0;    /* Number of option-flag toggling operations done */
/*
 * Read the cursor.
 */
      do {
	ierr = c_cursor(cp, 0, B_NORM, NULL, 0.0f, 0.0f, 1);
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
	ierr = c_plot(cp, C_REPLOT, 1, NULL) < 0;
      } else {
	switch (cp->cursor.key) {
	case KEY_CUR:      /* (un)flag baseline(s) of closest poin to cursor */
	  ierr = c_toggle(cp, &cp->cursor);
	  break;
	case KEY_CUT:     /* Flag points within a select-box */
	  ierr = c_box(cp, 1);
	  break;
	case KEY_REST:    /* Restore points inside a select-box */
	  ierr = c_box(cp, 0);
	  break;
	case KEY_MODE:     /* Toggle triangle editing mode */
	  ierr = c_newmode(cp, !cp->tri_ed, cp->if_ed, cp->ch_ed);
	  break;
	case KEY_IF:       /* Toggle IF editing mode */
	  ierr = c_newmode(cp, cp->tri_ed, !cp->if_ed, cp->ch_ed);
	  break;
	case KEY_CH:       /* Toggle channel editing mode */
	  ierr = c_newmode(cp, cp->tri_ed, cp->if_ed, !cp->ch_ed);
	  break;
	case KEY_DIS:      /* Re-display current plot */
	  ierr = c_plot(cp, C_REPLOT, 1, NULL) < 0;
	  break;
	case KEY_NXT:                     /* Plot next page */
	  ierr = c_plot(cp, cp->cursor.waslow ? C_NXT_TC:C_NXT_SUB, 1, NULL)<0;
	  break;
	case KEY_PRV:                     /* Plot previous page */
	  ierr = c_plot(cp, cp->cursor.waslow ? C_NXT_TC:C_NXT_SUB, 0, NULL)<0;
	  break;
	case KEY_TEL:                     /* New reference telescope */
/*
 * Get the user telescope/baseline request, and if valid, plot the result.
 */
	  {
	    Trispec *ts = read_Trispec(ob, NULL, NULL, cp->cplots[0].ts.isub);
	    ierr = ts && c_plot(cp, C_ALLNEW, 1, ts) < 0;
	  };
	  break;
	case KEY_UT:      /* Following cursor input selects UT ranges */
	  ierr = c_newut(cp);
	  break;
	case KEY_ZOOM:    /* Following cursor input selects phase range */
	  ierr = c_newphs(cp);
	  break;
	case KEY_NUMB:    /* Change number of plots per page */
	  ierr = c_newnum(cp);
	  break;
	case KEY_CROSS:   /* Toggle cross-hair cursor mode */
	  cp->docross = !cp->docross;
	  break;
	case KEY_ORDER:
	  cp->doall = !cp->doall;
	  if(c_plot(cp, C_RESET, 1, NULL) <= 0)
	    cp->doall = !cp->doall;
	  break;
	case KEY_PRVIF:
	case KEY_NXTIF:
	  {
	    int step = cp->cursor.key==KEY_NXTIF ? 1 : -1;
	    int cif = nextIF(ob, ob->stream.cif + step, 1, step);
	    ierr = cif >= 0 && (getIF(ob, cif) || c_redisp(cp));
	  };
	  break;
	case KEY_HELP:    /* Print usage info */
	  printf("Clsplot key bindings:\n");
	  printf(" %c - (right-mouse-button) exit clsplot.\n", KEY_QUIT);
	  printf(" %c - List key bindings.\n", KEY_HELP);
	  printf(" %c - (left-mouse-button) (un)flag baselines of closest point to cursor.\n", KEY_CUR);
	  printf(" %c - Flag all data inside a specified rectangular box.\n",
		 KEY_CUT);
	  printf(" %c - Restore data inside a specified rectangular box.\n",
		 KEY_REST);
	  printf(" %c - Toggle between baseline and triangle based editing.\n",
		 KEY_MODE);
	  printf(" %c - Toggle IF based editing.\n", KEY_IF);
	  printf(" %c - Toggle spectral-line channel based editing.\n",KEY_CH);
	  printf(" %c - Select UT range to be displayed.\n", KEY_UT);
	  printf(" %c - Zoom in or out on a selected phase range.\n",KEY_ZOOM);
	  printf(" %c - Redisplay current plot.\n", KEY_DIS);
	  printf(" %c - Plot the next set of triangles.\n", tolower(KEY_NXT));
	  printf(" %c - Plot the previous set of triangles.\n",
		 tolower(KEY_PRV));
	  printf(" %c - Skip to the next sub-array.\n", KEY_NXT);
	  printf(" %c - Skip to the preceding sub-array.\n", KEY_PRV);
	  printf(" %c - Plot from the next IF.\n", KEY_NXTIF);
	  printf(" %c - Plot from the preceding IF.\n", KEY_PRVIF);
	  printf(" %c - Toggle display of model data.\n", KEY_MOD);
	  printf(" %c - Toggle display of flagged data.\n", KEY_FLG);
	  printf(" %c - Toggle display of error bars.\n", KEY_ERR);
	  printf(" %c - Toggle whether to use a cross-hair cursor if available.\n", KEY_CROSS);
	  printf(" %c - Select number of sub-plots per page.\n", KEY_NUMB);
	  printf(" %c - Toggle between seeing all or just upper triangles.\n",
		 KEY_ORDER);
	  printf(" %c - Toggle breaking up plot into scans.\n", KEY_BRK);
	  printf(" %c - Prompt for a new set of closure triangles.\n", KEY_TEL);
	  printf("\n");
	  break;
	case KEY_QUIT:  /* Quit plotting session */
	  finished = 1;
	  break;
	default:
	  break;
	};
      };
    };
/*
 * Non-interactive plotting?
 */
  } else if(!ierr) {
/*
 * Plot as many pages as required, and keep the user informed.
 */
    while(!ierr && (npage <= 0 || cp->npage < npage)) {
      int nplotted = c_plot(cp, C_NXT_TRI, 1, NULL);
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
    *modified = cp->modified;
/*
 * Clean up.
 */
  cp = del_Clspar(cp);
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
 *  cp    Clspar *  The plot descriptor.
 *  key     char    The command key (upper case).
 *  waslow   int    If true the key was lower case.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Key unrecognized.
 */
static int c_flags(Clspar *cp, char key, int waslow)
{
  switch (key) {
  case KEY_MOD:      /* Toggle the display-model flag */
    cp->domod = !cp->domod;
    break;
  case KEY_FLG:      /* Toggle the display-flagged-data flag */
    cp->doflag = !cp->doflag;
    break;
  case KEY_ERR:
    cp->dobars = !cp->dobars;
    break;
  case KEY_BRK:
    cp->doscan = !cp->doscan;
    break;
  default:
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Create a new Clspar descriptor and contained Clssub array, and
 * assign given defaults.
 *
 * Input:
 *  ob Observation *  The descriptor of the observation.
 *  docurs     int    True if cursor control is required. If the current
 *                    device has no cursor, this will be ignored.
 *  doscan     int    Default for request to split data into scans.
 *  doflag     int    Default for request for flagged data to be plotted.
 *  domod      int    Default for request for model to be plotted.
 *  dobars     int    Default for request for error bars to be plotted.
 *  nrow       int    The initial number of sub-plots. If nrow<=0
 *                    the default sub->nstat-1 will be used.
 * Output:
 *  return  Clspar *  The new Clspar descriptor, or NULL on error.
 */
static Clspar *new_Clspar(Observation *ob, int docurs, int doscan, int doflag,
		    int domod, int dobars, int nrow)
{
  Clspar *cp;     /* Pointer to the new Clspar descriptor */
  int slen;       /* Length of string */
/*
 * Allocate the Clspar descriptor.
 */
  cp = (Clspar *) malloc(sizeof(Clspar));
  if(cp == NULL) {
    lprintf(stderr, "new_Clspar: Insufficient memory for plot descriptor.\n");
    return cp;
  };
/*
 * NULLify pointer members so that del_Clspar can be called before the
 * struct has been fully initialized.
 */
  cp->cplots = 0;
  cp->scans = 0;
/*
 * Record the descriptor of the observation.
 */
  cp->ob = ob;
  cp->sub = NULL;
/*
 * Assign passed defaults.
 */
  cp->doflag = doflag;
  cp->domod = domod;
  cp->dobars = dobars;
  cp->docross = 0;
  cp->doall = 1;
/*
 * Assign defaults to the rest of the members.
 */
  cp->nreq = nrow;
  cp->nrow = 0;
  cp->nplot = 0;
  cp->utmin = cp->utmax = cp->utsum = 0.0f;
  cp->phsmin = -pi;
  cp->phsmax = pi;
  cp->vxa = cp->vxb = cp->vya = cp->vyb = 0.0f;
  cp->uta = 0;
  cp->utb = 0;
  cp->utref = ob->date.ut;
  cp->nscan = 0;
  cp->doscan = doscan;
  cp->npage = 0;
/*
 * Mark the data as unmodified.
 */
  cp->modified = 0;
/*
 * Set defaults for IF and spectral-line channel editing.
 */
  cp->if_ed = 0;  /* Global IF editing */
  cp->ch_ed = 0;  /* Global channel editing */
  cp->tri_ed = 0;
/*
 * If cursor interaction is required, check if the device has a cursor.
 */
  if(docurs) {
    char answer[5];
    slen = sizeof(answer)-1;
    cpgqinf("CURSOR", answer, &slen);
    answer[3]='\0';
    docurs = (strncmp(answer, "YES" ,3) == 0);
  };
  cp->docurs = docurs;
  cp->cursor.key = KEY_NONE;
/*
 * Return the new descriptor.
 */
  return cp;
}

/*.......................................................................
 * Clspar (visibility plot descriptor) destructor function.
 *
 * Input:
 *  cp     Clspar *  Clspar pointer returned by new_Clspar().
 * Output:
 *  return Clspar *  Always NULL, so that you can write cp=del_Clspar(cp);
 */
static Clspar *del_Clspar(Clspar *cp)
{
  if(cp) {
    if(cp->cplots)
      free(cp->cplots);
    if(cp->scans)
      free(cp->scans);
    free(cp);
  };
  return NULL;
}

/*.......................................................................
 * Determine a new set of scans from a new time separator - and/or new
 * sub-array.
 *
 * Input:
 *  cp         Clspar *  The plot parameter container.
 *                        cp->doscan: If 0 the whole sub-array is
 *                                    treated as one scan.
 *                        cp->scans:  Must be NULL on first call.
 *                        cp->nscan:  Must be 0 on first call.
 * Input/Output:
 *  cp         Clspar *  The plot descriptor.
 *                       cp->scans will be allocated via malloc if
 *                       (cp->scans==NULL or cp->nscan==0).
 *                       Otherwise cp->scans will be realloc'd to
 *                       the new array size. cp->nscans will be set with
 *                       the number of scans initialized.
 *                       On error cp will be left unchanged.
 *  return        int    The number of scans initialized - this may be
 *                       0 if a memory-allocation failure occurs.
 */
static int get_Scans(Clspar *cp)
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
    lprintf(stderr, "get_Scans: NULL plot descriptor intercepted\n");
    return 0;
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
    lprintf(stderr, "get_Scans: Insufficient memory for new Scans\n");
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
 * Return the UT plot range for the ut range and plot options
 * in a passed Clspar descriptor.
 *
 * Input/Output:
 *  cp        Clspar * On entry this contains existing plotting attributes.
 *                     Currently only uta and utb need be initialized.
 *                     On output cp->utmin and cp->utmax will contain
 *                     the min and max UTs of the plot in seconds
 *                     wrt cp->utref, including margins.
 * Output:
 *  return      int    0 - OK.
 *                     On error -1 is returned and no changes are made
 *                     to *utmin or *utmax.
 */
static int c_utrange(Clspar *cp)
{
  float xa;   /* Start UT of range */
  float xb;   /* End UT of range */
  int scan;   /* The index of the scan being processed */
  Scans *sc;  /* The scan descriptor being processed */
/*
 * Check inputs.
 */
  if(cp==0) {
    lprintf(stderr, "c_utrange: NULL Clspar descriptor intercepted\n");
    return -1;
  };
/*
 * Valid uta and utb?
 */
  if(cp->uta < 0 || cp->uta>cp->utb || cp->utb >= cp->sub->ntime) {
    lprintf(stderr, "c_utrange: uta and utb are invalid.\n");
    return -1;
  };
/*
 * Determine the times corresponding to integrations uta and utb
 * with respect to the reference time vlb->ut (seconds).
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
 * Set up the viewport limits for the stack of plots leaving 4 char
 * heights on each side of plot for labelling.
 *
 * Input/Output:
 *  cp    Clspar *  The Cls-edit parameter struct. On input cp->nplot
 *                  must be set with the number of closure-phase plots
 *                  required. cp->cplots[] must have been pre-allocated
 *                  and on output the vxa,vxb,vya,vyb fields will be
 *                  initialized. All other fields are ignored.
 *  nrow     int    The number of sub-plot slots on the display
 *  nplot    int    The actual number of sub-plots to be plotted.
 *                  nplot <= nrow.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
static int c_cpwin(Clspar *cp, int nrow, int nplot)
{
  Clssub *cs;     /* The sub-plot being assigned a viewport. */
  float vxa,vxb;  /* X viewport limits enclosing whole stack */
  float vya,vyb;  /* Y viewport limits enclosing whole stack */
  float utsum;    /* Sum of scan UTs within current UT range */
  Scans *sc;      /* A scan descriptor from cp->scans[] */
  int scan;       /* Number of scan */
  int i;
/*
 * Check arguments.
 */
  if(nplot > nrow) {
    lprintf(stderr, "c_cpwin: Too many plots requested\n");
    return 1;
  };
  if(nplot <= 0) {
    lprintf(stderr, "c_cpwin: %d plots requested\?", nplot);
    return 1;
  };
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
  cp->nplot = nplot;
/*
 * Divide it into cp->nplot vertically adjacent viewports.
 */
  for(i=0; i<nplot; i++) {
    cs = &cp->cplots[i];
    cp->vxa = vxa;
    cp->vxb = vxb;
    cs->vyb = vyb - i*(vyb-vya)/nrow;
    cs->vya = cs->vyb - (vyb-vya)/nrow;
  };
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
/*
 * Scale the character height with the number of plots.
 */
  cpgsch(3.0f/cp->nplot);
  return 0;
}

/*.......................................................................
 * Draw axes for a given sub-plot.
 *
 * Input:
 *  cp      Clspar *  The plot descriptor.
 *  cs      Clssub *  The sub-plot descriptor.
 *  dotop      int    If true draw ticked axis along top axis of viewport.
 *  dobot      int    If true draw ticked and labelled axis along bottom
 *                    of the sub-plot viewport.
 *  erase      int    If true erase current axes instead of plotting.
 * Output:
 *  return     int    0 - OK.
 *                    Anything else if an error occured.
 */
static int c_plaxes(Clspar *cp, Clssub *cs, int dotop, int dobot, int erase)
{
  Subarray *sub; /* Local pointer to sub-array descriptor */
  Scans *sc;     /* The scan being labelled */
  float utmin;   /* Start UT + 1 day, in seconds since start of year */
  float utmax;   /* End UT + 1 day, in seconds since start of year */
  float ch;      /* Character height to use */
  int oldcol;    /* Color index on entry to function */
  int scan;      /* The scan index being processed */
  int first;     /* Index of first visible scan */
  int last;      /* Index of final visible scan */
  char label[25];/* Temporary string to compose baseline label in */
/*
 * Check arguments.
 */
  if(cp==NULL || cs==NULL) {
    lprintf(stderr, "c_plaxes: NULL %s descriptor intercepted\n",
	    (cp==NULL)?"plot":"sub-plot");
    return -1;
  };
/*
 * Get a local pointer to the current sub-array descriptor.
 */
  sub = cp->sub;
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
  ch = 1.0f/sqrt((double) cp->nplot);
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
  cpgsvp(cp->vxa, cp->vxb, cs->vya, cs->vyb);
  cpgswin(0.0f, 1.0f, cp->phsmin * rtod, cp->phsmax * rtod);
  cpgbox(" ", 0.0f, 0, "BCVNST", 0.0f, 0);
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
      cpgmove(sc->vxa, cs->vya);
      cpgdraw(sc->vxa, cs->vyb);
    };
    if(scan != last) {
      cpgmove(sc->vxb, cs->vya);
      cpgdraw(sc->vxb, cs->vyb);
    };
/*
 * Write numeric labels under the last plot.
 */
    cpgsvp(sc->vxa, sc->vxb, cs->vya, cs->vyb);
    cpgswin(utmin, utmax, 0.0f, 1.0f);
    cpgsch(dotop ? 0.7f : ch);
    cpgtbox("ZHCST", 0.0f, 0, " ", 0.0f, 0);
    cpgsch(dobot ? 0.7f : ch);
    cpgtbox(dobot?"ZHBNST":"ZHBST", 0.0f, 0, " ", 0.0f, 0);
  };
/*
 * Set viewport around whole sub-plot and write a closure-triangle label
 * inside the top right hand corner.
 */
  cpgsvp(cp->vxa, cp->vxb, cs->vya, cs->vyb);
  sprintf(label, "%.10s-%.10s-%.10s", sub->tel[cs->ts.ta].name,
	  sub->tel[cs->ts.tb].name, sub->tel[cs->ts.tc].name);
  cpgsch(0.5f);
  cpgmtxt("T", -1.5f, 0.99f, 1.0f, label);
/*
 * Restore entry color.
 */
  cpgsci(oldcol);
  return 0;
}

/*.......................................................................
 * Plot or erase closure-phase points.
 *
 * Input:
 *  cp       Clspar * The plot descriptor.
 *  cs       Clssub * The subplot descriptor.
 *  uta         int   Index of first integration to be plotted.
 *  utb         int   Index of second integration to be plotted.
 *  erase       int   If true erase points instead of plotting them.
 * Output:
 *  return      int   0 - OK.
 */
static int c_pldata(Clspar *cp, Clssub *cs, int uta, int utb, int erase)
{
  Subarray *sub;       /* Local pointer to the sub-array being displayed */
  Integration *integ;  /* The descriptor of the integration being displayed */
  Scans *sc;           /* The scan being plotted */
  static int oldcol;   /* Color index on entry to function */
  static int ut;       /* Index of integration being plotted */
  static int first;    /* True until the first point has been plotted */
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
 * Get the appropriate sub-array descriptor.
 */
  sub = cp->sub;
/*
 * Draw each point with the appropriate symbol and color for its
 * flag status.
 */
  sc = &cp->scans[0];
  first=1;
  for(ut=uta,integ = &sub->integ[ut]; ut<= utb; ut++,integ++) {
    Clphs *cphs = get_clphs(&cs->ts, integ->vis);
    float utval = integ->ut - cp->utref;
/*
 * Skip to the right scan for this point.
 */
    if(first || utval > sc->sutmax) {
      first = 0;
      while(utval > sc->sutmax)
	sc++;
      cpgsvp(sc->vxa, sc->vxb, cs->vya, cs->vyb);
      cpgswin(sc->utmin, sc->utmax, cp->phsmin, cp->phsmax);
    };
/*
 * Ignore deleted data.
 */
    if(!(cphs->bad & FLAG_CDEL)) {
/*
 * Should we plot this point and if so, which color and symbol should
 * be used?
 */
      if(!cphs->bad || cp->doflag) {
/*
 * Get the world coordinates to be plotted.
 */
	float phs = cphs->ophs;
	float phserr = 1.0f / sqrt(fabs(cphs->wt));
/*
 * Get the marker type and color to be used to plot the new point.
 */
	int isym, icol;  /* Plot marker and color */
	if(cphs->bad) {               /* Flagged correction */
	  if(cphs->bad & FLAG_CBAD) { /* Visibility flag */
	    isym = badsym;
	    icol = badcol;
	  } else {                    /* Selfcal correction flag */
	    isym = badcsym;
	    icol = badccol;
	  };
	} else {                      /* Not flagged */
	  isym = datsym;
	  icol = datcol;
	};
/*
 * Install the new color.
 */
	cpgsci(erase ? 0 : icol);
/*
 * Plot the point.
 */
	cpgpt(1, &utval, &phs, isym);
/*
 * Plot error bars if requested.
 */
	if(cp->dobars) {
	  cpgmove(utval, phs - phserr);
	  cpgdraw(utval, phs + phserr);
	};
      };/* End of  if not deleted */
    };  /* End of  if not flagged or allowed to plot flagged data */
  };    /* End of  loop over integrations */
/*
 * Restore entry color and terminate pgplot buffering.
 */
  cpgsci(oldcol);
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Plot or erase closure-phase model lines.
 *
 * Input:
 *  cp       Clspar * The plot descriptor.
 *  cs       Clssub * The sub-plot descriptor.
 *  erase       int   If true erase points instead of plotting them.
 * Output:
 *  return      int   0 - OK.
 */
static int c_plmodel(Clspar *cp, Clssub *cs, int erase)
{
  Subarray *sub;  /* The descriptor of the current sub-array */
  Integration *integ; /* The descriptor of the integration being displayed */
  Scans *sc;      /* The scan being plotted */
  float prevphs=0.0f; /* The previous plotted phase */
  float prevut =0.0f; /* The previous UT value plotted */
  int oldcol;     /* Color index on entry to function */
  int first;      /* Is 0 after first point has been plotted */
  int ut;         /* The index of the integration being plotted */
/*
 * Do nothing if no model exists or cp->domod==0.
 */
  if(!cp->ob->hasmod || !cp->domod)
    return 0;
/*
 * Get the descriptor of the sub-array being plotted.
 */
  sub = cp->sub;
/*
 * Store current color index, to be restored on return.
 */
  cpgqci(&oldcol);
/*
 * Set color for drawing or erasing.
 */
  cpgsci(erase ? 0 : modcol);
/*
 * Turn on pgplot buffering.
 */
  cpgbbuf();
/*
 * Plot the closure-phase model as a line.
 */
  sc = &cp->scans[0];
  first = 1;
  for(ut=cp->uta, integ = &sub->integ[ut]; ut<=cp->utb; ut++,integ++) {
    Clphs *cphs = get_clphs(&cs->ts, integ->vis);
    float utval = integ->ut - cp->utref;
/*
 * Ignore deleted data.
 */
    if(!(cphs->bad & FLAG_CDEL)) {
      float phs = cphs->mphs;
/*
 * If not in the current scan - skip to the right scan and position
 * for the start of the new model line.
 */
      if(first || utval > sc->sutmax || utval-prevut > sub->scangap) {
	while(utval > sc->sutmax)
	  sc++;
	cpgsvp(sc->vxa, sc->vxb, cs->vya, cs->vyb);
	cpgswin(sc->utmin, sc->utmax, cp->phsmin, cp->phsmax);
	first = 0;
	cpgmove(utval, phs);
      } else {
	float phsdif = phs - prevphs;  /* Phase excursion from previous phase */
/*
 * Because the closure phase is only known modulo 360 degrees, there are three
 * possible paths between the previous and current phase point within the
 * chosen -180 to 180 degree plot range. Take the shortest of these.
 */
	if(phsdif > pi) {
	  cpgdraw(utval, phs-twopi);
	  cpgmove(prevut, prevphs+twopi);
	  cpgdraw(utval, phs);
	} else if(phsdif < -pi) {
	  cpgdraw(utval, phs+twopi);
	  cpgmove(prevut, prevphs-twopi);
	  cpgdraw(utval, phs);
	} else {
	  cpgdraw(utval, phs); /* Continue the model line */
	};
      };
      prevut = utval;
      prevphs = phs;
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
 * cp->cursor.
 *
 * Input:
 *  cp     Clspar *  The plot descriptor.
 *  noout     int    If true then don't return until the cursor is
 *                   pressed inside a sub-plot.
 *  mode Bandmode    The desired type of cursor, from:
 *                     B_NORM  -  A single point is required - no banding.
 *                     B_LINE  -  Line band between vcref and the cursor.
 *                     B_RECT  -  Rectangular band between vcref and the
 *                               cursor.
 *                     B_YRNG - Two horizontal lines bracketing a Y-axis range.
 *                     B_XRNG - Two vertical lines bracketing an X-axis range.
 *                     B_YVAL - Vertical line through the cursor.
 *                     B_XVAL - Horizontal line through the cursor.
 *                     B_CROSS- Cross hair centered on cursor.
 *  csref  Clssub *  The descriptor of the sub-plot to which xref,yref refer.
 *                   This can be NULL if mode==B_NORM.
 *  xref    float    xref and yref denote the reference position for
 *  yref    float    band-type cursors. They are ignored if mode==B_NORM.
 *  ci        int    The color index to draw band cursors with.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int c_cursor(Clspar *cp, int noout, Bandmode mode,
		    Clssub *csref, float xref, float yref, int ci)
{
  static float xpos=0.5f; /* The X NDC position of the cursor */
  static float ypos=0.5f; /* The Y NDC position of the cursor */
  Clscurs *cc;       /* Pointer to cp->cursor */
  Clssub *cs=NULL;   /* The subplot being checked */
  int iplot=0;       /* The index of 'cs' in cp->cplots[] */
  Scans *sc;         /* The scan descriptor being checked */
  int scan;          /* The scan number being processed */
  char key;          /* The cursor selection key */
  int found=0;       /* True if the cursor selection was found in a sub-plot */
/*
 * Get the cursor descriptor.
 */
  cc = &cp->cursor;  
/*
 * Set the viewport around the whole viewsurface and make the world
 * coords the same as NDC so that the returned cursor position
 * is measured in NDC.
 */
  cpgsvp(0.0f, 1.0f, 0.0f, 1.0f);
  cpgswin(0.0f, 1.0f, 0.0f, 1.0f);
/*
 * If this is the first call of the plot session initialize the position
 * at which to bring up the cursor. Otherwise use the values retained from
 * the previous call in xpos and ypos.
 */
  if(cc->key == KEY_NONE) {
    xpos = 0.5f;
    ypos = 0.5f;
  };
/*
 * Initialize the return descriptor.
 */
  cc->key = KEY_NONE;
  cc->waslow = 0;
  cc->cs = NULL;
  cc->sc = NULL;
  cc->utval = 0.0f;
  cc->clphs = 0.0f;
/*
 * Substitute the cross-hair cursor for the normal cursor is requested.
 */
  if(mode == B_NORM && cp->docross)
    mode = B_CROSS;
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
 * Convert the reference UT and phase to the equivalent NDC position.
 */
      xref = sc->vxa + (xref - sc->utmin) * (sc->vxb - sc->vxa) /
	(sc->utmax - sc->utmin);
/*
 * Get the Y-axis reference value.
 */
      if(csref==NULL) {
	yref = 0.0;
      } else {
	yref = csref->vya   + (yref - cp->phsmin) * (csref->vyb - csref->vya) /
	  (cp->phsmax - cp->phsmin);
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
    cc->waslow = islower((int)key);
    cc->key = cc->waslow ? toupper((int)key) : key;
/*
 * See if the point is in any sub-plot.
 */
    if(xpos >= cp->vxa && xpos <= cp->vxb) {
      cs = cp->cplots;
      for(iplot=0; iplot<cp->nplot && (ypos < cs->vya || ypos > cs->vyb); iplot++,cs++);
      found = iplot<cp->nplot;
    };
/*
 * Was the cursor in a subplot?
 */
    if(found) {
/*
 * Record the details of the sub-plot that the cursor was located in.
 */
      cc->cs = cs;
      cc->iplot = iplot;
/*
 * Record the phase under the cursor.
 */
      cc->clphs = cp->phsmin + (ypos - cs->vya)/(cs->vyb - cs->vya) *
	         (cp->phsmax - cp->phsmin);
/*
 * Identify the scan that the cursor was in and use this to
 * determine the selected UT value.
 */
      for(scan=0; scan<cp->nscan; scan++) {
	sc = &cp->scans[scan];
	if(xpos >= sc->vxa && xpos <= sc->vxb) {
	  cc->utval = sc->utmin + (xpos - sc->vxa)/(sc->vxb - sc->vxa) *
	    (sc->utmax - sc->utmin);
	  cc->sc = sc;
	  break;
	};
      };
    };
    if(!found && noout)
      printf("The cursor must be in one of the plots.\n");
  } while(!found && noout); /* Repeat if outside plots and noout is true */
  return 0;
}

/*.......................................................................
 * Write labels around the frame enclosing all sub-plots.
 *
 * Input:
 *  cp       Clspar *  The plot descriptor.
 * Output:
 *  return      int    0 - OK.
 */
static int c_label(Clspar *cp)
{
  Observation *ob;      /* The descriptor of the observation being plotted */
  Subarray *sub;        /* The descriptor of the sub-array being plotted */
  char awrk[81];        /* Work string for labelling */
  char bwrk[81];        /* Work string for labelling */
/*
 * Get the descriptors of the observation and sub-array being plotted.
 */
  ob = cp->ob;
  sub = cp->sub;
/*
 * Set the viewport around the plot grid.
 */
  cpgsvp(cp->vxa, cp->vxb, cp->cplots[cp->nplot-1].vya, cp->vyb);
/*
 * Compose and write main title.
 */
  cpgsci(1);
  cpgsch(1.0f);
/*
 * Start the title with the source name and date.
 */
  sprintf(awrk, "%s  %s", ob->source.name,
	  sutdate(ob->date.year, ob->date.ut, bwrk));
  cpgmtxt("T", 1.7f, 0.0f, 0.0f, awrk);
  sprintf(awrk, "%s triangles of ", cp->doall ? "Closure":"Upper closure");
/*
 * Describe the fixed telescopes of the triangles.
 */
  if(write_Trispec(ob, &cp->cplots[0].ts, cp->nref, 1, sizeof(bwrk), bwrk) < 0)
    bwrk[0] = '\0';
/*
 * Combine the source/date and telescope description strings if there is
 * room?
 */
  strcat(awrk, (strlen(awrk) + strlen(bwrk) < sizeof(awrk)-1) ? bwrk : "..");
/*
 * Write the IF index part of the title and append it to the title
 * if there is room.
 */
  sprintf(bwrk, " in IF %d", ob->stream.cif+1);
  if(strlen(awrk) + strlen(bwrk) < sizeof(awrk)-1)
    strcat(awrk, bwrk);
  cpgmtxt("T", 0.5f, 0.0f, 0.0f, awrk);
/*
 * In non-interactive mode tell the user what is being plotted.
 */
  if(!cp->docurs) {
    lprintf(stdout, "Page %02.2d: %s of %s\n", cp->npage,
	    cp->doall ? "Triangles":"Upper triangles", bwrk);
  };
/*
 * Write Y label.
 */
  cpgmtxt("L", 3.0f, 0.5f, 0.5f, "Closure phase  (degrees)");
/*
 * Write X labels.
 */
  cpgmtxt("B", 2.5f, 0.5f, 0.5f, "UT");
  return 0;
}

/*.......................................................................
 * Replot the current plots to reflect new attribute selections such as
 * a new UT range. This function should not be called until the first
 * succesful call to c_plot() with oper=C_ALLNEW has been made.
 *
 * Input:
 *  cp        Clspar *  The plot descriptor.
 * Output:
 *  return       int    0 - OK.
 */
static int c_redisp(Clspar *cp)
{
  Clssub *cs; /* Pointer to current sub-plot descriptor */
  int iplot;  /* Number of sub-plot being drawn */
  int ierr=0; /* True if an error occurs */
/*
 * Cursory check of arguments.
 */
  if(cp==NULL) {
    lprintf(stderr, "c_redisp: NULL plot descriptor intercepted\n");
    return -1;
  };
/*
 * Nothing to plot?
 */
  if(cp->nplot<=0) {
    lprintf(stderr, "c_redisp: No plot rows have been initialized.\n");
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
  ierr = ierr || get_Scans(cp)==0;
/*
 * Determine the UT plot range for all plots.
 */
  ierr = ierr || c_utrange(cp);
/*
 * Set up viewport slots for each sub-plot.
 */
  ierr = ierr || c_cpwin(cp, cp->nrow, cp->nplot);
/*
 * Plot each sub-plot.
 */
  for(iplot=0; iplot<cp->nplot && !ierr; iplot++) {
    cs = &cp->cplots[iplot];
    cpgbbuf();
    ierr = ierr || c_plaxes(cp, cs, iplot==0, iplot==cp->nplot-1, 0);
    ierr = ierr || c_pldata(cp, cs, cp->uta, cp->utb, 0);
    ierr = ierr || c_plmodel(cp, cs, 0);
    if(iplot==0) {
      ierr = ierr || c_label(cp);
      ierr = ierr || (cp->docurs && c_mlab(cp, 0));
    };
    cpgebuf();
  };
  return ierr;
}

/*.......................................................................
 * Display a new page of closure triangles.
 *
 * Input:
 *  cp       Clspar *  The plot descriptor.
 *  oper      Clsop    The action to take to display the next page.
 *                      C_ALLNEW  - Start new plot using the given Trispec
 *                                  basis closure-triangle descriptor.
 *                      C_REPLOT  - Re-initialize the current plot to account
 *                                  for attribute changes.
 *                      C_RESET   - Replot from the start/end of the current
 *                                  reference telescope selection.
 *                      C_NEWNUM  - Re-initialize the current plot after
 *                                  a call to  c_setnrow().
 *                      C_NXT_SUB - Skip to the start of the next sub-array.
 *                      C_NXT_TC  - Plot the next page of triangles.
 *                      C_NXT_TRI - As T_NXT_TC, but stop when all the
 *                                  selected triangles have been seen.
 *  forward     int    The order in which to search for plottable closure
 *                     triangles.
 *                       0 - Search in the direction of decreasing
 *                           triangle order.
 *                       1 - Search in the direction of increasing
 *                           triangle order.
 *  refts   Trispec *  If oper==C_ALLNEW this closure-triangle descriptor
 *                     is to be used as the basis for the new plot.
 *                     Otherwise it is ignored and may be NULL.
 * Output:
 *  return      int    The number of sub-plots plotted. Or -1 on
 *                     error.
 */
static int c_plot(Clspar *cp, Clsop oper, int forward, Trispec *refts)
{
  Observation *ob; /* The descriptor of the parent observation */
  Trispec ts;      /* The new trial triangle specification */
/*
 * Check arguments.
 */
  if(cp==NULL) {
    lprintf(stderr, "c_plot: NULL plot descriptor intercepted\n");
    return -1;
  };
/*
 * Get the descriptor of the parent observation.
 */
  ob = cp->ob;
/*
 * The first call must use operation C_ALLNEW.
 */
  if(cp->nplot < 1 && oper != C_ALLNEW) {
    lprintf(stderr, "c_plot: First call must specify C_ALLNEW.\n");
    return -1;
  };
/*
 * Get the descriptor of the first triangle to be plotted.
 */
  switch(oper) {
  case C_ALLNEW:
    if(!refts) {
      lprintf(stderr, "c_plot: NULL basis descriptor received.\n");
      return -1;
    };
    ts = *refts;
/*
 * Record the reference telescope-aggregate type.
 */
    if(ts.nfix < REF_TEL)
      cp->nref = REF_TEL;
    else if(ts.nfix > REF_TRI)
      cp->nref = REF_TRI;
    else
      cp->nref = (Nref) ts.nfix;
/*
 * Locate the first matching triangle.
 */
    if(next_tri(ob, FIND_FIRST, forward, cp->nref, cp->doall, 0, 1, &ts))
      return 0;
    break;
  case C_REPLOT:
    ts = cp->cplots[0].ts;
    break;
  case C_RESET:
    ts = cp->cplots[0].ts;
    if(next_tri(ob, FIND_FIRST, forward, cp->nref, cp->doall, 1, 1, &ts))
      return 0;
    break;
  case C_NXT_SUB:
    ts = cp->cplots[forward ? cp->nplot-1 : 0].ts;
    if(next_tri(ob, SKIP_SUB,forward, cp->nref, cp->doall, 0, 1, &ts))
      return 0;
    break;
  case C_NXT_TC:
    ts = cp->cplots[forward ? cp->nplot-1 : 0].ts;
    if(next_tri(ob, SKIP_TC, forward, cp->nref, cp->doall, 0, 0, &ts) &&
       next_tri(ob, SKIP_TB, forward, cp->nref, cp->doall, 0, 0, &ts) &&
       next_tri(ob, SKIP_TA, forward, cp->nref, cp->doall, 0, 0, &ts) &&
       next_tri(ob, SKIP_SUB,forward, cp->nref, cp->doall, 0, 1, &ts))
      return 0;
    break;
  case C_NXT_TRI:
    ts = cp->cplots[forward ? cp->nplot-1 : 0].ts;
    if(next_tri(ob, FIND_NEXT, forward, cp->nref, cp->doall, 0, 0, &ts))
      return 0;
    break;
  };
/*
 * New sub-array?
 */
  if(cp->sub != ob->sub + ts.isub) {
/*
 * Get the descriptor of the new sub-array.
 */
    cp->sub = &ob->sub[ts.isub];
/*
 * Set up for the full time range of the new sub-array.
 */
    cp->uta = 0;
    cp->utb = cp->sub->ntime-1;
/*
 * Work out the scans of the new sub-array.
 */
    if(get_Scans(cp)==0)
      return -1;
  };
/*
 * Reset up the number of plots per page.
 */
  c_setnrow(cp, cp->nreq);
/*
 * Locate the rest of the triangles sub-array and reference type.
 */
  cp->nplot = 0;
  do {
    cp->cplots[cp->nplot].ts = ts;     /* Record the sub-plot baseline */
    cp->nplot++;                       /* Record addition to sub-plot list */
  } while(cp->nplot < cp->nrow && next_tri(ob, FIND_NEXT, forward, cp->nref,
					   cp->doall, 1, 0, &ts)==0);
/*
 * If we were searching in reverse, the triangles will now be reversed
 * in cplots[]. Rearrange them into forward order.
 */
  if(!forward) {
    Clssub *csa = &cp->cplots[0];
    Clssub *csb = &cp->cplots[cp->nplot-1];
    for( ; csa<csb; csa++,csb--) {
      Trispec nts = csa->ts;
      csa->ts = csb->ts;
      csb->ts = nts;
    };
  };
/*
 * Display the new baselines.
 */
  if(c_redisp(cp))
    return -1;
/*
 * Return the number of baselines plotted.
 */
  return cp->nplot;
}

/*.......................................................................
 * Handle a user request for a new number of plot slots.
 *
 * If the new number is different from the current number, or cp->cplots
 * is NULL then a new array of sub-plot descriptors will be allocated
 * and the new number of rows assigned to cp->nrow.
 *
 * NB. Uses cp->sub.
 *
 * Input:
 *  cp      Clspar *  The plot descriptor.
 *  nreq       int    The requested number of plots per page.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int c_setnrow(Clspar *cp, int nreq)
{
  int ntri=0;/* The number of plottable triangles in the current sub-array */
  int nrow;  /* The new number of rows */
/*
 * If nreq < 1 then signal a request for the defaul number, by assigning
 * nreq=0.
 */
  if(nreq < 1)
    nreq = 0;
/*
 * Determine the number of indenpendant triangles of the current sub-array
 * that contain the given number of fixed telescopes.
 */
  switch(cp->nref) {
  case REF_TRI:
    ntri = 1;   /* There is only one plottable triangle */
    break;
  case REF_BAS:
    ntri = cp->sub->nstat - 2;
    break;
  case REF_TEL:
    {
      long ltmp = cp->sub->nstat - 1;
      ntri = ltmp * (ltmp - 1) / 2;
    };
    break;
  };
/*
 * Start by assuming the requested number is ok, with the reservation that
 * nreq==0 requests the default number 5.
 */
  nrow = nreq > 0 ? nreq : 5;
/*
 * Limit the new number of rows to between 1 and ntri.
 */
  nrow = (nrow > ntri) ? ntri : nrow;
/*
 * Record the new request.
 */
  cp->nreq = nreq;
/*
 * Does this request change the number of sub-plot descriptors allocated?
 */
  if(cp->cplots==NULL || cp->nrow != nrow) {
    size_t nbytes = sizeof(Clssub) * nrow;
/*
 * Attempt to (re-)allocate the sub-plot descriptor array.
 */
    Clssub *cs = cp->cplots ? realloc(cp->cplots, nbytes) : malloc(nbytes);
    if(cs) {
      cp->cplots = cs;
    } else {
      lprintf(stderr, "c_setnrow: Insufficient memory.\n");
      if(cp->cplots==NULL)
	cp->nrow = 0;
      return 1;
    };
/*
 * Record the new number of rows.
 */
    cp->nrow = nrow;
  };
  return 0;
}

/*.......................................................................
 * Determine scaling factors required to convert from world coordinates
 * to mm in the given partition of a given sub-plot.
 *
 * Input:
 *  cp   Clspar *  The plot descriptor.
 *  cs   Clssub *  The sub-plot descriptor.
 * Output:
 *  xtomm float *  ut_shift * *xtomm yields the physical size
 *                 of the ut shift in mm.
 *  ytomm float *  Amplitude_or_phase_shift * *xtomm yields the physical
 *                 size of the amplitude or phase shift in mm.
 *  return  int    0 - OK.
 */
static int c_scale(Clspar *cp, Clssub *cs, float *xtomm, float *ytomm)
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
    lprintf(stderr, "c_scale: No scans visible.\n");
    return -1;
  };
/*
 * Determine the size of the viewport in physical device coordinates
 * (millimeters).
 */
  cpgsvp(sc->vxa, sc->vxb, cs->vya, cs->vyb);
  cpgqvp(2, &xa, &xb, &ya, &yb);
/*
 * Calculate factors to convert world coords into mm.
 */
  *xtomm = fabs((xb-xa)/(sc->utmax-sc->utmin));
  *ytomm = fabs((yb-ya)/(cp->phsmax-cp->phsmin));
  return 0;
}

/*.......................................................................
 * Take a cursor position returned by c_cursor() and locate the index
 * of the closest plotted point.
 *
 * Input:
 *  cp       Clspar *   The plot descriptor.
 *  cc      Clscurs *   The cursor descriptor.
 * Output:
 *  return      int     The integration index of the nearest point or
 *                      -1 if there is no displayed data in the zone
 *                      where the cursor was pressed.
 */
static int c_find(Clspar *cp, Clscurs *cc)
{
  Clssub *cs;        /* Subplot descriptor selected */
  Subarray *sub;     /* The descriptor of the displayed sub-array */
  Integration *integ;/* Descriptor of an integration */
  int ut;            /* The integration index being checked */
  double vlbut;      /* The ut selected */
  float phs;         /* The closure phase under the cursor */
  float xtomm,ytomm; /* Conversion factor between world coords and mm */
  int bestut=0;      /* The ut index of the point closest to the cursor */
  float mindist=0.0f;/* The min value of 'dist' */
  int first=1;       /* True until end of first iteration of search loop */
/*
 * Get the descriptor of the sub-plot selected by the user.
 */
  cs = cc->cs;
/*
 * Cursor pressed outside of any plot?
 */
  if(cs==NULL)
    return -1;
/*
 * Determine conversion factors from world coords to mm.
 */
  if(c_scale(cp, cs, &xtomm, &ytomm))
    return -1;
/*
 * Calculate the ut corresponding to the cursor selected time.
 */
  vlbut = cc->utval + cp->utref;
  phs = cc->clphs;
/*
 * Get the descriptor of the currently displayed sub-array.
 */
  sub = cp->sub;
/*
 * Locate the nearest point.
 */
  for(ut=cp->uta, integ = &sub->integ[ut]; ut<=cp->utb; ut++,integ++) {
    Clphs *cphs = get_clphs(&cs->ts, integ->vis);
/*
 * Skip deleted points. Also skip flagged data if not displayed.
 */
    if(!(cphs->bad & FLAG_CDEL) && (!cphs->bad || cp->doflag)) {
      float xdif = xtomm * (integ->ut - vlbut);
      float ydif = ytomm * (phs - cphs->ophs);
      float dist = xdif*xdif + ydif*ydif;
/*
 * Compare the squared distance from this point with the that of
 * the previous closest point.
 */
      if(first || dist < mindist) {
	first = 0;
	bestut = ut;
	mindist = dist;
      };
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
 * Plot an extra mode label for editting sessions.
 *
 * Input:
 *  cp       Clspar *  The plot descriptor.
 *  erase       int    If true, erase existing mode label.
 * Output:
 *  return      int    0 - OK.
 */
static int c_mlab(Clspar *cp, int erase)
{
  Observation *ob;  /* The descriptor of the observation being plotted */
  int oldcol;       /* Temporary storage for entry color index */
  char label[81];   /* Temporary work string to compose mode label in */
/*
 * Get the descriptor of the observation.
 */
  ob = cp->ob;
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
  cpgsvp(cp->vxa, cp->vxb, cp->vya, cp->vyb);
/*
 * Compose the mode label.
 */
  sprintf(label, "%s editing of %s channels of %s.",
	  (cp->tri_ed || cp->nref==REF_TRI) ? "Triangle" :
	     (cp->nref==REF_BAS ? "Baseline" : "Station"),
	  cp->ch_ed ? "selected" : "all",
	  cp->if_ed ? "the displayed IF" : "all IFs");
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
 *  cp       Clspar *  The plot descriptor.
 *  tri_ed      int    Ask for triangle (3 baseline) editing.
 *  if_ed       int    Select IF based editing if true.
 *  ch_ed       int    Select channel based editing if true.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int c_newmode(Clspar *cp, int tri_ed, int if_ed, int ch_ed)
{
/*
 * Buffer until the new text has been plotted.
 */
  cpgbbuf();
/*
 * Erase the existing mode line.
 */
  c_mlab(cp, 1);
/*
 * Install the new editing modes.
 */
  cp->if_ed = if_ed;
  cp->ch_ed = ch_ed;
  cp->tri_ed = tri_ed;
/*
 * Only allow baseline-based editing when plotting against a reference
 * baseline.
 */
/*
 * Draw the new mode line.
 */
  c_mlab(cp, 0); /* Plot new mode line */
/*
 * reveal the changes.
 */
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Edit one point (by baseline or triangle) and redisplay on each
 * respective displayed sub-plot.
 *
 * Input:
 *  cp       Clspar *  The plot descriptor.
 *  cs       Clssub *  The sole triangle to be editted if cp->base_ed==0.
 *  flag        int    If true, flag unflagged data. If false, restore
 *                     flagged data.
 *  ut          int    The index of the integration to be editted.
 * Output:
 *  return      int    0 - OK.
 */
static int c_edit(Clspar *cp, Clssub *cs, int flag, int ut)
{
  Observation *ob;  /* The descriptor of the observation being edited */
  Subarray *sub;    /* The descriptor of the displayed sub-array */
  int iplot;        /* Sub-plot index */
  int cif;          /* The index of the IF being edited */
  int ierr=0;       /* Error status */
/*
 * Check descriptors.
 */
  if(cp==NULL || cs==NULL) {
    lprintf(stderr, "c_edit: NULL %s descriptor intercepted.\n",
	    cp==NULL?"plot":"sub-plot");
    return 1;
  };
/*
 * Get the descriptor of the observation and sub-array being edited.
 */
  ob = cp->ob;
  sub = cp->sub;
/*
 * Get the index of the IF being edited.
 */
  cif = ob->stream.cif;
/*
 * This function modifies the data.
 */
  cp->modified = 1;
/*
 * Start by erasing the given integration from all sub-plots.
 */
  cpgbbuf();  /* Buffer changes */
  for(iplot=0; !ierr && iplot<cp->nplot; iplot++)
    ierr = c_pldata(cp, &cp->cplots[iplot], ut, ut, 1);
/*
 * Now edit the point on all relevant baselines.
 */
  if(!ierr) {
    Trispec *ts = &cs->ts;
    Nref nref = cp->tri_ed ? REF_TRI : cp->nref;
/*
 * Edit the three triangle baselines?
 */
    switch(nref) {
    case REF_TRI:
      ierr = ierr || ed_integ(ob, sub, ut, cif, flag, 1, 0,
			      cp->ch_ed, cp->if_ed, ts->b[0].base);
      ierr = ierr || ed_integ(ob, sub, ut, cif, flag, 1, 0,
			      cp->ch_ed, cp->if_ed, ts->b[1].base);
      ierr = ierr || ed_integ(ob, sub, ut, cif, flag, 1, 0,
			      cp->ch_ed, cp->if_ed, ts->b[2].base);
      break;
/*
 * Edit the reference baseline?
 */
    case REF_BAS:
      ierr = ierr || ed_integ(ob, sub, ut, cif, flag, 1, 0,
			      cp->ch_ed, cp->if_ed, ts->b[0].base);
      break;
/*
 * Edit the reference station?
 */
    case REF_TEL:
      ierr = ierr || ed_integ(ob, sub, ut, cif, flag, 0, 1,
			      cp->ch_ed, cp->if_ed, ts->ta);
      break;
    };
  };
/*
 * Re-plot the given integration on all sub-plots.
 */
  for(iplot=0; !ierr && iplot<cp->nplot; iplot++)
    ierr = c_pldata(cp, &cp->cplots[iplot], ut, ut, 0);
/*
 * End pgplot buffering.
 */
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Given the values returned by c_cursor(), toggle the flagged status of
 * the integration closest to the cursor.
 *
 * Input:
 *  cp       Clspar *  The plot descriptor.
 *  cc      Clscurs *  The cursor descriptor.
 * Output:
 *  return      int    0 - OK.
 */
static int c_toggle(Clspar *cp, Clscurs *cc)
{
  int ut;           /* The index of the nearest integration */
  int flagged;      /* True if nearest visibility is flagged */
/*
 * Find the nearest integration on the specified baseline.
 */
  ut = c_find(cp, cc);
  if(ut < 0)
    return 0; /* Ignore failed locations */
/*
 * Is the closure phase flagged?
 */
  flagged = get_clphs(&cc->cs->ts, cp->sub->integ[ut].vis)->bad & FLAG_BAD;
/*
 * Toggle the status of the relevant baselines associated with the closure
 * sample and propogate effects into other effected subplots.
 */
  return c_edit(cp, cc->cs, !flagged, ut);
}

/*.......................................................................
 * Allow a range box to be selected and optionally either flag all points
 * above and below the box, or restore all points inside the box.
 *
 * Input:
 *  cp      Clspar *  The visibility plot attributes to supply to visplt().
 *  doflag     int    If true flag data within box. If false unflag
 *                    everything within the box.
 * Output:
 *  return     int    0 - OK. Anything else means fatal error.
 */
static int c_box(Clspar *cp, int doflag)
{
  Integration *integ; /* Descriptor of an integration in the select box */
  int ut;           /* UT index loop variable */
  Clssub *cs=NULL;  /* Plot in which cursor was pressed */
  double utmin=0.0; /* Start of the selected ut range */
  double utmax=0.0; /* End of the selected ut range */
  float utrefc=0.0f;  /* Cursor reference UT */
  float phsrefc=0.0f; /* Cursor reference phase */
  float minphs=0.0; /* Lower phase bound of the select box */
  float maxphs=0.0; /* Higher phase bound of the select box */
  int iter;         /* Iterate over getting two valid keypresses */
/*
 * Get the first cursor position.
 */
  for(iter=0; iter<2; iter++) {
    int accepted = 0;
    while(!accepted) {
      if(c_cursor(cp, 1, iter==0 ? B_NORM:B_RECT, cs, utrefc, phsrefc,
		  doflag ? cutcol:rescol))
	return 1;
      switch(cp->cursor.key) {
      case KEY_QUIT: case KEY_CAN:     /* Abort box selection */
	return 0;
	break;
      case KEY_CUR:             /* Accept the selected vertex */
	if(iter==0) {
	  utrefc = cp->cursor.utval;
	  phsrefc = cp->cursor.clphs;
	  utmin = utmax = utrefc;
	  minphs = maxphs = phsrefc;
	  cs = cp->cursor.cs;
	} else {
	  if(cp->cursor.cs != cs) {
	    fprintf(stderr, "Select box spans more than one plot.\n");
	    return 0;
	  };
	  if(cp->cursor.utval < utmin)
	    utmin = cp->cursor.utval;
	  else
	    utmax = cp->cursor.utval;
	  if(cp->cursor.clphs < minphs)
	    minphs = cp->cursor.clphs;
	  else
	    maxphs = cp->cursor.clphs;
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
    };
  };
/*
 * Get the ut values.
 */
  utmin += cp->utref;
  utmax += cp->utref;
/*
 * Buffer PGPLOT instructions while the points are redrawn.
 */
  cpgbbuf();
/*
 * Only flag points in the UT range selected.
 */
  for(ut=cp->uta,integ = &cp->sub->integ[ut];  ut<=cp->utb;  ut++,integ++) {
    double vlbut = integ->ut;
    if(vlbut >= utmin && vlbut <= utmax) {  /* Only flag in given UT range */
      Clphs *cphs = get_clphs(&cs->ts, integ->vis);
/*
 * Only consider visible points.
 */
      if(!(cphs->bad & FLAG_CDEL) && (!cphs->bad || cp->doflag)) {
/*
 * Determine whether the current point is inside or outside the
 * select-box.
 */
	int inside = cphs->ophs >= minphs && cphs->ophs <= maxphs;
/*
 * Edit and redisplay the point in all relevant sub-plots.
 */
	if(inside && c_edit(cp, cs, doflag, ut))
	  return 1;
      };
    };
  };
  cpgebuf(); /* Plotting complete - end plot buffering */
  return 0;
}

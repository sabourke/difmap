#ifndef vplot_h
#define vplot_h

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
 * Elements of the following type are used to associate sorted times,
 * in the currently selected time system, with integrations.
 */
typedef struct {
  Integration *integ; /* The integration of this time sample */
  float t;            /* The time in the form used in plotting the X-axis */
} TimeSample;

/* Declare a container for sub-plot information */

typedef struct {
  float ampmin;   /* Min amplitude */
  float ampmax;   /* Max amplitude */
  float vya,vyb;  /* Min,max NDC Y-coords of viewport */
  float vymid;    /* NDC Y-coord where amp and phase plots touch */
  int base;       /* The baseline in this plot. */
} Vissub;

/* Declare a container for cursor selection details */

typedef struct {
  int key;       /* The upper-case value of the key used to return the cursor */
  int waslow;    /* True if 'key' was lower case */
  int wasamp;    /* True if the cursor was pressed in an amplitude plot */
  Vissub *vs;    /* Descriptor of the plot in which the cursor was pressed */
  int iplot;     /* The index of 'vs' in vp->vplots[] */
  Scan *scan;    /* The descriptor of the scan under the cursor */
  float tval;    /* The time value under the cursor */
  float value;   /* The amplitude or phase value under the cursor */
} Vcurs;

typedef struct {
  double utref;      /* Reference UT of observation (seconds) */
  double stref;      /* The apparent sidereal time at utref (seconds) */
  Observation *ob;   /* The descriptor of the observation being plotted */
  Subarray *sub;     /* The descriptor of sub-array refsub in ob->sub[] */
  Basespec bs_beg;   /* Baseline spec of the first plotted baseline */
  Basespec bs_end;   /* Baseline spec of the last plotted baseline */
  TimeSample *times; /* The time samples in X-axis plot order */
  int times_stale;   /* True if the times in times[] or the scans in 'scans' */
                     /*  are out of date. */
  FreeList *scan_mem;/* Memory for allocating Scan elements */
  Scan *scans;       /* List of scans */
  int scans_stale;   /* True if the scan list needs to be updated */
  float wxa;         /* World min X coordinate (60th seconds wrt utref) */
  float wxb;         /* World max X coordinate (60th seconds wrt utref) */
  float phsmin;      /* Minimum phase to plot (radians) */
  float phsmax;      /* Maximum phase to plot (radians) */
  float ampmin;      /* Minimum amplitude to plot (only if ampmax>ampmin) */
  float ampmax;      /* Maximum amplitude to plot (only if ampmax>ampmin) */
  float vxa,vxb;     /* Viewport surrounding grid of sub-plots */
  float vya,vyb;
  int modified;      /* This remains 0 unless the data are edited */
  int stat_ed;       /* If true, editing is station based, else baseline */
                     /*  based. */
  int if_ed;         /* If true, edits are restricted to the current IF */
  int ch_ed;         /* If true, edits are restricted to current freq */
                     /*  channels. */
  int ta,tb;         /* The indexes of the first and last plotted */
                     /*  integrations. */
  int docurs;        /* True when cursor control is in effect */
  int doamp;         /* If true then the amplitude will be plotted */
  int dophs;         /* If true then the phase will be plotted */
  int doflag;        /* True if flagged points are or will be plotted */
  int domod;         /* True if model lines are or will be plotted */
  int dobars;        /* True if error bars are or will be plotted */
  int docross;       /* True to enable cross-hair mode */
  int doutc;         /* True to plot UTC along the X-axis, false for sidereal */
                     /*  time. */
  int doall;         /* True to see all baselines of each reference telescope */
  int showall;       /* When true display full range of data+model */
  int dodiff;        /* If true, show the difference between the data and the */
                     /*  model. */
  int nrow;          /* The number of sub-plot slots on the display */
  int nplot;         /* The actual number of sub-plots plotted */
  int nreq;          /* The requested number of plots per page */
  int maxplot;       /* Number of slots available in 'vplots' */
  Vissub *vplots;    /* Dynamic array of sub->nstat subplot descriptors */
  int doscan;
  Vcurs cursor;      /* The descriptor of the last cursor selection */
  int npage;         /* The sequential number of the page being plotted */
  int old_if;        /* The index of the IF to restore on exit */
} Vedpar;

Vedpar *new_Vedpar(Observation *ob, int cif, int docurs, int doscan, int doamp,
		   int dophs, int doflag, int domod, int dobars, int showall,
		   int nrow);
Vedpar *del_Vedpar(Vedpar *vp);

int v_update_times(Vedpar *vp);
int v_update_scans(Vedpar *vp);

int v_pldata(Vedpar *vp, Vissub *vs, int ta, int tb, int erase);
int v_plmodel(Vedpar *vp, Vissub *vs, int erase);

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

int v_cursor(Vedpar *vp, int noout, Bandmode mode, int isamp, Vissub *vsref,
		    float xref, float yref, int ci);

int v_label(Vedpar *vp);
int v_newmode(Vedpar *vp, int stat_ed, int if_ed, int ch_ed);
int v_redisp(Vedpar *vp);

/* Define an enum used to specify what v_plot should plot next */

typedef enum {
  V_ALLNEW,     /* Start new plot using the given reftel, base and isub */
  V_REPLOT,     /* Replot the current page */
  V_RESET,      /* Reset to the first/end baseline of the current telescope */
  V_NXTSUB,     /* Start plotting from the next sub-array */
  V_NXT_TA,     /* Start plotting for the next reference telescope */
  V_NXT_TB,     /* Start plotting the next page of baselines. */
  V_NEXT        /* As V_NXT_TB, but stop when all selected baselines */
                /* have been seen. */
} Vedop;

int v_plot(Vedpar *vp, Vedop oper, int forward, Basespec *init);
int v_scale(Vedpar *vp, Vissub *vs, int doamp, float *xtomm, float *ytomm);
int v_setnrow(Vedpar *vp, int nreq);
void v_data_point(Vedpar *vp, Visibility *vis, float *amp, float *phs);
int v_toggle_timesys(Vedpar *vp);

#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "logio.h"
#include "mapmem.h"
#include "mapwin.h"
#include "units.h"
#include "vlbutil.h"
#include "vlbconst.h"
#include "vlbmath.h"
#include "maplot.h"
#include "cpgplot.h"

/* Set the keys used during cursor interaction */
 
enum {
  KEY_CORN   = 'A', /* Select start or opposite end vertex of clean window */
  KEY_DEL    = 'D', /* Delete the window who's corner is nearest the cursor */
  KEY_FIDL   = 'F', /* Fiddle the color table contrast and brightness */
  KEY_QUIT   = 'X', /* Quit interactive display */
  KEY_STAT   = 'S', /* Show stats of the window who's corner is nearest the */
                    /*  cursor. */
  KEY_TRAN   = 'T', /* Change color-table transfer function log <-> linear */
  KEY_COL    = 'C', /* Show map in pseudo-color */
  KEY_GRY    = 'G', /* Show map in grey-scale */
  KEY_KEEP   = 'K', /* Keep display limits for next invokation of maplot */
  KEY_DISP   = 'L', /* Re-display the plot */
  KEY_MOD    = 'M', /* Toggle display of the model */
  KEY_CMP    = 'N', /* Model component selection prefix */
  KEY_REM    = 'R', /* Key to select component to be removed */
  KEY_VAL    = 'V', /* Report the value of the pixel under the cursor */
  KEY_ZOOM   = 'Z', /* Key to introduce subimage selection */
  KEY_HELP   = 'H', /* Display the list of available keys */
  KEY_CROSS  = '+', /* Toggle cross-hair cursor mode */
  KEY_UNMARK = 'U'  /* Remove the marker nearest the cursor */
};

/*
 * Define the gaps (measured in character heights) between labels on each axis.
 */
static const float topsep =0.7f; /* Separation of title line from top axis */
static const float primsep=2.5f; /* Separation of axis labels from axis */
static const float clevsep=0.5f; /* clevsep+sepinc is separation of contour */
                                 /* text from X-axis label */
static const float sepinc=0.3f;  /* Separation between extra anotation lines */
static const float margin=0.01f; /* Margin around anotation text in NDC */
                                 /* smallest dimension of view surface */
static const float wdginc=0.2f;  /* Gap above wedge */
static const float wdgsiz=3.5f;  /* Height of wedge in char's */
static const int wincol=10;      /* The color of the window rectangles */
static const int zoomcol=5;      /* The color of zoom rectangles */
static const float echsize=0.8f; /* Char height of extra anotation text */

#define PORTWID 40  /* Width of the viewport in characters */

/* Define a structure to contain contour plot parameters */

typedef struct {
  float peak;            /* Peak value in image to be contoured */
  float cmin,cmax;       /* Min/max values in sub-image for contouring */
  float cmul;            /* Contour level scale factor */
  float *levs;           /* Contour level array */
  int nlev;              /* Number of levels in 'levs' */  
  int plevs;             /* If true then anotate with percentage levels */
} Contour;

/* Define a structure to contain color-map false-color plot parameters */

typedef struct {
  Ctable *ctab;   /* The color table */
  float vmin;     /* Min absolute false-color level */
  float vmax;     /* Max absolute false-color level */
} Cmpar;

/* Define a plot descriptor */

typedef struct {
  Observation *ob;       /* UV data container/descriptor */
  MapBeam *mb;           /* Pointer to map/beam container */
  Mapwin *mw;            /* Clean window list */
  MaplotBeam *mpb;       /* Symbolic beam display parameters */
  MaplotVect *vect;      /* Polarization vector display attributes */
  float *box;            /* Pointer to array of input displayed area limits */
  Model *newmod;         /* Way-station container for new model components */
  int hard;              /* True if plot device is hard-copy */
  int mono;              /* True if the plot device is effectively monochrome */
  int page;              /* Current page number */
  int cursor;            /* True if plot device has a cursor */
  int docross;           /* True to enable cross-hair mode */  
  int dowin;             /* True if windows should be plotted */
  int domod;             /* If true, display the model components */
  int dovar;             /* If true, display just variable model components */
  int docont;            /* If true plot contours */
  int dovect;            /* If true plot polarization vectors */
  int domap;             /* If true, the base plot is of the map */
  float *image;          /* Pointer to image to be plotted */
  Contour cpar;          /* Contour plot parameters */
  Cmpar cmpar;           /* Color-map parameters */
  int   pxa,pxb,pya,pyb; /* Range of pixels to be displayed */
  float wxa,wxb,wya,wyb; /* Range of world-coordinates */
  float xtomm,ytomm;     /* Conversion factors world coord -> millimeters */
  float tr[6];           /* Coord mapping matrix for pggrey,pgimag and pgcont */
  MarkerList *markers;   /* A list of markers to be displayed in the plot */
} Maplot;

/* Type used to record cursor selection info */

typedef struct {
  float x,y;      /* Cursor selected position in world coords */
  int key;        /* Upper-case version of the selection key pressed */
  int waslow;     /* True if the enterred key was given in lower case */
} Keypos;

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

static Maplot *new_Maplot(Observation *ob, MapBeam *mb, Mapwin *mw,
			  MaplotBeam *mpb, MaplotVect *vect, int domap,
			  Ctable *ctab, int docont, int dovect, int domod,
			  float *levs, int nlev, float cmul, float *box,
			  MarkerList *markers);
static Maplot *del_Maplot(Maplot * mp);
static int bad_Maplot(char *name, Maplot *mp);
static int setarea(Maplot *mp, float xmin, float xmax, float ymin, float ymax);
static int setcont(Maplot *mp, float *levs, float cmul, int nlev);
static int setcmpar(Maplot *mp, Ctable *ctab);
static char *lev_text(Maplot *mp, int *ndone);
static int setport(Maplot *mp);
static int replot(Maplot *mp);
static int plcont(Maplot *mp);
static int plvect(Maplot *mp);
static int plimage(Maplot *mp);
static int pllabel(Maplot *mp);
static int plmodel(Maplot *mp);
static int plmarkers(Maplot *mp);
static void plsmvp(float lftgap, float rgtgap, float botgap, float topgap,
		    float  chhgt);
static void plqcd(int units, int ishoriz, float *chsize);
static int plwedge(char *side, float disp, float width, float fg, float bg,
		   char *label, int mode, int doimag);
static Subwin *findwin(Maplot *mp, float xpos, float ypos);
static Subwin *win_limits(Maplot *mp, Subwin *win);
static int draw_win(Maplot *mp, Subwin *win, int erase);
static int plwins(Maplot *mp);
static int get_curs(Maplot *mp, int first, Bandmode mode,
		    float xref, float yref, int ci, Keypos *kp);
static int zap_win(Maplot *mp, Keypos *kp);
static int set_win(Maplot *mp, Keypos *kp);
static int interact(Maplot *mp);
static int set_zoom(Maplot *mp);

typedef struct {
  float x,y;     /* The world coordinates of the nearest pixel center */
  float value;   /* The value held by the nearest pixel */
  float poli;    /* The polarized intensity if mp->dovect is true */
  float pola;    /* The polarized angle (radians) if mp->dovect is true */
} Pixval;

static Pixval *pix_val(Maplot *mp, Keypos *kp);
static int make_cmp(Maplot *mp);
static int keep_cmp(Maplot *mp);

typedef struct {
  Model *mod;    /* Model containing the closest component */
  Modcmp *prev;  /* Component that precedes the closest component in mod */
  Modcmp *cmp;   /* The descriptor (prev->next) of the closest component */
  float roff;    /* Radial offset of component in map coordinates */
} Cmpfnd;

Cmpfnd *fnd_cmp(Maplot *mp, Keypos *kp);
static void cf_search(Cmpfnd *cf, Keypos *kp, Model *mod, int dofix, int dovar);
static int zapcmp(Maplot *mp, Keypos *kp);
static MarkerNode *find_marker(Maplot *mp, Keypos *kp);
static int zapmark(Maplot *mp, Keypos *kp);
static int change_cmap(Maplot *mp, Cmclass class, int ask);
static int change_transfer(Maplot *mp, Cmtran tran, int ask);

#define BASELEV 17   /* Lowest color index that PGGRAY will use. */
#define MINLEVS 15   /* Min number of gray-scale levels used by PGGRAY */
#define MAXLEVS 127  /* Max number of gray-scale levels used by PGGRAY */

/*.......................................................................
 * Set up the parameters for subsequent contour plotting.
 * Store the results in mp->cpar.
 *
 * Input/Output:
 *  mp    Maplot *   The Plot descriptor.
 *                   On return, mp->cpar will have been initialized.
 * Input:
 *  levs   float *   Array of 'nlev' contour levels, use NULL for default.
 *  cmul   float     Contour multiplier in map units. If zero then
 *                   1% level is assumed.
 *  nlev     int     The number of levels in 'levs' array.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Error.
 */
static int setcont(Maplot *mp, float *levs, float cmul, int nlev)
{
  Contour *cpar;         /* Pointer to mp->cpar */
  int   pxa,pxb,pya,pyb; /* Range of map pixels */
  float immin,immax;     /* Min,max values in image array */
/*
 * Define the default percentage contour levels for the map array.
 */
  static float mlevs[]={-1.0f, 1.0f, 2.0f, 4.0f, 8.0f, 16.0f, 32.0f, 64.0f};
  static int mnum=sizeof(mlevs)/sizeof(float);
/*
 * Define the default percentage contour levels for the beam.
 */
  static float blevs[]={-64.0f, -16.0f, -4.0f, 4.0f, 16.0f, 32.0f, 64.0f};
  static int bnum=sizeof(blevs)/sizeof(float);
/*
 * Check the Maplot descriptor.
 */
  if(bad_Maplot("setcont", mp))
    return 1;
/*
 * Get a local pointer to the contour descriptor.
 */
  cpar = &mp->cpar;
/*
 * Determine the displayable image area.
 */
  pxa = mp->mb->nx/4;
  pya = mp->mb->ny/4;
  pxb = 3*pxa - 1;
  pyb = 3*pya - 1;
/*
 * Record the given or default percentage contour levels.
 */
  if(levs==NULL || nlev<=0) {
    if(mp->domap) {     /* Select appropriate defaults */
      cpar->levs = &mlevs[0];
      cpar->nlev = mnum;
    } else {
      cpar->levs = &blevs[0];
      cpar->nlev = bnum;
    };
  } else {                   /* Use given levels */
    cpar->levs = levs;
    cpar->nlev = nlev;
  };
/*
 * Determine the range of values in the image.
 */
  imran(mp->image, mp->mb->nx, mp->mb->ny, pxa, pxb, pya,
	pyb, &immin, &immax);
/*
 * Record the maximum.
 */
  cpar->peak = fabs(immax) > fabs(immin) ? immax : immin;
/*
 * If no contour multiplier was given for the image, determine the
 * percentage multiplier for the image.
 */
  cpar->plevs = cmul < 1.0e-10;   /* Percentage levels */
  cpar->cmul = cpar->plevs ? (cpar->peak/100.0f) : cmul;
/*
 * Until we know any better, assume that the whole map area is available
 * for contouring.
 */
  cpar->cmin = immin;
  cpar->cmax = immax;
  return 0;
}

/*.......................................................................
 * Determine the false-color levels for the whole area of the image to
 * be imaged. Store the results in mp->cmpar.
 *
 * Input/Output:
 *  mp      Maplot *   The Plot descriptor.
 *                     On return, mp->cmpar will have been initialized.
 * Input:
 *  ctab    Ctable *   The descriptor of the color-table and transfer
 *                     function.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Error.
 */
static int setcmpar(Maplot *mp, Ctable *ctab)
{
  Cmpar *cmpar;          /* Local copy of mp->cmpar */
  int   pxa,pxb,pya,pyb; /* Range of map pixels */
  float immin,immax;     /* Min,max values in map array */
/*
 * Check the Maplot descriptor.
 */
  if(bad_Maplot("setcmpar", mp))
    return 1;
/*
 * Get a pointer to the color-map descriptor.
 */
  cmpar = &mp->cmpar;
/*
 * Record the color-table descriptor.
 */
  if(ctab) {
    cmpar->ctab = ctab;
  } else {
    lprintf(stderr, "setcmpar: NULL color-table descriptor intercepted.\n");
    return 1;
  };
/*
 * Determine the displayable map area.
 */
  pxa = mp->mb->nx/4;
  pya = mp->mb->ny/4;
  pxb = 3*pxa - 1;
  pyb = 3*pya - 1;
/*
 * Determine the range of values in the image.
 */
  imran(mp->image, mp->mb->nx, mp->mb->ny, pxa, pxb, pya,
	pyb, &immin, &immax);
/*
 * Determine absolute psuedo-color levels for the image.
 */
  if(fabs(ctab->vmin - ctab->vmax) < 1.0e-15) {
    cmpar->vmin = immin;
    cmpar->vmax = immax;
  } else {
    cmpar->vmin = ctab->vmin;
    cmpar->vmax = ctab->vmax;
  };
  return 0;
}

/*.......................................................................
 * Set the parameters in the plot descriptor that describe the displayed
 * area of the map/beam plot.
 * NB.  setcont() must be called before this function.
 *
 * The following members of the plot descriptor will be modified by this
 * function.
 *    pxa,pxb,pya,pyb  -  Range of pixels to be displayed.
 *    wxa,wxb,wya,wyb  -  Range of world-coordinates.
 *    tr               -  PGPLOT coordinate transformation matrix.
 *    cpar.cmin        -  Min value in contour sub-image.
 *    cpar.cmax        -  Max value in contour sub-image.
 *
 * Input:
 *  mp     Maplot *   The map-plot descriptor.
 *  xmin    float     The minimum X-axis position displayed (radians).
 *                    To select the whole X-axis plot area send xmin=xmax=0.0f.
 *  xmax    float     The maximum X-axis position displayed (radians).
 *  ymin    float     The minimum Y-axis position displayed (radians).
 *                    To select the whole Y-axis plot area send ymin=ymax=0.0f.
 *  ymax    float     The maximum Y-axis position displayed (radians).
 */
static int setarea(Maplot *mp, float xmin, float xmax, float ymin, float ymax)
{
  Contour *cpar; /* Contour plot descriptor */
  MapBeam *mb;   /* The map/beam descriptor. */
  int xcent;     /* The central X-axis pixel in the map/beam arrays */
  int ycent;     /* The central Y-axis pixel in the map/beam arrays */
  float xinc;    /* The X-axis pixel-increment measured in radians */
  float yinc;    /* The Y-axis pixel-increment measured in radians */
  float wxa;     /* Float pixel dist from map center to window X min */
  float wxb;     /* Float pixel dist from map center to window X max */
  float wya;     /* Float pixel dist from map center to window Y min */
  float wyb;     /* Float pixel dist from map center to window Y max */
  int xa,xb;     /* Start,end pixels along X-axis */
  int ya,yb;     /* Start,end pixels along Y-axis */
  int ixmin,ixmax; /* X-axis pixel limits imposed by map/beam array */
  int iymin,iymax; /* Y-axis pixel limits imposed by map/beam array */
  int incr;        /* Increment in pixel limits */
  float ftmp;
/*
 * Check the validity of the Maplot descriptor.
 */
  if(bad_Maplot("setarea", mp))
    return 1;
/*
 * Enforce xmin<xmax and ymin<ymax.
 */
  if(xmin>xmax) {ftmp=xmin; xmin=xmax; xmax=ftmp;};
  if(ymin>ymax) {ftmp=ymin; ymin=ymax; ymax=ftmp;};
/*
 * Get a pointer to the map/beam descriptor.
 */
  mb = mp->mb;
/*
 * Determine the central pixels of each axis.
 */
  xcent = mb->nx/2;
  ycent = mb->ny/2;
/*
 * Determine the pixel limits of the displayable map area.
 */
  ixmin = mb->nx/4;
  iymin = mb->ny/4;
  ixmax = 3*ixmin - 1;
  iymax = 3*iymin - 1;
/*
 * Use defaults?
 */
  if(xmin==xmax || ymin==ymax) {
    xa = ixmin;
    xb = ixmax;
    ya = iymin;
    yb = iymax;
  } else {
/*
 * Convert the window coordinates to numbers of pixels wrt the centre
 * of the map/beam.
 */
    wxa = xmin/mb->xinc;
    wxb = xmax/mb->xinc;
    wya = ymin/mb->yinc;
    wyb = ymax/mb->yinc;
/*
 * Convert the ranges from radians to map pixel element numbers in such
 * a way that all pixels whose centres are enclosed by the window
 * are included.
 */
    xa = xcent + (int) (wxa + ((wxa<0)?0.0:1.0));
    xb = xcent + (int) (wxb - ((wxb<0)?1.0:0.0));
    ya = ycent + (int) (wya + ((wya<0)?0.0:1.0));
    yb = ycent + (int) (wyb - ((wyb<0)?1.0:0.0));
/*
 * Enforce limits.
 */
    if(xa < ixmin)
      xa = ixmin;
    if(ya < iymin)
      ya = iymin;
    if(xb > ixmax)
      xb = ixmax;
    if(yb > iymax)
      yb = iymax;
/*
 * Symmetrically expand the range if less than two pixels were requested
 * along one or both axes. NB. It is possible to get xb<xa if no pixel
 * centers lay inside the selection window - the abs(xb-xa) in the increment
 * takes care of this.
 */
    if(xb-xa < 1) {
      incr = abs(xb-xa) + 1;  /* The increment to get a 3 pixel range */
      if(xa>ixmin)
	xa -= incr;
      if(xb<ixmax)
	xb += incr;
    };
    if(yb-ya < 1) {
      incr = abs(yb-ya) + 1;  /* The increment to get a 3 pixel range */
      if(ya>iymin)
	ya -= incr;
      if(yb<iymax)
	yb += incr;
    };
  };
/*
 * Record the new pixel limits.
 */
  mp->pxa = xa;
  mp->pxb = xb;
  mp->pya = ya;
  mp->pyb = yb;
/*
 * Determine the world-coordinate limits that these pixels imply.
 */
  xinc = mb->xinc;
  if(xinc > 0) {
    mp->wxa = (xa - xcent) * xinc;
    mp->wxb = (xb - xcent) * xinc;
  } else {
    mp->wxa = (xb - xcent) * xinc;
    mp->wxb = (xa - xcent) * xinc;
  };
  yinc = mb->yinc;
  if(yinc > 0) {
    mp->wya = (ya - ycent) * yinc;
    mp->wyb = (yb - ycent) * yinc;
  } else {
    mp->wya = (yb - ycent) * yinc;
    mp->wyb = (ya - ycent) * yinc;
  };
/*
 * Set up the PGPLOT coordinate transformation matrix.
 */
  mp->tr[0] = -xinc*(xcent+1);
  mp->tr[1] = xinc;
  mp->tr[2] = 0.0f;
  mp->tr[3] = -yinc*(ycent+1);
  mp->tr[4] = 0.0f;
  mp->tr[5] = yinc;
/*
 * Determine the min/max possible contour levels within the new area.
 */
  cpar = &mp->cpar;
  imran(mp->image, mp->mb->nx, mp->mb->ny, xa, xb, ya, yb,
	&cpar->cmin, &cpar->cmax);
  return 0;
}

/*.......................................................................
 * Test if a Maplot descriptor is NULL. If so, write an error message
 * of the form:
 *  lprintf(stderr, "%s: NULL Maplot descriptor intercepted\n", name);
 * and return true.
 *
 * Input:
 *  name   char *  The name of the calling funtion.
 *  mp   Maplot *  The descriptor pointer to be tested.
 * Output:
 *  return  int    mp==NULL.
 */
static int bad_Maplot(char *name, Maplot *mp)
{
  if(mp==NULL)
    lprintf(stderr, "%s: NULL Maplot descriptor intercepted\n", name);
  return mp==NULL;
}

/*.......................................................................
 * Provide interactive display of a map or beam.
 *
 * Input:
 *  ob     Observation *  The observation descriptor.
 *  mb         MapBeam *  Pointer to map/beam container.
 *  mw          Mapwin *  Clean window list.
 *  mpb     MaplotBeam *  The maplot beam position descriptor.
 *  vect    MaplotVect *  The object which describes how to draw polarization
 *                        vectors.
 *  domap          int    If true the map is to be plotted, else the beam.
 *  ctab        Ctable *  The color-table and transfer function to be used.
 *  docont         int    If true a contour plot is to be drawn.
 *  dovect         int    If true, and vect!=NULL, and domap is true,
 *                        polarization vectors will be drawn. The inner quarter
 *                        of the polarization intensity map is assumed to be in
 *                        the first quarter of mb->map[], and the polarization
 *                        angle map should be in the final quarter of mb->map[].
 *  levs         float *  Levels array to use for map - NULL selects the
 *                        default.
 *  nlev           int    The number of levels in levs. NB. levs, nlev and cmul
 *                        are only used when plotting the map.
 *  cmul         float    The multiplier to be applied to levs. If 0.0f,
 *                        percentage levels are selected.
 *  box          float[4] Array of the required X and Y limits of the plotted
 *                        region in radians, arranged as xa,xb,ya,yb.
 *                        If xa==xb or ya==yb then the maximum display area
 *                        will be plotted.
 *  markers MarkerList *  A list of markers to be displayed, or NULL if
 *                        not wanted.
 * Output:
 *  return         int    0 - OK.
 *                        1 - Error.
 */
int maplot(Observation *ob, MapBeam *mb, Mapwin *mw, MaplotBeam *mpb,
	   MaplotVect *vect, int domap, Ctable *ctab, int docont, int dovect,
	   int domod, float *levs, int nlev, float cmul, float *box,
	   MarkerList *markers)
{
  Maplot *mp;      /* Plot descriptor */
  int waserr = 0;  /* Error status */
/*
 * Allocate and initialize the plot descriptor.
 */
  mp = new_Maplot(ob, mb, mw, mpb, vect, domap, ctab, docont, dovect, domod,
		  levs, nlev, cmul, box, markers);
  if(mp==NULL)
    return 1;
/*
 * Interative session?
 */
  waserr = !mp->hard && mp->cursor && interact(mp);
/*
 * Transfer any newly created model components in mp->newmod to the
 * tentative model of the observation.
 */
  if(keep_cmp(mp))
    waserr = 1;
/*
 * Clean up.
 */
  mp = del_Maplot(mp);
  return waserr;
}

/*.......................................................................
 * Bring up the graphics cursor and wait for the user to press a key
 * to select a position on the plot. Return both the position selected
 * and the key pressed by the user.
 *
 * Input:
 *  mp     Maplot *  The plot descriptor.
 *  first     int    On the first call to this function for a given plot,
 *                   send this as true. The cursor will be positioned
 *                   in the center of the plot. Otherwise the cursor
 *                   will be positioned at its last selected position.
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
 * Input/Output:
 *  Keypos     kp *  The caller must send a pointer to a Keypos structure.
 *                   On output, the structure will contain the
 *                   position and key selected by the user.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int get_curs(Maplot *mp, int first, Bandmode mode,
		    float xref, float yref, int ci, Keypos *kp)
{
  static float xpos=0.0f; /* The X-world-coordinate of the cursor */
  static float ypos=0.0f; /* The Y-world-coordinate of the cursor */
  int waslow;             /* True if the key was lower case */
  char key;               /* The key that caused cpgcurs() to return */
/*
 * Check arguments.
 */
  if(mp==NULL || kp==NULL) {
    lprintf(stderr, "get_curs: NULL %s descriptor intercepted\n",
	    mp==NULL ? "Maplot":"Keypos");
    return 1;
  };
/*
 * Position the cursor in the center of the plot?
 */
  if(first) {
    xpos = (mp->wxa + mp->wxb) / 2.0f;
    ypos = (mp->wya + mp->wyb) / 2.0f;
  };
/*
 * Substitute the cross-hair cursor for the normal cursor is requested.
 */
  if(mode == B_NORM && mp->docross)
    mode = B_CROSS;
/*
 * Read a cursor position and the key that was pushed to return the
 * position.
 */
  cpgsci(ci);
  if(!cpgband((int) mode, 0, xref, yref, &xpos, &ypos, &key))
    return 1;
/*
 * Convert the selected key to upper case to simplify pending comparisons.
 */
  waslow = islower((int)key);
  if(waslow)
    key = toupper((int) key);
/*
 * Enforce plot bounds on the returned cursor position.
 */
  if(xpos < mp->wxa)
    xpos = mp->wxa;
  else if(xpos > mp->wxb)
    xpos = mp->wxb;
  if(ypos < mp->wya)
    ypos = mp->wya;
  else if(ypos > mp->wyb)
    ypos = mp->wyb;
/*
 * Record the results in the Keypos structure.
 */
  kp->x = xpos;
  kp->y = ypos;
  kp->key = key;
  kp->waslow = waslow;
  return 0;
}

/*.......................................................................
 * Delete and erase the CLEAN window, one of whose corners is the closest
 * window corner to the cursor.
 *
 * Input:
 *  mp    Maplot *  The plot descriptor.
 *  kp    Keypos *  The cursor selection.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static int zap_win(Maplot *mp, Keypos *kp)
{
  Subwin *win;
/*
 * Find the closest window.
 */
  win = findwin(mp, kp->x, kp->y);
/*
 * Erase the window from the display and delete it.
 */
  if(win) {
    draw_win(mp, win, 1);
    del_win(rem_win(mp->mw, NULL, win));
/*
 * Redraw all windows.
 */
    if(plwins(mp))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Allow the user to select a new CLEAN window using the cursor.
 *
 * Input:
 *  mp    Maplot *  The plot descriptor.
 *  kp    Keypos *  The initial cursor selection.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static int set_win(Maplot *mp, Keypos *kp)
{
  const int ptmark=1; /* The PGPLOT marker symbol to show the start vertex */
  Keypos kpb;         /* The second vertex selection */
  Subwin *win;        /* The new CLEAN window */
/*
 * Mark the start of the window selection with a dot.
 */
  cpgsci(wincol);
  cpgpt(1, &kp->x, &kp->y, ptmark);
/*
 * Get the second cursor selection.
 */
  for(;;) {
    kpb = *kp;
    if(get_curs(mp, 0, B_RECT, kp->x, kp->y, wincol, &kpb))
      return 1;
/*
 * Act on the key selection.
 */
    switch (kpb.key) {
    case KEY_CORN:
      if(kp->x!=kpb.x && kp->y!=kpb.y) {   /* Disallow 0-width windows */
/*
 * Add the new window to the list.
 */
	win = add_win(mp->mw, kp->x, kpb.x, kp->y, kpb.y);
	if(win == NULL)
	  return 1;
/*
 * Erase the start-vertex dot and draw the new window.
 */
	cpgsci(0);
	cpgpt(1, &kp->x, &kp->y, ptmark);
	cpgsci(1);
	draw_win(mp, win, 0);
	return 0;
      };
      break;
    case KEY_DEL:          /* Abort this window */
      cpgsci(0);
      cpgpt(1, &kp->x, &kp->y, ptmark);
      cpgsci(1);
      return 0;
      break;
    default:
      printf("You have selected one window corner - Use one of the following keys\n");
      printf(" %c - Select the opposite corner of the window you have started\n", KEY_CORN);
      printf(" %c - Discard the incomplete window\n", KEY_DEL);
    };
  };
}

/*.......................................................................
 * Allow the user to select a sub-image to be displayed, using the cursor.
 *
 * Input:
 *  mp    Maplot *   The plot descriptor.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
static int set_zoom(Maplot *mp)
{
  Keypos kp[2]; /* Vertex positions */
  int npts=0;   /* The number of vertices selected with the cursor */
/*
 * Tell user what to do.
 */
  printf("Select sub-image to be displayed - press %c for help\n", KEY_HELP);
/*
 * Get the two vertices of the sub-image.
 */
  while(npts<2) {
    if(get_curs(mp, 0, npts==0 ? B_NORM:B_RECT, kp[0].x, kp[0].y, zoomcol,
		&kp[npts]))
      return 1;
/*
 * Act on the key selection.
 */
    switch (kp[npts].key) {
    case KEY_ZOOM:   /* Select whole plot */
      kp[0].x = kp[0].y = kp[1].x = kp[1].y = 0.0f;
      npts = 2;
      break;
    case KEY_CORN:   /* Vertex selection */
      if(npts==0 || (kp[1].x != kp[0].x && kp[1].y != kp[0].y))
	npts++;
      break;
    case KEY_DEL:    /* Abort zoom */
      return 0;
      break;
    default:
      printf("You are currently in sub-image selection mode - please use keys:\n");
      printf(" %c - Select the %s of the required sub-image with this key\n",
	     KEY_CORN, npts==0 ? "two opposing corners":"opposite corner");
      printf(" %c - Select the whole %s\n", KEY_ZOOM, mp->domap?"map":"beam");
      printf(" %c - Abort selection\n", KEY_DEL);
      break;
    };
  };
/*
 * Initialze for the new area.
 */
  if(setarea(mp, kp[0].x, kp[1].x, kp[0].y, kp[1].y))
    return 1;
/*
 * Plot the requested area.
 */
  return replot(mp);
}

/*.......................................................................
 * Report on the value of the nearest pixel to the cursor.
 *
 * Input:
 *  mp     Maplot *  The plot descriptor.
 *  kp     Keypos *  The cursor selection.
 * Output:
 *  return Pixval *  A pointer to an internal static descriptor, containing
 *                   the details of the nearest pixel, or NULL on error.
 */
static Pixval *pix_val(Maplot *mp, Keypos *kp)
{
  static Pixval pv;/* The return descriptor */
  MapBeam *mb;     /* The map/beam descriptor. */
  int xcent;       /* The central X-axis pixel in the map/beam arrays */
  int ycent;       /* The central Y-axis pixel in the map/beam arrays */
  int ixmin,ixmax; /* X-axis pixel limits imposed by map/beam array */
  int iymin,iymax; /* Y-axis pixel limits imposed by map/beam array */
  int ix,iy;       /* Indexes of the pixel under the cursor */
/*
 * Check arguments.
 */
  if(mp==NULL || kp==NULL) {
    lprintf(stderr, "pix_val: NULL %s descriptor intercepted\n",
	    mp==NULL ? "Maplot":"Keypos");
    return NULL;
  };
/*
 * Get a pointer to the map/beam descriptor.
 */
  mb = mp->mb;
/*
 * Determine the central pixel on each axis.
 */
  xcent = mb->nx/2;
  ycent = mb->ny/2;
/*
 * Determine the pixel limits of the displayable map area.
 */
  ixmin = mb->nx/4;
  iymin = mb->ny/4;
  ixmax = 3*ixmin - 1;
  iymax = 3*iymin - 1;
/*
 * Find the pixel under to the cursor.
 */
  ix = xcent + (int) floor(kp->x/mb->xinc+0.5);
  iy = ycent + (int) floor(kp->y/mb->yinc+0.5);
/*
 * Enforce limits at edge of plot area.
 */
  if(ix<ixmin || ix>ixmax || iy<iymin || iy>iymax) {
    lprintf(stderr, "pix_val: Cursor out of plot bounds\n");
    return NULL;
  };
/*
 * Determine the world coords at the pixel center and the value of the pixel.
 */
  pv.x = radtoxy((ix-xcent) * mb->xinc);
  pv.y = radtoxy((iy-ycent) * mb->yinc);
  pv.value = mp->image[ix+iy*mb->nx];
  if(mp->dovect) {
    pv.poli = mp->image[(iy-iymin) * mb->nx/2 + (ix-ixmin)];
    pv.pola = mp->image[3*mb->ny/4 * mb->nx + (iy-iymin)*mb->nx/2+(ix-ixmin)];
  } else {
    pv.poli = pv.pola = 0.0f;
  };
/*
 * Return the pixel descriptor.
 */
  return &pv;
}

/*.......................................................................
 * Having alreay plotted the initial map/beam plot, start a session of
 * cursor directed interation.
 *
 * Input:
 *  mp    Maplot *   The plot descriptor.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
static int interact(Maplot *mp)
{
  Observation *ob;  /* The parent observation of the map */
  Ctable *ctab = mp->cmpar.ctab;  /* Color table */
  Keypos kp;        /* Cursor position and key-selection */
  Subwin *win;      /* Window closest to cursor */
  int first=1;      /* True after a new plot has been plotted */
  char buf[81];     /* Buffer to compose a label in */
/*
 * Tell user how to obtain key listing.
 */
  lprintf(stdout, "\nMove the cursor into the plot window and press \'%c\' for help\n",
	 KEY_HELP);
/*
 * Get a local pointer to the parent observation.
 */
  ob = mp->ob;
/*
 * Start the interactive editing loop.
 */
  do {
/*
 * Have the user select a position with the cursor. Get the cursor position
 * and the key that was pressed.
 */
    if(get_curs(mp, first, B_NORM, 0.0f, 0.0f, 1, &kp))
      return 1;
    first = 0;
/*
 * Act on the returned key selection.
 */
    switch(kp.key) {
    case KEY_CORN:           /* Set a new CLEAN window */
      if(set_win(mp, &kp))
	return 1;
      break;
    case KEY_DEL:            /* Delete the window with the closest side */
      if(zap_win(mp, &kp))
	return 1;
      break;
    case KEY_STAT:       /* Show stats on nearest window */
      win = findwin(mp, kp.x, kp.y);
      if(win)
	winstats(mp->mb, mp->domap, win, 0);
      break;
    case KEY_TRAN:
      if(change_transfer(mp, TR_LINEAR, 1))
	return 1;
      break;
    case KEY_DISP:       /* Re-display the plot */
      if(replot(mp))
	return 1;
      break;
    case KEY_ZOOM:       /* Allow selection of sub-image to be displayed */
      if(set_zoom(mp))
	return 1;
      first = 1;
      break;
    case KEY_MOD:        /* Toggle model component plotting */
/*
 * Clear currently plotted models?
 */
      {
	int doclr = (mp->domod || mp->dovar);
	if(kp.waslow) {
	  mp->domod = !mp->domod;
	  mp->dovar = 0;
	} else {
	  mp->dovar = !mp->dovar;
	  mp->domod = 0;
	};
	if((doclr && replot(mp)) || (!doclr && plmodel(mp)))
	  return 1;
      };
      break;
    case KEY_CMP:        /* Add a new cursor selected model component */
      if(make_cmp(mp))
	return 1;
      break;
    case KEY_REM:        /* Delete the nearest component to the cursor */
      if(zapcmp(mp, &kp))
	return 1;
      break;
    case KEY_VAL:        /* Report on the pixel under the cursor */
      {
	Pixval *pv = pix_val(mp, &kp);
	if(pv) {
	  double ra,dec;
	  lprintf(stdout, "Pixel value at x=%.3g y=%.3g (%s) is %.3g%s\n",
		  pv->x, pv->y, mapunits(U_PLAB), pv->value,
		  mp->domap ? " Jy/beam" : "");
/*
 * Display the polarization properties of a polarized pixel?
 */
	  if(mp->dovect) {
	    lprintf(stdout,
		    " Polarized flux=%.3g Jy/Beam, angle=%.4g degrees, P/%s=",
		    pv->poli, pv->pola * rtod,
		    Stokes_name(ob->stream.pol.type));
	    if(pv->value==0.0)
	      lprintf(stdout, "Infinity\n");
	    else
	      lprintf(stdout, "%.2g\n", pv->poli / pv->value);
	  };
/*
 * Display the Right Ascension and Declination of the selected point.
 */
	  ra = lmtora(ob->source.ra, ob->source.dec,
		      -ob->geom.east + xytorad(pv->x),
		      -ob->geom.north + xytorad(pv->y), ob->proj);
	  dec = lmtodec(ob->source.ra, ob->source.dec,
			-ob->geom.east + xytorad(pv->x),
			-ob->geom.north + xytorad(pv->y), ob->proj);
	  lprintf(stdout, " RA = %s,  ", sradhms(ra, 3, 0, buf));
	  lprintf(stdout, "Dec = %s (%.1f)\n", sraddms(dec, 3, 0, buf),
		  ob->source.epoch);
	};
      };
      break;
    case KEY_KEEP:       /* Save limits of displayed area for next time */
      mp->box[0] = mp->wxa;
      mp->box[1] = mp->wxb;
      mp->box[2] = mp->wya;
      mp->box[3] = mp->wyb;
      lprintf(stdout, "The displayed area limits have been saved for the next use of mapplot.\n");
      break;
    case KEY_COL:                  /* Change from grey-scale to psuedo-color */
      if(change_cmap(mp, CM_COLOR, !kp.waslow))
	return 1;
      break;
    case KEY_GRY:                  /* Change from psuedo-color to grey-scale */
      if(change_cmap(mp, CM_GREY, 0))
	return 1;
      break;
    case KEY_FIDL:                 /* Adjust colormap contrast and brightness */
      if(!mp->mono && ctab->cmap->class != CM_NONE) {
	if(kp.waslow) {
	  ctab->contra = 5.0 * kp.y / (kp.y < 0.0f ? mp->wya : -mp->wyb);
	  ctab->bright = (kp.x - mp->wxb)/(mp->wxa - mp->wxb);
	} else {
	  ctab->contra = 1.0f;
	  ctab->bright = 0.5f;
	  printf("Contrast and brightness reset.\n");
	};
	recolor(ctab->cmap, ctab->contra, ctab->bright);
      };
      break;
    case KEY_UNMARK:  /* Delete the nearest marker */
      if(zapmark(mp, &kp))
	return 1;
      break;
    case KEY_CROSS:   /* Toggle cross-hair cursor mode */
      mp->docross = !mp->docross;
      break;
    case KEY_HELP:
      printf("The following keys may be selected when the cursor is in the plot\n");
      printf(" %c - Quit this session\n", KEY_QUIT);
      printf(" %c - Select the two opposite corners of a new clean window.\n",
	     KEY_CORN);
      printf(" %c - Delete the window with a corner closest to the cursor.\n",
	     KEY_DEL);
      printf(" %c - Describe the area of the window with a corner closest to the cursor.\n", KEY_STAT);
      printf(" %c - Report the value of the pixel under the cursor.\n",KEY_VAL);
      printf(" %c - Fiddle the colormap contrast and brightness.\n",
	     tolower(KEY_FIDL));
      printf(" %c - Reset the colormap contrast and brightness to 1, 0.5.\n",
	     KEY_FIDL);
      printf(" %c - Re-display the plot.\n", KEY_DISP);
      printf(" %c - Install the default gray-scale color map.\n", KEY_GRY);
      printf(" %c - Install the default pseudo-color color map.\n",
	     tolower(KEY_COL));
      printf(" %c - Install a color map named at the keyboard.\n", KEY_COL);
      printf(" %c - Re-display with a different transfer function.\n",
	     KEY_TRAN);
      printf(" %c - Select a sub-image to be displayed.\n", KEY_ZOOM);
      printf(" %c - Retain the current sub-image limits for subsequent mapplot's\n",
	     KEY_KEEP);
      printf(" %c - Toggle display of the model.\n", tolower(KEY_MOD));
      printf(" %c - Toggle display of just the variable part of the model.\n",
	     KEY_MOD);
      printf(" %c - Initiate the description of a new model component.\n",
	     KEY_CMP);
      printf(" %c - Remove the model component closest to the cursor.\n",
	     KEY_REM);
      printf(" %c - Remove the marker closest to the cursor.\n", KEY_UNMARK);
      printf(" %c - Toggle whether to use a cross-hair cursor if available.\n",
	     KEY_CROSS);
      printf(" %c - List key bindings.\n", KEY_HELP);
      break;
    default:
      break;
    };
  } while(kp.key != KEY_QUIT);
  return 0;
}

/*.......................................................................
 * Replot the map/beam plot using the latest attributes in the Maplot
 * descriptor.
 *
 * Input:
 *  mp   Maplot *  The plot descriptor.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
static int replot(Maplot *mp)
{
  int ierr = 0;   /* Error return code */
/*
 * Sanity check.
 */
  if(bad_Maplot("replot", mp))
    return 1;
/*
 * Buffer output.
 */
  cpgbbuf();
/*
 * Start a new page.
 */
  cpgpage();
  mp->page++;
/*
 * Set up the color map for the first page.
 */
  ierr = mp->page==1 && change_cmap(mp, mp->cmpar.ctab->cmap->class, 0);
/*
 * Set up the viewport.
 */
  ierr = ierr || setport(mp);
/*
 * Display the false-color plot.
 */
  ierr = ierr || plimage(mp);
/*
 * Display the contour plot.
 */
  ierr = ierr || (mp->docont && plcont(mp));
/*
 * Display the polarization vectors.
 */
  ierr = ierr || (mp->dovect && plvect(mp));
/*
 * Display the CLEAN windows.
 */
  ierr = ierr || (mp->dowin && plwins(mp));
/*
 * Transfer the new model components in mp->newmod to the tentative
 * model of the observation.
 */
  ierr = ierr || keep_cmp(mp);
/*
 * Display model components.
 */
  ierr = ierr || plmodel(mp);
/*
 * Display markers.
 */
  ierr = ierr || plmarkers(mp);
/*
 * Label the plot.
 */
  ierr = ierr || pllabel(mp);
/*
 * Reveal the results.
 */
  cpgebuf();
  return ierr;
}

#define LWID 30  /* The max no. of char's used to represent a single level */
/*.......................................................................
 * Convert an array of contour levels into a textual listing, one line
 * per call. Only te levels visible in the current sub-image will be
 * anotated.
 *
 * Input:
 *  mp    Maplot *   The fully initialized plot descriptor.
 * Input/Output:
 *  ndone    int *   *ndone should hold the number of levels listed in
 *                   previous calls and on return will contain the
 *                   incremented count. On the first call this must be
 *                   pre-set to 0.
 * Output:
 *  return   int     Pointer to the returned label line, or NULL
 *                   on error.
 */
static char *lev_text(Maplot *mp, int *ndone)
{
  static char levtxt[81]; /* Buffer for the returned label */
  Contour *cpar;  /* Pointer to contour descriptor */
  int nchar=0;    /* The number of characters used in the output string */
  int trylen;     /* The no. of characters required to print the next level */
  char trytxt[LWID+1]; /* A string in which to find 'trylen'. */
  int ilev;       /* The next element of 'levs' to be processed */
  int nused=0;    /* Number of levels written */
  float scale;    /* Scale factor to convert to percentage levels. */
  float newlev;   /* The latest contour level */
  int maxwid;     /* Max numbers of chars to place in levtxt[] */
/*
 * Check arguments.
 */
  if(bad_Maplot("lev_text", mp))
    return NULL;
/*
 * Get a pointer to the contour plot descriptor.
 */
  cpar = &mp->cpar;
/*
 * Level array exhausted?
 */
  if(*ndone >= cpar->nlev)  /* All completed on previous call? */
    return NULL;
/*
 * Illegal peak value?
 */
  if(cpar->peak == 0.0)
    return NULL;
/*
 * What is the max number of chars that can go in levtxt[]?
 */
  maxwid = PORTWID/echsize;
  if(maxwid>=sizeof(levtxt)-1)
    maxwid = sizeof(levtxt)-1;
/*
 * Percentage levels?
 */
  if(cpar->plevs) {
/*
 * Write the line prefix.
 */
    sprintf(levtxt,"Contours %%:");
    nchar = strlen(levtxt);
/*
 * Determine the conversion factor to get percentage levels.
 */
    scale = 100.0/cpar->peak;
/*
 * Write as many levels into the return string as will fit.
 */
    for(ilev = *ndone; ilev<cpar->nlev; ilev++) {
      newlev = cpar->cmul * cpar->levs[ilev];
/*
 * Only anotate levels that are visible in the current sub-image.
 */
      if(newlev > cpar->cmin && newlev < cpar->cmax) {
/*
 * Write the next level into a temporary string to see if it will be
 * too long to append to the existing 'nchar' number of characters
 * in levtxt.
 */
	sprintf(trytxt, " %.3g", newlev * scale);
	trylen = strlen(trytxt);
/*
 * Stop if there is insufficient room.
 */
	if(nchar + trylen > maxwid)
	  break;
/*
 * Append the text for the new level.
 */
	strcpy(&levtxt[nchar], trytxt);
	nchar += trylen;
	nused++;
      };
    };
  }
/*
 * Absolute levels?
 */
  else {
/*
 * The first line starts with the multiplier.
 */
    if(*ndone==0)
      sprintf(levtxt, "Contours: %.3g %s x (",  cpar->cmul,
	      mp->domap ? "Jy/beam" : "/beam");
    else
      sprintf(levtxt, "Contours: ");
    nchar = strlen(levtxt);
/*
 * Write as many levels into the return string as will fit.
 */
    for(ilev = *ndone; ilev<cpar->nlev; ilev++) {
      newlev = cpar->cmul * cpar->levs[ilev];
/*
 * Only anotate levels that are visible in the current sub-image.
 */
      if(newlev > cpar->cmin && newlev < cpar->cmax) {
/*
 * Write the next level into a temporary string to see if it will be
 * too long to append to the existing 'nchar' number of characters
 * in levtxt.
 */
	sprintf(trytxt, "%.3g%c", cpar->levs[ilev],
		(ilev+1<cpar->nlev) ? ' ':')' );
	trylen = strlen(trytxt);
/*
 * Stop if there is insufficient room.
 */
	if(nchar + trylen > maxwid)
	  break;
/*
 * Append the text for the new level.
 */
	strcpy(&levtxt[nchar], trytxt);
	nchar += trylen;
	nused++;
      } else if(ilev+1==cpar->nlev) {  /* Ensure parenthesis closure */
	strcpy(&levtxt[nchar++], ")");
      };
    };
  };
/*
 * Unable to write a single level?
 */
  if(nused==0)
    return NULL;
  *ndone = ilev;
  return &levtxt[0];
}

/*.......................................................................
 * Determine how much room is required for labels and set the viewport
 * and world-coordinates for the image plot within this region.
 *
 * Input:
 *  mp   Maplot *  The fully initialized plot descriptor.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
static int setport(Maplot *mp)
{
  Contour *cpar;   /* Pointer to contour plot descriptor */
  Cmpar *cmpar;    /* Pointer to color-map plot descriptor */
  int nlines=0;    /* The number of anotation lines below primary X label */
  float botgap;    /* Gap to leave below viewport in character heights */
  float topgap;    /* Gap to leave above viewport in character heights */
  float lftgap;    /* Gap to leave left of viewport in character heights */
  float rgtgap;    /* Gap to leave right of viewport in character heights */
  float vxa,vxb,vya,vyb; /* Physical coordinates of viewport */
/*
 * Sanity check.
 */
  if(bad_Maplot("setport", mp))
    return 1;
/*
 * Get pointers to the contour and color-map plot descriptors.
 */
  cpar = &mp->cpar;
  cmpar = &mp->cmpar;
/*
 * Center RA/DEC anotation for map.
 */
  if(mp->domap)
    nlines++;
/*
 * Peak level anotation.
 */
  nlines++;
/*
 * Count the number of lines devoted to contour level anotation.
 */
  if(mp->docont) {
    int ndone=0;
    while(lev_text(mp, &ndone)!=NULL)
      nlines++;
  };
/*
 * Arrange to leave room for a beam description line?
 */
  if(mp->mb->ncmp)
    nlines++;
/*
 * Determine what gaps are to be left around the viewport to accomodate
 * text labels.
 */
  topgap = topsep + 1.0f + sepinc + 1.0f;
  lftgap = primsep + 1.0f;
  rgtgap = 0.0f;
  botgap = primsep + clevsep + echsize*nlines*(sepinc+1.0f) +
    (cmpar->ctab->cmap->class!=CM_NONE ? (wdginc+wdgsiz) : 0.0f);
/*
 * Set up the viewport leaving room for the labels (wedge) etc..
 * Also leave a margin around the labels.
 */
  cpgsvp(margin, 1.0f-margin, margin, 1.0f-margin);
  plsmvp(lftgap, rgtgap, botgap, topgap, 1.0f);
/*
 * Set up the window bounds (with RA going +ve to -ve).
 */
  cpgwnad(mp->wxb, mp->wxa, mp->wya, mp->wyb);
/*
 * Determine the physical aspect ratio of the plot window in
 * millimeters and use it to calculate the xtomm and ytomm conversions
 * from world coordinates to millimiters.
 */
  cpgqvp(2, &vxa, &vxb, &vya, &vyb);
  mp->xtomm = fabs((vxb - vxa)/(mp->wxb - mp->wxa));
  mp->ytomm = fabs((vyb - vya)/(mp->wyb - mp->wya));
  return 0;
}

/*.......................................................................
 * Display false-color plot of underlying image.
 * NB. setport() must be called before this function.
 *
 * Input:
 *  mp  Maplot *  The fully initialized plot descriptor.
 * Output:
 *  return int    0 - OK.
 *                1 - Error.
 */
static int plimage(Maplot *mp)
{
  Cmpar *cmpar;        /* Pointer to the color-map descriptor */
  int xdim;            /* X-axis dimension of image */
  int ydim;            /* Y-axis dimension of image */
  int pxa,pxb,pya,pyb; /* FORTRAN indices delimiting displayed area of array */
/*
 * Sanity check.
 */
  if(bad_Maplot("plimage", mp))
    return 1;
/*
 * Get a pointer to the color-map descriptor.
 */
  cmpar = &mp->cmpar;
/*
 * False-Color plot wanted?
 */
  if(cmpar->ctab->cmap->class != CM_NONE) {
/*
 * Get the dimensions of the image array.
 */
    xdim = mp->mb->nx;
    ydim = mp->mb->ny;
/*
 * Get FORTRAN array indices of plotted area of image array.
 */
    pxa = mp->pxa + 1;
    pxb = mp->pxb + 1;
    pya = mp->pya + 1;
    pyb = mp->pyb + 1;
/*
 * Plot the false-color image.
 */
    cpgsitf(cmpar->ctab->tran);
    if(!mp->mono)
      cpgimag(mp->image, xdim, ydim, pxa, pxb, pya, pyb,
	      cmpar->vmin, cmpar->vmax, mp->tr);
    else
      cpggray(mp->image, xdim, ydim, pxa, pxb, pya, pyb,
	      cmpar->vmin, cmpar->vmax, mp->tr);
  };
  return 0;
}

/*.......................................................................
 * Display contour plot of the image.
 * NB. setport() must be called before this function.
 *
 * Input:
 *  mp  Maplot *  The fully initialized plot descriptor.
 * Output:
 *  return int    0 - OK.
 *                1 - Error.
 */
static int plcont(Maplot *mp)
{
  const int poscol=1;  /* Color of positive contours */
  const int negcol=2;  /* Color of negative contours */
  Contour *cpar;       /* Pointer to the contour plot descriptor */
  float newlev;        /* The next contour level to be plotted */
  int xdim;            /* X-axis dimension of image */
  int ydim;            /* Y-axis dimension of image */
  int pxa,pxb,pya,pyb; /* FORTRAN indices delimiting displayed area of array */
  int i;
/*
 * Sanity check.
 */
  if(bad_Maplot("plcont", mp))
    return 1;
/*
 * Get a pointer to the contour plot descriptor.
 */
  cpar = &mp->cpar;
/*
 * Get the dimensions of the image array.
 */
  xdim = mp->mb->nx;
  ydim = mp->mb->ny;
/*
 * Get FORTRAN array indices of plotted area of image array.
 */
  pxa = mp->pxa + 1;
  pxb = mp->pxb + 1;
  pya = mp->pya + 1;
  pyb = mp->pyb + 1;
/*
 * Plot the contours - one at a time.
 */
  cpgbbuf();
  for(i=0; i<cpar->nlev; i++) {
    newlev = cpar->levs[i] * cpar->cmul;
/*
 * Don't attempt to plot contours unless they are going to be visible
 * in the current sub-image.
 */
    if(newlev > cpar->cmin && newlev < cpar->cmax) {
      cpgsci((newlev>=0.0f) ? poscol:negcol);
      cpgcont(mp->image, xdim, ydim, pxa, pxb, pya, pyb, &newlev, 1, mp->tr);
    };
  };
  cpgebuf();
  return 0;
}

/*.......................................................................
 * Write labels around the current viewport.
 *
 * Input:
 *  mp   Maplot *  The plot descriptor.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
static int pllabel(Maplot *mp)
{
  Observation *ob;      /* The descriptor of the observation */
  MapBeam *mb;          /* Map/beam descriptor */
  const int labcol=10;  /* Color of frame, labels etc.. */
  int nchar=0;          /* Number of characters used in label buffer */
  char buf1[81];        /* Buffer to compose a label in */
  char buf2[81];        /* Buffer to compose a label in */
  float xlabsep;        /* Incremental label separation from X axis */
/*
 * Sanity check.
 */
  if(bad_Maplot("pllabel", mp))
    return 1;
/*
 * Get a pointer to the map/beam descriptor.
 */
  mb = mp->mb;
/*
 * Get the observation descriptor.
 */
  ob = mp->ob;
/*
 * Plot axes.
 */
  cpgsch(1.0f);
  cpgsci(labcol);
/*
 * Temporarily change the units of the world-coordinates to those required
 * for labelling.
 */
  cpgswin(radtoxy(mp->wxb), radtoxy(mp->wxa),
	  radtoxy(mp->wya), radtoxy(mp->wyb));
  cpgbox("BCNST", 0.0f, 0, "BCNST", 0.0f, 0);
/*
 * Re-instate the normal world-coordinates.
 */
  cpgswin(mp->wxb, mp->wxa, mp->wya, mp->wyb);
/*
 * Compose a title.
 */
  sprintf(buf1, "%.16s\\fr at \\fn%.3f GHz %s", ob->source.name,
          getfreq(ob, -1)/1.0e9, sutdate(ob->date.year, ob->date.ut, buf2));
/*
 * Plot it.
 */
  cpgmtxt("T", topsep, 0.0f, 0.0f, buf1);
/*
 * Compose a miscellaneous title above the main title.
 */
  if(mp->domap)
    sprintf(buf1, "\\fr%s %s map. ", mb->ncmp ? "Clean" : "Residual",
	    Stokes_name(ob->stream.pol.type));
  else
    sprintf(buf1, "\\frDirty %s beam. ", Stokes_name(ob->stream.pol.type));
/*
 * Stations label or array label.
 */
  nchar = strlen(buf1);
  stnstr(ob, buf2, sizeof(buf2)-1);
  sprintf(&buf1[nchar], " Array: \\fn%.*s",
	  (int)(sizeof(buf1)-nchar-12), buf2);
/*
 * Display the new label.
 */
  cpgmtxt("T", topsep+1.0f+sepinc, 0.0f, 0.0f, buf1);
/*
 * Write the X and Y axis labels.
 */
  sprintf(buf1, "Relative Declination  (%s)", mapunits(U_PLAB));
  cpgmtxt("L", primsep, 0.5f, 0.5f, buf1);
  sprintf(buf1, "Right Ascension  (%s)", mapunits(U_PLAB));
  cpgmtxt("B", primsep, 0.5f, 0.5f, buf1);
/*
 * Proceed with the extra labels beneath the plot.
 */
  cpgsch(echsize);
  xlabsep = (primsep+clevsep)/echsize;
/*
 * Write the [shifted] center RA/DEC.
 */
  if(mp->domap) {
    double ra = lmtora(ob->source.ra, ob->source.dec, -ob->geom.east, -ob->geom.north, ob->proj);
    double dec = lmtodec(ob->source.ra, ob->source.dec, -ob->geom.east, -ob->geom.north, ob->proj);
    strcpy(buf1, "Map center:  ");
    sprintf(&buf1[strlen(buf1)], "RA: %s,  ", sradhms(ra, 3, 0, buf2));
    sprintf(&buf1[strlen(buf1)], "Dec: %s (%.1f)", sraddms(dec, 3, 0, buf2),
	    ob->source.epoch);
    xlabsep += 1.0f;
    cpgmtxt("B", xlabsep, 0.0f, 0.0f, buf1);
    xlabsep += sepinc;
  };
/*
 * Anotate peak level of map/beam.
 */
  if(mp->domap && mp->mb->ncmp)
    sprintf(buf1, "Map peak: %.3g Jy/beam", mp->cpar.peak);
  else
    sprintf(buf1, "Displayed range: %.3g \\(732) %.3g %s",
	    mp->cpar.cmin, mp->cpar.cmax, mp->domap ? "Jy/beam":" ");
  xlabsep += 1.0f;
  cpgmtxt("B", xlabsep, 0.0f, 0.0f, buf1);
  xlabsep += sepinc;
/*
 * Plot contour anotation if contours were drawn.
 */
  if(mp->docont) {
    int ndone=0;
    char *levtxt;
    while( (levtxt=lev_text(mp, &ndone)) != NULL) {
      xlabsep += 1.0f;
      cpgmtxt("B", xlabsep, 0.0f, 0.0f, levtxt);
      xlabsep += sepinc;
    };
  };
/*
 * Draw the restoring beam if this is a restored map and describe it
 * below the plot.
 */
  if(mb->ncmp) {
    MaplotBeam *mpb = mp->mpb;
    if(mpb) {
      float xmin = (mp->wxb - mp->wxa) * mpb->minsize;
      float xmax = (mp->wxb - mp->wxa) * mpb->maxsize;
      float ymin = (mp->wyb - mp->wya) * mpb->minsize;
      float ymax = (mp->wyb - mp->wya) * mpb->maxsize;
      plbeam(mb->bmin, mb->bmaj, mb->bpa, mpb->xc, mpb->yc,
	     xmin, xmax, ymin, ymax);
    };
    sprintf(buf1, "Beam FWHM: %.3g x %.3g (%s) at %.3g\\uo", radtoxy(mb->bmaj),
	    radtoxy(mb->bmin), mapunits(U_PLAB), mb->bpa*rtod);
    xlabsep += 1.0f;
    cpgmtxt("B", xlabsep, 0.0f, 0.0f, buf1);
    xlabsep += sepinc;
  };
/*
 * Draw a false-color wedge if required.
 */
  cpgsch(1.0f);
  xlabsep *= echsize;
  if(mp->cmpar.ctab->cmap->class != CM_NONE) {
    xlabsep += wdginc;
    plwedge("B", xlabsep, wdgsiz, mp->cmpar.vmax, mp->cmpar.vmin,
      mp->domap ? "Jy/beam" : "PSF response", mp->cmpar.ctab->tran, !mp->mono);
    xlabsep += wdgsiz;
  };
  return 0;
}

/*.......................................................................
 * Shrink the current viewport to allow for individually specified
 * margins around each edge, measured in character heights.
 *
 * Input:
 *  lftgap float  The no. of characters, to leave around left edge.
 *  rgtgap float  The no. of characters, to leave around right edge.
 *  topgap float  The no. of characters, to leave around bottom edge.
 *  botgap float  The no. of characters, to leave around top edge.
 *  chhgt  float  PGPLOT character height to use.
 */
static void plsmvp(float lftgap, float rgtgap, float botgap, float topgap,
	     float  chhgt)
{
  float vxa,vxb; /* Physical X limits of current viewport (mm) */
  float vya,vyb; /* Physical Y limits of current viewport (mm) */
  float xa,xb;   /* Physical X limits of the whole view-surface (mm) */
  float ya,yb;   /* Physical Y limits of the whole view-surface (mm) */
  float xwid_mm; /* Width of view-surface (mm). */
  float ywid_mm; /* Width of view-surface (mm). */
  float chsiz;   /* Physical character height (mm). */
/*
 * Find the physical dimensions of the current viewport (mm).
 */
  cpgqvp(2,&vxa,&vxb,&vya,&vyb);
/*
 * Find the physical dimensions of the whole view-surface (mm).
 */
  cpgsvp(0.0f,1.0f,0.0f,1.0f);
  cpgqvp(2,&xa,&xb,&ya,&yb);
  xwid_mm = xb-xa;
  ywid_mm = yb-ya;
/*
 * Determine the physical character height as set by PGPLOT.
 */
  plqcd(2, 1, &chsiz);
/*
 * Adjust the enclosing viewport to exclude the margins.
 */
  vxa += lftgap * chhgt * chsiz;
  vxb -= rgtgap * chhgt * chsiz;
  vya += botgap * chhgt * chsiz;
  vyb -= topgap * chhgt * chsiz;
/*
 * Set the adjusted viewport.
 */
  cpgsvp(vxa/xwid_mm, vxb/xwid_mm, vya/ywid_mm, vyb/ywid_mm);
  return;
}

/*.......................................................................
 * Return the default PGPLOT character height in a variety of units.
 * ie. the character height corresponding to cpgsch(1.0).
 *
 * Input:
 *   units     int    0 - Normalized device coordinates.
 *                    1 - inches.
 *                    2 - millimeters.
 *                    Other values are treated as 0.
 *   ishoriz   int    Only relevant when units=0 .
 *                    1 - Horizontally written characters.
 *                    0 - Vertically written characters.
 * Output:
 *   chdim   float *  The character height in the requested units.
 */
static void plqcd(int units, int ishoriz, float *chsize)
{
  float xa,xb;   /* Physical X limits of the whole view-surface (mm) */
  float ya,yb;   /* Physical Y limits of the whole view-surface (mm) */
  float xwid_mm; /* Width of view-surface (mm). */
  float ywid_mm; /* Width of view-surface (mm). */
  float vxa,vxb,vya,vyb; /* Temporary record of viewport coordinates */
/*
 * Store the current viewport coordinates.
 */
  cpgqvp(0,&vxa,&vxb,&vya,&vyb);
/*
 * Find the physical dimensions of the whole view-surface (mm).
 */
  cpgsvp(0.0f,1.0f,0.0f,1.0f);
  cpgqvp(2,&xa,&xb,&ya,&yb);
  xwid_mm = xb-xa;
  ywid_mm = yb-ya;
/*
 * Determine the physical character height as set by PGPLOT.
 * This chhgt/40 x the size of the shortest edge of the viewsurface.
 */
  if(ywid_mm/xwid_mm < 1.0f)
    *chsize = ywid_mm/40.0f;
  else
    *chsize = xwid_mm/40.0f;
/*
 * Return this size in the appropriate units.
 */
  switch (units) {
  case 1:
    *chsize /= 25.4f; /* Convert millimeters to inches */
    break;
  case 2:            /* Already in millimeters */
    break;
  default:           /* Convert to NDC */
    if(ishoriz)
      *chsize /= ywid_mm;
    else
      *chsize /= xwid_mm;
    break;
  };
/*
 * Reset to the original viewport.
 */
  cpgsvp(vxa,vxb,vya,vyb);
  return;
}

/*.......................................................................
 * Plot an anotated false-color calibration wedge alongside, beneath or
 * above the current viewport.
 * NB. The viewport and world coords will be restored to their original
 * values on return.
 *
 * Input:
 *   side    char * A single char string specifying the side of the
 *                  viewport parallel to which the wedge will be drawn.
 *                   "B" - Bottom edge.
 *                   "T" - Top edge.
 *                   "L" - Left edge.
 *                   "R" - Right edge.
 *   disp   float   Displacement outwards from the specified side in
 *                  units of character heights.
 *   width  float   Width of wedge + anotation in char's.
 *   fg     float   The label to be given to the foreground color.
 *   bg     float   The label to be given to the background color.
 *   label   char * The units label (if NULL no label will be drawn).
 *   mode     int   0 - Use linear transfer function.
 *                  1 - Use logarithmic transfer function.
 *   doimag   int   If true and PGIMAG is available, use PGIMAG.
 *                  Otherwise use PGGREY.
 * Output:
 *   return   int   0 - OK.
 *                  1 - Illegal side requested.
 */
static int plwedge(char *side, float disp, float width, float fg, float bg,
		   char *label, int mode, int doimag)
{
  float vxa,vxb,vya,vyb; /* New viewport coordinates */
  float wxa,wxb,wya,wyb; /* Temporary storage of world coordinates */
  float xa,xb,ya,yb;     /* Temporary storage of viewport coordinates */
  float old_ch;          /* Temporary storage of entry character height */
  float wdginc;          /* False-color increment per wedge element */
  float tr[6];           /* Gray-scale element mapping array */
  float chsize;          /* Default PGPLOT character height (mm) */
  float newhgt;          /* Anotation character height */
  float labwid;          /* Width of anotation (char's) */
  int ishoriz;           /* True if wedge to be drawn horizontally */
  const float txtfrc=0.6f;/* Fraction of 'width' used for anotation */
  const float margin=0.1f;/* Margin around all edges of viewport */
  const float txtsep=2.2f;/* Separation between units text and axis */
                          /* (character heights). */
#define WDGPIX 100
  float wdgarr[WDGPIX];   /* Temporary array to hold wedge */
  int i;
/*
 * Check that a valid 'side' was cited.
 */
  switch (*side) {
  case 'B': case 'T': case 'L': case 'R':
    break;
  default:
    lprintf(stderr, "plwedge: Illegal side requested\n");
    return 1;
  };
/*
 * Store the current world and viewport coordinates and the character height.
 */
  cpgqwin(&wxa, &wxb, &wya, &wyb);
  cpgqvp(0,&xa,&xb,&ya,&yb);
  cpgqch(&old_ch);
/*
 * Determine the default character height in NDC coords.
 */
  ishoriz = *side=='T' || *side=='B';
  plqcd(0, ishoriz, &chsize);
/*
 * Convert 'width' and 'disp' into viewport units.
 */
  width *= chsize * old_ch;
  disp  *= chsize * old_ch;
/*
 * Use these to determine viewport coordinates for the wedge + annotation.
 */
  vxa = xa;
  vxb = xb;
  vya = ya;
  vyb = yb;
  switch (*side) {
  case 'B':
    vyb = ya - disp;
    vya = vyb - width;
    break;
  case 'T':
    vya = yb + disp;
    vyb = vya + width;
    break;
  case 'L':
    vxb = xa - disp;
    vxa = vxb - width;
    break;
  case 'R':
    vxa = xb + disp;
    vxb = vxa + width;
    break;
  };
/*
 * Set the viewport.
 */
  cpgsvp(vxa, vxb, vya, vyb);
/*
 * Determine and set the character height required to fit the wedge
 * anotation text within the area allowed for it.
 */
  newhgt = txtfrc*width / ((txtsep+1.0f)*chsize); 
  cpgsch(newhgt);
/*
 * Adjust the viewport to enclose just the wedge, leaving room for the
 * numeric and units labels.
 */
  labwid = txtsep+((label!=0)?1.0f:0.0f);
  plsmvp(margin+((*side=='L')?labwid : 0.0f),
	 margin+((*side=='R')?labwid : 0.0f),
	 margin+((*side=='B')?labwid : 0.0f),
	 margin+((*side=='T')?labwid : 0.0f), newhgt);
/*
 * Set up for drawing the false-color wedge.
 */
  tr[0]=0.0f;
  tr[1]=1.0f;
  tr[2]=0.0f;
  tr[3]=0.0f;
  tr[4]=0.0f;
  tr[5]=1.0f;
  wdginc = (fg-bg)/(WDGPIX-1);
  for(i=0; i<WDGPIX; i++)
    wdgarr[i] = bg + i * wdginc;
/*
 * Draw the wedge then change the world coordinates for labelling.
 */
  if(ishoriz) {
    cpgswin(1.0f, (float) WDGPIX, 0.9f, 1.1f);
    cpgsitf(mode);
    if(doimag)
      cpgimag(wdgarr,WDGPIX,1, 1,WDGPIX, 1,1, bg,fg,tr);
    else
      cpggray(wdgarr,WDGPIX,1, 1,WDGPIX, 1,1, bg,fg,tr);
    cpgswin(bg,fg,0.0f,1.0f);
  } else {
    cpgswin(0.9f, 1.1f, 1.0f, (float) WDGPIX);
    cpgsitf(mode);
    if(doimag)
      cpgimag(wdgarr,1,WDGPIX, 1,1, 1,WDGPIX, bg,fg,tr);
    else
      cpggray(wdgarr,1,WDGPIX, 1,1, 1,WDGPIX, bg,fg,tr);
    cpgswin(0.0f,1.0f,bg,fg);
  };
/*
 * Set the world coordinates for labelling, then draw a labelled box.
 */
  switch (*side) {
  case 'B':
    cpgbox("BCNST",0.0f,0,"BC",0.0f,0);
    break;
  case 'T':
    cpgbox("BCMST",0.0f,0,"BC",0.0f,0);
    break;
  case 'L':
    cpgbox("BC",0.0f,0,"BCNST",0.0f,0);
    break;
  case 'R':
    cpgbox("BC",0.0f,0,"BCMST",0.0f,0);
    break;
  };
/*
 * Write the units label.
 */
  if(label != 0)
    cpgmtxt(side,txtsep,1.0f,1.0f,label);
/*
 * Reset the original viewport and world coordinates.
 */
  cpgsvp(xa,xb,ya,yb);
  cpgswin(wxa,wxb,wya,wyb);
  cpgsch(old_ch);
  return 0;
}

/*.......................................................................
 * Determine which plotted window has the nearest vertex to a given
 * cursor selected X/Y coordinate.
 *
 * Input:
 *  mp     Maplot *  The plot descriptor.
 *  xpos    float    The X-coordindate of the point (radians).
 *  ypos    float    The Y-coordindate of the point (radians).
 * Output:
 *  return Subwin    The window with the nearest vertex to xpos,ypos.
 *                   (or NULL if mwin==NULL or mwin->nwin==0).
 */
static Subwin *findwin(Maplot *mp, float xpos, float ypos)
{
  Mapwin *mwin;        /* The list of CLEAN windows */
  Subwin *w;           /* Displayed limits of latest window */
  Subwin *win;         /* The latest window being considered */
  Subwin *minwin=NULL; /* The closest window so far */
  float xadif,xbdif;   /* Squared X-dist from a window vertex to the cursor */
  float yadif,ybdif;   /* Squared Y-dist from a window vertex to the cursor */
  float rnew;          /* Min absolute dist from xpos,ypos to window corners */
  float rmin;          /* Minimum of 'rnew' over all windows */
/*
 * Sanity check.
 */
  if(bad_Maplot("findwin", mp))
    return NULL;
/*
 * Get the window list.
 */
  mwin = mp->mw;
/*
 * No window to return?.
 */
  if(mwin==0 || mwin->nwin==0)
    return NULL;
/*
 * Iterate through the list of windows.
 */
  rmin = -1.0f;  /* Impossible absolute distance */
  for(win=mwin->head; win != NULL; win = win->next) {
    w = win_limits(mp, win);
    if(w) {
/*
 * Find the physical (mm) distance (squared) of the vertex of
 * the current window that is closest to the cursor position xpos,ypos.
 * floatmin() is defined in minmax.c and returns its minimum float argument.
 */
      xadif = (w->xmin - xpos) * mp->xtomm; xadif *= xadif;
      xbdif = (w->xmax - xpos) * mp->xtomm; xbdif *= xbdif;
      yadif = (w->ymin - ypos) * mp->ytomm; yadif *= yadif;
      ybdif = (w->ymax - ypos) * mp->ytomm; ybdif *= ybdif;
      rnew = floatmin(xadif,xbdif) + floatmin(yadif,ybdif);
/*
 * If the closest side is closer than any for previous windows then
 * record the previous window since this gives the id of the current
 * window and that of the previous window (used for unlinking).
 */
      if(rnew < rmin || rmin < 0.0f) {
	rmin = rnew;
	minwin = win;
      };
    };
  };
  return minwin;
}

/*.......................................................................
 * Return the displayed limits of a window.
 *
 * Input:
 *  mp     Maplot *  The plot descriptor.
 *  win    Subwin *  The window.
 * Output:
 *  return Subwin *  Pointer to an internal static repository for the
 *                   returned limits, or NULL if the window was wholy
 *                   outside the plot.
 */
static Subwin *win_limits(Maplot *mp, Subwin *win)
{
  static Subwin w;  /* The displayed limits of the window */
/*
 * No window or display?
 */
  if(bad_Maplot("win_limits", mp))
    return NULL;
  if(win==NULL) {
    lprintf(stderr, "win_limits: Syserror: intercepted NULL window.\n");
    return NULL;
  };
/*
 * Shrink the window edges to fit the display.
 */
  w.xmin = (win->xmin < mp->wxa) ? mp->wxa : win->xmin;
  w.xmax = (win->xmax > mp->wxb) ? mp->wxb : win->xmax;
  w.ymin = (win->ymin < mp->wya) ? mp->wya : win->ymin;
  w.ymax = (win->ymax > mp->wyb) ? mp->wyb : win->ymax;
/*
 * Return w unless the window is wholly outside the displayed area.
 */
  if(w.xmin > mp->wxb || w.xmax < mp->wxa ||
     w.ymin > mp->wyb || w.ymax < mp->wya)
    return NULL;
  else {
    w.next = NULL;
    return &w;
  };
}

/*.......................................................................
 * Draw a single CLEAN window. Only the part (if any) that appears on the
 * display will be drawn.
 *
 * Input:
 *  mp     Maplot *  The plot descriptor.
 *  win    Subwin *  The window to be drawn.
 *  erase     int    If true then erase the window.
 */
static int draw_win(Maplot *mp, Subwin *win, int erase)
{
  Subwin *wlims;        /* *wlims Will contain window limits */
/*
 * No window to display?
 */
  if(win==NULL || mp==NULL) {
    lprintf(stderr, "draw_win: Syserror: intercepted NULL %s.\n",
	    win==NULL ? "window" : "plot descriptor");
    return 1;
  };
/*
 * Determine the visible area of the window.
 */
  wlims = win_limits(mp, win);
/*
 * Display or erase the window if it is visible.
 */
  if(wlims) {
    cpgsfs(2);
    cpgsci(erase ? 0 : wincol);
    cpgrect(wlims->xmin, wlims->xmax, wlims->ymin, wlims->ymax);
    cpgsci(wincol);
  };
  return 0;
}

/*.......................................................................
 * Display all the visible CLEAN windows.
 *
 * Input:
 *  mp   Maplot *  The plot descriptor.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
static int plwins(Maplot *mp)
{
  Mapwin *mw = mp->mw;  /* The window list container */
  Subwin *win;          /* A window from the list */
  int ierr = 0;         /* Error return value */
/*
 * Display the visible windows from the list.
 */
  if(mw) {
    cpgbbuf();
    for(win=mw->head; win!=NULL && !ierr; win = win->next)
      ierr = draw_win(mp, win, 0);
    cpgebuf();
  }
  return ierr;
}

/*.......................................................................
 * Display model components symbolically.
 *
 * Input:
 *  mp    Maplot *   The plot descriptor.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
static int plmodel(Maplot *mp)
{
/*
 * Should we plot anything from the established and tentative models?
 */
  if(mp->domod || mp->dovar) {
    int dovar = mp->domod || mp->dovar;
    int dofix = mp->domod;
    int nhidden = 0;
/*
 * Plot the components of the established model.
 */
    nhidden += modplot(mp->ob->model, dofix, dovar, mp->wxa, mp->wxb,
		       mp->wya, mp->wyb);
/*
 * Plot the components of the tentative model.
 */
    nhidden += modplot(mp->ob->newmod, dofix, dovar, mp->wxa, mp->wxb,
		       mp->wya, mp->wyb);
/*
 * If some components were not plotted because they lay outside the plot area, 
 * report their number.
 */
    if(nhidden) {
      lprintf(stdout, "Note that %d components lie outside the plot.\n",
	      nhidden);
    };
  };
  return 0;
}

/*.......................................................................
 * Allow the user to make a new model component with the cursor.
 *
 * Input:
 *  mp    Maplot *   The plot descriptor.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
static int make_cmp(Maplot *mp)
{
  const int ptmark=1; /* The PGPLOT marker symbol to show tentative vertex */
  const int tmpcol=11;/* Color used to display work lines */
  Keypos newkp;    /* The latest cursor selection */
  int npts=0;      /* The number of vertices selected with the cursor */
  int cancelled=0; /* True when the selection is cancelled */
  int completed=0; /* True when the selection is completed */
  int oldcol;      /* Color to be reinstated before return */
/*
 * Define default component values, to be updated as more parameters
 * are described.
 */
  int type=M_DELT;    /* The component type */
  float flux;         /* The flux of the component */
  float x=0.0f;       /* The X position of the component (default = 0) */
  float y=0.0f;       /* The Y position of the component (default = 0) */
  float major=0.0f;   /* The FWHM major axis of the component */
  float minor=0.0f;   /* The FWHM minor axis of the component */
  float ratio=1.0f;   /* The axial ratio minor/major */
  float phi=0.0f;     /* The position angle of the major axis (radians) */
  int freepar=M_FLUX|M_CENT; /* Bitmap of parameters to be made variable */
  float xmajor=0.0;   /* X-axis position of the end of the semi-major axis */
  float ymajor=0.0;   /* Y-axis position of the end of the semi-major axis */
  float xminor=0.0;   /* X-axis position of the end of the semi-major axis */
  float yminor=0.0;   /* Y-axis position of the end of the semi-major axis */
/*
 * Tell user what to do.
 */
  printf("Describe a new component (press '%c' for help, '%c' to cancel).\n",
	 KEY_HELP, KEY_DEL);
/*
 * Keep a copy of the entry color.
 */
  cpgqci(&oldcol);
/*
 * Give the work lines a distinct color.
 */
  cpgsci(tmpcol);
/*
 * Get up to three vertices and a terminator.
 */
  while(!completed && !cancelled) {
    int showhelp = 0;
/*
 * Get cursor input from the user.
 */
    if(get_curs(mp, 0, npts==1||npts==2 ? B_LINE:B_NORM, x, y, tmpcol, &newkp))
      return 1;
/*
 * Act on the key selection.
 */
    switch (newkp.key) {
    case KEY_CORN:     /* Select new vertex */
/*
 * If the point is acceptable, plot a point to mark its location,
 * record the key characteristics and increment the recorded number of
 * vertices selected.
 */
      if(npts < 3) {
	switch(++npts) {
	case 1:
	  freepar |= M_CENT;  /* Define component centroid */
	  x = newkp.x;
	  y = newkp.y;
	  cpgpt(1, &x, &y, ptmark); /* Plot central position */
	  break;
	case 2:
	  freepar |= M_MAJOR; /* Define the major axis */
	  {
	    float xdif = newkp.x - x;
	    float ydif = newkp.y - y;
	    major = 2.0 * sqrt(xdif * xdif + ydif * ydif);
	    phi = (xdif==0.0 && ydif==0.0) ? 0.0 : atan2(xdif, ydif);
	  };
/*
 * Determine the coordinates of the end of the semi-major axis.
 */
	  xmajor = x + 0.5 * major * sin(phi);
	  ymajor = y + 0.5 * major * cos(phi);
/*
 * Draw semi-major axis.
 */
	  cpgmove(x,y);
	  cpgdraw(xmajor, ymajor);
	  break;
	case 3:
	  freepar |= (M_RATIO | M_PHI);
	  {
	    float xdif = newkp.x - x;
	    float ydif = newkp.y - y;
/*
 * Get the angle corresponding to the semi-minor axis on the side of the major
 * axis in which the new selection was made.
 */
	    float xyphi = (xdif==0.0 && ydif==0.0) ? 0.0 : atan2(xdif, ydif);
	    float posphi = fmod(phi+pi,pi);
	    float minor_pa = (xyphi > 0.0 && xyphi < posphi) ||
	                     (xyphi > posphi-pi && xyphi < 0.0) ?
		       (posphi - halfpi) : (posphi + halfpi);
/*
 * Determine the FWHM minor axis.
 */
	    minor = 2.0 * sqrt(xdif * xdif + ydif * ydif);
/*
 * Determine the coordinates of the end of the semi-minor axis.
 */
	    xminor = x + 0.5 * minor * sin(minor_pa);
	    yminor = y + 0.5 * minor * cos(minor_pa);
	  };
/*
 * Draw semi-minor axis.
 */
	  cpgmove(x,y);
	  cpgdraw(xminor, yminor);
	  break;
	};
      } else {
        showhelp = 1;  /* No more vertices to select - show help info */
      };
      break;
    case KEY_DEL:      /* Cancel the selection */
      cancelled = 1;
      break;
    case KEY_CMP:      /* Complete the selection */
      completed = 1;
      break;
    default:           /* Un-bound key or KEY_HELP - display help info */
      showhelp = 1;
      break;
    };
/*
 * If requested, display help information pertinent to the number of
 * vertices so far selected.
 */
    if(showhelp) {
      printf("You are currently creating a new model component - use keys:\n");
      printf(" %c - Abort the component selection.\n", KEY_DEL);
      switch(npts) {
      case 0:
	printf(" %c - Select the component center.\n", KEY_CORN);
	printf(" %c - Install a delta component at the map center.\n", KEY_CMP);
	break;
      case 1:
	printf(" %c - Terminate the major axis radius vector.\n", KEY_CORN);
	printf(" %c - Install a delta component.\n", KEY_CMP);
	break;
      case 2:
	printf(" %c - Terminate the minor axis radius length.\n", KEY_CORN);
	printf(" %c - Install a circular gaussian component.\n", KEY_CMP);
	break;
      default:
	printf(" %c - Install an elliptical gaussian component.\n", KEY_CMP);
	break;
      };
    };
  };
/*
 * Erase the center marker and axes that marked the chosen selections
 * as they were drawn. (Note the intentional case fallthroughs).
 */
  cpgsci(0);
  switch(npts) {
  case 3:
    cpgmove(x,y);
    cpgdraw(xminor, yminor);
  case 2:
    cpgmove(x,y);
    cpgdraw(xmajor, ymajor);
  default:
    cpgpt(1, &x, &y, ptmark);
  };
/*
 * Re-instate the entry color.
 */
  cpgsci(oldcol);
/*
 * Install the new component unless cancelled.
 */
  if(!cancelled) {
/*
 * Make sure that the major axis is the longest of the axes.
 */
    if(major < minor) {
      float new_major = minor;
      minor = major;
      major = new_major;
      phi = fmod(phi + halfpi, pi);
    };
/*
 * What type does the component have?
 */
    type = (npts < 2 || major==0.0f) ? M_DELT : M_GAUS;
/*
 * If the component is elliptical, determine the axial ratio from
 * the selected major and minor axis extents.
 */
    if(type==M_GAUS && npts > 2)
      ratio = minor / major;
/*
 * Set the flux equal to the flux of the pixel at the center of the
 * component. A better method would be to determine the flux of the
 * delta function or 2-D gaussian (convolved with the clean when known)
 * that would produce the observed flux/density at the component center.
 * This would be complicated, and since this is only an estimated of the
 * flux, probably unwarranted, especially since we also aren't taking into
 * account other components that have been added during this mapplot session.
 */
    {
      Pixval *pv;
      newkp.x = x;
      newkp.y = y;
      pv = pix_val(mp, &newkp);
      if(pv==NULL)
	return 1;
      flux = pv->value;
    };
/*
 * Install the new component in the tentative model.
 */
    {
      Modcmp *cmp = add_xycmp(mp->newmod, 1, freepar, flux, x, y, major,
			    ratio, phi, type, getfreq(mp->ob, -1), 0.0);
      if(cmp==NULL)
	return 1;
/*
 * Display the new component.
 */
      cmpplot(cmp, mp->wxa, mp->wxb, mp->wya, mp->wyb, 0);
    };
  };
  return 0;
}

/*.......................................................................
 * Transfer the most recently created components (now in mp->newmod) to
 * the tentative model in the observation. This should be done each time
 * before the plot is re-displayed, and before return from this session.
 *
 * Input:
 *  mp    Maplot *  The plot descriptor.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static int keep_cmp(Maplot *mp)
{
/*
 * Are there any new components to be transfered?
 */
  if(mp->newmod->ncmp > 0) {
/*
 * The addition of a new component invalidates the current map.
 */
    mp->mb->domap = 1;
/*
 * Transfer the components from mp->newmod to the tentative model of the
 * observation.
 */
    if(obaddmod(mp->ob, mp->newmod, 0, 0, 1))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Create a map/beam plot descriptor to retain state info for the duration
 * of a plotting session. Call del_Maplot() to delete it at the end of the
 * session.
 *
 * Input:
 *  ob     Observation *  The observation descriptor.
 *  mb         MapBeam *  Pointer to map/beam container.
 *  mw          Mapwin *  Clean window list.
 *  mpb     MaplotBeam *  The maplot clean beam location descriptor.
 *  vect    MaplotVect *  The object which describes how to draw polarization
 *                        vectors.
 *  domap          int    If true the map is to be plotted, else the beam.
 *  ctab        Ctable *  The color-table and transfer function to use.
 *  docont         int    If true a contour plot is to be drawn.
 *  dovect         int    If true, and vect!=NULL, and domap is true,
 *                        polarization vectors will be drawn. The inner quarter
 *                        of the polarization intensity map is assumed to be in
 *                        the first quarter of mb->map[], and the polarization
 *                        angle map should be in the final quarter of mb->map[].
 *  levs         float *  Levels array to use for map - NULL selects the
 *                        default.
 *  nlev           int    The number of levels in levs. NB. levs, nlev and cmul
 *                        are only used when plotting the map.
 *  cmul         float    The multiplier to be applied to levs. If 0.0f,
 *                        percentage levels are selected.
 *  box          float[4] Array of the required X and Y limits of the plotted
 *                        region in radians, arranged as xa,xb,ya,yb.
 *                        If xa==xb or ya==yb then the maximum display area
 *                        will be plotted.
 *  markers MarkerList *  A list of markers to be displayed, or NULL if
 *                        not wanted.
 * Output:
 *  return      Maplot *  The descriptor, or NULL on error.
 */
static Maplot *new_Maplot(Observation *ob, MapBeam *mb, Mapwin *mw,
			  MaplotBeam *mpb, MaplotVect *vect, int domap,
			  Ctable *ctab, int docont, int dovect, int domod,
			  float *levs, int nlev, float cmul, float *box,
			  MarkerList *markers)
{
  Maplot *mp;      /* Plot descriptor */
  char answer[5];  /* String for reply from PGPLOT inquiry function */
  int slen;        /* Length of answer to PGPLOT inquiry */
  static float defbox[4]={0.0f,0.0f,0.0f,0.0f}; /* Default box */
/*
 * Check the required descriptors.
 */
  if(!ob_ready(ob, OB_SELECT, "mapplot"))
    return NULL;
  if(mb==NULL) {
    lprintf(stderr, "mapplot: NULL MapBeam descriptor intercepted.\n");
    return NULL;
  };
/*
 * Is there a PGPLOT device currently open?
 */
  slen = sizeof(answer)-1;
  cpgqinf("OPEN", answer, &slen);
  if(strncmp(answer, "NO", 2)==0) {
    lprintf(stderr, "new_Maplot: No PGPLOT device active\n");
    return NULL;
  };
/*
 * Allocate the plot descriptor.
 */
  mp = (Maplot *) malloc(sizeof(Maplot));
  if(mp==NULL) {
    lprintf(stderr, "new_Maplot: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize all pointers that will point to dynamically allocated memory
 * to NULL so that the descriptor can hereafter be safely sent to del_Maplot().
 */
  mp->newmod = NULL;
/*
 * Initialize the rest of the descriptor.
 */
  mp->ob = ob;
  mp->mb = mb;
  mp->mw = mw;
  mp->mpb = mpb;
  mp->vect = vect;
  mp->box = box ? box : &defbox[0];
/*
 * Allocate an empty model to be used as a way station for newly created
 * model components.
 */
  mp->newmod = new_Model();
  if(mp->newmod==NULL)
    return del_Maplot(mp);
/*
 * Is the plot device a hard-copy device?
 */
  slen = sizeof(answer)-1;
  cpgqinf("HARDCOPY", answer, &slen);
  mp->hard = strncmp(answer, "YES", 3)==0;
/*
 * Has the device got enough colors to sensibly use PGIMAG?
 */
  {
    int minind, maxind;  /* Min/max useable color index */
    cpgqcir(&minind, &maxind);
    mp->mono = maxind-minind+1 < MINLEVS;
  };
/*
 * Is the plot device a hard-copy device?
 */
  slen = sizeof(answer)-1;
  cpgqinf("CURSOR", answer, &slen);
  mp->cursor = strncmp(answer, "YES", 3)==0;
  mp->docross = 0;
/*
 * On restored maps, only plot windows on interactive devices.
 */
  mp->dowin = mw!=NULL && (!mp->mb->ncmp || !mp->hard);
/*
 * Should we plot the model?
 */
  mp->domod = domod && ob->model!=NULL;
  mp->dovar = 0;
/*
 * Image plotting options.
 */
  mp->docont = docont;
  mp->dovect = dovect && vect != NULL && domap && vect->scale > 0.0;
  mp->domap = domap;
  mp->image = domap ? mb->map : mb->beam;
/*
 * If neither false-color nor contours were requested, then signal an error.
 */
  if(!mp->docont && ctab->cmap->class == CM_NONE) {
    lprintf(stderr, "mapplot: Neither false-color nor contours requested\n");
    return del_Maplot(mp);
  };
/*
 * Initialize the contour plot parameters for subsequent plots.
 */
  if(setcont(mp, levs, cmul, nlev))
    return del_Maplot(mp);
/*
 * Initialize the color-map plot parameters for subsequent plots.
 */
  if(setcmpar(mp, ctab))
    return del_Maplot(mp);
/*
 * Set up to plot the given area of the map/beam.
 */
  if(setarea(mp, mp->box[0], mp->box[1], mp->box[2], mp->box[3]))
    return del_Maplot(mp);
/*
 * Record the optional list of markers to be displayed in maps.
 */
  mp->markers = domap ? markers : NULL;
/*
 * No pages yet.
 */
  mp->page = 0;
/*
 * Display the initial plot.
 */
  if(replot(mp))
    return del_Maplot(mp);
  return mp;
}

/*.......................................................................
 * Delete a Maplot descriptor.
 *
 * Input:
 *  mp     Maplot *   The descriptor to be deleted.
 * Output:
 *  return Maplot *   NULL.
 */
static Maplot *del_Maplot(Maplot * mp)
{
  if(mp) {
    del_Model(mp->newmod);
    free(mp);
  };
  return NULL;
}

/*.......................................................................
 * Return details of the closest model component to the given cursor
 * position.
 *
 * Input:
 *  mp      Maplot *   The plot descriptor.
 *  kp      Keypos *   The descriptor of the cursor position selection.
 * Output:
 *  return  Cmpfnd *   Pointer to static internal descriptor detailing
 *                     the nearest component, or NULL if no component
 *                     was found.
 */
Cmpfnd *fnd_cmp(Maplot *mp, Keypos *kp)
{
  static Cmpfnd cf;  /* Returned descriptor */
/*
 * Nothing found yet.
 */
  cf.mod = NULL;
  cf.prev = NULL;
  cf.cmp = NULL;
/*
 * Only search displayed models.
 */
  if(mp->domod || mp->dovar) {
    int dofix=mp->domod;              /* Are fixed components displayed? */
    int dovar=mp->domod || mp->dovar; /* Are variable components displayed? */
/*
 * Search the established model.
 */
    cf_search(&cf, kp, mp->ob->model, dofix, dovar);
/*
 * Search the tentative model.
 */
    cf_search(&cf, kp, mp->ob->newmod, dofix, dovar);
  };
/*
 * Search the list of displayed new model components.
 */
  cf_search(&cf, kp, mp->newmod, 1, 1);
/*
 * Anything found?
 */
  return cf.cmp ? &cf : NULL;
}

/*.......................................................................
 * Private work function of fnd_cmp() used to search a model for the
 * nearest component. If two components are equally close to the cursor
 * position, the later one is always chosen.
 *
 * Input:
 *  cf   Cmpfnd *  The search state descriptor, describing the current
 *                 state of the search.
 *  kp   Keypos *  The cursor selection descriptor.
 *  mod   Model *  The model to be searched.
 *  dofix   int    Only consider fixed components if dofix is true.
 *  dovar   int    Only consider variable components if dovar is true.
 */
static void cf_search(Cmpfnd *cf, Keypos *kp, Model *mod, int dofix, int dovar)
{
  Modcmp *cmp;       /* The latest component to be looked at */
  Modcmp *prev;      /* The component preceding cmp in the model list */
/*
 * Search the model for a closer component than cf->cmp.
 */
  prev = NULL;
  for(cmp=mod->head; cmp; prev=cmp, cmp = prev->next) {
    if((dofix && !cmp->freepar) || (dovar && cmp->freepar)) {
      float xoff = cmp->x - kp->x;
      float yoff = cmp->y - kp->y;
      float roff = sqrt(xoff * xoff + yoff * yoff);
      if(cf->cmp==NULL || roff <= cf->roff) {
	cf->roff = roff;
	cf->mod = mod;
	cf->prev = prev;
	cf->cmp = cmp;
      };
    };
  };
  return;
}

/*.......................................................................
 * Given a cursor selection descriptor, delete the model component displayed
 * closest to the cursor.
 *
 * Input:
 *  mp   Maplot *   The plot descriptor.
 *  kp   Keypos *   The cursor selection descriptor.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Fatal error.
 */
static int zapcmp(Maplot *mp, Keypos *kp)
{
  Cmpfnd *cf;  /* Pointer to a description of the nearest component */
  Modcmp *cmp; /* The component removed from the model */
/*
 * Locate the nearest displayed component to the cursor position described
 * in kp.
 */
  cf = fnd_cmp(mp, kp);
  if(cf==NULL)
    return 0;
/*
 * Extract a model from the scratch model?
 */
  if(cf->mod == mp->newmod) {
    cmp = rem_cmp(cf->mod, cf->prev, cf->cmp);
  } else {
    cmp = obremcmp(mp->ob, cf->cmp, cf->mod==mp->ob->model);
  };
/*
 * obremcmp() returns NULL if the component is not found or the model is
 * not in a modifiable state.
 */
  if(cmp) {
/*
 * Erase the component from the display.
 */
    cmpplot(cf->cmp, mp->wxa, mp->wxb, mp->wya, mp->wyb, 1);
/*
 * Delete the component from the given model.
 */
    cf->cmp = del_cmp(cmp);
/*
 * Mark the map as invalid unless the component was in mp->newmod.
 */
    if(cf->mod != mp->newmod)
      mp->mb->domap = 1;
  };
  return 0;
}

/*.......................................................................
 * Find the nearest marker to the cursor.
 *
 * Input:
 *  mp   Maplot *   The plot descriptor.
 *  kp   Keypos *   The cursor selection descriptor.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Fatal error.
 */
static MarkerNode *find_marker(Maplot *mp, Keypos *kp)
{
  MarkerNode *marker;         /* The marker being checked */
  MarkerNode *nearest = NULL; /* The nearest marker to the cursor */
  float minrsqr = 0.0;        /* The distance of the closest marker to kp */
/*
 * Get the parent observation.
 */
  Observation *ob = mp->ob;
/*
 * Get the source description container.
 */
  Source *src = &ob->source;
/*
 * Get the RA, Dec coordinate projection type.
 */
  Proj proj = ob->proj;
/*
 * Get the container that records the current offset of the map center.
 */
  UVgeom *geom = &ob->geom;
/*
 * No markers?
 */
  if(!mp->markers)
    return NULL;
/*
 * Scan the list of map markers looking for the closest.
 */
  for(marker=mp->markers->head; marker; marker=marker->next) {
/*
 * Get the XY coordinates of the marker wrt the map center.
 */
    float x = geom->east +
      radec_to_l(src->ra, src->dec, marker->ra, marker->dec, proj);
    float y = geom->north +
      radec_to_m(src->ra, src->dec, marker->ra, marker->dec, proj);
/*
 * Work out the X and Y offsets of the marker from the cursor.
 */
    float dx = x - kp->x;
    float dy = y - kp->y;
/*
 * Work out the square of the distance of the marker from the cursor position.
 */
    float rsqr = dx * dx + dy * dy;
/*
 * Is this the nearest marker to the cursor so far?
 */
    if(!nearest || rsqr < minrsqr) {
      nearest = marker;
      minrsqr = rsqr;
    };
  };
  return nearest;
}

/*.......................................................................
 * Given a cursor selection descriptor, delete the marker displayed
 * closest to the cursor.
 *
 * Input:
 *  mp   Maplot *   The plot descriptor.
 *  kp   Keypos *   The cursor selection descriptor.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Fatal error.
 */
static int zapmark(Maplot *mp, Keypos *kp)
{
/*
 * Find the nearest marker to the cursor.
 */
  MarkerNode *nearest = find_marker(mp, kp);
/*
 * Delete it from the marker list.
 */
  nearest = del_MarkerNode(mp->markers, nearest);
/*
 * Redraw the map.
 */
  return replot(mp);
}

/*.......................................................................
 * Change the installed colormap.
 *
 * Input:
 *  mp     Maplot *   The plot descriptor.
 *  class Cmclass     The class of colormap to change to:
 *                      CM_COLOR  -  A psuedo-color colormap.
 *                      CM_GREY   -  A grey-scale colormap.
 *  ask       int     If true, prompt the user for the name of the
 *                    colormap to be installed.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Fatal error.
 */
static int change_cmap(Maplot *mp, Cmclass class, int ask)
{
  char answer[80];  /* Buffer for user reply */
  Ctable *ctab;     /* Color table descriptor */
  int doflip;       /* If true, flip the sense of the color table */
/*
 * Get the color table descriptor.
 */
  ctab = mp->cmpar.ctab;
/*
 * The colormap is only relevant when already plotting grey-scale
 * or psuedo-color.
 */
  if(mp->mono || ctab->cmap->class == CM_NONE)
    return 0;
/*
 * Ask the user?
 */
  if(ask) {
    printf("Enter the name of a color map: ");
    if(fgets(answer, sizeof(answer), stdin)) {
/*
 * Strip leading and trailing white space and lookup the named colormap.
 */
      stripstr(answer, sizeof(answer));
      if(*answer == '\0' || get_Cmap(ctab, answer) == NULL)
	return 0;
    } else {
      return 0;
    };
/*
 * Override the input colormap class with that of the selected colormap.
 */
    class = ctab->cmap->class;
  };
/*
 * Install the new colormap if necessary.
 */
  if(class != ctab->cmap->class) {
    switch(class) {
    case CM_GREY:
      get_Cmap(ctab, "grey");
      break;
    case CM_COLOR:
      get_Cmap(ctab, "color");
      break;
    case CM_NONE:
      return 0;
      break;
    };
  };
/*
 * For grey-scale, reverse the ramp on grey-scale devices. This stops
 * one from getting a sheet of paper covered in black.
 */
  doflip = mp->hard && ctab->cmap->class==CM_GREY;
/*
 * Display the new colormap.
 */
  if(doflip)
    recolor(ctab->cmap, -ctab->contra, 1.0-ctab->bright);
  else
    recolor(ctab->cmap, ctab->contra, ctab->bright);
  return 0;
}

/*.......................................................................
 * Re-display the current image with a different transfer function.
 *
 * Input:
 *  mp     Maplot *   The plot descriptor.
 *  tran   Cmtran     The type of transfer function.
 *                      TR_LINEAR - Linear transfer function.
 *                      TR_LOG    - Logarithmic transfer function.
 *                      TR_SQRT   - Square root transfer function.
 *  ask       int     If true, prompt the user for the name of the
 *                    transfer function to be used.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Fatal error.
 */
static int change_transfer(Maplot *mp, Cmtran tran, int ask)
{
  char answer[80];  /* Buffer for user reply */
  Ctable *ctab;     /* Color table descriptor */
/*
 * Get the color table descriptor.
 */
  ctab = mp->cmpar.ctab;
/*
 * The transfer function is only relevant when plotting grey-scale
 * or psuedo-color.
 */
  if(ctab->cmap->class == CM_NONE)
    return 0;
/*
 * Ask the user?
 */
  if(ask) {
    printf("Enter the name of a transfer function: ");
    if(fgets(answer, sizeof(answer), stdin)) {
/*
 * Strip leading and trailing white space and lookup the named
 * transfer function.
 */
      stripstr(answer, sizeof(answer));
      if(*answer == '\0')
	return 0;
      ctab->tran = get_Cmtran(answer);
    } else {
      return 0;
    };
  } else {
    ctab->tran = tran;
  };
/*
 * Re-display with the new transfer function.
 */
  return replot(mp);
}

/*.......................................................................
 * Display a polarization vector plot from the polarized intensity and
 * angle maps that precede and follow the normal map.
 * NB. setport() must be called before this function.
 *
 * Input:
 *  mp  Maplot *  The fully initialized plot descriptor.
 * Output:
 *  return int    0 - OK.
 *                1 - Error.
 */
static int plvect(Maplot *mp)
{
  MaplotVect *vect;    /* Local pointer to *mp->vect */
  MapBeam *mb;         /* The map-beam container object */
  float *magptr;       /* A pointer into the polarization magnitude map */
  float *angptr;       /* A pointer into the polarization angle map */
  float *mapptr;       /* A pointer into the unpolarized map */
  int xa,xb,ya,yb;     /* The range of pixels within the maps to be displayed */
  int ix,iy;           /* Map indexes along x and y */
  int xskip;           /* The number of pixels not being plotted in x */
  int dx,dy;           /* The increments between x,y pixels in which to plot */
/*
 * Sanity check.
 */
  if(bad_Maplot("plcont", mp))
    return 1;
/*
 * Get a pointer to the polarization vector attributes description object.
 */
  vect = mp->vect;
/*
 * Get the increment between pixels at which to plot vectors.
 */
  dx = vect->dx > 0 ? vect->dx : 1;
  dy = vect->dy > 0 ? vect->dy : 1;
/*
 * Get a pointer to the map/beam container.
 */
  mb = mp->mb;
/*
 * Don't display the vectors until all have been drawn.
 */
  cpgbbuf();
/*
 * Compute the range of pixels to be displayed within the map.
 */
  xa = mp->pxa - mb->nx/4;
  xb = mp->pxb - mb->nx/4;
  ya = mp->pya - mb->ny/4;
  yb = mp->pyb - mb->ny/4;
/*
 * Compute the number of X-axis pixels of the inner quarter of the
 * maps to skip.
 */
  xskip = mb->nx/2 - (xb - xa + 1);
/*
 * Get pointers to the first pixels to be drawn from the polarized intensity
 * and magnitude maps.
 */
  magptr = mb->map + xa + ya * mb->nx/2;
  angptr = mb->map + xa + 3*mb->ny/4 * mb->nx + ya * mb->nx/2;
/*
 * Also get a pointer to the first pixel of the normal map.
 */
  mapptr = mb->map + mp->pxa + mp->pya * mb->nx;
/*
 * Draw the vectors.
 */
  for(iy=0; iy<=yb-ya; iy++, mapptr += mb->nx/2 + xskip, magptr += xskip,
      angptr += xskip) {
/*
 * Determine the world Y-axis position of the pixel.
 */
    float y = mp->wya + iy * mb->yinc;  /* World y-coordinate of pixel */
/*
 * See if it is ok to plot any vectors in this Y-axis pixel.
 */
    int yok = iy % dy == 0;
/*
 * Draw vectors within the current line of pixels.
 */
    for(ix=0; ix<=xb-xa; ix++) {
      float mag = *magptr++;            /* Polarized magnitude */
      float ang = *angptr++;            /* Polarized angle */
      float map = *mapptr++;            /* Unpolarized intensity */
      if(yok && ix % dx == 0 && mag > vect->pcut && fabs(map) > vect->icut) {
/*
 * Precompute the sine and cosine of the polarization angle.
 */
	float sin_ang = sin(ang);
	float cos_ang = cos(ang);
/*
 * Work out the scaled length of the vector divided by two.
 */
	float half_len = mag * vect->scale / 2.0;
/*
 * Get the world x-axis position of the pixel.
 */
	float x = mp->wxa + ix * mb->xinc;
/*
 * Draw the vector with its mid point at the middle of the pixel.
 */
	cpgmove(x - half_len * sin_ang, y - half_len * cos_ang);
	cpgdraw(x + half_len * sin_ang, y + half_len * cos_ang);
      };
    };
  };
/*
 * Reveal the vector plot.
 */
  cpgebuf();
  return 0;
}

/*.......................................................................
 * If provided, plot the list of markers in mp->markers.
 *
 * Input:
 *  mp    Maplot *  The resource object of the mapplot command.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static int plmarkers(Maplot *mp)
{
/*
 * Get the parent observation.
 */
  Observation *ob = mp->ob;
/*
 * Get the source description container.
 */
  Source *src = &ob->source;
/*
 * Get the RA, Dec coordinate projection type.
 */
  Proj proj = ob->proj;
/*
 * Get the container that records the current offset of the map center.
 */
  UVgeom *geom = &ob->geom;
/*
 * Are there any markers to be plotted?
 */
  if(mp->markers) {
    MarkerNode *marker;  /* The marker being plotted */
/*
 * Buffer output until the markers have all been drawn.
 */
    cpgbbuf();
/*
 * Traverse the list of markers, plotting them one at a time.
 */
    for(marker=mp->markers->head; marker; marker=marker->next) {
      float xbox[4], ybox[4]; /* The bounding box of displayed text */
      float xt,yt;  /* The coordinates of the start of the annotation string */
      float w, h;   /* The width and height of the text */
      float wc, hc; /* The width and height of a representative character */
      float chgap;  /* The horizontal gap between the text and the marker, */
                    /*  measured in characters. */
      int nchar;    /* The number of characters in the annotation text */
/*
 * Work out the offset of the marker RA, Dec from the center of
 * the map.
 */
      float x = geom->east +
	radec_to_l(src->ra, src->dec, marker->ra, marker->dec, proj);
      float y = geom->north +
	radec_to_m(src->ra, src->dec, marker->ra, marker->dec, proj);
/*
 * Set the requested character size.
 */
      cpgsch(marker->size);
/*
 * Work out the drawn width of a represetative character.
 */
      cpgqtxt(x, y, 0.0, 0.0, "X", xbox, ybox);
      wc = (xbox[2] - xbox[0]);
      hc = (ybox[1] - ybox[0]);
/*
 * Work out the drawn width and height of the annotation string.
 */
      if(marker->text) {
	cpgqtxt(x, y, 0.0, 0.0, marker->text, xbox, ybox);
	w = (xbox[2] - xbox[0]);
	h = (ybox[1] - ybox[0]);
	nchar = strlen(marker->text);
      } else {
	w = wc;
	h = hc;
	nchar = 0;
      };
/*
 * Work out the coordinates at which the start of the string should be
 * drawn.
 */
      xt = x - w * marker->just + wc * marker->xpos;
      yt = y - h * 0.5 + hc * marker->ypos;
/*
 * Work out the horizontal gap between the text and the marker, measured
 * in characters. If there is no gap, record the gap as zero.
 */
      chgap = marker->xpos - marker->just * nchar;
      if(chgap <= 0.0 && chgap >= -nchar)
	chgap = 0.0;
      else if(chgap < -nchar)
	chgap += nchar;
/*
 * Plot the symbol.
 */
      cpgsci(marker->color);
      switch(marker->sym) {
      case MK_ARROW:
	if(chgap==0.0 || fabs(marker->ypos) > 2.0) {
	  if(marker->ypos < 0.0)
	    cpgarro(xt + w/2.0, yt + hc, x, y);
	  else
	    cpgarro(xt + w/2.0, yt - hc/2.0, x, y);
	} else if(chgap < 0.0) {
	  cpgarro(xt + w + wc/2.0, yt + h/2.0, x, y);
	} else {
	  cpgarro(xt - wc/2.0, yt + h/2.0, x, y);
	};
	break;
      default:
	cpgpt(1, &x, &y, marker->sym);
      };
/*
 * Plot the annotation text, if there is any.
 */
      if(marker->text)
	cpgptxt(xt, yt, 0.0, 0.0, marker->text);
    };
/*
 * Reveal the markers.
 */
    cpgebuf();
  };
  return 0;
}

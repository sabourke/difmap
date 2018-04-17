/*
  This file contains user functions which call pgplot for display
  of data and user accessible variables concerned with these functions.
*/

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <math.h>

#include "sphere.h"
#include "helpdir.h"
#include "utils.h"
#include "cpgplot.h"

/*
  Declare variables that are to be aliased as user variables below. Only
  float, integer, char and logical variables are supported.
  NB. character strings must NOT be initialised here unless they are marked
  as R_ONLY parameters. This is to allow variable length strings where
  the previous string is often free'd first on the assumption that the
  memory for the string was allocated using malloc(), not by the
  compiler).
*/

extern float pi;
static char plot_open=0;
static char has_cursor=0;
static int ndiv=1;
static float xrange[2];
static float yrange[2];
static float zrange[2];
static float longlat[3];

static Descriptor plotv_type[] = {
  {'l' , '0' ,R_ONLY ,1, {1,1,1}, &plot_open},
  {'l' , '0' ,R_ONLY ,1, {1,1,1}, &has_cursor},
  {'i' , '0' ,NO_DEL ,1, {1,1,1}, &ndiv},
  {'f' , '1' ,NO_DEL ,2, {2,1,1}, &xrange},
  {'f' , '1' ,NO_DEL ,2, {2,1,1}, &yrange},
  {'f' , '1' ,NO_DEL ,2, {2,1,1}, &zrange},
  {'f' , '1' ,R_ONLY ,3, {3,1,1}, &longlat}
};

/*
  In the same order as the above array of types, define the array
  of user names for the arrays.
*/

static char *plotv_name[] = {
  "plot_open",
  "has_cursor",
  "ndiv",
  "xrange",
  "yrange",
  "zrange",
  "longlat"
};

/*
  Declare the user functions here.
*/

static Template(opdev_fn);
static Template(pgbox_fn);
static Template(advance_fn);
static Template(pgpap_fn);
static Template(pgdraw_fn);
static Template(pgline_fn);
static Template(pgerrx_fn);
static Template(pgerry_fn);
static Template(pggray_fn);
static Template(pgmove_fn);
static Template(pgpt_fn);
static Template(pglab_fn);
static Template(contour_fn);
static Template(grey_fn);
static Template(pgsci_fn);
static Template(pghist_fn);
static Template(pgbbuf_fn);
static Template(pgebuf_fn);
static Template(cursor_fn);
static Template(cursran_fn);
static Template(xyz_plot);
static Template(pgslw_fn);
static Template(pgtext_fn);
static Template(pgptext_fn);
static Template(pgmtxt_fn);
static Template(pgrect_fn);
static Template(pgsch_fn);
static Template(pgsls_fn);
static Template(window_fn);
static Template(axes_fn);
static Template(lgraph_fn);
static Template(pgraph_fn);
static Template(tvflag_fn);
static Template(pgscr_fn);
static Template(pgshls_fn);
static Template(pgscf_fn);
static Template(pgsave_fn);
static Template(pgunsa_fn);
static Template(pgarro_fn);
static Template(pgask_fn);
static Template(pgband_fn);
static Template(pgcirc_fn);
static Template(pgcont_fn);
static Template(pgcurs_fn);
static Template(pgend_fn);
static Template(pgenv_fn);
static Template(pgeras_fn);
static Template(pgerrb_fn);
static Template(pgscir_fn);
static Template(pgscrn_fn);
static Template(pgsfs_fn);
static Template(pgshs_fn);
static Template(pgswin_fn);
static Template(pgsvp_fn);
static Template(pgvstd_fn);
static Template(pgpage_fn);
static Template(pgpoly_fn);
static Template(pgwnad_fn);

static Template(pgqah_fn);
static Template(pgqcf_fn);
static Template(pgqch_fn);
static Template(pgqci_fn);
static Template(pgqcir_fn);
static Template(pgqcol_fn);
static Template(pgqcr_fn);
static Template(pgqcs_fn);
static Template(pgqfs_fn);
static Template(pgqhs_fn);
static Template(pgqid_fn);
static Template(pgqitf_fn);
static Template(pgqls_fn);
static Template(pgqlw_fn);
static Template(pgqpos_fn);
static Template(pgqvp_fn);
static Template(pgqvsz_fn);
static Template(pgqwin_fn);

/*
  Declare the function types below.
*/

static Functype plotf_type[] = {
   {opdev_fn,  NORM   , 0,3,    " Cf",  " 00",    " vv",  1 },
   {advance_fn,NORM   , 0,0,    " ",    " ",      " ",    1 },
   {pgpap_fn,  NORM   , 2,2,    " ff",  " 00",    " vv",  1 },
   {pgdraw_fn, NORM   , 2,2,    " ff",  " 00",    " vv",  1 },
   {pgline_fn, NORM   , 2,2,    " ff",  " 11",    " vv",  1 },
   {pgmove_fn, NORM   , 2,2,    " ff",  " 00",    " vv",  1 },
   {pgpt_fn,   NORM   , 2,3,    " ffi", " 110",   " vvv", 1 },
   {pglab_fn,  NORM   , 3,3,    " ccc", " 000",   " vvv", 1 },
   {contour_fn,NORM   , 2,6,    " fff", " 211",   " vvv",  1 },
   {grey_fn,   NORM   , 1,7,    " ffff", " 2001"," vvvv", 1 },
   {pgsci_fn,  NORM   , 1,1,     " i",    " 0",     " v",   1 },
   {pghist_fn, NORM   , 4,5,    " fffil"," 30000",  " vvvvv",1 },
   {pgbbuf_fn, NORM   , 0,0,    " ",    " ",      " ",    1 },
   {pgebuf_fn, NORM   , 0,0,    " ",    " ",      " ",    1 },
   {cursor_fn, NORM   , 2,2,    " ff",  " 00",    " rr",  1 },
   {cursran_fn,NORM   , 2,2,    " fi",  " 11",    " vN",  1 },
   {tvflag_fn, NORM   , 2,3,    "lffl", "1110",   "vvvv", 1 },
   {xyz_plot,  NORM   , 3,6,    " ffff"," 1110",  " vvvv",1 },
   {window_fn, NORM   , 2,5,     " fffff"," 11100", " vvvvv",1 },
   {axes_fn,   NORM   , 2,5,     " fffff"," 11100", " vvvvv",1 },
   {lgraph_fn, NORM   , 2,2,     " ff",   " 11",    " vv",1 },
   {pgraph_fn, NORM   , 2,3,     " ffi",  " 110",   " vvv",1 },
   {pgarro_fn, NORM   , 4,4,     " ffff", " 0000",  " vvvv",1 },
   {pgask_fn,  NORM   , 1,1,     " l",    " 0",     " v",   1 },
   {pgband_fn, NORM   , 7,7,     "iiiffffc", "00000000", "?vvvvrrr",   1 },
   {pgbbuf_fn, NORM   , 0,0,     " ",     " ",      " ",    1 },
   {pgbox_fn,  NORM   , 0,6,     " CfiCfi"," 000000"," vvvvvv",  1 },
   {pgcirc_fn, NORM   , 3,3,     " fff",  " 000",   " vvv", 1 },
   {pgcont_fn, NORM   , 7,7,     " fiiiiff"," 2000011"," vvvvvvv", 1},
   {pgcurs_fn, NORM   , 3,3,     "iffc",  "0000",   "?rrr", 1 },
   {pgdraw_fn, NORM   , 2,2,     " ff",   " 00",    " vv",  1 },
   {pgebuf_fn, NORM   , 0,0,     " ",     " ",      " ",    1 },
   {pgend_fn,  NORM   , 0,0,     " ",     " ",      " ",    1 },
   {pgenv_fn,  NORM   , 6,6,     " ffffii"," 000000", " vvvvvv", 1 },
   {pgeras_fn, NORM   , 0,0,     " ",     " ",      " ",    1 },
   {pgerrb_fn, NORM   , 4,5,     " iffff"," 01110",  " vvvvv",1 },
   {pgerrx_fn, NORM   , 3,4,     " ffff"," 1110",  " vvvv",1 },
   {pgerry_fn, NORM   , 3,4,     " ffff"," 1110",  " vvvv",1 },
   {pggray_fn, NORM   , 8,8,     " fiiiifff", " 20000001", " vvvvvvvv", 1 },
   {pghist_fn, NORM   , 4,5,     " fffil"," 30000",  " vvvvv",1 },
   {pglab_fn,  NORM   , 3,3,     " ccc", " 000",   " vvv", 1 },
   {pgline_fn, NORM   , 2,2,     " ff",  " 11",    " vv",  1 },
   {pgmove_fn, NORM   , 2,2,     " ff",  " 00",    " vv",  1 },
   {pgmtxt_fn, NORM   , 5,5,     " Cfffc"," 00000"," vvvvv", 1 },
   {pgmtxt_fn, NORM   , 5,5,     " Cfffc"," 00000"," vvvvv", 1 },
   {pgpage_fn, NORM   , 0,0,     " ",     " ",      " ",    1 },
   {pgpap_fn,  NORM   , 2,2,     " ff",  " 00",    " vv",  1 },
   {pgpoly_fn, NORM   , 2,2,     " ff",   " 11",   " vv",  1 },
   {pgpt_fn,   NORM   , 2,3,     " ffi", " 110",   " vvv", 1 },
   {pgptext_fn,NORM   , 5,5,     " ffffc"," 00000"," vvvvv", 1 },
   {pgptext_fn,NORM   , 5,5,     " ffffc"," 00000"," vvvvv", 1 },
   {pgrect_fn, NORM   , 4,4,     " ffff"," 0000", " vvvv",1 },
   {pgsave_fn, NORM   , 0,0,     " ",     " ",      " ",    1 },
   {pgunsa_fn, NORM   , 0,0,     " ",     " ",      " ",    1 },
   {pgscf_fn,  NORM   , 1,1,     " i",    " 0",     " v",    1 },
   {pgsch_fn,  NORM   , 1,1,     " f",   " 0",    " v",   1 },
   {pgsci_fn,  NORM   , 1,1,     " i",   " 0",     " v",   1 },
   {pgscir_fn, NORM   , 2,2,     " ii",  " 00",    " vv",  1 },
   {pgscr_fn,  NORM   , 4,4,     " ifff", " 0000",  " vvvv", 1 },
   {pgscrn_fn, NORM   , 2,3,     " iCi", " 000",   " vvv", 1 },
   {pgsfs_fn,  NORM   , 1,1,     " i",   " 0",     " v",   1 },
   {pgshls_fn, NORM   , 4,4,     " ifff", " 0000",  " vvvv", 1 },
   {pgshs_fn,  NORM   , 3,3,     " fff", " 000",   " vvv", 1 },
   {pgsls_fn,  NORM   , 1,1,     " i",   " 0",    " v",   1 },
   {pgslw_fn,  NORM   , 1,1,     " i",   " 0",     " v",   1 },
   {pgsvp_fn,  NORM   , 4,4,     " ffff"," 0000",  " vvvv",1 },
   {pgswin_fn, NORM   , 4,4,     " ffff"," 0000",  " vvvv",1 },
   {pgtext_fn, NORM   , 3,3,     " ffc"," 000",   " vvv", 1 },
   {pgvstd_fn, NORM   , 0,0,     " ",    " ",      " ",    1 },
   {pgwnad_fn, NORM   , 4,4,     " ffff",  " 0000",  " vvvv", 1},
   {pgqah_fn,  NORM   , 3,3,     " iff",   " 000",   " rrr",  1},
   {pgqcf_fn,  NORM   , 1,1,     " i",     " 0",     " r",    1},
   {pgqch_fn,  NORM   , 1,1,     " f",     " 0",     " r",    1},
   {pgqci_fn,  NORM   , 1,1,     " i",     " 0",     " r",    1},
   {pgqcir_fn, NORM   , 2,2,     " ii",    " 00",    " rr",   1},
   {pgqcol_fn, NORM   , 2,2,     " ii",    " 00",    " rr",   1},
   {pgqcr_fn,  NORM   , 4,4,     " ifff",  " 0000",  " vrrr", 1},
   {pgqcs_fn,  NORM   , 3,3,     " iff",   " 000",   " vrr",  1},
   {pgqfs_fn,  NORM   , 1,1,     " i",     " 0",     " r",    1},
   {pgqhs_fn,  NORM   , 3,3,     " fff",   " 000",   " rrr",  1},
   {pgqid_fn,  NORM   , 1,1,     " i",     " 0",     " r",    1},
   {pgqitf_fn, NORM   , 1,1,     " i",     " 0",     " r",    1},
   {pgqls_fn,  NORM   , 1,1,     " i",     " 0",     " r",    1},
   {pgqlw_fn,  NORM   , 1,1,     " i",     " 0",     " r",    1},
   {pgqpos_fn, NORM   , 2,2,     " ff",    " 00",    " rr",   1},
   {pgqvp_fn,  NORM   , 5,5,     " iffff", " 00000", " vrrrr",1},
   {pgqvsz_fn, NORM   , 5,5,     " iffff", " 00000", " vrrrr",1},
   {pgqwin_fn, NORM   , 4,4,     " ffff",  " 0000",  " rrrr", 1},
};

/*
  In the same order as the above array of types, define the array
  of user names for the functions.
*/

static char *plotf_name[] = {
   "device",
   "advance",
   "paper",
   "draw",
   "line",
   "move",
   "point",
   "label",
   "contour",
   "grey",
   "colour",
   "histogram",
   "bbuf",
   "ebuf",
   "cursor",
   "cursor_range",
   "tvflag",
   "xyz_plot",
   "pgwindow",
   "box",
   "lplot",
   "pplot",
   "pgarro",
   "pgask",
   "pgband",
   "pgbbuf",
   "pgbox",
   "pgcirc",
   "pgcont",
   "pgcurs",
   "pgdraw",
   "pgebuf",
   "pgend",
   "pgenv",
   "pgeras",
   "pgerrb",
   "pgerrx",
   "pgerry",
   "pggray",
   "pghist",
   "pglab",
   "pgline",
   "pgmove",
   "pgmtext",
   "pgmtxt",
   "pgpage",
   "pgpap",
   "pgpoly",
   "pgpt",
   "pgptext",
   "pgptxt",
   "pgrect",
   "pgsave",
   "pgunsa",
   "pgscf",
   "pgsch",
   "pgsci",
   "pgscir",
   "pgscr",
   "pgscrn",
   "pgsfs",
   "pgshls",
   "pgshs",
   "pgsls",
   "pgslw",
   "pgsvp",
   "pgswin",
   "pgtext",
   "pgvstd",
   "pgwnad",
   "pgqah",
   "pgqcf",
   "pgqch",
   "pgqci",
   "pgqcir",
   "pgqcol",
   "pgqcr",
   "pgqcs",
   "pgqfs",
   "pgqhs",
   "pgqid",
   "pgqitf",
   "pgqls",
   "pgqlw",
   "pgqpos",
   "pgqvp",
   "pgqvsz",
   "pgqwin",
};

static EXITFN(plot_end);  /* Module closedown function */

/*
  Record the above declarations etc for this module in a global
  structure for use when building the main symbol table.
*/

Module m_graphics = {
  "graphics",
  HELP_DIR,
  NULL, 0,
  plotv_type, plotv_name, COUNT(plotv_name),
  plotf_type, plotf_name, COUNT(plotf_name),
  0,
  plot_end
};

/*
  Local variables and functions.
*/

static int plot_limits(Descriptor *invals[],  int npar, float limits[4], float x_gap, float y_gap);

static int pgplot_newdev(char *name, int xnum, int ynum);

int make_open(void);
int check_open(void);
static int check_cursor(void);

/*.......................................................................
 * Plotlib closedown function.
 *
 * Input:
 *  code    Exitcode   DO_QUIT  -  Minimal fast closedown.
 *                     DO_EXIT  -  Full closedown.
 */
static EXITFN(plot_end)
{
  cpgend();    /* Close PGPLOT */
}

/*.......................................................................
  Open the display device named by the user, or ask pgplot to provide
  a selection of devices etc..
*/
static Template(opdev_fn)
{
        int xnum,ynum;  /* Number of plot subdivisions in X and Y */
	char *name;     /* Pointer to device name string */
/*
  Get the number of plot divisions if specified by the user.
  Default to one if not specified.
*/
	xnum = (npar>1) ? (int) *FLTPTR(invals[1]) : 1;
	ynum = (npar>2) ? (int) *FLTPTR(invals[2]) : 1;
/*
  Get a pointer to the device name to be used.
*/
	if(npar != 0)
	  name = *STRPTR(invals[0]);
	else
	  name = NULL;
/*
  Attempt to open the device and return the appropriate error code.
*/
	return pgplot_newdev(name, xnum, ynum);
}

/*.......................................................................
 * Call pgbox to draw plot axes. The user may specify optional
 * specifications for the axes, otherwise use defaults.
 */
static Template(pgbox_fn)
{
  char *xopt, *yopt;  /* The X and Y axis configuration options */
  float xtic,ytic;    /* The axis tick-mark interval */
  int nxsub,nysub;    /* The number of major ticks per axis */
/*
 * Device open?
 */
	if(make_open() == -1)
	  return -1;
/*
 *  Get the axis specifiers if given.
 */
  xopt = (npar > 0) ? *STRPTR(invals[0]) : "BCNST";
  xtic = (npar > 1) ? *FLTPTR(invals[1]) : 0.0;
  nxsub= (npar > 2) ? *INTPTR(invals[2]) : 0;
  yopt = (npar > 3) ? *STRPTR(invals[3]) : "BCNST";
  ytic = (npar > 4) ? *FLTPTR(invals[4]) : 0.0;
  nysub= (npar > 5) ? *INTPTR(invals[5]) : 0;
/*
 * Draw the box.
 */
  cpgbox(xopt,xtic,nxsub,yopt,ytic,nysub);
  return no_error;
}

/*.......................................................................
  Clear the plot display.
*/
static Template(advance_fn)
{
/*
  Device open?
*/
	if(make_open() == -1)
	  return -1;
/*
 * Start a new page.
 */
        cpgpage();
	cpgvstd();
	return no_error;
}

/*.......................................................................
 * Set the width and aspect ratio of the display surface.
 */
static Template(pgpap_fn)
{
  if(make_open() == -1)
    return -1;
  cpgpap(*FLTPTR(invals[0]), *FLTPTR(invals[1]));
  return no_error;
}

/*.......................................................................
 * Draw a line from the current cursor position to the co-ordinates
 * specified by the user.
 */
static Template(pgdraw_fn)
{
  if(check_open() == -1)
    return -1;
  cpgdraw(*FLTPTR(invals[0]), *FLTPTR(invals[1]));
  return no_error;
}

/*.......................................................................
 * Draw a line through the x-value and y-value arrays sent by the user.
 */
static Template(pgline_fn)
{
  int nvals;
/*
 * Device open?
 */
  if(check_open() == -1)
    return -1;
/*
 * Get the lengths of the two coordinate arrays.
 */
  nvals = invals[0]->adim[0];
  if(invals[1]->adim[0] != nvals) {
    lprintf(stderr, "The X and Y arrays differ in length\n");
    return -1;
  };
  cpgline(nvals, FLTPTR(invals[0]), FLTPTR(invals[1]));
  return no_error;
}

/*.......................................................................
 * Move the current cursor position to the co-ordinates specified by the
 * user.
 */
static Template(pgmove_fn)
{
  if(check_open() == -1)
    return -1;
  cpgmove(*FLTPTR(invals[0]), *FLTPTR(invals[1]));
  return no_error;
}

/*.......................................................................
 * Draw a point at each of the arrays of x and y coordinates.
 * The optional third argument will be used to specify the marker symbol.
 */
static Template(pgpt_fn)
{
  int nvals;    /* The number of points to be plotted */
  int marker;   /* The PGPLOT marker symbol to be plotted */
/*
 * Device open?
 */
  if(check_open() == -1)
    return -1;
/*
 * Set the number of points to the minimum of the dimensions of the
 * X and Y coordinate arrays.
 */
  nvals = (invals[0]->adim[0] > invals[1]->adim[0]) ? invals[0]->adim[0] : invals[1]->adim[0];
  marker = (npar > 2) ? *INTPTR(invals[2]) : 2;
  cpgpt(nvals, FLTPTR(invals[0]), FLTPTR(invals[1]), marker);
  return no_error;
}

/*.......................................................................
 * Write user specified axis labels and a title around the current plot.
 */
static Template(pglab_fn)
{
/*
 * Device open?
 */
  if(check_open() == -1)
    return -1;
  cpglab(*STRPTR(invals[0]), *STRPTR(invals[1]), *STRPTR(invals[2]));
  return no_error;
}

/*.......................................................................
  Draw a contour plot of a two dimensional data array, using the
  contour levels in the 1-D array of the second argument. Four more optional
  arguments determine the values to be associated with the 1st and last
  elements along each axis of the 2D data array. These arguments have
  the same options as window_fn, allowing one to directly specify bounds
  or indirect bounds from min/max values of 1D array arguments. If these
  arguments are not given then the default is to set the axes up to
  show element numbers.
*/
static Template(contour_fn)
{
        int nlevs, xnum, ynum;
	float tr[6]={0.0,1.0,0.0, 0.0,0.0,1.0};
/*
  Device open?
*/
	if(check_open() == -1)
	  return -1;
/*
  Determine the size of the user's 2-D array.
*/
	xnum=invals[0]->adim[0];
	ynum=invals[0]->adim[1];
/*
  And the number of contour levels specified.
*/
	nlevs=invals[1]->adim[0];
/*
  If the user has supplied arguments for associating axis limits
  with the grey scale, use them.
*/
	if(npar > 2) {
	  float limits[4];  /* Axis limits of contour plot */
	  if(plot_limits(&invals[2], npar-2, limits, 0.0, 0.0) == -1)
	    return -1;
/*
  Turn the limits into the transformation matrix required by pgplot.
  This matrix describes linear ramps in X and Y as:
  X(IX,IY) = tr[0] + IX * tr[1] + IY * tr[2];
  Y(IX,IY) = tr[3] + IX * tr[4] + IY * tr[5];
  Where I and J are fortran indexes into the arrays (starting at 1!!).
  We won't require the cross terms tr[2] and tr[4] - these remain zero.
*/
	  tr[1] = (limits[1] - limits[0]) / (xnum-1);
	  tr[0] = limits[0] - tr[1];
	  tr[5] = (limits[3] - limits[2]) / (ynum-1);
	  tr[3] = limits[2] - tr[5];
	};
/*
  Now contour the array.
*/
	cpgcont(FLTPTR(invals[0]), xnum, ynum, 1, xnum, 1, ynum,
		FLTPTR(invals[1]), nlevs, tr);
	return no_error;
}

/*.......................................................................
  Draw a grey-scale image of a two dimensional data array. Two optional
  arguments set the levels to be denoted by black and by white.
  If these are not given or identical then the minimum and maximum
  values in the data array will be used instead. Four more optional
  arguments determine the values to be associated with the 1st and last
  elements along each axis of the 2D data array. These arguments have
  the same options as window_fn, allowing one to directly specify bounds
  or indirect bounds from min/max values of 1D array arguments. If these
  arguments are not given then the default is to set the axes up to
  show element numbers.
*/
static Template(grey_fn)
{
        int xnum, ynum;	/* Number of points along each axis */
	int xynum;	/* Total number of points in array */
	float black;	/* Level to be shaded in black */
	float white;	/* Level to be shaded in white */
	float *inptr;	/* Pointer to start of data array */
	float limits[4];/* Axis limits to be associated with grey scale */
	float tr[6]={0.0,1.0,0.0, 0.0,0.0,1.0};
	int i;
/*
  Device open?
*/
	if(check_open() == -1)
	  return -1;
/*
  Determine the size of the user's 2-D array.
*/
	xnum = invals[0]->adim[0];
	ynum = invals[0]->adim[1];
	xynum = xnum * ynum;
/*
  Get a pointer to the start of the data array.
*/
	inptr = FLTPTR(invals[0]);
/*
  If user specified levels to be used for black and white, get them.
  Either both or neither levels must be supplied.
*/
	if(npar > 1 && npar < 3) {
	  lprintf(stderr, "Incomplete grey scale levels provided\n");
	  return -1;
	};
/*
  Get the levels.
*/
	black = white = 0.0;
	if(npar > 1) {
	  black = *FLTPTR(invals[1]);
	  white = *FLTPTR(invals[2]);
	};
/*
  If levels weren't provided or the levels that were, were equal, then
  use the min/max values in the data array.
*/
	if(npar == 1 || black == white) {
	  black = white = *inptr;
	  for(i=0; i<xynum; i++) {
	    if(*inptr > white)
	      white = *inptr;
	    else if(*inptr < black)
	      black = *inptr;
	    inptr++;
	  };
	};
/*
  If the data array is empty then black will still equal white.
*/
	if(black == white) {
	  lprintf(stderr, "Data array is uniform - auto-ranging failed\n");
	  return -1;
	};
/*
  If the user has supplied arguments for associating axis limits
  with the grey scale, use them.
*/
	if(npar > 3) {
	  if(plot_limits(&invals[3], npar-3, limits, 0.0, 0.0) == -1)
	    return -1;
/*
  Turn the limits into the transformation matrix required by pgplot.
  This matrix describes linear ramps in X and Y as:
  X(IX,IY) = tr[0] + IX * tr[1] + IY * tr[2];
  Y(IX,IY) = tr[3] + IX * tr[4] + IY * tr[5];
  Where I and J are fortran indexes into the arrays (starting at 1!!).
  We won't require the cross terms tr[2] and tr[4] - these remain zero.
*/
	  tr[1] = (limits[1] - limits[0]) / (xnum-1);
	  tr[0] = limits[0] - tr[1];
	  tr[5] = (limits[3] - limits[2]) / (ynum-1);
	  tr[3] = limits[2] - tr[5];
	};
/*
  Now (Finally!) grey-scale the array.
*/
	cpggray(FLTPTR(invals[0]), xnum, ynum, 1, xnum, 1, ynum, white, black, tr);
	return no_error;
}

/*.......................................................................
 * Change the current plotting colour to that specified by the user.
 */
static Template(pgsci_fn)
{
  int ci;
/*
 * Device open?
 */
  if(make_open() == -1)
    return -1;
/*
 * Extract the required colour index and check that it is within range.
 */
  ci = *INTPTR(invals[0]);
  if(ci > 255 || ci < 0) {
    lprintf(stderr, "Illegal colour: %d (0 - 255)\n",ci);
    return -1;
  };
/*
 * Apply the colour change.
 */
  cpgsci(ci);
  return no_error;
}

/*.......................................................................
 * Draw a histogram of a user array. The array is the first user argument,
 * and is binned by cpghist() into the number of bins specified via the
 * user's 4th argument. The min and max data values to be considered are
 * given by the user as the second and third arguments.
 */
static Template(pghist_fn)
{
  int nbins;           /* The number of histogram bins to plot */
  int nvals;           /* The number of data-points to bin */
  int no_clear;        /* If true, display the histogram in the current plot */
  float maxval,minval; /* The range of the histogram */
/*
 * Device open?
 */
  if(make_open() == -1)
    return -1;
/*
 * Find out the number of elements in the user array.
 */
  nvals = invals[0]->adim[0] * invals[0]->adim[1] * invals[0]->adim[2];
/*
 * Determine the min and max data values to be binned and check them.
 */
  minval = *FLTPTR(invals[1]);
  maxval = *FLTPTR(invals[2]);
  if(minval >= maxval) {
    lprintf(stderr,"Bad min=%f, max=%f given to histogram()\n", minval, maxval);
    return -1;
  };
/*
 * Determine the number of bins specified.
 */
  nbins = *INTPTR(invals[3]);
/*
 * Determine whether the user wants the histogram to overlay the current plot.
 */
  no_clear = (npar > 4) ? *LOGPTR(invals[4]) : 0;
/*
 * Plot the histogram.
 */
  cpghist(nvals, FLTPTR(invals[0]), minval, maxval, nbins, no_clear);
  return no_error;
}

/*.......................................................................
 * Turn pgplot buffering on.
 */
static Template(pgbbuf_fn)
{
/*
 * Device open?
 */
  if(check_open() == -1)
    return -1;
  cpgbbuf();
  return no_error;
}

/*.......................................................................
 * Turn pgplot buffering off.
 */
static Template(pgebuf_fn)
{
  if(check_open() == -1)
    return -1;
  cpgebuf();
  return no_error;
}

/*.......................................................................
  Given two arguments specifying the initial position of the cursor,
  bring the pgplot cursor up on the screen and return the position
  that the user specifies with a button press.
*/
static Template(cursor_fn)
{
        char ch[2];
	ch[1]='\0';
/*
  Device open?
*/
	if(check_open() == -1)
	  return -1;
/*
  Does the device have a cursor?
*/
	if(check_cursor() == -1)
	  return -1;
	if(cpgcurs(FLTPTR(invals[0]), FLTPTR(invals[1]), ch) == 0) {
	  lprintf(stderr, "Error getting cursor position.\n");
	  return -1;
	};
	return no_error;
}

/*.......................................................................
  Allow the user to delimit channel ranges along the x-axis of a plot
  by means of the pgplot cursor. The pairs of start and end ranges
  will be returned as a 1-D array. The user must send the x-array that
  was used to draw the current plot, and a reference to the return array.
*/
static Template(cursran_fn)
{
        int npts, nrange, endit, i,j,col,*outptr;
	float *inptr, ybot,ytop;
	float xmax,xmin,ymin,ymax,ycent,stem,bar,lastx,xpos,ypos;
	signed int bar_mult;
	char cdum[2];
	cdum[1]='\0';
/*
  Device open?
*/
	if(check_open() == -1)
	  return -1;
/*
  Does the current device have a cursor?
*/
	if(check_cursor() == -1)
	  return -1;
/*
  Find out the current pen colour.
*/
	cpgqci(&col);
/*
  Find the current plot's world co-ordinate limits and determine
  the y-position at which delimiting brackets will be plotted.
*/
	cpgqwin(&xmin, &xmax, &ymin, &ymax);
	ycent = (ymin+ymax)/2.0;
/*
  Determine the number of points in the x-array.
*/
	npts = invals[0]->adim[0];
/*
  Check to see if the data array is in ascending order.
*/
	inptr = FLTPTR(invals[0]);
	for(i=1; i<npts; i++)
	  if(inptr[i] < inptr[i-1]) {
	    lprintf(stderr, "x-array is not in ascending order.\n");
	    return -1;
	  };
/*
  Determine the max number of ranges that the return array can accomodate
  without redeclaration.
*/
	outptr = INTPTR(invals[1]);
	nrange = invals[1]->adim[0]/2;
	if(nrange == 0) {
	  lprintf(stderr, "No room in range return array.\n");
	  return -1;
	};
/*
  Set up the lengths of the square delimiter brackets in the plot
  world co-ordinates.
*/
	bar = (xmax-xmin)/200.0;
	stem = (ymax-ymin)/20.0;
	ytop = ycent + stem;
	ybot = ycent - stem;
/*
  Set the initial position of the cursor at the left side of the plot.
*/
	xpos = xmin;
	ypos = ycent;
	lastx = xmin-0.1;
	bar_mult = 1;
/*
  Get the cursor delimited ranges.
*/
	endit=0;
	i=0;
	do {
/*
  Get the next cursor position.
*/
	  if(cpgcurs(&xpos, &ypos, cdum) == 0) {
	    lprintf(stderr, "Error getting cursor position.\n");
	    return -1;
	  };
/*
  If the new position precedes the last position then delete the last
  delimiter and get ready to re-do that delimiter.
*/
	  if(xpos <= lastx && i>0) {
	    i--;
	    bar_mult *= -1;
/*
  Erase the delimiter by drawing over it in black.
*/
	    cpgsci(0);
	    cpgmove(lastx+bar * bar_mult, ybot);
	    cpgdraw(lastx, ybot);
	    cpgdraw(lastx, ytop);
	    cpgdraw(lastx+bar * bar_mult, ytop);
	    cpgsci(col);
	    lastx = (i>0) ? inptr[outptr[i-1]-1] : xmin-0.1;
	    continue;
	  };
/*
  If the cursor was pressed beyond the left edge of the window then
  default to that edge.
*/
	  if(xpos <= xmin)
	    xpos = xmin;
/*
  If the cursor was pressed beyond the right edge of the window then
  finish. If an unclosed range exists then use the value at the right edge
  for it.
*/
	  if(xpos >= xmax) {
	    xpos = xmax;
	    if(bar_mult == 1)
	      break;
	    endit=1;
	  };
/*
  Find the element of the x-array that lies closest to the cursor position.
*/
	  for(j=0; j<npts; j++) {
	    if(inptr[j] > xpos) {
	      if(j>0 && xpos - inptr[j-1] < inptr[j] - xpos)
		j--;
	      break;
	    };
	  };
	  if(j==npts) j--;
	  outptr[i] = j+1;
	  lastx = (j<npts) ? inptr[j]:inptr[npts-1];
/*
  Draw the range deliniter at the prescribed point.
*/
	  cpgmove(lastx+bar * bar_mult, ybot);
	  cpgdraw(lastx, ybot);
	  cpgdraw(lastx, ytop);
	  cpgdraw(lastx+bar * bar_mult, ytop);
	  bar_mult *= -1;
	  i++;
	} while(i < nrange * 2 && !endit);
/*
  If no range was set then report an error.
*/
	if(i<2) {
	  lprintf(stderr, "No limits set!\n");
	  return -1;
	};
/*
  Modify the dimension specification of the return array to reflect the
  number of elements actually used.
*/
	invals[1]->adim[0] = i;
	return no_error;
}

/*.......................................................................
  Given three 1-D arrays of x,y,z co-ordinates and three rotation
  angles about the x,y and z axes, display the data points rotated by the
  specified amounts, with perpective added to accentuate their 3-D
  distribution. Also change the size of the plotted points in relation to
  their post-rotation z-distance. The user must specify the rotation
  angles via the iterate variable, as {xmin,xmax,nxiter},{ymin...},{..}.
*/
static Template(xyz_plot)
{
        double x_angle=0.0, y_angle=0.0, z_angle=0.0;
	double cosx, cosy, cosz, sinx, siny, sinz;
	float xpos, ypos, zpos, ftmp, rmax, rmin, zinf, shrink;
	float xmid,ymid,zmid,zran, x_div, y_div, z_div;
	float *xvals, *yvals, *zvals, *xrot, *yrot, *zrot;
	extern float xrange[2], yrange[2], zrange[2];
	extern int ndiv;
	int i, j, num_pt, level, old_ci;
	int nlev=16;
	const float depth = 1.0;
	const double twopi = 6.2831853;
	const float rad_to_deg = 57.29578;  /* Conversion factor: radians to degrees */
/*
  Device open?
*/
	if(make_open() == -1)
	  return -1;
/*
  Get pointers to the data arrays.
*/
	xvals = FLTPTR(invals[0]);
	yvals = FLTPTR(invals[1]);
	zvals = FLTPTR(invals[2]);
/*
  Find the maximum number of data points of each user array argument. 
*/
	num_pt = invals[0]->adim[0];
	num_pt = (invals[1]->adim[0] > num_pt) ? invals[1]->adim[0] : num_pt;
	num_pt = (invals[2]->adim[0] > num_pt) ? invals[2]->adim[0] : num_pt;
/*
  Find the average x,y,z values using a running mean.
*/
	xmid=xvals[0]; ymid=yvals[0]; zmid=zvals[0];
	for(i=0; i<num_pt; i++) {
	  ftmp = 1.0/((float) i+1);
	  xmid += (xvals[i]-xmid) * ftmp;
	  ymid += (yvals[i]-ymid) * ftmp;
	  zmid += (zvals[i]-zmid) * ftmp;
	};
/*
  If the user specified their own xmid,ymid or zmid then substitute
  them for those calculated above.
*/
	xmid = (npar > 3) ? *FLTPTR(invals[3]) : xmid;
	ymid = (npar > 4) ? *FLTPTR(invals[4]) : ymid;
	zmid = (npar > 5) ? *FLTPTR(invals[5]) : zmid;
/*
  Find the minimum and maximum radius of the x,y,z points relative
  to their average x,y,z values.
*/
	rmax = rmin = 0.0;
	for(i=0;i<num_pt;i++) {
	  xpos = xvals[i]-xmid;
	  ypos = yvals[i]-ymid;
	  zpos = zvals[i]-zmid;
	  ftmp = xpos*xpos + ypos*ypos + zpos*zpos;
	  rmax = (ftmp > rmax) ? ftmp : rmax;
	  rmin = (ftmp < rmin) ? ftmp : rmin;
	};
	rmax=sqrt(rmax);
	rmin=sqrt(rmin);
/*
  Enforce a mimimum value of 1.0e-4 on rmax.
*/
	rmax = (rmax > 1.0e-4) ? rmax : 1.0e-4;
/*
  Compute the vanishing point to be a constant times the maximum radius.
  Then compute the reciprocal of the maximum distance possible between
  a data point and the vanishing point.
*/
	zinf = depth * rmax;
	zran = 1.0/(zinf+rmax);
/*
  Set up a plot box.
*/
        cpgenv(-rmax,rmax,-rmax, rmax, 1, 0);
/*
  Enforce a mimimum of one iteration of angle increments.
*/
	ndiv = (ndiv > 1) ? ndiv : 1.0;
/*
  Compute the repeat interval for the rotation angles around each axis.
*/
	x_div = (xrange[1]-xrange[0])/ndiv;
	y_div = (yrange[1]-yrange[0])/ndiv;
	z_div = (zrange[1]-zrange[0])/ndiv;
/*
  Allocate two work arrays which will hold the array of rotated
  x,y co-ordinates.
*/
	if( (xrot = (float *) calloc(num_pt+1, sizeof(float))) == NULL) {
	  lprintf(stderr, "Insufficient memory for work arrays in xyz_plot().\n");
	  return -1;
	};
	if( (yrot = (float *) calloc(num_pt+1, sizeof(float))) == NULL) {
	  lprintf(stderr, "Insufficient memory for work arrays in xyz_plot().\n");
	  free(xrot);
	  return -1;
	};
	if( (zrot = (float *) calloc(num_pt+1, sizeof(float))) == NULL) {
	  lprintf(stderr, "Insufficient memory for work arrays in xyz_plot().\n");
	  free(xrot);
	  free(yrot);
	  return -1;
	};
/*
  Re-define pgplot colour indices 16 to 31 for a linear ramp in grey scale.
*/
	for(i=0;i<nlev;i++) {
	  ftmp = 0.25 + i * 0.75/(nlev-1.0);
	  cpgscr(i+16, ftmp,ftmp,ftmp);
	};
/*
  Keep a record of the current pgplot colour index.
*/
	cpgqci(&old_ci);
/*
  Turn pgplot buffering on.
*/
	cpgbbuf();
/*
  Repeat the plot over each of the angles over the ranges specified
  by the user via xrange,yrange and zrange.
*/
	for(j=0; j<ndiv; j++) {
/*
  Get the user angle specifications.
*/
	  x_angle = xrange[0] + j * x_div;
	  y_angle = yrange[0] + j * y_div;
	  z_angle = zrange[0] + j * z_div;
/*
  Pre-calculate the cos()'s and sin()'s.
*/
	  cosx = cos(x_angle);
	  cosy = cos(y_angle);
	  cosz = cos(z_angle);
	  sinx = sin(x_angle);
	  siny = sin(y_angle);
	  sinz = sin(z_angle);
/*
  For each point apply the rotation matrix and plot the point.
*/
	  for(i=0; i<num_pt; i++) {
/*
  Get the raw x,y,z co-ordinates for the latest point and reference
  them to their mean value.
*/
	    xpos = xvals[i]-xmid;
	    ypos = yvals[i]-ymid;
	    zpos = zvals[i]-zmid;
/*
  Get to the required logitude with a rotation about the y-axis
  then get to the required latitude with a rotation about the
  new x-axis.
  Rotate about the y-axis.
*/
	    ftmp = xpos * cosy - zpos * siny;
	    zpos = zpos * cosy + xpos * siny;
	    xpos = ftmp;
/*
  Then about the x-axis.
*/
	    ftmp = zpos * cosx - ypos * sinx;
	    ypos = ypos * cosx + zpos * sinx;
	    zpos = ftmp;
/*
  Then about the z-axis.
*/
	    ftmp = xpos * cosz - ypos * sinz;
	    ypos = ypos * cosz + xpos * sinz;
	    xpos = ftmp;
	    xrot[i] = xpos;
	    yrot[i] = ypos;
	    zrot[i] = zpos;
	  };
/*
  Plot the data points - first in white.
*/
	  for(i=0; i<num_pt; i++) {
/*
  Compute the shrinkage factor for perspective visualisation.
*/
	    shrink = (zrot[i]+zinf) * zran;
	    level = (int) (16.0 + shrink * 15.0);
	    level = (level > 31) ? 31:level;
/*
  Use white for the closest points, light grey for intermediate
  distance points and dark grey for the distant points.
*/
	    cpgsci(level);
/*
  Shrink the x and y co-ordinates by this amount.
*/
	    xpos = xrot[i];
	    ypos = yrot[i];
	    cpgmove(0.9*xpos,0.9*ypos);
	    cpgdraw(xpos,ypos);
	  };
/*
  Then plot them in black - turn bufferring on such that
  the erase that this effects won't appear until buffering
  is turned off. ie when the next round of points have been
  plotted.
*/
	  cpgebuf();
/*
  Before erasing the current iteration - check that the user
  hasn't pressed ctrl-c or some other interrupt signal.
*/
	  if(no_error) break;
/*
  Don't erase the final plot!
*/
	  if(j < ndiv-1) {
	    cpgbbuf();
	    cpgsci(0);
	    for(i=0; i<num_pt; i++) {
/*
  Compute the shrinkage factor for perspective visualisation.
*/
/*	      shrink = (zrot[i]+zinf) * zran; */
/*
  Shrink the x and y co-ordinates by this amount.
*/
	      xpos = xrot[i];
	      ypos = yrot[i];
	      cpgmove(0.9*xpos,0.9*ypos);
	      cpgdraw(xpos,ypos);
	    };
	  };
	};
/*
  Delete the work arrays.
*/
	free(xrot);
	free(yrot);
	free(zrot);
/*
  Reset the colour index to its entry value.
*/
	cpgsci(old_ci);
/*
  Record the logitude and latitude of the last rotation displayed for the
  user to examine via the user parameter variable longlat.
*/
	longlat[0] = rad_to_deg * (float) fmod(y_angle, twopi);
	longlat[1] = rad_to_deg * (float) fmod(x_angle, twopi);
	longlat[2] = rad_to_deg * (float) fmod(z_angle, twopi);
	return no_error;
}

/*.......................................................................
 * Change the line width attribute of pgplot.
 */
static Template(pgslw_fn)
{
  int lw;
/*
 * Device open?
 */
  if(make_open() == -1)
    return -1;
/*
 * Extract the required line width and enforce limits of 1-200.
 */
  lw = *INTPTR(invals[0]);
  lw = (lw > 0) ? ((lw<200) ? lw:200) : 1;
/*
 * Apply the line width change.
 */
  cpgslw(lw);
  return no_error;
}

/*.......................................................................
 * Write text horizontally at given world-coordinates on the current plot.
 */
static Template(pgtext_fn)
{
  if(check_open() == -1)
    return -1;
  cpgtext(*FLTPTR(invals[0]),*FLTPTR(invals[1]),*STRPTR(invals[2]));
  return no_error;
}

/*.......................................................................
 * Write text at a given angle, world-coordinates and with a given
 * justification, on the current plot.
 */
static Template(pgptext_fn)
{
  if(check_open() == -1)
    return -1;
  cpgptxt(*FLTPTR(invals[0]), *FLTPTR(invals[1]), *FLTPTR(invals[2]),
	  *FLTPTR(invals[3]), *STRPTR(invals[4]));
  return no_error;
}

/*.......................................................................
 * Write text at a given angle, world-coordinates and with a given
 * justification, on the current plot.
 */
static Template(pgmtxt_fn)
{
  char side;
/*
 * Device open?
 */
  if(check_open() == -1)
    return -1;
/*
 * Make sure that the side designator is one of those accepted by pgmtext.
 */
  side = *STRPTR(invals[0])[0];
  switch(side) {
  case 'B': case 'L': case 'T': case 'R':
    break;
  default:
    lprintf(stderr,"Side option \'%c\' not one of B L T R\n",side);
    return -1;
    break;
  };
  cpgmtxt(*STRPTR(invals[0]), *FLTPTR(invals[1]), *FLTPTR(invals[2]),
	  *FLTPTR(invals[3]), *STRPTR(invals[4]));
  return no_error;
}

/*.......................................................................
 *  Draw a filled rectangle on the current plot given x1,x2,y1,y2.
 */
static Template(pgrect_fn)
{
  if(check_open() == -1)
    return -1;
  cpgrect(*FLTPTR(invals[0]), *FLTPTR(invals[1]),
	  *FLTPTR(invals[2]), *FLTPTR(invals[3]));
  return no_error;
}

/*.......................................................................
 * Set the character height on the current plot.
 */
static Template(pgsch_fn)
{
  if(make_open() == -1)
    return -1;
  cpgsch(*FLTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Change the line style attribute of pgplot.
 */
static Template(pgsls_fn)
{
  int ls;
/*
 * Device open?
 */
  if(make_open() == -1)
    return -1;
/*
 * Extract the required line style and enforce limits of 1-5.
 */
  ls = *INTPTR(invals[0]);
  ls = (ls > 0) ? ((ls<6) ? ls:5) : 1;
/*
 * Apply the line style change.
 */
  cpgsls(ls);
  return no_error;
}

/*.......................................................................
  This function is used to call pgwindow to set the world co-ordinate
  system of subsequent plots. The limits for the x-axis come first. These
  are specified either directly via two scalar arguments, or indirectly
  via an array whose min and max values are determined and used. The
  y-axis limits are specified next with the same options. Thus this
  function takes between two and four arguments.
*/
static Template(window_fn)
{
	float limits[4]; 	/* Will hold the plot limits */
	int used_args; /* Number of arguments used in determining limits */
	const float x_gap=0.05; /* Fractional margin around plot along x */
	const float y_gap=0.05; /* Fractional margin around plot along y */
/*
  Device open?
*/
	if(make_open() == -1)
	  return -1;
/*
  Get the plot limits.
*/
	used_args = plot_limits(invals, npar, limits, x_gap, y_gap);
	if(used_args == -1)
	  return -1;
/*
  Apply the limits to set up the world co-ordinate system.
  If there is one remaining argument which is greater than 0
  set a window and justify the viewport such that the scales
  are equal on the view surface.
*/
	if(npar > used_args && *FLTPTR(invals[used_args]) > 0)
	  cpgwnad(limits[0],limits[1],limits[2],limits[3]);
	else
	  cpgswin(limits[0],limits[1],limits[2],limits[3]);
	return no_error;
}

/*.......................................................................
  This function calls both window_fn and pgbox_fn in order to plot the
  axes of a graph.
*/
static Template(axes_fn)
{
  if(advance_fn(invals, 0, outvals) == -1 ||
     window_fn(invals, npar, outvals) == -1 ||
     pgbox_fn(invals,0,outvals) == -1)
    return -1;
  return no_error;
}

/*.......................................................................
  This function calls window_fn, pgbox_fn and pgline_fn in order to plot the
  axes for and a line through a set of data points.
*/
static Template(lgraph_fn)
{
  if(advance_fn(invals, 0, outvals) == -1 ||
     window_fn(invals, 2, outvals) == -1 ||
     pgbox_fn(invals,0,outvals) == -1 ||
     pgline_fn(invals,2,outvals) == -1)
    return -1;
  return no_error;
}

/*.......................................................................
  This function calls window_fn, pgbox_fn and pgpt_fn in order to plot the
  axes for and points at a set of data points.
*/
static Template(pgraph_fn)
{
  if(advance_fn(invals, 0, outvals) == -1 ||
     window_fn(invals, 2, outvals) == -1 ||
     pgbox_fn(invals,0,outvals) == -1 ||
     pgpt_fn(invals,npar,outvals) == -1)
    return -1;
  return no_error;
}

/*.......................................................................
  Allow the user to select channels from the x-axis of a plot by means
  of the pgplot cursor. The selected channels will be returned as a 1-D
  array of flags. The user must send the x and y arrays to be plotted
  If an optional third argument is set to ommitted or set to false then
  the graph will be drawn, otherwise the current graph will be used.
*/
static Template(tvflag_fn)
{
        int npts, i, col;
	float xmax,xmin,ymin,ymax,stem,xpos,ypos;
	float ydif,xdif,xnorm,ynorm,dist_min,xydist;
	float *xptr,*yptr;
	char *flagptr;
	int el_min;
	char cdum[2];
	cdum[1]='\0';
/*
  Device open?
*/
	if(check_open() == -1)
	  return -1;
/*
  Determine the number of points in the x/y-arrays.
*/
	npts = invals[0]->adim[0];
	if(invals[1]->adim[0] != npts) {
	  lprintf(stderr, "Unequal numbers of points in the x and y arrays\n");
	  return -1;
	};
/*
  If the third argument was ommitted or is false, then plot the points.
*/
	if(npar < 3 || *LOGPTR(invals[2]) == 0) {
	  if(pgraph_fn(invals,2,outvals) == -1)
	    return -1;
	};
/*
  Find out the current pen colour.
*/
	cpgqci(&col);
/*
  Find the current plot's world co-ordinate limits.
*/
	cpgqwin(&xmin, &xmax, &ymin, &ymax);
/*
  xnorm and ynorm are used to scale the x and y distances from
  world coordinates to normalized coordinates.
*/
	xnorm = 1.0/(xmax-xmin);
	ynorm = 1.0/(ymax-ymin);
/*
  Set up the length of the select bars in the world co-ordinates
  of the plot.
*/
	stem = (ymax-ymin)/20.0;
/*
  Get pointers to the x and y data arrays.
*/
	xptr = FLTPTR(invals[0]);
	yptr = FLTPTR(invals[1]);
/*
  Allocate the return flag array.
*/
	if( (VOIDPTR(outvals) = valof_alloc(npts,'l')) == NULL)
	  return -1;
/*
  Fill in the return-descriptor items.
*/
	outvals->adim[0] = outvals->num_el = npts;
/*
  Get a pointer to the return array.
*/
	flagptr = VOIDPTR(outvals);
/*
  Set the initial position of the cursor at the left side of the plot.
*/
	xpos = xmin;
	ypos = (ymax-ymin)/2.0;
/*
  If the current plot device has no cursor, inform the user but
  don't treat this as an error - treat it as though the user had
  chosen to flag no elements. The flag array currently contains
  no flagged elements.
*/
	if(!has_cursor) {
	  lprintf(stderr, "tvflag: No cursor available - no points flagged - continuing.\n");
	  return no_error;
	};
/*
  Get the cursor delimited ranges.
*/
	for(;;) {
/*
  Get the next cursor position.
*/
	  if(cpgcurs(&xpos, &ypos, cdum) == 0) {
	    lprintf(stderr, "cursor_sel: Error getting cursor position.\n");
	    return -1;
	  };
/*
  If the cursor lies outside the plot then quit.
*/
	  if(xpos < xmin || xpos > xmax || ypos < ymin || ypos > ymax)
	    return no_error;
/*
  Search the x and y data arrays for the point closest to the cursor.
*/
	  dist_min = 10.0;
	  el_min = 0;
	  for(i=0; i<npts; i++) {
	    xdif = (xpos - xptr[i]) * xnorm;
	    ydif = (ypos - yptr[i]) * ynorm;
	    xydist = xdif * xdif + ydif * ydif;
	    if(xydist < dist_min) {
	      dist_min = xydist;
	      el_min = i;
	    };
	  };
/*
  Switch the state of the flag for that element.
*/
	  flagptr[el_min] = !flagptr[el_min];
/*
  If the flag is now on then draw a vertical bar in red,
  otherwise draw it in black to erase the bar that was there.
*/
	  if(flagptr[el_min])
	    cpgsci(2);
	  else
	    cpgsci(0);
/*
  Draw the bar through the selected point.
*/
	  cpgmove(xptr[el_min], yptr[el_min]-stem);
	  cpgdraw(xptr[el_min], yptr[el_min]+stem);
	  cpgsci(col);
	};
}

/*.......................................................................
  This function is not called directly by the user. It is used by
  functions that are, however and is passed up to four user arguments in
  the normal way. The arguments are used to determine plot limits.
  The limits for the x-axis come first. These
  are specified either directly via two scalar arguments, or indirectly
  via an array whose min and max values are determined and used. The
  y-axis limits are specified next with the same options.
  Margins of (x_gap * x_range) and (y_gap * y_range) will be added to the
  limits. The limits are returned in the 4 element array:
   limits[0:3]=xmin,xmax,ymin,ymax.
  If an error occurs, the function return value is -1, otherwise the return
  vaue will be the number of user arguments actually used.
*/
static int plot_limits(Descriptor *invals[],  int npar, float limits[4], float x_gap, float y_gap)
{
        int npts;	/* Number of points in an array */
	int arg;        /* Number of next unused argument */
	float vmin,vmax;/* Used while finding min/max values in an array */
	float *inptr;	/* Pointer into a user array */
	int i;
/*
  X-limits.
*/
	npts = invals[0]->adim[0];
/*
  Two scalars?
*/
	if(npts == 1) {
	  if(invals[1]->adim[0] != 1) {
	    lprintf(stderr, "First x-limit scalar but not second.\n");
	    return -1;
	  };
	  vmin = *FLTPTR(invals[0]);
	  vmax = *FLTPTR(invals[1]);
	  arg=2;
	}
/*
  One array - find min and max values in array.
*/
	else {
	  inptr=FLTPTR(invals[0]);
	  vmin=vmax=*inptr;
	  for(i=0; i<npts; i++) {
/*
  Update the min and max for the x-array.
*/
	    if(vmin > *inptr)
	      vmin = *inptr;
	    else if(vmax < *inptr)
	      vmax = *inptr;
	    inptr++;
	  };
	  arg=1;
	};
/*
  Add a margin to the x-limits equal to the fraction (x_gap) of the
  total x-range.
*/
	limits[0] = vmin - (vmax-vmin) * x_gap;
	limits[1] = vmax + (vmax-vmin) * x_gap;
/*
  Now the y-limits.
*/
	if(npar < arg+1) {
	  lprintf(stderr, "No y-limits given\n");
	  return -1;
	};
	npts = invals[arg]->adim[0];
/*
  Two scalars?
*/
	if(npar > arg+1 && npts == 1) {
	  if(invals[arg+1]->adim[0] != 1) {
	    lprintf(stderr, "First x-limit scalar but not second.\n");
	    return -1;
	  };
	  vmin = *FLTPTR(invals[arg]);
	  vmax = *FLTPTR(invals[arg+1]);
	  arg += 2;
	}
/*
  One array - find min and max values in array.
*/
	else {
	  inptr=FLTPTR(invals[arg]);
	  vmin = vmax = *inptr;
	  for(i=0; i<npts; i++) {
/*
  Update the min and max for the x-array.
*/
	    if(vmin > *inptr)
	      vmin = *inptr;
	    else if(vmax < *inptr)
	      vmax = *inptr;
	    inptr++;
	  };
	  arg++;
	};
/*
  Add a margin to the y-limits equal to the fraction (y_gap) of the
  total y-range.
*/
	limits[2] = vmin - (vmax-vmin) * y_gap;
	limits[3] = vmax + (vmax-vmin) * y_gap;
/*
  Check the limits.
*/
	if(limits[0] == limits[1] || limits[2] == limits[3]) {
	  lprintf(stderr, "Illegal limits: %f,%f,%f,%f\n",
		  limits[0],limits[1],limits[2],limits[3]);
	  return -1;
	};
/*
  Return with a record of the number of user arguments used.
*/
	return arg;
}

/*.......................................................................
  Open a new PGPLOT device and determine some of the device
  characteristics, storing them in user parameters for access both by the
  user and the programmer. If the device name sent is NULL or "?" then
  list available pgplot devices and tell pgplot to prompt for a device
  name. xnum and ynum specify the number of sub-divisions of the plot
  surface. On error -1 is returned otherwise 0.
*/
static int pgplot_newdev(char *name, int xnum, int ynum)
{
        char *dev_name;  /* Pointer to device name string */
	char answer[10]; /* String for answer from pgplot inquiry routine */
	int slen;        /* Length of strig returned in answer[] */
/*
  Check that the xnum and ynum are valid and enforce defaults.
*/
	xnum = (xnum==0) ? 1:xnum;
	ynum = (ynum==0) ? 1:ynum;
/*
  Get a pointer to the device name required. This is "?" if no name
  was supplied, otherwise it is the string pointed by name.
  If name is NULL or contains "?" set up
  a pointer to a string consisting of "?" to get pgplot
  to prompt for the device name.
*/
	if(name == NULL || name[0] == '?') {
	  dev_name = "?";
	} else {
	  dev_name = name;
	  lprintf(stdout, "Attempting to open device: '%s'\n", name);
	};
/*
  Attempt to open the device.
*/
	if(cpgbeg(0, dev_name, xnum, ynum) != 1) {
/*
  No device is open - record this fact for examination by
  functions that require one to be and for examination by the user.
*/
	  plot_open = has_cursor = 0;
	  return -1;
	} else {
/*
  Flag success in the user parameter alias plot_open.
*/
	  plot_open = 1;
/*
  Ask PGPLOT if the device has a cursor.
*/
	  slen = 4;
	  cpgqinf("CURSOR", answer, &slen);
	  has_cursor = strncmp(answer,"YES",slen) == 0;
/*
 * Turn off PGPLOT prompting for next page.
 */
	  cpgask(0);
	};
/*
  Set up the standard view port.
*/
	cpgvstd();
	printf("\n");
	return no_error;
}

/*.......................................................................
  If a pgplot device is open then this routine returns immediately.
  Otherwise the user is prompted for a device name and an attempt is made
  to open it. The function returns -1 if the open fails.
*/
int make_open(void)
{
        if(plot_open) return no_error;
	return pgplot_newdev(NULL,1,1);
}

/*.......................................................................
  If a pgplot device has been opened then this function returns
  immediately. Otherwise an error is signalled by returning -1.
*/
int check_open(void)
{
        if(plot_open) return no_error;
	lprintf(stderr, "No PGPLOT device active\n");
	return -1;
}

/*.......................................................................
  If the current pgplot device has a cursor then this function returns
  immediately. Otherwise an error is signalled by returning -1.
*/
static int check_cursor(void)
{
        if(has_cursor) return no_error;
	lprintf(stderr, "The current plot device has no cursor\n");
	return -1;
}

/*.......................................................................
 * Set the color representation of a given color index using RGB levels.
 *
 * Input:
 *  ci   int   Color index who's representation is to be changed.
 *  cr float   0 -> 1  The fraction of saturation red color to use.
 *  cg float   0 -> 1  The fraction of saturation green color to use.
 *  cb float   0 -> 1  The fraction of saturation blue color to use.
 */
static Template(pgscr_fn)
{
  int ci = *INTPTR(invals[0]);
  float cr = *FLTPTR(invals[1]);
  float cg = *FLTPTR(invals[2]);
  float cb = *FLTPTR(invals[3]);
/*
 * Device open?
 */
  if(make_open() == -1)
    return -1;
/*
 * Check inputs.
 */
  if(ci < 0) {
    lprintf(stderr, "pgscr: Color index \'%d\' out of range\n", ci);
    return -1;
  };
  if(cr<0.0f || cr > 1.0f || cg<0.0f || cg > 1.0f || cb<0.0f || cb > 1.0f) {
    lprintf(stderr, "pgscr: Color representation %g,%g,%g out of range\n",
	    cr,cg,cb);
    return -1;
  };
/*
 * Set the color representation.
 */
  cpgscr(ci,cr,cg,cb);
  return no_error;
}

/*.......................................................................
 * Set the color representation of a given color index using
 * Hue-Lightness-Saturation levels.
 *
 * Input:
 *  ci   int   Color index who's representation is to be changed.
 *  ch float   0 -> 1  The hue fraction.
 *  cl float   0 -> 1  The lightness fraction.
 *  cs float   0 -> 1  The saturation fraction.
 */
static Template(pgshls_fn)
{
  int ci = *INTPTR(invals[0]);
  float ch = *FLTPTR(invals[1]);
  float cl = *FLTPTR(invals[2]);
  float cs = *FLTPTR(invals[3]);
/*
 * Device open?
 */
  if(make_open() == -1)
    return -1;
/*
 * Check inputs.
 */
  if(ci < 0) {
    lprintf(stderr, "pgshls: Color index \'%d\' out of range\n", ci);
    return -1;
  };
  if(ch<0.0f || ch > 360.0f || cl<0.0f || cl > 1.0f || cs<0.0f || cs > 1.0f) {
    lprintf(stderr, "pgshls: Color representation %g,%g,%g out of range\n",
	    ch,cl,cs);
    return -1;
  };
/*
 * Set the color representation.
 */
  cpgshls(ci,ch,cl,cs);
  return no_error;
}

/*.......................................................................
 * Set the PGPLOT character font.
 *
 * Input:
 *  font   int   The enumerated ID of the selected PGPLOT font.
 */
static Template(pgscf_fn)
{
  int font = *INTPTR(invals[0]);
/*
 * Device open?
 */
  if(make_open() == -1)
    return -1;
/*
 * Check inputs.
 */
  if(font < 1 || font > 4) {
    lprintf(stderr, "pgscf: Unknown font ID (%d).\n", font);
    return -1;
  };
/*
 * Set the font.
 */
  cpgscf(font);
  return no_error;
}

/*.......................................................................
 * Save the current PGPLOT attributes for later retrieval with pgunsa.
 */
static Template(pgsave_fn)
{
/*
 * Device open?
 */
  if(make_open() == -1)
    return -1;
  cpgsave();
  return no_error;
}

/*.......................................................................
 * Restore PGPLOT attributes from a previous call to pgsave.
 */
static Template(pgunsa_fn)
{
/*
 * Device open?
 */
  if(make_open() == -1)
    return -1;
  cpgunsa();
  return no_error;
}

/*.......................................................................
 * Draw a line with an arrow at one end. This is a wrapper around
 * the PGPLOT pgarro() subroutine.
 */
static Template(pgarro_fn)
{
  if(check_open() == -1)
    return -1;
  cpgarro(*FLTPTR(invals[0]), *FLTPTR(invals[1]), *FLTPTR(invals[2]),
	  *FLTPTR(invals[3]));
  return no_error;
}

/*.......................................................................
 * Toggle new-page prompting. This is a wrapper around the PGPLOT pgask()
 * subroutine.
 */
static Template(pgask_fn)
{
  if(check_open() == -1)
    return -1;
  cpgask(*LOGPTR(invals[0])!=0);
  return no_error;
}

/*.......................................................................
 * Display the PGPLOT cursor and read back its position and optionally the
 * key pressed when the user selects a point with it. This is a wrapper
 * around the PGPLOT pgband() function.
 */
static Template(pgband_fn)
{
  int mode;  /* The type of cursor to display */
  int posn;  /* Whether to position the cursor when it is first displayed */
  float xref,yref; /* The reference point of the cursor, where relevant */
  float x,y; /* The position that the user selected */
  char ch;   /* The character that the user selected */
  int ok;    /* True if the call sucedes */
/*
 * Make sure that PGPLOT is open.
 */
  if(check_open() == -1)
    return -1;
/*
 * Get the function arguments.
 */
  mode = *INTPTR(invals[0]);
  posn = *INTPTR(invals[1]);
  xref = *FLTPTR(invals[2]);
  yref = *FLTPTR(invals[3]);
  x = *FLTPTR(invals[4]);
  y = *FLTPTR(invals[5]);
  ch = '\0';
/*
 * Call pgband.
 */
  ok = cpgband(mode, posn, xref, yref, &x, &y, &ch);
/*
 * Set up the function return value if this is being used
 * as a function. If not then translate the pgband error return
 * into an exception.
 */
  if(outvals)
    *INTPTR(outvals) = ok;
  else if(ok != 1) {
    lprintf(stderr, "cpgband: PGPLOT cpgband() returned an error.\n");
    return -1;
  };
/*
 * Set up the return arguments.
 */
  if(ok) {
    *FLTPTR(invals[4]) = x;
    *FLTPTR(invals[5]) = y;
/*
 * Only allocate a new string if necessary.
 */
    if(strlen(*STRPTR(invals[6])) < 1) {
      char *ptr = stralloc(1);
      if(!ptr)
	return -1;
      free(*STRPTR(invals[6]));
      *STRPTR(invals[6]) = ptr;
    };
    (*STRPTR(invals[6]))[0] = ch;
    (*STRPTR(invals[6]))[1] = '\0';
  };
  return no_error;
}

/*.......................................................................
 * Draw a circle. This is a wrapper around the PGPLOT pgcirc() subroutine.
 */
static Template(pgcirc_fn)
{
  if(check_open() == -1)
    return -1;
  cpgcirc(*FLTPTR(invals[0]), *FLTPTR(invals[1]), *FLTPTR(invals[2]));
  return no_error;
}

/*.......................................................................
 * This is a wrapper around the pgcont() contour plotting function.
 */
static Template(pgcont_fn)
{
  float *array;   /* The 2D array to be contoured */
  int idim,jdim;  /* The dimensions of the 2D array */
  int i1,i2;      /* The ranges of the first index of the array */
  int j1,j2;      /* The ranges of the second index of the array */
  float *levels;  /* The array of contour levels */
  int nlev;       /* The number of contour levels */
  float *tr;      /* The transformation matrix */
/*
 * Make sure that we have a pgplot device open.
 */
  if(check_open() == -1)
    return -1;
/*
 * Get the 2D array and its dimensions.
 */
  array = FLTPTR(invals[0]);
  idim=invals[0]->adim[0];
  jdim=invals[0]->adim[1];
/*
 * Get the ranges of indexes which delimit the sub-array to be contoured.
 */
  i1 = *INTPTR(invals[1]);
  i2 = *INTPTR(invals[2]);
  j1 = *INTPTR(invals[3]);
  j2 = *INTPTR(invals[4]);
/*
 * Check that they make sense.
 */
  if(i1 < 1 || i1 > idim || i2 < 1 || i2 > idim) {
    lprintf(stderr, "pgcont: i indexes %d-%d out of range %d-%d\n", i1, i2,
	    1, idim);
    return -1;
  };
  if(j1 < 1 || j1 > jdim || j2 < 1 || j2 > jdim) {
    lprintf(stderr, "pgcont: j indexes %d-%d out of range %d-%d\n", j1, j2,
	    1, jdim);
    return -1;
  };
/*
 * Get the array of contour levels and their number.
 */
  levels = FLTPTR(invals[5]);
  nlev = invals[5]->adim[0];
/*
 * Get the transformation matrix, and check that it has the required
 * dimensions.
 */
  tr = FLTPTR(invals[6]);
  if(invals[6]->adim[0] < 6) {
    lprintf(stderr,
	 "pgcont: The tr argument has less than the necessary 6 elements.\n");
    return -1;
  };
/*
 * Hand over to pgplot.
 */
  cpgcont(array, idim, jdim, i1, i2, j1, j2, levels, nlev, tr);
  return no_error;
}

/*.......................................................................
 * Display the PGPLOT cursor and read back its position and optionally the
 * key pressed when the user selects a point with it. This is a wrapper
 * around the PGPLOT pgcurs() function.
 */
static Template(pgcurs_fn)
{
  float x,y; /* The position that the user selected */
  char ch;   /* The character that the user selected */
  int ok;    /* True if the call sucedes */
/*
 * Make sure that PGPLOT is open.
 */
  if(check_open() == -1)
    return -1;
/*
 * Get the function arguments.
 */
  x = *FLTPTR(invals[0]);
  y = *FLTPTR(invals[1]);
  ch = '\0';
/*
 * Call pgcurs.
 */
  ok = cpgcurs(&x, &y, &ch);
/*
 * Set up the function return value if this is being used
 * as a function. If not then translate the pgcurs error return
 * into an exception.
 */
  if(outvals)
    *INTPTR(outvals) = ok;
  else if(ok != 1) {
    lprintf(stderr, "cpgcurs: PGPLOT cpgcurs() returned an error.\n");
    return -1;
  };
/*
 * Set up the return arguments.
 */
  if(ok) {
    *FLTPTR(invals[0]) = x;
    *FLTPTR(invals[1]) = y;
/*
 * Only allocate a new string if necessary.
 */
    if(strlen(*STRPTR(invals[2])) < 1) {
      char *ptr = stralloc(1);
      if(!ptr)
	return -1;
      free(*STRPTR(invals[2]));
      *STRPTR(invals[2]) = ptr;
    };
    (*STRPTR(invals[2]))[0] = ch;
    (*STRPTR(invals[2]))[1] = '\0';
  };
  return no_error;
}

/*.......................................................................
 * Close all PGPLOT devices. This is a wrapper function for the PGPLOT
 * pgend() subroutine.
 */
static Template(pgend_fn)
{
  cpgend();
  plot_open = 0;
  return no_error;
}

/*.......................................................................
 * Draw a line with an arrow at one end. This is a wrapper around
 * the PGPLOT pgenv() subroutine.
 */
static Template(pgenv_fn)
{
  if(check_open() == -1)
    return -1;
  cpgenv(*FLTPTR(invals[0]), *FLTPTR(invals[1]), *FLTPTR(invals[2]),
	 *FLTPTR(invals[3]), *INTPTR(invals[4]),*INTPTR(invals[5]));
  return no_error;
}

/*.......................................................................
 * Erase all PGPLOT graphics. This is a wrapper for the PGPLOT
 * pgeras() subroutine.
 */
static Template(pgeras_fn)
{
  if(make_open() == -1)
    return -1;
  cpgeras();
  return no_error;
}

/*.......................................................................
 * Plot one or more vertical error bars horizontally or vertically.
 *
 * Input:
 *  dir     int    The direction of the error-bar.
 *  x     float *  The X-axis coordinate of each error bar.
 *  y1    float *  The lower Y coordinate of each error bar.
 *  y2    float *  The upper Y coordinate of each error bar.
 *  size  float    The size of the error bar terminal wrt the default.
 */
static Template(pgerrb_fn)
{
  int nvals;         /* The number of values to be plotted */
  float size = 1.0f; /* The normalized size of the error bar terminals */
  int i;
/*
 * Device open?
 */
  if(check_open() == -1)
    return -1;
/*
 * Find the minimum of the sizes of each of the provided arrays.
 */
  nvals = invals[1]->adim[0];
  for(i=2; i<4; i++) {
    if(invals[i]->adim[0] < nvals)
      nvals = invals[i]->adim[0];
  };
/*
 * See if a error-bar terminal size was provided.
 */
  if(npar > 4)
    size = *FLTPTR(invals[4]);
/*
 * Pass the results to PGPLOT.
 */
  cpgerrb(*INTPTR(invals[0]), nvals, FLTPTR(invals[1]), FLTPTR(invals[2]),
	  FLTPTR(invals[3]), size);
  return no_error;
}

/*.......................................................................
 * Plot one or more horizontal error bars.
 *
 * Input:
 *  x1    float *  The leftmost X coordinate of each error bar.
 *  x2    float *  The rightmost X coordinate of each error bar.
 *  y     float *  The Y-axis coordinate of each error bar.
 *  size  float    The size of the error bar terminal wrt the default.
 */
static Template(pgerrx_fn)
{
  int nvals;         /* The number of values to be plotted */
  float size = 1.0f; /* The normalized size of the error bar terminals */
  int i;
/*
 * Device open?
 */
  if(check_open() == -1)
    return -1;
/*
 * Find the minimum of the sizes of each of the provided arrays.
 */
  nvals = invals[0]->adim[0];
  for(i=1; i<3; i++) {
    if(invals[i]->adim[0] < nvals)
      nvals = invals[i]->adim[0];
  };
/*
 * See if a error-bar terminal size was provided.
 */
  if(npar > 3)
    size = *FLTPTR(invals[3]);
/*
 * Pass the results to PGPLOT.
 */
  cpgerrx(nvals, FLTPTR(invals[0]), FLTPTR(invals[1]), FLTPTR(invals[2]),
	  size);
  return no_error;
}


/*.......................................................................
 * Plot one or more vertical error bars.
 *
 * Input:
 *  x     float *  The X-axis coordinate of each error bar.
 *  y1    float *  The lower Y coordinate of each error bar.
 *  y2    float *  The upper Y coordinate of each error bar.
 *  size  float    The size of the error bar terminal wrt the default.
 */
static Template(pgerry_fn)
{
  int nvals;         /* The number of values to be plotted */
  float size = 1.0f; /* The normalized size of the error bar terminals */
  int i;
/*
 * Device open?
 */
  if(check_open() == -1)
    return -1;
/*
 * Find the minimum of the sizes of each of the provided arrays.
 */
  nvals = invals[0]->adim[0];
  for(i=1; i<3; i++) {
    if(invals[i]->adim[0] < nvals)
      nvals = invals[i]->adim[0];
  };
/*
 * See if a error-bar terminal size was provided.
 */
  if(npar > 3)
    size = *FLTPTR(invals[3]);
/*
 * Pass the results to PGPLOT.
 */
  cpgerry(nvals, FLTPTR(invals[0]), FLTPTR(invals[1]), FLTPTR(invals[2]),
	  size);
  return no_error;
}

/*.......................................................................
 * This is a wrapper around the pggray() gray-scale imaging function.
 */
static Template(pggray_fn)
{
  float *array;   /* The 2D array to be displayed */
  int idim,jdim;  /* The dimensions of the 2D array */
  int i1,i2;      /* The ranges of the first index of the array */
  int j1,j2;      /* The ranges of the second index of the array */
  float bg,fg;    /* The background and foreground levels of the array */
  float *tr;      /* The transformation matrix */
/*
 * Make sure that we have a pgplot device open.
 */
  if(check_open() == -1)
    return -1;
/*
 * Get the 2D array and its dimensions.
 */
  array = FLTPTR(invals[0]);
  idim=invals[0]->adim[0];
  jdim=invals[0]->adim[1];
/*
 * Get the ranges of indexes which delimit the sub-array to be gray-scaled.
 */
  i1 = *INTPTR(invals[1]);
  i2 = *INTPTR(invals[2]);
  j1 = *INTPTR(invals[3]);
  j2 = *INTPTR(invals[4]);
/*
 * Check that they make sense.
 */
  if(i1 < 1 || i1 > idim || i2 < 1 || i2 > idim) {
    lprintf(stderr, "pggray: i indexes %d-%d out of range %d-%d\n", i1, i2,
	    1, idim);
    return -1;
  };
  if(j1 < 1 || j1 > jdim || j2 < 1 || j2 > jdim) {
    lprintf(stderr, "pggray: j indexes %d-%d out of range %d-%d\n", j1, j2,
	    1, jdim);
    return -1;
  };
/*
 * Get the foreground and background levels.
 */
  fg = *FLTPTR(invals[5]);
  bg = *FLTPTR(invals[6]);
/*
 * Get the transformation matrix, and check that it has the required
 * dimensions.
 */
  tr = FLTPTR(invals[7]);
  if(invals[7]->adim[0] < 6) {
    lprintf(stderr,
	 "pggray: The tr argument has less than the necessary 6 elements.\n");
    return -1;
  };
/*
 * Hand over to pgplot.
 */
  cpggray(array, idim, jdim, i1, i2, j1, j2, fg, bg, tr);
  return no_error;
}

/*.......................................................................
 * Set the color representation of a PGPLOT color index by name.
 * This is a wrapper function for the PGPLOT pgscir() subroutine.
 */
static Template(pgscir_fn)
{
  cpgscir(*INTPTR(invals[0]), *INTPTR(invals[1]));
  return no_error;
}

/*.......................................................................
 * Set the color representation of a PGPLOT color index by name.
 * This is a wrapper function for the PGPLOT pgscrn() subroutine.
 */
static Template(pgscrn_fn)
{
  int ier;   /* The error-status returned by PGPLOT */
  cpgscrn(*INTPTR(invals[0]), *STRPTR(invals[1]), &ier);
  if(npar > 2)
    *INTPTR(invals[2]) = ier;
  else if(ier == 1)
    return -1;
  return no_error;
}

/*.......................................................................
 * Set the PGPLOT fill-area style for subsequent graphics.
 *
 * Input:
 *  font   int   The enumerated ID of the style.
 */
static Template(pgsfs_fn)
{
  if(make_open() == -1)
    return -1;
  cpgsfs(*INTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Set the PGPLOT hatching style for use with fill-area style 3. This
 * is a wrapper function around the PGPLOT pgshs() subroutine.
 */
static Template(pgshs_fn)
{
  if(make_open() == -1)
    return -1;
  cpgshs(*FLTPTR(invals[0]), *FLTPTR(invals[1]), *FLTPTR(invals[2]));
  return no_error;
}


/*.......................................................................
 * Set the current world coordinates. This is a wrapper function around
 * the PGPLOT pgswin() subroutine.
 */
static Template(pgswin_fn)
{
  cpgswin(*FLTPTR(invals[0]), *FLTPTR(invals[1]),
	  *FLTPTR(invals[2]), *FLTPTR(invals[3]));
  return no_error;
}

/*.......................................................................
 * Set the current viewport. This is a wrapper function around
 * the PGPLOT pgsvp() subroutine.
 */
static Template(pgsvp_fn)
{
  cpgsvp(*FLTPTR(invals[0]), *FLTPTR(invals[1]),
	 *FLTPTR(invals[2]), *FLTPTR(invals[3]));
  return no_error;
}

/*.......................................................................
 * Set the standard viewport. This is a wrapper function around
 * the PGPLOT pgsvp() subroutine.
 */
static Template(pgvstd_fn)
{
  cpgvstd();
  return no_error;
}

/*.......................................................................
 * Start a new page on the current PGPLOT device. This is a wrapper for
 * the PGPLOT pgpage() subroutine.
 */
static Template(pgpage_fn)
{
  if(make_open() == -1)
    return -1;
  cpgpage();
  return no_error;
}

/*.......................................................................
 * Set the current world coordinates and adjust the viewport to the
 * same aspect. This is a wrapper function around the PGPLOT pgwnad()
 * subroutine.
 */
static Template(pgwnad_fn)
{
  cpgwnad(*FLTPTR(invals[0]), *FLTPTR(invals[1]),
	  *FLTPTR(invals[2]), *FLTPTR(invals[3]));
  return no_error;
}

/*.......................................................................
 * Shade a polygonal area. This is a wrapper around the PGPLOT pgpoly()
 * subroutine.
 */
static Template(pgpoly_fn)
{
  int nvals;
  if(check_open() == -1)
    return -1;
/*
 * Determine the size of the vertex arrays.
 */
  nvals = invals[0]->adim[0];
  if(invals[1]->adim[0] != nvals) {
    lprintf(stderr, "pgpoly: The X and Y arrays differ in length\n");
    return -1;
  };
  cpgpoly(nvals, FLTPTR(invals[0]), FLTPTR(invals[1]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqah().
 */
static Template(pgqah_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqah(INTPTR(invals[0]), FLTPTR(invals[1]), FLTPTR(invals[2]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqcf().
 */
static Template(pgqcf_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqcf(INTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqch().
 */
static Template(pgqch_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqch(FLTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqci().
 */
static Template(pgqci_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqci(INTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqcir().
 */
static Template(pgqcir_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqcir(INTPTR(invals[0]), INTPTR(invals[1]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqcol().
 */
static Template(pgqcol_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqcol(INTPTR(invals[0]), INTPTR(invals[1]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqcr().
 */
static Template(pgqcr_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqcr(*INTPTR(invals[0]), FLTPTR(invals[1]), FLTPTR(invals[2]),
	 FLTPTR(invals[3]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqcs().
 */
static Template(pgqcs_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqcs(*INTPTR(invals[0]), FLTPTR(invals[1]), FLTPTR(invals[2]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqfs().
 */
static Template(pgqfs_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqfs(INTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqhs().
 */
static Template(pgqhs_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqhs(FLTPTR(invals[0]), FLTPTR(invals[1]), FLTPTR(invals[2]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqid().
 */
static Template(pgqid_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqid(INTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqitf().
 */
static Template(pgqitf_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqitf(INTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqls().
 */
static Template(pgqls_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqls(INTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqlw().
 */
static Template(pgqlw_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqlw(INTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqpos().
 */
static Template(pgqpos_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqpos(FLTPTR(invals[0]), FLTPTR(invals[1]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgvp().
 */
static Template(pgqvp_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqvp(*INTPTR(invals[0]),
	 FLTPTR(invals[1]), FLTPTR(invals[2]),
	 FLTPTR(invals[3]), FLTPTR(invals[4]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgvsz().
 */
static Template(pgqvsz_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqvsz(*INTPTR(invals[0]),
	  FLTPTR(invals[1]), FLTPTR(invals[2]),
	  FLTPTR(invals[3]), FLTPTR(invals[4]));
  return no_error;
}

/*.......................................................................
 * Provide a wrapper around pgqwin().
 */
static Template(pgqwin_fn)
{
  if(check_open() == -1)
    return -1;
  cpgqwin(FLTPTR(invals[0]), FLTPTR(invals[1]),
	  FLTPTR(invals[2]), FLTPTR(invals[3]));
  return no_error;
}



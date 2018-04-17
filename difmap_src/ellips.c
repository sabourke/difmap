#include <math.h>
#include <stdio.h>

#include "vlbconst.h"
#include "logio.h"
#include "cpgplot.h"

#include "ellips.h"

static int el_bad(Ellipse *el, char *fname);

/*.......................................................................
 * Fill an Ellipse descriptor, via the dimensions of the required ellipse.
 * This descriptor may then be used in transactions with other ellipse
 * method functions.
 *
 * Input:
 *  el    Ellipse * The pointer to the Ellipse descriptor to be filled.
 *  minor   float   The minor axis diameter (world coordinates).
 *  major   float   The major axis diameter (world coordinates).
 *  pa      float   The clockwise angle which the major axis subtends to
 *                  the +ve Y-axis (radians).
 *  xc      float   The X coordinate of the ellipse center.
 *  yc      float   The Y coordinate of the ellipse center.
 *                  This may be changed later via el_move().
 */
void el_define(Ellipse *el, float minor, float major, float pa,
	       float xc, float yc)
{
  float xang;  /* The position angle of the locus with the largest X position */
  float yang;  /* The position angle of the locus with the largest Y position */
/*
 * Invalid descriptor?
 */
  if(el_bad(el, "el_define"))
    return;
/*
 * Enforce positivity on the axis extents.
 */
  if(minor < 0.0)
    minor = -minor;
  if(major < 0.0)
    major = -major;
/*
 * Record the ellipse chracteristics.
 */
  if(major >= minor) {
    el->minor = minor;
    el->major = major;
    el->pa = pa;
  } else {              /* Swap major and minor if major < minor */
    el->minor = major;
    el->major = minor;
    el->pa = pa - pi/2; /* Adjust the pa for the new major axis */
  };
/*
 * Record the center position of the ellipse.
 */
  el->xc = xc;
  el->yc = yc;
/*
 * Calculate the maximum x-extent of the ellipse.  By differentiating the
 * equation for X, vs clockwise angle and setting the result to zero,
 * one finds that this occurs at angle xang = atan(minor/major/tan(pa).
 * This equation produces overflows near 0 and halfpi radians, so work
 * out these special cases first.
 */
  if(fabs(el->pa) < 0.01 || el->major==0.0)
    xang = halfpi;
  else if(fabs(fabs(el->pa)-halfpi) < 0.01)
    xang = 0.0f;
  else
    xang = atan((1.0f/tan(el->pa)) * el->minor/el->major);
/*
 * Determine the X-axis extent from the X-position on the ellipse
 * at xang.
 */
  {
    float x;
    el_locus(el, xang, &x, NULL);
    el->xwid = 2.0 * fabs(x - el->xc);
  };
/*
 * The max Y-extent is at angle yang = -atan(minor/major*tan(pa)).
 */
  if(fabs(fabs(el->pa)-halfpi) < 0.01 || el->major==0.0)
    yang = -halfpi;
  else
    yang = -atan(tan(el->pa)*el->minor/el->major);
/*
 * Determine the Y-axis extent from the Y-position on the ellipse
 * at yang.
 */
  {
    float y;
    el_locus(el, yang, NULL, &y);
    el->ywid = 2.0 * fabs(y - el->yc);
  };
  return;
}

/*.......................................................................
 * Define a new center for a given ellipse.
 *
 * Input:
 *  el   Ellipse *  An ellipse descriptor previously initialized via
 *                  el_define().
 *  xc     float    The X-coordinate of the new centroid (world coordinates).
 *  yc     float    The Y-coordinate of the new centroid (world coordinates).
 */
void el_move(Ellipse *el, float xc, float yc)
{
  if(el_bad(el, "el_move"))
    return;
  el->xc = xc;
  el->yc = yc;
  return;
}

/*.......................................................................
 * Return the X and Y axis positions that correspond to the locus of
 * a given ellipse at a given position angle.
 *
 * Input:
 *  el   Ellipse *  An ellipse descriptor previously initialized via
 *                  el_define().
 *  theta  float    The clockwise polar position angle wrt the +ve Y-axis,
 *                  defining the required position on the ellipse (radians).
 * Input/Output:
 *  x      float *  If x!=NULL, *x will be assigned the X-axis position of
 *                  the locus of the ellipse at polar angle 'theta'.
 *  y      float *  If y!=NULL, *y will be assigned the Y-axis position of
 *                  the locus of the ellipse at polar angle 'theta'.
 */
void el_locus(Ellipse *el, float theta, float *x, float *y)
{
  float sinpa; /* sin(pa) */
  float cospa; /* cos(pa) */
  float minax; /* The distance of the point along the minor axis */
  float majax; /* The distance of the point along the major axis */
/*
 * Calculate the unrotated minor and major axis coordinates of the point.
 */
  minax = el->minor * sin(theta)/2.0;
  majax = el->major * cos(theta)/2.0;
/*
 * Pre-calculate sin(pa) and cos(pa) since they are each used twice.
 */
  sinpa = sin(el->pa);
  cospa = cos(el->pa);
/*
 * Return the position, rotated clockwise by 'pa' radians, and offset
 * by el->xc,el->yc.
 */
  if(x) *x = el->xc + minax * cospa + majax * sinpa;
  if(y) *y = el->yc + majax * cospa - minax * sinpa;
  return;
}

/*.......................................................................
 * Report whether a given ellipse lies partially or fully within the
 * confines of a specified rectangular area.
 *
 * Input:
 *  el   Ellipse *  An ellipse descriptor previously initialized via
 *                  el_define().
 *  xa     float    The X-axis location of one edge of the area (world coords).
 *  xb     float    The X-axis location of the opposite edge of the area to xa.
 *  ya     float    The Y-axis location of one edge of the area (world coords).
 *  yb     float    The Y-axis location of the opposite edge of the area to ya.
 *  state Elstat    Specify what is considered to be sufficiently visible.
 *                   EL_FULL - Fully visible.
 *                   EL_PART - Partially visible.
 *                   EL_CENT - Center visible.
 * Output:
 *  return   int    0 - The ellipse is visible.
 *                  1 - The ellipse is not visible.
 */
int el_visible(Ellipse *el, float xa, float xb, float ya, float yb,
	       Elstat state)
{
  float exa,exb;  /* The min,max X coordinates of the ellipse limits */
  float eya,eyb;  /* The min,max Y coordinates of the ellipse limits */
/*
 * Swap xa and xb such that xa < xb.
 */
  if(xa > xb) {
    float x = xa;
    xa = xb;
    xb = x;
  };
/*
 * Swap ya and yb such that ya < yb.
 */
  if(ya > yb) {
    float y = ya;
    ya = yb;
    yb = y;
  };
/*
 * Determine the limits of the beam wrt the current ellipse center.
 */
  exa = el->xc - el->xwid/2.0;
  exb = exa + el->xwid;
  eya = el->yc - el->ywid/2.0;
  eyb = eya + el->ywid;
/*
 * Determine the requested visibility status.
 */
  switch(state) {
  case EL_FULL:
    return (exa >= xa && exb <= xb) && (eya >= ya && eyb <= yb);
    break;
  case EL_PART:
    return !((exb < xa || exa > xb) || (eyb < ya || eya > yb));
    break;
  case EL_CENT:
    return (el->xc >= xa && el->xc <= xb) && (el->xc >= xa && el->xc <= xb);
    break;
  default:
    lprintf(stderr, "el_visible: Unknown visibility designation.\n");
    return 0;
  };
}

/*.......................................................................
 * Plot a given ellipse, optionally filled and/or outlined.
 *
 * Input:
 *  el   Ellipse *  An ellipse descriptor previously initialized via
 *                  el_define().
 *  outline  int    The color index to use to outline the ellipse,
 *                  or -1 if no outline is wanted.
 *  fill     int    The color index to use to fill the area inside the
 *                  ellipse, or -1 if no filling is wanted.
 *  cross    int    Line style to plot major and minor axes with, or 0
 *                  if the axes are not required.
 *  nmax     int    The maximum number of points to be used to draw the
 *                  ellipse (an internal maximum is also set).
 */
void el_plot(Ellipse *el, int outline, int fill, int cross, int nmax)
{
  enum {ELMAX=50};
  static float xp[ELMAX];   /* The array for the array of X coordinates */
  static float yp[ELMAX];   /* The array for the array of Y coordinates */
  float step;               /* The angle step size between points (radians) */
  int oldcol, oldfil, oldls;/* PGPLOT attributes to be restored on return */
  int i;
/*
 * Unuseable descriptor?
 */
  if(el_bad(el, "el_plot"))
    return;
/*
 * Record the current color index, fill-style and line style attributes.
 */
  cpgqci(&oldcol);
  cpgqfs(&oldfil);
  cpgqls(&oldls);
/*
 * Limit the requested maximum number of points to be used to draw the
 * ellipse.
 */
  if(nmax <= 2 || nmax > ELMAX)
    nmax = ELMAX;
/*
 * Determine the angle step increment required to plot a complete ellipse
 * with 'nmax' points.
 */
  step = twopi / (nmax-1);
/*
 * Get the X and Y coordinates of the ellipse loci.
 */
  for(i=0; i<nmax; i++)
    el_locus(el, i*step, &xp[i], &yp[i]);
/*
 * Solid line for ellipse.
 */
  cpgsls(1);
/*
 * Fill the ellipse if requested.
 */
  if(fill >= 0) {
    cpgsci(fill);
    cpgsfs(1);
    cpgpoly(nmax, xp, yp);
  };
/*
 * Outline the ellipse if requested.
 */
  if(outline >= 0) {
    cpgsci(outline);
    cpgline(nmax, xp, yp);
  };
/*
 * Draw axis lines if required.
 */
  if(cross > 0 && cross < 6) {
    float xoff,yoff;
    cpgsls(cross);
/*
 * First the major axis.
 */
    xoff = 0.5 * el->major * sin(el->pa);
    yoff = 0.5 * el->major * cos(el->pa);
    cpgmove(el->xc - xoff, el->yc - yoff);
    cpgdraw(el->xc + xoff, el->yc + yoff);
/*
 * Then the minor axis.
 */
    xoff =  0.5 * el->minor * cos(el->pa);
    yoff = -0.5 * el->minor * sin(el->pa);
    cpgmove(el->xc - xoff, el->yc - yoff);
    cpgdraw(el->xc + xoff, el->yc + yoff);
  };
/*
 * Restore the color index, fill-area and line styles to the input values.
 */
  cpgsci(oldcol);
  cpgsfs(oldfil);
  cpgsls(oldls);
  return;
}

/*.......................................................................
 * Return true if the descriptor is invalid.
 *
 * Input:
 *  el    Ellipse *  The descriptor to be checked.
 *  fname    char *  The name of the calling function.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Descriptor invalid.
 */
static int el_bad(Ellipse *el, char *fname)
{
  if(el==NULL) {
    lprintf(stderr, "%s: NULL Ellipse descriptor received.\n",
	    fname ? fname : "(unknown)");
    return 1;
  };
  return 0;
}

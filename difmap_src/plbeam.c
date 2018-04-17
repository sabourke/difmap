#include <math.h>
#include <stdio.h>

#include "logio.h"
#include "vlbmath.h"
#include "vlbconst.h"
#include "ellips.h"
#include "cpgplot.h"

#define NBEAM 30    /* Number of vectors in beam ellipse */

/*.......................................................................
 * Draw a restoring beam as an ellipse at a given position within the
 * current plot. The ellipse is shaded by drawing many parallel lines
 * within the ellipse, perpendicular to the major axis.
 * If either 'xpos' or 'ypos' are not between 0 and 1, inclusive, then no
 * beam will be drawn. Otherwise the beam will be positioned at the
 * position given, with the exception that if the position would place
 * part of the beam outside the map it will be moved inwards, suitably
 * spaced from the edge. If either bmin or bmaj are <= 0 then no beam
 * will be plotted.
 *
 * NB. The map must be plotted on the current PGPLOT device before
 * calling this function.
 *
 * Inputs:
 *   bmin   float    Minor axis length of beam, in world coordinate units.
 *   bmaj   float    Major axis length of beam, in world coordinate units.
 *   bpa    float    Beam position angle in radians.
 *   xpos   float    The position to draw the centre of the beam at
 *                   0 -> 1 where 0=left edge, 1=right edeg of plot.
 *   ypos   float    The position to draw the centre of the beam at.
 *                   0 -> 1 where 0=bottom edge, 1=top edeg of plot.
 *   xmin   float    If xmin>0.0f and the beam X-extent is less than xmin
 *                   the beam will not be plotted and 1 will be returned.
 *   xmax   float    If xmax>0.0f and the beam X-extent is greater than xmax
 *                   the beam will not be plotted and 1 will be returned.
 *   ymin   float    If ymin>0.0f and the beam Y-extent is less than ymin
 *                   the beam will not be plotted and 1 will be returned.
 *   ymax   float    If ymax>0.0f and the beam Y-extent is greater than ymax
 *                   the beam will not be plotted and 1 will be returned.
 *
 * Output:
 *   return int      0 - Beam was valid and plotted.
 *                   1 - no beam was plotted.
 */
int plbeam(float bmin, float bmaj, float bpa, float xpos, float ypos,
	   float xmin, float xmax, float ymin, float ymax)
{
  const float margin=0.05f; /* The min dist between beam and plot edge */
                            /* in normalised viewport coordinates */
  Ellipse el;               /* The descriptor of the beam ellipse */
  float xa,xb,ya,yb;        /* Current world coordinate limits of plot */
  float xwid, ywid;         /* abs(xb-xa) and abs(yb-ya) respectively */
/*
 * Check whether the beam position is outside the plot.
 */
  if(xpos < 0.0f || xpos > 1.0f || ypos < 0.0f || ypos > 1.0f)
    return 1;
/*
 * Parameterise the beam ellipse.
 */
  el_define(&el, bmin, bmaj, bpa, 0.0f, 0.0f);
/*
 * Determine the current world coordinates of the plotted map.
 */
  cpgqwin(&xa, &xb, &ya, &yb);
  xwid = fabs(xb-xa);
  ywid = fabs(yb-ya);
/*
 * Is the beam too big to be plotted?
 */
  if((xmax>0.0f && el.xwid > xmax) || (ymax>0.0f && el.ywid > ymax) ||
     (xmin>0.0f && el.xwid < xmin) || (ymin>0.0f && el.ywid < ymin))
    return 1;
/*
 * If necessary move the beam centre in normalised world coords, to
 * enforce a minimum distance between the beam edge and the plot edges.
 */
  xpos = floatmin(floatmax(xpos, el.xwid/2.0/xwid+margin),
	      1.0-el.xwid/2.0/xwid-margin);
  ypos = floatmin(floatmax(ypos, el.ywid/2.0/ywid+margin),
	      1.0-el.ywid/2.0/ywid-margin);
/*
 * Find the world coordinate equivalent values of xpos and ypos.
 */
  xpos = xa + xpos * (xb-xa);
  ypos = ya + ypos * (yb-ya);
/*
 * Plot the beam as a filled and outlined ellipse.
 */
  el_move(&el, xpos, ypos);
  el_plot(&el, 5, 14, 0, 0);
  return 0;
}


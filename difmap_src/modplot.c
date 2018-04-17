#include <stdlib.h>
#include <stdio.h>

#include "logio.h"
#include "vlbconst.h"
#include "model.h"
#include "ellips.h"
#include "cpgplot.h"

/*.......................................................................
 * Plot a single model component if visible within the specified axis
 * limits. If you are calling this function repeatedly to plot a block
 * of components, bracket the block of calls with cpgbbuf() and cpgebuf()
 * calls to buffer the PGPLOT output. This will greatly speed up the plotting.
 * 
 * Input:
 *  cmp   Modcmp *  The descriptor of the model component to be plotted.
 *  xa     float    The X-axis world coordinate of one edge of the plot area.
 *  xb     float    The X-axis world coordinate of the opposite edge (to xa)
 *                  of the plot area.
 *  ya     float    The Y-axis world coordinate of one edge of the plot area.
 *  yb     float    The Y-axis world coordinate of the opposite edge (to ya)
 *                  of the plot area.
 *  erase    int    If true, erase the component.
 * Output:
 *  return   int    0 - The center of the component was outside the plot area.
 *                  1 - The center of the component was inside the plot area.
 */
int cmpplot(Modcmp *cmp, float xa, float xb, float ya, float yb, int erase)
{
  const int fpcol=10;  /* The color to plot positive fixed components. */
  const int fncol=2;   /* The color to plot negative fixed components. */
  const int vpcol=7;   /* The color to plot positive variable components. */
  const int vncol=8;   /* The color to plot negative variable components. */
  const int delta_pt=2;/* The PGPLOT marker symbol for delta components */
  int visible;  /* True if the center of the component is visible */
  int color;    /* The color index to plot the component with */
  float xc,yc;  /* The coordinates of the center of the component */
/*
 * Sanity check the component.
 */
  if(cmp==NULL) {
    lprintf(stderr, "cmpplot: NULL component intercepted.\n");
    return 0;
  };
/*
 * Make xa<xb && ya<yb.
 */
  if(xa > xb) {float xtmp = xa; xa = xb; xb = xtmp;};
  if(ya > yb) {float ytmp = ya; ya = yb; yb = ytmp;};
/*
 * Get the coordinates of the center of the pixel.
 */
  xc = cmp->x;
  yc = cmp->y;
/*
 * Buffer the plot.
 */
  cpgbbuf();
/*
 * Set the color index appropriate to the sign of the component flux,
 * whether it is fixed or variable, and whether it is to be erased or
 * drawn.
 */
  if(erase)
    color = 0;
  else if(cmp->freepar)
    color = cmp->flux>0.0f ? vpcol:vncol;
  else 
    color = cmp->flux>0.0f ? fpcol:fncol;
/*
 * Is the center of the component visible?
 */
  visible = (xc>=xa && xc<=xb && yc>=ya && yc<=yb);
/*
 * Plot the delta function components with '+' symbols and the rest
 * with ellipses showing the dimensions of the component.
 */
  if(cmp->type==M_DELT) {
    if(visible) {
      cpgsci(color);
      cpgpt(1, &xc, &yc, delta_pt);
    };
  } else {
    Ellipse el;  /* Descriptor of the elliptical aspect of the component */
/*
 * Describe the elliptical characteristics of the component.
 */
    el_define(&el, cmp->major * cmp->ratio, cmp->major, cmp->phi, xc, yc);
/*
 * If the ellipse is at least partially visible, plot its outline.
 */
    if(el_visible(&el, xa, xb, ya, yb, EL_PART))
      el_plot(&el, color, -1, 1, 0);
  };
/*
 * End pgplot buffering.
 */
  cpgebuf();
  return visible;
}

/*.......................................................................
 * Plot either or both of the fixed and variable parts of the established
 * and tentative models of an observation.
 *
 * Input:
 *  mod       Model * The descriptor of the model to be plotted, (NULL is ok).
 *  dofix       int   If true plot the fixed model.
 *  dovar       int   If true plot the variable model.
 *  xa        float   The X-axis world coordinate of one edge of the plot area.
 *  xb        float   The X-axis world coordinate of the edge opposite to xa.
 *  ya        float   The Y-axis world coordinate of one edge of the plot area.
 *  yb        float   The Y-axis world coordinate of the edge opposite to ya.
 * Output:
 *  return      int   The number of components outside the plot area.
 */
int modplot(Model *mod, int dofix, int dovar, float xa, float xb, float ya, float yb)
{
  Modcmp *cmp;         /* The component being plotted */
  int nhidden=0;       /* Number of components outside of the plot */
  int oldcol;          /* The color on entry to this routine */
/*
 * No model to be plotted?
 */
  if(mod==NULL || mod->ncmp < 1)
    return 0;
/*
 * Arrange for xa<xb && ya<yb.
 */
  if(xa > xb) {float xtmp = xa; xa = xb; xb = xtmp;};
  if(ya > yb) {float ytmp = ya; ya = yb; yb = ytmp;};
/*
 * Store original PGPLOT attributes.
 */
  cpgqci(&oldcol);
/*
 * Plot the components.
 */
  cpgbbuf();
  for(cmp=mod->head; cmp!=NULL; cmp=cmp->next) {
    if((!cmp->freepar && dofix) || (cmp->freepar && dovar)) {
      if(!cmpplot(cmp, xa, xb, ya, yb, 0))
	nhidden++;
    };
  };
  cpgebuf();
/*
 * Restore the original color index.
 */
  cpgsci(oldcol);
/*
 * Return the number of points that lay beyond the plot bounds.
 */
  return nhidden;
}

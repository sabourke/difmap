#include "vlbmath.h"
/*.......................................................................
 * Find the min/max value in a patch of a 2D array. NB the x and y bounds
 * of the patch will be checked against each other and xdim and ydim
 * and enforced to lie within the array.
 *
 * Input:
 *   map    float *  Pointer to first element of 1D array of xdim*ydim
 *                   elements.
 *   xdim   int      The number of float elements along the X axis.
 *   ydim   int      The number of float elements along the Y axis.
 *   xa     int      The smallest X element of the patch.
 *   xb     int      The highest X element of the patch.
 *   ya     int      The smallest Y element of the patch.
 *   yb     int      The highest Y element of the patch.
 * Output:
 *   mapmin float *  Pointer to variable into which the min value of the
 *                   array is deposited.
 *   mapmax float *  Pointer to variable into which the max value of the
 *                   array is deposited.
 */
void imran(float *map, int xdim, int ydim, int xa, int xb, int ya, int yb,
	    float *mapmin, float *mapmax)
{
  float *mapptr;  /* Pointer into 'map' */
  float vmin,vmax;/* Internal versions of mapmin,mapmax */
  int xskip;      /* Offset from end of one row to start of next */
  int xwid;       /* Number of elements across patch */
  int ywid;       /* Number of elements down patch */
  int itmp;       /* Used as intermediary to swap patch bounds if wrong */
  int ix,iy;      /* Pixel coordinates in 'map' */
/*
 * Swap xa with xb and/or ya with yb if necessary.
 */
  if(xa > xb) {
    itmp = xa;
    xa = xb;
    xb = itmp;
  };
  if(ya > yb) {
    itmp = ya;
    ya = yb;
    yb = itmp;
  };
/*
 * Enforce bounds on xa,xb and ya,yb.
 */
  if(xa < 0) xa = 0;
  if(xa >= xdim) xa = xdim - 1;
  if(xb < 0) xb = 0;
  if(xb >= xdim) xb = xdim - 1;
  if(ya < 0) ya = 0;
  if(ya >= ydim) ya = ydim - 1;
  if(yb < 0) yb = 0;
  if(yb >= ydim) yb = ydim - 1;
/*
 * Determine the width and height of the patch in pixels.
 */
  xwid = xb - xa + 1;
  ywid = yb - ya + 1;
/*
 * Find the pointer to the start of the patch (at X=xa,Y=ya).
 */
  mapptr = map + xa + ya * xdim;
/*
 * Work out the increment between the end of one row of the patch and
 * the start of the next.
 */
  xskip = xdim - xwid;
/*
 * Initialise min/max using first element of the patch.
 */
  vmin = vmax = *mapptr;
/*
 * Step through the array searching for the min/max values.
 */
  for(iy=ya; iy<=yb; iy++) {
    for(ix=xa; ix<=xb; ix++) {
      if(*mapptr < vmin)
	vmin = *mapptr;
      if(*mapptr > vmax)
	vmax = *mapptr;
      mapptr++;
    };
    mapptr += xskip;
  };
/*
 * Return the results.
 */
  *mapmin = vmin;
  *mapmax = vmax;
  return;
}

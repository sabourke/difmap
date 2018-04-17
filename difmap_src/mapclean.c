#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "logio.h"
#include "vlbconst.h"
#include "vlbmath.h"
#include "mapmem.h"
#include "mapwin.h"
#include "model.h"
#include "mapcln.h"

static float *absmax(float *map, int nx, int ny, Winran *wins, int nwin);

static void subcc(MapBeam *mb, float *cmpptr, float cmpval,
		  int ixmin, int ixmax, int iymin, int iymax);

/*.......................................................................
 * Clean a map using the beam and windows given and with optional
 * flux and iteration limits.
 *
 * Input:
 *  ob  Observation *  The observation responsible for the map.
 * Input/Ouput:
 *  mb      MapBeam *  The Map and beam. On output the map will be the
 *                     residual map left after the subtraction of all
 *                     components.
 * Input:
 *  win      Subwin *  The container of a list of clean windows. If this
 *                     is NULL or contains no windows then the whole area
 *                     defined in 'mb' will be searched for components.
 *  mod       Model *  This may be NULL, in which case a new model will be
 *                     created. Otherwise, new model components will
 *                     be appended to the existing model in 'mod'. Whichever
 *                     model is used will be returned.
 *  maxcmp      int    The absolute value 'maxcmp' specifies the max number
 *                     of components to be found. If this is negative then
 *                     cleaning will stop before the given limit as soon
 *                     as a -ve component is seen - the -ve component won't
 *                     be subtracted.
 *  cutoff    float    The residual flux to stop cleaning at. (Jy/Beam).
 *  gain      float    The CLEAN loop gain < 1.0 .
 *  docomp      int    If true then compress the model components.
 * Output:

 * return     Model *  The mew Model components that were subtracted
 *                     or NULL on error. If an input model was
 *                     specifed then the input model with the new
 *                     components appended will be returned.
 *                     Otherwise, it is a model that was allocated
 *                     internally (also deleted internally on error)
 *                     and should be deleted using function
 *                     del_Model() when no longer required.
 */
Model *mapclean(Observation *ob, MapBeam *mb, Mapwin *mw, int maxcmp,
		float cutoff, float gain, int docomp)
{
  Subwin *win;          /* A window in the list contained in 'mw' */
  Winran *wins;         /* Array of window indexes */
  Winran *swin;         /* A window in 'wins' */
  Model *mod;           /* The new clean model */
  float *absptr;        /* Contains pointer to the pixel in map with the */
                        /* max absolute value in the map */
  float maxval;         /* Value of latest clean component */
  float xval,yval;      /* Coordinate of latest clean component (radians) */
  int nwin;             /* Number of clean windows */
  int ixmin, ixmax;     /* First and last X-pixel of CLEAN area */
  int iymin, iymax;     /* First and last Y-pixel of CLEAN area */
  int xcent,ycent;      /* Pixel coordinate of centre of map grid */
  int cntr;             /* Offset to centre of map and beam arrays */
  float bmax;           /* Value of beam at beam centre */
  float ccsum=0.0f;     /* Sum of clean-component fluxes */
  int noneg;            /* If true then stop at the first negative component */
  int niter;            /* Number of clean iterations completed */
/*
 * Warn user if the map is marked as restored.
 */
  if(mb->ncmp) {
    lprintf(stderr,
	    "mapclean: Warning: You appear to be cleaning a restored map\n");
  };
/*
 * Trap bad clean gain.
 */
  if(gain<=0.0 || gain > 1.0) {
    lprintf(stderr, "mapclean: Ridiculous clean gain: %g\n", gain);
    return NULL;
  };
/*
 * Set up the Model container to receive CLEAN components.
 */
  mod = new_Model();
  if(!mod)
    return NULL;
/*
 * Determine the area of the map that can be cleaned. This must be the
 * inner half of the array so that the beam array is twice as big.
 * Also determine the coordinate of the centre of the map grid - this
 * corresponds to map coordinate 0,0.
 */
  ixmin = mb->nx/4;
  iymin = mb->ny/4;
  ixmax = mb->nx - ixmin-1;
  iymax = mb->ny - iymin-1;
  xcent = mb->nx/2;
  ycent = mb->ny/2;
/*
 * Allocate as many Winran structures as windows in 'mw' or if 'mw'
 * is empty 1 Winran structure to delimit the whole clean area
 * as a window.
 */
  nwin = (mw == NULL || mw->nwin == 0) ? 1:mw->nwin;
  wins = (Winran *) malloc(nwin * sizeof(Winran));
  if(wins == NULL) {
    lprintf(stderr, "Insufficient memory to CLEAN map\n");
    return del_Model(mod);
  };
/*
 * If no windows were specified set one up to select the whole CLEAN
 * area.
 */
  if(mw == NULL || mw->nwin == 0) {
    wins->xa = ixmin;
    wins->xb = ixmax;
    wins->ya = iymin;
    wins->yb = iymax;
  } else {
/*
 * Convert the window limits to element ranges and enforce bounds.
 */
    swin = wins;
    for(win=mw->head; win != NULL; win=win->next) {
      if(win_pix(win, mb, ixmin, ixmax, iymin, iymax, swin) == 0)
	swin++;
      else
	nwin--;
    };
/*
 * Were there any windows in the clean area?
 */
    if(nwin==0) {
      lprintf(stderr,
	      "clean: All your CLEAN windows lie outside the CLEAN area\n");
      lprintf(stderr, "clean: No CLEANing performed.\n");
      free(wins);
      return del_Model(mod);
    };
  };
/*
 * The cutoff is for comparison to absolute values so enforce positivity.
 */
  cutoff = fabs(cutoff);
/*
 * Process maxcmp to see if to stop at the first negative.
 */
  noneg = maxcmp < 0;
  if(noneg) maxcmp = -maxcmp;
/*
 * Determine the offset in floats to the centre of the beam and map arrays.
 */
  cntr = mb->nx/2 + mb->nx * mb->ny/2;
/*
 * Determine the value at the centre of the beam - this will be used
 * to scale pixel values to Jy/Beam.
 */
  bmax = mb->beam[cntr];
  if(bmax == 0.0f) {
    lprintf(stderr, "clean: invalid dirty beam supplied - try using invert\n");
    free(wins);
    return del_Model(mod);
  };
/*
 * CLEAN loop.
 */
  for(niter=0; niter<maxcmp;) {
/*
 * Search the CLEAN windows for the next largest absolute pixel value
 * in the residual map.
 */
    absptr = absmax(mb->map, mb->nx, mb->ny, wins, nwin);
    if(absptr == NULL) {
      lprintf(stderr, "clean: No flux left in map - finishing early\n");
      break;
    };
/*
 * Copy the max value and convert to Jy/beam.
 */
    maxval = *absptr/bmax;    
/*
 * See if the CLEAN loop has converged.
 */
    if(fabs(maxval) <= cutoff) {
      lprintf(stdout, "Clean target residual flux of %g Jy/beam attained\n",
	      cutoff);
      break;
    };
/*
 * Stop if a negative component is detected and negative components
 * were prohibited by the user.
 */
    if(noneg && maxval < 0.0f) {
      lprintf(stdout, "Clean halted at first negative component\n");
      break;
    };
/*
 * Apply the CLEAN gain.
 */
    maxval *= gain;
/*
 * Subtract the component from the map.
 */
    subcc(mb, absptr, maxval, ixmin, ixmax, iymin, iymax);
/*
 * Keep a record of the flux subtracted so far and report it every 50
 * components.
 */
    niter++;
    ccsum += maxval;
    if(niter && niter % 50 == 0)
      lprintf(stdout, "Component: %3.3d  -  total flux cleaned = %g Jy\n",
	     niter, ccsum);
/*
 * Append the new model component to the model.
 */
    xval = ((absptr - mb->map) % mb->nx - xcent) * mb->xinc;
    yval = ((absptr - mb->map) / mb->nx - ycent) * mb->yinc;
    if(add_xycmp(mod, docomp, 0, maxval, xval, yval, 0.0, 0.0, 0.0, M_DELT,
		 0.0, 0.0) == NULL) {
      lprintf(stderr, "Leaving CLEAN early due to memory problems\n");
      break;
    };
  };
/*
 * Determine and record map statistics.
 */
  mapstats(ob, mb);
/*
 * Report some statistics of the cleaned map.
 */
  lprintf(stdout, "Total flux subtracted in %d components = %g Jy\n", niter,
	 ccsum);
  lprintf(stdout, "Clean residual min=%f max=%f Jy/beam\n", mb->minpix.value,
	  mb->maxpix.value);
  lprintf(stdout, "Clean residual mean=%f rms=%f Jy/beam\n", mb->mapmean,
	  mb->maprms);
/*
 * Mark the map as un-restored.
 */
  mb->ncmp = 0;
  free(wins);
  return mod;
}

/*.......................................................................
 * Find and return a pointer to the pixel containing the max absolute
 * value within windowed areas of a 2D array.
 *
 * Input:
 *   map  float *   A pointer to the 1st element of the 2D array.
 *   nx   int       Number of elements along the X-axis.
 *   ny   int       Number of elements along the Y-axis.
 *   wins Winran *  Pointer to first element of an array of window
 *                  ranges.
 *   nwin int       The number of elements in Winran.
 * Output:
 *   return float * The pointer to the max absolute value in the windowed
 *                  region. If this is NULL then no
 *                  non-zero points were found within the windows.
 */
static float *absmax(float *map, int nx, int ny, Winran *wins, int nwin)
{
  float maxabs; /* Max absolute value yet found */
  float *maxptr;/* Pointer to */
  int ix,iy;   /* Pixel being examined */
  float *fptr; /* Pointer into map array */
  int xwid;    /* Number of pixels along X in current window */
  int ywid;    /* Number of pixels along Y in current window */
  int xskip;   /* The number of floats between a right and next left edge */
  int iwin;    /* Number of current CLEAN window */
/*
 * Loop through windows and find and record the max absolute value and its
 * pointer. 
 */
  maxabs = 0.0f;
  maxptr = NULL;
  for(iwin=0; iwin < nwin; iwin++, wins++) {
    xwid = wins->xb - wins->xa + 1;
    ywid = wins->yb - wins->ya + 1;
/*
 * Determine the number of pixels between the right edge of the current
 * window and the left edge of the same window on the next row.
 */
    xskip = nx - xwid;
/*
 * Pointer to corner of current window with lowest address.
 */
    fptr = map + wins->xa + wins->ya * nx;
    for(iy=0; iy < ywid; iy++,fptr += xskip) {
      for(ix=0; ix < xwid; ix++,fptr++) {
	if(*fptr < -maxabs || *fptr > maxabs) {
	  maxabs = (*fptr < 0.0f) ? -(*fptr):*fptr;
	  maxptr = fptr;
	};
      };
    };
  };
/*
 * Return the pointer to the pixel of the max absolute value.
 */
  return maxptr;
}

/*.......................................................................
 * Subtract a model component from the residual map.
 *
 * Input/Output:
 *  mb  MapBeam *  The residual map and beam etc..
 * Input:
 *  cmpptr float * The address of the component pixel in the map.
 *  cmpval float   The component value (Jy/Beam).
 *  ixmin  int     The X pixel at the left edge of the CLEAN area.
 *  ixmax  int     The X pixel at the right edge of the CLEAN area.
 *  iymin  int     The Y pixel at the lower edge of the CLEAN area.
 *  iymax  int     The Y pixel at the upper edge of the CLEAN area.
 */
static void subcc(MapBeam *mb, float *cmpptr, float cmpval,
		  int ixmin, int ixmax, int iymin, int iymax)
{
  int xwid;    /* Width of CLEAN area in pixels along X axis */
  int ywid;    /* Width of CLEAN area in pixels along Y axis */
  int xskip;   /* No. of pixels between left and right edges of CLEAN area */
  int ix,iy;   /* Used to record position within xwid and ywid */
  int cntr;    /* Offset in floats to centre of map and beam arrays */
  float *mptr; /* Pointer into the map */
  float *bptr; /* Pointer into the beam array */
/*
 * Determine the size of the clean area in pixels.
 */
  xwid = ixmax - ixmin + 1;
  ywid = iymax - iymin + 1;
/*
 * Determine the number of pixels between the right edge of the CLEAN area
 * and the left edge of the same window on the next row.
 */
  xskip = mb->nx - xwid;
/*
 * Determine the offset in floats to the centre of the beam and map arrays.
 */
  cntr = mb->nx/2 + mb->nx * mb->ny/2;
/*
 * Pointer to corner of CLEAN area with lowest address.
 */
  mptr = mb->map + ixmin + iymin * mb->nx;
/*
 * Get the pointer to the pixel in the beam array that should correspond
 * to this point.
 */
  bptr = mb->beam + cntr + (mptr - cmpptr);
/*
 * Subtract the component.
 */
  for(iy=0; iy < ywid; iy++) {
    for(ix=0; ix < xwid; ix++)
      *(mptr++) -= *(bptr++) * cmpval;
    mptr += xskip;
    bptr += xskip;
  };
  return;
}

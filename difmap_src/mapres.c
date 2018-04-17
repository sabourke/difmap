#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "logio.h"
#include "vlbconst.h"
#include "vlbutil.h"
#include "mapmem.h"
#include "model.h"
#include "mapres.h"

#define ETSIZ 1024  /* Size of exponential lookup table */
static float exptab[ETSIZ];
static int ready=0; /* Set to 1 after the table has been initialised */

static void gauconv(float maj_a, float min_a, float ang_a,
	     float maj_b, float min_b, float ang_b,
	     float *major, float *minor, float *angle);
static void res_smooth(float *map, int nx, int ny);

/*.......................................................................
 * Restore a clean map from the residual map and model generated by
 * mapclean(). Only delta-function components are supported presently.
 *
 * Input
 *  ob  Observation *
 *  mb      MapBeam *   On input, the Residual map and beam. On output
 *                      the restored map. The beam is not used.
 *  mod       Model *   The model produced by CLEAN.
 *  clnmap    float *   The pointer to the first element of an array
 *                      of length mb->nx * mb->ny to be filled with
 *                      the restored map. There are 3 options for
 *                      this:
 *                       1. NULL - In this case the array will be allocated
 *                          internally and returned as the function
 *                          return value. It is the caller's
 *                          responsibility to free this array later.
 *                       2. mb->map - The residual map will be
 *                          replaced by the clean map.
 *                       3. An externally allocated array.
 *  bmaj      float     The major axis of the restoring beam (radians).
 *  bmin      float     The minor axis of the restoring beam (radians).
 *  bpa       float     The position angle of the major axis of the
 *                      restoring beam (radians).
 *  dosub       int     If this is true then the model will be subtracted
 *                      from the input map rather than added. ie this
 *                      will return a restored map to the residual map.
 *  noresid     int     If true the residual map will be zeroed before
 *                      restoring the model. The resulting map will just
 *                      be that of the model.
 *  dosmth       int    If true then smooth the residual map before
 *                      restoring.
 *  freq       float    The frequency at which to restore the model.
 * Output:
 *  return     float *  The pointer to the first element of the clean map
 *                      or NULL on error.
 */
float *mapres(Observation *ob, MapBeam *mb, Model *mod, float *clnmp,
	      float bmaj, float bmin, float bpa, int dosub, int noresid,
	      int dosmth, float freq)
{
  const float nsigma=4.5; /* Number of sigma in convolution beam patch */
  Modcmp *cmp;   /* A single model component */
  float minfac;  /* Twice the minor axis gaussian reciprocal variance */
  float majfac;  /* Twice the major axis gaussian reciprocal variance */
  float xminor;  /* X-axis pixel contribution along minor axis of gaussian */
  float yminor;  /* Y-axis pixel contribution along minor axis of gaussian */
  float xmajor;  /* X-axis pixel contribution along major axis of gaussian */
  float ymajor;  /* Y-axis pixel contribution along major axis of gaussian */
  float minor;   /* Distance along minor axis of clean beam gaussian */
  float major;   /* Distance along major axis of clean beam gaussian */
  float arg;     /* Exponential argument in computing gaussian beam */
  int iarg;      /* Index into exptab[] array for exp(arg) */
  int ix,iy;     /* X,Y pixel positions */
  float modx,mody;/* Model X,Y position in map (Pixels) */
  int imodx,imody;/* Integer versions of modx,mody */
  float fx,fy;   /* Position in beam array relative to its centre */
  float *rmptr;  /* Pointer into residual map */
  float *maptr;  /* Pointer into map array */
  int mpinc;     /* Elemental increment in map array */
  float flux;    /* Flux of a clean component */
  int xa,xb,ya,yb;/* Indexes of start and end X and Y elements in map */
  float cmin,cmaj,cpa; /* The convolved dimensions of a gaussian component */
  int nxpix;     /* The number of X-axis pixels out to 2 sigma of gaussian */
  int nypix;     /* The number of X-axis pixels out to 2 sigma of gaussian */
  float bfac;    /* Conversion factor FWHM(radians) -> std-dev (pixels) */
  float expconv; /* Conversion factor from exponent argument to exptab index */
  float ftmp;
/*
 * Compute the exponential lookup table if not already initialised.
 */
  expconv = ETSIZ/(0.5*nsigma*nsigma);
  if(!ready) {
    ready=1;
    for(xa=0; xa<ETSIZ; xa++)
      exptab[xa] = exp((double) -xa/expconv);
  };
/*
 * Swap bmin and bmaj if bmin>bmaj.
 */
  if(bmin>bmaj) {ftmp=bmin; bmin=bmaj; bmaj=ftmp;};
/*
 * Compute the conversion factor between gaussian FWHM and standard deviation.
 */
  bfac = 1.0/sqrt(log(256.0));
/*
 * Zero the residuals if requested.
 */
  if(noresid) {
    maptr = clnmp;
    for(ix=0; ix<mb->nx*mb->ny; ix++)
      *(maptr++) = 0.0f;
  } else {
/*
 * Copy the residual map into the clean map. (Unless the clean map
 * and residual map are the same array).
 */
    if(clnmp != mb->map) {
      rmptr = mb->map;
      maptr = clnmp;
      for(ix=0; ix<mb->nx*mb->ny; ix++)
	*(maptr++) = *(rmptr++);
    };
/*
 * Pre-smooth the residuals in the inner cleaned part of the residual map.
 */
    if(dosmth)
      res_smooth(clnmp, mb->nx, mb->ny);
  };
/*
 * Convolve each model component with the clean beam.
 */
  for(cmp=mod->head; cmp != NULL; cmp=cmp->next) {
/*
 * Get the size of the gaussian to be added into the map.
 */
    switch (cmp->type) {
    case M_DELT:
      cmin = bmin;
      cmaj = bmaj;
      cpa  = bpa;
      break;
    case M_GAUS:
/*
 * Find the dimensions of the gaussian resulting from convolving the beam
 * with the component gaussian.
 */
      gauconv(bmin, bmaj, bpa,
	      cmp->ratio*cmp->major, cmp->major, cmp->phi,
	      &cmin, &cmaj, &cpa);
      break;
/*
 * Only gaussians and delta functions are supported.
 */
    default:
      lprintf(stderr, "mapres: Non delta/gaussian function component not supported\n");
      continue;
    };
/*
 * Determine the scale factor to place in front of the exponential part of
 * the gaussian to yield Janskies per beam.
 */
    flux = cmp->flux*bmaj*bmin/(cmin*cmaj);
/*
 * If the component has a spectral index, modify the flux accordingly.
 */
    if(cmp->spcind != 0.0)
      flux *= pow(freq/cmp->freq0, cmp->spcind);
/*
 * Convert the major and minor axis FWHMs to gaussian standard deviations
 * and from them determine the max number of pixels on each axis, required
 * to sample the gaussian out to 'nsigma' sigma (assuming no beam rotation),
 * and the scale factor in front of the gaussian to yield Janskies per beam.
 */
    cmin *= bfac;
    cmaj *= bfac;
    nxpix = nsigma*cmaj/mb->xinc;
    nypix = nsigma*cmaj/mb->yinc;
/*
 * Determine 1/(2.sigma^2).
 */
    minfac = 0.5f/(cmin*cmin);
    majfac = 0.5f/(cmaj*cmaj);
/*
 * Pre-compute the contribution of X and Y pixel distances along the
 * major and minor axes in the gaussian, measured in radians.
 */
    xminor =  mb->xinc * cos(cpa);
    yminor = -mb->yinc * sin(cpa);
    xmajor =  mb->xinc * sin(cpa);
    ymajor =  mb->yinc * cos(cpa);
/*
 * Determine the nearest pixel to the gaussian center and the first/last
 * pixels on each axis to be processed.
 */
    imodx = modx = mb->nx/2 + cmp->x/mb->xinc;
    imody = mody = mb->ny/2 + cmp->y/mb->yinc;
    xa = (nxpix > imodx) ? 0 : (imodx - nxpix);
    xb = imodx + nxpix;
    if(xb >= mb->nx) xb = mb->nx;
    ya = (nypix > imody) ? 0 : (imody - nypix);
    yb = imody + nypix;
    if(yb >= mb->ny) yb = mb->ny;
/*
 * Add in the gaussian.
 */
    maptr = clnmp + xa + ya * mb->nx;
    mpinc = mb->nx - (xb-xa+1);
    for(iy=ya; iy<=yb; iy++) {
      fy = mody-iy;
      for(ix=xa; ix<=xb; ix++) {
	fx = modx-ix;
	minor = xminor * fx + yminor * fy;
	major = xmajor * fx + ymajor * fy;
	arg = minfac*minor*minor+majfac*major*major;
	iarg = arg * expconv;
	if(iarg < ETSIZ) {
	  if(!dosub)
	    *maptr += flux * exptab[iarg];
	  else
	    *maptr -= flux * exptab[iarg];
	};
	maptr++;
      };
      maptr += mpinc;
    };
  };
/*
 * Determine and record map statistics.
 */
  mapstats(ob, mb);
/*
 * Report the flux range of the restored map.
 */
  lprintf(stdout, "Clean map  min=%.5g  max=%.5g Jy/beam\n",
	  mb->minpix.value, mb->maxpix.value);
/*
 * Mark the map as being restored and record the beam parameters.
 */
  mb->ncmp += mod->ncmp;
  mb->bmin = bmin;
  mb->bmaj = bmaj;
  mb->bpa = bpa;
  return clnmp;
}

#define TINY 1.0e-10
/*.......................................................................
 * Given the parameters of two gaussians, return the parameters for a
 * new gaussian consisting of the convolution of the originals.
 * Derived from Wild J. P., Aust. J. Phys., 1970, 23, 113-115.
 *
 * Input:
 *  min_a   float    The minor axis length of the first gaussian.
 *  maj_a   float    The major axis length of the first gaussian.
 *  ang_a   float    The position angle between major axis and vertical.
 *                   (radians).
 *  min_b   float    The minor axis length of the second gaussian.
 *  maj_b   float    The major axis length of the second gaussian.
 *  ang_b   float    The position angle between major axis and vertical.
 *                   (radians).
 * Output:
 *  minor   float *  The minor axis length of the first gaussian.
 *  major   float *  The major axis length of the first gaussian.
 *  angle   float *  The position angle between major axis and vertical.
 *                   (radians).
 */
static void gauconv(float min_a, float maj_a, float ang_a,
		    float min_b, float maj_b, float ang_b,
		    float *minor, float *major, float *angle)
{
  float sum7;  /* The sum on the RHS of equation 7, Wild (1970) */
  float sum8;  /* The sum on the RHS of equation 8, Wild (1970) */
  float sum9;  /* The sum on the RHS of equation 9, Wild (1970) */
  float sumvar;/* Sum of variances of the output gaussian */  
/*
 * All variances are added and subtracted in quadrature so square the major
 * and minor axes.
 */
  maj_a *= maj_a;
  min_a *= min_a;
  maj_b *= maj_b;
  min_b *= min_b;
/*
 * Determine the RHS sums of equations 7,8,9 of Wild (1970).
 */
  sum7 = (maj_a - min_a) * sin(2.0*ang_a) + (maj_b - min_b) * sin(2.0*ang_b);
  sum8 = (maj_a + min_a) + (maj_b + min_b);
  sum9 = (maj_a - min_a) * cos(2.0*ang_a) + (maj_b - min_b) * cos(2.0*ang_b);
/*
 * Determine the angle of the convolved gaussian.
 */
  *angle = (fabs(sum7)==0.0 && fabs(sum9)==0.0) ? 0.0:(0.5*atan2(sum7,sum9));
/*
 * The sum of the squares and difference of the squares of the major and
 * minor axis lengths are given by sqrt(sum7^2+sum9^2) and sum8 respectively.
 */
  sumvar = sqrt(sum7*sum7 + sum9*sum9);
/*
 * Determine the new major and minor axis lengths.
 */
  *major = sqrt(0.5f * (sum8 + sumvar));
  *minor = sqrt(fabs(0.5f * (sum8 - sumvar)));
  return;
}


/*.......................................................................
 * Smooth the central nx/2 x ny/2 area of the residual map using
 * a fixed 3x3 mask.
 *
 * Input:
 *  map   float *   Map array of size nx x ny.
 *  nx      int     The dimension of the X-axis.
 *  ny      int     The dimension of the Y-axis.
 */
static void res_smooth(float *map, int nx, int ny)
{
  float *mapptr;  /* Points to start of map area in 'map' */
  float *endptr;  /* Points to last pixel of map in 'map' */
  float *dest_ptr;/* Points to the destination element in the 'map' array */
  float *tmp_ptr; /* Temporary pointer into 'map' array */
  size_t yinc;    /* No. elements between end of one row and start of next */
  size_t yshift;  /* The no. of elements by which the map array is shifted */
  size_t mask_offset; /* Map offset between corner and center of mask */
  int ix,iy;  /* X and Y array indexes */
  int mx,my;  /* X and Y indexes into 'mask' */
  int xa,xb;  /* Start and end X indexes of map within 'map' */
  int ya,yb;  /* Start and end Y indexes of map within 'map' */
  float sum;  /* Sum of contributions to one smoothed pixel */
/*
 * Define the smoothing mask.
 * Note that N_MASK must be an odd number.
 */
#define N_MASK 3
  float mask[N_MASK][N_MASK]={
    {0.0625,0.125,0.0625},
    {0.125, 0.25, 0.125},
    {0.0625,0.125,0.0625}
  };
  int nhalf=N_MASK/2;
/*
 * Determine the start and end indexes of the central map.
 */
  xa = nx/4;
  ya = ny/4;
  xb = 3*xa - 1;
  yb = 3*ya - 1;
/*
 * Work out pointers to the start and end pixels of the central map.
 */
  mapptr = map + xa + ya * nx;
  endptr = map + xb + yb * nx;
/*
 * Work out the increment in elements to get from the end of a row
 * to the start of the next row.
 */
  yinc = nx/2;
/*
 * First shift the central map patch by yinc elements "right"
 * into the unused margin of the map grid so that the map can be
 * smoothed in-place.
 */
  yshift = yinc;
  dest_ptr = mapptr + yshift;
  for(iy=ya; iy<=yb; iy++, dest_ptr+=yinc) {
    for(ix=xa; ix<=xb; ix++, dest_ptr++)
      *dest_ptr = *(dest_ptr-yshift);
  };
/*
 * Now smooth the map back into its original position.
 * Note that we will not attempt to smooth a margin of width 'nhalf'
 * around the edges of the map.
 *
 * First determine the map element offset between the center of the mask
 * and the first corner of the mask.
 */
  mask_offset = nhalf + nx * nhalf;
  dest_ptr = mapptr + mask_offset;  /* Start at [nhalf][nhalf] */
  yinc = nx/2 + 2*nhalf;
  for(iy=ya+nhalf; iy<=yb-nhalf; iy++, dest_ptr+=yinc) {
    for(ix=xa+nhalf; ix<=xb-nhalf; ix++, dest_ptr++) {
/*
 * Find the first (corner) element of the patch of the shifted image to
 * be used in smoothing
 */
      tmp_ptr = dest_ptr + yshift - mask_offset;
/*
 * Sum mask-weighted neighbouring elements of the map.
 */
      sum=0.0f;
      for(my=0; my<N_MASK; my++,tmp_ptr+=nx-N_MASK)
	for(mx=0; mx<N_MASK; mx++,tmp_ptr++)
	  sum += *tmp_ptr * mask[my][mx];
      *dest_ptr = sum;
    };
  };
  return;
}

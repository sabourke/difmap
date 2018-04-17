#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "logio.h"
#include "mapmem.h"

static int ispow2(int n);
static UVbin *new_UVbin(MapBeam *mb, int nu, int nv);
static UVbin *del_UVbin(MapBeam *mb);

/*.......................................................................
 * Allocate a MapBeam struct and map and beam arrays of a specified
 * size to be contained by it.
 *
 * Input:
 *   oldmap MapBeam *  Optional previous MapBeam to be reused if
 *                     nx and ny haven't changed and deleted otherwise.
 *                     If not required, send NULL.
 *   nx         int    The number of pixels along the X axis of the map
 *                     and beam arrays. The actual number allocated will
 *                     be nx+2, to accomodate the use of the grid for a
 *                     half-conjugate real FFT from the UV plane.
 *   xinc     float    The width of a cell along the X-axis of the map
 *                     and beam arrays (radians).
 *   ny         int    The number of pixels along the Y axis of the map
 *                     and beam arrays.
 *   yinc     float    The width of a cell along the Y-axis of the map
 *                     and beam arrays (radians).
 * Output:
 *   return MapBeam *  Pointer to the initialised container struct or
 *                     NULL if there is a memory allocation failure.
 */
MapBeam *new_MapBeam(MapBeam *oldmap, int nx, float xinc, int ny, float yinc)
{
  MapBeam *mb;    /* Pointer to new container */
  size_t nbytes;  /* Number of bytes required for the map/beam arrays */
  float *mptr;    /* A pointer into the map array */
  float *bptr;    /* A pointer into the beam array */
  size_t i;
/*
 * Check arguments.
 */
  if(nx<=32 || !ispow2(nx) || ny<=32 || !ispow2(ny)) {
    lprintf(stderr, "new_MapBeam: map grid size must be a power of 2 > 32.\n");
    return del_MapBeam(oldmap);
  };
  if(xinc <= 0.0f || yinc <= 0.0f) {
    lprintf(stderr, "new_MapBeam: cell-size must be finite and positive.\n");
    return del_MapBeam(oldmap);
  };
/*
 * Allocate the container and the map and beam arrays if the map
 * size has changed or never been previously allocated.
 */
  if(oldmap != NULL && oldmap->nx == nx && oldmap->ny == ny) {
    mb = oldmap;
  } else {
    del_MapBeam(oldmap);
/*
 * Get the container.
 */
    mb = (MapBeam *) malloc(sizeof(MapBeam));
    if(mb == NULL) {
      lprintf(stderr, "Insufficient memory for a new Map and Beam\n");
      return mb;
    };
/*
 * NULL all pointers such that del_MapBeam() can hereafter be called.
 */
    mb->map = NULL;
    mb->beam = NULL;
    mb->rxft = NULL;
    mb->ryft = NULL;
    mb->bin = NULL;
/*
 * We require two 2D arrays of size (nx+2)*ny for the half of the
 * conjugate symmetric UV array used in gridding. We also need extra
 * work arrays of length 'nx' and 'ny' to hold the normalized reciprocal
 * fourier transforms along the respective axes, of the gridding
 * convolution function, for use in uvtrans().
 */
    nbytes = (nx+2)*(ny+1)*sizeof(float);
    mb->map = (float *) malloc(nbytes);
    mb->beam = (float *) malloc(nbytes);
    mb->rxft = (float *) malloc((nx+1) * sizeof(float));
    mb->ryft = (float *) malloc((ny+1) * sizeof(float));
    if(mb->map==NULL || mb->beam==NULL || mb->rxft==NULL || mb->ryft==NULL) {
      lprintf(stderr,"Insufficient memory for map/beam of size %d,%d\n", nx,ny);
      return del_MapBeam(mb);
    };
  };
/*
 * Zero the map and beam arrays.
 */
  mptr = mb->map;
  bptr = mb->beam;
  for(i=0; i<(nx+2)*ny; i++) {
    *(mptr++) = 0.0f;
    *(bptr++) = 0.0f;
  };
/*
 * Initialize the recorded map statistics.
 */
  {
    const Mappix nopix = {0.0,0.0,0.0,0.0,0.00,0};
    mb->maxpix = mb->minpix = nopix;
    mb->maprms = 0.0f;
    mb->mapmean = 0.0f;
    mb->mapflux = 0.0f;
    mb->noise = 0.0f;
  };
/*
 * Install other parameters.
 */
  mb->xinc = xinc;
  mb->yinc = yinc;
  mb->uinc = 1.0/(xinc*nx);
  mb->vinc = 1.0/(yinc*ny);
  mb->nx = nx;
  mb->ny = ny;
  mb->ncmp = 0;
  mb->domap = mb->dobeam = 1;
  mb->bmin = mb->bmaj = mb->bpa = 0.0f;
  mb->e_bmin = mb->e_bmaj = mb->e_bpa = 0.0f;
  mb->maparea.ixmin = mb->nx/4;
  mb->maparea.ixmax = mb->nx - mb->nx/4 - 1;
  mb->maparea.iymin = mb->ny/4;
  mb->maparea.iymax = mb->ny - mb->ny/4 - 1;
/*
 * Allocate and/or initialize the bin array container.
 */
  if(new_UVbin(mb, nx/4, ny/2)==NULL)
    return del_MapBeam(mb);
  return mb;
}

/*.......................................................................
 * Delete a MapBeam container and its contained map and beam arrays.
 *
 * Input:
 *  mb      MapBeam *   The container to be deleted.
 * Output:
 *  return  MapBeam *   Always NULL, ie (MapBeam *) 0 .
 */
MapBeam *del_MapBeam(MapBeam *mb)
{
  if(mb) {
    if(mb->map)
      free(mb->map);
    if(mb->beam)
      free(mb->beam);
    if(mb->rxft)
      free(mb->rxft);
    if(mb->ryft)
      free(mb->ryft);
    del_UVbin(mb);
    free(mb);
  };
  return NULL;
}

/*.......................................................................
 * Create/initialize a UVbin container in a given MapBeam structure.
 *
 * Input:
 *  mb     MapBeam *  The container of the UVbin array. mb->bin must be
 *                    NULL unless it has been previously allocated.
 *  nu         int    The number of bins required on the U axis.
 *                    This should be equal to mb->nx/4.
 *  nv         int    The number of bins required on the V axis.
 *                    This should be equal to mb->ny/2.
 * Output:
 *  return   UVbin *  The initialized container, or NULL on error.
 *                    mb->bin will be deleted on error.
 */
static UVbin *new_UVbin(MapBeam *mb, int nu, int nv)
{
  UVbin *bin;   /* The new descriptor */
  Bincell *bc;  /* Pointer to bin->bins[] */
/*
 * Allocate a new descriptor?
 */
  if(mb->bin==NULL) {
    mb->bin = malloc(sizeof(UVbin));
    if(mb->bin==NULL) {
      lprintf(stderr, "Insuficient memory to allocate UV bin array.\n");
      return del_UVbin(mb);
    };
/*
 * Default initialize all contained members.
 */
    bin = mb->bin;
    bin->bins = NULL;
    bin->nu = bin->nv = bin->nbin = 0;
    bin->utopix = bin->vtopix = 0.0f;
  };
/*
 * Get a pointer to the new descriptor.
 */
  bin = mb->bin;
/*
 * Determine and record the dimensions of the required bin array.
 */
  bin->nu = nu;
  bin->nv = nv;
  bin->nbin = bin->nu * bin->nv;
/*
 * Allocate or re-size the bin array in the UVbin container.
 */
  if(bin->bins==NULL)
    bc = malloc(sizeof(Bincell) * bin->nbin);
  else
    bc = realloc(bin->bins, sizeof(Bincell) * bin->nbin);
  if(bc==NULL) {
    lprintf(stderr, "Insufficient memory for new UV bin array.\n");
    return del_UVbin(mb);
  };
/*
 * Install the new/re-sized array.
 */
  bin->bins = bc;
/*
 * Assign 0 to the conversion factors from U and V in wavelengths to
 * bin cell widths. These will are re-assigned in uvbin().
 */
  bin->utopix = 0.0f;
  bin->vtopix = 0.0f;
/*
 * Return the initialized container.
 */
  return bin;
}

/*.......................................................................
 * Delete the UVbin descriptor contained in a MapBeam container.
 *
 * Input:
 *  mb   MapBeam *  The container of the descriptor to be deleted.
 * Output:
 *  return UVbin *  The deleted descriptor ie. NULL.
 */
static UVbin *del_UVbin(MapBeam *mb)
{
  if(mb && mb->bin) {
    UVbin *bin = mb->bin;
/*
 * Delete the bin array.
 */
    if(bin->bins)
      free(bin->bins);
/*
 * Free the container.
 */
    free(bin);
    mb->bin = NULL;
  };
  return NULL;
}

/*.......................................................................
 * Private function of new_MapBeam() used to determine whether a given
 * grid-size is a power of 2. This method is a lot faster and more
 * robust than using logs etc..
 *
 * Input:
 *  n        int   The grid size.
 * Output:
 *  return   int   0 - n is not a power of 2.
 *                 1 - n is a power of 2.
 */
static int ispow2(int n)
{
/*
 * Negative numbers are not powers of 2.
 */
  if(n<=0)
    return 0;
/*
 * Divide by 2 until the remainder of dividing by two is non-zero, or
 * 1 is reached. 1 is only reached if the number is a power of two.
 */
  while(n>1 && n%2==0)
    n >>= 1;
  return n==1;
}

/*.......................................................................
 * Determine the statistics of the cleanable part of a map,
 * and record them in mb->maxpix, mb->minpix, mb->maprms and mb->mapmean.
 * This function should be called whenever the map is changed.
 *
 * Input:
 *  ob  Observation *  The observation responsible for the current map.
 *  mb      MapBeam *  The map and beam descriptor/container.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int mapstats(Observation *ob, MapBeam *mb)
{
  Mappix pmin;    /* Details of minimum valued pixel in map */
  Mappix pmax;    /* Details of maximum valued pixel in map */
  int xa,xb;      /* Indexes of the first and last pixels along the X-axis */
  int ya,yb;      /* Indexes of the first and last pixels along the Y-axis */
  int ix,iy;      /* The indexes of the pixel being checked */
  int xskip;      /* Number of margin pixels to skip after scanning X-axis */
  float mean_flux;/* Running mean of residual flux */
  float mean_sqr; /* Running mean square offset from the mean flux */
  size_t nsum;    /* Number of pixels in running means */
  float *mapptr;  /* Pointer into the map array */
/*
 * Check the descriptor.
 */
  if(!mb || !ob) {
    lprintf(stderr, "mapstats: NULL argument(s).\n");
    return 1;
  };
/*
 * Determine the range of array indexes required to cover the cleanable
 * area of the map.
 */
  xa = mb->maparea.ixmin;
  ya = mb->maparea.iymin;
  xb = mb->maparea.ixmax;
  yb = mb->maparea.iymax;
/*
 * Find the pointer to the start of the patch (at X=xa,Y=ya).
 */
  mapptr = mb->map + xa + ya * mb->nx;
/*
 * Work out the increment between the end of one row of the patch and
 * the start of the next.
 */
  xskip = mb->nx - (xb-xa+1);
/*
 * Initialise min/max using first element of the patch.
 */
  pmin.value = pmax.value = *mapptr;
  pmin.ix = pmax.ix = xa;
  pmin.iy = pmax.iy = ya;
/*
 * Initialize the mean flux.
 */
  mean_flux = 0.0f;
  nsum = 0;
/*
 * Step through the array searching for the min/max values and accumulate
 * the running mean pixel flux.
 */
  for(iy=ya; iy<=yb; iy++) {
    for(ix=xa; ix<=xb; ix++) {
      float value = *mapptr++;
/*
 * Update the min/max pixels.
 */
      if(value > pmax.value)  {
	pmax.value = value;
	pmax.ix = ix;
	pmax.iy = iy;
      } else if(value < pmin.value) {
	pmin.value = value;
	pmin.ix = ix;
	pmin.iy = iy;
      };
/*
 * Update the running mean flux.
 */
      nsum++;
      mean_flux += (value - mean_flux) / nsum;
    };
    mapptr += xskip;
  };
/*
 * Go through the map a second time to find the mean square difference
 * between pixel fluxes and the mean pixel flux.
 */
  mean_sqr = 0.0f;
  nsum = 0;
  mapptr = mb->map + xa + ya * mb->nx;
  for(iy=ya; iy<=yb; iy++) {
    for(ix=xa; ix<=xb; ix++) {
      float value = *mapptr++ - mean_flux;
/*
 * Update the running mean.
 */
      nsum++;
      mean_sqr += (value * value - mean_sqr) / nsum;
    };
    mapptr += xskip;
  };
/*
 * Convert from position units of pixel indexes to units of radians
 * wrt the map center.
 */
  pmin.xpos = (pmin.ix - mb->nx/2) * mb->xinc;
  pmin.ypos = (pmin.iy - mb->ny/2) * mb->yinc;
  pmax.xpos = (pmax.ix - mb->nx/2) * mb->xinc;
  pmax.ypos = (pmax.iy - mb->ny/2) * mb->yinc;
/*
 * Calculate the RA and Dec equivalents of the pixel coordinates.
 */
  pmin.ra = lmtora(ob->source.ra, ob->source.dec,
		   -ob->geom.east + map_x_pixel_to_coord(mb, pmin.ix),
		   -ob->geom.north + map_y_pixel_to_coord(mb, pmin.iy),
		   ob->proj);
  pmin.dec = lmtodec(ob->source.ra, ob->source.dec,
		     -ob->geom.east + map_x_pixel_to_coord(mb, pmin.ix),
		     -ob->geom.north + map_y_pixel_to_coord(mb, pmin.iy),
		     ob->proj);
  pmax.ra = lmtora(ob->source.ra, ob->source.dec,
		   -ob->geom.east + map_x_pixel_to_coord(mb, pmax.ix),
		   -ob->geom.north + map_y_pixel_to_coord(mb, pmax.iy),
		   ob->proj);
  pmax.dec = lmtodec(ob->source.ra, ob->source.dec,
		     -ob->geom.east + map_x_pixel_to_coord(mb, pmax.ix),
		     -ob->geom.north + map_y_pixel_to_coord(mb, pmax.iy),
		     ob->proj);
/*
 * Record the results.
 */
  mb->minpix = pmin;
  mb->maxpix = pmax;
  mb->maprms = sqrt(mean_sqr);
  mb->mapmean = mean_flux;
  mb->mapflux = mean_flux * nsum;
  return 0;
}

/*.......................................................................
 * Convert a map x coordinate (radians) wrt the map center to the
 * corresponding 2D map array X-index of the nearest pixel.
 *
 * Input:
 *  mb  MapBeam * The map/beam container.
 *  x     float   The X-axis coordinate (radians).
 * Output:
 *  return  int   The index of the nearest pixel.
 */
int map_x_coord_to_pixel(MapBeam *mb, float x)
{
  return mb ? (mb->nx/2 + (int) floor(x/mb->xinc + 0.5)) : 0;
}

/*.......................................................................
 * Convert a map y coordinate (radians) wrt the map center to the
 * corresponding 2D map array Y-index of the nearest pixel.
 *
 * Input:
 *  mb  MapBeam * The map/beam container.
 *  y     float   The Y-axis coordinate (radians).
 * Output:
 *  return  int   The index of the nearest pixel.
 */
int map_y_coord_to_pixel(MapBeam *mb, float y)
{
  return mb ? (mb->ny/2 + (int) floor(y/mb->yinc + 0.5)) : 0;
}

/*.......................................................................
 * Convert a 2D map-array X-index to the corresponding x coordinate
 * (radians) wrt the map center.
 *
 * Input:
 *  mb  MapBeam * The map/beam container.
 *  ix      int   The X-axis array index.
 * Output:
 *  return  int   The X-coordinate of the pixel.
 */
float map_x_pixel_to_coord(MapBeam *mb, int ix)
{
  return mb ? (ix - mb->nx/2) * mb->xinc : 0.0f;
}

/*.......................................................................
 * Convert a 2D map-array Y-index to the corresponding y coordinate
 * (radians) wrt the map center.
 *
 * Input:
 *  mb  MapBeam * The map/beam container.
 *  iy      int   The Y-axis array index.
 * Output:
 *  return  int   The Y-coordinate of the pixel.
 */
float map_y_pixel_to_coord(MapBeam *mb, int iy)
{
  return mb ? (iy - mb->ny/2) * mb->yinc : 0.0f;
}


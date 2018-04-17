/* Container class for a map and its associated dirty beam */

#ifndef mapmem_h
#define mapmem_h

#include "obs.h"

/*
 * Define a container for a uniform-weighting bin array and associated
 * parameters. This retains information set by a call to uvbin()
 * such that it may be used in subsequent calls to getuvbin() in the
 * gridding function uvgrid().
 */
typedef int Bincell;
typedef struct {
  Bincell *bins;   /* Array of nbin = nu U-coord * nv V-coord UV bins */
  int nu,nv;       /* Dimension of U and V bin array axes */
  int nbin;        /* The total number of elements in bins = nu*nv */
  float utopix;    /* Conversion factor from U (wavelengths) to bin index */
  float vtopix;    /* Conversion factor from V (wavelengths) to bin index */
} UVbin;

/* Define a type used to record details about a map pixel */

typedef struct {
  float value;     /* Value of pixel */
  float xpos,ypos; /* Coordinates of pixel wrt map center (radians) */
  double ra, dec;  /* The Right Ascension and declination of the pixel */
  int ix,iy;       /* Array coordinate indexes of pixel */
} Mappix;

/*
 * The following container records the 2D array indexes of the corners
 * of the map area within the map array.
 */
typedef struct {
  int ixmin,ixmax;       /* Min/max X indexes of map area */
  int iymin,iymax;       /* Min/max Y indexes of map area */
} MapArea;

typedef struct {
  float *map;       /* Map array */
  float *beam;      /* Beam array */
  Mappix maxpix;    /* Details of max value'd pixel in map */
  Mappix minpix;    /* Details of min value'd pixel in map */
  float maprms;     /* RMS flux in map */
  float mapmean;    /* Mean flux in map */
  float mapflux;    /* Total flux in map */
  float noise;      /* The rms of the map noise predicted from the weights */
  int domap;        /* If 1 then the map is out of date wrt the UV data */
  int dobeam;       /* If 1 then the beam is out of date wrt the UV data */
  int nx,ny;        /* Map and Beam have (nx+2)*ny elements */
  int ncmp;         /* The number of model components restored in map[] */
  float xinc;       /* Map cellsize along X-axis (radians) */
  float yinc;       /* Map cellsize along Y-axis (radians) */
  float uinc;       /* UV grid U axis cell-size (wavelengths) */
  float vinc;       /* UV grid V axis cell-size (wavelengths) */
  float bmin;       /* The beam minor axis (radians) last used by restore() */
  float bmaj;       /* The beam major axis (radians) last used by restore() */
  float bpa;        /* Beam position angle (radians) last used by restore() */  
  float e_bmin;     /* Estimated beam minor axis (radians) from grid-weights */
  float e_bmaj;     /* Estimated beam major axis (radians) from grid-weights */
  float e_bpa;      /* Estimated beam posn. angle (radians) from grid-weights */
  float *rxft;      /* Reciprocal Fourier transform of X-axis gridding fn. */
  float *ryft;      /* Reciprocal Fourier transform of Y-axis gridding fn. */
  MapArea maparea;  /* The 2D pixel bounds of the map within 'mb->map' */
  UVbin *bin;       /* Uniform-weighting bin array */
} MapBeam;

MapBeam *new_MapBeam(MapBeam *oldmap, int nx, float xinc, int ny, float yinc);
MapBeam *del_MapBeam(MapBeam *mb);

/*
 * Record map statistics in mb->maxpix and mb->minpix,
 * and mb->maprms, mb->mapmean.
 */

int mapstats(Observation *ob, MapBeam *mb);

/*
 * Center relative radian coordinate <-> 2D array index converters.
 */

int map_x_coord_to_pixel(MapBeam *mb, float x);
int map_y_coord_to_pixel(MapBeam *mb, float y);
float map_x_pixel_to_coord(MapBeam *mb, int ix);
float map_y_pixel_to_coord(MapBeam *mb, int iy);

#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "logio.h"
#include "obs.h"
#include "mapmem.h"
#include "mapcor.h"

/*
 * Specify how many samples of the primary beam to compute along the range
 * of radial distances from the pointing center covered by the map.
 */
#define PB_COR_NSAMPLE 512

/*.......................................................................
 * Divide a given map by the combined primary beam of all baselines.
 *
 * Input:
 *  ob   Observation *  The observation that was used to generate the map.
 *  mb       MapBeam *  The container of the map to be corrected.
 *  cutoff     float    The primary beam value below which pixels are
 *                      replaced with zero flux.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int pb_cor_map(Observation *ob, MapBeam *mb, float cutoff)
{
  float *pb;         /* An array of npix samples of the primary beam */
  float factor;      /* A primary beam factor */
  float px, py;      /* The pixel location of the pointing center */
  float r1,r2,r3,r4; /* The radial distances of each of the four corners */
                     /* of the map from the pointing center. */
  float rmin,rmax;   /* The minimum and maximum of r1,r2,r3,r4 */
  float dr;          /* The step size in radial distance */
  float *mapptr;     /* A pointer into the map array */
  int xa,xb;         /* Indexes of the first and last pixels along the X-axis */
  int ya,yb;         /* Indexes of the first and last pixels along the Y-axis */
  int ix,iy;         /* The indexes of the pixel being checked */
  int xskip;         /* Number of margin pixels to skip after scanning X-axis */
  int pbmax;         /* The last index of the pb[] array that is filled */
  float max_radius;  /* The maximum radial distance at which the primary */
                     /*  beam is known not to have fallen below 1/max_cor */
  int i;
/*
 * Check the arguments.
 */
  if(!ob_ready(ob, OB_SELECT, "pb_cor_map"))
    return 1;
  if(!mb) {
    lprintf(stderr, "pb_cor_map: NULL argument(s).\n");
    return 1;
  };
/*
 * If no antenna beams are currently defined, there is nothing to be done.
 */
  if(count_antenna_beams(ob->ab) <= 0)
    return 0;
/*
 * Update the per-baseline weights in all IFs.
 */
  if(update_baseline_weights(ob, -1))
    return 1;
/*
 * Determine the range of array indexes required to cover the cleanable
 * area of the map.
 */
  xa = mb->maparea.ixmin;
  ya = mb->maparea.iymin;
  xb = mb->maparea.ixmax;
  yb = mb->maparea.iymax;
/*
 * Get the location of the pointing center, in map pixels.
 */
  px = map_x_coord_to_pixel(mb, ob->geom.east - ob->source.east);
  py = map_y_coord_to_pixel(mb, ob->geom.north - ob->source.north);
/*
 * Compute the radial distances of the four corners of the map area.
 */
  r1 = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xa),
			    map_y_pixel_to_coord(mb, ya));
  r2 = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xa),
			    map_y_pixel_to_coord(mb, yb));
  r3 = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xb),
			    map_y_pixel_to_coord(mb, ya));
  r4 = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xb),
			    map_y_pixel_to_coord(mb, yb));
/*
 * Get the maximum of the above radii.
 */
  rmax = r1;
  if(r2 > rmax)
    rmax = r2;
  if(r3 > rmax)
    rmax = r3;
  if(r4 > rmax)
    rmax = r4;
/*
 * Find the minimum radial distance of any map pixel from the pointing
 * center.
 */
  if(px >= xa && px <= xb) {
/*
 * Is the pointing center also within the y-axis range of the map?
 */
    if(py >= ya && py <= yb) {
      rmin = 0.0;
    } else if(py < ya) {
      rmin = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, px),
				  map_y_pixel_to_coord(mb, ya));
    } else {              /* py > yb */
      rmin = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, px),
				  map_y_pixel_to_coord(mb, yb));
    };
  } else if(px < xa) {
    if(py >= ya && py <= yb) {
      rmin = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xa),
				  map_y_pixel_to_coord(mb, py));
    } else if(py < ya) {
      rmin = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xa),
				  map_y_pixel_to_coord(mb, ya));
    } else {              /* py > yb */
      rmin = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xb),
				  map_y_pixel_to_coord(mb, ya));
    };
  } else {                 /* px > xb */
    if(py >= ya && py <= yb) {
      rmin = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xb),
				  map_y_pixel_to_coord(mb, py));
    } else if(py < ya) {
      rmin = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xb),
				  map_y_pixel_to_coord(mb, ya));
    } else {              /* py > yb */
      rmin = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, xb),
				  map_y_pixel_to_coord(mb, yb));
    };
  };
/*
 * Create an array in which we can compute the circularly symmetric primary
 * beam 
 */
  pb = (float *) malloc(sizeof(*pb) * PB_COR_NSAMPLE);
  if(!pb) {
    lprintf(stderr,
	    "pb_cor_map: Insufficient memory to compute primary beam.\n");
    return 1;
  };
/*
 * Compute the step-size between elements of the primary beam array.
 */
  dr = (rmax - rmin) / (PB_COR_NSAMPLE - 1);
/*
 * Sample the primary beam at each element of pb, covering the range
 * rmin to rmax in radial distance.
 */
  max_radius = rmin;
  for(i=0; i<PB_COR_NSAMPLE; i++) {
    float radius = rmin + i*dr;
    if(pb_scale_factor(ob, radius, &factor)) {
      free(pb);
      return 1;
    };
/*
 * Assign the primary beam factor to the array. Note that this has
 * to be done before the check for the cutoff, because we need one
 * value above the cutoff radius, to be able to interpolate up to
 * max_radius.
 */
    pb[i] = factor;
    pbmax = i;
/*
 * Have we passed the primary beam cutoff?
 */
    if(factor < cutoff)
      break;
/*
 * We are still below the cutoff radius.
 */
    max_radius = radius;
  };
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
 * Use the above array to correct each of the pixels of the map.
 */
  for(iy=ya; iy<=yb; iy++, mapptr += xskip) {
    for(ix=xa; ix<=xb; ix++,mapptr++) {
/*
 * Get the radial distance of the latest pixel from the pointing center.
 */
      float radius = calc_pointing_offset(ob, map_x_pixel_to_coord(mb, ix),
					  map_y_pixel_to_coord(mb, iy));
/*
 * Fill pixels beyond the primary beam cutoff with zeroes.
 */
      if(radius > max_radius) {
	*mapptr = 0.0;
      } else {
/*
 * Determine the two elements of the pb[] array which bracket the above
 * radius.
 */
	float p = radius / dr;
	int ia = floor(p);
	int ib = ceil(p);
/*
 * The following check probably isn't necessary, but it does no harm
 * to be a little paranoid.
 */
	if(ib >= pbmax)
	  ia = ib = pbmax - 1;
/*
 * Interpolate for the primary beam factor.
 */
	factor = pb[ia] + (p - ia) * (pb[ib] - pb[ia]);
/*
 * Correct the map pixel by the primary beam factor.
 */
	*mapptr /= factor;
      };
    };
  };
/*
 * Clean up.
 */
  free(pb);
  return 0;
}

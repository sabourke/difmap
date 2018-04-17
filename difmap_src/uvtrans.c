#include <stdio.h>
#include <stdlib.h>
#include "obs.h"
#include "vlbinv.h"
#include "vlbfft.h"

/*.......................................................................
 * Take the UV grid returned by uvgrid(), phase shift it and inverse
 * transform it to produce a dirty map or beam. The input array should
 * be a half conjugate symmetric array as returned by uvgrid() with
 * nx/2+1,ny complex pairs of floats and having U=0,V=0 at array
 * element 0,0. The returned array will be the square map or beam in the
 * first nx*ny floats. The extra floats at the end will be zeroed.
 *
 * Input/Output:
 *  mb    MapBeam *  The map/beam container, containing the gridded UV
 *                   data to be inverted, and sensitivity
 *                   deconvolution functions.
 * Input:
 *  domap     int    If true, invert the map grid. If false, invert the
 *                   beam grid.
 */
void uvtrans(MapBeam *mb, int domap)
{
  float *image;  /* Pointer to the grid to be inverted */
  float *rxptr;  /* Pointer into rtrans for current X pixel coordinate */
  float *ryptr;  /* Pointer into rtrans for current Y pixel coordinate */
  float *imptr;  /* Pointer into 'image' */
  int ix,iy;     /* Indexes x,y in image */
/*
 * Get a pointer to the map or beam.
 */
  image = domap ? mb->map : mb->beam;
/*
 * Apply phase shifts to make the map centre appear at the centre of
 * the map grid.
 */
  cnj_shift(image, mb->nx, mb->ny);
/*
 * Inverse transform the UV grid.
 */
  newfft(image, mb->nx/2, mb->ny, -1, 1, 0);
/*
 * Multiply the image throughout by the sensitivity function to remove the
 * gridding convolution function.
 */
  imptr = image;
  ryptr = mb->ryft;
  for(iy=0; iy<mb->ny; iy++) {
    rxptr = mb->rxft;
    for(ix=0; ix<mb->nx; ix++)
      *(imptr++) *= *(rxptr++) * *ryptr;
    ryptr++;
  };
  return;
}



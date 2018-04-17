#include "vlbmath.h"
/*.......................................................................
 * Find the min/max value of a 1D array.
 *
 * Input:
 *   vec    float *  Pointer to first element of 1D array of dimension
 *                   vecdim.
 *   vecdim int      The number of float elements vec.
 * Output:
 *   vecmin float *  Pointer to variable into which the min value of the
 *                   array is deposited.
 *   vecmax float *  Pointer to variable into which the max value of the
 *                   array is deposited.
 */
void frange(float *vec, int vecdim, float *vecmin, float *vecmax)
{
  float *vecptr;  /* Pointer into 'vec' */
  float *endptr;  /* Pointer to last element in 'vec' */
  float vmin,vmax;/* Internal versions of vecmin,vecmax */
/*
 * Find the end element of the array.
 */
  endptr = vec+vecdim-1;
/*
 * Initialise mn/max using first element of the array.
 */
  vecptr = vec;
  vmin = vmax = *vecptr;
/*
 * Step through the array searching for the min/max values.
 */
  while(++vecptr <= endptr) {
    if(*vecptr < vmin)
      vmin = *vecptr;
    if(*vecptr > vmax)
      vmax = *vecptr;
  };
/*
 * Return the results.
 */
  *vecmin = vmin;
  *vecmax = vmax;
  return;
}

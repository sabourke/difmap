#include <math.h>
#include "vlbmath.h"
#include "vlbconst.h"

/*.......................................................................
 * Perform a discrete cosine transform between a psuedo continuous
 * interpolation function containing the 1st half of an even function and
 * an output array of the dimension of the grid that the interpolation
 * function corresponds to. The 0 of the output array will be shifted
 * to pixel nout/2. Trigonometric recurrence relations are used to
 * increment cos() in the inner loop. This speeds the algorithm up
 * impressively.
 * 
 * Input:
 *  inparr   float *   The input 1D array of 'ninp' elements.
 *  ninp     int       The number of elements in 'inparr'.
 *  inwid    float     Number of elements that the interpolation
 *                     array corresponds to in the input plane.
 * Output:
 *  outarr   float *   The output 1D array of 'nout' elements.
 *  nout     int       The number of elements in 'outarr'.
 */
void costran(float *inparr, int ninp, float inwid, float *outarr, int nout)
{
  int inp;     /* An element in the input array */
  int out;     /* An element in the output array */
  int icent;   /* Middle element of output array */
  float theta; /* Argument to scale pixel positions in cos(2piux) */
  float *inpptr; /* Pointer into input array */
  float *outptr; /* Pointer into output array */
  double sininc; /* Sine of an angle increment */
  double cosinc; /* Cosine of an angle increment */
  double newcos; /* Current value of cosine */
  double newsin; /* Current value of sine */
  double wtmp;   /* Double precision calculation intermediary */
  icent = nout/2;
/*
 * Determine the 2.pi.u factor for the COS.
 */
  theta = twopi * inwid/ninp/nout;
  outptr = outarr;
  for(out=0; out<=icent; out++) {
    inpptr = inparr;
    *outptr = 0.0f;
    cosinc = cos(theta*(out-icent));
    sininc = sin(theta*(out-icent));
    newcos = 1.0;
    newsin = 0.0;
    for(inp=0; inp<ninp; inp++) {
      *outptr += *(inpptr++) * newcos; /* cos((out-icent)*inp*theta); */
      wtmp = newcos;
      newcos = wtmp*cosinc - newsin*sininc;
      newsin = wtmp*sininc + newsin*cosinc;
    };
    outptr++;
  };
/*
 * Copy the first half of the even transform as calculated above
 * into the second half of the array to complete the transform.
 */
  inpptr = outarr+icent-1;
  outptr = outarr+icent+1;
  for(inp=0; inp<icent-1; inp++)
    *(outptr++) = *(inpptr--);
  return;
}

#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include "utils.h"

extern float pi;
extern int no_error;

#define MAX_ORDER 20

/*.......................................................................
  This routine decomposes the fourier series from data supplied by the
  user. The x and y data arrays are passed to the function along with
  the number of points in them and the period of the fundamental. The
  max_order, amplitudes and phases are returned in arrays asum[] and
  bsum[] (asum[] and bsum[] are also used to accumulate the fourier
  A and B coefficients from which the amplitudes and phases are then
  determined. The data need not be sampled upon a regular grid.
  To ensure orthogonality, the integrations are performed assuming that
  straight lines connect each data point. Thus for instance to integrate
  data(x=XA) to data(x=XB) multiplied by sin(x) the integral of
  (line_connecting_XA_to_XB)*sin(x) is evaluated from XA to XB.
  The output of the routine goes to the user variables AMP(n) and
  PHASE(n), each element of said holding the fourier n'th component.
  The user may then evaluate the Fourier series using these components
  via the user function FVAL(x). (Subroutine FEVAL).
*/
int fourier_series(float x_data[], float y_data[], int npts, float period, float asum[], float bsum[], int max_order)
{
/*
  Local declarations.
*/
	float xa, xb, ya, yb, yend,pgrad,rl, omega,xstart;
	double ph_a, ph_b;
        const float twopi = 6.2831853;
        int last,k,l;
/*
  Check the PERIOD value.
*/
        if(period <= 0.0) {
	  lprintf(stderr, "four: Illegal period given: %f\n", period);
	  return -1;
	};
/*
  Check that the data are in increasing order.
*/
	for(k=1; k<npts; k++) {
	  if(x_data[k] < x_data[k-1]) {
	    lprintf(stderr, "four: The x-array is not in increasing order.\n");
	    return -1;
	  };
	};
/*
  Find the position of the last point within one period of the first point.
  last will actually point at the point following this one.
  If the data don't cover a whole period then last will be the element
  number just off then end of the array.
*/
	for(last=0; last < npts && x_data[last]-x_data[0] < period;)
	  last++;
/*
  The integrals must span a whole period, whereas it is unlikely that
  there is a point with a time of exactly one period from the first
  point. Make the last point have the same value as the first.
*/
	yend = y_data[0];
	xstart = x_data[0];
/*
  Convert the period to a frequency.
*/
        omega = twopi/period;
/*
  Perform the integrations for each harmonic.
*/
	for(l=0; l<max_order; l++) {
/*
  Generate a float version of L.
*/
          rl = (float) l;
/*
  Initialise the A & B coefficient sums to zero.
*/
          asum[l]=0.0;
          bsum[l]=0.0;
/*
  Get the time of the first point, converted into a new form that
  goes between 0 and 2 x pi. Also get the first y-value.
*/
	  xb = 0.0;
	  yb = *y_data;
/*
  Integrate over all of the points up to the end of the required
  period to gain the sum's for the Fourier A(L) & B(L) coefficients.
*/
          for(k=0; k < last;) {
/*
  The second point of the last pair of points is the first point
  of the next pair.
*/
	    xa = xb;
	    ya = yb;
/*
  Get the second (converted) time co-ordinate and y value.
  The second time value must be sufficiently different from
  the first time that an overflow won't occur when the
  gradient of the line connecting the two points is
  evaluated. When they are too close, skip to the next
  point and try again. Use the interpolated/extrapolated
  end point for the last point.
*/
	    while(xb-xa < 1.0e-7 && k < last) {
	      if(++k != last) {
		xb = omega * (x_data[k] - xstart);
		yb = y_data[k];
	      }
	      else {
		xb = twopi;
		yb = yend;
	      };
	    };
/*
  Compute the gradient of the line connecting the current and next
  point (in the above defined coordinate system).
*/
            pgrad = (yb-ya)/(xb-xa);
/*
  For every coefficient except the 0th Fourier coeff evaluate the
  interpolated integral as follows.
*/
            if(l != 0) {
/*
  Pre-generate the phase arguments for the various SIN's and COS's.
*/
	      ph_a = (double) xa * rl;
	      ph_b = (double) xb * rl;
/*
  Increment sum for the Fourier A(L) coefficient.
*/
              asum[l] += (yb * sin(ph_b)  -
			  ya * sin(ph_a)  +
			  pgrad/rl * (cos(ph_b)-cos(ph_a)) )/rl;
/*
  Increment sum for the Fourier B(L) coefficient.
*/
              bsum[l] -= (yb * cos(ph_b)  -
			  ya * cos(ph_a)  -
			  pgrad/rl * (sin(ph_b)-sin(ph_a)) )/rl;
/*
  The zeroth coefficient will be handled in a simpler manner - Simply
  perform a trapezium rule integration between the current points.
*/
	    }
            else {
              asum[0] += (xb-xa) * (ya + (xb-xa) * pgrad/2.0);
	    };
	  };
/*
  The Fourier A & B coefficients are now complete for coeff L, apart
  from a scaling factor of 1/PI.
*/
          asum[l] /= pi;
          bsum[l] /= pi;
/*
  Evaluate the new harmonic contribution and subtract it from the
  previous residual in the work array. The integration for the
  next harmonic will be performed upon this residual. This reduces
  precision related errors.
          if(l == 0) {
	    tmp = asum[0]/2.0;
	    for(k=0; k<last; k++)
              y_data[k] -= tmp;
	    yend -= tmp;
	  }
          else {
	    for(k=0; k<last; k++) {
              tmp = rl * omega * (x_data[k]-x_data[0]);
	      y_data[k] -= asum[l] * cos(tmp) - bsum[l] * sin(tmp);
	    };
	    yend -= asum[l];
	  };
*/
/*
  If any more harmonics have been requested then go back and
  integrate to find the contribution of the next harmonic.
*/
	};
/*
  The decomposition has been made so now turn the A and B coefficients
  into amplitudes and phases.  Place the zero frequency term in the
  0'th coefficient.  The phases are in the time units of the input
  array relative to the start time and actually specify the time
  of the first zero on the y-axis of the a sine wave of the given
  harmonic. Obviously this makes no sense in terms of the zero frequency
  term so set its phase equal to start time of the data.
*/
        asum[0] /= 2.0;
        bsum[0] = xstart;
/*
  For each harmonic convert the cosine+sine representation into
  amplitude versus phase representation, in the same units as
  the input data.
*/
	for(k=1; k < max_order; k++) {
/*
  Use xa and xb to hold temporary copies of asum and bsum.
*/
	  xa = asum[k];
	  xb = bsum[k];
/*
  Start with the amplitude.
*/
          asum[k] = sqrt(xb * xb + xa * xa);
/*
  Then get the phase, in the same units as the input period,
  ie the first zero phase time following the start time of the data.
*/
          if(xa != 0 || xb != 0) {
            bsum[k] = atan2(xb,xa) * period/(k * twopi) + xstart;
	  }
          else {
            bsum[k] = 0;
	  };
	};
	return no_error;
}

/*.......................................................................
  Evaluate the fourier series fit of a past invokation of subroutine
  FOUR. It is assumed that the amplitudes and phases of each order
  of the fit are held in the user variables AMPLITUDE & PHASE.
  The current values of the user variables, PERIOD,ORDER,
  FILTER and FTYPE are also used. Parsed to the routine are the time
  value XXVAL at which to evaluate the series, and TSHFT to specify
  any shift between the required model and the original fit. YYVAL is
  returned containing the series value at XXVAL (shifted be TSHFT).
  ORDER specifies the maximum harmonic to be included. PERIOD specifies
  the fourier period, as in FOUR. The fourier component amplitudes
  used are multiplied by the user variable FILTER(n), such that the
  user can low-pass high-pass, etc.. filter the series. FTYPE
  determines what type of curve to be returned. For values, 0 <=-1,
  >=1, the series value, the integral, or the first derivatives are
  returned respectively.
*/
int fourier_series_value(float xval, float *yval, int differential_order, float period, float amp[], float phase[], float filter[], int max_order)
{
/*
  Local declarations.
*/
        static float omega,tmp,harm,ampl;
	static int l;
        const float twopi = 6.2831853;
/*
  Check the period value.
*/
        if(period <= 0.0) {
	  lprintf(stderr, "fsval: Illegal period given: %f\n", period);
	  return -1;
	};
/*
  Evaluate the fourier series on the input time grid.
*/
        omega = twopi/period;
/*
  Start with the contribution from the DC (freq=0) term.
  Get the relevant value depending upon the diffential order.
  -ve orders mean integrate.
*/
	switch (differential_order) {
/*
  0'th order. DC term contributes directly.
*/
	case 0:
          *yval = amp[0] * filter[0];
	  break;
/*
  For derivatives the zero frequency term makes no contribution.
*/
	case 1: case 2:
          *yval = 0.0;
	  break;
/*
  In the integral the DC term integrates to amplitude * time.
*/
	case -1:
          *yval = amp[0] * filter[0] * (xval-phase[0]);
	  break;
	default:
	  lprintf(stderr, "fsval: Unsupported differential order: %d\n", differential_order);
	  return -1;
	  break;
	};
/*
  Add in the subsequent harmonics one at a time.
*/
	for(l=1; l<max_order; l++) {
/*
  Apply the user harmonic filter.
*/
          ampl = amp[l] * filter[l];
/*
  Calculate the angular frequency of the current harmonic.
*/
          harm = l * omega;
/*
  Calculate the (nwt) term for the current harmonic.
*/
          tmp = harm * (xval-phase[l]);
/*
  Add in the contribution from the current harmonic, L at the
  time XVAL. This requires different expressions for
  different derivatives.
*/
	  switch (differential_order) {
	  case 0:
            *yval += ampl * cos(tmp);
	    break;
/*
  First derivative.
*/
	  case 1:
            *yval -= ampl * harm * sin(tmp);
	    break;
/*
  Second derivative.
*/
	  case 2:
            *yval -= ampl * harm * harm * cos(tmp);
	    break;
/*
  Integral.
*/
	  case -1:
            *yval += ampl * sin(tmp)/harm;
	    break;
	  };
	};
	return no_error;
}

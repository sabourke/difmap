#include <math.h>

#include "besj.h"

/*.......................................................................
 * Return the bessel function Jo(x). Translated to C from FORTRAN
 * numerical recipes.
 *
 * Input:
 *  x       float   The value to compute Jo at.
 * Output:
 *  return  float   The value Jo(x).
 */
float c_besj0(float x)
{
/*
 * Declare bessel constants.
 */
  static double p1 =  1.0;
  static double p2 = -1.098628627e-3;
  static double p3 =  2.734510407e-5;
  static double p4 = -2.073370639e-6;
  static double p5 =  2.093887211e-7;

  static double q1 = -1.562499995e-2;
  static double q2 =  1.430488765e-4;
  static double q3 = -6.911147651e-6;
  static double q4 =  7.621095161e-7;
  static double q5 = -9.34945152e-8;

  static double r1 =  57568490574.0;
  static double r2 = -13362590354.0;
  static double r3 =  651619640.7;
  static double r4 = -11214424.18;
  static double r5 =  77392.33017;
  static double r6 = -184.9052456;

  static double s1 = 57568490411.0;
  static double s2 = 1029532985.0;
  static double s3 = 9494680.718;
  static double s4 = 59272.64853;
  static double s5 = 267.8532712;
  static double s6 = 1.0;
  double y,xx,z;
  double abs_x;    /* Absolute version of x */
/*
 * Use rational-function approximation if x is small.
 */
  abs_x = fabs(x);
  if(abs_x < 8.0) {
    y = x * x;
    return (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) /
           (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));
  } else {
    z = 8.0/abs_x;
    y = z * z;
    xx = abs_x-0.785398164;
    return sqrt(0.636619772/abs_x) *
      ( cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5)))) - 
      z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))) );
  };
}

/*.......................................................................
 * Return the bessel function J1(x). Translated to C from FORTRAN
 * numerical recipes.
 *
 * Input:
 *  x       float   The value to compute J1 at.
 * Output:
 *  return  float   The value J1(x).
 */
float c_besj1(float x)
{
/*
 * Declare bessel constants.
 */
  static double p1 =  1.0;
  static double p2 =  1.83105e-3;
  static double p3 = -3.516396496e-5;
  static double p4 =  2.457520174e-6;
  static double p5 = -2.40337019e-7;

  static double q1 =  0.04687499995;
  static double q2 = -2.002690873e-4;
  static double q3 =  8.449199096e-6;
  static double q4 = -8.8228987e-7;
  static double q5 =  1.05787412e-7;

  static double r1 =  72362614232.0;
  static double r2 = -7895059235.0;
  static double r3 =  242396853.1; 
  static double r4 = -2972611.439;
  static double r5 =  15704.48260;
  static double r6 = -30.16036606;

  static double s1 = 144725228442.0;
  static double s2 = 2300535178.0;
  static double s3 = 18583304.74;
  static double s4 = 99447.43394;
  static double s5 = 376.9991397;
  static double s6 = 1.0;
  double y,xx,z;
  double abs_x;    /* Absolute version of x */
/*
 * Use rational-function approximation if x is small.
 */
  abs_x = fabs(x);
  if(abs_x < 8.0) {
    y = x * x;
    return x * (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) /
               (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));
  } else {
    z = 8.0/abs_x;
    y = z * z;
    xx = abs_x - 2.356194491;
    return (x>=0.0 ? 1.0:-1.0) * sqrt(0.636619772/abs_x) *
      ( cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5)))) - 
      z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))) );
  };
}

/*.......................................................................
 * Return the bessel function J2(x). Based on the numerical recipes
 * function BESSJ().
 *
 * Input:
 *  x       float   The value to compute J2 at.
 * Output:
 *  return  float   The value J2(x).
 */
float c_besj2(float x)
{
  float retval;   /* The return value of the function */
/*
 * The function is even, so substitute the absolute value of x.
 */
  if(x < 0.0f)
    x = -x;
/*
 * Special case for x=0.
 */
  if(x==0.0f)
    return 0.0f;
/*
 * If x exceeds the Bessel function order n, (ie. 2), use the
 * recurrence relation:
 *  Jn+1(x) = 2*n*Jn(x)/x - Jn-1(x).
 */
  if(x > 2.0f) {
    retval = 2.0 * c_besj1(x)/x - c_besj0(x);
  }
/*
 * The recurrence relation is unstable for values of x <= the order.
 * In such cases revert to downward recurrence using Miller's algorithm.
 * To do this choose arbitrary values for the start even and odd
 * bessel functions for a given starting order (n>>2), accumulate the
 * sum J0 + 2(J2 + J4 + J6 etc...) which is known to aymptotically
 * approach 1.0, preserve the arbirarily scaled J2 value when the
 * downward recurrence reaches it, and then when J0 has been reached
 * re-normalize it with the constant needed to make the above sum
 * evaluate to 1.0.
 */
  else {
    const float large=1.0e10f; /* Trigger value for overflow prevention */
    int start = 10;    /* Start order from Numerical Recipes for order=2 */
    float bjpp=0.0f;   /* Jn+2 (previous previous Jn) */
    float bjp=1.0f;    /* Jn+1 (previous Jn) */
    float bj;          /* Jn */
    float normsum=0.0f;/* J0 + 2(J2 + J4 + J6 etc...) sum for renormalization */
    int order;         /* The order being computed */
/*
 * Start the reverse recurrence: Jn(x) = 2*(n+1)*Jn+1(x)/x - Jn+2(x).
 */
    float recfac = 2.0f / x;  /* Jn(x) = recfac*(n+1)*Jn+1(x) - Jn+2(x). */
    retval = 0.0;
    for(order=start; order>1; order--) {
      bj = recfac * (order+1) * bjp - bjpp;
/*
 * Shift queue of bj values.
 */
      bjpp = bjp;
      bjp = bj;
/*
 * The scale factor in the recurrence is arbitrary so re-normalize
 * if the values get too large.
 */
      if(bj > large) {
	bj /= large;
	bjp /= large;
	bjpp /= large;
	retval /= large;
      };
/*
 * Accumulate the normalization sum: J0 + 2(J2 + J4 + J6 etc...).
 */
      if(order % 2 == 0)
	normsum += (order ? 2.0 : 1.0) * bj;
/*
 * Record the required order when it is computed.
 */
      if(order==2)
	retval = bj;
    };
/*
 * Re-normalize the return value, from the knowledge that normsum if
 * properly scaled would equal 1.0.
 */
    retval /= normsum;
  };
/*
 * Return the bessel function value.
 */
  return retval;
}

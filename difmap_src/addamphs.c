#include <math.h>
#include <stdio.h>
#include "vlbutil.h"

/*.......................................................................
 * Given two amplitude,phase pairs, add them in the complex plane
 * and return the result as amplitude and phase.
 *
 * Input/Output:
 *  amp   float *   The amplitude of the coordinate pair to be added to.
 *  phs   float *   The phase of the coordinate pair to be added to.
 * Input:
 *  addamp float    The amplitude of the coordinate pair to be added.
 *  addphs float    The phase of the coordinate pair to be added.
 */
void addamphs(float *amp, float *phs, float addamp, float addphs)
{
  float x;  /* Real part of complex sum */
  float y;  /* Imaginary part of complex sum */
/*
 * Add the complex versions of the polar numbers.
 */
  x = *amp * cos(*phs) + addamp * cos(addphs);
  y = *amp * sin(*phs) + addamp * sin(addphs);
/*
 * Convert back to amplitude and phase.
 */
  *amp = sqrt(x*x+y*y);
  *phs = (x==0.0f && y==0.0f) ? 0.0f : atan2(y,x);
  return;
}

/*.......................................................................
 * Add a complex value expressed as real and imaginary parts to another
 * complex lvalue expressed in amplitude and phase.
 *
 * Input/Output:
 *  amp    float *  The amplitude of the coordinate pair to be added to.
 *  phs    float *  The phase of the coordinate pair to be added to.
 * Input:
 *  re,im  float    The real and imaginary parts to be added.
 */
void add_cart_to_polar(float *amp, float *phs, float re, float im)
{
  float x;  /* Real part of complex sum */
  float y;  /* Imaginary part of complex sum */
/*
 * Convert the amplitude and phase to real,imaginary to perform the
 * addition.
 */
  x = *amp * cos(*phs) + re;
  y = *amp * sin(*phs) + im;
/*
 * Convert back to amplitude and phase.
 */
  *amp = sqrt(x*x+y*y);
  *phs = (x==0.0f && y==0.0f) ? 0.0f : atan2(y,x);
  return;
}

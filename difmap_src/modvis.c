#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "modvis.h"
#include "model.h"
#include "besj.h"
#include "vlbconst.h"
#include "obs.h"
#include "logio.h"

/*
 * The following declarations implement a cosine lookup table.
 */
#define CTSIZ 4096                 /* Size of lookup table */

static float costab[CTSIZ+2];      /* Cosine lookup table */
static float *cos_tab=&costab[1];  /* Allow elements -1 to CTSIZ */
static int ready=0;                /* True after table is initialized */
static const int soffset=CTSIZ+(CTSIZ/4); /* Index to cos(360+90 degrees) */

static void add_cmpvis_to_modvis(float amp, float phs, float *re, float *im);

/*.......................................................................
 * Compute the visibility amplitude and phase at a given U,V coordinate
 * corresponding to a single model component.
 *
 * Input:
 *  cmp    Modcmp *  The descriptor of the model component.
 *  sub  Subarray *  The parent subarray of the model visibility.
 *  base      int    The index of the parent baseline.
 *  freq    float    The frequency for which to compute the component (Hz).
 *  uu      float    The U coordinate in the UV plane (wavelengths).
 *  vv      float    The V coordinate in the UV plane (wavelengths).
 * Output:
 *  amp     float *  The corresponding amplitude.
 *  phs     float *  The corresponding phase (radians).
 */
void cmpvis(Modcmp *cmp, Subarray *sub, int base, float freq, float uu,
	    float vv, float *amp, float *phs)
{
  float cmpamp;  /* The amplitude of the component */
/*
 * Compute the primary beam factor.
 */
  float pb = pb_bl_factor(sub, base, freq,
			  calc_pointing_offset(sub->ob, cmp->x, cmp->y));
/*
 * Compute the spectral index factor.
 */
  float spec = cmp->spcind==0.0 ? 1.0 : pow(freq/cmp->freq0, cmp->spcind);
/*
 * Compute the flux of the component.
 */
  float flux = cmp->flux * spec * pb;
/*
 * All the components are even functions so the phase is just the
 * fourier component phase at the centroid of the component.
 */
  float cmpphs = (float) twopi * (uu*cmp->x + vv*cmp->y);
/*
 * The component amplitude is different for each component type.
 * Delta components are the easiest and commonest.
 */
  if(cmp->type == M_DELT) {
    cmpamp = flux;
  }
/*
 * The other types are a little more complicated!
 */
  else {
/*
 * Pre-compute parameters.
 */
    double sinphi = sin(cmp->phi);
    double cosphi = cos(cmp->phi);
    double tmpa = (vv*cosphi+uu*sinphi);
    double tmpb = (cmp->ratio*(uu*cosphi-vv*sinphi));
    double tmpc = pi * cmp->major * sqrt(tmpa*tmpa + tmpb*tmpb);
/*
 * Limit tmpc to sensible values to prevent underflow,overflow and
 * divide-by-zero errors.
 */
    if(tmpc < 1.0e-9)
      tmpc = 1.0e-9;
/*
 * See the "Introduction to Caltech VLBI programs" for details of
 * the type-specific calculations below.
 */
    switch(cmp->type) {
    case M_GAUS:
      cmpamp = flux * (tmpc<12.0 ? exp(-0.3606737602 * tmpc*tmpc):0.0);
      break;
    case M_DISK:
      cmpamp = 2.0f * flux * c_besj1(tmpc)/tmpc;
      break;
    case M_ELLI:
      cmpamp = 3.0f * flux * (sin(tmpc)-tmpc*cos(tmpc))/(tmpc*tmpc*tmpc);
      break;
    case M_RING:
      cmpamp = flux * c_besj0(tmpc);
      break;
    case M_RECT:
      tmpa = pi * cmp->major * (uu*sinphi+vv*cosphi);
      cmpamp = flux * (fabs(tmpa) > 0.001 ? sin(tmpa)/tmpa : 1.0);
      break;
    case M_SZ:
      cmpamp = flux * (tmpc < 50.0 ? exp(-tmpc) : 0.0) / tmpc;
      break;
    default:
      lprintf(stderr, "Ignoring unknown model component type: %d\n", cmp->type);
      cmpamp = 0.0f;
      cmpphs = 0.0f;
    };
  };
/*
 * Assign the return values.
 */
  *amp = cmpamp;
  *phs = cmpphs;
  return;
}

/*.......................................................................
 * Add the contribution of a model component to a given model visibility.
 *
 * Input:
 *  cmp     Modcmp *   The component who's contribution is to be added.
 *  sub   Subarray *   The parent subarray of the model visibility.
 *  base       int     The index of the parent baseline.
 *  freq     float     The frequency for which to compute the component (Hz).
 *  uu       float     The U coordinate in the UV plane (wavelengths).
 *  vv       float     The V coordinate in the UV plane (wavelengths).
 * Input/Output:
 *  re, im   float *   The real and imaginary parts of the model component
 *                     visibility will be added to *re and *im, respectively.
 */
void add_cmp_to_modvis(Modcmp *cmp, Subarray *sub, int base, float freq,
		       float uu, float vv, float *re, float *im)
{
  float amp, phs;   /* The amplitude and phase of the component visibility */
/*
 * Compute the amplitude and phase of the model component visibility.
 */
  cmpvis(cmp, sub, base, freq, uu, vv, &amp, &phs);
/*
 * Add the amplitude and phase to the model visibility.
 */
  add_cmpvis_to_modvis(amp, phs, re, im);
  return;
}

/*.......................................................................
 * This is a private function of add_cmp_to_modvis(). It converts the
 * specified amplitude and phase to the correspondong real and
 * imaginary parts, and adds these to the given model imaginary and
 * phase.
 *
 * Input:
 *   amp     float    The amplitude of the model component visibility.
 *   phs     float    The phase of the model component visibility.
 * Input/Output:
 *  re, im   float *  The real and imaginary parts of the model visibility
 *                    to which to add those of the new component.
 */
static void add_cmpvis_to_modvis(float amp, float phs, float *re, float *im)
{
/*
 * Divide the phase by 2*pi.
 */
  float off = phs / twopi;
/*
 * Get the sign of this phase.
 */
  int isign = (off<0.0f)?-1:1;
/*
 * Work out the index of the normalized phase in the cosine lookup table.
 */
  float ftmp = off*(isign*CTSIZ);
/*
 * Convert this to an integer.
 */
  int cos_indx = (int) ftmp;
/*
 * Get the non-integral part of the index.
 */
  float err_indx = ftmp - cos_indx;
/*
 * Fold the table index to the size of the table, then get a pointer
 * to the corresponding table value.
 */
  float *cos_ptr = &cos_tab[cos_indx %= CTSIZ]; /* Element of folded index */
/*
 * Initialise the cosine lookup table on first entry to this function.
 */
  if(!ready) {
    int itab;
    for(itab=-1; itab<CTSIZ+1; itab++)
      cos_tab[itab]=cos(twopi*itab/CTSIZ);
    ready = 1;
  };
/*
 * Iterpolate between cos_ptr[0] and cos_ptr[1] to determine cos(twopi*off),
 * and use this to compute the real part of the model visibility. Add this
 * to the overall model visibility.
 */
  *re += amp * (*cos_ptr + err_indx*(cos_ptr[1] - *cos_ptr));
/*
 * Now work out the integer part of the index of sin(twopi*off) in
 * the cosine table, noting that 'soffset' offsets the index by
 * 360+(90-angle).
 */
  cos_indx = (soffset - isign*cos_indx) % CTSIZ;
/*
 * Get a pointer to the corresponding table entry.
 */
  cos_ptr  = &cos_tab[cos_indx];
/*
 * Iterpolate between cos_ptr[0] and cos_ptr[-isign] to determine
 * cos(twopi*off), and use this to compute the imaginary part of the model
 * visibility. Add this to the overall model visibility.
 */
  *im += amp * (*cos_ptr + err_indx*(cos_ptr[-isign] - *cos_ptr));
  return;
}

/*
 * Copyright (c) 1997 by the California Institute of Technology.
 * All rights reserved.
 *
 * Original author: Martin Shepherd  (mcs@astro.caltech.edu)
 *
 * This file is copyrighted but available free for non-commercial
 * purposes under the understanding that neither the author, nor
 * Caltech shall be held responsible for damages resulting from
 * its use. It can also be modified under the understanding that
 * this copyright notice not be removed, and with the proviso
 * that the author of subsequent modifications appends note to this
 * notice to say that the file has been modified by them.
 */

#include <math.h>
#include <stdio.h>
#include "vlbconst.h"
#include "vlbfft.h"

static void bitswap(float *image, int curdim, int curinc, int othdim,
		    int othinc);
static void fixreal(float *image, int adim, int bdim, int isign);
static void opn_fft(float *image, int adim, int bdim);
static void cls_fft(float *image, int adim, int bdim);

/*.......................................................................
 * Perform a 2D FFT.
 *
 * Input:
 *  image  float *      The image to be transformed.
 *
 *                      For Complex transforms, this must be a 1D
 *                      array of length 2*adim*bdim, where the 2
 *                      reflects the fact that real and imaginary
 *                      parts are stored in consecutive float's. The
 *                      'image' must be organised as 'bdim'
 *                      consecutive vectors of 'adim' elements. 
 *
 *                      For Real <=> Conjugate-Symmetric transforms
 *                      image[] must be an array of 2*(adim+1)*bdim
 *                      floats. For the forward transform, real
 *                      data are packed into the first 2*adim*bdim
 *                      elements. The resulting transform yields one
 *                      half of a conjugate symmetric array of
 *                      (adim+1)*bdim complex numbers.
 *                      
 *  adim     int        The number of complex points in a row of
 *                      'image'. This must be a power of 2.
 *  bdim     int        The number of rows 'image'. This must be a power
 *                      of 2. (bdim may be 1).
 *  
 *  isign    int        +1 => do the forward transform.
 *                      -1 => do the inverse transform. The inverse
 *                            transform will be correctly scaled.
 *  isreal   int        For a real transform to a half conjugate symmetric
 *                      array, set isreal=1.
 *  rescale  int        If true, scale the transform by 1/(adim*bdim).
 *                      This should be done either on forward transforms
 *                      or on the inverse transforms, but not both.
 *                      Which one you choose is largely a matter of
 *                      convention and the nature of your data.
 */
void newfft(float *image, int adim, int bdim, int isign, int isreal,
	    int rescale)
{
  int ini_dim; /* Number of first dimension to be transformed */
  int end_dim; /* Number of last dimension to be transformed */
  int inc_dim; /* Increment between ini_dim and end_dim */
  int curdim=0;/* Number of complex elements on current axis */
  int othdim=0;/* Number of elements on the other of the two axes */
  int curinc=0;/* Number of floats between complex elements on current axis */
  int othinc=0;/* Number of floats between complex elements on other axis */
  int slot;    /* An element number on current axis */
  int ipos;    /* Index into row or column on the "other" axis */
  float *vecptr;/* Pointer to first element of a row or column vector */
  float *kptr;  /* Pointer to k'th element within a sub-transform */
  float *erptr; /* Pointer to the real part of an element in even tranform */
  float *eiptr; /* Pointer to the imag part of an element in even transform */
  float *orptr; /* Pointer to the real part of an element in odd tranform */
  float *oiptr; /* Pointer to the imag part of an element in odd transform */
  double wr;    /* Real part of W(k) = exp(2.pi.i.k/N) */
  double wi;    /* Imag part of W(k) = exp(2.pi.i.k/N) */
  double omega; /* +/-2.pi/N - the increment angle to get W(k+1) from W(k) */
  double sininc;/* sin(omega) for use in recurrence relations */
  double cosinc;/* cos(omega) for use in recurrence relations */
  double wtmp;  /* Intermediary for double precision calculations */
  float fwr,fwi;/* Float versions of wr and wi. */
  float wreal;  /* Real part of W(k)*F(k_odd) */
  float wimag;  /* Imaginary part of W(k)*F(k_odd) */
  int ntrans;   /* Number of points in transform being built */
  int nprev;    /* Width of even/odd transforms to be combined to ntrans */
  int traninc;  /* Increment from one 'ntrans' transform to next */
  int previnc;  /* Increment from one half-transform to next */
  int k;        /* k'th element in a transform */
  int i;
/*
 * Loop over dimensions for the bit-reversal phase. When the input array
 * is a real-array being forward transformed then the first dimension
 * (rows) must be tranformed first, whereas where the input array is
 * a half-conjugate-symmetric array being inverse transformed, the
 * first dimension must be tranformed last. For normal tranforms the
 * order is irrelavant.
 */
  if(isign == 1) {
    ini_dim=0;
    end_dim=2;
    inc_dim=1;
  } else {
    ini_dim=1;
    end_dim=-1;
    inc_dim=-1;
  };
  for(i=ini_dim; i!=end_dim; i += inc_dim) {
/*
 * Set up element increments and dimensions for the current axis.
 * Also handle conversion of the rows of real/conjugate-symmetric
 * transform pairs - this must be done after all other dimensions
 * have been transformed and for this reason the 1st dimension of
 * the image will be transformed last.
 */
    switch(i) {
    case 0:
      curdim = adim;   /* Number of elements on current axis (ncol) */
      othdim = bdim;   /* Number of elements on other axis (nrow) */
      curinc = 2;      /* Two floats between elements on current axis */
      othinc = 2*adim; /* Number of floats between elements on other axis */
      if(isreal && isign == -1)
	fixreal(image,adim,bdim,isign);
      break;
    case 1:
      curdim = bdim;
      othdim = adim;
      curinc = 2*adim;
      othinc = 2;
/*
 * Real tranforms actually have an extra complex column to be
 * transformed.
 */
      if(isreal) {
	if(isign == 1)
	  fixreal(image,adim,bdim,isign);
	othdim = adim+1;
	curinc = 2*(adim+1);
      };
      break;
    };
/*
 * Swap columns/rows in which the integer array indexes are mirrored
 * in their bit representation.
 */
    bitswap(image, curdim, curinc, othdim, othinc);
/*
 * Now combine all the one-point transformed that resulted from the
 * above swapping, first into 2-point transforms, then 4-point transforms
 * etc... 'ntrans' is the number of complex points per transform at
 * each step.
 */
    vecptr = image;
    for(ipos=0; ipos < othdim; ipos++) {
      for(ntrans=2; ntrans <= curdim; ntrans<<=1) {
	nprev = ntrans/2;  /* Length of odd and even parts being combined */
/*
 * Determine the increment required to step over one completed transform
 * of 'ntrans' points. Also the increment over one of the sub-transforms.
 */
	traninc = ntrans * curinc;
	previnc = nprev * curinc;
/*
 * Pre-calculate part of the argument of exp(2.pi.i/ntrans) as
 * required in the Danielson and Lanczos formula. The 'isign'
 * term decides if the transform is forward or reverse.
 */
	omega = isign*twopi/ntrans;
/*
 * wr and wi are the real and imaginary parts of exp(i.omega.k).
 * Where 'k' is the element of the transform being calculated.
 * For k=0 they are 1 and 0 respectively.
 */
	fwr = (float) (wr = 1.0);
	fwi = (float) (wi = 0.0);
/*
 * As k increments by 1 within each iteration of the following loop,
 * we can simply calculate wr and wi by simple trig angle addition
 * equations. For these we will require cos() and sin() of omega.
 */
	cosinc = cos(omega);
	sininc = sin(omega);
/*
* Combine the adjacent even and odd transforms. NB. Both transforms
* are period in 1/ntrans which is why we can store each in just
* 1/ntrans elements. Since we need each element of the two transforms
* twice, once for k=n (n=0->ntrans/2-1) and then again for k=n+ntrans/2,
* The calculation for these two elements will be performed in one go.
* In order to do this, note that if W(k)=exp(2.pi.i.k/ntrans) then
* W(k+ntrans/2) = -W(k).
*/
	kptr=vecptr;  /* Start of current row or column in image */
	for(k=0; k<nprev; k++) {
/*
 * Loop for each tranform.
 */
	  erptr = kptr;         /* Pointer to element k of even transform */
	  orptr = kptr+previnc; /* Pointer to element k of odd transform */
	  for(slot=0; slot<curdim; slot += ntrans) {
	    eiptr = erptr+1;
	    oiptr = orptr+1;
	    wreal = *orptr*fwr - *oiptr*fwi; /* W(k)*F(k_odd) - real */
	    wimag = *orptr*fwi + *oiptr*fwr; /* W(k)*F(k_odd) - imaginary */
/*
 * First evaluate F(kp) where kp=k+nprev and W(kp) is -W(k).
 * Thus F(kp)=F(k_even)-W(k)*F(k_odd).
 */
	    *orptr = *erptr - wreal;
	    *oiptr = *eiptr - wimag;
/*
 * Now evaluate F(k)=F(k_even)+W(k)*F(k_odd).
 */
	    *erptr += wreal;
	    *eiptr += wimag;
/*
 * Move to element k of the next even and odd transforms.
 */
	    erptr += traninc;
	    orptr += traninc;
	  };
/*
 * Since in the next iteration k goes up by one, compute sin() and
 * cos() of omega plus the current angle.
 */
	  wtmp = wr;
	  wr = wtmp*cosinc - wi*sininc;
	  wi = wtmp*sininc + wi*cosinc;
	  fwr = (float) wr;  /* Single precision is sufficient */
	  fwi = (float) wi;  /* outside the recurrence */
/*
 * Next element of target 'ntrans' transform.
 */
	  kptr += curinc;
	};
      };
/*
 * Next row or column of image to be transformed.
 */
      vecptr += othinc;
    };
  };
/*
 * On the reverse transform divide throughout by the number of complex data
 * points.
 */
  if(rescale) {
    int nfloat=(isreal ? (adim+1):adim)*bdim*2;
    float ncomplex = (isreal ? adim*2 : adim) * bdim;
    for(i=0; i<nfloat; i++)
      image[i] /= ncomplex;
  };
  return;
}

/*.......................................................................
 * Given an input image of dimensions (adim,bdim) apply the phase shift
 * specified in the shift theorem, necessary to move the centre of the
 * fourier transform of the input image either from the array centre
 * to 0,0 or from 0,0 to the array centre. The array centre is taken as
 * being element adim/2,bdim/2. The phase shift for an adim/2,bdim/2
 * coordinate shift is exp(-pi*i*(adim*u+bdim*v)), which is simply
 * equal to -1 when (adim*u+bdim*v) is odd and +1 when it is even.
 *
 * Input/Output:
 *  image  float *   A 1D array with adim*bdim*2 elements, where adim
 *                   and bdim denote the dimensions of the equivalent
 *                   2D array and must be integral powers of two and
 *                   the factor of 2 indicates that each element of this
 *                   array is made up of two consecutive floats, the
 *                   first is the real part and the second, the imaginary
 *                   part.
 * Input:
 *  adim   int       The first dimension (fastest changing index) of the
 *                   2D array equivalent of 'image'.
 *  bdim   int       The second dimension of the 2D array equivalent of
 *                   'image'.
 */
void fft_shift(float *image, int adim, int bdim)
{
  int ia;            /* Step count along 1st dimension */
  int ib;            /* Step count along 2nd dimension */
  float *fptr=image; /* Pointer to next element to be scaled by -1 */
/*
 * Loop over pairs of even and odd rows.
 */
  for(ib=0; ib<bdim; ib+=2) {
/*
 * Even rows.
 */
    fptr += 2;
    for(ia=0; ia<adim; ia+=2, fptr += 4) {
      *fptr *= -1.0f;     /* Real part */
      *(fptr+1) *= -1.0f; /* Imaginary part */
    };
/*
 * The odd row following the above even row.
 */
    fptr -= 2;
    for(ia=0; ia<adim; ia+=2, fptr += 4) {
      *fptr *= -1.0f;
      *(fptr+1) *= -1.0f;
    };
  };
  return;
}

/*.......................................................................
 * Swap rows or columns etc of an FFT array where the index of each
 * swapped row or column is the bit-reversed value of the other.
 *
 * Input:
 *   float *image;    Pointer to image being FFT'd.
 *   int curdim;      The number of elements along the current axis.
 *   int curinc;      The increment (in floats) between elements on the
 *                    current axis.
 *   int othdim;      The number of elements on the perpendicular 'other'
 *                    axis.
 *   int othinc;      The increment (in floats) between elements of the other
 *                    axis.
 *
 */
static void bitswap(float *image, int curdim, int curinc, int othdim,
		    int othinc)
{
  static int slot;    /* An element number on current axis */
  static int ipos;    /* An element number on other axis axis */
  static int idim;    /* Temporary copy of curdim */
  static int swapd;   /* The index of one of two rows/columns to be swapped */
  static int orig;    /* The index of one of two rows/columns to be swapped */
  static float ftmp;  /* Intermediary for swapping element */
  static float *aptr; /* Pointer to be stepped along a column or row */
  static float *bptr; /* Pointer to be stepped along a column or row */
/*
 * Swap each vector that is perpendicular to the current dimension at the
 * current element (slot) on the current dimension.
 */
  for(slot=0; slot<curdim; slot++) {
/*
 * Find the bit-reversed value of 'slot'. Build it in 'swapd'.
 */
      idim = curdim;
      orig = slot; swapd = 0;
      while(idim>>=1) {
	swapd <<= 1;        /* Make room for the new bit */
	swapd |= orig & 01; /* Add the new bit to the vacated slot */
	orig  >>= 1;        /* Demote the next msb of 'slot' to be its lsb */
      };
/*
 * Clause to prevent columns/rows being swapped twice or to themselves.
 * If it is met then swap column/row 'slot' for 'swapd'.
 */
      if(swapd < slot) {
	aptr = image+slot*curinc;   /* Pointer to start of 1st column/row */
	bptr = image+swapd*curinc;  /* Pointer to start of 2nd column/row */
	for(ipos=0; ipos < othdim; ipos++) {
	  ftmp  = *aptr;          /* Swap the real part */
	  *aptr = *bptr;
	  *bptr = ftmp;
	  ftmp    = aptr[1];      /* Swap the imaginary part */
	  aptr[1] = bptr[1];
	  bptr[1] = ftmp;
	  aptr += othinc;         /* Go to next element in each column/row */
	  bptr += othinc;
	};
      };
    };
    return;
}

/*.......................................................................
 * This routine is provided for:
 * 1) forward transforming real 2D images of
 *    dimensions 2*adim*bdim floats to half its conjugate symmetric
 *    tranform of (adim+1)*bdim complex numbers (pairs of floats).
 * 2) performing the inverse of the above.
 * The array must have space for 2*(adim+1)*bdim floats.
 *
 * Input/Output:
 *  image   float *   The image of dimensions adim*2 by bdim floats.
 *                    And half of the conjugate symmetric tranform in
 *                    (adim+1)*bdim complex float pairs.
 * Input:
 *  adim    int       See under image. adim must be a power of 2.
 *  bdim    int       See under image. bdim must be a power of 2.
 *  isign   int        1 - forward transform real image to complex.
 *                    -1 - inverse tranform complex to real image.
 */
static void fixreal(float *image, int adim, int bdim, int isign)
{
  int irow;        /* Row number */
  int icol;        /* Column number */
  int rowinc;      /* Number of floats between one row and next */
  float *rowptr;   /* Pointer to beiginning of a row */
  float *cnjptr;   /* Pointer to pixel symmetrically opposite to rowptr */
  double theta;    /* Increment angle for complex exponential */
  double sininc;   /* sin(theta) */
  double cosinc;   /* cos(theta) */
  double wr,wi;    /* Real/imag values of exp(i.n.theta) */
  float fwr,fwi;   /* Single precision versions of fwr,fwi */
  float rsum_a;    /* Intermediate summation (real part) */
  float isum_a;    /* Intermediate summation (imag part) */
  float rsum_b;    /* Intermediate summation (real part) */
  float isum_b;    /* Intermediate summation (imag part) */
  float *rnptr;    /* Pointer to real part of complex_row[n] */
  float *inptr;    /* Pointer to imag part of complex_row[n] */
  float *rmptr;    /* Pointer to real part of complex_row[adim-n] */
  float *imptr;    /* Pointer to imag part of complex_row[adim-n] */
  double wtmp;     /* Double precision calculation intermediary */
  float scal;      /* Transform direction-sensitive scale factor */
  float norm;      /* Normalization factor to account for half-size array */
/*
 * Forward transformed real array.
 * Make way for the two floats required at x=N/2.
 */
  if(isign == 1)
    opn_fft(image, adim, bdim);
/*
 * One of the algorithm's scale factors is -0.5 for forward
 * transforms and +0.5 on reverse transforms. Precompute this.
 * Also compute the re-normalization factor needed to account for
 * the fact that we are simulating an array of twice the length.
 */
  if(isign == 1) {
    scal = -0.5f;
    norm = 1.0f;
  } else {
    scal = 0.5f;
    norm = 2.0f;
  };
/*
 * We need sin() and cos() of pi/adim for calculating successive
 * real/imag parts of exp(-pi.i.n/adim) n=0,...adim-1.
 */
  theta = isign*pi/adim;
  sininc=sin(theta);
  cosinc=cos(theta);
/*
 * Unpack each row.
 */
  rowptr = image;
  rowinc = 2*(adim+1);
  cnjptr = image + rowinc-2;
  for(irow=0; irow<bdim; irow++, rowptr+=rowinc, cnjptr+=rowinc) {
    rnptr = rowptr;   /* Pointer to real part of x=0,y=irow */
    inptr = rnptr+1;  /* Pointer to imag part of x=0,y=irow */
    rmptr = cnjptr;   /* Pointer to real part of x=adim-0,y=bdim-irow */
    imptr = rmptr+1;  /* Pointer to imag part of x=adim-0,y=bdim-irow */
/*
 * Process column pair 0 and N outside the main loop.
 */
    if(isign == 1) {
      *rmptr =  *rnptr - *inptr;
      *rnptr += *inptr;
      *inptr = *imptr = 0.0f;
    } else {
      *inptr = norm * 0.5f * (*rnptr - *rmptr);
      *rnptr = norm * 0.5f * (*rnptr + *rmptr);
    };
/*
 * Initial value of exp(-pi.i.n/adim) for n=1.
 */
    fwr = (float) (wr = cosinc);   /* Real part */
    fwi = (float) (wi = sininc);   /* Imag part */
/*
 * Process column pairs n,N-n from n=1 to n=N/2-1.
 */
    rnptr+=2;
    inptr+=2;
    rmptr-=2;
    imptr-=2;
    for(icol=1;icol<=adim/2; icol++) {
/*
 * sum_a(n) = 0.5( rowptr[n] + conjugate(rowptr[adim-n]) ).
 */
      rsum_a = 0.5f * (*rnptr + *rmptr);
      isum_a = 0.5f * (*inptr - *imptr);
/*
 * sum_b(n) = 0.5i*exp(pi.i.n/adim)*(rowptr[n]-conjugate(rowptr[adim-n]))
 */
      rsum_b = -scal * (*inptr + *imptr);
      isum_b =  scal * (*rnptr - *rmptr);
/*
 * The new data points are data(n)=sum_a(n)+sum_b(n)*exp(-pi.i.n/adim)
 */
      *rnptr = norm * (rsum_a + fwr * rsum_b - fwi * isum_b);
      *inptr = norm * (isum_a + fwr * isum_b + fwi * rsum_b);
/*
 * SImilarly, data(adim-n)=conjugate(sum_a(n)) - conjugate(sum_b(n)) *
 *                         conjugate(exp(-pi.i.n/adim))
 */
      *rmptr = norm * ( rsum_a - fwr * rsum_b + fwi * isum_b);
      *imptr = norm * (-isum_a + fwr * isum_b + fwi * rsum_b);
/*
 * Compute sin() and cos() of theta plus the current angle.
 */
      wtmp = wr;
      wr = wtmp*cosinc - wi*sininc;
      wi = wtmp*sininc + wi*cosinc;
      fwr = (float) wr;  /* Single precision is sufficient */
      fwi = (float) wi;  /* outside the recurrence */
/*
 * Next even/odd pair of elements inwards from each end of the array.
 */
      rnptr+=2;
      inptr+=2;
      rmptr-=2;
      imptr-=2;
    };
  };
/*
 * Remove the redundant extra column at n_x=N_x/2 before the inverse fft.
 */
  if(isign == -1)
    cls_fft(image, adim, bdim);
  return;
}

/*.......................................................................
 * Internal realfft() function to take an adim*bdim complex transform
 * and insert an extra complex element at the end of each row
 * thus turning the array into an (adim+1)*bdim complex array as required
 * for the conjugate symmetric array. The extra element is effectively
 * n_x = N_x/2 and because of the N/2 periodicity in the underlying
 * transform is copied from n_x=0.
 */
static void opn_fft(float *image, int adim, int bdim)
{
  int oldinc=2*adim;  /* Number of floats per row in original array */
  float *origptr;     /* Pointer to a float of original array to be copied */
  float *destptr;     /* Pointer to destination of *origptr original float */
  float *finptr;      /* Pointer to the final float to be moved */
  float *prevptr;     /* Points to last float of the previous row */
/*
 * In order not to overwrite originals prematurely, start at the end of
 * the two arrays.
 */
  origptr =image+2*adim*bdim-1; 
  destptr = image+2*(adim+1)*bdim-3;
/*
 * The two extra floats at the end of each row should be copied
 * from the first two elements of that row. ie. f(N/2)=f(0).
 * In the loop below this assignment is performed at the end of
 * the loop, for the next row back so do the current row here.
 */
  *(destptr+1) = *(origptr-oldinc+1);
  *(destptr+2) = *(origptr-oldinc+2);
/*
 * The first row is already in place so set the finish pointer to the
 * float following it last float.
 */
  finptr = image+oldinc;
/*
 * Copy a row at a time.
 */
  while(origptr > finptr) {
    prevptr = origptr-oldinc;
    while(origptr != prevptr)
      *(destptr--) = *(origptr--);
/*
 * Insert the two extra floats at the end of the next complex row.
 * These are copies of the first two floats on that row.
 */
    *(destptr--) = *(origptr-oldinc+2);  /* Imag part */
    *(destptr--) = *(origptr-oldinc+1);  /* Real part */
  };
  return;
}

/*.......................................................................
 * Internal realfft() function to take an (adim+1)*bdim complex array
 * and remove the redundant (0,0) complex element at the end of each row
 * thus turning the array into an adim*bdim complex array suitable for
 * fft'ing.
 */
static void cls_fft(float *image, int adim, int bdim)
{
  int oldinc=2*(adim+1);/* Number of floats per row in original array */
  int newinc=2*adim;  /* Number of floats per row in new array */
  float *origptr;     /* Pointer to a float of original array to be copied */
  float *destptr;     /* Pointer to destination of *origptr original float */
  float *finptr;      /* Pointer to the final float to be moved */
  float *nextptr;     /* Points to first float in the next row */
/*
 * Start at the start of the second row - the first row doesn't need
 * moving.
 */
  origptr = image+oldinc;
  destptr = origptr-2;
/*
 * Determine the final pointer.
 */
  finptr = image+bdim*oldinc;
/*
 * Copy a row at a time.
 */
  while(origptr < finptr) {
    nextptr = destptr+newinc;
    while(destptr < nextptr)
      *(destptr++) = *(origptr++);
/*
 * Skip the following redundant complex pair of floats.
 */
    origptr += 2;
  };
/*
 * Zero the unused elements at the end of the array - destptr already
 * points to the first of these and there are bdim*2 of them.
 */
  finptr = destptr + 2*bdim;
  while(destptr < finptr)
    *(destptr++) = 0.0f;
  return;
}

/*.......................................................................
 * This is a modified version of fft_shift() for applying a phase
 * shift to the half conjugate symmetric UV data array taken by
 * realfft() (inverse transform only) in order to make the resulting
 * image have its center at N/2,N/2. For the forward transform use
 * fft_shift(). Note that the dimensions should be the same as
 * required by realfft() ie. adim+2,bdim floats arranged as a 1D
 * array. The reason that you must use this rather than fft_shift for
 * the inverse transform, is simply the extra pair of floats on the
 * first dimension. The phase shift specified by the shift theorem,
 * for an adim/2,bdim/2 coordinate shift is
 * exp(-pi*i*(adim*u+bdim*v)), which is simply equal to -1 when
 * (adim*u+bdim*v) is odd and +1 when it is even.
 *
 * Input/Output:
 *  image  float *   A 1D array with 2(adim/2+1)*bdim elements, where adim+2
 *                   and bdim denote the dimensions of the equivalent
 *                   2D array.  The factor of 2 indicates that each
 *                   element of this array is made up of two
 *                   consecutive floats, the first is the real part
 *                   and the second, the imaginary part.  adim and
 *                   bdim must be integral powers of two
 * Input:
 *  adim   int       There should be adim+2 floats along the first
 *                   dimension (fastest changing index) of the
 *                   2D array equivalent of 'image'.
 *  bdim   int       The second dimension of the 2D array equivalent of
 *                   'image'.
 */
void cnj_shift(float *image, int adim, int bdim)
{
  int ib;            /* Step count along 2nd dimension */
  float *fptr=image; /* Pointer to next element to be scaled by -1 */
/*
 * Multiply every second complex element by -1.
 */
  fptr += 2;
  for(ib=0; ib<bdim*(adim+2)/4; ib++, fptr += 4) {
    *fptr *= -1.0f;     /* Real part */
    *(fptr+1) *= -1.0f; /* Imaginary part */
  };
  return;
}

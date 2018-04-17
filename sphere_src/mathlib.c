/*
  This file contains an assortment of user functions and user
  accessible variables concerned with these functions.
*/

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <math.h>

#include "sphere.h"
#include "utils.h"
#include "helpdir.h"

/*
 * The conversion factor to convert from the FWHM of a Guassian
 * to its standard deviation [ie. 1/sqrt(8*ln(2))]
 */
#define FWHM_TO_STDDEV 0.4246609001

/*
  Declare variables that are to be aliased as user variables below. Only
  float, integer, char and logical variables are supported.
  NB. character strings must NOT be initialised here unless they are marked
  as R_ONLY parameters. This is to allow variable length strings where
  the previous string is often free'd first on the assumption that the
  memory for the string was allocated using malloc(), not by the
  compiler).
*/

float pi = 3.14159265;
static float period=0.0;
static float grad, yint, yinterr,graderr;

static Descriptor filter    = {'f' , '1' ,RWD    ,20,{20,1,1}, NULL};
static Descriptor amplitude = {'f' , '1' ,RWD    ,20,{20,1,1}, NULL};
static Descriptor phase     = {'f' , '1' ,RWD    ,20,{20,1,1}, NULL};

static Descriptor mathv_type[] = {
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &pi},
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &period},
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &grad},
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &yint},
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &yinterr},
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &graderr},
   {'D' , '0' ,NO_DEL ,1, {1,1,1}, &filter},
   {'D' , '0' ,NO_DEL ,1, {1,1,1}, &amplitude},
   {'D' , '0' ,NO_DEL ,1, {1,1,1}, &phase}
};

/*
  In the same order as the above array of types, define the array
  of user names for the arrays.
*/

static char *mathv_name[] = {
   "pi",
   "period",
   "gradient",
   "yintercept",
   "yinterr",
   "graderr",
   "filter",
   "amplitude",
   "phase"
};

/*
  Declare the user functions here.
*/

static Template(sin_fn);
static Template(cos_fn);
static Template(tan_fn);
static Template(asin_fn);
static Template(acos_fn);
static Template(atan_fn);
static Template(sqrt_fn);
static Template(abs_fn);
static Template(ln_fn);
static Template(log_fn);
static Template(int_fn);
static Template(nint_fn);
static Template(real_fn);
static Template(exp_fn);
static Template(atan2_fn);
static Template(mod_fn);
static Template(gran_fn);
static Template(uran_fn);
static Template(min_fn);
static Template(max_fn);
static Template(seed_fn);
static Template(mean_fn);
static Template(sum_fn);
static Template(rms_fn);
static Template(ramp_fn);
static Template(fht_fn);
static Template(smooth_fn);
static Template(minmax_fn);
static Template(four_fn);
static Template(fsval_fn);
static Template(trans_fn);
static Template(sort_fn);
static Template(fold_fn);
static Template(integ_fn);
static Template(median_fn);
static Template(correl_fn);
static Template(fitline_fn);
static Template(grid_fn);
static Template(garray_fn);

/*
  Declare the function types below.
*/

static Functype mathf_type[] = {
   {sin_fn,   NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {cos_fn,   NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {tan_fn,   NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {asin_fn,  NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {acos_fn,  NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {atan_fn,  NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {sqrt_fn,  NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {abs_fn,   NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {ln_fn,    NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {log_fn,   NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {int_fn,   NORM    , 1,1,    "if",   "00",     "vv",   0 },
   {nint_fn,  NORM    , 1,1,    "if",   "00",     "vv",   0 },
   {real_fn,  NORM    , 1,1,    "fn",   "00",     "vv",   0 },
   {exp_fn,   NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {atan2_fn, NORM    , 2,2,    "fff",  "000",    "vvv",  0 },
   {mod_fn,   NORM    , 2,2,    "fff",  "000",    "vvv",  0 },
   {gran_fn,  NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {uran_fn,  NORM    , 1,1,    "ff",   "00",     "vv",   0 },
   {min_fn,   NORM    , 1,MAXARG,"ff",  "00",     "vv",   0 },
   {max_fn,   NORM    , 1,MAXARG,"ff",  "00",     "vv",   0 },
   {seed_fn,  NORM    , 1,1,    " i",   " 0",     " v",   1 },
   {mean_fn,  NORM    , 1,1,    "ff",   "0*",     "vv",   1 },
   {sum_fn,   NORM    , 1,1,    "ff",   "0*",     "vv",   1 },
   {rms_fn,   NORM    , 1,1,    "ff",   "0*",     "vv",   1 },
   {ramp_fn,   NORM   , 2,3,    "fff",  "100",    "vvv",  1 },
   {fht_fn,    NORM   , 1,1,    "ff",   "*2",     "vv",   1 },
   {smooth_fn, NORM   , 2,3,    "ffff", "*200",   "vvvv", 1 },
   {minmax_fn, NORM   , 1,1,    "ff",   "23",     "vv",   1 },
   {four_fn,   NORM   , 2,2,     " ff", " 11",    " vv",  1 },
   {fsval_fn,  NORM   , 1,2,     "ffi", "000",    "vvv",  0 },
   {trans_fn,  NORM   , 2,2,     "ffC", "*30",    "vvv",  1 },
   {sort_fn,   NORM   , 2,MAXARG," fff"," 01*",   " vrr", 1 },
   {fold_fn,   NORM   , 3,MAXARG," ffff"," 001*", " vvrr",1 },
   {integ_fn,  NORM   , 3,4,     "fffff","301*0", "vvvvv",1 },
   {median_fn, NORM   , 1,1,     "ff",  "0*",     "vv",   1 },
   {correl_fn, NORM   , 2,2,     "fff", "0**",    "vvv",  1 },
   {fitline_fn,NORM   , 2,3,     " fff"," 111",   " vvv", 1 },
   {grid_fn,   NORM   , 5,5,     " iffff", " 01*1*",  " vvvrr",1 },
   {garray_fn, NORM   , 3,3,      "ffff", "1111", "vvvv", 1 },
};

/*
  In the same order as the above array of types, define the array
  of user names for the functions.
*/

static char *mathf_name[] = {
   "sin",
   "cos",
   "tan",
   "asin",
   "acos",
   "atan",
   "sqrt",
   "abs",
   "ln",
   "log",
   "int",
   "nint",
   "real",
   "exp",
   "atan2",
   "mod",
   "gauss_rand",
   "uniform_rand",
   "min",
   "max",
   "seed_rand",
   "mean",
   "sum",
   "rms",
   "ramp",
   "hartley",
   "smooth",
   "minmax",
   "fourier",
   "fsval",
   "transpose",
   "sort",
   "fold",
   "trap_int",
   "median",
   "correl",
   "fit_line",
   "grid",
   "gauss_array",
};

/*
  Record the above declarations etc for this module in a global
  structure for use when building the main symbol table.
*/

Module m_maths = {
  "maths",
  HELP_DIR,
  NULL, 0,
  mathv_type, mathv_name, COUNT(mathv_name),
  mathf_type, mathf_name, COUNT(mathf_name)
};

/*
  Local variables and functions.
*/

static float fnum_a, fnum_b, fnum;

/*.......................................................................
  Take the trigonomeric sin() of a single number (in radians).
*/
static Template(sin_fn)
{
        *FLTPTR(outvals) = (float) sin((double) *FLTPTR(invals[0]));
        return no_error;
}

/*.......................................................................
  Take the trigonomeric cos() of a single number (in radians).
*/
static Template(cos_fn)
{
        *FLTPTR(outvals) = (float) cos((double) *FLTPTR(invals[0]));
        return no_error;
}

/*.......................................................................
  Take the trigonomeric tan() of a single number (in radians).
*/
static Template(tan_fn)
{
        *FLTPTR(outvals) = (float) tan((double) *FLTPTR(invals[0]));
        return no_error;
}

/*.......................................................................
  Take the trigonomeric asin() of a single number - returns radians.
*/
static Template(asin_fn)
{
        fnum = *FLTPTR(invals[0]);
        if( fnum > 1.0 || fnum < -1.0) {
          lprintf(stderr,"Illegal operand value: asin(%f)\n",fnum);
          return -1;
        };
        *FLTPTR(outvals) = (float) asin((double) fnum);
        return no_error;
}

/*.......................................................................
  Take the trigonomeric acos() of a single number - returns radians.
*/
static Template(acos_fn)
{
        fnum = *FLTPTR(invals[0]);
        if( fnum > 1.0 || fnum < -1.0) {
          lprintf(stderr,"Illegal operand value: acos(%f)\n",fnum);
          return -1;
        };
        *FLTPTR(outvals) = (float) acos((double) fnum);
        return no_error;
}

/*.......................................................................
  Take the trigonomeric atan() of a single number - returns radians.
*/
static Template(atan_fn)
{
        fnum = *FLTPTR(invals[0]);
        *FLTPTR(outvals) = (float) atan((double) fnum);
        return no_error;
}

/*.......................................................................
  Take the trigonomeric atan2(x,y) two operands - returns radians.
*/
static Template(atan2_fn)
{
        *FLTPTR(outvals) = (float)
         atan2((double) *FLTPTR(invals[0]), (double) *FLTPTR(invals[1]));
        return no_error;
}

/*.......................................................................
  Take the square-root of a single number.
*/
static Template(sqrt_fn)
{
        fnum = *FLTPTR(invals[0]);
        if( fnum < 0.0) {
          lprintf(stderr,"Illegal operand value: sqrt(%f)\n",fnum);
          return -1;
        };
        *FLTPTR(outvals) = (float) sqrt((double) fnum);
        return no_error;
}

/*.......................................................................
  Return the absolute value of a single number.
*/
static Template(abs_fn)
{
        fnum = *FLTPTR(invals[0]);
        *FLTPTR(outvals) = (float) fabs((double) fnum);
        return no_error;
}

/*.......................................................................
  Return the fractional remainder after a division.
*/
static Template(mod_fn)
{
        fnum_a = *FLTPTR(invals[0]);
        fnum_b = *FLTPTR(invals[1]);
        if(fnum_b == 0.0) {
          lprintf(stderr,"Divide by zero error: mod(%f,%f)\n",fnum_a,fnum_b);
          return -1;
        };
        *FLTPTR(outvals) = (float) fmod((double) fnum_a, (double) fnum_b);
        return no_error;
}


/*.......................................................................
  Return a random number from a gaussian distrubution with standard
  deviation equal to the single argument.
*/
static Template(gran_fn)
{
        fnum = *FLTPTR(invals[0]);
        *FLTPTR(outvals) = gauss_rand(fnum);
        return no_error;
}


/*.......................................................................
  Return a random number from a uniform probability distribution, between
  - and + the number argument.
*/
static Template(uran_fn)
{
        fnum = *FLTPTR(invals[0]);
        *FLTPTR(outvals) = uniform_rand(fnum);
        return no_error;
}

/*.......................................................................
  Find and return the minimum value of the user's scalar arguments.
*/
static Template(min_fn)
{
        static int i;
	static float minval;
        minval = *FLTPTR(invals[0]);
	for(i=1;i<npar; i++) {
	  if( (fnum = *FLTPTR(invals[i])) < minval) minval=fnum;
	};
	*FLTPTR(outvals) = minval;
	return no_error;
}

/*.......................................................................
  Find and return the maximum value of the user's scalar arguments.
*/
static Template(max_fn)
{
        static int i;
	static float maxval;
        maxval = *FLTPTR(invals[0]);
	for(i=1;i<npar; i++)
	  if( (fnum= *FLTPTR(invals[i])) > maxval) maxval=fnum;
	*FLTPTR(outvals) = maxval;
	return no_error;
}

/*.......................................................................
  Re-seed the random number generator using the single argument as the
  seed.
*/
static Template(seed_fn)
{
        float frand(unsigned int iseed);
	frand((unsigned int) *INTPTR(invals[0]));
	return no_error;
}

/*.......................................................................
  Take the natural log of a single number.
*/
static Template(ln_fn)
{
        fnum = *FLTPTR(invals[0]);
        if( fnum <= 0.0) {
          lprintf(stderr,"Illegal operand value: ln(%f)\n",fnum);
          return -1;
        };
        *FLTPTR(outvals) = (float) log((double) fnum);
        return no_error;
}

/*.......................................................................
  Take the natural log to the base 10, of a single number.
*/
static Template(log_fn)
{
        fnum = *FLTPTR(invals[0]);
        if( fnum <= 0.0) {
          lprintf(stderr,"Illegal operand value: log(%f)\n",fnum);
          return -1;
        };
        *FLTPTR(outvals) = (float) log10((double) fnum);
        return no_error;
}


/*.......................................................................
  Return the integer value of a single number.
*/
static Template(int_fn)
{
        *INTPTR(outvals) = *FLTPTR(invals[0]);
        return no_error;
}

/*.......................................................................
 * Return the nearest integer to a given floating point number.
 */
static Template(nint_fn)
{
  *INTPTR(outvals) = floor(*FLTPTR(invals[0]) + 0.5f);
  return no_error;
}

/*.......................................................................
  Return the float representation of an integer argument.
*/
static Template(real_fn)
{
  switch (invals[0]->atyp) {
  case 'i':
    *FLTPTR(outvals) = *INTPTR(invals[0]);
    break;
  case 'f':
    *FLTPTR(outvals) = *FLTPTR(invals[0]);
    break;
  default:
    *FLTPTR(outvals) = 0.0f;
    lprintf(stderr, "real(): Unrecognised type\n");
    break;
  };
  return no_error;
}

/*.......................................................................
  Return the value of e raised to the single operand. ie e^x.
*/
static Template(exp_fn)
{
        *FLTPTR(outvals) = (float) exp((double) *FLTPTR(invals[0]));
        return no_error;
}

/*.......................................................................
  Return the mean of a multi-dimensional array.
*/
static Template(mean_fn)
{
        int num_el,i;
        float *fptr;
	float runmean,nmean;
/*
  Determine the number of elements in the array.
*/
        num_el = 1;
        for(i=0;i<3;i++)
          num_el *= invals[0]->adim[i];
/*
  Get the pointer to the first element.
*/
        fptr = FLTPTR(invals[0]);
/*
  Sum the array. To avoid precision problems due to adding small
  numbers to an ever increasing sum, form the running mean.
  To get the sum this is then multiplied by NL.
*/
        nmean=runmean=0.0f;
        for(i=0; i<num_el; i++)
          runmean += (*(fptr++)-runmean) / ++nmean;
/*
  Return the result.
*/
        *FLTPTR(outvals) = runmean;
        return no_error;
}

/*.......................................................................
  Return the mean of a multi-dimensional array.
*/
static Template(sum_fn)
{
        int num_el,i;
/*
  Determine the number of elements in the array.
*/
        num_el = 1;
        for(i=0;i<3;i++)
          num_el *= invals[0]->adim[i];
/*
  Find the mean of the array.
*/
        mean_fn(invals, 1, outvals);
/*
  Turn the running mean into the required sum.
*/
        *FLTPTR(outvals) =  *FLTPTR(outvals) * num_el;
        return no_error;
}

/*.......................................................................
  Return the standard deviation of a multi-dimensional array.
*/
static Template(rms_fn)
{
        int num_el,i;
        float *fptr,mean_val,temp;
/*
  Determine the number of elements in the array.
*/
        num_el = 1;
        for(i=0;i<3;i++)
          num_el *= invals[0]->adim[i];
/*
  Find the mean of the array.
*/
        mean_fn(invals, 1, outvals);
	mean_val = *FLTPTR(outvals);
/*
  Get the pointer to the first element.
*/
        fptr = FLTPTR(invals[0]);
/*
  Find the mean square difference of the data array values and the mean.
  To avoid precision problems due to adding small
  numbers to an ever increasing sum, form the running mean.
*/
        fnum = 0;
        for(i=1; i<=num_el; i++) {
	  temp = *(fptr++)-mean_val;
          fnum += (temp*temp-fnum)/((float) i);
        };
/*
  Square root the mean square difference to get the standard deviation.
  Also correct for the fact that one degree of freedom was used to form
  mean. This bias is removed by scaling by n/(n-1).
*/
        *FLTPTR(outvals) = sqrt(fnum*num_el/((num_el > 1) ? num_el-1:1));
        return no_error;
}


/*.......................................................................
  A test array-return function. It takes three arguments and returns
  an array of numbers that increments over the range specified in the
  first two arguments, with the increment (from element to element)
  given by the third argument.
*/
static Template(ramp_fn)
{
        float start_val,end_val,inc_val;
        int nvals,i;
/*
  Copy the user specifications for the ramp into local variables.
*/
        start_val=*FLTPTR(invals[0]);
        end_val=*FLTPTR(invals[1]);
/*
  Specification of the increment is optional.
*/
        if(npar == 3)
          inc_val=*FLTPTR(invals[2]);
        else
          inc_val= (end_val > start_val) ? 1.0 : -1.0;
/*
  Test the sign of the increment against that of the range.
*/
        if(inc_val == 0 || (end_val > start_val && inc_val < 0) ||
        (end_val < start_val && inc_val > 0)) {
          lprintf(stderr,"Illegal increment value in ramp(%f,%f,%f)",
           start_val,end_val,inc_val);
          return -1;
        };
/*
  Determine the number of elements that the return array will require.
*/
        nvals = 1+(end_val-start_val)/inc_val;
/*
  Allocate memory for the return array.
*/
        if((VOIDPTR(outvals)=valof_alloc(nvals,'f')) == NULL)
          return -1;
/*
  Build the ramp.
*/
        for(i=0;i<nvals;i++)
          *(FLTPTR(outvals)+i)=start_val+i*inc_val;
/*
  Fill in the dimensional aspects of the return descriptor.
*/
        outvals->num_el=nvals;
        outvals->adim[0] = nvals;
        return no_error;
}

/*.......................................................................
  Take the Fast Hartley Transform (FHT) of an input 2-D array, returning
  the result as the return value of the function.
*/
static Template(fht_fn)
{
        int xnum,ynum;
/*
  Determine the number of elements on each dimension of the input
  array.
*/
	xnum = invals[0]->adim[0];
	ynum = invals[0]->adim[1];
/*
  Check that both dimensions have a number of elements which is
  a power of 2. This check is also made in the function that performs
  the transform, so it is slightly superfluous. However, by checking
  here, we avoid an inefficient allocation and de-allocation of the
  large array that will contain the return value.
*/
	if(!is_pow_of_two(xnum) || !is_pow_of_two(ynum)) {
	  lprintf(stderr, "Illegal array size (%d,%d) - not a power of two - sent to\nthe Fast Hartley Transform function.\n", xnum, ynum);
	  return -1;
	};
/*
  The return descriptor will be identical to the input descriptor
  except that it will be given a different valof pointer.
*/
        *outvals = *invals[0];
/*
  Allocate memory for the return array.
*/
        if((VOIDPTR(outvals)=valof_alloc(invals[0]->num_el,'f')) == NULL)
          return -1;
/*
  Copy the input array into the return array.
*/
	memcpy(FLTPTR(outvals), FLTPTR(invals[0]), xnum*ynum*sizeof(float));
/*
  Perform the Hartley transform.
*/
	if(two_dim_FHT(FLTPTR(outvals), xnum, ynum, 1) == -1)
	  return -1;
        return no_error;
}

/*.......................................................................
  Smooth a 2-D array via the Hartley plane. Convolution by a gaussian
  is performed by multiplication in the Hartley plane. The Gaussian
  FWHM widths along each dimension are given by the user, in units
  of 1 channel.
*/
static Template(smooth_fn)
{
        int i,j,xnum,ynum;
	double xwid,ywid,mul_fac;
	const double two_pi_sq = 19.739209;
	float *pos_ptr, *neg_ptr, *fptr;
/*
  Take the Hartley transform of the input data array. The hartley
  transform is returned in outvals.
*/
	if(fht_fn(invals,1,outvals) == -1)
	  return -1;
/*
  Determine the number of elements on each dimension of the output
  array.
*/
	xnum = outvals->adim[0];
	ynum = outvals->adim[1];
/*
  Ascertain the Gaussian smoothing widths.
*/
	xwid = (double) *FLTPTR(invals[1]);
	ywid = (npar == 3) ? (double) *FLTPTR(invals[2]) : 0.0;
/*
  Convert these FWHM Gaussian widths to the equivalent standard
  deviation. sigma = FWHM/2.sqrt(2ln2).
*/
	xwid = xwid*0.424661;
	ywid = ywid*0.424661;
/*
  The Fourier transform of the probability gaussian:
   1/(sqrt(2.pi).sigma) exp(-x^2/(2.sigma^2)) <==> e^(2.(pi.sigma.f)^2
  Turn xwid and ywid into their fourier plane (scaling) equivalents.
  NB. Take into acount the conversion factor between elements
  and frequency in the Fourier plane. ie. df=di/N where N is
  the total number of elements in the transform.
*/
	xwid = -two_pi_sq * xwid * xwid/(xnum * xnum);
	ywid = -two_pi_sq * ywid * ywid/(ynum * ynum);
/*
  Multiply the first dimension by its smoothing Gaussian.
  The FHT of a gaussian is wholly real so we need to multiply
  both the negative and positive frequencies by the same
  factors. The pointers, pos_ptr and neg_ptr will
  be used to step along the positive and negative parts of the
  FHT. Do the first dimension first. In order bot to have to
  evaluate the same exponential separately for every element
  along the second dimension it is advisable to treat the two
  dimensions separately.
*/
	fptr = FLTPTR(outvals);
	for(i=1; i<=xnum/2; i++) {
/*
  Get pointers to the positive and negative frequencies, i and -i.
*/
	  pos_ptr=fptr+i;
	  neg_ptr=fptr+xnum-i;
/*
  Evaluate the Gaussian FFT value for this frequency.
*/
	  mul_fac = exp(xwid*i*i);
/*
  Multiply through the +ve and -ve frequency stripes (parallel
  to the y-axis).
*/
	  for(j=0; j<ynum; j++) {
	    *pos_ptr = *pos_ptr * mul_fac;
/*
  The last frequency is both -ve and +ve! - don't scale it a second time.
*/
	    if(neg_ptr != pos_ptr) *neg_ptr = *neg_ptr * mul_fac;
/*
  Step the pointers on to the next y-axis positions.
*/
	    pos_ptr += xnum;
	    neg_ptr += xnum;
	  };
	};
/*
  Now the second dimension.
*/
	for(i=1; i<=ynum/2; i++) {
/*
  Get pointers to the positive and negative frequencies, i and -i.
*/
	  pos_ptr=fptr+i*xnum;
	  neg_ptr=fptr+(ynum-i)*xnum;
/*
  Evaluate the Gaussian FFT value for this frequency.
*/
	  mul_fac = exp(ywid*i*i);
/*
  Multiply through the +ve and -ve frequency stripes (parallel
  to the y-axis).
*/
	  for(j=0; j<xnum; j++) {
	    *pos_ptr = *pos_ptr * mul_fac;
/*
  The last frequency is both -ve and +ve! - don't scale it a second time.
*/
	    if(neg_ptr != pos_ptr) *neg_ptr = *neg_ptr * mul_fac;
/*
  Step the pointers on to the next y-axis positions.
*/
	    pos_ptr++;
	    neg_ptr++;
	  };
	};
/*
  Transform the multiplied array back to the image plane.
*/
	if(two_dim_FHT(FLTPTR(outvals), xnum, ynum, 0) == -1)
	  return -1;
        return no_error;
}

/*.......................................................................
  Take an array of up to three dimensions and return a 4x2 element array
  containing the minimum value of the array and its x,y,z position
  in the array and the same for the maximum.
*/
static Template(minmax_fn)
{
        int x,y,z, xmax,ymax,zmax, xmin,ymin,zmin;
	float maxval, minval,*fptr,temp;
/*
  Allocate memory for the return array.
*/
        if((VOIDPTR(outvals)=valof_alloc(8,'f')) == NULL)
          return -1;
/*
  Fill in the dimensional aspects of the return descriptor.
*/
        outvals->adim[0] = 4;
	outvals->adim[1] = 2;
	outvals->adim[2] = 1;
	outvals->num_el = 8;
/*
  Search the array for its maximum.
*/
	fptr = FLTPTR(invals[0]);
	xmax = ymax = zmax = xmin = ymin = zmin = 0;
	maxval = minval = *fptr;
	for(z=0; z<invals[0]->adim[2]; z++) {
	  for(y=0; y<invals[0]->adim[1]; y++) {
	    for(x=0; x<invals[0]->adim[0]; x++) {
	      temp = *fptr;
	      if(temp > maxval) {
		maxval = temp;
		xmax=x; ymax=y; zmax=z;
	      }
	      else if(temp < minval) {
		minval= temp;
		xmin=x; ymin=y; zmin=z;
	      };
	      fptr++;
	    };
	  };
	};
/*
  Record the results in the return array.
*/
	FLTPTR(outvals)[0] = minval;
	FLTPTR(outvals)[1] = xmin+1;
	FLTPTR(outvals)[2] = ymin+1;
	FLTPTR(outvals)[3] = zmin+1;
	FLTPTR(outvals)[4] = maxval;
	FLTPTR(outvals)[5] = xmax+1;
	FLTPTR(outvals)[6] = ymax+1;
	FLTPTR(outvals)[7] = zmax+1;
        return no_error;
}

/*.......................................................................
  Decompose the time and data arrays into fourier series components.
  The function deposits the resulting amplitudes and phases in the global
  user variables of the same name. It also uses the user variable, period,
  to set the fundamental period to be used.
*/
static Template(four_fn)
{
int fourier_series(float x_data[], float y_data[], int npts, float period, float asum[], float bsum[], int max_order);
        int npts, max_order;
	float *amp, *ph;
/*
  Get pointers to the amplitude and phase arrays.
*/
	amp = FLTPTR(&amplitude);
	ph  = FLTPTR(&phase);
/*
  Set max_order equal to the minimum number of elements in the two arrays.
*/
	if(amplitude.adim[0] < phase.adim[0])
	  max_order = amplitude.adim[0];
	else
	  max_order = phase.adim[0];
/*
  Similarly for the number of data points in the time and data arrays.
*/
	if(invals[0]->adim[0] < invals[1]->adim[0])
	  npts = invals[0]->adim[0];
	else
	  npts = invals[1]->adim[0];
/*
  Evaluate the fourier series.
*/
	return fourier_series(FLTPTR(invals[0]), FLTPTR(invals[1]),
			      npts, period, amp, ph, max_order);
}

/*.......................................................................
  Return the fourier series value at the time given as argument.
  The function uses the amplitudes and phases and the filter values
  in the global user variables of the same name. It also uses the user
  variable, period, to set the fundamental period to be used. If a second
  argument was given, it will be interpretted as the differential order.
*/
static Template(fsval_fn)
{
        int fourier_series_value(float xval, float *yval, int differential_order, float period, float amp[], float phase[], float filter[], int max_order);
        int max_order, diff_order;
        float *amp, *ph, *filt;
/*
  Get pointers to the amplitude and phase arrays.
*/
	amp = FLTPTR(&amplitude);
	ph  = FLTPTR(&phase);
	filt= FLTPTR(&filter);
/*
  Set max_order equal to the minimum number of elements in the three arrays.
*/
	if(amplitude.adim[0] < phase.adim[0]) {
	  if(amplitude.adim[0] < filter.adim[0])
	    max_order = amplitude.adim[0];
	  else
	    max_order = filter.adim[0];
	}
	else {
	  if(phase.adim[0] < filter.adim[0])
	    max_order = phase.adim[0];
	  else
	    max_order = filter.adim[0];
	};
/*
  Get the diffential order - default to zero order.
*/
	if(npar > 1)
	  diff_order = *INTPTR(invals[1]);
	else
	  diff_order = 0;
/*
  Evaluate the fourier series.
*/
	return fourier_series_value(*FLTPTR(invals[0]),
				    FLTPTR(outvals), diff_order,
				    period, amp, ph, filt, max_order);
}

/*.......................................................................
  Transpose a user n-D array using the specification code in the string
  first argument. The user function returns the transposed array.
*/
static Template(trans_fn)
{
        char *spec;
	int olddim[3], newdim[3], axis[3], newaxis[3], add[3], dim, pos, i,j,k;
	float *inptr, *outptr;
/*
  Get the user specification string.
*/
	spec = *STRPTR(invals[1]);
/*
  Get the dimensions of the original array.
*/
	for(i=0; i<3; i++)
	  olddim[i] = invals[0]->adim[i];
	dim = invals[0]->dim - '0';
/*
  Compare to the length of the specification string.
*/
	if(strlen(spec) < dim) {
	  lprintf(stderr, "trans(): Specification string has too few items.\n");
	  return -1;
	};
/*
  Parse the user specification string. The array axis[i] will contain the
  number of the old axis to be placed as axis number, i, of the new array.
*/
	for(i=0; i<3; i++) {
	  if(i<dim) {
	    axis[i] = pos = spec[i] - '0';
	    if(pos < 0 || pos >= dim) {
	      lprintf(stderr, "trans: Specifier-string item out of range.\n");
	      return -1;
	    };
	  }
	  else
	    axis[i] = pos = i;
/*
  Check that the specified axis was not previously specified as a previous
  item.
*/
	  for(j=0; j<i; j++) {
	    if(pos == axis[j]) {
	      lprintf(stderr, "trans: Duplicate specifier-string item.\n");
	      return -1;
	    };
	  };
	};
/*
  Determine the total number of channels in the input array.
*/
	for(j=1,i=0; i<3; i++)
	  j *= invals[0]->adim[i];
/*
  Allocate memory for the return array.
*/
	if( (VOIDPTR(outvals) = valof_alloc(j,invals[0]->atyp)) == NULL)
	  return -1;
/*
  Fill in the return descriptor items.
*/
	outvals->num_el = j;
	for(i=0; i<3; i++)
	  outvals->adim[i] = newdim[i] = olddim[axis[i]];
/*
  Transpose the arrays.
*/
	inptr = FLTPTR(invals[0]);
	outptr = FLTPTR(outvals);
/*
  Determine the order in which the axes of the transposed array
  should be stepped through.
*/
	for(i=0; i<3; i++)
	  newaxis[axis[i]] = i;
/*
  Calculate the increment in elements for each axis of the new array
  corresponding to an increment of 1 position along the original array
  axis.
*/
	get_increments(newaxis,newdim,add);
/*
  Transpose.
*/
	for(i=0; i<olddim[2]; i++) {
	  for(j=0; j<olddim[1]; j++) {
	    for(k=0; k<olddim[0]; k++) {
	      *outptr = *inptr;
	      inptr++;
	      outptr += add[0];
	    };
	    outptr += add[1];
	  };
	  outptr += add[2];
	};
	return no_error;
}

/*.......................................................................
  Given an array to be used to produce an index array, and at least one
  data array, sort the data arrays. All the arrays must have the same
  number of elements.
*/
static Template(sort_fn)
{
        int nskip, npts,i,j,k, arg, axis[3], ndim[3], add[3], sort_axis;
	float *work, *work_ptr, *data_ptr, *outptr;
	int *index, *index_ptr;
/*
  Find out the number of the axis to be sorted.
*/
	sort_axis = (int) *FLTPTR(invals[0]);
/*
  Check its legality.
*/
	if(sort_axis < 0 || sort_axis > 2) {
	  lprintf(stderr, "sort(): Axis specification (%d) out of bounds.\n", sort_axis);
	  return -1;
	};
/*
  Get the number of points in the index source array.
*/
	npts = invals[1]->adim[0];
/*
  Check that all the input arrays have the same number of elements
  along their y-axes.
*/
	for(i=2; i<npar; i++) {
	  if(invals[i]->adim[sort_axis] != npts) {
	    lprintf(stderr, "sort(): The argument arrays have differing numbers of elements along the requested axis.\n");
	    return -1;
	  };
	};
/*
  One point? don't bother sorting!
*/
	if(npts==1)
	  return no_error;
/*
  Get the index array with which to sort the arrays.
*/
	if(indexx(npts, FLTPTR(invals[1]), &index) == -1)
	  return -1;
/*
  Allocate memory for a work array into which each array will be
  sorted before being copied back to the original.
*/
	if( (work = (float *) calloc(npts, sizeof(float))) == NULL) {
	  lprintf(stderr, "sort: Memory allocation of work array failed.\n");
	  return -1;
	};
/*
  Rearrange the indexing array separately since it is to be re-orderred
  along the x-axis, whereas the others are to be orderred along an
  arbitrary.
*/
	data_ptr=FLTPTR(invals[1]);
	work_ptr=work;
	index_ptr=index;
	for(j=0; j<npts; j++, index_ptr++, work_ptr++)
	  *work_ptr = data_ptr[*index_ptr];
/*
  Copy the sorted array into the original array.
*/
	work_ptr=work;
	for(j=0; j<npts; j++, data_ptr++, work_ptr++)
	  *data_ptr = *work_ptr;
/*
  Assemble an array specifying the order in which to step through the
  array. ie. the inner loop to step through the array should step
  along the sort axis.
*/
	axis[0] = sort_axis;
	for(i=1,j=0; i<3; i++, j++)
	  axis[i] = (j==sort_axis) ? ++j : j;
/*
  Use the index array to sort each of the data arrays in turn
  into the work array and then copy back to overwrite the original.
*/
	for(arg=2; arg<npar; arg++) {
/*
  Determine the dimensions of the current array.
*/
	  for(j=0; j<3; j++)
	    ndim[j] = invals[arg]->adim[j];
/*
  Determine the element increments required to step along the
  required axis and to the start of the next occurence of that
  axis.
*/
	  get_increments(axis,ndim,add);
/*
  Determine the number of elements to be skipped in order to step along
  the axis being sorted.
*/
	  nskip = 1;
	  for(i=0; i<axis[0]; i++)
	    nskip *= ndim[i]; 
/*
  Perform the rearrangements of each array along the required axis.
*/
	  outptr = FLTPTR(invals[arg]);
	  for(i=0; i<ndim[axis[2]]; i++) {
	    for(j=0; j<ndim[axis[1]]; j++) {
/*
  Sort the current sub-array into the work array.
*/
	      data_ptr = outptr;
	      index_ptr = index;
	      work_ptr = work;
	      for(k=0; k<ndim[axis[0]]; k++) {
		*work_ptr = data_ptr[*index_ptr * nskip];
		work_ptr++;
		index_ptr++;
	      };
/*
  Copy the sorted array back into the original array.
*/
	      work_ptr = work;
	      for(k=0; k<ndim[axis[0]]; k++) {
		*outptr = *work_ptr;
		work_ptr++;
		outptr += add[0];
	      };
	      outptr += add[1];
	    };
	    outptr += add[2];
	  };
	};
/*
  Return the memory allocated to work arrays.
*/
	free(index);
	free(work);
	return no_error;
}

/*.......................................................................
  Given a folding period, a 1D array of times, and associated data arrays,
  fold the time array into one period, re-ordering the data arrays to
  keep elemental coincidence.
*/
static Template(fold_fn)
{
	float *time_ptr, start;
	int i,fold_axis, npts;
/*
  Get and check the folding period.
*/
	period = *FLTPTR(invals[0]);
	if(period <= 0.0) {
	  lprintf(stderr, "fold(): Unphysical period: %f\n",period);
	  return -1;
	};
/*
  Find out the number of the axis to be sorted.
*/
	fold_axis = (int) *FLTPTR(invals[1]);
/*
  Check its legality.
*/
	if(fold_axis < 0 || fold_axis > 2) {
	  lprintf(stderr, "fold(): Axis specification (%d) out of bounds.\n", fold_axis);
	  return -1;
	};
/*
  Get the number of points in the time array.
*/
	npts = invals[2]->adim[0];
/*
  Check that all the input arrays have the same number of elements along
  the specified axis.
*/
	for(i=3; i<npar; i++) {
	  if(invals[i]->adim[fold_axis] != npts) {
	    lprintf(stderr, "fold(): The argument arrays have differing numbers of elements along the specified axis.\n");
	    return -1;
	  };
	};
/*
  Find the remainder of dividing by the folding period for each element
  of the time array.
*/
	time_ptr = FLTPTR(invals[2]);
	start = *time_ptr;
	for(i=0; i<npts; i++,time_ptr++)
	  *time_ptr = (float) fmod((double) (*time_ptr-start), period);
/*
  Sort the time array and re-order the rest of the data arrays to
  maintain elemental coincidence.
*/
	return sort_fn(&invals[1],npar-1,outvals);
}

/*.......................................................................
  Integrate an n-D array along a specified axis. Return the array result.
  If an optional 4th (float) argument is specified, then the integral
  is of a periodic function with a period given by this argument. In this
  case the data must stretch less than a period and the integral will be
  over exactly one period. It is unlikely that it will cover exactly one
  period, so the final point in the integral, will be the first point
  placed at 1 period from the start. The result will be divided by 1
  period to give the zero-offset of the periodic function.
*/
static Template(integ_fn)
{
	int indim[3], outdim[3], inaxis[3], inadd[3];
	int integ_axis, i,j,k, npts;
	float *inptr, *outptr, *xptr;
	float wrap_per=0.0f; /* Wrap-around period for periodic functions */
	float wrap_len=0.0f; /* Distance from final point to wrap_per */
	float first_y;	     /* First data value in each integral */
	int do_wrap=0;  /* If a wrap period was given do a periodic integral */
/*
  Find out the number of the axis to be sorted.
*/
	integ_axis = (int) *FLTPTR(invals[0]);
/*
  Check its legality.
*/
	if(integ_axis < 0 || integ_axis > 2) {
	  lprintf(stderr, "integ(): Axis specification (%d) out of bounds.\n", integ_axis);
	  return -1;
	};
/*
  Get the number of points in the 1D array that is to be integrated
  against.
*/
	npts = invals[1]->adim[0];
/*
  One point? Illegal.
*/
	if(npts == 1) {
	  lprintf(stderr, "integ(): Illegal request for integral of 1 point.\n");
	  return -1;
	};
/*
  Check that the data array has the same number of elements
  along the axis to be integrated.
*/
	if(invals[2]->adim[integ_axis] != npts) {
	  lprintf(stderr, "integ(): The number of elements in the x and y arrays differ.\n");
	  return -1;
	};
/*
  Make sure that the x-axis array is in ascending order.
*/
	inptr = FLTPTR(invals[1]);
	for(i=1; i<npts; i++) {
	  if(inptr[i] - inptr[i-1] < 0) {
	    lprintf(stderr, "integ(): x-array not in ascending order.\n");
	    return -1;
	  };
	};
/*
  Also check to see if its first and last values are significantly
  different.
*/
	if(inptr[npts-1] - inptr[0] < 1e-20) {
	  lprintf(stderr,"integ(): x-array not in ascending order.\n");
	  return -1;
	};
/*
  See if the user wants a wrapped periodic integral.
*/
	if(npar > 3) {
	  do_wrap=1;
	  wrap_per = *FLTPTR(invals[3]);
	  if(wrap_per <= 0.0) {
	    lprintf(stderr, "Invalid period: %f\n", wrap_per);
	    return -1;
	  };
/*
  Determine the x-value of the end of the final period following
  the last point, and the difference between this and the last point.
*/
	  wrap_len = wrap_per - inptr[npts-1];
	  if(wrap_len < 0) {
	    lprintf(stderr,"Data covers %f more than one period\n", -wrap_len);
	    return -1;
	  };
	};
/*
  Determine the number of elements per axis of the return array.
  NB. The return array lacks the integrated axis.
  Also get the dimensions of the input data array.
*/
	for(i=0,j=0; i<3; i++) {
	  indim[i] = invals[2]->adim[i];
	  if(i != integ_axis)
	    outdim[j++] = indim[i];
	};
	outdim[j]=1;
/*
  Determine the total number of channels in the return array.
*/
	for(j=1,i=0; i<3; i++)
	    j *= outdim[i];
/*
  Allocate memory for the return array.
*/
	if( (VOIDPTR(outvals) = valof_alloc(j,invals[2]->atyp)) == NULL)
	  return -1;
/*
  Fill in the return-descriptor items.
*/
	outvals->num_el = j;
	for(i=0; i<3; i++)
	  outvals->adim[i] = outdim[i];
/*
  Work out the order in which the axes of the input and output
  arrays should be stepped through. The axis of the input array
  should be stepped through in the innermost loop.
*/
	inaxis[0] = integ_axis;
	for(i=0,j=1;i<3;i++) {
	  if(i != integ_axis)
	    inaxis[j++] = i;
	};
/*
  With these axis orders, determine the element increments required
  to step through the axes in the specified order.
*/
	get_increments(inaxis,indim,inadd);
/*
  Get pointers to the x-axis array, the input data array and the output
  data array.
*/
	xptr = FLTPTR(invals[1]);
	inptr  = FLTPTR(invals[2]);
	outptr = FLTPTR(outvals);
/*
  Perform the integrations.
*/
	for(i=0; i<indim[inaxis[2]]; i++) {
	  for(j=0; j<indim[inaxis[1]]; j++) {
/*
  Keep a record of the 1st data point in case this is to be a wrapped
  integral and initialize the integral output value to zero.
*/
	    first_y = *inptr;
	    *outptr = 0;
/*
  Inside the following loop, the value of *inptr will be that of the first
  of succesive pairs of points.
*/
	    for(k=1; k<npts; k++) {
	      *outptr += 0.5 * (*inptr + *(inptr+inadd[0])) * (xptr[k] - xptr[k-1]);
	      inptr += inadd[0];
	    };
/*
  Complete as a wrapped integral if requested.
*/
	    if(do_wrap) {
	      *outptr += 0.5 * (*inptr + first_y) * wrap_len;
	      *outptr /= wrap_per;
	    };
/*
  Next integral.
*/
	    outptr++;
	    inptr += inadd[0] + inadd[1];
	  };
	  inptr += inadd[2];
	};
	return no_error;
}

/*.......................................................................
  User function to return the median value in an array.
*/
static Template(median_fn)
{
        int i,npts,*indx;
/*
  Determine the total number of elements in the input array.
*/
	npts=1;
	for(i=0;i<3;i++)
	  npts *= invals[0]->adim[i];
/*
  The standard procedure used to find the median is to sort
  the array into ascending order and then pick out the value
  in the middle element of the resulting array. In our case
  we will not both going as far as producing a sorted array,
  but can simply form the index array that would be used for
  the actual sort and deference the required value through
  its middle element.
*/
	if(indexx(npts, FLTPTR(invals[0]), &indx) == -1)
	  return -1;
/*
  Pick out the middle value.
*/
	*FLTPTR(outvals) = FLTPTR(invals[0])[indx[npts/2]];
	free(indx);
	return no_error;
}

/*.......................................................................
  Evaluate the correlation coefficient between two arrays.
*/
static Template(correl_fn)
{
        float x_sdev,y_sdev, x_mean, y_mean, xy_cov,temp,num;
	float *xptr, *yptr;
	int i, npts;
/*
  Signal an error if the two arrays have different dimensions.
*/
	npts=1;
	for(i=0; i<3; i++) {
	  npts *= invals[0]->adim[i];
	  if(invals[0]->adim[i] != invals[1]->adim[i]) {
	    lprintf(stderr, "correl(): Differing array dimensions.\n");
	    return -1;
	  };
	};
/*
  Determine the means of the two arrays separately.
*/
	mean_fn(invals, 1, outvals);
	x_mean = *FLTPTR(outvals);
	mean_fn(&invals[1], 1, outvals);
	y_mean = *FLTPTR(outvals);
/*
  Get the pointer to the first element.
*/
        xptr = FLTPTR(invals[0]);
        yptr = FLTPTR(invals[1]);
/*
  Determine the standard deviation of each array individually.
  Use a running mean to avoid precision problems.
*/
        x_sdev = y_sdev = 0;
	num = 0.0;
        for(i=1; i<=npts; i++) {
	  num += 1;
	  temp = *xptr - x_mean;
          x_sdev += (temp*temp-x_sdev)/num;
	  xptr++;
	  temp = *yptr - y_mean;
          y_sdev += (temp*temp-y_sdev)/num;
	  yptr++;
        };
/*
  Square root the mean square difference to get the standard deviation.
*/
        x_sdev = sqrt(x_sdev);
        y_sdev = sqrt(y_sdev);
/*
  If either standard deviation is zero signal an error.
*/
	if(x_sdev == 0.0 || y_sdev == 0.0) {
	  lprintf(stderr, "correl(): Zero standard deviation encounterred.\n");
	  return -1;
	};
/*
  Now determine the covariance.
  Get the pointer to the first element.
*/
        xptr = FLTPTR(invals[0]);
        yptr = FLTPTR(invals[1]);
/*
  Use a running mean to avoid precision problems.
*/
        xy_cov = 0;
	num = 0.0;
        for(i=1; i<=npts; i++) {
	  num += 1.0;
	  xy_cov += ((*xptr - x_mean) * (*yptr - y_mean) - xy_cov)/num;
	  xptr++; yptr++;
        };
/*
  The covariance divided by the product of the standard deviations
  is the required correlation coefficient.
*/
	*FLTPTR(outvals) = xy_cov / (x_sdev * y_sdev);
	return no_error;
}

/*.......................................................................
  Given two 1D data arrays and an optional weight array, perform a least
  squares fit for the straight line that they represent.
*/
static Template(fitline_fn)
{
        int npts,i;
	float mx,my,mxy,mxx,weight,wsum,nsum,scale,xval,yval;
	float *xptr, *yptr, *wptr;
/*
  Make sure that all the arrays have the same number of elements.
*/
	npts = invals[0]->adim[0];
	for(i=1; i<npar; i++)
	  if(invals[i]->adim[0] != npts) {
	    lprintf(stderr, "fit_line: Differing input array sizes.\n");
	    return -1;
	  };
/*
  Get pointers to the data arrays and the weight array.
*/
	xptr = FLTPTR(invals[0]);
	yptr = FLTPTR(invals[1]);
	wptr = (npar > 2) ? FLTPTR(invals[0]) : NULL;
/*
  We require weighted mean values of xy,x,y and xx.
  Start by finding running means of wxy,wx,wy, and wxx
  and then scale them by nsum/wsum to recover the weighted
  means.
*/
	mx = my = mxy = mxx = wsum = nsum = 0.0;
	for(i=0; i<npts; i++) {
	  weight = (wptr == NULL)? 1.0 : *(wptr++);
	  xval = *(xptr++);
	  yval = *(yptr++);
	  wsum += weight;
	  nsum += 1.0;
	  mx  += (weight*xval - mx)/nsum;
	  my  += (weight*yval - my)/nsum;
	  mxx += (weight*xval*xval - mxx)/nsum;
	  mxy += (weight*xval*yval - mxy)/nsum;
	};
/*
  Recover weighted means.
*/
	scale = nsum/wsum;
	mx *= scale;
	my *= scale;
	mxx *= scale;
	mxy *= scale;
/*
  Check for inifinite gradient.
*/
	if(mxx - mx*mx <= 1e-30) {
	  lprintf(stderr, "fit_line: infinite gradient found\n");
	  return -1;
	};
/*
  Determine the required gradient.
*/
	grad = (mxy - mx * my)/(mxx - mx * mx);
/*
  Determine the y-intercept.
*/
	yint = my - grad * mx;
/*
  Determine the standard deviation of the gradient and intercept.
*/
	graderr = sqrt(fabs(1/(wsum * (mxx - mx * mx))));
	yinterr = graderr * mxx;
	return no_error;
}

/*.......................................................................
  Regrid a data array onto a new regular coordinate grid. The original
  coordinate array is sent as the first argument. The array to be
  interpolated onto the new grid is sent as the third, preceded by a number
  to specify the axis of the array to be interpolated. The new number of
  channels is specified in the forth argument. The regridded array
  is returned as the value of the function.
*/
static Template(grid_fn)
{
	int axis;	/* The axis to be regridded */
	int nold;	/* The required number of channels on the old grid */
	int nnew;	/* The required number of channels on the new grid */
	int olddim[3];	/* No. of elements on each axis of input array */
	int newdim[3];	/* No. of elements on each axis of return array */
	int order[3];	/* Axis order in which to step through in/out arrays */
	int oldadd[3];	/* Element increments to step through input array */
	int newadd[3];	/* Element increments to step through return array */
	int nj,nk;	/* Nested for-loop variables for stepping along axes */
	float *oldptr;	/* Pointer to original array */
	float *newptr;	/* Pointer to new array */
	float *old_grid;	/* Pointer to orignal time grid array */
	float *new_grid;	/* Pointer to new time grid array */
	float start,inc;	/* New grid parameters */
	float newpos;	/* A new grid position. */
	int ach,bch;    /* Element numbers of elements that straddle newpos */
	float frac;     /* Fractional displacement of newpos between two
			   two stradling old grid positions */
	float *aptr,*bptr; /* Pointers to elements ach and bch of old array */
	int axis_inc;   /* Object increment per element along 'axis' */
	int i,j,k;
/*
  Get the axis designation number - must be between 0 and 2.
*/
	axis = *INTPTR(invals[0]);
	if(axis < 0 || axis > 2) {
	  lprintf(stderr, "Non-existent axis: %d\n", axis);
	  return -1;
	};
/*
  Get the original coordinate grid and make sure it has more than one element.
*/
	old_grid = FLTPTR(invals[1]);
	nold = invals[1]->adim[0];
	if(nold == 1) {
	  lprintf(stderr, "Can\'t re-grid a one element array!\n");
	  return -1;
	};
/*
  Check that the orginal coordinate grid is in increasing order.
*/
	for(i=1; i<nold; i++) {
	  if(old_grid[i] - old_grid[i-1] < 1e-20) {
	    lprintf(stderr,"Old grid not in ascending order.\n");
	    return -1;
	  };
	};
/*
  Get a pointer to the original data array.
*/
	oldptr = FLTPTR(invals[2]);
/*
  Get a copy of the input array dimensions in olddim[].
*/
	for(i=0; i<3; i++)
	  olddim[i] = invals[2]->adim[i];
/*
  Check that along axis 'axis' it has the same number of elements as
  the coordinate array.
*/
	if(invals[2]->adim[axis] != nold) {
	  lprintf(stderr, "The original axis grid and data arrays have conflicting sizes\n");
	  return -1;
	};
/*
  Get a pointer to the array for the new grid. Also get the required
  number of elements in the new grid and check it.
*/
	new_grid = FLTPTR(invals[3]);
	nnew = invals[3]->adim[0];
	if(nnew < 2) {
	  lprintf(stderr, "Can't interpolate onto a %d point grid\n",nnew);
	  return -1;
	};
/*
  Since old_grid is in ascending order, the minimum and maximum values
  of that grid are at the start and end of old_grid. The new grid
  will also have these limits - set up start and inc (increment) to
  reflect the start value and increment from point to point in the new
  grid.  
*/
	start = old_grid[0];
	inc = (old_grid[nold-1] - start) / (nnew-1);
/*
  Finally get a pointer to the start of the array to receive the re-gridded
  array and check that it matches the size of the original array along
  all axes except 'axis'. Also check the number of points along the
  re-grid axis equals that in the new grid array.
*/
	newptr = FLTPTR(invals[4]);
	for(i=0; i<3; i++) {
	  newdim[i] = invals[4]->adim[i];
	  if(i==axis) {
	    if(newdim[i] != nnew) {
	      lprintf(stderr, "New grid array and data arrays differ in size\n");
	      return -1;
	    };
	  } else if(newdim[i] != olddim[i]) {
	    lprintf(stderr, "New and old arrays differ in shape\n");
	    return -1;
	  };
	};
/*
  All OK. In order to treat every point in the data array we shall employ
  nested for-loops one per axis. The outer loop will be that for the
  interpolation axis since we only want to have to search the original
  grid array once per new grid position to find the two points straddling
  the new position. Record this order as axis numbers in array order[3].
*/
	for(i=0,j=0; i<3; i++)
	  if(i != axis) order[j++] = i;
	order[j]=axis;
/*
  Get the increment in floats required to step by one element along each
  axis in the order described in order[] at the end of each iteration of the
  respective nested for-loop. Do this for both the input and return arrays.
*/
	get_increments(order,olddim,oldadd);
	get_increments(order,newdim,newadd);
/*
  Also determine the increment required to step by one element along the
  interpolation axis of the data array without including the affects
  of the inner nested for-loops. This is neccessary for absolute
  accessing of elements rather than relative to the last iteration
  of a for-loop.
*/
	axis_inc = 1;
	for(i=0; i<axis; i++)
	  axis_inc *= olddim[i];
/*
  Perform the interpolated regridding. For each of the new grid positions,
  the original grid will be searched for the two points that straddle
  that position. Since the grids are in scending order the search will
  not start from the beginning of the grid each time but from the last
  two points. These will be held in ach and bch. Set them up for the
  first point.
*/
	ach=0;
	bch=1;
/*
  Pre-calculate the upper limits of the inner nested for-loops.
*/
	nj = newdim[order[1]];
	nk = newdim[order[0]];
	for(i=0; i<nnew; i++) {
/*
  Determine position of the next of the new grid points.
*/
	  newpos = start + inc * i;
	  new_grid[i] = newpos;
/*
  Search for the element numbers, ach & bch of the two points that
  straddle it in the old grid. NB.
*/
	  for(;;) {
	    if( (old_grid[ach] <= newpos && old_grid[bch] >= newpos) ||
	       bch == nold-1)
	      break;
	    ach++;
	    bch++;
	  };
/*
  Find the fraction of the distance between the two old grid points
  to the new grid point.
*/
	  frac = (newpos-old_grid[ach]) / (old_grid[bch]-old_grid[ach]);
/*
  Set aptr and bptr to point to the 1st elements of the old array
  at channels ach and bch on axis 'axis'.
*/
	  aptr = oldptr + ach * axis_inc;
	  bptr = oldptr + bch * axis_inc;
/*
  Loop for each element along other axes at the current channels ach and bch.
*/
	  for(j=0; j<nj; j++) {
/*
  Inside the following loop, the value of *inptr will be that of the second
  point of succesive pairs of points - so pre-skip to the second point.
*/
	    for(k=0; k<nk; k++) {
/*
  Inner loop - perform interpolation.
*/
	      *newptr = *aptr + (*bptr - *aptr) * frac;
	      newptr += newadd[0];
	      aptr += oldadd[0];
	      bptr += oldadd[0];
	    };
	    newptr += newadd[1];
	    aptr += oldadd[1];
	    bptr += oldadd[1];
	  };
	  newptr += newadd[2];
	};
	return no_error;
}

/*.......................................................................
 * Return an array which samples a guassian out to a given number of
 * sigma.
 *
 * Input:
 *  fwhm    float  The full-width-half-maximum of the gaussian.
 *  nsigma  float  The number of sigma out to which to sample the
 *                 gaussian.
 *  step    float  The x-axis increment between samples.
 */
static Template(garray_fn)
{
  float fwhm;    /* The FWHM of the Gaussian */
  float nsigma;  /* The number of sigma out to which to sample the gaussian */
  float step;    /* The increment between samples. */
  float sigma;   /* The standard deviation of the Gaussian */
  int nvals;     /* The number of elements in the returned array */
  int i;
/*
 * Copy the user arguments into local variables.
 */
  fwhm = *FLTPTR(invals[0]);
  nsigma = *FLTPTR(invals[1]);
  step = *FLTPTR(invals[2]);
/*
 * Validate the arguments.
 */
  if(fwhm <= 0.0 || nsigma <= 0.0 || step <= 0.0) {
    lprintf(stderr, "gaussian_array: Invalid negative value(s).\n");
    return -1;
  };
/*
 * Convert the FWHM to a standard deviation.
 */
  sigma = FWHM_TO_STDDEV * fwhm;
/*
 * Determine the number of elements that the return array will require.
 */
  nvals = (sigma * nsigma) / step + 1;
/*
 * Allocate memory for the return array.
 */
  if((VOIDPTR(outvals)=valof_alloc(nvals,'f')) == NULL)
    return -1;
/*
 * Fill in the dimensional aspects of the return descriptor.
 */
  outvals->num_el  = nvals;
  outvals->adim[0] = nvals;
/*
 * Compute the guassian at each of the x-axis values needed by the
 * sample array.
 */
  for(i=0; i<nvals; i++) {
    float x = (step * i) / sigma;
    FLTPTR(outvals)[i] = exp(-0.5 * x*x);
  };
  return no_error;
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "obs.h"
#include "logio.h"
#include "visaver.h"

static Scatsum *new_Scatsum(long nvis);               /* Constructor */
static Scatsum *del_Scatsum(Scatsum *scatsum);        /* Destructor */
static void ini_Scatsum(Scatsum *scatsum, long nvis); /* Initialize sums */

static Basesum *new_Basesum(int nbmax);            /* Constructor */
static Basesum *del_Basesum(Basesum *bsum);        /* Destructor */
static void ini_Basesum(Basesum *bsum, int nbmax); /* Initialize sums */


/*.......................................................................
 * Allocate and return a Visaver class.
 *
 *  dp         Dpage *  The descriptor of the output file for the averaged
 *                      visibilities.
 *  avtime    double    The solution bin width (seconds).
 *  scatter      int    If true, allocate an extra buffer to use in finding
 *                      uncertainties from the scatter of the visibilities.
 * Output:
 *  return   Visaver *  The new visaver class instance, or NULL on error.
 */
Visaver *new_Visaver(Dpage *dp, double avtime, int scatter)
{
  Visaver *av = NULL;  /* The pointer to the new visaver instance */
/*
 * Allocate the container.
 */
  av = malloc(sizeof(Visaver));
  if(av==NULL) {
    lprintf(stderr, "new_Visaver: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize at least up to the point at which del_Visaver() can safely be
 * called.
 */
  av->nvis = dp->nvis;
  av->nbmax = dp->nbase;
  av->nbase = 0;
  av->vis = NULL;
  av->scatsum = NULL;
  av->basesum = NULL;
  av->dp = dp;
/*
 * Allocate the scatter sum array if required.
 */
  if(scatter) {
    av->scatsum = new_Scatsum(av->nvis);
    if(av->scatsum==NULL)
      return del_Visaver(av);
  };
/*
 * Allocate the baseline sum array.
 */
  av->basesum = new_Basesum(av->nbmax);
  if(av->basesum==NULL)
    return del_Visaver(av);
/*
 * Return the new visaver instance.
 */
  return av;
}

/*.......................................................................
 * Delete a visaver class instance.
 *
 * Input:
 *  av     Visaver *  The visaver instance to be deleted.
 * Output:
 *  return Visaver *  The deleted instance, ie NULL. Use like:
 *                    av = del_Visaver(av);
 */
Visaver *del_Visaver(Visaver *av)
{
  if(av) {
/*
 * Delete the contents.
 */
    del_Scatsum(av->scatsum);
    del_Basesum(av->basesum);
/*
 * Delete the container.
 */
    free(av);
  };
  return NULL;
}

/*.......................................................................
 * Allocate an array of container classes used to store intermediate
 * sums used in determining visibility weights from their scatter within
 * each averaging solution bin. ini_Scatsum() should be called before
 * each use of the returned array.
 *
 * Input:
 *  nvis     long    The max number of visibilities in an integration
 *                   record.
 * Output:
 *  return Scatsum * The pointer to the first elements of the new array,
 *                   or NULL on error.
 */
static Scatsum *new_Scatsum(long nvis)
{
  Scatsum *ssum;   /* The pointer to the array to be returned */
/*
 * Sanity checks.
 */
  if(nvis <= 0) {
    lprintf(stderr, "new_Scatsum: No visibilities.\n");
    return NULL;
  };
/*
 * Allocate the array.
 */
  ssum = malloc(sizeof(Scatsum) * nvis);
  if(ssum==NULL) {
    lprintf(stderr, "new_Scatsum: Insufficient memory for sums.\n");
    return NULL;
  };
/*
 * Return the array.
 */
  return ssum;
}

/*.......................................................................
 * Delete an array of scatter sums.
 *
 * Input:
 *  ssum    Scatsum *  The pointer to the first element of the array to
 *                      be deleted.
 * Output:
 *  return  Scatsum *  The deleted array, ie. NULL.
 */
static Scatsum *del_Scatsum(Scatsum *ssum)
{
  if(ssum)
    free(ssum);
  return NULL;
}

/*.......................................................................
 * Initialize the elements of a scatter estimate array.
 *
 * Input:
 *  ssum  Scatsum *  The pointer to the first element of the array to
 *                   be initialized.
 *  nvis     long    The number of elements to be initialized.
 */
static void ini_Scatsum(Scatsum *ssum, long nvis)
{
  static const Scatsum clr_stat = {0, 0.0f};
  long i;
/*
 * Initialize the elements of the array.
 */
  if(ssum) {
    for(i=0; i<nvis; i++,ssum++)
      *ssum = clr_stat;
  };
  return;
}

/*.......................................................................
 * Allocate an array of container classes used to store intermediate
 * sums used in determining weighted baseline U,V,W and UT.
 * ini_Basesum() should be called before each use of the returned array.
 *
 * Input:
 *  nbmax      int   The max number of baselines in any sub-array.
 * Output:
 *  return Basesum * The pointer to the first elements of the new array,
 *                    or NULL on error.
 */
static Basesum *new_Basesum(int nbmax)
{
  Basesum *bsum;   /* The pointer to the array to be returned */
/*
 * Sanity checks.
 */
  if(nbmax <= 0) {
    lprintf(stderr, "new_Basesum: Zero baselines?\n");
    return NULL;
  };
/*
 * Allocate the array.
 */
  bsum = malloc(sizeof(Basesum) * nbmax);
  if(bsum==NULL) {
    lprintf(stderr, "new_Basesum: Insufficient memory for sums.\n");
    return NULL;
  };
/*
 * Return the array.
 */
  return bsum;
}

/*.......................................................................
 * Delete an array of baseline sums.
 *
 * Input:
 *  bsum    Basesum *  The pointer to the first element of the array to
 *                     be deleted.
 * Output:
 *  return  Basesum *  The deleted array, ie. NULL.
 */
static Basesum *del_Basesum(Basesum *bsum)
{
  if(bsum)
    free(bsum);
  return NULL;
}

/*.......................................................................
 * Initialize the elements of a baseline sum array.
 *
 * Input:
 *  bsum  Basesum *  The pointer to the first element of the array to
 *                   be initialized.
 *  nbmax     int    The number of elements to be initialized.
 */
static void ini_Basesum(Basesum *bsum, int nbmax)
{
  static const Basesum clr_bsum = {0.0f};
  long i;
/*
 * Initialize the elements of the array.
 */
  for(i=0; i<nbmax; i++,bsum++)
    *bsum = clr_bsum;
  return;
}

/*.......................................................................
 * Reset output visibilities and averaging sums perparatory to starting
 * a new averaging bin.
 *
 * Input:
 *  av      Visaver *  The visibility average descriptor.
 *  vis  Visibility *  The array of nbase baselines visibilities in
 *                     which running mean U,V,W are to be accumulated.
 *                     (No other members will be touched).
 *  nbase       int    The number of baselines in 'vis[]'.
 *  irec       long    The output integration record number.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
int av_newint(Visaver *av, Visibility *vis, int nbase, long irec)
{
/*
 * Check args.
 */
  if(av==NULL || vis==NULL) {
    lprintf(stderr, "av_newint: NULL %s intercepted.\n",
	    av==NULL ? "Visaver" : "Visibility");
    return 1;
  };
/*
 * Clear the accumulated visibility sums of the last solution in preparation
 * for the next solution.
 */
  ini_Scatsum(av->scatsum, av->nvis);
  ini_Basesum(av->basesum, av->nbmax);
/*
 * Clear the output buffer.
 */
  dp_clear(av->dp, irec);
/*
 * Install the new visibility array.
 */
  av->vis = vis;
  av->nbase = nbase;
  return 0;
}

/*.......................................................................
 * Complete the job of averaging an output integration started by
 * av_newint(). This currently does two things:
 *
 * 1. It clears the U,V,W coordinates of unsampled visibilities in the
 *    array of visibility descriptors previously provided to av_newint().
 *
 * 2. It uses the scatter statistics that were accumulated to calculate
 *    and record the weights for each output (dp->cvis[]) visibility.
 *
 * Input:
 *  av      Visaver *  The visibility average descriptor.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
int av_endint(Visaver *av)
{
/*
 * Check args.
 */
  if(av==NULL) {
    lprintf(stderr, "av_endint: NULL Visaver intercepted.\n");
    return 1;
  };
/*
 * Clear the U,V,W coordinates of un-sampled visibilities.
 */
  if(av->vis) {
    Visibility *vis = av->vis;
    Basesum *bsum = av->basesum;
    int base;
    for(base=0; base<av->nbase; base++,vis++,bsum++) {
      if(bsum->wtsum==0.0f)
	vis->u = vis->v = vis->w = 0.0f;
    };
  };
/*
 * Calculate and record the output weights deduced from the scatter
 * sums.
 */
  if(av->scatsum) {
    const float tiny = 1.0e-16;     /* Minimum variance allowed */
    Cvis *cvis = av->dp->cvis;      /* Output integration record */
    Scatsum *scatsum = av->scatsum; /* The accumulated scatter stats */
    long iv;                        /* Index into cvis[] */
    for(iv=0; iv<av->nvis; iv++,cvis++,scatsum++) {
/*
 * Ignore deleted visibilities.
 */
      if(cvis->wt != 0.0f) {
	float new_wt;   /* The new weight to be assigned to the visibility */
/*
 * If there were insufficient points to provide a scatter estimate,
 * leave the original weight in place but make it -ve to flag the point.
 */
	if(scatsum->nsum < 2) {
	  new_wt = -fabs(cvis->wt);
	} else {
	  float im = cvis->im;     /* Imaginary part of mean visibility */
	  float re = cvis->re;     /* Real part of mean visibility */
	  float n = scatsum->nsum; /* The number of visibilities in the sums */
/*
 * Determine the amplitude variance from the data scatter. Note that this
 * is the variance of the mean, which is defined as the variance of the
 * sampled data divided by the number of points in the sample. Furthermore
 * since we want an unbiased estimate of the data variance, we correct
 * the variance (derived from weighted means) by a factor n/(n-1).
 */
	  float variance = 0.5f * (scatsum->sqr_mean - im*im - re*re) / (n-1.0);
/*
 * The weight to give the visibility is the reciprocal of the amplitude
 * variance. Take care to prevent divide by zero.
 */
	  new_wt = 1.0f / (variance>tiny ? variance:tiny);
	};
/*
 * Take care not to un-flag a flagged visibility when assigning the new weight.
 */
	cvis->wt = cvis->wt > 0.0f ? new_wt : -fabs(new_wt);
      };
    };
  };
/*
 * Vector averaging can yield 0 amplitudes from good input visibilities.
 * If this happens, flag the affected visibility.
 */
  {
    Cvis *cvis = av->dp->cvis;      /* Output integration record */
    long iv;                        /* Index into cvis[] */
    for(iv=0; iv<av->nvis; iv++,cvis++) {
      if(cvis->wt!=0.0f && (cvis->re==0.0f && cvis->im==0.0f))
	cvis->wt = 0.0f;
    };
  };
  return 0;
}

/*.......................................................................
 * Include a new visibility in the running averages of the output
 * integration record.
 *
 * To preserve precision, running weighted means will be used. Flagged
 * data will only be included in these means until supplanted by the
 * un-flagged visibilities.
 *
 * Input:
 *  av   Visaver *   The visibility average descriptor.
 *  re     float     The real part of the input visibility.
 *  im     float     The imaginary part of the input visibility.
 *  wt     float     The weight of the input visibility.
 *  ivis    long     The index of the visibility in the output
 *                   (Dpage *dp->cvis[]) integration record.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int av_dp(Visaver *av, float re, float im, float wt, long ivis)
{
/*
 * Consider all except deleted visibilities.
 */
  if(wt!=0.0f) {
    Cvis *cvis = &av->dp->cvis[ivis];
/*
 * Initialize the output running means when the first flagged or
 * unflagged visibility is encountered and again when the first un-flagged
 * visibility is encountered.
 */
    if(cvis->wt==0.0f || (cvis->wt<0.0f && wt>0.0f)) {
      cvis->wt = wt;
      cvis->im = im;
      cvis->re = re;
/*
 * Optionally accumulate running means for use in estimates of the data scatter.
 */
      if(av->scatsum) {
	Scatsum *scatsum = &av->scatsum[ivis];
	scatsum->sqr_mean = (re*re+im*im);
	scatsum->nsum = 1;
      };
    }
/*
 * Accumulate the running mean real and imaginary parts of the visibility,
 * using flagged visibilities only until the flagged means are reset
 * (see above) to un-flagged by the first good visibility.
 */
    else if(wt > 0.0f || cvis->wt < 0.0f) {
      float runwt = wt / (cvis->wt += wt);
      cvis->im += runwt * (im - cvis->im);
      cvis->re += runwt * (re - cvis->re);
/*
 * Optionally accumulate running means for use in estimates of the data scatter.
 */
      if(av->scatsum) {
	Scatsum *scatsum = &av->scatsum[ivis];
	scatsum->sqr_mean += runwt * (re*re+im*im - scatsum->sqr_mean);
	scatsum->nsum++;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Include a new visibility in the running mean of baseline U,V,W
 * coordinates, and increment the integration time accordingly.
 *
 * To preserve precision, running weighted means will be used. Flagged
 * data will only be included in these means until supplanted by
 * un-flagged visibilities.
 *
 * Input:
 *  av   Visaver *   The visibility average descriptor.
 *  uu     float     The input U coordinate.
 *  vv     float     The input V coordinate.
 *  ww     float     The input W coordinate.
 *  wt     float     The weight of the input visibility.
 *  inttim float     The integration time of the visibility, or 0.0 if
 *                   unknown.
 *  base     int     The index of the baseline to accumulate the means
 *                   for.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int av_uvwt(Visaver *av, float uu, float vv, float ww, float wt, float inttim,
	    int base)
{
/*
 * Consider all except deleted visibilities.
 */
  if(wt!=0.0f) {
    Basesum *bsum = &av->basesum[base];
    Visibility *vis = &av->vis[base];
/*
 * Initialize the output U,V,W running means when the first flagged or
 * unflagged visibility is encountered and again when the first un-flagged
 * visibility is encountered.
 */
    if(bsum->wtsum==0.0f || (bsum->wtsum<0.0f && wt>0.0f)) {
      bsum->wtsum = wt;
      vis->u = uu;
      vis->v = vv;
      vis->w = ww;
      vis->dt = inttim;
    }
/*
 * Accumulate the running mean U,V,W, using flagged visibilities only until
 * the flagged means are reset (see above) to un-flagged by the first good
 * visibility.
 */
    else if(wt > 0.0f || bsum->wtsum < 0.0f) {
      float runwt = wt / (bsum->wtsum += wt);
      vis->u += runwt * (uu - vis->u);
      vis->v += runwt * (vv - vis->v);
      vis->w += runwt * (ww - vis->w);
      vis->dt += inttim;
    };
  };
  return 0;
}

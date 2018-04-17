#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "logio.h"
#include "obs.h"
#include "vlbconst.h"
#include "telspec.h"
#include "baselist.h"
#include "spectra.h"
#include "obedit.h"

static int dp_sumspec(Observation *ob, Spectrum *spec);
static int bad_Spectrum(Spectrum *spec, char *fname);
static int bad_Spectra(Spectra *spectra, char *fname);

/*.......................................................................
 * Create a spectrum list container for a given observation.
 *
 * Input:
 *  ob   Observation *   The observation to which the list pertains.
 * Output:
 *  return   Spectra *   The container of the list of spectra, or NULL
 *                       on error.
 */
Spectra *new_Spectra(Observation *ob)
{
  Spectra *spectra;  /* The container to be returned */
/*
 * We need baseline indexes, integration time-stamps, UV radii etc..
 */
  if(!ob_ready(ob, OB_INDEX, "new_Spectra"))
    return NULL;
/*
 * Allocate the new container.
 */
  spectra = (Spectra *) malloc(sizeof(Spectra));
  if(spectra == NULL) {
    lprintf(stderr, "new_Spectra: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize the container at least up to the point at which it can
 * safely be sent to del_Spectra().
 */
  spectra->ob = ob;
  spectra->head = NULL;
  spectra->tail = NULL;
  return spectra;
}

/*.......................................................................
 * Delete a list of spectra previously allocated with new_Spectra().
 *
 * Input:
 *  spectra   Spectra *  The list to be deleted.
 * Output:
 *  return    Spectra *  Allways NULL.
 */
Spectra *del_Spectra(Spectra *spectra)
{
  if(spectra) {
    if(spectra->head) {
      Spectrum *head = NULL;
      Spectrum *next = spectra->head;
      while(next) {
	head = next;
	next = head->next;
	head = del_Spectrum(head);
      };
    };
    free(spectra);
  };
  return NULL;
}

/*.......................................................................
 * Create a new spectrum container and add it to a spectrum list.
 * Note that the spectrum itself will not be filled until
 * get_Spectra(spectra) is called.
 *
 * Input:
 *  spectra    Spectra *  The spectrum list to add to.
 *  dovector       int    0 - Scalar averaged real_complex=amplitude spectrum.
 *                        1 - Vector averaged complex spectrum.
 *  stokes      Stokes    The polarization of the spectrum.
 *  uta            int    The ob->rec[] index of the first integration to
 *                        sample.
 *  utb            int    The ob->rec[] index of the last integration to
 *                        sample.
 *  uvmin        float    The minimum UV radius to sample.
 *  uvmax        float    The maximum UV radius to sample. For the
 *                        whole UV range, set uvmin=uvmax=0.
 *  bgrp       Basegrp *  A baseline group, or NULL if
 *                        this is to be defered to a later call to
 *                        spc_set_bgrp().
 * Output:
 *  return    Spectrum *  The new spectrum container, or NULL on
 *                        error.
 */
Spectrum *add_Spectrum(Spectra *spectra, int dovector, Stokes stokes,
		       int uta, int utb, float uvmin, float uvmax,
		       Basegrp *bgrp)
{
  Observation *ob;  /* The observation of the spectrum */
  Spectrum *spec;   /* The container to be returned */
  int cif;          /* An IF index */
  int isub;         /* Sub-array index */
  int nbtotal;      /* The total number of baselines selected */
/*
 * Check that we have a list to add to.
 */
  if(!spectra) {
    lprintf(stderr, "add_Spectrum: No list of spectra to add to.\n");
    return NULL;
  };
/*
 * Check that we have a list of baselines to be plotted.
 */
  if(!bgrp) {
    lprintf(stderr, "add_Spectrum: No baseline specification list provided.\n");
    return NULL;
  };
/*
 * Get the descriptor of the observation associated with the spectrum.
 */
  ob = spectra->ob;
/*
 * Create a new spectrum container.
 */
  spec = (Spectrum *) malloc(sizeof(Spectrum));
  if(!spec) {    
    lprintf(stderr, "add_Spectrum: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize the container at least to the point at which it is safe
 * to call del_Spectrum(spec) before doing anything that might fail.
 */
  spec->uta = uta;
  spec->utb = utb;
  spec->uvmin = 0.0f;
  spec->uvmax = 0.0f;
  spec->ssub = NULL;
  spec->nsub = ob->nsub;
  spec->ifs = NULL;
  spec->nif = ob->nif;
  spec->dovector = dovector;
  spec->next = NULL;
/*
 * Assign the required polarization.
 */
  if(spc_set_pol(spectra, spec, stokes))
    return del_Spectrum(spec);
/*
 * Set the start and end integration indexes.
 */
  if(spc_set_ut(spectra, spec, uta, utb))
    return del_Spectrum(spec);
/*
 * Set the intial UV range to encompass all UV radii.
 */
  if(spc_set_uvrange(spectra, spec, uvmin, uvmax))
    return del_Spectrum(spec);
/*
 * Allocate the array of sub-array specific lists of baselines.
 */
  spec->ssub = (Specsub *) malloc(sizeof(Specsub) * spec->nsub);
  if(!spec->ssub) {
    lprintf(stderr, "add_Spectrum: Insufficient memory.\n");
    return del_Spectrum(spec);
  };
/*
 * Initialize the array up to the point at which spec is again a viable
 * argument to del_Spectrum().
 */
  nbtotal = 0;
  for(isub=0; isub<spec->nsub; isub++) {
    Specsub *ssub = spec->ssub + isub;
    ssub->nbmax = ob->sub[isub].nbase;
    ssub->nbase = 0;
    ssub->baselines = NULL;
    nbtotal += ssub->nbmax;  /* Count the total number of baselines */
  };
/*
 * Allocate all the arrays of baselines as a single array, anchored in
 * the first sub-array container.
 */
  spec->ssub->baselines = (int *) malloc(sizeof(int) * nbtotal);
  if(!spec->ssub->baselines) {
    lprintf(stderr, "add_Spectrum: Insufficient memory.\n");
    return del_Spectrum(spec);
  };
/*
 * Distribute the baseline array between sub-arrays.
 */
  for(isub=0; isub < spec->nsub - 1; isub++) {
    Specsub *ssub = spec->ssub + isub;
    ssub[1].baselines = ssub->baselines + ssub->nbmax;
  };
/*
 * Install the initial list of baselines.
 */
  if(spc_set_bgrp(spectra, spec, bgrp))
    return del_Spectrum(spec);
/*
 * Allocate the container of IF spectra.
 */
  spec->ifs = (IFspec *) malloc(sizeof(IFspec) * spec->nif);
  if(!spec->ifs) {
    lprintf(stderr, "add_Spectrum: Insufficient memory.\n");
    return del_Spectrum(spec);
  };
/*
 * Initialize each IF container at least to the point at which the
 * array can be sent to del_Spectrum() before doing anything else that
 * might fail.
 */
  for(cif=0; cif<spec->nif; cif++) {
    IFspec *ifs = spec->ifs + cif;
    ifs->chan = NULL;
    ifs->nchan = ob->nchan;
  };
/*
 * Now allocate separate spectrum arrays for each IF.
 */
  for(cif=0; cif<spec->nif; cif++) {
    IFspec *ifs = spec->ifs + cif;
    ifs->chan = (Cvis *) malloc(sizeof(Cvis) * ifs->nchan);
    if(!ifs->chan) {
      lprintf(stderr, "Insufficient memory.\n");
      return del_Spectrum(spec);
    };
  };
/*
 * Append the spectrum to the spectrum list.
 */
  if(spectra->tail==NULL) {
    spectra->head = spectra->tail = spec;
  } else {
    spectra->tail->next = spec;
    spectra->tail = spec;
  };
  return spec;
}

/*.......................................................................
 * Report an error if the given spectrum descriptor is NULL.
 *
 * Input:
 *  spec    Spectrum *  The spectrum to check.
 *  fname       char *  The name of the calling function.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Bad spectrum.
 */
static int bad_Spectrum(Spectrum *spec, char *fname)
{
  if(spec==NULL) {
    lprintf(stderr, "%s: NULL Spectrum descriptor.\n",
	    fname ? fname : "bad_Spectrum");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Report an error if the given spectrum list descriptor is NULL.
 *
 * Input:
 *  spectra Spectra *  The spectrum list to check.
 *  fname      char *  The name of the calling function.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Bad spectrum list.
 */
static int bad_Spectra(Spectra *spectra, char *fname)
{
  if(spectra==NULL) {
    lprintf(stderr, "%s: NULL Spectrum list descriptor.\n",
	    fname ? fname : "bad_Spectra");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Change the current baseline selection list of a given spectrum
 * that was previously allocated by add_Spectrum().
 *
 * Input:
 *  spectra Spectra *  The spectrum list containing spec.
 *  spec   Spectrum *  The spectrum whose baseline selection list is to
 *                     be changed.
 *  bgrp    Basegrp *  The new baseline group selection to be installed.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error (spec is left unchanged).
 */
int spc_set_bgrp(Spectra *spectra, Spectrum *spec, Basegrp *bgrp)
{
  Observation *ob;  /* The descriptor of the observation of the spectra */
  int isub;         /* Sub-array index */
/*
 * Check arguments.
 */
  if(bad_Spectra(spectra, "spc_set_bgrp") ||
     bad_Spectrum(spec, "spc_set_bgrp"))
    return 1;
/*
 * Get the descriptor of the observation.
 */
  ob = spectra->ob;
/*
 * Compile the new baseline selection list.
 */
  for(isub=0; isub<spec->nsub; isub++) {
    Specsub *ssub = spec->ssub + isub;
    Subarray *sub = ob->sub + isub;
    ssub->nbase = 0;
    if(bgrp && bgrp->nnode > 0) {
      int base;
      for(base=0; base < sub->nbase; base++) {
	if(in_Basegrp(ob, isub, base, bgrp))
	  ssub->baselines[ssub->nbase++] = base;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Change the current polarization selection of a given spectrum
 * that was previously allocated by add_Spectrum().
 *
 * Input:
 *  spectra Spectra *  The spectrum list containing spec.
 *  spec   Spectrum *  The spectrum whose polarization selection is to
 *                     be changed.
 *  pol      Stokes    The new stokes parameter.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error (spec is left unchanged).
 */
int spc_set_pol(Spectra *spectra, Spectrum *spec, Stokes pol)
{
  Obpol obpol;
/*
 * Check arguments.
 */
  if(bad_Spectra(spectra, "spc_set_pol") ||
     bad_Spectrum(spec, "spc_set_pol"))
    return 1;
  if(get_Obpol(spectra->ob, pol, 1, &obpol))
    return 1;
  spec->obpol = obpol;
  return 0;
}

/*.......................................................................
 * Change the current polarization selection of a given spectrum
 * that was previously allocated by add_Spectrum().
 *
 * Input:
 *  spectra Spectra *  The spectrum list containing spec.
 *  spec   Spectrum *  The spectrum whose polarization selection is to
 *                     be changed.
 *  uta         int    The ob->rec[] index of the first integration to
 *                     sample.
 *  utb         int    The ob->rec[] index of the last integration to
 *                     sample.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error (spec is left unchanged).
 */
int spc_set_ut(Spectra *spectra, Spectrum *spec, int uta, int utb)
{
/*
 * Check arguments.
 */
  if(bad_Spectra(spectra, "spc_set_ut") ||
     bad_Spectrum(spec, "spc_set_ut"))
    return 1;
/*
 * Swap uta and utb if necessary.
 */
  if(uta > utb) {
    int uttmp = uta;
    uta = utb;
    utb = uttmp;
  };
/*
 * Check the indexes.
 */
  if(uta < 0)
    uta = 0;
  if(utb > spectra->ob->nrec-1)
    utb = spectra->ob->nrec-1;
  if(uta > utb) {
    lprintf(stderr, "spc_set_ut: Bad time range.\n");
    return 1;
  };
/*
 * Record the new integration range.
 */
  spec->uta = uta;
  spec->utb = utb;
  return 0;
}

/*.......................................................................
 * Change the visibility time-averaging mode of a spectrum.
 *
 * Input:
 *  spectra Spectra *  The spectrum list containing spec.
 *  spec   Spectrum *  The spectrum whose polarization selection is to
 *                     be changed.
 *  dovector    int    0 - Scalar average.
 *                     1 - Vector average.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error (dovector is left unchanged).
 */
int spc_set_avmode(Spectra *spectra, Spectrum *spec, int dovector)
{
/*
 * Check arguments.
 */
  if(bad_Spectra(spectra, "spc_set_avmode") ||
     bad_Spectrum(spec, "spc_set_avmode"))
    return 1;
/*
 * Install the new mode.
 */
  spec->dovector = dovector;
  return 0;
}

/*.......................................................................
 * Change the current UV range selection of a given spectrum that was
 * previously allocated by add_Spectrum().
 *
 * Input:
 *  spectra Spectra *  The spectrum list containing spec.
 *  spec   Spectrum *  The spectrum whose polarization selection is to
 *                     be changed.
 *  uvmin     float    The minimum UV radius to use (wavelengths).
 *  uvmax     float    The maximum UV radius to use. If uvmin==uvmax==0
 *                     then all UV radii will be included.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error (spec is left unchanged).
 */
int spc_set_uvrange(Spectra *spectra, Spectrum *spec, float uvmin, float uvmax)
{
/*
 * Check arguments.
 */
  if(bad_Spectra(spectra, "spc_set_uvrange") ||
     bad_Spectrum(spec, "spc_set_uvrange"))
    return 1;
/*
 * Check the values.
 */
  if(uvmin < 0.0f)
    uvmin = 0.0f;
  if(uvmax < 0.0f)
    uvmax = 0.0f;
/*
 * Swap uvmin and uvmax if necessary.
 */
  if(uvmin > uvmax) {
    int uvtmp = uvmin;
    uvmin = uvmax;
    uvmax = uvtmp;
  };
/*
 * Record the new integration range.
 */
  spec->uvmin = uvmin;
  spec->uvmax = uvmax;
  return 0;
}

/*.......................................................................
 * Remove a given spectrum from a list of spectra.
 *
 * Input:
 *  spectra  Spectra *   The list of spectra to remove the spectrum
 *                       from.
 *  spec    Spectrum *   The spectrum to be removed.
 * Output:
 *  return  Spectrum *   The removed spectrum, or NULL if not found.
 */
Spectrum *rem_Spectrum(Spectra *spectra, Spectrum *spec)
{
  Spectrum *prev; /* The spectrum preceding the spectrum being compared */
  Spectrum *next; /* The spectrum being compared against spec */
  if(!spectra) {
    lprintf(stderr, "rem_Spectrum: No spectrum list.\n");
    return NULL;
  };
/*
 * Search for the spectrum in the list.
 */
  prev = NULL;
  next = spectra->head;
  while(next && next!=spec) {
    prev = next;
    next = next->next;
  };
/*
 * Not found?
 */
  if(!next) {
    lprintf(stderr, "rem_Spectrum: Spectrum not in list.\n");
    return NULL;
  };
/*
 * Remove the spectrum from the list and re-link the list around it.
 */
  if(prev)
    prev->next = spec->next;
  else
    spectra->head = spec->next;
  if(spectra->tail == spec)
    spectra->tail = prev;
  spec->next = NULL;
  return spec;
}

/*.......................................................................
 * Delete a spectrum previously allocated by add_Spectrum(). Note that
 * the spectrum must first have been removed from the list with
 * rem_Spectrum() before calling this function.
 *
 * Input:
 *  spec   Spectrum *  The spectrum to be deleted (NULL is allowed).
 * Output:
 *  return Spectrum *  Allways NULL.
 */
Spectrum *del_Spectrum(Spectrum *spec)
{
  if(spec) {
/*
 * Delete any sub-array baseline lists.
 */
    if(spec->ssub) {
      if(spec->ssub->baselines)
	free(spec->ssub->baselines);
      free(spec->ssub);
    };
/*
 * Delete any per-IF spectra and their container.
 */
    if(spec->ifs) {
      int cif;
      for(cif=0; cif < spec->nif; cif++) {
	IFspec *ifs = spec->ifs + cif;
	if(ifs->chan)
	  free(ifs->chan);
      };
      free(spec->ifs);
    };
  };
  return NULL;
}

/*.......................................................................
 * Clear the per-IF spectrum arrays of a given spectrum container.
 *
 * Input:
 *  spec   Spectrum *   The spectrum whose arrays are to be zeroed.
 */
void clr_Spectrum(Spectrum *spec)
{
  Cvis zero_cvis = {0.0f, 0.0f, 0.0f};
  if(spec) {
    int fc;  /* Frequency channel index */
    int cif; /* The index of the IF being processed */
/*
 * Clear the spectrum of each IF.
 */
    for(cif=0; cif<spec->nif; cif++) {
      IFspec *ifs = spec->ifs + cif;
      for(fc=0; fc<ifs->nchan; fc++)
	ifs->chan[fc] = zero_cvis;
    };
  };
  return;
}

/*.......................................................................
 * Construct given spectra from the raw data in an ob->dp paging file.
 *
 * Input:
 *  spectra  Spectra *  The list of spectra to fill.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int get_Spectra(Spectra *spectra)
{
  Observation *ob;  /* The observation descriptor of the spectra */
  Spectrum *spec;   /* The spectrum being processed */
  int irec;         /* The index of the integration record being processed */
/*
 * NULL container?
 */
  if(bad_Spectra(spectra, "get_Spectra"))
    return 1;
/*
 * Are there no spectra to be built?
 */
  if(!spectra->head)
    return 0;
/*
 * Get the observation descriptor for the spectra.
 */
  ob = spectra->ob;
/*
 * Flush cached edits so that the UV data scratch file is up to date.
 */
  if(ed_flush(ob))
    return 1;
/*
 * Clear the spectrum arrays.
 */
  for(spec=spectra->head; spec; spec=spec->next)
    clr_Spectrum(spec);
/*
 * Initialize to read whole integrations.
 */
  if(dp_crange(ob->dp, 0, ob->nchan-1) ||
     dp_irange(ob->dp, 0, ob->nif-1)   ||
     dp_brange(ob->dp, 0, ob->nbmax-1) ||
     dp_srange(ob->dp, 0, ob->npol-1))
    return 1;
/*
 * Read sampled integration records and construct their spectra.
 */
  irec = 0;
  do {
/*
 * Find the spectrum having the earliest start integration record index
 * that is also sampled at or beyond the current record. Record its record
 * index extent in uta and utb.
 */
    int found=0;
    int uta=0;
    int utb=0;
    for(spec=spectra->head; spec; spec=spec->next) {
      if(spec->utb >= irec) {
	if(!found || spec->uta < uta) {
	  found = 1;
	  uta = spec->uta;
	  utb = spec->utb;
	};
      };
    };
/*
 * Stop when no further spectra exist.
 */
    if(!found)
      break;
/*
 * Process the located spectrum over the new index range, along with
 * any other spectra that happen to be sampled over the same interval.
 */
    for( ; irec<=utb; irec++) {
/*
 * Read the next integration of raw visibilities from the uvdata.scr
 * paging file.
 */
      if(dp_read(ob->dp, irec))
	return 1;
/*
 * Calibrate the raw data.
 */
      if(dp_cal(ob))
	return 1;
/*
 * Apply the current stream position shift, if any.
 */
      if(dp_shift(ob))
	return 1;
/*
 * Add to the weighted sum spectra from this integration.
 */
      for(spec=spectra->head; spec; spec=spec->next) {
	if(dp_sumspec(ob, spec))
	  return 1;
      };
    };
  } while(irec < ob->nrec);
/*
 * Process all listed spectra together, one integration at a time.
 */
/*
 * Turn the weighted sum spectra into weighted mean spectra.
 */
  for(spec=spectra->head; spec; spec=spec->next) {
/*
 * Each IF has its own spectrum.
 */
    int cif;
    for(cif=0; cif<spec->nif; cif++) {
      IFspec *ifs = spec->ifs + cif;
      Cvis *cvis;
/*
 * Get the weighted mean of each channel in the current IF.
 * Also apply the current weight scale factor to the final weights.
 */
      for(cvis=ifs->chan; cvis < ifs->chan + ifs->nchan; cvis++) {
	if(cvis->wt != 0.0f) {    /* Ignore deleted visibilities */
	  cvis->re /= cvis->wt;
	  cvis->im /= cvis->wt;
	  cvis->wt *= ob->geom.wtscale;
	};
      };
    };
  };
  return 0;
}

/*.......................................................................
 * A private function of get_Spectra() used to Extract the spectrum of
 * given baselines from the current sub-array and integration in the
 * ob->dp I/O buffer, and add it to the spectrum in 'spec'. 
 * Note that 'spec' must have been created for 'ob' with add_Spectrum().
 *
 * Input:
 *  ob   Observation *   The observation descriptor.
 *  spec    Spectrum *   The spectrum container to add to.
 * Output:
 *  return       int     0 - OK.
 *                       1 - Error.
 */
static int dp_sumspec(Observation *ob, Spectrum *spec)
{
  Integration *integ;/* The descriptor of the integration in the I/O buffer */
  Specsub *ssub;     /* Baseline set for the current sub-array */
  Dpage *dp;         /* The descriptor of the uvdata.scr file */
  Dif *dif;          /* Pointer to the visibilities of the IF being processed */
  IFspec *ifs;       /* The spectrum of a single IF */
  If *ifp;           /* The observation descriptor of a single IF */
  int cif;           /* The index of the IF being processed */
  int fc;            /* The index of the channel being processed */
  int bmax;          /* The maximum allowable baseline index */
/*
 * Sanity check the arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "dp_sumspec"))
    return 1;
/*
 * Get the uvdata file descriptor.
 */
  dp = ob->dp;
  if(dp->ut < 0 || dp->ut >= ob->nrec) {
    lprintf(stderr, "dp_sumspec: Invalid integration.\n");
    return 1;
  };
/*
 * Get the associated integration descriptor.
 */
  integ = ob->rec[dp->ut].integ;
/*
 * Is the current integration within the range of the given
 * spectrum?
 */
  if(dp->ut < spec->uta || dp->ut > spec->utb)
    return 0;
/*
 * Get the list of baselines to be sampled, corresponding to the
 * sub-array that is currently in the I/O buffer.
 */
  ssub = spec->ssub + (integ->sub - ob->sub);
/*
 * Are no baselines required for this sub-array?
 */
  if(ssub->nbase <= 0 || ssub->baselines==NULL)
    return 0;
/*
 * Check that the spectrum container is compatible with this observation.
 */
  if(spec->nif != ob->nif) {
    lprintf(stderr, "dp_sumspec: Incompatible IF structure.\n");
    return 1;
  };
/*
 * Get the index of the last baseline that appears both in the current
 * sub-array and the dp I/O buffer.
 */
  bmax = integ->sub->nbase - 1;
  if(bmax > dp->bb)
    bmax = dp->bb;
/*
 * Are there no baselines to be processed?
 */
  if(bmax < dp->ba)
    return 0;
/*
 * Loop through each IF in the I/O buffer.
 */
  for(cif=dp->ia, dif=dp->ifs+cif, ifs=spec->ifs+cif, ifp=ob->ifs+cif;
      cif<=dp->ib;
      cif++, dif++, ifs++, ifp++) {
    float f;           /* Frequency of channel 0 */
    float df;          /* Frequency increment per channel */
    Cvis *svis;        /* Spectrum visibility */
    Dchan *dchan;      /* Channel pointer in the I/O buffer */
/*
 * Check that the IF holds the expected number of channels.
 */
    if(ifs->nchan != ob->nchan) {
      lprintf(stderr, "dp_sumspec: Incompatible channel structure.\n");
      return 1;
    };
/*
 * Determine the frequency of channel 0 and the frequency increment per
 * channel.
 */
    f = ifp->freq;
    df = ifp->df;
/*
 * Form sums for each spectral line channel in the IF.
 */
    for(fc=dp->ca, dchan=dif->chan+fc, svis=ifs->chan+fc;
	fc<=dp->cb; fc++, dchan++, svis++) {
      Dbase *dbase = dchan->base;
/*
 * Iterate over all requested baselines.
 */
      int *bptr = ssub->baselines;
      for( ; bptr < ssub->baselines+ssub->nbase; bptr++) {
	int base = *bptr;
/*
 * Only process baselines that are within the range last read into the
 * I/O buffer.
 */
	if(base >= dp->ba && base <= bmax) {
/*
 * Check whether the UV radius of the baseline is in the required range.
 */
	  Visibility *vis = integ->vis + base;
	  float freq = f + fc * df;
	  float u = vis->u * freq;
	  float v = vis->v * freq;
	  float uvrad = sqrt(u*u+v*v);
/*
 * Is the visibility within the specified range of UV radii?
 */
	  if(spec->uvmax <= 0.0f ||
	     (uvrad >= spec->uvmin && uvrad <= spec->uvmax)) {
/*
 * Get the appropriate polarization.
 */
	    Cvis cvis;
	    spec->obpol.getpol(&spec->obpol, dbase[*bptr].pol, &cvis);
/*
 * Accumulate the visibility spectrum as a weighted sum of good visibilities.
 */
	    if(cvis.wt > 0.0f) {
	      if(spec->dovector) {               /* Vector average */
		svis->re += cvis.wt * cvis.re;
		svis->im += cvis.wt * cvis.im;
	      } else {                           /* Scalar average */
		svis->re += cvis.wt * sqrt(cvis.re*cvis.re + cvis.im*cvis.im);
		svis->im += cvis.wt * ((cvis.re!=0.0f || cvis.im!=0.0f) ?
		  atan2(cvis.im, cvis.re) : 0.0f);
	      };
	      svis->wt += cvis.wt;
	    };
	  };
	};
      };
    };
  };
  return 0;
}


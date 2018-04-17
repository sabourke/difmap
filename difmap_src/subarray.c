#include <stdio.h>
#include <stdlib.h>

#include "logio.h"
#include "obs.h"
#include "scans.h"

/*
 * This file contains Subarray memory management functions.
 * All new_*() functions here-in are capable of re-sizing the arrays
 * of descriptors that they are responsible for if the number of
 * stations,integrations, baselines and/or IFs are changed. To this
 * end they assume that the given sub-array descriptor contains counts
 * of the existing number of stations, baselines, IFs and
 * integrations, (all zero if not previously intialized) and that the
 * newly required sizes are provided as function arguments.
 * New_Subarray() ensures that this is true when allocating new elements.
 *
 * Also, all functions initialize structure members (initially pointers are
 * assigned NULL) to ensure that at any failure point del_Subarray() and
 * its minions may be called safely.
 */

static Subarray *clr_Subarray(Subarray *sub);
static Station *new_Stations(Subarray *sub, int nstat);
static Station *del_Stations(Subarray *sub);
static Baseline *new_Baselines(Subarray *sub, int nif, int nbase);
static Baseline *del_Baselines(Subarray *sub);
static Basmem *new_Basmem(Subarray *sub, int nif, int nbase);
static Basmem *del_Basmem(Subarray *sub);
static Bascor *new_Bascor(Basmem *bm, int nif, int nbase);
static Bascor *del_Bascor(Basmem *bm);
static double *new_p_diff(Subarray *sub, int nif);
static double *del_p_diff(Subarray *sub);

static Integration *new_Integrations(Subarray *sub, int nif, int nbase,
				     int nstat, int ntime);
static Integration *del_Integrations(Subarray *sub);
static Intmem *new_Intmem(Subarray *sub, int nif, int nbase, int nstat,
			  int ntime);
static Intmem *del_Intmem(Subarray *sub);
static Visibility *new_Visibilities(Intmem *imem, int nbase, int ntime);
static Visibility *del_Visibilities(Intmem *imem);
static Intcor *new_Intcor(Intmem *imem, int nif, int ntime);
static Intcor *del_Intcor(Intmem *imem);
static Telcor *new_Telcor(Intmem *imem, int nif, int nstat, int ntime);
static Telcor *del_Telcor(Intmem *imem);
static Baswt *new_Baswt(Basmem *bm, int nif, int nbase);
static Baswt *del_Baswt(Basmem *bm);

/*.......................................................................
 * Allocate and/or (re-)intialize an array of Subarray descriptors.
 * This function only allocates the array of nsub Subarrays. For each
 * element of this array ini_Subarray() must be called to allocate
 * its contents.
 *
 * Input:
 *  ob  Observation * The Observation containing the sub-arrays.
 *  nsub        int   The required number of sub-array descriptors.
 * Output:
 *  return Subarray * The descriptor of the subarray, or NULL onm error.
 */
Subarray *new_Subarray(Observation *ob, int nsub)
{
  Subarray *sub;  /* Pointer into the revised array of sub-array descriptors */
  int isub;       /* The index of a sub-array descriptor in ob->sub[] */
/*
 * Valid Observation descriptor?
 */
  if(ob==NULL) {
    lprintf(stderr, "new_Subarray: NULL Observation descriptor intercepted.\n");
    return NULL;
  };
/*
 * Sanity check the rest of the arguments.
 */
  if(nsub<=0) {
    lprintf(stderr, "new_Subarray: Item count 0 or -ve.\n");
    return NULL;
  };
/*
 * Do we need to allocate an array of Subarrays?
 */
  if(ob->sub==NULL) {
    ob->nsub = 0;
    ob->sub = malloc(sizeof(Subarray) * nsub);
    if(ob->sub==NULL) {
      lprintf(stderr, "new_Subarray: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size the array of sub-array descriptors?
 */
  else if(nsub!=ob->nsub) {
/*
 * If the array is to be shrunk then delete the contents of members that are
 * no longer required before shrinking.
 */
    if(nsub<ob->nsub) {
      sub = &ob->sub[nsub];
      for(isub=nsub; isub<ob->nsub; isub++,sub++)
	clr_Subarray(sub);
    };
/*
 * Re-size the array.
 */
    sub = realloc(ob->sub, sizeof(Subarray) * nsub);
    if(sub)
      ob->sub = sub;
    else if(nsub>ob->nsub) {
      lprintf(stderr, "new_Subarray: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Initialize any new elements.
 */
  if(nsub > ob->nsub) {
    sub = &ob->sub[ob->nsub];
    for(isub=ob->nsub; isub<nsub; isub++,sub++) {
      sub->scangap = DEFGAP;
      sub->datutc = 0.0;
      sub->nif = 0;
      sub->ntime = 0;
      sub->nstat = 0;
      sub->nbase = 0;
      sub->tel = NULL;
      sub->base = NULL;
      sub->binan = NULL;
      sub->p_refant = -1;
      sub->p_diff = NULL;
      sub->integ = NULL;
      sub->ob = ob;
      sub->imem = NULL;
      sub->bmem = NULL;
    };
  };
/*
 * Record the new number of sub-arrays descriptors.
 */
  ob->nsub = nsub;
/*
 * Return the initialized array.
 */
  return ob->sub;
}

/*.......................................................................
 * Delete the array of Subarray descriptors in a given Observation.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation containing the
 *                     sub-array array to be deleted.
 * Output:
 *  return Subarray *  The deleted array (ie. allways NULL).
 */
Subarray *del_Subarray(Observation *ob)
{
  if(ob && ob->sub) {
    int isub;
    for(isub=0; isub<ob->nsub; isub++)
      clr_Subarray(&ob->sub[isub]);
/*
 * Delete the array of descriptors.
 */
    free(ob->sub);
    ob->sub = NULL;
  };
  return NULL;
}

/*.......................................................................
 * Initialize (allocate/resize) the contents of a Subarray descriptor
 *
 * Input:
 *  sub  Subarray *  The pointer to an existing descriptor to be
 *                   (re)-initialized.
 *  nif       int    The number of IFs in the observation.
 *  nbase     int    The number of baselines in the subarray.
 *  nstat     int    The number of stations in the subarray.
 *  ntime     int    The number of integrations in the using this subarray.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
int ini_Subarray(Subarray *sub, int nif, int nbase, int nstat, int ntime)
{
/*
 * Valid Subarray descriptor?
 */
  if(sub==NULL) {
    lprintf(stderr, "ini_Subarray: NULL Subarray descriptor intercepted.\n");
    return 1;
  };
/*
 * Sanity check the arguments.
 */
  if(nif<=0 || nbase<=0 || nstat<=0 || ntime<=0) {
    lprintf(stderr, "ini_Subarray: Item count 0 or -ve.\n");
    return 1;
  };
/*
 * (Re-)allocate/initialize new Station descriptor array.
 */
  if(new_Stations(sub, nstat)==NULL)
    return 1;
/*
 * (Re-)allocate/initialize new Station descriptor array.
 */
  if(new_Baselines(sub, nif, nbase)==NULL)
    return 1;
/*
 * If a Binan descriptor already exists and needs to be re-sized,
 * do so.
 */
  if(sub->binan && nstat!=sub->nstat &&
     new_Binan(sub, nstat, sub->binan->nopcal, sub->binan->numorb)==NULL)
    return 1;
/*
 * (Re-)allocate/initialize a new R-L phase difference array.
 */
  if(new_p_diff(sub, nif) == NULL)
    return 1;
/*
 * (Re-)allocate/initialize new Integration descriptor array.
 */
  if(new_Integrations(sub, nif, nbase, nstat, ntime)==NULL)
    return 1;
/*
 * Record the new numbers of baselines, stations and integrations.
 */
  sub->nif = nif;
  sub->nbase = nbase;
  sub->nstat = nstat;
  sub->ntime = ntime;
/*
 * Return the initialized container.
 */
  return 0;
}

/*.......................................................................
 * Report an error and return true if the given sub-array pointer is
 * NULL. This is the minimal check that functions that take Subarray
 * pointers are expected to make. The error message is:
 *
 *   lprintf(stderr, "%s: NULL Subarray descriptor intercepted.\n", fn);
 *
 * Input:
 *  sub   Subarray *  The pointer to check.
 *  fn        char *  The name of the calling function.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
int sub_bad(Subarray *sub, const char *fn)
{
  if(!sub) {
    lprintf(stderr, "%s: NULL Subarray descriptor intercepted.\n", fn);
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Clear a Subarray container by deleting its contents. Don't delete the
 * container.
 *
 * Input:
 *  sub    Subarray *   The descriptor to be cleared.
 * Output:
 *  return Subarray *   The cleared descriptor == sub.
 */
static Subarray *clr_Subarray(Subarray *sub)
{
  if(sub) {
/*
 * Delete the Station descriptor array.
 */
    del_Stations(sub);
/*
 * Delete the Baseline descriptor array.
 */
    del_Baselines(sub);
/*
 * Delete the AIPS binary AN-table descriptor.
 */
    del_Binan(sub);
/*
 * Delete the array of R-L phase differences.
 */
    del_p_diff(sub);
/*
 * Delete the array of integrations.
 */
    del_Integrations(sub);
/*
 * Zero the index counts.
 */
    sub->nif = 0;
    sub->ntime = 0;
    sub->nstat = 0;
    sub->nbase = 0;
  };
/*
 * Return the empty container.
 */
  return sub;
}

/*.......................................................................
 * (Re-)allocate/initialize an array of Station descriptors for a given
 * Subarray. 
 *
 * Input:
 *  sub   Subarray *  The descriptor of the subarray containing the
 *                    station descriptor array.
 *  nstat      int    The number of stations in the sub-array.
 * Output:
 *  return Station *  Pointer to the array of descriptors, or NULL on
 *                    error.
 */
static Station *new_Stations(Subarray *sub, int nstat)
{
  Station *tel;  /* Pointer into the new station descriptor array */
  int itel;      /* Index of telescope */
/*
 * Allocate a new array?
 */
  if(sub->tel==NULL) {
    sub->tel = malloc(sizeof(Station) * nstat);
    if(sub->tel==NULL) {
      lprintf(stderr, "new_Stations: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing array?
 */
  else if(sub->nstat != nstat) {
    tel = realloc(sub->tel, sizeof(Station) * nstat);
    if(tel)
      sub->tel = tel;
    else if(nstat > sub->nstat) {
      lprintf(stderr, "new_Stations: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Default initialize any new members of the array.
 */
  if(nstat > sub->nstat) {
    tel = &sub->tel[sub->nstat];
    for(itel=sub->nstat; itel<nstat; itel++,tel++) {
      tel->name[0] = '\0';
      tel->antno = 0;
      tel->antfix = 0;
      tel->antwt = 1.0f;
      tel->type = GROUND;
      tel->geo.gnd.x = tel->geo.gnd.y = tel->geo.gnd.z = 0.0f;
      tel->vb = NULL;
    };
  };
/*
 * Return the revised array.
 */
  return sub->tel;
}

/*.......................................................................
 * Delete the array of station descriptors in a given subarray.
 *
 * Input:
 *  sub   Subarray *  The descriptor of the subarray containing the
 *                    station array to be deleted.
 * Output:
 *  return Station *  The deleted array - ie. allways NULL.
 */
static Station *del_Stations(Subarray *sub)
{
  int i;
  if(sub && sub->tel) {
    for(i=0; i<sub->nstat; i++) {
      Station *tel = sub->tel + i;
      tel->vb = del_VoltageBeam(tel->vb);
    };
    free(sub->tel);
    sub->tel = NULL;
  };
  return NULL;
}

/*.......................................................................
 * (Re-)allocate/initialize an array of Baseline descriptors for a given
 * Subarray. 
 *
 * Input:
 *  sub    Subarray *  The descriptor of the subarray containing the
 *                     baseline descriptor array.
 *  nif         int    The new number of IFs for which baseline gain
 *                     corrections are required.
 *  nbase       int    The new number of baselines in the sub-array.
 * Output:
 *  return Baseline *  Pointer to the array of descriptors, or NULL on
 *                     error.
 */
static Baseline *new_Baselines(Subarray *sub, int nif, int nbase)
{
  Baseline *base; /* Pointer into the new baseline descriptor array */
  int ibase;      /* Index of baseline */
/*
 * Allocate a new array?
 */
  if(sub->base==NULL) {
    sub->base = malloc(sizeof(Baseline) * nbase);
    if(sub->base==NULL) {
      lprintf(stderr, "new_Baselines: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing array?
 */
  else if(sub->nbase != nbase) {
    base = realloc(sub->base, sizeof(Baseline) * nbase);
    if(base)
      sub->base = base;
    else if(nbase > sub->nbase) {
      lprintf(stderr, "new_Baselines: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Default initialize any new members of the array.
 */
  if(nbase > sub->nbase) {
    base = &sub->base[sub->nbase];
    for(ibase=sub->nbase; ibase<nbase; ibase++,base++) {
      base->tel_a = base->tel_b = 0;
      base->boff = base->bxy = base->bz = 0.0;
    };
  };
/*
 * Allocate arrays that are to be distributed between each baseline
 * descriptor. Also thread them.
 */
  if(new_Basmem(sub, nif, nbase)==NULL)
    return NULL;
/*
 * Return the revised array.
 */
  return sub->base;
}

/*.......................................................................
 * Delete the array of baseline descriptors in a given subarray.
 *
 * Input:
 *  sub    Subarray * The descriptor of the subarray containing the
 *                    baseline array to be deleted.
 * Output:
 *  return Baseline * The deleted array - ie. allways NULL.
 */
static Baseline *del_Baselines(Subarray *sub)
{
  if(sub && sub->base) {
    free(sub->base);
    sub->base = NULL;
  };
/*
 * Delete the threaded arrays used in the baseline descriptors.
 */
  del_Basmem(sub);
  return NULL;
}

/*.......................................................................
 * Revise the container of arrays to be distributed between baselines.
 * The distributing of the arrays is also performed via this function,
 * so the baseline array must have been revised BEFORE calling this function.
 *
 * Input:
 *  sub    Subarray *  The descriptor of sub-array containing the
 *                     Basmem and Baseline members being updated.
 *  nif         int    The new number of IFs to cater for.
 *  nbase       int    The new number of baselines to cater for.
 * Output:
 *  return   Basmem *  The revised container, or NULL on error.
 */
static Basmem *new_Basmem(Subarray *sub, int nif, int nbase)
{
  Basmem *bm;  /* Pointer to sub->bmem */
/*
 * Create a new container?
 */
  if(sub->bmem==NULL) {
    sub->bmem = malloc(sizeof(Basmem));
    if(sub->bmem == NULL) {
      lprintf(stderr, "new_Basmem: Insufficient memory.\n");
      return NULL;
    };
/*
 * Initialize the container.
 */
    bm = sub->bmem;
    bm->nbase = bm->nif = 0;
    bm->bcor = NULL;
    bm->bwt = NULL;
  };
/*
 * Get a local pointer to the container.
 */
  bm = sub->bmem;
/*
 * Revise the array of baseline corrections.
 */
  if(new_Bascor(bm, nif, nbase)==NULL)
    return NULL;
/*
 * Revise the array of baseline weights.
 */
  if(new_Baswt(bm, nif, nbase)==NULL)
    return NULL;
/*
 * Record the new array sizes.
 */
  bm->nif = nif;
  bm->nbase = nbase;
/*
 * Thread the members into the sub-array baseline array.
 */
  {
    int b;     /* Index of an baseline */
    Baseline *base;
    Bascor *bcor = bm->bcor;
    Baswt *bwt = bm->bwt;
/*
 * Assign nif baseline corrections per baseline.
 */
    for(b=0,base=sub->base; b<nbase; b++,base++,bcor+=nif)
      base->bcor = bcor;
/*
 * Assign nif baseline weights per baseline.
 */
    for(b=0,base=sub->base; b<nbase; b++,base++,bwt+=nif)
      base->bwt = bwt;
  };
/*
 * Return the revised container.
 */
  return sub->bmem;
}

/*.......................................................................
 * Delete the Basmem container of a sub-array, and its contents.
 *
 * Input:
 *  sub   Subarray *  The sub-array containing the container.
 * Output:
 *  return  Basmem *  The deleted container, ie. allways NULL.
 */
static Basmem *del_Basmem(Subarray *sub)
{
  if(sub && sub->bmem) {
    Basmem *bm = sub->bmem;
/*
 * Delete the contents of the container.
 */
    del_Bascor(bm);
    del_Baswt(bm);
/*
 * Delete the container.
 */
    free(sub->bmem);
    sub->bmem = NULL;
  };
  return NULL;
}

/*.......................................................................
 * (Re-)allocate/initialize new a resoff time-invariant baseline correction
 * array in a given IF (and subarray).
 *
 * Input:
 *  bm       Basmem *  The container of the array.
 *  nif         int    The new number of IFs to cater for.
 *  nbase       int    The new number of baselines to cater for.
 * Output:
 *  return   Bascor *  The pointer to the correction array, or NULL on error.
 */
static Bascor *new_Bascor(Basmem *bm, int nif, int nbase)
{
  Bascor *bcor;   /* Pointer into the new correction array */
  int base;       /* Index of a baseline */
  size_t nold;    /* The existing number of elements in the array */
  size_t nnew;    /* The new required number of elements in the array */
/*
 * Determine the exiting and required sizes of the array.
 */
  nold = bm->nif * bm->nbase;
  nnew = nif * nbase;
/*
 * Allocate a new Bascor array?
 */
  if(bm->bcor==NULL) {
    bm->bcor = malloc(sizeof(Bascor) * nnew);
    if(bm->bcor==NULL) {
      lprintf(stderr, "new_Bascor: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing Bascor array?
 */
  else if(nnew != nold) {
    bcor = realloc(bm->bcor, sizeof(Bascor) * nnew);
    if(bcor)
      bm->bcor = bcor;
    else if(nnew > nold) {
      lprintf(stderr, "new_Bascor: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Initialize the new elements of the array.
 */
  if(nnew > nold) {
    bcor = &bm->bcor[nold];
    for(base=nold; base<nnew; base++,bcor++) {
      bcor->amp_cor = 1.0f;
      bcor->phs_cor = 0.0f;
    };
  };
/*
 * Return the new array.
 */
  return bm->bcor;
}

/*.......................................................................
 * Delete the baseline correction array of a given Basmem container.
 *
 * Input:
 *  bm   Basmem *  The descriptor whose Bascor array is to be deleted.
 * Output:
 *  return Bascor *  The deleted array - ie. allways NULL.
 */
static Bascor *del_Bascor(Basmem *bm)
{
  if(bm && bm->bcor) {
    free(bm->bcor);
    bm->bcor = NULL;
  };
  return NULL;
}

/*.......................................................................
 * (Re-)allocate/initialize a time-invariant baseline weight array.
 *
 * Input:
 *  bm       Basmem *  The container of the array.
 *  nif         int    The new number of IFs to cater for.
 *  nbase       int    The new number of baselines to cater for.
 * Output:
 *  return   Baswt *  The pointer to the correction array, or NULL on error.
 */
static Baswt *new_Baswt(Basmem *bm, int nif, int nbase)
{
  Baswt *bwt;     /* Pointer into the new weight array */
  int base;       /* Index of a baseline */
  size_t nold;    /* The existing number of elements in the array */
  size_t nnew;    /* The new required number of elements in the array */
/*
 * Determine the exiting and required sizes of the array.
 */
  nold = bm->nif * bm->nbase;
  nnew = nif * nbase;
/*
 * Allocate a new Baswt array?
 */
  if(bm->bwt==NULL) {
    bm->bwt = malloc(sizeof(Baswt) * nnew);
    if(bm->bwt==NULL) {
      lprintf(stderr, "new_Baswt: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing Baswt array?
 */
  else if(nnew != nold) {
    bwt = realloc(bm->bwt, sizeof(Baswt) * nnew);
    if(bwt)
      bm->bwt = bwt;
    else if(nnew > nold) {
      lprintf(stderr, "new_Baswt: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Initialize the new elements of the array.
 */
  if(nnew > nold) {
    bwt = &bm->bwt[nold];
    for(base=nold; base<nnew; base++,bwt++)
      bwt->wtsum = 0.0f;
  };
/*
 * Return the new array.
 */
  return bm->bwt;
}

/*.......................................................................
 * Delete the baseline weight array of a given Basmem container.
 *
 * Input:
 *  bm    Basmem *  The descriptor whose Baswt array is to be deleted.
 * Output:
 *  return Baswt *  The deleted array - ie. allways NULL.
 */
static Baswt *del_Baswt(Basmem *bm)
{
  if(bm && bm->bwt) {
    free(bm->bwt);
    bm->bwt = NULL;
  };
  return NULL;
}

/*.......................................................................
 * (Re-)allocate/initialize the array of integrations associated with
 * a given sub-array.
 *
 * Input:
 *  sub       Subarray * The descriptor of the sub-array containing the
 *                       integration array.
 *  nif            int   The new number of IFs for which corrections are
 *                       required.
 *  nbase          int   The new number of baselines in the sub-array.
 *  nstat          int   The new number of stations for which corrections
 *                       are required.
 *  ntime          int   The new number of integrations required.
 * Output:
 *  return Integration * The revised array of integrations, or NULL
 *                       on error.
 */
static Integration *new_Integrations(Subarray *sub, int nif, int nbase,
				     int nstat, int ntime)
{
  Integration *integ;  /* Pointer into the array of integrations */
  int ut;              /* The index of an integration */
/*
 * Allocate a new array?
 */
  if(sub->integ==NULL) {
    sub->integ = malloc(sizeof(Integration) * ntime);
    if(sub->integ==NULL) {
      lprintf(stderr, "new_Integrations: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing array?
 */
  else if(sub->ntime != ntime) {
    integ = realloc(sub->integ, sizeof(Integration) * ntime);
    if(integ)
      sub->integ = integ;
    else if(ntime > sub->ntime) {
      lprintf(stderr, "new_Integrations: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Default initialize the new members of the array.
 */
  if(ntime > sub->ntime) {
    integ = &sub->integ[sub->ntime];
    for(ut=sub->ntime; ut<ntime; ut++,integ++) {
      integ->ut = 0.0;
      integ->irec = 0;
      integ->sub = sub;
      integ->vis = NULL;
      integ->icor = NULL;
      integ->edlist = NULL;
    };
  };
/*
 * Revise the arrays that are distributed between integrations.
 */
  if(new_Intmem(sub, nif, nbase, nstat, ntime)==NULL)
    return NULL;
/*
 * Return the revised integration array.
 */
  return sub->integ;
}

/*.......................................................................
 * Delete the array of integrations of a given sub-array.
 *
 * Input:
 *  sub       Subarray * The subarray containing the array of integrations
 *                       to be deleted.
 * Output:
 *  return Integration * The deleted array - ie. allways NULL.
 */
static Integration *del_Integrations(Subarray *sub)
{
  if(sub && sub->integ) {
/*
 * Delete the integration descriptor.
 */
    free(sub->integ);
    sub->integ = NULL;
  };
/*
 * Delete the distributed integration member arrays in sub->imem.
 */
  del_Intmem(sub);
  return NULL;
}

/*.......................................................................
 * Create a new or revised container of arrays to be distributed between
 * integrations. The distributing of the arrays is also performed
 * via this function, so the integration array must have been revised
 * BEFORE calling this function.
 *
 * Input:
 *  sub    Subarray *  The descriptor of sub-array containing the
 *                     Intmem and Integration members being updated.
 *  nif         int    The new number of IFs to cater for.
 *  nbase       int    The new number of baselines to cater for.
 *  nstat       int    The new number of stations to cater for.
 *  ntime       int    The new number of integrations to cater for.
 * Output:
 *  return   Intmem *  The revised container, or NULL on error.
 */
static Intmem *new_Intmem(Subarray *sub, int nif, int nbase, int nstat,
			  int ntime)
{
  Intmem *im;  /* Pointer to sub->imem */
/*
 * Create a new container?
 */
  if(sub->imem==NULL) {
    sub->imem = malloc(sizeof(Intmem));
    if(sub->imem == NULL) {
      lprintf(stderr, "new_Intmem: Insufficient memory.\n");
      return NULL;
    };
/*
 * Initialize the container.
 */
    im = sub->imem;
    im->ntime = im->nbase = im->nstat = im->nif = 0;
    im->vis = NULL;
    im->icor = NULL;
    im->tcor = NULL;
  };
/*
 * Get a local pointer to sub->imem.
 */
  im = sub->imem;
/*
 * Revise the array of visibilities.
 */
  if(new_Visibilities(im, nbase, ntime)==NULL)
    return NULL;
/*
 * Revise the array of Integration/IF corrections.
 */
  if(new_Intcor(im, nif, ntime)==NULL)
    return NULL;
/*
 * Revise the array of telescope gain corrections.
 */
  if(new_Telcor(im, nif, nstat, ntime)==NULL)
    return NULL;
/*
 * Thread the members into the sub-array integration array.
 */
  {
    int ut;    /* Index of an integration */
    int cif;   /* Index of an IF */
    Integration *integ = sub->integ;
    Visibility *vis = im->vis;
    Intcor *icor = im->icor;
    Telcor *tcor = im->tcor;
/*
 * Thread nbase visibilities and nif integration corrections per integration.
 */
    for(ut=0; ut<ntime; ut++,integ++,vis+=nbase) {
      integ->vis = vis;
      integ->icor = icor;
/*
 * Thread nstat telescope corrections into each icor element.
 */
      for(cif=0; cif<nif; cif++,icor++,tcor+=nstat)
	icor->tcor = tcor;
    };
  };
/*
 * Record the revised sizes of the contained arrays.
 */
  im->ntime = ntime;
  im->nbase = nbase;
  im->nstat = nstat;
  im->nif   = nif;
/*
 * Return the revised container.
 */
  return sub->imem;
}

/*.......................................................................
 * Delete the Intmem container of a sub-array, and its contents.
 *
 * Input:
 *  sub   Subarray *  The sub-array containing the container.
 * Output:
 *  return  Intmem *  The deleted container, ie. allways NULL.
 */
static Intmem * del_Intmem(Subarray *sub)
{
  if(sub && sub->imem) {
    Intmem *im = sub->imem;
/*
 * Delete the contents of the container.
 */
    del_Visibilities(im);
    del_Intcor(im);
    del_Telcor(im);
/*
 * Delete the container.
 */
    free(sub->imem);
    sub->imem = NULL;
  };
  return NULL;
}

/*.......................................................................
 * (Re-)allocate/initialize an array of ntime * nbase visibilities.
 *
 * Input:
 *  im         Intmem * The container of the visibility array.
 *  nbase         int   The new number of baselines to cater for.
 *  ntime         int   The new number of integrations to cater for.
 * Output:
 *  return Visibility * The revised array, or NULL on error.
 */
static Visibility *new_Visibilities(Intmem *im, int nbase, int ntime)
{
  Visibility *vis; /* Pointer into the visibility array */
  size_t nold;     /* The existing number of elements in the array */
  size_t nnew;     /* The new required number of elements in the array */
  size_t i;
/*
 * Determine the existing and new number of elements in the array.
 */
  nold = im->ntime * im->nbase;
  nnew = ntime * nbase;
/*
 * Allocate a new array?
 */
  if(im->vis==NULL) {
    im->vis = malloc(sizeof(Visibility) * nnew);
    if(im->vis==NULL) {
      lprintf(stderr, "new_Visibilities: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing array?
 */
  else if(nnew != nold) {
/*
 * Re-size the array.
 */
    vis = realloc(im->vis, sizeof(Visibility) * nnew);
    if(vis)
      im->vis = vis;
    else if(nnew > nold) {
      lprintf(stderr, "new_Visibilities: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Default initialize any new elements of the array.
 */
  if(nnew > nold) {
    vis = &im->vis[nold];
    for(i=nold; i<nnew; i++,vis++) {
      vis->amp = vis->modamp = 0.0f;
      vis->phs = vis->modphs = 0.0f;
      vis->wt = 0.0f;
      vis->u = vis->v = vis->w = 0.0f;
      vis->dt = 0.0f;
      vis->bad = FLAG_DEL;
    };
  };
/*
 * Return the revised array.
 */
  return im->vis;
}

/*.......................................................................
 * Delete the array of visibilities of a sub-array.
 *
 * Input:
 *  im         Intmem *  The container of the array to be deleted.
 * Output:
 *  return Visibility *  The deleted array - ie. allways NULL.
 */
static Visibility *del_Visibilities(Intmem *im)
{
  if(im && im->vis) {
    free(im->vis);
    im->vis = NULL;
  };
  return NULL;
}

/*.......................................................................
 * Revise an array of nif * ntime containers of arrays of nstat telescope
 * corrections. 
 *
 * Input:
 *  im      Intmem *  The container of the Intcor array.
 *  nif        int    The new number of IFs.
 *  ntime      int    The new number of integrations.
 * Output:
 *  return  Intcor *  The array of containers, or NULL on error.
 */
static Intcor *new_Intcor(Intmem *im, int nif, int ntime)
{
  Intcor *icor;    /* Pointer into the Intcor array */
  size_t nold;     /* The existing number of elements in the array */
  size_t nnew;     /* The new required number of elements in the array */
  size_t i;
/*
 * Determine the existing and new number of elements in the array.
 */
  nold = im->ntime * im->nif;
  nnew = ntime * nif;
/*
 * Allocate a new array?
 */
  if(im->icor==NULL) {
    im->icor = malloc(sizeof(Intcor) * nnew);
    if(im->icor==NULL) {
      lprintf(stderr, "new_Intcor: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing array?
 */
  else if(nnew != nold) {
/*
 * Re-size the array.
 */
    icor = realloc(im->icor, sizeof(Intcor) * nnew);
    if(icor)
      im->icor = icor;
    else if(nnew > nold) {
      lprintf(stderr, "new_Intcor: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Default initialize any new elements of the array.
 */
  if(nnew > nold) {
    icor = &im->icor[nold];
    for(i=nold; i<nnew; i++,icor++)
      icor->tcor = NULL;
  };
/*
 * Return the revised array.
 */
  return im->icor;
}

/*.......................................................................
 * Delete an array of nif * ntime containers of arrays of nstat telescope
 * corrections. 
 *
 * Input:
 *  im      Intmem *  The container of the Intcor array.
 * Output:
 *  return  Intcor *  The deleted array, ie. allways NULL.
 */
static Intcor *del_Intcor(Intmem *im)
{
  if(im && im->icor) {
    free(im->icor);
    im->icor = NULL;
  };
  return NULL;
}

/*.......................................................................
 * (Re-)allocate/initialize an array of ntime * nif * nstat telescope
 * corrections.
 *
 * Input:
 *  im     Intmem *  The container of the array.
 *  nif       int    The new number of IFs to cater for.
 *  nstat     int    The new number of stations to cater for.
 *  ntime     int    The new number of integrations to cater for.
 * Output:
 *  return Telcor *  The revised array, or NULL on error.
 */
static Telcor *new_Telcor(Intmem *im, int nif, int nstat, int ntime)
{
  Telcor *tcor; /* Pointer into the telescope correction array */
  size_t nold;  /* The existing number of elements in the array */
  size_t nnew;  /* The new required number of elements in the array */
  size_t i;
/*
 * Determine the existing and new number of elements in the array.
 */
  nold = (size_t) im->ntime * im->nif * im->nstat;
  nnew = (size_t) ntime * nif * nstat;
/*
 * Allocate a new array?
 */
  if(im->tcor==NULL) {
    im->tcor = malloc(sizeof(Telcor) * nnew);
    if(im->tcor==NULL) {
      lprintf(stderr, "new_Telcor: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing array?
 */
  else if(nnew != nold) {
/*
 * Re-size the array.
 */
    tcor = realloc(im->tcor, sizeof(Telcor) * nnew);
    if(tcor)
      im->tcor = tcor;
    else if(nnew > nold) {
      lprintf(stderr, "new_Telcor: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Default initialize any new elements of the array.
 */
  if(nnew > nold) {
    tcor = &im->tcor[nold];
    for(i=nold; i<nnew; i++,tcor++) {
      tcor->amp_cor = 1.0f;
      tcor->phs_cor = 0.0f;
      tcor->bad = 0;
    };
  };
/*
 * Return the revised array.
 */
  return im->tcor;
}

/*.......................................................................
 * Delete an array of ntime * nif * nstat telescope corrections.
 *
 * Input:
 *  im     Intmem *  The container of the array.
 * Output:
 *  return Telcor *  The deleted array - ie. allways NULL.
 */
static Telcor *del_Telcor(Intmem *im)
{
  if(im && im->tcor) {
    free(im->tcor);
    im->tcor = NULL;
  };
  return NULL;
}

/*.......................................................................
 * (Re-)allocate/initialize an array of R-L phase differences.
 *
 * Input:
 *  sub   Subarray *  The descriptor of the subarray containing the
 *                    station descriptor array.
 *  nif        int    The number of IFs to allocate elements for.
 * Output:
 *  return  double *  Pointer to the array of 'nif' phase differences,
 *                    or NULL on error.
 */
static double *new_p_diff(Subarray *sub, int nif)
{
  double *p_diff; /* Pointer into the new array */
  int cif;        /* Index of IF */
/*
 * Allocate a new array?
 */
  if(sub->p_diff==NULL) {
    sub->p_diff = malloc(sizeof(double) * nif);
    if(sub->p_diff==NULL) {
      lprintf(stderr, "new_p_diff: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing array?
 */
  else if(sub->nif != nif) {
    p_diff = realloc(sub->p_diff, sizeof(double) * nif);
    if(p_diff)
      sub->p_diff = p_diff;
    else if(nif > sub->nif) {
      lprintf(stderr, "new_p_diff: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Default initialize any new members of the array.
 */
  if(nif > sub->nif) {
    p_diff = &sub->p_diff[sub->nif];
    for(cif=sub->nif; cif<nif; cif++,p_diff++)
      *p_diff = 0.0f;
  };
/*
 * Return the revised array.
 */
  return sub->p_diff;
}

/*.......................................................................
 * Delete a sub->p_diff array of phase differences.
 *
 * Input:
 *  sub    Subarray *  The descriptor of the sub-array to delete from.
 * Output:
 *  return    double *  Allways NULL.
 */
static double *del_p_diff(Subarray *sub)
{
  if(sub && sub->p_diff)
    free(sub->p_diff);
  return NULL;
}

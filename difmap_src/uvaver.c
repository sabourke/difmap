#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "logio.h"
#include "obs.h"
#include "visaver.h"
#include "utbin.h"
#include "vlbconst.h"
#include "obedit.h"

/*
 * Class to record the details of the next un-averaged solution bin in a
 * given sub-array.
 */

typedef struct Solbin {
  double ut;            /* Time stamp of the pending averaged integration */
  Subarray *sub;        /* The descriptor of the associated sub-array */
  Integration *integ;   /* Pointer to first un-used integration in sub-array */
  Integration *aver;    /* Pointer to output averaged integration <= integ */
  int ibin;             /* The index of the current bin (0 relative) */
  int ntime;            /* Number of sub-array integrations in solution bin */
  int nleft;            /* The number of un-processed integrations remaining */
  struct Solbin *next;  /* Pointer to the container with next highest UT */
} Solbin;

/*
 * Descriptor of solution bin iterator.
 */

typedef struct {
  double avtime;    /* Width of solution bin (seconds) */
  double origin;    /* The origin time of the solution bins */
  long irec;        /* Record number in output file */
  int nsub;         /* The number of sub-arrays */
  Solbin *binmem;   /* Dynamic array of nsub un-ordered Solbin containers */
  Solbin *head;     /* Head of list of pending solution bins for each */
                    /* Subarray, in order of increasing solution UT */
} Biniter;

static Biniter *new_Biniter(Observation *ob, double avtime);
static Biniter *del_Biniter(Biniter *iter);
static Solbin *nextbin(Biniter *iter);        /* Solution iterator */
static void newbin(Biniter *iter, Solbin *sbin);

typedef struct {
  long nrec;         /* The total number of averaged integrations */
  Biniter *iter;     /* The solution bin iterator */
  Visaver *av;       /* The visibility averaging container */
  Dpage *dp;         /* The output uvdata.scr file */
} UVaver;

static UVaver *new_UVaver(Observation *ob, double avtime, int scatter);
static UVaver *del_UVaver(UVaver *av);

static int dp_aver(Observation *ob, UVaver *av, Solbin *sbin);

static Observation *uvend(Observation *ob, UVaver *av);

/*.......................................................................
 * Perform a coherent average of a UV data set and shrink the data-set
 * to reclaim the newly released memory and scratch file disk space.
 *
 * Input:
 *  ob     Observation *  The UV data set to be averaged.
 *  avtime      double    The averaging time in seconds.
 *  scatter        int    If true, attempt to compute uncertainties from
 *                        the scatter of the input visibilities.
 * Output:
 *  return Observation *  The averaged 'ob' or NULL if 'ob' had to be
 *                        deleted (if all data flagged).
 */
Observation *uvaver(Observation *ob, float avtime, int scatter)
{
  UVaver *av;        /* Container for intermediate objects */
  Solbin *sbin;      /* Solution bin descriptor */
  int isub;          /* The index of the sub-array being processed */
  int was_select;    /* If true a selection existed before averaging */
/*
 * Sanity check the observation descriptor.
 */
  if(!ob_ready(ob, OB_INDEX, "uvaver"))
    return del_Observation(ob);
/*
 * See if their is a current selection, so that we know whether to
 * restore one.
 */
  was_select = ob_ready(ob, OB_SELECT, NULL);
/*
 * Flush cached edits.
 */
  if(ed_flush(ob))
    return ob;
/*
 * Allocate intermediate objects, including the output file.
 */
  av = new_UVaver(ob, avtime, scatter);
  if(av==NULL)
    return ob;    /* No harm done yet - ob is still valid */
/*
 * Update the user.
 */
  lprintf(stdout, "Averaging into %g second bins.\n", avtime);
/*
 * Up to this point no changes have been made to the Observation structure,
 * and error recovery was painless. If an error occurs from here onwards
 * the observation structure must be deleted.
 *
 * Compose the averaged visibilities of each solution bin from the existing
 * ob->dp uvdata.scr file and write them to the output av->dp uvdata.scr
 * file.
 */
  while((sbin=nextbin(av->iter)) != NULL) {
    if(dp_aver(ob, av, sbin))
      return uvend(del_Observation(ob), av);
  };
/*
 * Clear self-cal and resoff correction records. Such corrections have now been
 * frozen into the averaged data.
 */
  uncalib(ob, 1, 1, 1, 1);
  ini_bcor(ob, 1, 1, 1);
/*
 * If the scatter method was selected, such that weights were deduced
 * directly from the data, reset the weight scale factor to 1.
 */
  if(scatter)
    ob->geom.wtscale = 1.0;
/*
 * Delete the un-averaged uvdata.scr paging file and replace it with
 * the new scratch file.
 */
  ob->dp = del_Dpage(ob->dp);
  ob->dp = av->dp;
/*
 * The new uvdata.scr file is now referenced by ob, and must not be
 * redundantly deleted when del_UVaver() is called, so unreference it here.
 */
  av->dp = NULL;
/*
 * Shrink wrap the Observation structure and create the new IF and UV model
 * scratch files.
 */
  ob = Obs_alloc(ob, av->nrec, ob->nbmax, ob->nsub, ob->nif, ob->npol, ob->nchan);
  if(ob==NULL)
    return uvend(del_Observation(ob), av);
/*
 * Shrink-wrap the sub-array descriptors.
 */
  for(isub=0; isub<ob->nsub; isub++) {
    Subarray *sub = &ob->sub[isub];
    if(ini_Subarray(sub, sub->nif, sub->nbase, sub->nstat,
		    av->iter->binmem[isub].ibin))
      return uvend(del_Observation(ob), av);
  };
/*
 * We have data in the new ob->dp paging file.
 */
  ob->state = OB_DATA;
/*
 * Index the new integrations in ob->rec[].
 */
  if(ini_Intrec(ob))
    return uvend(del_Observation(ob), av);
/*
 * Re-instate the previous stream selection if one was selected
 * before uvaver was called.
 */
  if(was_select && ob_select(ob, 1, ob->stream.cl, ob->stream.pol.type))
    return uvend(del_Observation(ob), av);
/*
 * Return the succesfully averaged data.
 */
  return uvend(ob, av);
}

/*.......................................................................
 * Private cleanup function of uvaver().
 */
static Observation *uvend(Observation *ob, UVaver *av)
{
  del_UVaver(av);
  return ob;
}

/*.......................................................................
 * Create a new solution bin iterator.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation being averaged.
 *  avtime    double    The solution bin width (seconds).
 * Output:
 *  return   Biniter *  The solution bin iterator.
 */
static Biniter *new_Biniter(Observation *ob, double avtime)
{
  Biniter *iter;   /* The new iterator descriptor */
/*
 * Sanity checks.
 */
  if(avtime < 1.0) {
    lprintf(stderr, "new_Biniter: Solution bin width less than 1 second.\n");
    return NULL;
  };
/*
 * Allocate the iterator descriptor.
 */
  iter = malloc(sizeof(Biniter));
  if(iter==NULL) {
    lprintf(stderr, "new_Biniter: Insufficient memory for iterator.\n");
    return NULL;
  };
/*
 * Initialize the descriptor to at least the point at which it can be sent
 * to del_Biniter() on error.
 */
  iter->avtime = avtime;
/*
 * The origin of the binning grid is UT=0 on first day of the observation.
 */
  iter->origin = ob->date.ut;
  iter->irec = -1;
  iter->nsub = ob->nsub;
  iter->binmem = NULL;
  iter->head = NULL;
/*
 * Allocate a solution bin descriptor for each sub-array.
 */
  iter->binmem = malloc(sizeof(Solbin) * iter->nsub);
  if(iter->binmem==NULL) {
    lprintf(stderr, "new_Biniter: Insufficient memory for iterator.\n");
    return del_Biniter(iter);
  };
/*
 * Initialize each element of binmem - only initialize members that
 * are not reset at the end of each cycle through the iterator.
 * The first call to nextbin() will do the rest.
 */
  {
    int isub;
    Subarray *sub = ob->sub;
    Solbin *bin = iter->binmem;
    for(isub=0; isub<ob->nsub; isub++,bin++,sub++)
      bin->sub = sub;
  };
/*
 * Initialize the list of solution bins.
 */
  iter->head = NULL;
/*
 * Return the new iterator.
 */
  return iter;
}

/*.......................................................................
 * Delete a solution bin iterator.
 *
 * Input:
 *  iter    Biniter *  The iterator to be deleted.
 * Output:
 *  return  Biniter *  The deleted iterator (allways NULL).
 *                     Use like: iter = del_Biniter(iter);
 */
static Biniter *del_Biniter(Biniter *iter)
{
  if(iter) {
    if(iter->binmem)
      free(iter->binmem);
    free(iter);
  };
  return NULL;
}

/*.......................................................................
 * Get the descriptor of the next solution bin to be processed.
 *
 * If there are no further solution bins, the function returns NULL, and
 * the iterator is re-started on the next call.
 *
 * Input:
 *  iter     Biniter *  An iterator descriptor returned by new_Biniter().
 * Output:
 *  return   Solbin *  The descriptor of the solution bin to be processed
 *                      next, or NULL if there are no more bins.
 */
static Solbin *nextbin(Biniter *iter)
{
  Solbin *sbin;   /* The descriptor of the solution bin being processed */
/*
 * Initialize the iterator if the iterator list is empty. This
 * is the case both on the first call to the iterator, and also when the
 * previous call ended a cycle through the iterator.
 */
  if(iter->head == NULL) {
    int isub;      /* The index of the sub-array being processed */
/*
 * Initialize the solution bin list with the first solution of
 * each sub-array.
 */
    sbin = iter->binmem;
    for(isub=0; isub<iter->nsub; isub++,sbin++) {
      sbin->integ = sbin->aver = sbin->sub->integ;
      sbin->nleft = sbin->sub->ntime;
      sbin->ntime = 0;
      sbin->ibin = 0;
/*
 * Initialize the solution-bin to the first solution of the associated
 * sub-array and insert the result into the iterator list.
 */
      newbin(iter, sbin);
    };
/*
 * Rewind output record number.
 */
    iter->irec = 0;
  }
/*
 * Iterate.
 */
  else {
/*
 * Remove the descriptor of the completed solution bin from the head of the
 * solution bin list, but maintain a pointer to it in sbin.
 */
    sbin = iter->head;
    iter->head = iter->head->next;
/*
 * Update the sub-array bin index and the output-file record number.
 */
    iter->irec++;
    sbin->ibin++;
/*
 * Advance over the integrations of the completed solution bin.
 */
    sbin->nleft -= sbin->ntime;
    sbin->integ += sbin->ntime;
/*
 * Advance the output integration pointer.
 */
    sbin->aver++;
/*
 * Re-initialize the descriptor of the completed solution bin with the
 * next solution in its sub-array if there is one and re-insert it in the
 * iterator list. If there is no further solution, then the solution bin
 * will be left out of the list since it is no longer required.
 */
    newbin(iter, sbin);
  };
/*
 * Return the solution bin at the head of the iterator list.
 */
  return iter->head;
}

/*.......................................................................
 * Get the details of the next solution bin and insert the result in the
 * solution bin list. If there are no further solution bins then
 * the obsolete descriptor will not be re-inserted in the list.
 *
 * Input:
 *  iter   Biniter *  The descriptor of the iterator.
 * Input/Output:
 *  sbin    Solbin *  The descriptor of the solution bin.
 */
static void newbin(Biniter *iter, Solbin *sbin)
{
/*
 * Are there any integrations left to form solution bins from?
 */
  if(sbin->nleft > 0) {
    Integration *integ; /* Integration being checked for inclusion in bin */
    int nleft;          /* Number of integrations left to be considered */
/*
 * Determine the time stamp of start and end of the new solution bin in such
 * a way that integration bins are centered on a regular grid of interval
 * iter->avtime wrt a bin grid origin of iter->origin.
 */
    UTbin *utbin = bintime(iter->origin, sbin->integ->ut, iter->avtime);
/*
 * Find the last integration that lies in the solution bin.
 */
    integ = sbin->integ;
    nleft = sbin->nleft;
    while(nleft>0 && integ->ut<=utbin->end_ut) {
      nleft--;
      integ++;
    };
/*
 * How many integrations are there in the new bin?
 */
    sbin->ntime = integ - sbin->integ;
/*
 * Record the time-stamp of the center of the bin.
 */
    sbin->ut = utbin->mid_ut;
/*
 * Insert the new solution bin at the position in the list that will keep
 * the list in order of time-stamp.
 */
    {
      Solbin *node;       /* Pointer to node being checked */
      Solbin *prev=NULL;  /* Pointer to node that precedes 'node' */
      for(node=iter->head; node!=NULL && sbin->ut > node->ut; node=node->next)
	prev = node;
/*
 * Insert at head of list?
 */
      if(prev==NULL) {
	sbin->next = node;
	iter->head = sbin;
      } else {
	sbin->next = node;
	prev->next = sbin;
      };
    };
  };
  return;
}

/*.......................................................................
 * Take records from the existing ob->dp uvdata.scr file, average
 * self-cal'd versions of them into solution bins and place the results
 * in a new uvdata.scr file av->dp.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation being averaged.
 *  av       UVaver *  The average state descriptor.
 *  sbin     Solbin *  The solution bin to be averaged.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int dp_aver(Observation *ob, UVaver *av, Solbin *sbin)
{
  Integration *integ;   /* Pointer into sbin->integ[0..sbin->ntime-1] */
  int nbase;            /* The number of baselines in the current sub-array */
  int i;
/*
 * Initialize to read whole integrations.
 */
  if(dp_crange(ob->dp, 0, ob->nchan-1) ||
     dp_irange(ob->dp, 0, ob->nif-1)   ||
     dp_brange(ob->dp, 0, ob->nbmax-1) ||
     dp_srange(ob->dp, 0, ob->npol-1))
    return 1;
/*
 * Determine the number of baselines in the sub-array of the latest solution
 * bin.
 */
  nbase = sbin->sub->nbase;
/*
 * Clear the accumulated visibility sums of the last solution in preparation
 * for the next solution.
 */
  if(av_newint(av->av, sbin->aver->vis, nbase, av->iter->irec))
    return 1;
/*
 * Process each integration of the solution bin.
 */
  integ = sbin->integ;
  for(i=0; i<sbin->ntime; i++,integ++) {
/*
 * Read the next un-averaged integration to be included in the bin and
 * apply self-cal and resoff corrections.
 */
    if(dp_read(ob->dp, integ->irec) || dp_cal(ob))
      return 1;
/*
 * Accumulate the weighted mean U,V,W for each baseline.
 */
    {
/*
 * Loop through each IF in the input integration record.
 */
      Dif *ifs = ob->dp->ifs;
      int cif;
      for(cif=0; cif<ob->nif; cif++,ifs++) {
/*
 * Loop through all spectral-line channels in the current IF.
 */
	Dchan *dchan = ifs->chan;
	int fc;
	for(fc=0; fc<ob->nchan; fc++,dchan++) {
/*
 * Loop through each baseline, and for each baseline accumulate the
 * weighted mean U,V,W coordinates and the sum of weights.
 */
	  Visibility *vis = integ->vis;  /* Input visibility */
	  Dbase *dbase = dchan->base;
	  int base;
	  for(base=0; base<nbase; base++,dbase++,vis++) {
/*
 * Loop through each polarization recorded on the current baseline.
 */
	    Cvis *cvis = dbase->pol;
	    long ivis = cvis - ob->dp->cvis; /* Index into integration record */
	    int pol;
	    for(pol=0; pol<ob->npol; pol++,cvis++,ivis++) {
/*
 * Accumulate running mean data values in the output buffer.
 */
	      if(av_dp(av->av, cvis->re, cvis->im, cvis->wt, ivis))
		return 1;
/*
 * Accumulate running mean U,V,W values in the output visibility descriptors.
 */
	      if(av_uvwt(av->av, vis->u, vis->v, vis->w, cvis->wt, vis->dt,
			 base))
		return 1;
	    };
	  };
	};
      };
    };
  };
/*
 * Set U,V and W to zero on baselines which were totally un-sampled, and
 * if requested, replace output weights with those deduced from the data
 * scatter.
 */
  if(av_endint(av->av))
    return 1;
/*
 * Write the averaged solution bin to the output uvdata.scr file.
 */
  if(dp_write(av->dp, av->iter->irec))
    return 1;
/*
 * Assign the output record index and time-stamp to the new integration.
 */
  sbin->aver->irec = av->iter->irec;
  sbin->aver->ut = sbin->ut;
/*
 * Ready for next bin.
 */
  return 0;
}

/*.......................................................................
 * Allocate and return a UVaver class.
 *
 *  ob   Observation *  The descriptor of the observation being averaged.
 *  avtime    double    The solution bin width (seconds).
 *  scatter      int    If true, allocate an extra buffer to use in finding
 *                      uncertainties from the scatter of the visibilities.
 * Output:
 *  return    UVaver *  The new uvaver class instance, or NULL on error.
 */
static UVaver *new_UVaver(Observation *ob, double avtime, int scatter)
{
  UVaver *av = NULL;  /* The pointer to the new uvaver instance */
/*
 * Allocate the container.
 */
  av = malloc(sizeof(UVaver));
  if(av==NULL) {
    lprintf(stderr, "new_UVaver: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize at least up to the point at which del_UVaver() can safely be
 * called.
 */
  av->av = NULL;
  av->iter = NULL;
  av->dp = NULL;
/*
 * Allocate a solution-bin iterator.
 */
  av->iter = new_Biniter(ob, avtime);
  if(av->iter==NULL)
    return del_UVaver(av);
/*
 * Loop through the solution bin iterator to count the total number of
 * output integration records required.
 */
  for(av->nrec=0; nextbin(av->iter) != NULL; av->nrec++);
/*
 * Open the output uvdata.scr file.
 */
  av->dp = new_Dpage(av->nrec, ob->nbmax, ob->nchan, ob->nif, ob->npol);
  if(av->dp==NULL)
    return del_UVaver(av);
/*
 * Initialize to write whole integrations.
 */
  if(dp_crange(av->dp, 0, ob->nchan-1) ||
     dp_irange(av->dp, 0, ob->nif-1)   ||
     dp_brange(av->dp, 0, ob->nbmax-1) ||
     dp_srange(av->dp, 0, ob->npol-1))
    return del_UVaver(av);
/*
 * Allocate the visibility averaging container.
 */
  av->av = new_Visaver(av->dp, avtime, scatter);
  if(av->av==NULL)
    return del_UVaver(av);
/*
 * Return the new uvaver instance.
 */
  return av;
}

/*.......................................................................
 * Delete a uvaver class instance.
 *
 * Input:
 *  av     UVaver *  The uvaver instance to be deleted.
 * Output:
 *  return UVaver *  The deleted instance, ie NULL. Use like:
 *                   av = del_UVaver(av);
 */
static UVaver *del_UVaver(UVaver *av)
{
  if(av) {
/*
 * Delete the contents.
 */
    del_Biniter(av->iter);
    del_Dpage(av->dp);
    del_Visaver(av->av);
/*
 * Delete the container.
 */
    free(av);
  };
  return NULL;
}


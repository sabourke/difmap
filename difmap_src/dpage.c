#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "logio.h"
#include "recio.h"
#include "dpage.h"
#include "obs.h"
#include "vlbconst.h"

/*
 * Define a Cvis structure used to initialize visibilities in the
 * integration buffer.
 */
static Cvis zero_vis = {0.0f,0.0f,0.0f};

static Dpage *dpmemerr(Dpage *dp);
static int dp_error(Dpage *dp, const char *fname);

/*.......................................................................
 * Allocate and intialize a Data-paging file descriptor. This is used
 * for buffered reading and writing of up to an integration at a time
 * from the data paging file.
 *
 * Input:
 *  ntime    int   The number of integrations in the observation.
 *  nbase    int   The number of baselines per integration.
 *  nchan    int   The number of spectral line channels per integration.
 *  nif      int   The number of IFs per integration.
 *  npol     int   The number of stokes parameters per integration.
 * Output:
 *  return Dpage * The allocated and initialized descriptor, or NULL
 *                 on error.
 */
Dpage *new_Dpage(int ntime, int nbase, int nchan, int nif, int npol)
{
  Dpage *dp;    /* Pointer to the new descriptor */
  size_t node;  /* The index of the visibility tree node being intialized */
  size_t nnode; /* The number of nodes to be initialized */
  Dif *dif;     /* Array of IF nodes */
  Dchan *dchan; /* Array of spectral-line channel nodes */
  Dbase *dbase; /* Array of baseline nodes */
/*
 * Attempt to allocate the new descriptor.
 */
  dp = malloc(sizeof(*dp));
  if(dp==NULL)
    return dpmemerr(dp);
/*
 * Intialize dp at least up to the point at which it can safely be passed to
 * del_Dpage().
 */
  dp->rio = NULL;
  dp->cvis = NULL;
  dp->ifs = NULL;
  dp->nvis = 0;
  dp->first = 0;
  dp->nbuff = 0;
  dp->ioerr = 0;
/*
 * Install input parameters.
 */
  dp->ntime = ntime;
  dp->nbase = nbase;
  dp->nchan = nchan;
  dp->nif = nif;
  dp->npol = npol;
/*
 * Determine the indexing offsets between each of the elements of the
 * integration record.
 */
  dp->soff = 1;
  dp->boff = dp->soff * dp->npol;
  dp->coff = dp->boff * dp->nbase;
  dp->ioff = dp->coff * dp->nchan;
/*
 * Initialize the read/write visibility ranges.
 */
  dp->ca = 0;  dp->cb = nchan-1;
  dp->ia = 0;  dp->ib = nif-1;
  dp->sa = 0;  dp->sb = npol-1;
  dp->ba = 0;  dp->bb = nbase-1;
/*
 * Determine the number of visibilities per integration.
 */
  dp->nbuff = dp->nvis = (size_t) nbase * nchan * nif * npol;
/*
 * Open the paging scratch file.
 */
  dp->rio = new_Recio("uvdata.scr", IS_SCR, 0, dp->nvis * sizeof(Cvis));
  if(dp->rio==NULL)
    return del_Dpage(dp);
/*
 * Allocate a visibility buffer of this size.
 */
  dp->cvis = malloc(dp->nvis * sizeof(Cvis));
  if(dp->cvis==NULL)
    return dpmemerr(dp);
/*
 * Initialize the buffer.
 */
  if(dp_clear(dp, -1))
    return del_Dpage(dp);
/*
 * Now set up the hierarchical tree of visibility levels under dp->ifs.
 * First allocate an array of dp->nif nodes conataining pointers to arrays of
 * spetral-line channel nodes.
 */
  dif = dp->ifs = malloc(dp->nif * sizeof(Dif));
  if(dif==NULL)
    return dpmemerr(dp);
/*
 * Allocate an array of dp->nif * dp->nchan spectral-line channel nodes.
 */
  dchan = dif->chan = malloc(dp->nif * dp->nchan * sizeof(Dchan));
  if(dchan == NULL)
    return dpmemerr(dp);
/*
 * Apportion the array of spectral-line channel nodes between IF nodes.
 */
  nnode = dp->nif;
  for(node=1; node<nnode; node++)
    dif[node].chan = dif[node-1].chan + dp->nchan;
/*
 * Allocate an array of dp->nif * dp->nchan * dp->nbase baseline nodes.
 */
  dbase = dchan->base = malloc(sizeof(Dbase) * dp->nif * dp->nchan * dp->nbase);
  if(dbase==NULL)
    return dpmemerr(dp);
/*
 * Apportion the array of baseline nodes between spectral-line channel nodes.
 */
  nnode = dp->nif * dp->nchan;
  for(node=1; node<nnode; node++)
    dchan[node].base = dchan[node-1].base + dp->nbase;
/*
 * The leaves of the tree are the visibilities in the dp->cvis[]
 * integration. Point each baseline node at the corresponding array of
 * dp->npol visibilities in the integration record buffer.
 */
  dbase->pol = dp->cvis;
  nnode = dp->nif * dp->nchan * dp->nbase;
  for(node=1; node<nnode; node++)
    dbase[node].pol = dbase[node-1].pol + dp->npol;
/*
 * Return the intialized descriptor.
 */
  return dp;
}

/*.......................................................................
 * Private cleanup function of new_Dpage() for memory allocation failures.
 *
 * Input:
 *  dp      Dpage *   The partially initialized Dpage descriptor.
 * Output:
 *  return  Dpage *   Allways NULL.
 */
static Dpage *dpmemerr(Dpage *dp)
{
  lprintf(stderr, "new_Dpage: Insufficient memory.\n");
  return del_Dpage(dp);
}

/*.......................................................................
 * Delete a Dpage descriptor and its contents.
 *
 * Input:
 *  dp      Dpage *   A descriptor originally returned by new_Dpage().
 * Output:
 *  return  Dpage *   Always NULL. Use like dp=del_Dpage(dp);
 */
Dpage *del_Dpage(Dpage *dp)
{
  if(dp) {
/*
 * Close and delete the scratch file.
 */
    if(dp->rio)
      dp->rio = del_Recio(dp->rio);
/*
 * Release memory from the visibility tree.
 */
    if(dp->ifs) {
      Dchan *dchan = dp->ifs->chan;
      if(dchan) {
	Dbase *dbase = dchan->base;
	if(dbase)
	  free(dbase);
	free(dchan);
      };
      free(dp->ifs);
    };
/*
 * Release memory allocated to the visibility I/O buffer.
 */
    if(dp->cvis)
      free(dp->cvis);
/*
 * Delete the container.
 */
    free(dp);
  };
  return NULL;
}

/*.......................................................................
 * Write the previously read portion of the visibility buffer dp->cvis[]
 * to the scratch file.
 *
 * Input:
 *  dp     Dpage *  The data paging descriptor.
 *  ut      long    The (0-relative) index of the integration to be
 *                  written.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int dp_write(Dpage *dp, long ut)
{
  if(dp_error(dp, "dp_write"))
    return 1;
/*
 * Check the integration index.
 */
  if(ut < 0 || ut >= dp->ntime) {
    lprintf(stderr, "dp_write: Integration index out of range.\n");
    return 1;
  };
  dp->ut = ut;
/*
 * Position the file if necessary.
 */
  if(rec_seek(dp->rio, ut, dp->first * sizeof(Cvis))) {
    dp->ioerr = 1;
    return 1;
  };
/*
 * Write the buffer to the scratch file.
 */
  if(rec_write(dp->rio, dp->nbuff, sizeof(Cvis), &dp->cvis[dp->first]) < dp->nbuff) {
    lprintf(stderr, "dp_write: Error writing to scratch file.\n");
    dp->ioerr = 1;
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Read a given portion of an integration into the corresponding
 * portion of the visibility buffer dp->cvis[].
 * The range of visibilities read will be that set by dp->first and
 * dp->nbuff.
 *
 * Input:
 *  dp     Dpage *  The data paging descriptor.
 *  ut      long    The (0-relative) index of the integration to be read.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int dp_read(Dpage *dp, long ut)
{
/*
 * Check validity of dp.
 */
  if(dp_error(dp, "dp_read"))
    return 1;
/*
 * Check integration index.
 */
  if(ut < 0 || ut >= dp->ntime) {
    lprintf(stderr, "dp_read: Integration index out of range.\n");
    return 1;
  };
  dp->ut = ut;
/*
 * Position the file if necessary.
 */
  if(rec_seek(dp->rio, ut, dp->first * sizeof(Cvis))) {
    dp->ioerr = 1;
    return 1;
  };
/*
 * Read from the scratch file into the buffer.
 */
  if(rec_read(dp->rio, dp->nbuff, sizeof(Cvis), &dp->cvis[dp->first]) < dp->nbuff) {
    lprintf(stderr, "dp_read: Error reading from scratch file.\n");
    dp->ioerr = 1;
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Set the range of spectral-line channels whose visibilities are to be
 * read and written in subsequent reads and writes.
 *
 * Input:
 *  dp   Dpage *   The paging file descriptor.
 *  ca     int     The index of the first spectral-line channel.
 *  cb     int     The index of the last spectral-line channel.
 * Output:
 *  return int     0 - OK.
 *                 1 - Error.
 */
int dp_crange(Dpage *dp, int ca, int cb)
{
  long ainc;   /* Increment to index of start visibility */
  long binc;   /* Increment to index of end visibility */
/*
 * Quietly ignore range setting if dp is NULL.
 */
  if(dp==NULL)
    return 0;
/*
 * Check the validity of dp.
 */
  if(dp_error(dp, "dp_crange"))
    return 1;
/*
 * Ensure that all ranges increase from a to b.
 */
  if(ca>cb) {int tmp=ca;ca=cb;cb=tmp;};
/*
 * Sanity check the given range.
 */
  if(ca<0 || cb>=dp->nchan) {
    lprintf(stderr, "dp_crange: Out of range spectral-line channel indexes.\n");
    return 1;
  };
/*
 * Determine the required increments to the indexes of the start and end
 * of the visibility range.
 */
  ainc = (ca - dp->ca) * dp->coff;
  binc = (cb - dp->cb) * dp->coff;
/*
 * Modify the read/write visibility range.
 */
  dp->first += ainc;
  dp->nbuff += binc - ainc;
/*
 * Record the changed channel range.
 */
  dp->ca = ca;
  dp->cb = cb;
  return 0;
}

/*.......................................................................
 * Set the range of IFs whose visibilities are to be
 * read and written in subsequent reads and writes.
 *
 * Input:
 *  dp   Dpage *   The paging file descriptor.
 *  ifa     int    The index of the first IF.
 *  ifb     int    The index of the last IF.
 * Output:
 *  return int     0 - OK.
 *                 1 - Error.
 */
int dp_irange(Dpage *dp, int ifa, int ifb)
{
  long ainc;   /* Increment to index of start visibility */
  long binc;   /* Increment to index of end visibility */
/*
 * Quietly ignore range setting if dp is NULL.
 */
  if(dp==NULL)
    return 0;
/*
 * Check the validity of dp.
 */
  if(dp_error(dp, "dp_irange"))
    return 1;
/*
 * Ensure that all ranges increase from a to b.
 */
  if(ifa>ifb) {int tmp=ifa;ifa=ifb;ifb=tmp;};
/*
 * Sanity check the given range.
 */
  if(ifa<0 || ifb>=dp->nif) {
    lprintf(stderr, "dp_irange: Out of range IF indexes.\n");
    return 1;
  };
/*
 * Determine the required increments to the indexes of the start and end
 * of the visibility range.
 */
  ainc = (ifa - dp->ia) * dp->ioff;
  binc = (ifb - dp->ib) * dp->ioff;
/*
 * Modify the read/write visibility range.
 */
  dp->first += ainc;
  dp->nbuff += binc - ainc;
/*
 * Record the changed IF range.
 */
  dp->ia = ifa;
  dp->ib = ifb;
  return 0;
}

/*.......................................................................
 * Set the range of stokes parameters whose visibilities are to be
 * read and written in subsequent reads and writes.
 *
 * Input:
 *  dp   Dpage *   The paging file descriptor.
 *  sa     int     The index of the first IF.
 *  sb     int     The index of the last IF.
 * Output:
 *  return int     0 - OK.
 *                 1 - Error.
 */
int dp_srange(Dpage *dp, int sa, int sb)
{
  long ainc;   /* Increment to index of start visibility */
  long binc;   /* Increment to index of end visibility */
/*
 * Quietly ignore range setting if dp is NULL.
 */
  if(dp==NULL)
    return 0;
/*
 * Check the validity of dp.
 */
  if(dp_error(dp, "dp_srange"))
    return 1;
/*
 * Ensure that all ranges increase from a to b.
 */
  if(sa>sb) {int tmp=sa;sa=sb;sb=tmp;};
/*
 * Sanity check the given range.
 */
  if(sa<0 || sb>=dp->npol) {
    lprintf(stderr, "dp_srange: Out of range stokes indexes.\n");
    return 1;
  };
/*
 * Determine the required increments to the indexes of the start and end
 * of the visibility range.
 */
  ainc = (sa - dp->sa) * dp->soff;
  binc = (sb - dp->sb) * dp->soff;
/*
 * Modify the read/write visibility range.
 */
  dp->first += ainc;
  dp->nbuff += binc - ainc;
/*
 * Record the changed stokes range.
 */
  dp->sa = sa;
  dp->sb = sb;
  return 0;
}

/*.......................................................................
 * Set the range of baselines whose visibilities are to be read and
 * written in subsequent reads and writes.
 *
 * Input:
 *  dp   Dpage *   The paging file descriptor.
 *  ba     int     The index of the first IF.
 *  bb     int     The index of the last IF.
 * Output:
 *  return int     0 - OK.
 *                 1 - Error.
 */
int dp_brange(Dpage *dp, int ba, int bb)
{
  long ainc;   /* Increment to index of start visibility */
  long binc;   /* Increment to index of end visibility */
/*
 * Quietly ignore range setting if dp is NULL.
 */
  if(dp==NULL)
    return 0;
/*
 * Check the validity of dp.
 */
  if(dp_error(dp, "dp_brange"))
    return 1;
/*
 * Ensure that all ranges increase from a to b.
 */
  if(ba>bb) {int tmp=ba;ba=bb;bb=tmp;};
/*
 * Sanity check the given range.
 */
  if(ba<0 || bb>=dp->nbase) {
    lprintf(stderr, "dp_brange: Out of range baseline indexes.\n");
    return 1;
  };
/*
 * Determine the required increments to the indexes of the start and end
 * of the visibility range.
 */
  ainc = (ba - dp->ba) * dp->boff;
  binc = (bb - dp->bb) * dp->boff;
/*
 * Modify the read/write visibility range.
 */
  dp->first += ainc;
  dp->nbuff += binc - ainc;
/*
 * Record the changed baseline range.
 */
  dp->ba = ba;
  dp->bb = bb;
  return 0;
}

/*.......................................................................
 * Check the validity of the Dpage descriptor.
 *
 * Input:
 *  dp     Dpage *   The data paging file descriptor.
 *  fname   char *   The name of the calling function.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
static int dp_error(Dpage *dp, const char *fname)
{
/*
 * No descriptor?
 */
  if(dp==NULL) {
    lprintf(stderr, "%s: Intercepted NULL Dpage descriptor.\n", fname);
    return 1;
  };
/*
 * Previous I/O error?
 */
  if(dp->ioerr)
    return 1;
  return 0;
}

/*.......................................................................
 * Clear the whole input/output buffer to prepare for a new integration.
 *
 * Input:
 *  dp    Dpage *  The data paging file descriptor.
 *  ut      int    If you are clearing the buffer to start building a new
 *                 integration, specify the index of that integration here.
 *                 Otherwise specify -1 to say that there is no integration
 *                 associated with the data in dp->cvis[].
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int dp_clear(Dpage *dp, int ut)
{
  Cvis *cvis;   /* Pointer to array of dp->nvis visibilities */
  long nvis;    /* Number of visibilities left to be cleared */
/*
 * Quietly ignore this setup call if dp is NULL.
 */
  if(dp==NULL)
    return 0;
/*
 * Check the validity of dp.
 */
  if(dp_error(dp, "dp_clear"))
    return 1;
/*
 * Clear the buffer.
 */
  nvis = dp->nvis;
  cvis = dp->cvis;
  while(nvis--)
    *cvis++ = zero_vis;
  dp->ut = ut;
  return 0;
}

/*.......................................................................
 * Make sure that a data paging file is up to date by flushing all
 * I/O.
 *
 * Input:
 *  dp    Dpage *  The data paging file descriptor.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int dp_flush(Dpage *dp)
{
  return dp ? rec_flush(dp->rio) : 0;
}

/*.......................................................................
 * Apply self-cal and resoff corrections to the intergration record
 * currently in the ob->dp I/O buffer. Corrections are only applied to
 * the last region read. This is determined by looking at the ranges
 * and integration record number recorded in ob->dp, along with the
 * number of baselines in the associated sub-array of *ob. Note that
 * dp->ut must be set to the index of the integration being processed.
 * dp_read(), dp_write() and dp_clear() initialize dp->ut correctly.
 *
 * Input:
 *  ob    Observation *   The descriptor of the observation.
 * Output:
 *  return        int     0 - OK.
 *                        1 - Error.
 */
int dp_cal(Observation *ob)
{
  Integration *integ;/* The descriptor of the integration being processed */
  Subarray *sub;     /* The descriptor of the sub-array of the integration */
  Dpage *dp;         /* The descriptor of the uvdata.scr file */
  Dif *ifs;          /* Pointer to the visibilities of the IF being processed */
  int cif;           /* The index of the IF being processed */
  int fc;            /* The index of the channel being processed */
  int base;          /* The index of the baseline being processed */
  int pol;           /* The index of the polarization being processed */
  int bb;            /* The index of the last baseline to be calibrated */
/*
 * Sanity check the arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "dp_cal"))
    return 1;
/*
 * Get the uvdata file descriptor.
 */
  dp = ob->dp;
  if(dp->ut < 0) {
    lprintf(stderr, "dp_cal: Invalid integration.\n");
    return 1;
  };
/*
 * Get the descriptors of the integration and the sub-array from which
 * it was taken.
 */
  integ = ob->rec[dp->ut].integ;
  sub = integ->sub;
/*
 * The Dpage I/O buffer is dimensioned to hold the number of baselines
 * in the biggest sub-array. Restrict calibrations to the actual number
 * of baselines in the current sub-array.
 */
  bb = dp->bb < sub->nbase ? dp->bb : sub->nbase-1;
/*
 * Are there any baselines to calibrate?
 */
  if(dp->ba > bb)
    return 0;
/*
 * Loop through each IF in the I/O buffer.
 */
  for(cif=dp->ia, ifs=dp->ifs+cif; cif<=dp->ib; cif++,ifs++) {
    Telcor *tcor = integ->icor[cif].tcor; /* IF telescope corrections */
    Baseline *bptr;                       /* Baseline descriptor */
/*
 * All polarizations and spectral-line channels receive the same corrections
 * on a given baseline, so do one baseline at a time.
 */
    for(base=dp->ba, bptr=sub->base+base; base<=bb; base++,bptr++) {
      Bascor *bcor = &bptr->bcor[cif];    /* IF baseline corrections */
      int ta = bptr->tel_a;  /* Index of first telescope on baseline */
      int tb = bptr->tel_b;  /* Index of second telescope on baseline */
/*
 * Combine selfcal and resoff contributions to the amp and phase corrections.
 */
      float amp_cor = tcor[ta].amp_cor * tcor[tb].amp_cor * bcor->amp_cor;
      float phs_cor = tcor[ta].phs_cor - tcor[tb].phs_cor + bcor->phs_cor;
/*
 * Pre-compute the cos and sin of the phase correction for use in
 * correcting the complex representation of the visibilities.
 */
      float cosphi = cos(phs_cor);
      float sinphi = sin(phs_cor);
/*
 * Determine whether the correction for this baseline is flagged.
 */
      int bad_cor = tcor[ta].bad || tcor[tb].bad;
      Dchan *dchan;
/*
 * The amplitude correction must be finite.
 */
      if(amp_cor <= 0.0f)
	amp_cor = 1.0f;
/*
 * Apply corrections to all polarizations of all spectral-line channels.
 */
      for(fc=dp->ca, dchan=ifs->chan+fc; fc<=dp->cb; fc++,dchan++) {
	Cvis *cvis;
	for(pol=dp->sa, cvis=dchan->base[base].pol+pol; pol<=dp->sb;
	    pol++, cvis++) {
	  float re = cvis->re;
	  float im = cvis->im;
	  float wt = cvis->wt;
/*
 * The complex correction is: (x+iy) * amp_cor * exp(i.phs_cor).
 */
	  cvis->re = amp_cor * (re * cosphi - im * sinphi); /* Real */
	  cvis->im = amp_cor * (re * sinphi + im * cosphi); /* Imaginary */
	  cvis->wt = (bad_cor && wt>0.0f ? -wt:wt) / (amp_cor * amp_cor);
	             /* Weight = 1/amp_err^2 */
	};
      };
    };
  };
/*
 * The corrected visibilities are now in the ob->dp I/O buffer.
 */
  return 0;
}

/*.......................................................................
 * Apply the current stream shift to the data in the ob->dp I/O buffer.
 *
 * Input:
 *  ob    Observation *   The descriptor of the observation.
 * Output:
 *  return        int     0 - OK.
 *                        1 - Error.
 */
int dp_shift(Observation *ob)
{
  Integration *integ;/* The descriptor of the integration being processed */
  Subarray *sub;     /* The descriptor of the sub-array of the integration */
  Dpage *dp;         /* The descriptor of the uvdata.scr file */
  Dif *dif;          /* Pointer to the visibilities of the IF being processed */
  If *ifp;           /* The IF descriptor of the IF being processed */
  int cif;           /* The index of the IF being processed */
  int fc;            /* The index of the channel being processed */
  int base;          /* The index of the baseline being processed */
  int pol;           /* The index of the polarization being processed */
  int bb;            /* The index of the last baseline to be calibrated */
/*
 * Sanity check the arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "dp_shift"))
    return 1;
/*
 * Get the uvdata file descriptor.
 */
  dp = ob->dp;
  if(dp->ut < 0) {
    lprintf(stderr, "dp_shift: Invalid integration.\n");
    return 1;
  };
/*
 * Get the descriptors of the integration and the sub-array from which
 * it was taken.
 */
  integ = ob->rec[dp->ut].integ;
  sub = integ->sub;
/*
 * The Dpage I/O buffer is dimensioned to hold the number of baselines
 * in the biggest sub-array. Restrict processing to the actual number
 * of baselines in the current sub-array.
 */
  bb = dp->bb < sub->nbase ? dp->bb : sub->nbase-1;
/*
 * Are there any baselines to shift?
 */
  if(dp->ba > bb)
    return 0;
/*
 * If there is no shift to be applied, avoid the overhead.
 */
  if(ob->geom.east == 0.0f && ob->geom.north == 0.0f)
    return 0;
/*
 * Loop through each IF in the I/O buffer.
 */
  for(cif=dp->ia, dif=dp->ifs+cif, ifp=ob->ifs+cif; cif<=dp->ib;
      cif++, dif++, ifp++) {
/*
 * Determine the frequency of the start channel and the increment between
 * channels.
 */
    float freq = ifp->freq + dp->ca * ifp->df;
    float df = ifp->df;
/*
 * Baselines.
 */
    Visibility *vis;                      /* Visibility of baseline, for UVW */
    for(base=dp->ba, vis=integ->vis+base; base<=bb; base++, vis++) {
/*
 * To perform the shift, we need cos() and sin() of the phase shift for
 * each frequency channel. Fortunately the phase shift increments by
 * a fixed amount (for a given U,V coordinate) between frequency channels
 * so we can use a recurrence relation to determine the sin() and cos() for
 * all channels from the sin() and cos() of the phase shift for the
 * start frequency channel, and the phase shift increment between channels.
 */
      float shift = (ob->geom.east * twopi * vis->u +
		     ob->geom.north * twopi * vis->v);
      float sinphs = sin(freq * shift);
      float cosphs = cos(freq * shift);
      float sininc = sin(df * shift);
      float cosinc = cos(df * shift);
/*
 * Frequency channels.
 */
      Dchan *dchan;
      for(fc=dp->ca, dchan=dif->chan+fc; fc<=dp->cb; fc++,dchan++) {
	Cvis *cvis;
/*
 * Stokes.
 */
	for(pol=dp->sa, cvis=dchan->base[base].pol+pol; pol<=dp->sb;
	    pol++, cvis++) {
	  if(cvis->wt != 0.0f) {
	    float re = cvis->re;
	    float im = cvis->im;
	    cvis->re = re * cosphs - im * sinphs;
	    cvis->im = re * sinphs + im * cosphs;
	  };
	};
/*
 * Determine sin(phase+increment) and cos(phase+increment) for the next channel.
 */
	{
	  float cphs = cosphs*cosinc - sinphs*sininc;
	  float sphs = cosphs*sininc + sinphs*cosinc;
	  cosphs = cphs;
	  sinphs = sphs;
	};
      };
    };
  };
/*
 * The shifted visibilities are now in the ob->dp I/O buffer.
 */
  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "logio.h"
#include "obs.h"
#include "obedit.h"
#include "modeltab.h"

static Observation *obalerr(Observation *ob);
static Stokes *del_Stokes(Observation *ob);
static Stokes *new_Stokes(Observation *ob, int npol);

/*
 * Set the size of the model hash table. This should be a prime number.
 */
#define MTAB_SIZE 113
/*
 * Set the blocking size for allocating model nodes in the model hash
 * table. This specifies how many model nodes should be allocated at
 * once when expanding the table.
 */
#define MTAB_BLK 50

/*.......................................................................
 * Read the details of an observation from a file, and record them
 * in an Observation structure.
 *
 * Input:
 *  name    const char *  The name of the input file.
 *  binwid      double    The integration time to bin the data into (seconds).
 *                        No binning will be performed if binwid < 1 second.
 *  scatter        int    If true replace data weights with weights derived
 *                        from the data scatter.
 *  keepant        int    If true, allocate space for all antennas and
 *                        associated baselines. If false, discard all
 *                        antennas and baselines that don't have
 *                        visibilities associated with them.
 *  cl          Chlist *  The required list of initial channel ranges, or
 *                        NULL for the default.
 *  stokes      Stokes    The initial stokes parameter, or NO_POL for
 *                        the default (I,RR,LL or the first polarization).
 * Output:
 *  return Observation *  The descriptor of the observation, or NULL
 *                        on error.
 */
Observation *new_Observation(const char *name, double binwid, int scatter,
			     int keepant, Chlist *cl, Stokes stokes)
{
  Observation *ob;  /* The descriptor to return */
/*
 * Assume that the file is a UV-FITS file.
 */
  ob = uvf_read(name, binwid, scatter, keepant);
/*
 * Select the initial processing stream.
 */
  if(ob) {
/*
 * We have now acquired the raw data.
 */
    ob->state = OB_DATA;
/*
 * Index the sub-array integrations, in time order in ob->rec.
 */
    if(ini_Intrec(ob))
      return del_Observation(ob);
/*
 * Clear the model.
 */
    clrmod(ob, 1, 1, 1);
/*
 * Automatically make a selection only if, either a non-default channel range
 * and polarization is given, or if there is only one possible selection.
 */
    if((cl!=NULL && stokes!=NO_POL) || (ob->npol==1 && ob->nctotal==1))
      ob_select(ob, 0, cl, stokes);
  };
  return ob;
}

/*.......................................................................
 * Allocate an uninitialised Observation structure and its constituent
 * arrays. The passed parameters will be installed in their equivalent
 * slots in the container.
 *
 * The observation returned may be sent to del_Observation() without
 * further initialisation.
 *
 * IMPORTANT: Before attempting to use the returned Observation descriptor,
 * be sure to call ini_Subarray() to fill each sub-array and then call
 * ini_Intrec() to index all sub-array integrations in ob->rec.
 *
 * Input:
 *  ob Observation * An existing descriptor to be re-sized, or NULL.
 *  nrec       int   The total number of integrations in the observation.
 *  nbmax      int   The max number of used baselines in any sub-array.
 *  nsub       int   The number of telescope sub-arrays in the observation.
 *  nif        int   The number of IFs in the observation.
 *  npol       int   The number of polarizations or stokes parameters.
 *  nchan      int   The number of spectral-line channels.
 * Output:
 *  return Observation *  The allocated observation container plus its
 *                        constituent arrays, or NULL on error.
 */
Observation *Obs_alloc(Observation *ob, int nrec, int nbmax, int nsub,
		       int nif, int npol, int nchan)
{
/*
 * Attempt to allocate the Observation container.
 */
  if(ob==NULL) {
    ob = (Observation *) malloc(sizeof(Observation));
    if(ob == NULL)
      return obalerr(ob);
/*
 * Mark the container as allocated.
 */
    ob->state = OB_ALLOC;
/*
 * Clear the container.
 */
    ob->nhist = 0;
    ob->nsub  = 0;
    ob->nrec  = 0;
    ob->nif   = 0;
    ob->npol  = 0;
    ob->nchan = 0;
    ob->nbmax = 0;
    ob->nctotal = 0;
    ob->have_inttim = 0;
/*
 * Clear all date fields.
 */
    ob->date.year = 0;
    ob->date.utc_ref = 0.0;
    ob->date.ut = 0.0;
    ob->date.app_st = 0.0;
    ob->date.cav_tim = 0.0;
    ob->date.iav_tim = 0.0;
/*
 * We know of no velocity info yet.
 */
    ob->vel.velref = 0;
    ob->vel.altrval = 0.0;
    ob->vel.altrpix = 0.0;
    ob->vel.restfreq = 0.0;
/*
 * Default to use the SIN projection.
 */
    ob->proj = PRJ_SIN;
/*
 * Initialize the stream record for a single channel, single IF, single
 * polarization data set.
 */
    ob->stream.cif = 0;
    ob->stream.cl = NULL;
    ob->stream.pol.type = NO_POL;
    ob->stream.pol.pa = ob->stream.pol.pb = -1;
    ob->stream.pol.getpol = 0;
    ob->stream.uvscale = 1.0f;
/*
 * Initialize the UV geometry transformations for no initial tranformation.
 */
    ob->geom.east = 0.0f;
    ob->geom.north = 0.0f;
    ob->geom.uvangle = 0.0f;
    ob->geom.wtscale = 1.0f;
/*
 * Initialize the zero-spacing flux to undefined.
 */
    ob->uvzero.amp = 0.0f;
    ob->uvzero.modamp = 0.0f;
    ob->uvzero.wt = 0.0f;
/*
 * Clear the source container.
 */
    ob->source.name[0] = '\0';
    ob->source.epoch = 0.0;
    ob->source.ra = 0.0;
    ob->source.dec = 0.0;
    ob->source.app_ra = 0.0;
    ob->source.app_dec = 0.0;
    ob->source.tot_flux = 0.0;
    ob->source.have_obs = 0;
    ob->source.obsra = 0.0;
    ob->source.obsdec = 0.0;
    ob->source.east = 0.0;
    ob->source.north = 0.0;
/*
 * No model visibilities yet.
 */
    ob->hasmod=0;
/*
 * Make all pointer members in the container NULL so that hereafter
 * the container can be sent to del_Observation().
 */
    ob->sub  = 0;
    ob->rec  = 0;
    ob->his  = 0;
    ob->ifs  = 0;
    ob->pols = 0;
    ob->dp   = 0;
    ob->ip   = 0;
    ob->uvp  = 0;
    ob->model= 0;
    ob->newmod = 0;
    ob->cmodel= 0;
    ob->cnewmod = 0;
    ob->mtab = NULL;
    ob->obed = 0;
    ob->ab = NULL;
    ini_Obhead(ob, NULL, NULL, NULL, NULL, NULL, NULL, 0.0);
/*
 * Create a scratch file to store FITS history in.
 */
    ob->his = new_Recio("history.scr", IS_SCR, 0, 80L);
    if(ob->his==NULL)
      return del_Observation(ob);
/*
 * Allocate and initialize the container of the deferred scratch file
 * editing free-list.
 */
    ob->obed = new_Obedit(ob);
    if(ob->obed == NULL)
      return del_Observation(ob);
  };
/*
 * Mark the descriptor as allocated.
 */
  ob->state = OB_ALLOC;
/*
 * (Re-)allocate an array of nsub sub-array descriptors.
 * Each new element must be externally filled via ini_Subarray() before
 * the Observation descriptor can be used.
 */
  if(new_Subarray(ob, nsub)==NULL)
    return obalerr(ob);
/*
 * (Re-)allocate a new Intrec array.
 */
  if(new_Intrec(ob, nrec)==NULL)
    return obalerr(ob);
/*
 * (Re-)allocate and initialize an array of nif IF descriptors and
 * their contents.
 */
  if(new_If(ob, nif)==NULL)
    return obalerr(ob);
/*
 * Allocate an array of npol stokes descriptors.
 */
  if(new_Stokes(ob, npol)==NULL)
    return obalerr(ob);
/*
 * Create a UVDATA paging file to hold the raw data if one does not already
 * exist.
 */
  if(ob->dp==NULL) {
    ob->dp = new_Dpage(nrec, nbmax, nchan, nif, npol);
    if(ob->dp==NULL)
      return del_Observation(ob);
/*
 * If a UVDATA paging file does exist then it must have the same dimensions
 * as now required.
 */
  } else {
    Dpage *dp = ob->dp;
    if(dp->ntime != nrec || dp->nbase != nbmax || dp->nchan != nchan ||
       dp->nif != nif || dp->npol != npol) {
      lprintf(stderr, "Obs_alloc: Can't re-size uvdata.scr files.\n");
      return del_Observation(ob);
    };
  };
/*
 * Delete obsolete IF paging file if not of the required size.
 */
  if(ob->ip) {
    IFpage *ip = ob->ip;
    if(ip->nif != nif || ip->nbase != nbmax || ip->ntime != nrec)
      ob->ip = del_IFpage(ob->ip);
  };
/*
 * Create an IF scratch file, but only if required and one does not
 * already exist.
 */
  if(ob->ip==NULL && nif>1) {
    ob->ip = new_IFpage(nif, nbmax, nrec);
    if(ob->ip==NULL)
      return del_Observation(ob);
  };
/*
 * Delete obsolete UV-model paging file if not of the required size.
 */
  if(ob->uvp) {
    UVpage *uvp = ob->uvp;
    if(uvp->nif != nif || uvp->nbase != nbmax || uvp->ntime != nrec)
      ob->uvp = del_UVpage(ob->uvp);
  };
/*
 * Create a scratch file to store UV-model representations for each
 * IF, but only if required and one does not already exist.
 */
  if(ob->uvp==NULL && nif>1) {
    ob->uvp = new_UVpage(nrec, nbmax , nif);
    if(ob->uvp==NULL)
      return del_Observation(ob);
  };
/*
 * Allocate empty established and tentative models if necessary.
 */
  if(ob->model==NULL) {
    ob->model = new_Model();
    if(ob->model==NULL)
      return del_Observation(ob);
  };
  if(ob->newmod==NULL) {
    ob->newmod = new_Model();
    if(ob->newmod==NULL)
      return del_Observation(ob);
  };
/*
 * Allocate empty established and tentative continuum models if necessary.
 */
  if(ob->cmodel==NULL) {
    ob->cmodel = new_Model();
    if(ob->cmodel==NULL)
      return del_Observation(ob);
  };
  if(ob->cnewmod==NULL) {
    ob->cnewmod = new_Model();
    if(ob->cnewmod==NULL)
      return del_Observation(ob);
  };
/*
 * Allocate a new model table?
 */
  if(!ob->mtab) {
    ob->mtab = new_ModelTable(MTAB_SIZE, MTAB_BLK);
    if(!ob->mtab)
      return del_Observation(ob);
  };
/*
 * Allocate a new container object for antenna and primary beams?
 */
  if(!ob->ab) {
    ob->ab = new_AntennaBeams();
    if(!ob->ab)
      return del_Observation(ob);
  };
/*
 * Install the new dimensions.
 */
  ob->nsub = nsub;
  ob->nrec = nrec;
  ob->nif  = nif;
  ob->npol  = npol;
  ob->nchan = nchan;
  ob->nbmax = nbmax;
  ob->nctotal = nif * nchan;
/*
 * Mark the descriptor as successfully allocated.
 */
  ob->state = OB_ALLOC;
/*
 * Return the initialized observation descriptor.
 */
  return ob;
}

/*.......................................................................
 * Private function of Obs_alloc, used to report memory allocation
 * failures, delete the observation and return NULL. It may thus be used
 * as a return value.
 */
static Observation *obalerr(Observation *ob)
{
  lprintf(stderr, "Insufficient memory for new observation.\n");
  return del_Observation(ob);
}

/*.......................................................................
 * Destroy an Observation data structure. This mainly involves
 * de-allocating memory allocated in function new_Observation().
 * This function always returns a NULL Observation pointer.
 */
Observation *del_Observation(Observation *ob)
{
  if(ob) {
/*
 * Delete the index array of file record integrations.
 */
    del_Intrec(ob);
/*
 * Free the array of sub-array descriptors and their contents.
 */
    del_Subarray(ob);
/*
 * Free IF descriptors.
 */
    del_If(ob);
/*
 * Free polarization descriptors.
 */
    del_Stokes(ob);
/*
 * Delete the uvdata scratch file.
 */
    ob->dp = del_Dpage(ob->dp);
/*
 * Delete the IF scratch file.
 */
    ob->ip = del_IFpage(ob->ip);
/*
 * Delete the history scratch file.
 */
    ob->his = del_Recio(ob->his);
/*
 * Delete the UV model scratch file.
 */
    ob->uvp = del_UVpage(ob->uvp);
/*
 * Delete the lists of map-plane model components.
 */
    ob->model = del_Model(ob->model);
    ob->newmod = del_Model(ob->newmod);
    ob->cmodel = del_Model(ob->cmodel);
    ob->cnewmod = del_Model(ob->cnewmod);
/*
 * Delete the table of models.
 */
    ob->mtab = del_ModelTable(ob->mtab);
/*
 * Delete the deferred editing free-list.
 */
    ob->obed = del_Obedit(ob);
/*
 * Delete misc FITS header strings.
 */
    clr_Obhead(ob);
/*
 * Delete the antenna/primary beam container object.
 * Note that this has to be done after calling del_Subarray(),
 * since the station structures of each subarray contain references
 * to objects which are deleted by this function.
 */
    ob->ab = del_AntennaBeams(ob->ab);
/*
 * Finally free the observation descriptor itself.
 */
    free(ob);
  };
  return NULL;
}

/*.......................................................................
 * (Re-)allocate an array of npol polarization descriptors.
 *
 * Input:
 *  ob  Observation *  The descriptor of the containing observation.
 *  npol        int    The new number of polarizations required.
 * Output:
 *  return   Stokes *  The revised stokes array, or NULL on error.
 */
static Stokes *new_Stokes(Observation *ob, int npol)
{
  Stokes *pol;   /* Pointer into ob->pols[] */
/*
 * Allocate a new array?
 */
  if(ob->pols==NULL) {
    ob->pols = malloc(sizeof(Stokes) * npol);
    if(ob->pols == NULL) {
      lprintf(stderr, "new_Stokes: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing array?
 */
  else if(ob->npol != npol) {
    pol = realloc(ob->pols, sizeof(Stokes) * npol);
    if(ob->pols == NULL) {
      lprintf(stderr, "new_Stokes: Insufficient memory.\n");
      return NULL;
    };    
  };
/*
 * Return the revised array.
 */
  return ob->pols;
}

/*.......................................................................
 * Delete the array of polarization descriptors of an observation.
 *
 * Input:
 *  ob  Observation *   The descriptor of the observation containing the
 *                      array to be deleted.
 * Output:
 *  return   Stokes *   The deleted array - ie. allways NULL.
 */
static Stokes *del_Stokes(Observation *ob)
{
  if(ob && ob->pols) {
    free(ob->pols);
    ob->pols = NULL;
  };
  return NULL;
}

/*.......................................................................
 * Clear the current history by reseting ob->nhist to 0.
 *
 * Input:
 *  ob     Observation *  The descriptor of the observation whos history is
 *                        to be cleared.
 * Output:
 *  return         int    0 - OK.
 *                        1 - Error.
 */
int clr_hist(Observation *ob)
{
  if(!ob_ready(ob, OB_ALLOC, "clr_hist"))
    return 1;
/*
 * Reset the recorded number of history records to zero.
 */
  ob->nhist = 0;
  return 0;
}

/*.......................................................................
 * To allow primary beam calculations in observations that don't contain
 * OBSRA and OBSDEC keywords, specify the Right Ascension and Declination
 * of the original pointing center of the observation.
 *
 * Input:
 *  ob   Observation *  The container of the observation.
 *  obsra     double    The center Right Ascension of the pointing.
 *                      (radians)
 *  obsdec    double    The center Declination of the pointing.
 *                      (radians)
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int set_obs_radec(Observation *ob, double obsra, double obsdec)
{
/*
 * Check the arguments.
 */
  if(!ob_ready(ob, OB_ALLOC, "set_obs_radec"))
    return 1;
/*
 * Register the new pointing center.
 */
  ob->source.have_obs = 1;
  ob->source.obsra = obsra;
  ob->source.obsdec = obsdec;
/*
 * Work out the sin-projection offset of the pointing center from
 * the recorded source position.
 */
  ob->source.east = radec_to_l(ob->source.ra, ob->source.dec,
			       obsra, obsdec, ob->proj);
  ob->source.north = radec_to_m(ob->source.ra, ob->source.dec,
				obsra, obsdec, ob->proj);
  return 0;
}

/*.......................................................................
 * Get the offset of a given sky-projected position to the pointing center.
 *
 * Input:
 *  ob        Observation *  The container of the observation.
 *  x, y            float    The east,north map-projection offset of the
 *                           position from the current map center (radians).
 * Output:
 *  return          float    The distance of the specified position from
 *                           the pointing center (radians).
 */
float calc_pointing_offset(Observation *ob, float x, float y)
{
  float east, north;   /* The east and north offsets of the position from */
                       /*  the pointing center */
/*
 * Check the arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "calc_pointing_offset"))
    return 1;
/*
 * Sum up the various offsets from the pointing center.
 */
  east = ob->source.east + x - ob->geom.east;
  north = ob->source.north + y - ob->geom.north;
/*
 * Compute the corresponding radial offset from the pointing center.
 */
  return sqrt(east*east + north*north);
}


/* Observation "method" file.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include "logio.h"
#include "obs.h"
#include "obedit.h"
#include "modeltab.h"

static int ob_get_select(Observation *ob, Chlist *cl, Stokes pol);
static int bad_ob_chlist(Observation *ob, Chlist *def_cl, Chlist **if_cl);
static int iniIF(Observation *ob, int cif);

/*.......................................................................
 * Add a new line of history to observation 'ob'.
 * This involves writing the history line to the history scratch file and
 * updating the count of history lines in ob->nhist.
 *
 * Input:
 *   ob   Observation *  The observation to which the history pertains.
 *   hist        char *  The new line of history to be added. Only the
 *                       first 80 characters will be used.
 * Output:
 *   return       int    0 - OK.
 *                       1 - Error.
 */
int add_hist(Observation *ob, char *hisrec)
{
  char newhis[81];  /* The history line to be written */
  int hislen;       /* Length of input string */
/*
 * Determine the length of the input string.
 */
  hislen = strlen(hisrec);
/*
 * Copy up to 80 characters to newhis[].
 */
  strncpy(newhis, hisrec, 80);
/*
 * If the history line contains less than 80 characters pad with
 * spaces.
 */
  while(hislen < 80)
    newhis[hislen++] = ' ';
  newhis[80] = '\0';
/*
 * Append the history line to the history scratch file.
 */
  if(rec_seek(ob->his, ob->nhist, 0L) ||
     rec_write(ob->his, 80, sizeof(char), newhis) < 80)
    return 1;
/*
 * Increment the count of history lines.
 */
  ob->nhist++;
  return 0;
}

/*.......................................................................
 * Return the FITS name of the given spherical coordinate projection.
 *
 * Input:
 *  proj     Proj   The enumerator to be named.
 * Output:
 *  return   char * A pointer to an internal static string naming the
 *                  projection. "   " is returned if not recognised.
 */
char *Proj_name(Proj proj)
{
  switch(proj) {
  case PRJ_SIN:
    return "SIN";
    break;
  case PRJ_NCP:
    return "NCP";
    break;
  default:
    break;
  };
  return "   ";
}

/*.......................................................................
 * Return the enumerator associated with a given FITS coordinate
 * projection name.
 *
 * Input:
 *  name     char *   The projection name to look up. This should be
 *                    a 3 character all upper-case name.
 * Output:
 *  return   Proj     The enumerator found, or PRJ_NON if not found.
 */
Proj name_Proj(char *name)
{
  if(name) {
    if(strcmp(name, "SIN")==0)
      return PRJ_SIN;
    else if(strcmp(name, "NCP")==0)
      return PRJ_NCP;
  };
  return PRJ_NON;
}

/*.......................................................................
 * Make a given IF the current IF by reading it from the IF paging file
 * into its parent Observation and applying accumulated corrections.
 * Also read the associated UV model.
 *
 * In the special case of the requested IF already being the the current
 * IF, no operations will be performed and the function will return
 * sucessfully.
 *
 * If there is no IF paging file or no UV model paging file, then this
 * function may be called only if the requested IF is the one in memory.
 *
 * Input:
 *  ob    Observation *   The observation to read into.
 *  cif           int     The index of the IF to be read.
 * Output:
 *  return        int     0 - OK (ob->state=OB_GETIF).
 *                        1 - Error.
 */
int getIF(Observation *ob, int cif)
{
  Intrec *rec; /* Pointer into integration record array ob->rec. */
  IFpage *ip;  /* The IF paging descriptor */
  int ut;      /* The index of the current integration */
  int base;    /* The index of the current baseline */
  int nodata;  /* True if no channels of the IF are selected */
/*
 * Check validity of the descriptor.
 */
  if(!ob_ready(ob, OB_SELECT, "getIF"))
    return 1;
/*
 * Check that the requested IF is in range.
 */
  if(cif < 0 || cif >= ob->nif) {
    lprintf(stderr, "getIF: IF %d is unavailable.\n", cif);
    return 1;
  };
/*
 * If the requested IF is already in memory, do nothing further.
 * Note that no IF paging file is created for single IF data-sets,
 * so their visibilities are always in memory, but ob->state can
 * still get demoted from OB_GETIF.
 */
  if((ob_ready(ob, OB_GETIF, NULL) && cif==ob->stream.cif) || ob->nif==1) {
    ob->state = OB_GETIF;
    return 0;
  };
/*
 * Are any channels selected in the requested IF.
 */
  nodata = ob->ifs[cif].cl == NULL;
/*
 * Get the IF paging file descriptor.
 */
  ip = ob->ip;
/*
 * If there is no IF paging file then only the current IF is accessible.
 */
  if(ip==NULL || ob->uvp==NULL) {
    lprintf(stderr,
	    "getIF: There is no %s paging file to retrieve IF %d from.\n",
	    ip==NULL ? "IF":"UV model", cif);
    return 1;
  };
/*
 * Check the IF paging file descriptor for previous I/O errors.
 */
  if(ip_error(ip, "getIF"))
    return 1;
/*
 * Select the visibility range and IF to be read.
 */
  if(ip_range(ip, cif, 0, ob->nbmax-1))
    return 1;
/*
 * Set the observation state to selection status until the new IF
 * has been succesfully been acquired.
 */
  ob->state = OB_SELECT;
/*
 * Read each integration from the IF paging file.
 */
  rec = ob->rec;
  for(ut=0; ut<ob->nrec; ut++,rec++) {
    Integration *integ = rec->integ;
    Visibility *vis = integ->vis;
    int nbase = integ->sub->nbase;
    Dvis *dvis = ip->dvis;
/*
 * Read the next integration, or simply clear it if no channels have
 * been selected.
 */
    if(nodata ? ip_clear(ip) : ip_read(ip, ut))
      return 1;
/*
 * Copy each of its visibilities into the corresponding integration
 * of the Observation structure.
 */
    for(base=0; base<nbase; base++,vis++,dvis++) {
      vis->amp = dvis->amp;
      vis->phs = dvis->phs;
      if(dvis->wt > 0.0f) {        /* +ve weight means a good visibility */
	vis->wt = dvis->wt;
	vis->bad = 0;
      } else if(dvis->wt < 0.0f) { /* -ve weight means a flagged visibility */
	vis->wt = -dvis->wt;
	vis->bad = FLAG_BAD;
      } else {                     /* Zero weight means a deleted visibility */
	vis->wt = 0.0f;
	vis->bad = FLAG_DEL;
      };
    };
  };
/*
 * Read the associated UV model from the UV model paging file.
 */
  if(getmodel(ob, cif))
    return 1;
/*
 * Apply corrections, geometric shifts, weight scales etc to the data
 * and upgrade the observation state to OB_GETIF on success.
 */
  if(iniIF(ob, cif)) {
    ob->state = OB_SELECT;
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Private function of nextIF() and ob_select(), used to apply corrections
 * to visibilities that have just been read from the uvdata.scr scratch
 * file and record details of the new IF selection. The observation state
 * will be upgraded to OB_GETIF on success.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation.
 *  cif          int    The IF whose visibilities are to be corrected.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
static int iniIF(Observation *ob, int cif)
{
/*
 * We now have raw visibilities in memory.
 */
  ob->state = OB_RAWIF;
/*
 * Record the index of the new IF.
 */
  ob->stream.cif = cif;
/*
 * Set the factor used to scale UVW coords to wavelengths at the
 * frequency of the new IF.
 */
  ob->stream.uvscale = getuvscale(ob, cif);
/*
 * Apply corrections if there is anything to correct.
 */
  if(ob->ifs[cif].cl) {
/*
 * Apply self-cal corrections to the new data in the Observation structure.
 */
    if(app_Telcor(ob, cif))
      return 1;
/*
 * Apply resoff corrections to the new data in the Observation structure.
 */
    if(app_bcor(ob, cif))
      return 1;
/*
 * Apply pending edits.
 */
    if(app_Obedit(ob, cif))
      return 1;
/*
 * Shift the phase center of the data in the Observation structure.
 */
    if(ob->geom.east != 0.0f || ob->geom.north != 0.0f) {
      if(uvshift(ob, ob->geom.east, ob->geom.north))
	return 1;
    };
/*
 * Scale the weights?
 */
    if(ob->geom.wtscale != 1.0f) {
      float wtscale = ob->geom.wtscale;  /* Get the scale to be applied */
      int ut;                            /* Index of the next integration */
      int base;                          /* Index of the next baseline */
      int isub;                          /* Index of the next sub-array */
      Subarray *sub;                     /* Pointer to next sub-array */
/*
 * Loop over the integration arrays in each sub-array.
 */
      sub = ob->sub;
      for(isub=0; isub<ob->nsub; isub++,sub++) {
	Integration *integ = sub->integ;
	for(ut=0; ut<sub->ntime; ut++,integ++) {
	  Visibility *vis = integ->vis;
/*
 * Apply the scale factor to the visibilities of all baselines in the
 * current integration.
 */
	  for(base=0; base<sub->nbase; base++,vis++)
	    vis->wt *= wtscale;
	};
      };
    };
  };
/*
 * Record the new observation readiness state.
 */
  ob->state = OB_GETIF;
  return 0;
}

/*.......................................................................
 * Read the UV model of an IF into an Observation structure from the
 * UV model paging file.
 *
 * If the model is already in memory, nothing will be done, otherwise
 * if an IF is currently in memory it will be invalid after this call
 * and ob->state will be reverted to OB_SELECT to reflect this.
 *
 * Input:
 *  ob    Observation *   The observation to read into.
 *  cif           int     The index of the IF to be read.
 * Output:
 *  return        int     0 - OK.
 *                        1 - Error.
 */
int getmodel(Observation *ob, int cif)
{
  UVpage *uvp; /* The UV model paging descriptor */
  int ut;      /* The index of the current integration */
  int base;    /* The index of the current baseline */
  Intrec *rec; /* Pointer into integration record array ob->rec */
/*
 * Check validity of the descriptor.
 */
  if(!ob_ready(ob, OB_SELECT, "getmodel"))
    return 1;
/*
 * Check that the requested IF is in range.
 */
  if(cif < 0 || cif >= ob->nif) {
    lprintf(stderr, "getmodel: IF %d is unavailable.\n", cif);
    return 1;
  };
/*
 * If the model is already in memory, don't do anything.
 * Note that no model paging file is created for single IF data-sets,
 * so their model visibilities are always in memory, but ob->state can
 * still get demoted from OB_GETIF.
 */
  if((ob_ready(ob, OB_GETIF, NULL) && cif==ob->stream.cif) || ob->nif==1)
    return 0;
/*
 * Get the descriptor of the uvmodel scratch file.
 */
  uvp = ob->uvp;
/*
 * No uvmodel.scr file?
 */
  if(uvp==NULL) {
    lprintf(stderr,
	"getmodel: There is no UV model paging file to retrieve IF %d from.\n",
	cif);
    return 1;
  };
/*
 * Check the paging file descriptor for previous I/O errors.
 */
  if(uvp_error(uvp, "getmodel"))
    return 1;
/*
 * If an IF is currently in memory then we are about to over-write its
 * model visibilities so mark the IF as invalid.
 */
  ob->state = OB_SELECT;
/*
 * Read each integration from the UV model paging file.
 */
  rec = ob->rec;
  for(ut=0; ut<ob->nrec; ut++,rec++) {
    Integration *integ = rec->integ;
    Visibility *vis = integ->vis;
    int nbase = integ->sub->nbase;
    Mvis *mvis = uvp->mvis;
/*
 * Read the model for the next integration.
 */
    if(uvp_read(uvp, ut, cif))
      return 1;
/*
 * Copy each of its visibilities into the corresponding integration
 * of the Observation structure.
 */
    for(base=0; base<nbase; base++,vis++,mvis++) {
      vis->modamp = mvis->amp;
      vis->modphs = mvis->phs;
    };
  };
/*
 * Model OK.
 */
  return 0;
}

/*.......................................................................
 * Write the UV model of the current IF from the Observation descriptor
 * to the model paging file. If there is no model paging file
 * (ob->uvp==NULL) then this function does nothing, but returns as
 * though successful. 
 *
 * Input:
 *  ob    Observation *   The observation to write from.
 *  cif           int     The index of the IF to be written.
 * Output:
 *  return        int     0 - OK.
 *                        1 - Error.
 */
int putmodel(Observation *ob, int cif)
{
  Intrec *rec; /* Pointer into integration record array ob->rec */
  UVpage *uvp;       /* The UV model paging descriptor */
  int ut;            /* The index of the current integration */
  int base;          /* The index of the baseline being copied */
/*
 * Check validity of the descriptor.
 */
  if(!ob_ready(ob, OB_INDEX, "putmodel"))
    return 1;
/*
 * Check that the requested IF is in range.
 */
  if(cif < 0 || cif >= ob->nif) {
    lprintf(stderr, "putmodel: IF %d does not exist.\n", cif);
    return 1;
  };
/*
 * No model file is created for single IF data-sets, so do nothing
 * in such cases.
 */
  if(ob->nif==1)
    return 0;
/*
 * Get the UV paging descriptor.
 */
  uvp = ob->uvp;
  if(uvp==NULL) {
    lprintf(stderr, "putmodel: There's no model scratch file to write to.\n");
    return 1;
  };
/*
 * Check the IF paging descriptor for I/O errors.
 */
  if(uvp_error(uvp, "putmodel"))
    return 1;
/*
 * Write the model of each integration of the Observation structure to
 * the UV model paging file.
 */
  rec = ob->rec;
  for(ut=0; ut<ob->nrec; ut++,rec++) {
    Integration *integ = rec->integ;
    Visibility *vis = integ->vis;
    Mvis *mvis = uvp->mvis;
    int nbase = integ->sub->nbase;
/*
 * Copy visibilities from the current integration in the Observation
 * structure to the output buffer.
 */
    for(base=0; base<nbase; base++,vis++,mvis++) {
      mvis->amp = vis->modamp;
      mvis->phs = vis->modphs;
    };
/*
 * Write the integration to the IF paging file.
 */
    if(uvp_write(uvp, ut, cif))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Select a new UV data stream, compose it from the uvdata scratch file
 * and write the IF scratch file. Note that the UV representation of the
 * established model is always cleared, but the model components may
 * optionally be preserved in the tentative model.
 *
 * Input:
 *  ob     Observation * The descriptor of the observation.
 *  keep           int   If true retain the model of the current selection
 *                       for use in the new selection. Otherwise record
 *                       the current model for use when the current selection
 *                       is next requested, and reinstall the model that
 *                       applied to the new selection the last time that
 *                       it was selected (this will be an empty model if
 *                       the new selection hasn't been requested before).
 *                       In either case, the model will be entirely contained
 *                       in the tentative model, since the visibilities of
 *                       the established model will have to be recomputed
 *                       for the new selection.
 *  cl          Chlist * To select a new list of channel ranges send a
 *                       container returned by new_chlist(). Alternatively,
 *                       If this container contains 0 ranges, or cl==NULL
 *                       the current channels will be retained. The ranges
 *                       will be correctly truncated if the channel ranges
 *                       request channels beyond the highest channel. Do not
 *                       attempt to use cl after calling this function. It
 *                       will be deleted internally.
 *  stokes      Stokes   The required stokes parameter, or NO_POL for
 *                       the default (I or the first polarization).
 * Output:
 *  return         int   0 - OK (ob->state = OB_SELECT).
 *                       1 - Error (ob->state will reflect the state
 *                           that ob has been left in. Note that this
 *                           will still be OB_SELECT if the previous
 *                           selection is still valid).
 */
int ob_select(Observation *ob, int keep, Chlist *cl, Stokes stokes)
{
  Intrec *rec;   /* Pointer into integration record array ob->rec */
  Obpol *obpol;  /* Stream polarization descriptor */
  int ut;        /* The index of the integration being processed */
  int cif;       /* The index of the IF bieng processed */
  int base;      /* The index of the baseline being processed */
  int cr;        /* The channel range being processed */
  int chan;      /* The index of the spectral-line channel being processed */
/*
 * Check validity of arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "ob_select"))
    return 1;
/*
 * Ensure that all deferred edits have been applied.
 */
  if(ed_flush(ob))
    return 1;
/*
 * Preserve the current established model in the tentative model.
 */
  mergemod(ob, 0);
/*
 * Unless we are keeping the same model for all selections,
 * make a record of the current model for use the next time that
 * the current selection is requested.
 */
  if(!keep && ob_ready(ob, OB_SELECT, NULL) && ob_record_select_model(ob)) {
    cl = del_Chlist(cl);
    return 1;
  };
/*
 * Install the new channel list and polarization selections if valid.
 */
  if(ob_get_select(ob, cl, stokes)) {
    cl = del_Chlist(cl);
    return 1;
  };
/*
 * Extract pointers to the new selection descriptions.
 */
  cl = ob->stream.cl;
  obpol = &ob->stream.pol;
/*
 * Mark all per-baseline sums of weights as out of date.
 */
  flag_baseline_weights(ob, -1);
/*
 * If we are keeping different models for different channel/polarization
 * selections, restore any model that was previously made for this
 * selection. If no previously made model exists, clear the current model,
 * so as to start afresh.
 */
  if(!keep && ob_install_select_model(ob))
    return 1;
/*
 * Report the stream selection being used.
 */
  lprintf(stdout, "Selecting polarization: %s,  channels:",
	  Stokes_name(ob->stream.pol.type));
  for(cr=0; cr < cl->nrange; cr++)
    lprintf(stdout, " %d..%d", cl->range[cr].ca+1, cl->range[cr].cb+1);
  lprintf(stdout, "\n");
/*
 * For each IF separately read the uvdata.scr file, compose the new stream
 * in ob and then write it to the IF scratch file.
 */
  for(cif=0; cif<ob->nif; cif++) {
    Dif *dif = &ob->dp->ifs[cif];
    If *ifp = ob->ifs + cif;
    Chlist *if_cl = ifp->cl;
/*
 * Select the next IF in the output ifdata.scr file.
 */
    if(ip_range(ob->ip, cif, 0, ob->nbmax-1))
      return 1;
/*
 * Start the report of the number of channels being read from each IF.
 */
    lprintf(stdout, "Reading IF %d channels:", cif+1);
/*
 * Are any channels of this IF sampled?
 */
    if(if_cl) {
/*
 * List the channel ranges that will be read.
 */
      for(cr=0; cr < if_cl->nrange; cr++) {
	lprintf(stdout, " %d..%d",
		ifp->coff + if_cl->range[cr].ca+1,
		ifp->coff + if_cl->range[cr].cb+1);
      };
      lprintf(stdout, "\n");
/*
 * Set required uvdata-file paging ranges.
 */
      if(dp_crange(ob->dp, if_cl->ca, if_cl->cb) ||
	 dp_srange(ob->dp, 0, ob->npol-1)  ||
	 dp_brange(ob->dp, 0, ob->nbmax-1) ||
	 dp_irange(ob->dp, cif, cif))
	return 1;
/*
 * Read an integration at a time from the uvdata file and combine the
 * required spectral line channels and stokes parameters separately
 * for each IF. Place the results in ob.
 */
      rec = ob->rec;
      for(ut=0; ut<ob->nrec; ut++,rec++) {
	Integration *integ = rec->integ;
	Visibility *vis = integ->vis;
	int nbase = integ->sub->nbase;
/*
 * Read the next integration of the current IF from the uvdata scratch file
 * and construct the required polarization in the I/O buffer.
 */
	if(dp_read(ob->dp, ut))
	  return 1;
/*
 * Accumulate one visibility for each baseline.
 */
	for(base=0; base<nbase; base++,vis++) {
	  int deleted=0;  /* True if the target visibility should be deleted */
	  int flagged=0;  /* True if the target visibility should be flagged */
	  int npts=0;     /* The number of visibilities in the mean */
/*
 * Accumulate the weighted complex sum of visibilities over the
 * required spectral-line channels.
 */
	  Cvis sumvis={0.0f,0.0f,0.0f};
	  for(cr=0; cr < if_cl->nrange && !deleted; cr++) {
	    for(chan=if_cl->range[cr].ca; chan<=if_cl->range[cr].cb; chan++) {
	      Cvis curvis;
	      obpol->getpol(obpol, dif->chan[chan].base[base].pol, &curvis);
/*
 * Ascertain the usability of the new visibility.
 */
	      if(curvis.wt == 0.0f) {
		deleted = 1;
		break;
	      } else {
		if(curvis.wt < 0.0f) {
		  flagged = 1;
		  curvis.wt = -curvis.wt;
		};
/*
 * Accumulate the unweighted sum of spectral line channel visibilities.
 */
		npts++;
		sumvis.re += curvis.re;
		sumvis.im += curvis.im;
		sumvis.wt += 1.0f/curvis.wt; /* Accumulate variance sum */
	      };
	    };
	  };
/*
 * Convert the visibility sum into a mean.
 */
	  if(deleted || sumvis.wt==0.0f || npts==0) {
	    sumvis.re = sumvis.im = sumvis.wt = 0.0f; /* Deleted */
	  } else {
	    sumvis.re /= npts;
	    sumvis.im /= npts;
	    sumvis.wt = npts * npts / sumvis.wt;
	  };
/*
 * Record the visibility on the latest baseline in the respective
 * visibility in the current integration of ob. Note that if there
 * is only one IF in the observation, there will be no IF scratch
 * file.
 */
	  if(deleted || (sumvis.re==0.0f && sumvis.im==0.0f)) {
	    vis->amp = vis->phs = vis->wt = 0.0f;
	    vis->bad = FLAG_DEL;
	  } else {
	    vis->amp = sqrt(sumvis.re * sumvis.re + sumvis.im * sumvis.im);
	    vis->phs = atan2(sumvis.im, sumvis.re);
	    vis->wt = sumvis.wt;
	    vis->bad = flagged ? FLAG_BAD : 0;
	  };
	};
      };
/*
 * Having cached a whole IF in the observation structure, copy it to
 * the IF scratch file if there is one.
 */
      if(ob->ip) {
	rec = ob->rec;
	for(ut=0; ut<ob->nrec; ut++,rec++) {
	  Integration *integ = rec->integ;
	  Visibility *vis = integ->vis;
	  int nbase = integ->sub->nbase;
	  Dvis *dvis = ob->ip->dvis;
/*
 * Copy visibilities from the current integration in the Observation
 * structure to the output buffer.
 */
	  for(base=0; base<nbase; base++,vis++,dvis++) {
	    dvis->amp = vis->amp;
	    dvis->phs = vis->phs;
	    if(!(vis->bad & (FLAG_DEL | FLAG_BAD))) {
	      dvis->wt = vis->wt;   /* +ve weight means a good visibility */
	    } else if(vis->bad & FLAG_BAD) {
	      dvis->wt = -vis->wt;  /* -ve weight means a flagged visibility */
	    } else {
	      dvis->wt = 0.0f;      /* Zero weight means a deleted visibility */
	    };
	  };
/*
 * Write the latest integration to the output IF scratch file.
 */
	  if(ip_write(ob->ip, ut))
	    return 1;
	};
      };
/*
 * If there are no channels selected in the current IF, simply clear
 * output buffer and write all integrations.
 */
    } else {
/*
 * Report that no channels are being read from this IF.
 */
      lprintf(stdout, " (none)\n");
/*
 * Zero fill the IF scratch file records for the empty IF.
 */
      if(ob->ip) {
	ip_clear(ob->ip);
	for(ut=0; ut<ob->nrec; ut++,rec++) {
	  if(ip_write(ob->ip, ut))
	    return 1;
	};
      };
    };
  };
/*
 * Record the new observation state.
 */
  ob->state = OB_SELECT;
/*
 * If only one IF exists its visibilities are currently in memory, so
 * apply corrections to them and place the observation in OB_GETIF state.
 */
  if(ob->nif==1 && iniIF(ob, 0)) {
    ob->state = OB_SELECT;
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Call this function to check whether an observation desciptor is in
 * an appropraite state. This includes checking to see if the descriptor
 * is NULL, then checking whether the descriptor status recorded in
 * ob->state is at least as high as the state of readiness requested.
 * An error message is generated if the descriptor does not meet the
 * requested criteria, and 1 is returned.
 *
 * Input:
 *  ob    Observation *   The descriptor to be checked.
 *  state     Obstate     The required minimum descriptor state.
 *                        See obs.h for the list of valid values.
 *  name         char *   The name of the calling function, for use
 *                        in error message reporting.
 *                        Send NULL If you want to know the status
 *                        without an error message being reported.
 * Output:
 *  return        int     1 - Descriptor is ready for use.
 *                        0 - Descriptor is not ready at the requested level.
 */
int ob_ready(Observation *ob, Obstate state, const char *name)
{
  char *message=NULL; /* The error message to display */
/*
 * Check the Observation structure.
 */
  if(ob==NULL) {
    message = "Intercepted NULL Observation";
  } else if((int) state > (int) ob->state) {
    switch(ob->state) {
    default:
    case OB_BAD:
      message = "Observation corrupt";
      break;
    case OB_ALLOC:
      message = "No data read yet";
      break;
    case OB_DATA:
      message = "Integrations have not yet been indexed";
      break;
    case OB_INDEX:
      message = "No data stream selected yet";
      break;
    case OB_SELECT:
      message = "No IF in memory";
      break;
    };
  };
/*
 * Was an error detected?
 */
  if(message) {
/*
 * Only display the error message if a function name was given.
 */
    if(name)
      lprintf(stderr, "%s: %s.\n", name, message);
    return 0;
  };
/*
 * The observation appears to be OK.
 */
  return 1;
}

/*.......................................................................
 * Return the approriate UV scale factor for the given IF.
 *
 *  ob    Observation *   The descriptor of the observation.
 *  cif           int     The index of the IF to parameterise.
 * Output:
 *  return      float     The scale factor to use for the specified
 *                        IF. To get UVW coords in wavelengths multiply
 *                        by this factor. On error, 0.0f is returned.
 */
float getuvscale(Observation *ob, int cif)
{
/*
 * Check validity of the descriptor.
 */
  if(!ob_ready(ob, OB_INDEX, "getuvscale"))
    return 0.0f;
/*
 * Check that the requested IF is in range.
 */
  if(cif < 0 || cif >= ob->nif) {
    lprintf(stderr, "getuvscale: IF %d does not exist.\n", cif+1);
    return 0.0f;
  };
/*
 * The scale factor is the frequency of cited IF.
 */
  return (float) getfreq(ob, cif);
}

/*.......................................................................
 * Return the mean frequency pertaining to the currently selected
 * range of channels in one or all IFs.
 *
 * Input:
 *  ob   Observation *   The descriptor of the observation.
 *  cif          int     The index of the IF, or -1 for all IFs.
 * Output:
 *  return    double     The mean frequency of the selected channels of
 *                       the requested IFs.
 */
double getfreq(Observation *ob, int cif)
{
  double w_f_sum;/* Channel-width weighted sum of channel center frequencies */
  double w_sum;  /* Sum of channel widths (ie. overall bandwidth) */
  int bif,eif;   /* Start and end IF indexes */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_SELECT, "getfreq"))
    return 0.0;
/*
 * Get the indexes of the first and last IFs to be processed.
 */
  if(cif == -1) {
    bif = 0;
    eif = ob->nif - 1;
  } else if(cif>=0 && cif < ob->nif) {
    bif = eif = cif;
  } else {
    lprintf(stderr, "getfreq: IF index out of range.\n");
    return 0.0;
  };
/*
 * Sum the frequencies of the selected channels of each IF.
 */
  w_f_sum = 0.0f;
  w_sum = 0.0f;
  for(cif=bif; cif<=eif; cif++) {
/*
 * Get the descriptor of the IF.
 */
    If *ifp = &ob->ifs[cif];
/*
 * Get the list of selected channels.
 */
    Chlist *cl = ifp->cl;
    if(cl) {
      int sc=0;      /* Twice sum of selected channel indexes */
      int nc=0;      /* The number of channels in sc */
      Chans *cr;     /* The channel range being processed */
      for(cr=cl->range; cr < cl->range + cl->nrange; cr++) {
	int ca = cr->ca;
	int cb = cr->cb;
	int n = cb - ca + 1;
	nc += n;
	sc += n * (ca + cb);
      };
      if(nc) {
/*
 * Accumulate the sum of channel center frequencies, weighted by the
 * channel bandwidth. Note that sc = 2 * sum_i(c[i]).
 */
	w_f_sum += fabs(ifp->df) * (nc * ifp->freq + 0.5 * sc * ifp->df);
/*
 * Accumulate the total bandwidth covered by all selected channels.
 */
	w_sum += nc * fabs(ifp->df);
      };
    };
  };
  return (w_sum > 0.0) ? (w_f_sum / w_sum) : (ob->ifs[(bif+eif)/2].freq);
}

/*.......................................................................
 * Return the total bandwidth pertaining to the currently selected
 * range of channels in one or all IFs.
 *
 * Input:
 *  ob   Observation *   The descriptor of the observation.
 *  cif          int     The index of the IF, or -1 for all IFs.
 * Output:
 *  return    double     The bandwidth covered by the selected channel
 *                       range in the given IF. On error 0.0 is returned.
 */
double getbw(Observation *ob, int cif)
{
  int bif,eif;         /* Start and end IF indexes */
  double bw_sum = 0.0; /* Sum of channel bandwidths */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_SELECT, "getbw"))
    return 0.0;
/*
 * Get the indexes of the first and last IFs to be processed.
 */
  if(cif == -1) {
    bif = 0;
    eif = ob->nif - 1;
  } else if(cif>=0 && cif < ob->nif) {
    bif = eif = cif;
  } else {
    lprintf(stderr, "getbw: IF index out of range.\n");
    return 0.0;
  };
/*
 * Sum the channel bandwidths of all selected channels of IFs
 * bif to eif.
 */
  for(cif=bif; cif<=eif; cif++) {
/*
 * Get the descriptor of the IF.
 */
    If *ifp = &ob->ifs[cif];
/*
 * Get the list of selected channels.
 */
    Chlist *cl = ifp->cl;
    if(cl) {
      int nc=0;      /* The number of channels in all channel ranges of cl */
      Chans *cr;     /* The channel range being processed */
      for(cr=cl->range; cr < cl->range + cl->nrange; cr++)
	nc += cr->cb - cr->ca + 1;
      if(nc)
	bw_sum += nc * fabs(ifp->df);
    };
  };
/*
 * Return the bandwidth.
 */
  return bw_sum;
}

/*.......................................................................
 * Find the ob->rec[] index of the nearest integration who's
 * time-stamp matches a given relational test with respect to a given
 * time-stamp.
 *
 * Input:
 *  ob   Observation *   The observation to be searched.
 *  ut        double     The time-stamp to compare against.
 *                       (seconds from the start of year ob->year).
 *  op        UTfind     The relational test, from:
 *                        UT_LT   -   time < ut.
 *                        UT_LE   -   time <= ut.
 *                        UT_NR   -   time ~= ut (nearest UT).
 *                        UT_GE   -   time >= ut.
 *                        UT_GT   -   time
 * Output:
 *  return       int     The index of a matching integration record, or
 *                       -1 if not found.
 */
int ob_find_ut(Observation *ob, double ut, UTfind op)
{
  int lo;         /* Index of integration with time < ut */
  int best;       /* Index of integration with time nearest ut */
  int hi;         /* Index of integration with time > ut */
  int slot;       /* The returned slot, or -1 if not found */
  double lo_ut;   /* Time-stamp corresponding to 'lo' */
  double best_ut; /* Time-stamp corresponding to 'best' */
  double hi_ut;   /* Time-stamp corresponding to 'hi' */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "ob_find_ut"))
    return -1;
/*
 * Perform a binary search for the closest match to the given time-stamp.
 */
  lo = 0;
  hi = ob->nrec-1;
  while(lo <= hi) {
    int mid = (lo+hi)/2;
    double utmid = ob->rec[mid].integ->ut;
    if(ut < utmid)
      hi = mid - 1;
    else
      lo = mid + 1;
  };
/*
 * Currently hi > lo. Swap them, to ease comprehension.
 */
  {int itmp=lo; lo=hi; hi=itmp;}
/*
 * Limit the upper and lower search indexes.
 */
  if(lo < 0)
    lo = 0;
  if(hi >= ob->nrec)
    hi = ob->nrec - 1;
/*
 * Get the time-stamps associated with lo and hi.
 */
  lo_ut = ob->rec[lo].integ->ut;
  hi_ut = ob->rec[hi].integ->ut;
/*
 * Arrange for lo<best<hi.
 */
  if(ut - lo_ut < hi_ut - ut) {
    best = lo;
    best_ut = lo_ut;
  } else {
    best = hi;
    best_ut = hi_ut;
  };
  if(lo==best && lo>0) {
    lo--;
    lo_ut = ob->rec[lo].integ->ut;
  };
  if(hi==best && hi<ob->nrec-1) {
    hi++;
    hi_ut = ob->rec[hi].integ->ut;
  };
/*
 * Find the output slot according to the given relational operator.
 */
  slot = -1;
  switch(op) {
  case UT_LT:
    if(best_ut < ut)
      slot = best;
    else if(lo_ut < ut)
      slot = lo;
    break;
  case UT_LE:
    if(best_ut <= ut)
      slot = best;
    else if(lo_ut <= ut)
      slot = lo;
    break;
  case UT_NR:
    slot = best;
    break;
  case UT_GE:
    if(best_ut >= ut)
      slot = best;
    else if(hi_ut >= ut)
      slot = hi;
    break;
  case UT_GT:
    if(best_ut > ut)
      slot = best;
    else if(hi_ut > ut)
      slot = hi;
    break;
  default:
    lprintf(stderr, "ob_find_ut: Unrecognized relational operator.\n");
    slot = -1;
  };
  return slot;
}

/*.......................................................................
 * Find the sub->integ[] index of the nearest integration who's
 * time-stamp matches a given relational test with respect to a given
 * time-stamp.
 *
 * Input:
 *  sub     Subarray *   The sub-array to be searched.
 *  ut        double     The time-stamp to compare against.
 *                       (seconds from the start of year ob->year).
 *  op        UTfind     The relational test, from:
 *                        UT_LT   -   time < ut.
 *                        UT_LE   -   time <= ut.
 *                        UT_NR   -   time ~= ut (nearest UT).
 *                        UT_GE   -   time >= ut.
 *                        UT_GT   -   time
 * Output:
 *  return       int     The index of a matching integration record, or
 *                       -1 if not found.
 */
int sub_find_ut(Subarray *sub, double ut, UTfind op)
{
  int lo;         /* Index of integration with time < ut */
  int best;       /* Index of integration with time nearest ut */
  int hi;         /* Index of integration with time > ut */
  int slot;       /* The returned slot, or -1 if not found */
  double lo_ut;   /* Time-stamp corresponding to 'lo' */
  double best_ut; /* Time-stamp corresponding to 'best' */
  double hi_ut;   /* Time-stamp corresponding to 'hi' */
/*
 * Check arguments.
 */
  if(sub_bad(sub, "sub_find_ut"))
    return -1;
/*
 * Perform a binary search for the closest match to the given time-stamp.
 */
  lo = 0;
  hi = sub->ntime-1;
  while(lo < hi) {
    int mid = (lo+hi)/2;
    double utmid = sub->integ[mid].ut;
    if(ut < utmid)
      hi = mid - 1;
    else
      lo = mid + 1;
  };
/*
 * Currently hi > lo. Swap them, to ease comprehension.
 */
  {int itmp=lo; lo=hi; hi=itmp;}
/*
 * Limit the upper and lower search indexes.
 */
  if(lo < 0)
    lo = 0;
  if(hi >= sub->ntime)
    hi = sub->ntime - 1;
/*
 * Get the time-stamps associated with lo and hi.
 */
  lo_ut = sub->integ[lo].ut;
  hi_ut = sub->integ[hi].ut;
/*
 * Arrange for lo<best<hi.
 */
  if(ut - lo_ut < hi_ut - ut) {
    best = lo;
    best_ut = lo_ut;
  } else {
    best = hi;
    best_ut = hi_ut;
  };
  if(lo==best && lo>0) {
    lo--;
    lo_ut = sub->integ[lo].ut;
  };
  if(hi==best && hi<sub->ntime-1) {
    hi++;
    hi_ut = sub->integ[hi].ut;
  };
/*
 * Find the output slot according to the given relational operator.
 */
  slot = -1;
  switch(op) {
  case UT_LT:
    if(best_ut < ut)
      slot = best;
    else if(lo_ut < ut)
      slot = lo;
    break;
  case UT_LE:
    if(best_ut <= ut)
      slot = best;
    else if(lo_ut <= ut)
      slot = lo;
    break;
  case UT_NR:
    slot = best;
    break;
  case UT_GE:
    if(best_ut >= ut)
      slot = best;
    else if(hi_ut >= ut)
      slot = hi;
    break;
  case UT_GT:
    if(best_ut > ut)
      slot = best;
    else if(hi_ut > ut)
      slot = hi;
    break;
  default:
    lprintf(stderr, "sub_find_ut: Unrecognized relational operator.\n");
    slot = -1;
  };
  return slot;
}

/*.......................................................................
 * Return the Right-Acsension corresponding to easterly and northerly
 * SIN-projection direction cosine offsets from a given reference RA
 * and Dec. See AIPS Memo 27 for details.
 *
 * Input:
 *  ra     double    The reference R.A., (radians).
 *  dec    double    The reference Dec, (radians).
 *  l      double    Easterly projected direction cosine offset -1...1.
 *  m      double    Northerly projected direction cosine offset -1...1.
 *  proj     Proj    The spherical projection of l,m.
 * Output:
 *  return double    The offset Right-Ascension (radians).
 */
double lmtora(double ra, double dec, double l, double m, Proj proj)
{
  double newra = 0.0;
  double rtmp;
/*
 * Direction cosines run between -1 and 1.
 */
  if(l <= 1.0 && l >= -1.0 && m <= 1.0 && m >= -1.0) {
    switch(proj) {
    case PRJ_SIN:
      rtmp = cos(dec)*sqrt(fabs(1.0-l*l-m*m)) - m*sin(dec);
      if(rtmp!=0.0)
	newra = ra + atan2(l, rtmp);
      break;
    case PRJ_NCP:
      rtmp = cos(dec) - m * sin(dec);
      if(rtmp!=0.0)
	newra = ra + atan2(l, rtmp);
      break;
    default:
      lprintf(stderr, "lmtora: Unrecognized projection (%s).", Proj_name(proj));
      break;
    };
  };
  return newra;
}

/*.......................................................................
 * Return the Declination corresponding to easterly and northerly
 * SIN-projection direction cosine offsets from a given reference RA
 * and Dec. See AIPS Memo 27 for details.
 *
 * Input:
 *  ra     double    The reference R.A., (radians).
 *  dec    double    The reference Dec, (radians).
 *  l      double    Easterly SIN-proj direction cosine offset -1..1.
 *  m      double    Northerly SIN-proj direction cosine offset -1..1.
 *  proj     Proj    The spherical projection of l,m.
 * Output:
 *  return double    The offset Right-Ascension (radians).
 */
double lmtodec(double ra, double dec, double l, double m, Proj proj)
{
  double newdec = 0.0;  /* The new declination to be returned */
  double dtmp;          /* Work variable */
/*
 * Direction cosines run between -1 and 1.
 */
  if(l <= 1.0 && l >= -1.0 && m <= 1.0 && m >= -1.0) {
/*
 * Convert the l,m projection offsets into a new declination.
 */
    switch(proj) {
    case PRJ_SIN:
      dtmp = m * cos(dec) + sin(dec) * sqrt(fabs(1.0 - l*l - m*m));
      if(fabs(dtmp) <= 1.0)
	newdec = asin(dtmp);
      break;
    case PRJ_NCP:
      dtmp = cos(dec) - m * sin(dec);
      if(dtmp != 0.0) {
	double dtmp1 = cos(atan2(l,dtmp));
	if(dtmp1 != 0.0) {
	  dtmp1 = dtmp / dtmp1;
	  if(fabs(dtmp1) <= 1.0)
	    newdec = (dec < 0.0 ? -1.0 : 1.0) * acos(dtmp1);
	};
      };
      break;
    default:
      lprintf(stderr, "lmtodec: Unrecognized projection (%s).",Proj_name(proj));
      break;
    };
  };
  return newdec;
}

/*.......................................................................
 * Given the RA,DEC of the pointing center of an observation
 * and the RA,DEC of a neighboring point, return the equivalent easterly
 * SIN-projection direction-cosine offset from the pointing center.
 * See AIPS Memo 27 for details.
 *
 * Input:
 *  ref_ra   double    The reference R.A., (radians).
 *  ref_dec  double    The reference Dec, (radians).
 *  ra, dec  double    The target RA,Dec (radians).
 *  proj       Proj    The spherical projection of the returned offset.
 * Output:
 *  return   double    Easterly projected direction cosine offset.
 */
double radec_to_l(double ref_ra, double ref_dec, double ra, double dec,
		  Proj proj)
{
/*
 * Compute the projected coordinate according to the type of projection
 * requested.
 */
  switch(proj) {
  case PRJ_SIN:
  case PRJ_NCP:
    return cos(dec) * sin(ra - ref_ra);
    break;
  default:
    lprintf(stderr, "radec_to_l: Unrecognized projection (%s).",
	    Proj_name(proj));
    return 0;
  };
}

/*.......................................................................
 * Given the RA,DEC of the pointing center of an observation
 * and the RA,DEC of a neighboring point, return the equivalent northerly
 * SIN-projection direction-cosine offset from the pointing center.
 * See AIPS Memo 27 for details.
 *
 * Input:
 *  ref_ra   double    The reference R.A., (radians).
 *  ref_dec  double    The reference Dec, (radians).
 *  ra, dec  double    The target RA,Dec (radians).
 *  proj       Proj    The spherical projection of the returned offset.
 * Output:
 *  return   double    Northerly projected direction cosine offset.
 */
double radec_to_m(double ref_ra, double ref_dec, double ra, double dec,
		  Proj proj)
{
  double tmp;
/*
 * Perform the computation in accordance with the requested projection.
 */
  switch(proj) {
  case PRJ_SIN:
    return sin(dec)*cos(ref_dec) - cos(dec) * sin(ref_dec) * cos(ra - ref_ra);
    break;
  case PRJ_NCP:
    tmp = sin(ref_dec);
    if(tmp==0.0) {
      lprintf(stderr, "radec_to_m: NCP projection isn't defined at dec=0.0.\n");
      return 0.0;
    };
    return (cos(ref_dec) - cos(dec)*cos(ra - ref_ra))/tmp;
    break;
  default:
    lprintf(stderr, "radec_to_m: Unrecognized projection (%s).",
	    Proj_name(proj));
    return 0.0;
  };
}

/*.......................................................................
 * Private function of ob_select() used to check and install new channel
 * and polarization selection descriptors. On success, ob->state is
 * demoted to OB_INDEX to reflect the fact that the selection has not
 * been established by ob_select().
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 *  cl       Chlist *  The list of channels to be installed.
 *                     Channels in the list should be in the domain
 *                     0..(ob->nif * ob->nchan - 1).
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error creating the new lists.
 *                         The observation is left as it was on entry.
 */    
static int ob_get_select(Observation *ob, Chlist *cl, Stokes pol)
{
  Obpol obpol;          /* The polarization descriptor */
  Chlist *def_cl=NULL;  /* Default channel list if allocated */
  Chlist **if_cl=NULL;  /* Array of ob->nif Chlist's */
  int cif;              /* The index of the IF being processed */
/*
 * Check the validity of the arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "ob_chlist"))
    return 1;
/*
 * See if the required polarization was observed or can be derived, and
 * get its descriptor.
 */
  if(get_Obpol(ob, pol, 1, &obpol) &&
     (pol!=NO_POL && get_Obpol(ob, NO_POL, 1, &obpol)))
    return 1;
/*
 * Do we need to install a new channel list?
 */
  if(!ob->stream.cl || (cl && cl != ob->stream.cl)) {
/*
 * If no channel list was given, substitute a default.
 */
    if(!cl) {
/*
 * Create a default list covering all channels in the observation.
 */
      def_cl = new_Chlist();
      if(!def_cl || add_crange(def_cl, 0, ob->nctotal-1))
	return bad_ob_chlist(ob, def_cl, if_cl);
/*
 * Substitute the default channel list.
 * Note that def_cl must still be deleted on error.
 */
      cl = def_cl;
    };
/*
 * Limit the new channel list to the range of channels available
 * in the observation and require at least one channel range to be left.
 */
    {
      int nleft = lim_Chlist(cl, ob->nctotal);
      if(nleft < 0)
	return bad_ob_chlist(ob, def_cl, if_cl);
      if(nleft < 1) {
	lprintf(stderr, "ob_chlist: No channels selected.\n");
	return bad_ob_chlist(ob, def_cl, if_cl);
      };
    };
/*
 * Allocate an array of ob->nif Chlist pointers so that we can construct
 * the new channel lists without touching the current lists.
 * If there is an error we can then delete the array without touching ob.
 */
    if_cl = (Chlist **) malloc(ob->nif * sizeof(Chlist *));
    if(!if_cl) {
      lprintf(stderr, "ob_chlist: Insufficient memory.\n");
      return bad_ob_chlist(ob, def_cl, if_cl);
    };
/*
 * Initialize the array so that we know which elements have been allocated.
 */
    for(cif=0; cif<ob->nif; cif++)
      if_cl[cif] = NULL;
/*
 * Construct separate channel lists for each IF.
 */
    for(cif=0; cif<ob->nif; cif++) {
      Chlist *clp = sub_Chlist(cl, ob->ifs[cif].coff, ob->nchan);
      if(!clp)
	return bad_ob_chlist(ob, def_cl, if_cl);
      else if(clp->nrange < 1)
	clp = del_Chlist(clp);
      else
	if_cl[cif] = clp;
    };
  };
/*
 * Record the fact that the new selection has not been established yet.
 */
  ob->state = OB_INDEX;
/*
 * Now that we know that both the polarization and channel lists are
 * ok, install them.
 */
  if(if_cl) {
    ob->stream.cl = del_Chlist(ob->stream.cl);
    ob->stream.cl = cl;
    for(cif=0; cif<ob->nif; cif++) {
      If *ifp = ob->ifs + cif;
      ifp->cl = del_Chlist(ifp->cl);
      ifp->cl = if_cl[cif];
    };
/*
 * Discard the temporary array of IF channel lists.
 */
    free(if_cl);
  };
/*
 * Store the new polarization selection.
 */
  ob->stream.pol = obpol;
  return 0;
}

/*.......................................................................
 * Private error cleanup and return function of ob_chlist().
 */
static int bad_ob_chlist(Observation *ob, Chlist *def_cl, Chlist **if_cl)
{
  int cif;
  def_cl = del_Chlist(def_cl);
  if(if_cl) {
    for(cif=0; cif < ob->nif; cif++)
      if_cl[cif] = del_Chlist(if_cl[cif]);
  };
  return 1;  /* Error return value of ob_chlist() */
}

/*.......................................................................
 * Return the index of the current stream IF, or -1 if not available.
 * This is meant to be paired with set_cif_state() to bracket code that
 * potentially changes the current IF.
 *
 * Input:
 *  ob    Observation *  The descriptor of the observation.
 * Output:
 *  return        int    The index of the current IF, or -1 if none
 *                       is currently in memory.
 */
int get_cif_state(Observation *ob)
{
  return ob_ready(ob, OB_GETIF, NULL) ? ob->stream.cif : -1;
}

/*.......................................................................
 * Take a value previously returned by get_cif_state() and use it to
 * restore the IF of that time if one was available then.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation.
 *  cif          int    A previous return value of get_cif_state().
 *                                 -1 - Do nothing.
 *                       0..ob->nif-1 - The index of the IF to be
 *                                      restored into memory.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int set_cif_state(Observation *ob, int cif)
{
  return cif == -1 ? 0 : getIF(ob, cif);
}

/*.......................................................................
 * Replace the current model with any model saved for the current
 * channel-range and polarization.
 *
 * Input:
 *  ob   Observation *  The observation in which to install the model.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int ob_install_select_model(Observation *ob)
{
  Model *newmod;   /* The new model to be installed */
/*
 * Check the arguments.
 */
  if(!ob) {
    lprintf(stderr, "ob_install_select_model: NULL argument(s).\n");
    return 1;
  };
/*
 * We need a channel-range and polarization to lookup the model for.
 * Note that we can't use ob_ready(OB_SELECT) to test for this, because
 * this function is called from ob_select() itself before it has
 * changed the observation state.
 */
  if(!ob->stream.cl || ob->stream.pol.type==NO_POL) {
    lprintf(stderr, "ob_install_select_model: No stream has been selected.\n");
    return 1;
  };
/*
 * Clear the now redundant tentative and established model.
 */
  if(clrmod(ob, 1, 1, 0))
    return 1;
/*
 * See if there is an entry for the given channel-range/polarization
 * selection.
 */
  newmod = rem_ModelEntry(ob->mtab, ob->stream.cl, ob->stream.pol.type,
			  ob->geom.east, ob->geom.north);
/*
 * If there is a model for this selection, replace the now empty
 * tentative model with the new model. Otherwise leave the tentative
 * model empty.
 */
  if(newmod) {
    ob->newmod = del_Model(ob->newmod);
    ob->newmod = newmod;
    if(newmod->ncmp > 0)
      lprintf(stdout, "Restored previously made model of latest selection.\n");
  };
  return 0;
}

/*.......................................................................
 * Record the current model in the model table, indexed by the
 * currently selected channel-range and polarization. Note that this
 * function has the side effect of unestablishing the established model.
 *
 * Input:
 *  ob     Observation *  The observation containing the model.
 * Output:
 *  return         int    0 - OK.
 *                        1 - Error.
 */
int ob_record_select_model(Observation *ob)
{
/*
 * We need a stream to save the model with.
 */
  if(!ob_ready(ob, OB_SELECT, "ob_record_select_model"))
    return 1;
/*
 * Merge the current established model into the tentative model.
 */
  mergemod(ob, 0);
/*
 * Record the tentative model in the model table, indexed by the
 * currently selected channel-range and polarization.
 */
  if(!add_ModelEntry(ob->mtab, ob->newmod, ob->stream.cl, ob->stream.pol.type,
		     ob->geom.east, ob->geom.north))
    return 1;
  return 0;
}

/*.......................................................................
 * Compute the sum of per-baseline weights of one or all IFs that are
 * flagged as being out of date.
 *
 * Input:
 *  ob    Observation *  The container of the observed data.
 *  cif            int   The index of the target IF, or -1 for all IFs.
 * Output:
 *  return         int   0 - OK.
 *                       1 - Error.
 */
int update_baseline_weights(Observation *ob, int cif)
{
  int bif,eif;      /* Start and end IF indexes */
  int old_if;       /* The selected IF index on entry to this function */
  int isub;         /* The index of the sub-array being traversed */
  int base;         /* The index of a baseline of the current subarray */
  int ut;           /* The index of the integration being processed */
/*
 * Do nothing if no observed data has been selected yet.
 */
  if(!ob_ready(ob, OB_SELECT, NULL))
    return 0;
/*
 * Get the indexes of the first and last IFs to be processed.
 */
  if(cif == -1) {
    bif = 0;
    eif = ob->nif - 1;
  } else if(cif >= 0 && cif < ob->nif) {
    bif = eif = cif;
  } else {
    lprintf(stderr, "update_baseline_weights: IF index out of range.\n");
    return 1;
  };
/*
 * Store the index of the currently selected IF, so that we can restore
 * it before returning.
 */
  old_if = get_cif_state(ob);
/*
 * Process the selected IFs.
 */
  for(cif=bif; (cif=nextIF(ob, cif, 1, 1)) >= 0 && cif<=eif; cif++) {
    If *ifp = ob->ifs + cif;
/*
 * Do the baseline sums need to be updated in this IF?
 */
    if(ifp->wtsum_bad) {
/*
 * Get the new IF.
 */
      if(getIF(ob, cif))
	return 1;
/*
 * Zero the per-baseline sums of the current IF.
 */
      for(isub=0; isub<ob->nsub; isub++) {
	Subarray *sub = ob->sub + isub;
	for(base=0; base<sub->nbase; base++)
	  sub->base[base].bwt[cif].wtsum = 0.0;
      };
/*
 * Add up per-baseline weights of unflagged visibilities of the current IF.
 */
      for(isub=0; isub<ob->nsub; isub++) {
	Subarray *sub = ob->sub + isub;
	for(ut=0; ut<sub->ntime; ut++) {
	  Integration *integ = sub->integ + ut;
	  for(base=0; base<sub->nbase; base++) {
	    Baseline *bptr = sub->base + base;
	    Visibility *vis = integ->vis + base;
	    if(!vis->bad)
	      bptr->bwt[cif].wtsum += vis->wt;
	  };
	};
      };
/*
 * Mark the sum as valid.
 */
      ifp->wtsum_bad = 0;
    };
  };
/*
 * Reinstate the IF of the caller.
 */
  if(set_cif_state(ob, old_if))
    return 1;
  return 0;
}

/*.......................................................................
 * This function should be called when the weights of any visibilities
 * of a given IF (or all IFs) are changed. It marks the per-baseline
 * sums of weights as needing recomputation for the specified IFs.
 *
 * Input:
 *  ob    Observation *  The container of the observed data.
 *  cif            int   The index of the target IF, or -1 for all IFs.
 * Output:
 *  return         int   0 - OK.
 *                       1 - Error.
 */
int flag_baseline_weights(Observation *ob, int cif)
{
  int bif,eif;      /* Start and end IF indexes */
/*
 * Do nothing if no observed data has been selected yet.
 */
  if(!ob_ready(ob, OB_SELECT, NULL))
    return 0;
/*
 * Get the indexes of the first and last IFs to be flagged.
 */
  if(cif == -1) {
    bif = 0;
    eif = ob->nif - 1;
  } else if(cif >= 0 && cif < ob->nif) {
    bif = eif = cif;
  } else {
    lprintf(stderr, "flag_baseline_weights: IF index out of range.\n");
    return 1;
  };
/*
 * Apply the flags.
 */
  for(cif=bif; cif<=eif; cif++)
    ob->ifs[cif].wtsum_bad = 1;
  return 0;
}

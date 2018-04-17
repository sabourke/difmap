#include <stdio.h>
#include <stdlib.h>

#include "modvis.h"
#include "obs.h"
#include "vlbutil.h"
#include "mapwin.h"
#include "winmod.h"
#include "obwin.h"
#include "logio.h"

/*
 * This module is dedicated to the components and UV representations
 * of models recorded in an Observation structure. All but the
 * UV representation can be handled with the observation in OB_INDEX
 * state. Changes to the UV representation of the established models
 * are handled through fixmod(). This function quietly ignores
 * requests when there is no selection to compute visibilities for.
 * The only function that will touch the model visibilities when not
 * in OB_SELECT state is clrmod(), which when asked to clear
 * the established model, will clear the visibilities in the UV model
 * scratch file, preparatory to making a new selection.
 *
 * Note that ob_select() moves all established model components into
 * the un-established Model containers. This ensures that if model
 * components get placed into the established model containers when no
 * selection is in effect, they will be retreived when the next
 * selection is made.
 *
 * Also see obshift.c.
 */

static int fixmod(Observation *ob, Model *mod, Modcmp *cmp, int doadd);
static int uvaddmod(Observation *ob, Model *mod);
static int uvsubmod(Observation *ob, Model *mod);

static int obvret(Model *tmpmod, int iret);

/*.......................................................................
 * Add a model to the established or tentative models of an Observation.
 *
 * Input:
 *  ob  Observation *  The observation to which the model is to be
 *                     added.
 *  mod       Model *  The model to be added.
 *                     On output the contents of this model will
 *                     have been transferred to ob->model, and the
 *                     container will be empty and should be deleted
 *                     if no longer required.
 *  keep        int    0 - Add to the tentative model.
 *                     1 - Add to the established model.
 *  docont      int    0 - Add to the normal model.
 *                     1 - Add to the continuum model.
 *  append      int    0 - Prepend components of 'mod' to the specified model.
 *                     1 - Append components of 'mod' to the specified model.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int obaddmod(Observation *ob, Model *mod, int keep, int docont, int append)
{
/*
 * Observation ok?
 */
  if(!ob_ready(ob, OB_INDEX, "obaddmod"))
    return 1;
/*
 * Anything to add?
 */
  if(mod && mod->ncmp>0) {
/*
 * Add to the established model?
 */
    if(keep) {
/*
 * Add the UV representation of the model.
 */
      if(uvaddmod(ob, mod))
	return 1;
/*
 * Add to the appropriate established model component list.
 */
      if(docont) {
	ob->cmodel = add_mod(ob->cmodel, mod, 1, append);
	lprintf(stdout,
		"The UV continuum model now contains %d components and %g Jy\n",
		ob->cmodel->ncmp, ob->cmodel->flux);
      } else {
	ob->model  = add_mod(ob->model,  mod, 1, append);
	lprintf(stdout,
		"The established model now contains %d components and %g Jy\n",
	      ob->model->ncmp, ob->model->flux);
      };
    }
/*
 * Add to the appropriate tentative model.
 */
    else {
      if(docont)
	ob->cnewmod = add_mod(ob->cnewmod, mod, 1, append);
      else
	ob->newmod  = add_mod(ob->newmod,  mod, 1, append);
    };
  };
  return 0;
}

/*.......................................................................
 * Clear the UV representation of the established model but preserve the
 * components of the established model in the tentative model.
 *
 * Input:
 *  ob  Observation *  The observation to which the model is to be
 *                     added.
 *  doold       int    If true establish the tentative model.
 *                     If false relegate the established model back
 *                     to the tentative model.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int mergemod(Observation *ob, int doold)
{
/*
 * Observation ok?
 */
  if(!ob_ready(ob, OB_INDEX, "mergemod"))
    return 1;
/*
 * Establish the tentative models?
 */
  if(doold) {
    if(obaddmod(ob, ob->newmod, 1, 0, 1) ||
       obaddmod(ob, ob->cnewmod, 1, 1, 1))
      return 1;
  } else {
/*
 * Prepend the model components of the established models to the tentative
 * models.
 */
    if(obaddmod(ob, ob->model, 0, 0, 0) ||
       obaddmod(ob, ob->cmodel, 0, 1, 0))
      return 1;
/*
 * Clear the UV representation of the now empty established model.
 */
    clrmod(ob, 1, 0, 0);
  };
  return 0;
}

/*.......................................................................
 * Combine the continuum and normal models to a combined continuum or
 * normal model.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation containing the
 *                      models to be combined.
 *  tocont       int    0 - Prepend the current continuum models to the
 *                          normal models.
 *                      1 - Append the current normal models to the
 *                          continuum models.
 * Output:
 *  return       int    0 - OK.
 */
int setcmod(Observation *ob, int tocont)
{
/*
 * Observation ok?
 */
  if(!ob_ready(ob, OB_INDEX, "setcmod"))
    return 1;
  if(tocont) {
    ob->cnewmod = add_mod(ob->cnewmod, ob->newmod, 1, 1);
    ob->cmodel  = add_mod(ob->cmodel,  ob->model, 1, 1);
  } else {
    ob->newmod = add_mod(ob->newmod, ob->cnewmod, 1, 0);
    ob->model  = add_mod(ob->model,  ob->cmodel, 1, 0);
  };
  return 0;
}

/*.......................................................................
 * Delete the components of the established and tentative models that lie
 * optionally either inside or outside a given list of map windows.
 *
 * Input:
 *  ob   Observation *  The observation whose established model is to be
 *                      modified.
 *  mw        Mapwin *  The list of map windows.
 *  doout        int    If 0, keep only components that lie within the
 *                      list of windows. Otherwise keep only those outside
 *                      the windows.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int obwinmod(Observation *ob, Mapwin *mw, int doout)
{
  Model *imod;    /* List of components lying within windows */
  Model *omod;    /* List of components lying outside all windows */
  int ncmp;       /* The original number of components in the model */
/*
 * Observation ok?
 */
  if(!ob_ready(ob, OB_INDEX, "obwinmod"))
    return 1;
/*
 * No clean windows to apply?
 */
  if(mw==NULL || mw->nwin<1)
    return 0;
/*
 * Window the established model.
 */
  if(ob->model->ncmp>0) {
/*
 * Record the existing size of the established model.
 */
    ncmp = ob->model->ncmp;
/*
 * Split the model into two lists of components, one containing the
 * components that lie outside the windows, and another containing the rest.
 */
    imod = win_mod(ob->model, mw, 1);
    if(imod==NULL)
      return 1;
    omod = ob->model;
/*
 * Delete the UV representation of the obsolete part of the established
 * model.
 */
    if(doout) {
      uvsubmod(ob, imod);
      imod = del_Model(imod);
      ob->model = omod;
    } else {
      uvsubmod(ob, omod);
      omod = del_Model(omod);
      ob->model = imod;
    };
/*
 * Report the state of the established model.
 */
    lprintf(stdout,
 "The established clean model now contains %d of the original %d components.\n",
	    ob->model->ncmp, ncmp);
  };
/*
 * Window the tentative clean model, keeping the required parts of the
 * model in 'ob->newmod' and discarding the unwanted parts.
 */
  if(ob->newmod->ncmp > 0) {
    ncmp = ob->newmod->ncmp;
/*
 * Split the model into two lists of components, one containing the
 * components that lie outside the windows, and another containing the rest.
 */
    imod = win_mod(ob->newmod, mw, 1);
    if(imod==NULL)
      return 1;
    omod = ob->newmod;
/*
 * Install the required part of the model and delete the the other part.
 */
    if(doout) {
      imod = del_Model(imod);
      ob->newmod = omod;
    } else {
      omod = del_Model(omod);
      ob->newmod = imod;
    };
    lprintf(stdout,
	 "The tentative model now contains %d of the original %d components.\n",
	    ob->newmod->ncmp, ncmp);
  };
/*
 * Report the combined flux in the latest and established clean models.
 */
  lprintf(stdout,
	  "Remaining flux in the tentative and established models = %g Jy\n",
	  ob->newmod->flux + ob->model->flux);  
  return 0;
}

/*.......................................................................
 * Add a model to the model visibilities of an observation.
 * 'mod' may be sent as NULL in which case the function will return
 * without doing anything.
 *
 * Input:
 *  ob  Observation *  The observation to which the model is to be
 *                     added.
 *  mod       Model *  The model to be added.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int uvaddmod(Observation *ob, Model *mod)
{
  return fixmod(ob, mod, NULL, 1);
}

/*.......................................................................
 * Subtract a model from the model visibilities of an observation.
 * 'mod' may be sent as NULL in which case the function will return
 * without doing anything.
 *
 * Input:
 *  ob  Observation *  The observation from which the model is to be
 *                     subtracted.
 *  mod       Model *  The model to be subtracted.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int uvsubmod(Observation *ob, Model *mod)
{
  return fixmod(ob, mod, NULL, 0);
}

/*.......................................................................
 * Add or subtract a model from the UV models in an Observation.
 *
 * Input:
 *  ob  Observation *  The observation to which the model is to be
 *                     added.
 *  mod       Model *  The model to be added.
 *  cmp      Modcmp *  If mod==NULL, component cmp will be used instead.
 *  doadd       int    If true add the model.
 *                     If false subtract the model.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int fixmod(Observation *ob, Model *mod, Modcmp *cmp, int doadd)
{
  float uvscale; /* The UVW coordinate scale factor. */
  int cif;       /* The index of the IF being processed */
  int isub;      /* The index of the sub-array being processed */
  int ut;        /* The index of the integration being processed */
  int base;      /* The index of baseline visibility being processed */
  int domod;     /* If true use mod argument, otherwise use cmp argument */
  int old_if;    /* State of current IF to be restored on exit */
/*
 * Quietly ignore this call if there is no selection to compute model
 * visibilities for.
 */
  if(!ob_ready(ob, OB_SELECT, NULL))
    return 0;
/*
 * Model or component?
 */
  domod = (mod != NULL);
/*
 * Nothing to add?
 */
  if((domod && mod->ncmp < 1) || (!domod && cmp==NULL))
    return 0;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Inform user.
 */
  if(doadd) {
    lprintf(stdout, "Adding %d model components to the UV plane model.\n",
	    domod ? mod->ncmp : 1);
  } else {
    lprintf(stdout, "Extracting %d model components from the UV plane model.\n",
	    domod ? mod->ncmp : 1);
  };
/*
 * Record the fact that model visibilities now exist in 'ob'.
 */
  ob->hasmod=1;
/*
 * Fix all sampled IFs.
 * Note that the model visibilities are particular to this stream and
 * will be discarded when another stream is selected, so there is no
 * need to calculate visibilities for unsampled IFs.  
 */
  for(cif=0; (cif=nextIF(ob,cif,1,1)) >= 0; cif++) {
    Subarray *sub = ob->sub;
/*
 * Get the frequency of the IF.
 */
    float freq = getfreq(ob, cif);
/*
 * Get the model of the next IF.
 */
    if(getIF(ob, cif))
      return 1;
/*
 * Get the factor required to convert the UVW coordinates
 * from light seconds to wavelengths.
 */
    uvscale = ob->stream.uvscale;
/*
 * Loop through all sub-arrays of the new IF.
 */
    for(isub=0; isub<ob->nsub; isub++, sub++) {
      Integration *integ = sub->integ;
/*
 * Loop through the integrations of the current sub-array.
 */
      for(ut=0; ut<sub->ntime; ut++,integ++) {
	Visibility *vis = integ->vis;
/*
 * Fix the model of each visibility in the integration.
 */
	for(base=0; base<sub->nbase; base++,vis++) {
	  float u = vis->u * uvscale;
	  float v = vis->v * uvscale;
/*
 * Only calculate the model for good and flagged visibilities of the
 * current stream. This is both an optimization and a way of avoiding
 * the garbage UVW coordinates that often accompany deleted visibilities.
 * Note that it _is_ necessary to calculate the model for flagged visibilities
 * in case the user unflags them at a later time.
 */
	  if(!(vis->bad & FLAG_DEL)) {
/*
 * Accumulate the real and imaginary parts of the model visibility.
 */
	    float re=0.0, im=0.0;
	    if(domod) {
	      for(cmp=mod->head; cmp; cmp=cmp->next)
		add_cmp_to_modvis(cmp, sub, base, freq, u, v, &re, &im);
	    } else {
	      add_cmp_to_modvis(cmp, sub, base, freq, u, v, &re, &im);
	    };
/*
 * If subtracting the model, negate the visibility.
 */
	    if(!doadd) {
	      re = -re;
	      im = -im;
	    };
/*
 * Add the new model visibilty to the existing one.
 */
	    add_cart_to_polar(&vis->modamp, &vis->modphs, re, im);
	  };
	};
      };
    };
/*
 * Store the modified model in the uvmodel.scr scratch file.
 */
    if(putmodel(ob, cif))
      return 1;
  };
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    return 1;
/*
 * Now update the zero-spacing model amplitude.
 * This is simply the sum of all model-component fluxes.
 */
  {
    float modamp = 0.0f;
    if(domod) {
      for(cmp=mod->head; cmp; cmp = cmp->next)
	modamp += cmp->flux;
    } else {
      modamp = cmp->flux;
    };
    ob->uvzero.modamp += doadd ? modamp : -modamp;
  };
  return 0;
}

/*.......................................................................
 * Re-arrange the established model and tentative models such
 * that the fixed components of both models become the new established
 * model, while the remaining components - those with free parameters -
 * become the new tentative model, for use in model-fitting related
 * applications. 
 *
 * Input:
 *  ob  Observation *  The observation from whos established model
 *                     the variable components are to be removed.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int obvarmod(Observation *ob)
{
  Model *tmpmod;       /* A temporary model */
/*
 * Observation ok?
 */
  if(!ob_ready(ob, OB_INDEX, "obvarmod"))
    return 1;
/*
 * Keep the user informed.
 */
  lprintf(stdout, "Partitioning the model into established and variable parts.\n");
/*
 * Create a temporary model.
 */
  tmpmod = new_Model();
  if(tmpmod==NULL)
    return 1;
/*
 * Establish the continuum model, since we won't be fitting this.
 */
  if(obaddmod(ob, ob->cnewmod, 1, 1, 1))
    return 1;
/*
 * Place the variable components of the tentative model in tmpmod.
 */
  if(ob->newmod->ncmp > 0) {
    var_mod(ob->newmod, tmpmod);
/*
 * Place the remaining fixed components of the tentative model in the
 * established model.
 */
    if(ob->newmod->ncmp > 0 && obaddmod(ob, ob->newmod, 1, 0, 1))
      return obvret(tmpmod, 1);
  };
/*
 * Place the variable components of the established model in the now emptied
 * tentative model.
 */
  if(ob->model->ncmp > 0) {
    var_mod(ob->model, ob->newmod);
/*
 * Subtract the UV representation of the extracted components.
 */
    if(ob->newmod->ncmp > 0 && uvsubmod(ob, ob->newmod))
      return obvret(tmpmod, 1);
  };
/*
 * Append the variable components of the original tentative model
 * back into the tentative model (Note that doing things in this order
 * ensures that components from the established model always appear first,
 * thus preserving the original time ordering).
 */
  if(tmpmod->ncmp > 0 && obaddmod(ob, tmpmod, 0, 0, 1))
    return obvret(tmpmod, 1);
/*
 * Report the stats of the fixed and variable parts of the model.
 */
  lprintf(stdout,
	  "The fixed established model contains %d components (%g Jy).\n",
	  ob->model->ncmp, ob->model->flux);
  lprintf(stdout,
	  "The variable part of the model contains %d components (%g Jy).\n",
	  ob->newmod->ncmp, ob->newmod->flux);
  return obvret(tmpmod, 0);
}

/*.......................................................................
 * Private cleanup/return function of obvarmod().
 */
static int obvret(Model *tmpmod, int iret)
{
  tmpmod = del_Model(tmpmod);
  return iret;
}

/*.......................................................................
 * Allow the user to edit the variable part or all of the established and
 * tentative models in an external editor. The edited model will be the
 * new tentative model ob->newmod. 
 *
 * Input:
 *  ob  Observation *  The observation from whos established model
 *                     the variable components are to be removed and
 *                     edited.
 *  dovar       int    If true edit just the variable components.
 *                     If false, edit all components.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int obedmod(Observation *ob, int dovar)
{
/*
 * Check that we have been given a valid observation.
 */
  if(!ob_ready(ob, OB_INDEX, "obedmod"))
    return 1;
/*
 * Rearrange the tentative and established model components such that
 * those components that are to be edited are placed in the tentative
 * model while the rest are placed in the established model.
 */
  if(dovar ? obvarmod(ob) : mergemod(ob, 0))
    return 1;
/*
 * Edit the variable part of the model.
 */
  ob->newmod = ed_model(ob->newmod);
  return 0;
}

/*.......................................................................
 * Add a model component to the established or tentative model of an
 * observation. If the model is to be added to the established model then
 * this includes computing its UV representation.
 *
 * Note that if you are adding more than one component then it is much
 * more efficient to collect the components into a model and use obaddmod().
 *
 * Input:
 *  ob   Observation * The descriptor of the observation.
 *  cmp       Modcmp * The component to be added.
 *  keep        int    0 - Add to the tentative model.
 *                     1 - Add to the established model.
 * Output:
 *  return    Modcmp * The destination model component. Note that because
 *                     of the possibility of component merging this may
 *                     not be cmp, in which case cmp will have been deleted
 *                     and should not be used further. If the specified
 *                     model is not in a modifiable state, NULL is
 *                     returned, (cmp is not deleted).
 */
Modcmp *obaddcmp(Observation *ob, Modcmp *cmp, int keep)
{
/*
 * Signal an error if the observation is not in an appropriate state to
 * accept modifications to the specified model.
 */
  if(!ob_ready(ob, OB_INDEX, "obaddcmp"))
    return NULL;
/*
 * Add to the established model?
 */
  if(keep) {
/*
 * Compute the UV representation of the component (Note that this must
 * be done before calling add_cmp() since add_cmp() may merge the component
 * with another component).
 */
    if(fixmod(ob, NULL, cmp, 1))
      return NULL;
/*
 * Add the component to the established model component list.
 */
    cmp = add_cmp(cmp, ob->model, 1);
  }
/*
 * Add to the tentative model.
 */
  else {
    cmp = add_cmp(cmp, ob->newmod, 1);
  };
  return cmp;
}

/*.......................................................................
 * Remove a model component from the established or tentative model of an
 * observation. If the model is to be removed from the established model
 * then this includes computing its UV representation.
 *
 * Note that doing this one component at a time is inefficient.
 *
 * Input:
 *  ob   Observation * The descriptor of the observation.
 *  cmp       Modcmp * The component to be added.
 *  keep        int    0 - Remove from the tentative model.
 *                     1 - Remove from the established model.
 * Output:
 *  return    Modcmp * The removed model component. This should be deleted
 *                     if no longer required. If the component is not found
 *                     or the specified model is not in a modifiable state,
 *                     NULL is returned.
 */
Modcmp *obremcmp(Observation *ob, Modcmp *cmp, int keep)
{
  Modcmp *next;  /* The next component to be checked */
  Modcmp *prev;  /* The prev component to be checked */
/*
 * Signal an error if the observation is not in an appropriate state to
 * accept modifications to the specified model.
 */
  if(!ob_ready(ob, OB_INDEX, "obremcmp"))
    return NULL;
/*
 * Remove from the established model?
 */
  if(keep) {
/*
 * Locate the component to be removed.
 */
    prev = NULL;
    for(next=ob->model->head; next && next!=cmp; prev=next,next=prev->next);
/*
 * Was the component found?
 */
    if(next && next==cmp) {
/*
 * Compute the UV representation of the component.
 */
      if(fixmod(ob, NULL, cmp, 0))
	return NULL;
/*
 * Remove the component from the established model component list.
 */
      cmp = rem_cmp(ob->model, prev, next);
    } else {
      lprintf(stderr, "obremcmp: Component not found.\n");
      return NULL;
    };
  }
/*
 * Remove from the tentative model.
 */
  else {
/*
 * Locate the component to be removed.
 */
    prev = NULL;
    for(next=ob->newmod->head; next && next!=cmp; prev=next,next=prev->next);
/*
 * Was the component found?
 */
    if(next && next==cmp) {
/*
 * Remove the component from the established model component list.
 */
      cmp = rem_cmp(ob->newmod, prev, next);
    } else {
      lprintf(stderr, "obremcmp: Component not found.\n");
      return NULL;
    };
  };
  return cmp;
}

/*.......................................................................
 * Clear the established and/or tentative and/or continuum models of an
 * observation.
 * If the established model is cleared, its UV representation in memory and
 * in the model paging file will be cleared and the hasmod flag will be
 * reset to 0. The model visibilities will always be cleared when doold
 * is true, regardless of the value of the 'hasmod' flag. This ensures that
 * there is a method to clear the UV model, regardless of state.
 *
 * Input:
 *  ob  Observation *  The observation whose model visibilities are to
 *                     be cleared.
 *  doold       int    If true, clear the established model and its
 *                     UV representation.
 *  donew       int    If true, clear the tentative model.
 *  docont      int    If true, clear the continuum models.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int clrmod(Observation *ob, int doold, int donew, int docont)
{
  int wasclr=0;     /* True if the pertinent models were already clear */
  int iret=0;       /* Return code - set to 1 on error */
/*
 * No observation to delete model from?
 */
  if(!ob_ready(ob, OB_INDEX, "clrmod"))
    return iret;
/*
 * Are the pertinent models already clear?
 */
  wasclr = !((doold && ob->model->ncmp > 0) ||
	     (donew && ob->newmod->ncmp > 0) ||
	     (docont && (ob->cnewmod->ncmp + ob->cmodel->ncmp) > 0));
/*
 * Always clear model visibilities if doold is true. Also clear model
 * visibilities if clearing the established continuum model when no
 * normal established model exists - this is quicker than using subtraction.
 */
  if(doold || (docont && ob->model->ncmp<1)) {
    int isub;         /* The index of the sub-array being processed */
    int ut;           /* The integration ut number */
    int base;         /* The baseline number */
/*
 * Clear the model visibilities of each Sub-array.
 */
    for(isub=0; isub<ob->nsub; isub++) {
      Subarray *sub = &ob->sub[isub];
      for(ut=0; ut<sub->ntime; ut++) {
	Visibility *vis = sub->integ[ut].vis;
	for(base=0; base<sub->nbase; base++,vis++)
	  vis->modamp = vis->modphs = 0.0f;
      };
    };
/*
 * If there is a UV model paging file, have the model cleared there as
 * well.
 */
    if(ob->uvp) {
      int cif;                      /* The index of the IF being corrected */
/*
 * Copy the cleared model to each IF in the UV model paging file.
 */
      for(cif=0; cif<ob->nif; cif++) {
	if(putmodel(ob, cif))
	  iret = 1;
      };
    };
/*
 * Clear the zero-baseline flux model.
 */
    ob->uvzero.modamp = 0.0f;
/*
 * Record the fact that there are now no model visibilities.
 */
    ob->hasmod = 0;
  }
/*
 * If we need to delete the established continuum model while there is
 * an established model, clear it by subtracting the continuum model
 * from the model visibilities.
 */
  else if(docont) {
    uvsubmod(ob, ob->cmodel);
  };
/*
 * Delete the contents of the continuum models?
 */
  if(docont) {
    ob->cmodel = clr_Model(ob->cmodel);
    ob->cnewmod = clr_Model(ob->cnewmod);
  };
/*
 * Delete the components of the established model?
 */
  if(doold) {
    ob->model = clr_Model(ob->model);
  };
/*
 * Clear the tentative model?
 */
  if(donew) {
    ob->newmod = clr_Model(ob->newmod);
  };
/*
 * If the model visibilities were cleared to delete the established model,
 * but the continuum model is not being deleted, prepend the established part
 * of the continuum model to ob->cnewmod.
 */
  if(!ob->hasmod && ob->cmodel->ncmp > 0)
    ob->cnewmod = add_mod(ob->cnewmod, ob->cmodel, 1, 0);
/*
 * Inform the user only if components were actually cleared.
 * This allows the function to be called to clear the model visibilities, or
 * to be called by multiple functions without producing redundant messages.
 */
  if(!wasclr) {
    int nclr = (doold!=0) + (donew!=0) + (docont!=0);
/*
 * Compose a single output line naming the models that have been cleared.
 */
    if(nclr > 0) {
      int nmore = nclr;
      lprintf(stdout, "clrmod: Cleared the");
      if(doold) {
	nmore--;
	lprintf(stdout," established%s", nmore==0 ? "":(nmore==1 ? " and":","));
      };
      if(donew) {
	nmore--;
	lprintf(stdout, " tentative%s", nmore==0 ? "":" and");
      };
      if(docont)
	lprintf(stdout, " continuum");
      lprintf(stdout, " model%s.\n", nclr>1 ? "s":"");
    };
  };
  return iret;
}


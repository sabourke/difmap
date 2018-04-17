#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logio.h"
#include "vlbconst.h"
#include "obs.h"

/*.......................................................................
 * Determine and apply a single correction amp+phase residual offset
 * to a given baseline in the current IF. The amplitude correction is
 * determined from the ratio of the weighted means of the observed and
 * established model amplitudes (after establishing the tentative model).
 * The phase correction is determined from the weighted mean phase offset
 * between the observed and model phases.
 *
 * Input/Output:
 *  ob  Observation *  The observation descriptor.
 *                     The resulting corrections are applied to the data
 *                     and stored in the baseline descriptors in *ob.
 * doall        int    If true correct all IFs. If 0, correct only the
 *                     current IF.
 *  base        int    The index of the baseline to be corrected.
 *  isub        int    The index of the sub-array to be corrected.
 * Output:
 *  return      int    0 - OK.
 */
int resoff(Observation *ob, int doall, int base, int isub)
{
  Visibility *vis; /* A single visibility */
  Subarray *sub;   /* The sub-array being corrected */
  float amp;       /* The weighted sum of observed amplitudes on a baseline */
  float modamp;    /* The weighted sum of model amplitudes on a baseline */
  float phsdif;    /* The difference between model and observed phase */
  float phsoff;    /* The weighted sum of phase differences. */
  float wt,wtsum;  /* The weight of a visibility and sum of weights */
  float ampcor;    /* The amplitude correction for a baseline */
  float phscor;    /* The amplitude correction for a baseline */
  int ut;          /* Indexes of baselines and integrations */
  int bif, eif;    /* The indexes of the first and last IF to be processed */
  int cif;         /* The index of the IF being processed */
  int old_if;      /* State of current IF to be restored on exit */
/*
 * Check inputs.
 */
  if(!ob_ready(ob, doall ? OB_SELECT : OB_GETIF, "resoff"))
    return 1;
  if(isub < 0 || isub >= ob->nsub) {
    lprintf(stderr, "resoff: Out of range sub-array index received.\n");
    return 1;
  };
/*
 * Get the subarray descriptor.
 */
  sub = &ob->sub[isub];
/*
 * Check requested baseline index.
 */
  if(base < 0 || base >= sub->nbase) {
    lprintf(stderr, "resoff: Out-of range baseline index received.\n");
    return 1;
  };
/*
 * Determine the range of IFs to be processed.
 */
  if(doall) {
    bif = 0;
    eif = ob->nif - 1;
  } else {
    bif = eif = ob->stream.cif;
  };  
/*
 * Establish the tentative model.
 */
  if(mergemod(ob, 1))
    return 1;
/*
 * Abort if there is no established model.
 */
  if(ob->model->ncmp + ob->cmodel->ncmp < 1) {
    lprintf(stderr, "resoff: There is no model to reference to.\n");
    return 1;
  };
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Process all sampled IFs over the range bif to eif.
 */
  for(cif=bif; (cif=nextIF(ob, cif, 1, 1)) >= 0 && cif<=eif; cif++) {
/*
 * Get the visibilities of the IF to be processed next.
 */
    if(getIF(ob, cif))
      return 1;
/*
 * Mark the per-baseline sums of weights of the IF as stale.
 */
    flag_baseline_weights(ob, cif);
/*
 * Determine the means for the given baseline.
 */
    amp = modamp = phsoff = wtsum = 0.0f;
    for(ut=0; ut<sub->ntime; ut++) {
      Integration *integ = &sub->integ[ut];
      vis = &integ->vis[base];
      if(!vis->bad) {
/*
 * Use the visibility weight to weight the means.
 */
	wt = vis->wt;
	wtsum += wt;
/*
 * Weighted sums of observed and model amplitudes.
 */
	amp += wt * vis->amp;
	modamp += wt * vis->modamp;
/*
 * Determine the phase offset of this visibility, wrapped into
 * the range -pi to pi.
 */
	phsdif = vis->modphs - vis->phs;
	phsdif -= twopi*floor(phsdif/twopi + 0.5);
/*
 * Accumulate weighted sum of phase offsets.
 */
	phsoff += wt * phsdif;
      };
    };
/*
 * Where there any unflagged visibilities?
 */
    if(wtsum > 0.0f) {
/*
 * Determine the phase correction.
 */
      phscor = phsoff / wtsum;
/*
 * Determine the amplitude correction.
 */
      ampcor = (modamp > 0.0f && amp > 0.0f) ? (modamp / amp) : 1.0f;
/*
 * Apply the corrections to all but the deleted visibilities.
 */
      for(ut=0; ut<sub->ntime; ut++) {
	Integration *integ = &sub->integ[ut];
	Visibility *vis = &integ->vis[base];
/*
 * Don't correct deleted visibilities.
 */
	if(!(vis->bad & FLAG_DEL)) {
	  vis->amp *= ampcor;
	  vis->wt /= ampcor * ampcor;   /* NB. wt = 1/amp_err^2 */
	  vis->phs += phscor;
	};
      };
/*
 * Store the modified baseline corrections.
 */
      {
	Bascor *bcor = &sub->base[base].bcor[cif];
	bcor->amp_cor *= ampcor;
	bcor->phs_cor += phscor;
      };
/*
 * Report the corrections.
 */
      lprintf(stdout,
        "IF %d %d:%.10s-%.10s: %g amplitude change, %g degrees phase offset.\n",
        cif+1, isub+1, sub->tel[sub->base[base].tel_a].name,
        sub->tel[sub->base[base].tel_b].name,
        ampcor, rtod*phscor);
    };
  };
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    return 1;
  return 0;
}

/*.......................................................................
 * Undo all baseline based calibrations of one or all IFs.
 * If the observation state is currently OB_GETIF and the resident IF
 * is one of the IFs that is being uncorrected, then the visibilities
 * in memory will be modified accordingly.
 *
 * Input:
 *  ob    Observation *  The observation descriptor.
 *  doall         int    If true process all IFs. If 0, process only the
 *                       current IF.
 *  doamp         int    If true, undo the baseline amplitude calibration.
 *  dophs         int    If true, undo the baseline phase calibration.
 * Output:
 *  return        int    0 - OK.
 */
int clroff(Observation *ob, int doall, int doamp, int dophs)
{
  int base,ut;     /* Baseline and integration indexes */
  int isub;        /* The index of the subarray being processed */
  float ampcor;    /* The amplitude correction for a baseline */
  float phscor;    /* The phase correction for a baseline */
/*
 * Check inputs.
 */
  if(!ob_ready(ob, doall ? OB_INDEX : OB_GETIF, "clroff"))
    return 1;
/*
 * Propogate changes to the current IF if there is one and it is
 * sampled.
 */
  if(ob_ready(ob, OB_GETIF, NULL) && ob->ifs[ob->stream.cif].cl) {
/*
 * Undo the corrections of the current IF one sub-array at a time.
 */
    for(isub=0; isub<ob->nsub; isub++) {
      Subarray *sub = &ob->sub[isub];
      Integration *integ = sub->integ;
/*
 * Undo the corrections in each integration of the current IF and sub-array.
 */
      for(ut=0; ut<sub->ntime; ut++,integ++) {
	Baseline *bptr = sub->base;
	Visibility *vis = integ->vis;
/*
 * Undo each baseline correction.
 */
	for(base=0; base<sub->nbase; base++,bptr++,vis++) {
	  Bascor *bcor = &bptr->bcor[ob->stream.cif];
	  if(!(vis->bad & FLAG_DEL)) {
	    ampcor = bcor->amp_cor;
	    phscor = bcor->phs_cor;
	    if(doamp && ampcor > 0.0f) {
	      vis->amp /= ampcor;
	      vis->wt *= ampcor * ampcor;
	    };
	    if(dophs)
	      vis->phs -= phscor;
	  };
	};
      };
    };
  };
/*
 * Reset the corrections for the specified IFs.
 */
  if(ini_bcor(ob, doall, doamp, dophs))
    return 1;
  return 0;
}

/*.......................................................................
 * Apply all the stored resoff corrections to the visibilities. This should
 * only be done when a new uncalibrated IF-stream is paged into the
 * Observation structure, via iniIF().
 *
 * Input:
 *  ob     Observation * The descriptor of the observation.
 *  cif            int   The index of the IF whose corrections are to be
 *                       applied.
 * Output:
 *  return         int   0 - OK.
 *                       1 - Error.
 */
int app_bcor(Observation *ob, int cif)
{
  int isub; /* The index of the subarray being processed */
  int ut;   /* The index of integration *integ */
  int base; /* The index of the baseline being corrected */
/*
 * Make sure that we are applying corrections to uncorrected visibilities.
 * obutil.c::iniIF() sets ob->state=OB_RAWIF until all corrections have
 * been applied.
 */
  if(!ob || ob->state != OB_RAWIF) {
    lprintf(stderr, "app_Telcor: No uncorrected visibilities to correct.\n");
    return 1;
  };
  if(cif<0 || cif>ob->nif-1) {
    lprintf(stderr, "app_bcor: Out of bounds IF index intercepted.\n");
    return 1;
  };
/*
 * Apply the corrections of the current IF one sub-array at a time.
 */
  for(isub=0; isub<ob->nsub; isub++) {
    Subarray *sub = &ob->sub[isub];
    Integration *integ = sub->integ;
/*
 * Undo the corrections in each integration of the current IF and sub-array.
 */
    for(ut=0; ut<sub->ntime; ut++,integ++) {
      Baseline *bptr = sub->base;
      Visibility *vis = integ->vis;
/*
 * Calibrate each baseline correction.
 */
      for(base=0; base<sub->nbase; base++,bptr++,vis++) {
	Bascor *bcor = &bptr->bcor[cif];
	if(!(vis->bad & FLAG_DEL)) {
	  float ampcor = bcor->amp_cor;
	  float phscor = bcor->phs_cor;
	  if(ampcor > 0.0f) {
	    vis->amp *= ampcor;
	    vis->wt /= ampcor * ampcor;   /* NB. wt = 1/amp_err^2 */
	  };
	  vis->phs += phscor;
	};
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Initialize all resoff corrections in all IFs and sub-arrays.
 *
 * Input:
 *  ob     Observation * The descriptor of the observation.
 *  doall          int   If true process all IFs. If 0, process only the
 *                       current IF.
 *  doamp          int   If true, clear the baseline amplitude calibrations.
 *  dophs          int   If true, clear the baseline phase calibrations.
 * Output:
 *  return         int   0 - OK.
 *                       1 - Error.
 */
int ini_bcor(Observation *ob, int doall, int doamp, int dophs)
{
  Subarray *sub; /* The sub-array being processed */
  int isub;      /* The index of the subarray being processed */
  int base;      /* The index of the baseline being processed */
  int bif, eif;  /* The indexes of the first and last IFs to be processed */
  int cif;       /* The index of the IF being processed */
/*
 * Sanity check.
 */
  if(!ob_ready(ob, doall ? OB_INDEX : OB_GETIF, "ini_bcor"))
    return 1;
/*
 * Determine the range of IFs to be processed.
 */
  if(doall) {
    bif = 0;
    eif = ob->nif - 1;
  } else {
    bif = eif = ob->stream.cif;
  };  
/*
 * Initialize the corrections one sub-array at a time.
 */
  sub = ob->sub;
  for(isub=0; isub<ob->nsub; isub++,sub++) {
    Baseline *bptr = sub->base;
/*
 * Reset the corrections in the baseline descriptors of each sub-array.
 */
    for(base=0; base<sub->nbase; base++,bptr++) {
      Bascor *bcor = bptr->bcor;
/*
 * Initialize the corrections of all the specified IFs.
 */
      for(cif=bif; cif<=eif; cif++,bcor++) {
/*
 * Mark the per-baseline sums of weights of the IF as stale.
 */
	flag_baseline_weights(ob, cif);
/*
 * Reset the specified calibrations.
 */
	if(doamp)
	  bcor->amp_cor = 1.0f;
	if(dophs)
	  bcor->phs_cor = 0.0f;
      };
    };
  };
  return 0;
}

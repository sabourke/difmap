#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "obs.h"
#include "logio.h"

/*.......................................................................
 * Undo all recorded telescope amplitude and/or phase corrections
 * and/or correction flags for all IFs. If an IF is currently in
 * memory, its visibilities will be modified accordingly.
 *
 * Input:
 *  ob     Observation *   The observation to be uncorrected.
 *  doamp          int     If true undo amplitude corrections.
 *  dophs          int     If true undo phase corrections.
 *  doflag         int     If true undo correction flags.
 *  doreset        int     If true reset the recorded corrections to
 *                         null values. This should only be false if the
 *                         caller needs to re-apply modified versions
 *                         of the corrections using 'recalib()'.
 */
void uncalib(Observation *ob, int doamp, int dophs, int doflag, int doreset)
{
  Subarray *sub;    /* The sub-array being processed */
  int ut;           /* The index of the Integration being processed */
  int base;         /* The index of the baseline being processed */
  int isub;         /* The index of the sub-array being processed */
  int itel;         /* The index of a telescope */
/*
 * Are there corrections to be reset?
 */
  if(!ob_ready(ob, OB_INDEX, "uncalib"))
    return;
/*
 * Is there anything to be done?
 */
  if(doamp || dophs || doflag) {
/*
 * If changes to either the weights or the flagging of any visibilities
 * have been requested, mark the per-baseline sums of weights as being
 * stale.
 */
    if(doamp || doflag)
      flag_baseline_weights(ob, -1);
/*
 * If an IF is in memory, remove the corrections.
 */
    if(ob_ready(ob, OB_GETIF, NULL)) {
/*
 * Get the index of the current IF.
 */
      int cif = ob->stream.cif;
/*
 * Uncalibrate each sub-array of the current IF.
 */
      sub = ob->sub;
      for(isub=0; isub<ob->nsub; isub++,sub++) {
	Integration *integ = sub->integ;
/*
 * Loop over each integration of the sub-array to be corrected.
 */
	for(ut=0; ut<sub->ntime; ut++,integ++) {
	  Telcor *tcor = integ->icor[cif].tcor;
	  Visibility *vis = integ->vis;
	  Baseline *bptr = sub->base;
/*
 * Undo the corrections on each baseline at the current ut.
 */
	  for(base=0; base<sub->nbase; base++,vis++,bptr++) {
	    Telcor *ta_cor = &tcor[bptr->tel_a];
	    Telcor *tb_cor = &tcor[bptr->tel_b];
/*
 * Undo visibility phase correction.
 */
	    if(dophs)
	      vis->phs -= ta_cor->phs_cor - tb_cor->phs_cor;
/*
 * Undo visibility amplitude correction.
 */
	    if(doamp) {
	      float gcor = ta_cor->amp_cor * tb_cor->amp_cor;
	      if(gcor > 0.0f) {
		vis->amp /= gcor;
		vis->wt *= gcor * gcor;   /* NB. wt = 1/amp_err^2 */
	      };
	    };
/*
 * Undo visibility correction flags.
 */
	    if(doflag)
	      vis->bad &= ~(FLAG_TA | FLAG_TB);
	  };
	};
      };
    };
/*
 * Reset the corrections in all IFs of each sub-array.
 */
    if(doreset) {
/*
 * Loop for the corrections of each sub-array.
 */
      sub = ob->sub;
      for(isub=0; isub<ob->nsub; isub++,sub++) {
	Integration *integ = sub->integ;
/*
 * Loop over integrations in the sub-array.
 */
	for(ut=0; ut<sub->ntime; ut++,integ++) {
	  Intcor *icor = integ->icor;
	  int cif;
/*
 * Loop for the corrections in each IF in the sub-array.
 */
	  for(cif=0; cif<ob->nif; cif++,icor++) {
	    Telcor *tcor = icor->tcor;
/*
 * Reset the correction pertaining to each telescope of the integration.
 */
	    for(itel=0; itel<sub->nstat; itel++,tcor++) {
/*
 * Reset the phase correction.
 */
	      if(dophs)
		tcor->phs_cor = 0.0;
/*
 * Reset the amplitude correction.
 */
	      if(doamp)
		tcor->amp_cor = 1.0;
/*
 * Reset the correction flag.
 */
	      if(doflag)
		tcor->bad = 0;
	    };
	  };
	};
      };
    };
  };
  return;
}

/*.......................................................................
 * Apply recorded telescope amplitude and phase corrections along with
 * correction flags, to the IF in memory. This should only be called
 * by iniIF().
 *
 * Input:
 *  ob     Observation *   The observation to be corrected.
 *  cif            int     The index of the IF whose corrections are
 *                         to be applied. No check is made to ensure
 *                         that this IF is in memory, because this
 *                         function is called by getIF itself.
 */
int app_Telcor(Observation *ob, int cif)
{
  Subarray *sub;    /* The sub-array being processed */
  int ut;           /* The index of the Integration being processed */
  int isub;         /* The index of the sub-array being processed */
  int base;         /* The index of the baseline being processed */
/*
 * Make sure that we are applying corrections to uncorrected visibilities.
 * obutil.c::iniIF() sets ob->state=OB_RAWIF until all corrections have
 * been applied.
 */
  if(!ob || ob->state != OB_RAWIF) {
    lprintf(stderr, "app_Telcor: No uncorrected visibilities to correct.\n");
    return 1;
  };
/*
 * Re-calibrate each sub-array of the given IF.
 */
  sub = ob->sub;
  for(isub=0; isub<ob->nsub; isub++,sub++) {
    Integration *integ = sub->integ;
/*
 * Loop over each integration of the sub-array to be corrected.
 */
    for(ut=0; ut<sub->ntime; ut++,integ++) {
      Visibility *vis = integ->vis;
      Baseline *bptr = sub->base;
      Telcor *tcor = integ->icor[cif].tcor;
/*
 * Apply the corrections on each baseline at the current ut.
 */
      for(base=0; base<sub->nbase; base++,bptr++,vis++) {
	Telcor *ta_cor = &tcor[bptr->tel_a];
	Telcor *tb_cor = &tcor[bptr->tel_b];
/*
 * Apply phase correction.
 */
	vis->phs += ta_cor->phs_cor - tb_cor->phs_cor;
/*
 * Apply amplitude correction.
 */
	{
	  float gcor = ta_cor->amp_cor * tb_cor->amp_cor;
	  if(gcor > 0.0f) {
	    vis->amp *= gcor;
	    vis->wt /= gcor * gcor;   /* NB. wt = 1/amp_err^2 */
	  };
	};
/*
 * Apply telescope correction flags.
 */
	if(ta_cor->bad)
	  vis->bad |= FLAG_TA;
	if(tb_cor->bad)
	  vis->bad |= FLAG_TB;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Flag or unflag a telescope correction and propagate the flags to all
 * effected visibilities.
 *
 * Input:
 *  ob     Observation *   The observation to be flagged.
 *  sub       Subarray *   The subarray to be flagged.
 *  cif            int     The index of the IF of the integration to be
 *                         flagged. If this IF is in memory, its flags will
 *                         be modified.
 *  ut             int     The index of the integration to be flagged.
 *  itel           int     The index of the telescope to be flagged.
 *  doflag         int     0 - Remove flags.
 *                         1 - Apply flags.
 * Output:
 *  return         int     0 - OK.
 *                         1 - Error.
 */
int ed_Telcor(Observation *ob, Subarray *sub, int cif, int ut, int itel, int doflag)
{
  Integration *integ;   /* The integration to be flagged */
  Baseline *bptr;       /* Pointer into sub->base[] */
  Visibility *vis;      /* Pointer into integ->visp[] */
  int base;             /* Baseline index */
/*
 * Get the integration to be flagged.
 */
  integ = &sub->integ[ut];
/*
 * Record the correction flag.
 */
  integ->icor[cif].tcor[itel].bad = doflag;
/*
 * Mark the per-baseline sums of weights as stale.
 */
  flag_baseline_weights(ob, cif);
/*
 * Flag or unflag all visibilities that lie on baselines of telescope 'itel'.
 */
  if(ob_ready(ob, OB_GETIF, NULL) && cif==ob->stream.cif) {
    vis = integ->vis;
    bptr = sub->base;
    if(doflag) {
      for(base=0; base<sub->nbase; base++,vis++,bptr++) {
	if(bptr->tel_a == itel)
	  vis->bad |= FLAG_TA;
	else if(bptr->tel_b == itel)
	  vis->bad |= FLAG_TB;
      };
    } else {
      for(base=0; base<sub->nbase; base++,vis++,bptr++) {
	if(bptr->tel_a == itel)
	  vis->bad &= ~FLAG_TA;
	else if(bptr->tel_b == itel)
	  vis->bad &= ~FLAG_TB;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Adjust a telescope correction and propagate the result to all
 * affected visibilities.
 *
 * Input:
 *  ob     Observation *   The observation to be flagged.
 *  sub       Subarray *   The subarray to be flagged.
 *  cif            int     The index of the IF to be corrected. If
 *                         this IF is in memory, its visibilities will
 *                         be corrected.
 *  ut             int     The index of the integration to be flagged.
 *  itel           int     The index of the telescope to be flagged.
 *  amp_cor      float     A scale-factor to apply to the amplitude correction.
 *  phs_cor      float     An offset to add to the phase correction (radians).
 * Output:
 *  return         int     0 - OK.
 *                         1 - Error.
 */
int adj_Telcor(Observation *ob, Subarray *sub, int cif, int ut, int itel,
	       float amp_cor, float phs_cor)
{
  Integration *integ;   /* The integration to be adjusted */
  Baseline *bptr;       /* Pointer into sub->base[] */
  Visibility *vis;      /* Pointer into integ->visp[] */
  Telcor *tcor;         /* The container of the corrections to be adjusted */
  int base;             /* Baseline index */
/*
 * Get the integration to be adjusted.
 */
  integ = &sub->integ[ut];
/*
 * Disallow zero and negative amplitude corrections.
 */
  if(amp_cor <= 0.0f)
    amp_cor = 1.0f;
/*
 * Get the telescope correction container.
 */
  tcor = &integ->icor[cif].tcor[itel];
/*
 * Adjust the recorded corrections.
 */
  tcor->amp_cor *= amp_cor;
  tcor->phs_cor += phs_cor;
/*
 * Mark the corresponding per-baseline sums of weights as stale.
 */
  flag_baseline_weights(ob, cif);
/*
 * Adjust the corrections applied to all visibilities that lie on
 * baselines of telescope 'itel'.
 */
  if(ob_ready(ob, OB_GETIF, NULL) && cif==ob->stream.cif) {
    vis = integ->vis;
    bptr = sub->base;
    for(base=0; base<sub->nbase; base++,vis++,bptr++) {
      if(bptr->tel_a == itel) {
	vis->phs += phs_cor;
	vis->amp *= amp_cor;
	vis->wt /= amp_cor * amp_cor;
      } else if(bptr->tel_b == itel) {
	vis->phs -= phs_cor;
	vis->amp *= amp_cor;
	vis->wt /= amp_cor * amp_cor;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Clear the recorded corrections of a given telescope, and propagate the
 * results to all effected visibilities.
 *
 * Input:
 *  ob     Observation *   The observation to be flagged.
 *  sub       Subarray *   The subarray to be flagged.
 *  cif            int     The index of the IF of the correction to be
 *                         cleared. If this IF is in memory, its visibilities
 *                         will be modified to appropriately.
 *  ut             int     The index of the integration to be flagged.
 *  itel           int     The index of the telescope to be flagged.
 * Output:
 *  return         int     0 - OK.
 *                         1 - Error.
 */
int clr_Telcor(Observation *ob, Subarray *sub, int cif, int ut, int itel)
{
  Integration *integ;   /* The integration to be flagged */
  Baseline *bptr;       /* Pointer into sub->base[] */
  Visibility *vis;      /* Pointer into integ->visp[] */
  Telcor *tcor;         /* The container of the corrections to be removed */
  int base;             /* Baseline index */
  float amp_cor,phs_cor;/* The amp and phase corrections to be removed */
/*
 * Get the integration to be uncorrected.
 */
  integ = &sub->integ[ut];
/*
 * Get the telescope correction container.
 */
  tcor = &integ->icor[cif].tcor[itel];
/*
 * Get local copies of the corrections to be removed.
 */
  amp_cor = tcor->amp_cor;
  phs_cor = tcor->phs_cor;
/*
 * Disallow zero and negative amplitude corrections.
 */
  if(amp_cor <= 0.0f)
    amp_cor = 1.0f;
/*
 * Clear the recorded corrections.
 */
  tcor->amp_cor = 1.0f;
  tcor->phs_cor = 0.0f;
/*
 * Mark the corresponding per-baseline sums of weights as stale.
 */
  flag_baseline_weights(ob, cif);
/*
 * Remove the corrections from all visibilities that lie on baselines of
 * telescope 'itel'.
 */
  if(ob_ready(ob, OB_GETIF, NULL) && cif==ob->stream.cif) {
    vis = integ->vis;
    bptr = sub->base;
    for(base=0; base<sub->nbase; base++,vis++,bptr++) {
      if(bptr->tel_a == itel) {
	vis->phs -= phs_cor;
	vis->amp /= amp_cor;
	vis->wt *= amp_cor * amp_cor;
      } else if(bptr->tel_b == itel) {
	vis->phs += phs_cor;
	vis->amp /= amp_cor;
	vis->wt *= amp_cor * amp_cor;
      };
    };
  };
  return 0;
}

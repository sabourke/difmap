#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vlbconst.h"
#include "logio.h"
#include "obs.h"
#include "slfcal.h"

/*
 * Types specific to the self-cal functions.
 */
typedef struct {  /* Weighted real-imag representation of complex number */
  float re;       /* Real part */
  float im;       /* Imaginary part */
  float wt;       /* Weight part */
} Scvis;

typedef struct {  /* Corrections for a single telescope and time */
  float amp_cor;
  float phs_cor;
  float weight;
} Cor;

typedef struct {  /* Container for all solutions at a specified time */
  double begut;   /* Start UT of bin (seconds) */
  double endut;   /* End UT of bin (seconds) */
  Cor *cors;      /* Array of nstat corrections */
} Solns;

/* Define a container for dynamic memory required by self-cal routines */

typedef struct {
  Scvis *memblk;  /* Memory for rows of 'nvis' */
  Scvis *gain;    /* 'nstat' complex gain corrections */
  Scvis *gnew;    /* 'nstat' intermediary new estimates gains */
  Scvis **nvis;   /* 2D array of weighted,model-normalized visibilities */
  Cor *cors;      /* 'nstat' Amplitude and phase corrections */
  Solns *solns;   /* Array of correction arrays for each solution bin */
  int *usable;    /* 'nbase' flags designating usable baselines */
  int *iwrk;      /* Int work array */
  int nstat;      /* Number of baselines */
  int nbase;      /* Number of stations */
  int nbin;       /* The number of solution bins */
} Scalmem;

static Scalmem *new_Scal(Subarray *sub, double utint, int doone);
static Scalmem *del_Scal(Scalmem *scal);   /* Scalmem destructor */
static float slfdif(Scvis **nvis, Scvis *gain, int nstat);
static void getgain(Subarray *sub, Scvis **nvis, Scvis *gain, Scvis *gnew,
		    int doamp, int dophs, float slfgain);
static int endbin(Subarray *sub, int uta, double utint);
static int count_bins(Subarray *sub, double utint);
static int get_cors(Subarray *sub, int isbad, int dophs, float maxphs,
		    int doamp, float maxamp, Scvis *gain, Cor *cors);
static void apply_cors(Subarray *sub, int cif, int uta, int utb,
		       int doamp, int dophs, Cor *cors);
static void rep_cors(Observation *ob, int isub, Cor *cors,
		     int doamp, int dophs);
static void apply_solns(Subarray *sub, Scalmem *scal, int cif, float solint,
			int doamp, int dophs);
static double get_area(double xa, double xb, double sigma);
static void sum_ratios(Subarray *sub, int ut, float gaufac, Scalmem *scal);
static int get_usable(Observation *ob, Subarray *sub, int ut,
		      float uvmin, float uvmax, int mintel, int doflag,
		      int *usable, int *telnum, int *nbadtel);
static int count_tel(Observation *ob, Subarray *sub, int ut, int doclose,
		     int doflag, int *usable, int *telnum);
static float norm_cors(Subarray *sub, int cif, Cor *cors);
static int slfsub(Observation *ob, int isub, float gauval, float gaurad,
		  float solint, int doamp, int dophs, int dofloat, int mintel,
		  int doflag, int doone, float maxamp, float maxphs,
		  float uvmin, float uvmax, int *flagged);

static Scvis czero={0.0f,0.0f,0.0f}; /* Used to initialize Scvis elements */

/*.......................................................................
 * Self-calibrate an observation with respect to the established and
 * model after establishing the tentative model.
 *
 * The self-cal technique used (applied separately to each sub-array of
 * each IF), uses least squares to minimize equation 9.5 of
 * chapter 9 ("SELF-CALIBRATION by T.Cornwell and E.B. Formalont) of
 * "Synthesis Imaging in Radio Astronomy", eds. R.A.Perley, F.R.Schwab,
 * A.H.Bridle, published in 1989 by the "Astronomical Society of the Pacific".
 *
 * Note that instead of summing equation 9.5 wrt time for each solution bin
 * a weighted mean of Xij(t) is formed, weighted by 1/(variance of
 * Xij(t)). Also note that over a given solution bin the complex
 * antenna gains are assumed to be constant. Smoothing is applied to
 * the resulting incremental complex gain solutions to remove the
 * stair-step sampling artefacts.
 *
 *
 * Input:
 *   ob   Observation *  The parent observation.
 *   isub         int    The index of the sub-array to be corrected, or
 *                       -1 if all sub-arrays are to be corrected.
 *   doall        int    If true correct all IFs. If 0, correct only the
 *                       current IF.
 *   gauval     float    The value of the weighting gaussian at UV radius
 *                       gaurad, between 0 and 1. If <=0 or >=1, no gaussian
 *                       taper is applied.
 *   gaurad     float    The radius (wavelengths) in the UV plane at which
 *                       the gaussian weighting function has value 'gauval'.
 *                       If <=0, no gaussian taper is applied.
 *   solint     float    The piece-wise amplitude solution interval (minutes).
 *                       solint <= 0.0 will result in separate solutions
 *                       for each integration.
 *   doamp        int    If true then do amplitude corrections.
 *   dophs        int    If true then do phase corrections.
 *   dofloat      int    Gain corrections are usually normalised, to stop
 *                       the flux scale from wandering. If dofloat is true
 *                       (dofloat != 0) then the normalisation will be
 *                       turned off to allow the flux scale to float.
 *   mintel       int    The minimum number of telescopes in closed arrays
 *                       required before attempting a solution.
 *   doflag       int    If true, flag un-correctable visibilities.
 *   doone        int    If true, ignore solint and find a single
 *                       overall correction for the whole time range.
 *   maxamp     float    If maxamp>1.0f and any amplitude correction
 *                       is > maxamp or < 1.0/maxamp, then the affected
 *                       solution interval will be left un-corrected.
 *   maxphs     float    If maxphs > 0.0f and any phase correction
 *                       is > maxphs or < -maxphs, then the affected solution
 *                       interval will be left un-corrected.
 *   uvmin      float    The minimum UV radius to take visibilities from.
 *   uvmax      float    The maximum UV radius to take visibilities from.
 * Input/Output:
 *   flagged      int *  If flagged!=NULL then *flagged will be assigned
 *                       the value 1 if any data are flagged, and 0 otherwise.
 * Output:
 *   return       int    0 - no errors.
 */
int slfcal(Observation *ob, int isub, int doall, float gauval, float gaurad,
	   float solint, int doamp, int dophs, int dofloat, int mintel,
	   int doflag, int doone, float maxamp, float maxphs,
	   float uvmin, float uvmax, int *flagged)
{
  Moddif before,after;/* The goodness of fit before and after selfcal */
  int keepif;         /* The index of the current IF */
  int cif;            /* The index of the IF being corrected */
  int is;             /* The index of the sub-array being corrected */
  int ifa, ifb;       /* The range of IF indexes to be corrected */
  int isa, isb;       /* The range of sub-array indexes to be corrected */
  int old_if;         /* State of current IF to be restored on exit */
/*
 * No data flagged yet.
 */
  if(flagged!=NULL)
    *flagged = 0;
/*
 * Check arguments.
 */
  if(!ob_ready(ob, doall ? OB_SELECT : OB_GETIF, "slfcal"))
    return 1;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Sub-array index OK?
 */
  if(isub<-1 || isub>ob->nsub-1) {
    lprintf(stderr, "slfcal: Sub-array index out of range.\n");
    return 1;
  };
/*
 * Determine the range of IF and sub-array indexes to be used.
 */
  ifa = doall ? 0 : ob->stream.cif;
  ifb = doall ? ob->nif-1 : ob->stream.cif;
  isa = isub>=0 ? isub : 0;
  isb = isub>=0 ? isub : ob->nsub-1;
/*
 * Keep a record of the current IF, so that it can be restored.
 */
  keepif = ob->stream.cif;
/*
 * Establish the tentative model.
 */
  if(mergemod(ob, 1))
    return 1;
/*
 * Determine the overall goodness of fit between the model and observed
 * visibilities before the self-cal.
 */
  if(moddif(ob, &before, uvmin, uvmax))
    return 1;
/*
 * Correct one IF at a time.
 */
  for(cif=ifa; (cif=nextIF(ob, cif, 0, 1)) >= 0 && cif<=ifb; cif++) {
/*
 * Ignore unsampled IFs.
 */
    if(ob->ifs[cif].cl==NULL) {
      lprintf(stdout, "\nNot correcting unselected IF %d.\n", cif+1);
    } else {
/*
 * Report progress.
 */
      lprintf(stdout, "\nCorrecting IF %d.\n", cif+1);
/*
 * Read the current IF.
 */
      if(getIF(ob, cif))
	return 1;
/*
 * Mark the per-baseline sums of weights of the IF as stale.
 */
      flag_baseline_weights(ob, cif);
/*
 * Correct one sub-array of the new IF at a time.
 */
      for(is=isa; is<=isb; is++) {
/*
 * Correct the sub-array.
 */
	if(slfsub(ob, is, gauval, gaurad, solint, doamp, dophs, dofloat,
		  mintel, doflag, doone, maxamp, maxphs, uvmin, uvmax,
		  flagged))
	  return 1;
      };
    };
  };
/*
 * Determine the overall goodness of fit between the model and observed
 * visibilities after the self-cal.
 */
  if(moddif(ob, &after, uvmin, uvmax))
    return 1;
/*
 * Report the change in the fit due to self-cal.
 */
  lprintf(stdout,"\n");
  lprintf(stdout,"Fit before self-cal, rms=%fJy  sigma=%f\n",
	  before.rms, sqrt(before.chisq/before.ndata));
  lprintf(stdout,"Fit after  self-cal, rms=%fJy  sigma=%f\n",
	  after.rms, sqrt(after.chisq/after.ndata));
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    return 1;
  return 0;
}

/*.......................................................................
 * Self-calibrate a sub-array of an observation with the UV-plane model,
 * stored along with the visibilities.
 *
 * The self-cal technique used, uses least squares to minimize equation
 * 9.5 of chapter 9 ("SELF-CALIBRATION by T.Cornwell and E.B. Formalont)
 * of "Synthesis Imaging in Radio Astronomy", eds. R.A.Perley, F.R.Schwab,
 * A.H.Bridle, published in 1989 by the "Astronomical Society of the Pacific".
 *
 * Note that instead of summing equation 9.5 wrt time for each
 * solution bin a weighted mean of Xij(t) is formed, weighted by
 * 1/(variance of Xij(t)). Also note that over a given solution bin
 * the complex antenna gains are assumed to be constant. Smoothing is
 * applied to the resulting incremental complex gain solutions to
 * remove the stair-step sampling artefacts.
 *
 *
 * Input:
 *   ob   Observation *  The parent observation.
 *   isub         int    The index of the sub-array to be corrected.
 *   gauval     float    The value of the weighting gaussian at UV radius
 *                       gaurad, between 0 and 1. If <=0 or >=1, no gaussian
 *                       taper is applied.
 *   gaurad     float    The radius (wavelengths) in the UV plane at which
 *                       the gaussian weighting function has value 'gauval'.
 *                       If <=0, no gaussian taper is applied.
 *   solint     float    The piece-wise amplitude solution interval (minutes).
 *                       solint <= 0.0 will result in separate solutions
 *                       for each integration.
 *   doamp        int    If true then do amplitude corrections.
 *   dophs        int    If true then do phase corrections.
 *   dofloat      int    Gain corrections are usually normalised, to stop
 *                       the flux scale from wandering. If dofloat is true
 *                       (dofloat != 0) then the normalisation will be
 *                       turned off to allow the flux scale to float.
 *   mintel       int    The minimum number of telescopes in closed arrays
 *                       required before attempting a solution.
 *   doflag       int    If true, flag un-correctable visibilities.
 *   doone        int    If true, ignore solint and find a single
 *                       overall correction for the whole time range.
 *   maxamp     float    If maxamp>1.0f and any amplitude correction
 *                       is > maxamp or < 1.0/maxamp, then the affected
 *                       solution interval will be left un-corrected.
 *   maxphs     float    If maxphs > 0.0f and any phase correction
 *                       is > maxphs or < -maxphs, then the affected solution
 *                       interval will be left un-corrected.
 *   uvmin      float    The minimum UV radius to take visibilities from.
 *   uvmax      float    The maximum UV radius to take visibilities from.
 * Input/Output:
 *   flagged      int *  If flagged!=NULL then *flagged will be assigned
 *                       the value 1 if any data are flagged, and 0 otherwise.
 * Output:
 *   return       int    0 - no errors.
 */
static int slfsub(Observation *ob, int isub, float gauval, float gaurad,
		  float solint, int doamp, int dophs, int dofloat, int mintel,
		  int doflag, int doone, float maxamp, float maxphs,
		  float uvmin, float uvmax, int *flagged)
{
  const int niter=100;         /* Max number of gradient-search iterations */
  const float slfgain=0.5f;    /* Loop gain of gradient search */
  const float epsilon=1.0e-6f; /* Acceptable relative change in residuals */
  Subarray *sub;      /* The descriptor of the sub-array being corrected */
  Solns *soln;        /* Pointer to solution bin container in scal->solns[] */
  int ut;             /* The index of the integration being processed */
  Scvis *gain;        /* A gain correction from scal->gain[] */
  Scvis *ctmp;        /* Element of 'scal->nvis' */
  int ita, itb;       /* The numbers of the telescopes on a given baseline */
  Scalmem *scal;      /* Container of dynamically allocated arrays */
  float ini_res;      /* Solution residual before any fitting */
  float old_res;      /* Residual from previous iteration */
  float new_res=0.0f; /* Residual after latest iteration */
  float gfac=0.0f;    /* The reciprocal variance of the gaussian taper */
  double utint;       /* The solution interval */
  double utmid;       /* Mid point of solution interval */
  int uta, utb;       /* The start and end integration indexes in a bin */
  int n_ut;           /* The number of usable integrations in a bin */
  int iter;           /* Current iteration number */
  int nbadtel=0;      /* The number of bad telescope corrections */
  int nbadsol=0;      /* The number of un-usable solution intervals */
/*
 * Get the descriptor of the specified sub-array.
 */
  if(isub<0 || isub>= ob->nsub) {
    lprintf(stderr, "slfsub: Sub-array index out of bounds.\n");
    return 1;
  };
  sub = &ob->sub[isub];
/*
 * Convert the requested solution interval to seconds.
 * 'utint' will double as a flag, in that if no UT interval was selected
 * (ie. solint <= 0.0) it will be made 0.
 */
  utint = solint * 60.0f;
  utint = (utint <= 1.0) ? 0.0 : utint;
/*
 * Allocate all memory required in the self-calibration.
 */
  scal=new_Scal(sub, utint, doone);
  if(scal == NULL)
    return -1;
/*
 * Record whether a gaussian taper was specified and if so work out the
 * -ve reciprocal of the gaussian variance.
 */
  if(gaurad > 0.0 && gauval > 0.0 && gauval < 1.0) {
/*
 * Convert from wavelengths to the units of the recorded UV coordinates.
 */
    gaurad /= ob->stream.uvscale;
    gfac = log(1.0f-gauval)/gaurad/gaurad;
  };
/*
 *---------------------------------------------------------------------
 * Solve for telescope complex gain errors over each solution interval.
 *---------------------------------------------------------------------
 */
  soln = scal->solns;
  for(uta=0; uta<sub->ntime; uta=utb+1,soln += utint>0.0?1:0) {
/*
 * Determine the end UT of the next solution bin.
 */
    utb = doone ? sub->ntime-1 : endbin(sub, uta, utint);
/*
 * Record the start and end times of the solution bin.
 * Take care when finding midut to precision losses.
 */
    utmid = sub->integ[uta].ut + (sub->integ[utb].ut - sub->integ[uta].ut)/2.0;
    soln->begut = utmid-(utint/2.0);
    soln->endut = utmid+(utint/2.0);
/*
 * Initialize the 'nvis' model-normalized visibility array to zero for the
 * new bin.
 * NB. The nvis 2-D array is formed from the 1D array scal->memblk, which
 * has sub->nstat * sub->nstat elements.
 */
    ctmp = scal->memblk;
    for(ita=0;ita<sub->nstat; ita++)
      for(itb=0;itb<sub->nstat; itb++)
	*ctmp++ = czero;
/*
 * Process UTs in the new bin.
 */
    n_ut = 0;
    for(ut=uta; ut <= utb; ut++) {
/*
 * Determine which visibilities in the current integration are usable
 * for self-calibration, and whether there is sufficient data for
 * a solution. Also flag un-correctable data if doflag is true.
 */
      if(get_usable(ob, sub, ut, uvmin, uvmax, mintel, doflag, scal->usable,
		    scal->iwrk, &nbadtel)==0) {
/*
 * Count usable integrations.
 */
	n_ut++;
/*
 * Accumulate the weighted sums of the ratio of observed and model
 * visibilities per baseline, over integrations in the current solution bin.
 */
	sum_ratios(sub, ut, gfac, scal);
      };
    };
/*
 * If non of the integrations within the latest bin were usable,
 * skip to the next bin.
 */
    if(!n_ut)
      continue;
/*
 * Convert weighted observed/model visibility sums into weighted means
 * by dividing by the sums of weights.
 */
    for(ita=0; ita<sub->nstat; ita++) {
      ctmp = &scal->nvis[ita][0];
      for(itb=0; itb<sub->nstat; itb++, ctmp++) {
	float wt = ctmp->wt;
	if(wt > 0.0f) {
	  ctmp->re /= wt;
	  ctmp->im /= wt;
	};
      };
    };
/*
 * Temporarily set gain corrections to 1+0i with 0 weight and determine
 * the starting residual.
 */
    gain = &scal->gain[0];
    for(ita=0; ita<sub->nstat; ita++, gain++) {
      gain->re = 1.0f;
      gain->im = 0.0f;
      gain->wt = 0.0f;
    };
    ini_res = slfdif(scal->nvis,scal->gain,sub->nstat);
/*
 * Initial estimate of telescope gain corrections comes from the weighted
 * mean normalized visibility over the respective telescope.
 */
    getgain(sub, scal->nvis, scal->gain, scal->gnew, doamp, dophs, 1.0f);
/*
 * Compute residual for this estimate.
 */
    old_res = slfdif(scal->nvis,scal->gain,sub->nstat);
/*
 * Iterate for better gain solutions.
 */
    for(iter=0; iter<niter; iter++) {
/*
 * Get better gain estimate gain[ita] for each telescope ita.
 */
      getgain(sub, scal->nvis, scal->gain, scal->gnew, doamp, dophs, slfgain);
/*
 * Determine the resulting residuals to the new fit.
 */
      new_res = slfdif(scal->nvis,scal->gain,sub->nstat);
/*
 * Compare change in residuals from this iteration with the
 * initial residual.
 */
      if(fabs(new_res-old_res) <= epsilon*ini_res)
	break;
      old_res = new_res;
    };
/*
 * Convert the complex reciprocal gains of scal->gain[] to amplitude and
 * phase corrections. Store these in soln->cors and check them against
 * any user-specified limits. Substitute zero-weight unit corrections,
 * if the solution turns out to be un-usable.
 */
    if(get_cors(sub, ini_res<new_res, dophs,maxphs, doamp,maxamp,
		scal->gain, soln->cors)) {
      nbadsol++;  /* Count the number of un-usable solution intervals */
    }
/*
 * The solutions were OK. Apply them now if no interpolation or
 * smoothing is required.
 */
    else {
/*
 * Apply corrections now if no finite solution interval was requested.
 */
      if(doone || utint<=0.0)
	apply_cors(sub, ob->stream.cif, uta, utb, doamp, dophs, soln->cors);
    };
  };          /* end of loop over UT bins. */
/*
 * If smoothed interpolated solutions are required, smooth and
 * interpolate them onto the olbservation time grid, and apply the
 * resulting corrections to the observed data.
 */
  if(utint>0.0 && !doone)
    apply_solns(sub, scal, ob->stream.cif, solint, doamp, dophs);
/*
 * Report flagging operations.
 */
  if(nbadtel>0) {
    lprintf(stdout,
        " A total of %d telescope corrections were %s in sub-array %d.\n",
	 nbadtel, doflag ? "flagged" : "ignored", isub+1);
  };
/*
 * Normalize absolute gain corrections?
 */
  if(doamp && !dofloat) {
    lprintf(stdout, " Amplitude normalization factor in sub-array %d: %g\n",
	    isub+1, norm_cors(sub, ob->stream.cif, scal->cors));
  };
/*
 * Report the corrections if single overall corrections were requested.
 */
  if(doone)
    rep_cors(ob, isub, soln->cors, doamp, dophs);
/*
 * Finished - free all dynamically allocated memory and return.
 */
  del_Scal(scal);
/*
 * Any data flagged?
 */
  if(flagged!=NULL && doflag && nbadtel > 0)
    *flagged = 1;
  return 0;
}

/*.......................................................................
 * Compute the residuals of a self-cal fit from the latest estimate of
 * telescope complex gains and normalized visibilities.
 *
 * Input:
 *  nvis   Scvis **  The model normalized weighted mean baseline-visibility
 *                   matrix. Dimension nvis[nstat][nstat] .
 *  gain   Scvis *   Complex telescope gains. Dimension gain[nstat].
 *  nstat    int     Number of stations.
 * Output:
 *  return float     The residual.
 */
static float slfdif(Scvis **nvis, Scvis *gain, int nstat)
{
  Scvis *ga, *gb;    /* Complex gains of antennas ita,itb */
  Scvis *ctmp;       /* Pointer to nvis[ita][itb] */
  int ita,itb;       /* Indexes of telescope pair */
  float resid=0.0f;  /* Residual of fit. */
  float wtsum=0.0f;  /* Sum of weights */
  float re;          /* Real part of complex temporary */
  float im;          /* Imaginary part of complex temporary */
/*
 * Residual calculated as weighted mean of:
 * sqrmodulus(GAIN_i * conjugate(GAIN_j) - NVIS_ij)
 * This is eqn. 9.5 of chapter 9 ("SELF-CALIBRATION by T.Cornwell and
 * E.B. Formalont) of "Synthesis Imaging in Radio Astronomy", eds.
 * R.A.Perley, F.R.Schwab, A.H.Bridle, published in 1989 by the
 * "Astronomical Society of the Pacific", abitrarily normalized by the
 * sum of the term outside the square modulus.
 */
  for(ita=0; ita<nstat; ita++) {
    ga = &gain[ita];
    gb = &gain[0];
    ctmp = &nvis[ita][0];
    for(itb=0; itb<nstat; itb++, ctmp++, gb++) {
/*
 * ga * cnj(gb) - ctmp.
 */
      re = ga->re*gb->re + ga->im*gb->im - ctmp->re;
      im = ga->im*gb->re - ga->re*gb->im - ctmp->im;
/*
 * Weighted sum of square modulus of the above.
 */
      resid += ctmp->wt * (re*re+im*im);
      wtsum += ctmp->wt;
    };
  };
/*
 * Convert weighted sum into weighted mean by dividing by sum of weights.
 */
  if(resid > 0.0f && wtsum > 0.0f)
    resid /= wtsum;
  else
    resid = 0.0f;
  return resid;
}

/*.......................................................................
 * Determine new estimates for the complex gain corrections.
 *
 * Input/Output:
 *  gain    Scvis *  Complex gain corrections.
 *                    On input : Previous estimates.
 *                    On Output: New estimates.
 * Input:
 *  sub  Subarray *  The descriptor of the sub-array being corrected.
 *  nvis    Scvis ** Array of weighted sum of model-normalized visibilities
 *  gnew    Scvis *  Work array for intermediary new estimates of gains.
 *  doamp     int    If true, then the new estimates will include
 *                   non-unity amplitude estimates.
 *  dophs     int    If true, then the new estimates will include
 *                   non-zero phase estimates.
 *  slfgain float    Iteration gain.
 */
static void getgain(Subarray *sub, Scvis **nvis, Scvis *gain, Scvis *gnew,
		    int doamp, int dophs, float slfgain)
{
  static Station *tel;  /* The descriptor of a telescope */
  static Scvis *ga, *gb;/* Telescope gains for telescope pair a,b */
  static Scvis *gn;     /* Pointer into gnew[] */
  static Scvis *ctmp;
  static Scvis top;     /* Complex temporary (numerator of quotient) */
  static float bot;     /* Temporary real (denominator of quotient) */
  static float amp;     /* An amplitude */
  static float wt;      /* A weight */
  static float wt_sum;  /* Sum of weights */
  static int ita,itb;   /* Indexes of baseline telescope pair a,b */
/*
 * Get better gain estimate gnew[ita] for each telescope ita.
 */
  ga = &gain[0];
  gn = &gnew[0];
  for(ita=0; ita<sub->nstat; ita++, ga++, gn++) {
    ctmp = &nvis[ita][0];
    gb = &gain[0];
/*
 * Determine a least-squares solution for GAIN_a to eqution:
 *
 *   SUM_a_b(WEIGHT_a_b * |NVIS_a_b - GAIN_a.GAIN_b*|^2).
 *
 * With 'a' fixed, and GAIN_[b!=a] fixed. This results in the equation:
 *
 *   GAIN_b = SUM_b(WEIGHT_a_b * GAIN_b * NVIS_a_b)
 *            -------------------------------------
 *               SUM_b(WEIGHT_a_b * |GAIN_b|^2)
 *
 * Perform the above numerator (top) and denominator (bot) sums.
 */
    top = czero;
    bot = wt_sum = 0.0f;
    for(itb=0; itb<sub->nstat; itb++,ctmp++,gb++) {
      wt = ctmp->wt;
      if(wt > 0.0f) {
	top.re += wt * (gb->re*ctmp->re - gb->im*ctmp->im);
	top.im += wt * (gb->re*ctmp->im + gb->im*ctmp->re);
	bot += wt * (gb->re*gb->re + gb->im*gb->im);
	wt_sum += wt;
      };
    };
/*
 * The new gain estimate is a weighted combination of the previous
 * and latest estimates.
 */
    if(bot > 0.0f) {
      gn->re = (1.0f-slfgain) * ga->re + slfgain * top.re/bot;
      gn->im = (1.0f-slfgain) * ga->im + slfgain * top.im/bot;
      gn->wt = wt_sum;
    };
/*
 * Replace bad solutions with the previous best gain estimate.
 */
    if(bot <= 0.0f || (gn->re==0.0f && gn->im==0.0f))
      *gn = *ga;
  };
/*
 * Copy the new gain estimates to gain[]. At the same time remove
 * phase corrections if dophs is false and/or amplitude corrections
 * if doamp is false.
 */
  gn = &gnew[0];
  gb = &gain[0];
  tel = sub->tel;
  for(ita=0; ita<sub->nstat; ita++, gn++, gb++, tel++) {
    if(gn->wt>0.0f) {
      amp = sqrt(gn->re*gn->re + gn->im*gn->im);
/*
 * Return telescope correction to 1+0i if changes are not allowed.
 */
      if(tel->antfix) {
	gn->re = 1.0f;
	gn->im = 0.0f;
      }
/*
 * Remove phase correction if requested.
 */
      else if(!dophs) {
	gn->re = amp;
	gn->im = 0.0f;
      }
/*
 * Remove amplitude correction if requested.
 */
      else if(!doamp) {
	gn->re /= amp;
	gn->im /= amp;
      };
    };
/*
 * Copy the new estimate.
 */
    *gb = *gn;
  };
  return;
}

/*.......................................................................
 * Scalmem constructor.
 *
 * Input:
 *  sub     Subarray *  The descriptor of the sub-array being corrected.
 *  utint     double    The solution interval (seconds).
 *  doone        int    If true, only one overall correction is to be
 *                      determined.
 * Output:
 *  return   Scalmem *  The initialised Scalmem container, or NULL on
 *                      error.
 */
static Scalmem *new_Scal(Subarray *sub, double utint, int doone)
{
  Scalmem *scal; /* The returned container */
  int i;
  static char *errmsg="Insufficient memory to self-cal\n";
/*
 * Allocate the container.
 */
  scal = (Scalmem *) malloc(sizeof(Scalmem));
  if(scal == NULL) {
    lprintf(stderr, errmsg);
    return scal;
  };
/*
 * Make all the member pointers null so that del_Scal() can be called
 * before the container has been fully initialised to clean up on
 * errors. NB calloc() can't be used to do this since some machines
 * don't represent pointers as zero bytes. Only assignment with zero is
 * guaranteed to work by the ANSI standard.
 */
  scal->gain = 0;
  scal->gnew = 0;
  scal->nvis = 0;
  scal->usable = 0;
  scal->iwrk = 0;
  scal->memblk = 0;
  scal->cors = 0;
  scal->solns = 0;
/*
 * Install a record of the number of stations and baselines in the
 * observation to be self-cal'd.
 */
  scal->nbase = sub->nbase;
  scal->nstat = sub->nstat;
/*
 * Now allocate the various arrays.
 * First those that depend on the number of baselines.
 */
  scal->usable = (int *) malloc(sub->nbase * sizeof(int));
/*
 * Pause to check for memory allocation errors.
 */
  if(scal->usable==NULL) {
    lprintf(stderr, errmsg);
    return del_Scal(scal);
  };
/*
 * Now vectors that depend on the number of stations.
 */
  scal->gain  = (Scvis *) malloc(sub->nstat * sizeof(Scvis));
  scal->gnew  = (Scvis *) malloc(sub->nstat * sizeof(Scvis));
  scal->cors  = (Cor *) malloc(sub->nstat * sizeof(Cor));
  scal->nvis  = (Scvis **) malloc(sub->nstat * sizeof(Scvis *));
  scal->iwrk  = (int *) malloc(sub->nstat * sizeof(int));
/*
 * Check for memory allocation errors.
 */
  if(scal->gain == NULL || scal->gnew == NULL || scal->cors == NULL ||
     scal->nvis == NULL || scal->iwrk == NULL) {
    lprintf(stderr, errmsg);
    return del_Scal(scal);
  };
/*
 * Now allocate an array of size sub->nstat*sub->nstat, to be
 * indexed into by the scal->nvis array of pointers to form
 * a square matrix.
 */
  scal->memblk = (Scvis *) malloc(sizeof(Scvis) * sub->nstat * sub->nstat);
  if(scal->memblk == NULL) {
    lprintf(stderr, errmsg);
    return del_Scal(scal);
  };
/*
 * Make the scal->nvis array of pointers to Scvis, point to successive
 * linear arrays in scal->memblk to form a 2D square matrix.
 */
  scal->nvis[0] = scal->memblk;
  for(i=1; i<sub->nstat; i++)
    scal->nvis[i] = scal->nvis[i-1]+sub->nstat;
/*
 * Determine the number of solution containers required.
 */
  scal->nbin = (utint>0.0 && !doone) ? count_bins(sub, utint) : 1;
/*
 * Allocate memory for 'scal->nbin' correction array containers.
 */
  scal->solns = (Solns *) malloc(scal->nbin * sizeof(Solns));
  if(scal->solns == NULL) {
    lprintf(stderr, errmsg);
    return del_Scal(scal);
  };
/*
 * Allocate sufficient memory for scal->nbin * nstat corrections.
 */
  scal->solns[0].cors = calloc((size_t) scal->nbin * sub->nstat, sizeof(Cor));
  if(scal->solns[0].cors == NULL) {
    lprintf(stderr, errmsg);
    return del_Scal(scal);
  };
/*
 * Assign blocks of nstat corrections to each solution bin container.
 */
  for(i=1; i<scal->nbin; i++)
    scal->solns[i].cors = scal->solns[i-1].cors + sub->nstat;
/*
 * Initialize all corrections to have zero weight (until used).
 */
  for(i=0; i<scal->nbin*sub->nstat; i++)
    scal->solns[0].cors[i].weight = 0.0f;
  return scal;
}

/*.......................................................................
 * Scalmem destructor. Release all memory contained by the Scalmem
 * container class and then free the container itself.
 *
 * Input:
 *  scalmem  Scalmem *  The Scalmem container (+ contained memory) to
 *                      be deleted.
 * Output:
 *  return Scalmem *    Always NULL.
 */
static Scalmem *del_Scal(Scalmem *scal)
{
  if(scal) {
    if(scal->memblk)
      free(scal->memblk);
    if(scal->gain)
      free(scal->gain);
    if(scal->gnew)
      free(scal->gnew);
    if(scal->cors)
      free(scal->cors);
    if(scal->nvis)
      free(scal->nvis);
    if(scal->usable)
      free(scal->usable);
    if(scal->iwrk)
      free(scal->iwrk);
    if(scal->solns) {
      if(scal->solns[0].cors)
	free(scal->solns[0].cors);
      free(scal->solns);
    };
    free(scal);
  };
  return NULL;
}

/*.......................................................................
 * Given the index of the first integration of a solution bin, work out
 * the index of the last integration that falls within the solution
 * interval.
 *
 * Input:
 *  sub    Subarray *  The descriptor of the sub-array being sampled.
 *  uta         int    The index of the first integration in the bin.
 *  utint    double    The solution interval in seconds, or 0.0 if no
 *                     binning is being used.
 * Output:
 *  return      int    The index of the last integration in the bin.
 */
static int endbin(Subarray *sub, int uta, double utint)
{
  Integration *integ; /* The descriptor of an integration in the bin */
  double begut; /* Start UT of bin */
  double endut; /* End UT of bin */
  int utb;      /* Index of last integration in bin */
/*
 * Get the descriptor of the start integration.
 */
  integ = &sub->integ[uta];
/*
 * Determine the start and end UT of the next solution bin.
 */
  if(utint>0.0) {
    begut = utint * floor(integ->ut/utint);
    endut = begut + utint;
/*
 * Locate the end integration index of the bin.
 */
    for(utb=uta; utb<sub->ntime && integ->ut <= endut; utb++,integ++);
    utb--;
  } else {
    utb=uta;
    begut = endut = integ->ut;
  };
  return utb;
}

/*.......................................................................
 * Count the number of sampled solution bins in a given sub-array.
 *
 * Input:
 *  sub    Subarray *  The descriptor of the sub-array being sampled.
 *  utint    double    The solution interval (seconds).
 * Output:
 *  return      int    The number of solution bins counted.
 */
static int count_bins(Subarray *sub, double utint)
{
  int uta;       /* The index of the start integration of a bin */
  int nbin=0;    /* The number of bins counted */
/*
 * Loop over all times sampled by the observation and work out which bins
 * they fall into.
 */
  for(uta=0; uta<sub->ntime; uta=endbin(sub, uta, utint)+1)
    nbin++;
  return nbin;
}

/*.......................................................................
 * Record and apply self-calibration corrections for one integration of
 * a given sub-array.
 *
 * Input:
 *  sub     Subarray *  The descriptor of the sub-array being corrected.
 *  cif          int    The index of the IF for which the corrections are
 *                      to be recorded.
 *  uta          int    The index of the first integration to be corrected.
 *  utb          int    The index of the last integration to be corrected.
 *  doamp        int    If true, apply amplitude corrections.
 *  dophs        int    If true, apply phase corrections.
 *  cors         Cor *  The array of sub->nstat corrections.
 */
static void apply_cors(Subarray *sub, int cif, int uta, int utb,
		       int doamp, int dophs, Cor *cors)
{
  Integration *integ;/* The descriptor of the integration being corrected */
  Telcor *ocor;     /* Pointer into the array of applied corrections */
  Cor *icor;        /* Pointer into the arrray of input corrections */
  int base;         /* The index of the baseline being corrected */
  int ita,itb;      /* The indexes of telescopes on a given baseline */
  int ut;           /* The index of the integration being corrected */
/*
 * Correct each integration.
 */
  integ = &sub->integ[uta];
  for(ut=uta; ut<=utb; ut++,integ++) {
/*
 * Get the array of visibilities that are to be corrected and the
 * parallel baseline array.
 */
    Visibility *vis = integ->vis;
    Baseline *bptr = sub->base;
/*
 * Correct the visibility associated with each baseline.
 */
    for(base=0; base<sub->nbase; base++,vis++,bptr++) {
/*
 * Get the indexes of the stations of the baseline that measured the visibility.
 */
      ita = bptr->tel_a;
      itb = bptr->tel_b;
/*
 * Amplitude corrections:
 */
      if(doamp) {
	float ftmp = cors[ita].amp_cor * cors[itb].amp_cor;
	vis->amp *= ftmp;
	vis->wt /= ftmp * ftmp;   /* NB. wt = 1/amp_err^2 */
      };
/*
 * Now the phase corrections.
 */
      if(dophs)
	vis->phs += cors[ita].phs_cor - cors[itb].phs_cor;
    };
/*
 * Record the corrections. Note that if the accumulated correction
 * is negative then the telescope has not been calibrated before
 * and should be unflagged now if a correction is being made.
 */
    icor = cors;
    ocor = integ->icor[cif].tcor;
    for(ita=0; ita<sub->nstat; ita++,ocor++,icor++) {
      if(dophs)
	ocor->phs_cor += icor->phs_cor;
      if(doamp) {
	ocor->amp_cor *= icor->amp_cor;
/*
 * Remove non-calib flag if a correction was made.
 */
	if(ocor->amp_cor < 0.0f && icor->weight > 0.0f)
	  ocor->amp_cor *= -1.0;
      };
    };
  };
  return;
}

/*.......................................................................
 * For each integration take the corrections contained in the scal->solns
 * solutions, smooth, interpolate and apply them to the observed data.
 *
 * Input:
 *  sub    Subarray *  The descriptor of the sub-array being corrected.
 *  scal    Scalmem *  The container of dynamic objects required by
 *  cif         int    The index of the IF for which the corrections are
 *                     to be recorded.
 *  solint    float    The solution interval in minutes.
 *  doamp       int    If true, apply amplitude corrections.
 *  dophs       int    If true, apply phase corrections.
 *                     self-cal.
 */
static void apply_solns(Subarray *sub, Scalmem *scal, int cif, float solint,
			int doamp, int dophs)
{
  Integration *integ; /* The descriptor of the integration being corrected */
  static const double nsigma=2.5;     /* Max number of standard deviations. */
  double sigma;   /* The standard deviation of the smoothing gaussian */
  double maxoff;  /* The time difference defining a solution as a neighbour */
  double utval;   /* The time of the latest integration */
  Solns *soln;    /* The solution within which the latest integration falls */
  Cor *ocor;      /* Pointer into scal->cors (output corrections array) */
  Cor *icor;      /* Pointer into one of the input corrections arrays */
  int ut;         /* The index of the integration to be corrected */
  int itel;       /* The index of the station correction */
  int sa;         /* First solution in range of the current integration */
  int sb;         /* Next solution in range of the current integration */
/*
 * Determine the standard deviation of the smoothing gaussian such that
 * the fourier transform of the gaussian has a half-width at half power of
 * 1/(2.utint) (Nyquist sampling). This comes out as:
 *  sigma = 2 * utint * sqrt(ln(sqrt(2))) / pi.
 */
  sigma = solint * 0.37478125;
/*
 * Calculate the boundary beyond which the gaussian tails are to be
 * considered insignificant.
 */
  maxoff = nsigma*sigma;
/*
 * Determine and apply corrections for each integration.
 * This involves interpolating the corrections from neighbouring solution
 * intervals, weighted by the individual correction weights and
 * by a smoothing Gaussian function.
 */
  sa = 0;
  integ = sub->integ;
  for(ut=0; ut<sub->ntime; ut++,integ++) {
    utval = integ->ut;
    soln = &scal->solns[sa];
/*
 * Initialize the output corrections array.
 */
    ocor = scal->cors;
    for(itel=0; itel<sub->nstat; itel++,ocor++) {
      ocor->amp_cor = 0.0f;
      ocor->phs_cor = 0.0f;
      ocor->weight = 0.0f;
    };
/*
 * Find the next solution whose solution interval lies less than
 * maxoff away from the new integration.
 */
    for(; sa<scal->nbin && (utval-soln->endut)*uttomin >= maxoff;
	sa++,soln++);
/*
 * For each solution whose borders lie within maxoff of the current
 * integration, determine the weight assigned to its corrections,
 * and sum the weighted corrections for each telescope.
 */
    for(sb=sa; sb<scal->nbin && (soln->begut-utval)*uttomin < maxoff;
	sb++,soln++) {
/*
 * Get the area under the smoothing gaussian (centered on utval) within the
 * solution bin.
 */
      double b_start = uttomin * (soln->begut - utval); /* Start time of bin */
      double b_end   = uttomin * (soln->endut - utval); /* End time of bin */
      float area;                     /* The sampled area under the gaussian */
/*
 * Limit the bounds of the binned area to be within maxoff of utval.
 */
      if(b_start < -maxoff)
	b_start = -maxoff;
      if(b_end > maxoff)
	b_end = maxoff;
/*
 * Calculate the area between b_start and b_end.
 */
      area = get_area(b_start, b_end, sigma);
/*
 * Sum corrections into scal->cor[0...nitel-1], using the correction weights
 * and the gaussian area under the solution bin as overall weights.
 */
      ocor = scal->cors;
      icor = soln->cors;
      for(itel=0; itel<sub->nstat; itel++,icor++,ocor++) {
	if(icor->weight > 0.0f) {  /* Usable correction? */
	  float wt = area * icor->weight;
	  ocor->amp_cor += wt * icor->amp_cor;
	  ocor->phs_cor += wt * icor->phs_cor;
	  ocor->weight += wt;
	};
      };
    };
/*
 * Form weighted means from the weighted correction sums.
 */
    ocor = scal->cors;
    for(itel=0; itel<sub->nstat; itel++,ocor++) {
      if(ocor->weight > 0.0f) {  /* Usable correction? */
	ocor->amp_cor /= ocor->weight;
	ocor->phs_cor /= ocor->weight;
      } else {
	ocor->amp_cor = 1.0f;
	ocor->phs_cor = 0.0f;
      };
    };
/*
 * Apply corrections to the latest integration.
 */
    apply_cors(sub, cif, ut, ut, doamp, dophs, scal->cors);
  };
  return;
}

/*.......................................................................
 * Return the approximate area under a unit area gaussian of standard
 * deviation 'sigma', from X=xa,xb. This is determined through
 * interpolating a coarsely tabulated sampling of a rational approximation
 * to the error function (divided by two).
 *
 * Input:
 *  xa     double   The X value of the gaussian, to start at.
 *  xb     double   The X value of the gaussian, to end at.
 *  sigma  double   The standard deviation of the gaussian.
 * Output:
 *  return double   The area between X=xa,xb of the gaussian.
 */
static double get_area(double xa, double xb, double sigma)
{
  enum {ERFSIZ=16};           /* The size of the interpolation table */
  static const double s2=1.4142136; /* sqrt(2.0) */
  static const double nsigma=2.5;   /* Max no. std deviations for gaussian */
  static double erfconv;       /* Maps between coordinate and table index */
  static double erftab[ERFSIZ+1]; /* The tabulated error function */
  static int done=0;           /* True once the table has been initialized */
  static double za,zb;         /* xa and xb converted to x/(2.sigma) */
  static double asgn,bsgn;     /* Unit signs of za and zb */
  static double apos,bpos;     /* Decimal index in lookup table for za and zb */
  static int aind,bind;     /* Integer index in lookup table for za and zb */
  static double a_area,b_area; /* Areas over 0,za and 0,zb */
/*
 * Initialize the error function table up to x=sigma, using a
 * rational approximation to the error function.
 */
  if(!done) {
    double z,t;
    int i;
/*
 * Mark the lookup table as initialized.
 */
    done=1;
/*
 * Determine the mapping between error function argument and integer
 * index into the lookup table.
 */
    erfconv = ((ERFSIZ-1) * s2) / nsigma;
/*
 * Fill the lookup table from rational approximation of error function/2.0.
 */
    for(i=0; i<=ERFSIZ; i++) { /* One extra value simplifies interpolation */
      z = i / erfconv;
      t = 1.0/(1.0 + 0.47047 * z);
      erftab[i] = 0.5 -
	(0.1740121*t*(1.0 + -0.2754975*t*(1.0 + -7.7999287*t))) * exp(-z*z);
    };
  };
/*
 * Change units of xa and xb to z=x/(sqrt(2).sigma).
 * In these units the gaussian becomes: 2/sqrt(pi) * e^(z^2) and
 * are the units used in determining error functions.
 */
  za = xa/(s2 * sigma);
  zb = xb/(s2 * sigma);
/*
 * Determine the area over z=0,za.
 */
  asgn = za < 0.0 ? -1 : 1;   /* Sign of za */
  apos = erfconv * asgn * za;  /* Decimal index into lookup table */
  aind = (int) apos;           /* Truncated integer index into lookup table */
  if(aind<ERFSIZ) {
    double a1 = erftab[aind];   /* First value for interpolation */
    double a2 = erftab[aind+1]; /* Second value for interpolation */
    a_area = a1 + (apos-aind) * (a2 - a1);  /* Interpolated area */
  } else {
    a_area = erftab[ERFSIZ];   /* NB there are actually ERFSIZ+1 entries */
  };
/*
 * Determine the area over z=0,zb.
 */
  bsgn = zb < 0.0 ? -1 : 1;   /* Sign of zb */
  bpos = erfconv * bsgn * zb;  /* Decimal index into lookup table */
  bind = (int) bpos;           /* Truncated integer index into lookup table */
  if(bind<ERFSIZ) {
    double b1 = erftab[bind];   /* First value for interpolation */
    double b2 = erftab[bind+1]; /* Second value for interpolation */
    b_area = b1 + (bpos-bind) * (b2 - b1);  /* Interpolated area */
  } else {
    b_area = erftab[ERFSIZ];   /* NB there are actually ERFSIZ+1 entries */
  };
/*
 * Work out the area between za and zb.
 */
  return fabs(asgn * a_area - bsgn * b_area);
}

/*.......................................................................
 * Report corrections to each telescope of a sub-array.
 *
 * Input:
 *  ob    Observation *   The observation descriptor.
 *  isub          int     The index of the sub-array to report about.
 *  cors          Cor *   The array of ob->sub[isub].nstat corrections to list.
 *  doamp         int     If true, report amplitude corrections.
 *  dophs         int     If true, report phase corrections.
 */
static void rep_cors(Observation *ob, int isub, Cor *cors, int doamp, int dophs)
{
  Subarray *sub; /* The descriptor of the sub-array */
  Cor *cor;      /* The correction being reported */
  int itel;      /* The index of the telescope correction to report */
/*
 * Get the descriptor of the sub-array.
 */
  sub = &ob->sub[isub];
/*
 * Report the individual corrections to the user's terminal.
 */
  if(doamp || dophs) {
    lprintf(stdout, " Telescope %s%s corrections in sub-array %d:",
	    doamp ? "amplitude" : "phase",
	    doamp && dophs ? " and phase" : "", isub+1);
/*
 * List telescope names followed by their corrections.
 */
    cor = cors;
    for(itel=0; itel<sub->nstat; itel++,cor++) {
/*
 * Start a new line when there is no room left on the current line.
 */
      if(itel % (doamp && dophs ? 3 : 4) == 0)
	lprintf(stdout, "\n  ");
/*
 * List the telescope name.
 */
      lprintf(stdout, "%-8s", sub->tel[itel].name);
/*
 * List the amplitude and/or phase corrections for this telescope.
 */
      if(doamp)
	lprintf(stdout, " %5.2f", cor->amp_cor);
      if(dophs)
	lprintf(stdout, "%c%5.2f", doamp ? ',':' ', cor->phs_cor);
/*
 * Mark the corrections of uncorrected telescopes with an asterix.
 */
      lprintf(stdout, "%c    ", cor->weight > 0.0f ? ' ':'*');
    };
/*
 * Terminate the listing.
 */
    lprintf(stdout, "\n\n");
  };
  return;
}

/*.......................................................................
 * Determine the weighted complex ratios of observed/model visibilities
 * of each baseline at a given UT, and sum these into the square
 * nstat x nstat matrix scal->nvis[][] to form a weighted sum of ratios
 * with those of other UTs in the current solution interval.
 *
 * Input:
 *  sub     Subarray *  The descriptor of the sub-array being self-calibrated.
 *  ut           int    The index of the integration to be used.
 *  gaufac     float    The factor to multiply an unscaled UV radius by to
 *                      apply a gaussian taper. Send 0.0f if no taper
 *                      is to be applied.
 * Input/Output:
 *  scal     Scalmem *  The self-cal dynamic memory container.
 *                      1. On input, scal->usable must be initialized
 *                         such that only baselines with usable visibilities
 *                         are flagged as true (!=0).
 *                      2. On output, scal->nvis[][] ratios will have been
 *                         updated.
 */
static void sum_ratios(Subarray *sub, int ut, float gaufac, Scalmem *scal)
{
  Visibility *vis;  /* Pointer into visibility array for this integration */
  Baseline *bptr;   /* Pointer into the baseline descriptor array */
  Station *tel;     /* The array of station descriptors in the sub-array */
  int *usable;      /* Pointer into scal->usable[] */
  Scvis *ctmp;      /* Pointer to nvis[ita][itb] */
  float amp;        /* The weighted obs/model amp ratio of a visibility */
  float phs;        /* The obs-model phase difference of a visibility */
  float wt;         /* The weight of a visibility ratio */
  float re;         /* Real part of complex temporary */
  float im;         /* Imaginary part of complex temporary */
  int base;         /* The index of the baseline being processed */
  int ita,itb;      /* The indexes of telescopes on a given baseline */
/*
 * Get a pointer to the array of telescope descriptors of this sub-array.
 */
  tel = sub->tel;
/*
 * Get a pointer to the array of visibilities for this integration, and
 * the associated baseline descriptor and usable arrays.
 */
  vis = sub->integ[ut].vis;
  bptr = sub->base;
  usable = scal->usable;
/*
 * Get the visibility ratios for each baseline.
 */
  for(base=0; base<sub->nbase; base++,vis++,bptr++,usable++) {
/*
 * Ignore baselines with unusable visibilities.
 */
    if(*usable && vis->modamp!=0.0f) {
/*
 * Get the indexes of the two telescopes in the current baseline.
 */
      ita = bptr->tel_a;
      itb = bptr->tel_b;
/*
 * Get the weight. This is 1/variance of the ratio being weighted, ie.
 * if Vobs is the observed complex visibility and Vmod is the
 * corresponding mode visibility:
 *
 *  Variance(Vobs/Vmod) = Variance(Vobs) / |Vmod|^2
 *
 *  Weight = |Vmod|^2 / Variance(Vobs).
 *
 * The definition of vis->wt is 1/Variance(Vobs), and |Vmod| is the
 * amplitude of the model, leading to the following equation for the weight.
 */
      wt = vis->wt * vis->modamp * vis->modamp;
/*
 * Gaussian taper to apply to the weights?
 */
      if(gaufac < 0.0f) {
	float uu = vis->u;
	float vv = vis->v;
	wt *= 1.0-exp(gaufac*(uu*uu+vv*vv));
      };
/*
 * Apply extra telescope weights.
 */
      wt *= fabs(tel[ita].antwt * tel[itb].antwt);
/*
 * Divide polar representation of visibility by model and multiply by
 * weight.
 */
      amp = wt * vis->amp / vis->modamp;
      phs = vis->phs - vis->modphs;
/*
 * Get the complex representation of amp and phase.
 */
      re = amp * cos(phs);
      im = amp * sin(phs);
/*
 * Sum into baseline visibility matrix.
 */
      ctmp = &scal->nvis[ita][itb];
      ctmp->re += re;
      ctmp->im += im;
      ctmp->wt += wt;
/*
 * Also add to conjugate baseline.
 */
      ctmp = &scal->nvis[itb][ita];
      ctmp->re += re;
      ctmp->im -= im;
      ctmp->wt += wt;
    };
  };
  return;
}

/*.......................................................................
 * Convert the complex reciprocal gains of gain[] to amplitude and
 * phase corrections. Store these in cors[] and check them against
 * any user-specified limits. If the corrections are marked as already
 * bad with the isbad argument, or they turn out to be bad when
 * compared against maxphs and/or maxamp, then the returned
 * corrections are set to zero weight and unit gain, and the function
 * returns an indication that the input corrections were un-usable.
 *
 * Input:
 *  sub  Subarray *  The descriptor of the sub-array being corrected.
 *  isbad     int    If true then the solutions in gain[] are
 *                   considered unusable (see above).
 *  dophs     int    If true, check phase limits.
 *  maxphs  float    If dophs is true and maxphs > 0.0f and any phase
 *                   correction is > maxphs or < -maxphs, then the
 *                   solutions in gain[] are considered un-usable (see
 *                   above).
 *  doamp     int    If true, check amplitude limits.
 *  maxamp  float    If doamp is true, and maxamp>1.0f and any
 *                   amplitude correction is > maxamp or < 1.0/maxamp,
 *                   then the solutions in gain[] are considered
 *                   un-usable (see above).
 *  gain    Scvis *  The array of sub->nstat reciprocal complex
 *                   telescope gain corrections to be translated from.
 * Input/Output:
 *  cors      Cor *  Send an empty array for sub->nstat corrections. On
 *                   output this will contain the telescope amplitude
 *                   and phase corrections corresponding to the
 *                   reciprocal gains in gain[].
 * Output:
 *  return    int    0 - OK.
 *                   1 - isbad was true or the corrections exceeded
 *                       user specified limits.
 */
static int get_cors(Subarray *sub, int isbad, int dophs, float maxphs,
		    int doamp, float maxamp, Scvis *gain, Cor *cors)
{
  Scvis *gptr;  /* Pointer into gain[] */
  Cor *cptr;    /* Pointer into cors[] */
  float minamp; /* Minimum acceptable amplitude correction */
  int doplim;   /* True if phase limits are to be checked */
  int doalim;   /* True if amplitude limits are to be applied */
  int itel;     /* The index of the telescope correction being processed */
/*
 * Which, if any correction limits should be checked?
 */
  doplim = dophs && maxphs > 0.0f;                /* Limit phases? */
  doalim = doamp && maxamp > 1.0f;                /* Limit amplitudes? */
/*
 * Get the minimum acceptable amplitude correction.
 */
  minamp = (doalim && maxamp!=0.0f) ? 1.0f/maxamp : 0.0f;
/*
 * Convert from real-imaginary notation to polar notation and change
 * the corrections from model multipliers to observed data corrections.
 */
  gptr = gain;
  cptr = cors;
  for(itel=0; itel<sub->nstat && !isbad; itel++, cptr++, gptr++) {
    if(gptr->re==0.0f && gptr->im==0.0f) {
      cptr->amp_cor = 1.0f;
      cptr->phs_cor = 0.0f;
      cptr->weight = 0.0f;
    } else {
      cptr->amp_cor = 1.0f/sqrt(gptr->re*gptr->re + gptr->im*gptr->im);
      cptr->phs_cor = -atan2(gptr->im,gptr->re);
      cptr->weight = gptr->wt;
/*
 * Check against user specified correction limits.
 */
      if((doplim && (cptr->phs_cor > maxphs || cptr->phs_cor < -maxphs)) ||
         (doalim && (cptr->amp_cor > maxamp || cptr->amp_cor <  minamp)) )
	isbad = 1;
    };
  };
/*
 * If the corrections are bad, set their weights to zero and assign unit
 * gains to them.
 */
  if(isbad) {
    for(cptr=cors,itel=0; itel<sub->nstat; itel++,cptr++) {
      cptr->amp_cor = 1.0f;
      cptr->phs_cor = 0.0f;
      cptr->weight = 0.0f;
    };
  };
  return isbad;
}

/*.......................................................................
 * Given an integration, optional UV radius limits, and a minimum for
 * for the number of telescopes in closed sub-arrays, determine which
 * visibilities in the integration are usable for self-calibration.
 * Also, if requested, flag either just the visibilities of baselines
 * that are not members of closed telescope sub-arrays, or flag all
 * visibilities if there turn out to be no usable visibities.
 *
 * Input:
 *  ob   Observation *  The descriptor of the parent observation.
 *  sub     Subarray *  The descriptor of the sub-array being corrected.
 *  ut           int    The index of the integration to be processed.
 *  uvmin      float    The min allowed UV radius (Wavelengths).
 *  uvmax      float    The max allowed UV radius (Wavelengths).
 *                      NB. the range will only be applied if the greater of
 *                      uvmin and uvmax > 0.0f;
 *  mintel       int    The minimum accetible number of telescopes in
 *                      closed sub-arrays. Normally this would be 3 for
 *                      phase self-cal and 4 for amplitude self-cal.
 *  doflag       int    Flag un-correctable baselines if true (!=0).
 * Input/Output:
 *  usable       int *  The caller must supply an un-initialized array of
 *                      sub->nbase int's (one element per visibility in the
 *                      current integration). On return usable
 *                      visibilities will be flagged as true (!=0) in
 *                      this array.
 *  telnum       int *  The caller must supply an un-initialized array
 *                      of nstat int's (One element per telescope).
 *                      On return this will contain a count of how
 *                      many connected baselines each telescope
 *                      belongs to.
 *  nbadtel      int *  If not sent as NULL, then the variable pointed to
 *                      by this argument will be assigned its current
 *                      value + the number of telescopes that are too
 *                      undersampled to be corrected.
 * Output:
 *  return       int    0 - Self-calibration may proceed if restricted
 *                          to the usable visibilities marked in usable[].
 *                      1 - Insufficient data for self-cal.
 */
static int get_usable(Observation *ob, Subarray *sub, int ut,
		      float uvmin, float uvmax, int mintel, int doflag,
		      int *usable, int *telnum, int *nbadtel)
{
  Integration *integ; /* The descriptor of the integration being checked */
  Telcor *tcor;       /* An existing telescope correction */
  int ntel;           /* The number of constrained telescopes */
  int itel;           /* A telescope index */
/*
 * Sanity checks.
 */
  if(sub_bad(sub, "get_usable"))
    return 1;
  if(usable==NULL || telnum==NULL) {
    lprintf(stderr, "get_usable: NULL work array(s) intercepted\n");
    return 1;
  };
/*
 * Get the descriptor of the specified integration.
 */
  integ = &sub->integ[ut];
/*
 * Ascertain which visibilities are usable, based solely on existing flags
 * and the given UV range.
 */
  visflags(ob, integ->vis, sub->nbase, uvmin, uvmax, usable);
/*
 * Now exclude all telescopes and their visibilities that lie on
 * baselines which are not part of closed sub-arrays
 * and count the number telescopes left in play.
 */
  ntel = count_tel(ob, sub, ut, mintel>2, doflag, usable, telnum);
/*
 * If there are too few telescopes arrange to ignore all baselines and
 * telescopes that are still marked as usable.
 */
  if(ntel < mintel) {
    int base;
    for(base=0; base<sub->nbase; base++)
      usable[base] = 0;
    for(itel=0; itel<sub->nstat; itel++)
      telnum[itel] = 0;
  };
/*
 * Count and optionally flag dropped telescope corrections, except where
 * corrections have not been requested.
 */
  tcor = &integ->icor[ob->stream.cif].tcor[0];
  for(itel=0; itel<sub->nstat; itel++) {
    if(telnum[itel]==0 && !sub->tel[itel].antfix) {
      if(doflag && !tcor[itel].bad) {
	(*nbadtel)++;
	ed_Telcor(ob, sub, ob->stream.cif, ut, itel, 1);
      };
    };
  };
/*
 * Return an indication of whether there is sufficient data for self-cal
 * to proceed.
 */
  return ntel < mintel ? 1 : 0;
}

/*.......................................................................
 * Where there is insufficient data to solve for the phase and gain
 * errors at a given telescope mark the baselines involving that
 * telescope as un-usable in usable[].
 *
 * Input:
 *  ob   Observation *  The descriptor of the parent observation.
 *  sub     Subarray *  The descriptor of the sub-array being corrected.
 *  ut           int    The index of the integration to be processed.
 *  doclose      int    If true, then eliminate visibilities that do
 *                      not lie on closed arrays of telescopes from the
 *                      usable[] array before counting used telescopes.
 *  doflag       int    If both doflag and doclose are true, flag
 *                      visibilities that are marked as un-usable by
 *                      this function.
 * Input/Output:
 *  usable       int *  The caller must supply an un-initialized array of
 *                      sub->nbase int's (one element per visibility in the
 *                      current integration). On return usable
 *                      visibilities will be flagged as true (!=0) in
 *                      this array.
 *  telnum       int *  The caller must supply an un-initialized array
 *                      of nstat int's (One element per telescope).
 *                      On return this will contain a count of how
 *                      many connected baselines each telescope
 *                      belongs to.
 * Output:
 *  return       int    The number of telescopes for which solutions can
 *                      be found.
 */
static int count_tel(Observation *ob, Subarray *sub, int ut, int doclose,
		     int doflag, int *usable, int *telnum)
{
  Baseline *bptr; /* Pointer into sub->base[] */
  int *uptr;      /* Pointer into the usable[] array */
  int base;       /* The index of the baseline being checked */
  int itel;       /* Checks each telescope in turn */
  int ntel;       /* The number of telescopes that can be solved for */
/*
 * Count the number of times each telescope is referenced on a
 * usable baseline.
 */
  for(itel=0;itel<sub->nstat;itel++)
    telnum[itel]=0;
  bptr = sub->base;
  uptr = usable;
  for(base=0; base<sub->nbase; base++,bptr++,uptr++) {
    if(*uptr) {
      telnum[bptr->tel_a]++;
      telnum[bptr->tel_b]++;
    };
  };
/*
 * Find telescopes that are only cited on one baseline and check
 * the telescope at the other end of the baseline...
 */
  if(doclose) {
    for(itel=0; itel<sub->nstat; itel++) {
      int newtel=itel;
/*
 * If the current telescope appears on just one baseline locate and
 * flag that baseline as unusable and see if the removal of that baseline
 * makes the other telescope on that baseline insoluble.
 */
      while(telnum[newtel]==1) {
/*
 * Search baselines affected by the removal of 'newtel'.
 */
	bptr = sub->base;
	uptr = usable;
	for(base=0; base<sub->nbase; base++,bptr++,uptr++) {
	  if(*uptr) {
	    int ita = bptr->tel_a;
	    int itb = bptr->tel_b;
/*
 * Locate the lone baseline.
 */
	    if(ita==newtel || itb==newtel) {
	      *uptr=0;
	      telnum[ita]--;
	      telnum[itb]--;
/*
 * Setup to check the other telescope on this baseline.
 */
	      newtel = (ita==newtel)?itb:ita;
	      break;
	    };
	  };
	};
      };
    };
  };
/*
 * Count the number of telescopes for which solutions can be found.
 */
  for(ntel=itel=0; itel<sub->nstat; itel++) {
    if(telnum[itel])
      ntel++;
  };
  return ntel;
}

/*.......................................................................
 * Normalize the established gain corrections to prevent the flux scale
 * of the sub-array from wandering over repeated iterations of
 * self-cal and CLEAN.
 *
 * Input:
 *  sub  Subarray *  The descriptor of the sub-array to be corrected.
 *  cif       int    The index of the IF for which the corrections are
 *                   to be recorded.
 *  cors      Cor *  A work array to compile corrections in.
 * Output:
 *  return  float    The telescope normalization factor applied.
 */
static float norm_cors(Subarray *sub, int cif, Cor *cors)
{
  Integration *integ;  /* The descriptor of an integration */
  Station *tel;        /* The descriptor of telescope 'itel' */
  double amp_sum=0.0;  /* The sum of amplitude corrections */
  Cor *icor;           /* Pointer into input correction array, cors[] */
  Telcor *ocor;        /* Pointer into output correction array */
  float amp_cor;       /* The mean amplitude correction */
  int namp=0;          /* The number of corrections in the sum */
  int itel;            /* Index of telescope correction being processed */
  int ut;              /* Index of integration being processed */
/*
 * Accumulate the sum of all non-fixed amplitude corrections.
 */
  integ = sub->integ;
  for(ut=0; ut<sub->ntime; ut++,integ++) {
    ocor = integ->icor[cif].tcor;
    tel = sub->tel;
    for(itel=0; itel<sub->nstat; itel++,ocor++,tel++) {
      if(ocor->amp_cor > 0.0f && !tel->antfix) {
	amp_sum += ocor->amp_cor;
	namp++;
      };
    };
  };
/*
 * No corrections in sum?
 */
  if(namp<1)
    return 1.0f;
/*
 * Turn the sum into the reciprocal of the mean amplitude correction.
 */
  amp_cor = namp/amp_sum;
/*
 * Apply this as a correction to all output corrections and data, but
 * take care not to mark un-corrected corrections as corrected, by
 * making the input correction weights agree with the current output
 * correction statuses.
 */
  integ = sub->integ;
  for(ut=0; ut<sub->ntime; ut++,integ++) {
/*
 * Compile the array of corrections to be applied - don't normalize
 * the corrections of fixed antennas (if any).
 */
    ocor = integ->icor[cif].tcor;
    icor = cors;
    tel = sub->tel;
    for(itel=0; itel<sub->nstat; itel++,icor++,ocor++,tel++) {
      icor->amp_cor = tel->antfix ? 1.0f : amp_cor;
      icor->weight = ocor->amp_cor > 0.0f ? 1.0f : 0.0f;
    };
/*
 * Apply the corrections.
 */
    apply_cors(sub, cif, ut, ut, 1, 0, cors);
  };
  return amp_cor;
}

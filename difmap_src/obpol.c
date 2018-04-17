#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "logio.h"
#include "obs.h"

static GETPOL_FN(get_pol);
static GETPOL_FN(get_ipol);
static GETPOL_FN(get_qpol);
static GETPOL_FN(get_upol);
static GETPOL_FN(get_vpol);
static GETPOL_FN(get_pi_pol);

static int find_stokes(Observation *ob, Stokes pol);

/*.......................................................................
 * Search for a way to create a given polarization from the stokes
 * parameters recorded in an observation. Record the results in a
 * given container along with a function that can use the contents of
 * the container to extract visibilities of the specified polarization
 * from the array ob->npol polarized visibilities contained in the Dbase 
 * level of an ob->dp I/O buffer.
 *
 * Input:
 *  ob   Observation *  The observation containing the polarized
 *                      visibilities.
 *  stokes    Stokes    The type of polarization required.
 *  report       int    0 - Don't report failed searches.
 *                      1 - Report if the requested polarization is
 *                          unavailable.
 * Input/Output:
 *  obpol      Obpol *  The container to be filled, or NULL if this
 *                      function is simply being used to check the validity
 *                      of a polarization. obpol will not be changed unless
 *                      a valid polarization is found.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Polarization not found.
 */
int get_Obpol(Observation *ob, Stokes stokes, int report, Obpol *obpol)
{
  Obpol ptmp;  /* Work container to be copied when complete */
  int found=0; /* True if the requested polarization is acquired */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "get_Obpol"))
    return 1;
/*
 * Initialize the work container.
 */
  ptmp.type = stokes;
  ptmp.pa = ptmp.pb = -1;
  ptmp.getpol = 0;
/*
 * Is the requested polarization directly recorded in the data?
 */
  ptmp.pa = find_stokes(ob, stokes);
/*
 * If the polarization was not recorded with the observation, see if
 * it can be derived from two other polarizations that were observed.
 */
  if(ptmp.pa >= 0) {
    ptmp.getpol = get_pol;
    found = 1;
  } else if(stokes == PI_POL) {
    ptmp.pa = find_stokes(ob, RR);
    ptmp.pb = find_stokes(ob, LL);
/*
 * If only one of RR or LL is available, make sure that pa points to it.
 */
    if(ptmp.pa < 0 && ptmp.pb >= 0) {
      ptmp.pa = ptmp.pb;
      ptmp.pb = -1;
    };
    ptmp.getpol = get_pi_pol;
/*
 * Was at least one or RR or LL found?
 */
    found = ptmp.pa >= 0;
  } else {
    switch(stokes) {
    case SI:
      ptmp.pa = find_stokes(ob, RR);
      ptmp.pb = find_stokes(ob, LL);
      ptmp.getpol = get_ipol;
      break;
    case SV:
      ptmp.pa = find_stokes(ob, RR);
      ptmp.pb = find_stokes(ob, LL);
      ptmp.getpol = get_vpol;
      break;
    case SQ:
      ptmp.pa = find_stokes(ob, RL);
      ptmp.pb = find_stokes(ob, LR);
      ptmp.getpol = get_qpol;
      break;
    case SU:
      ptmp.pa = find_stokes(ob, LR);
      ptmp.pb = find_stokes(ob, RL);
      ptmp.getpol = get_upol;
      break;
    case NO_POL:   /* Substitute default polarization */
      if(ob_ready(ob, OB_SELECT, NULL) && ob->stream.pol.type != NO_POL) {
	return get_Obpol(ob, ob->stream.pol.type, 0, obpol);
      } else if(get_Obpol(ob, SI, 0, obpol)==0) {
	return 0;
      } else if(ob->pols[0] != NO_POL) {
	return get_Obpol(ob, ob->pols[0], 1, obpol);
      };
      break;
    default:
      break;
    };
/*
 * Was the combined polarization available?
 */
    found = ptmp.pa >= 0 && ptmp.pb >= 0;
  };
/*
 * No such polarization?
 */
  if(!found) {
    if(report)
      lprintf(stderr, "Polarization %s is unavailable.\n", Stokes_name(stokes));
    return 1;
  };
/*
 * Return a copy of the work container.
 */
  if(obpol)
    *obpol = ptmp;
  return 0;
}

/*.......................................................................
 * The function used to extract a visibility of a directly recorded
 * polarization from an ob->dp->ifs[].chan[].base[].pol array of ob->npol
 * visibilities.
 *
 * Input:
 *  pol   Obpol *  The polarization descriptor containing the index of
 *                 the polarization to be extracted.
 *  pvis   Cvis *  A ob->dp->ifs[].chan[].base[].pol array of ob->npol
 *                 polarized visibilities to extract the required
 *                 visibility from.
 * Input/Output:
 *  out    Cvis *  The container to hold the extracted visibility of
 *                 the selected polarization.
 */
static GETPOL_FN(get_pol)
{
  *out = pvis[pol->pa];
  return;
}

/*.......................................................................
 * The function used to combine recorded RR and LL visibilities in an
 * ob->dp->ifs[].chan[].base[].pol array to extract a visibility of
 * stokes I = (RR+LL)/2.
 *
 * Input:
 *  pol   Obpol *  The polarization descriptor containing the index of
 *                 the polarization to be extracted.
 *  pvis   Cvis *  A ob->dp->ifs[].chan[].base[].pol array of ob->npol
 *                 polarized visibilities to extract the required
 *                 visibility from.
 * Input/Output:
 *  out    Cvis *  The container to hold the extracted visibility of
 *                 the selected polarization.
 */
static GETPOL_FN(get_ipol)
{
  Cvis *avis = pvis + pol->pa;  /* RR */
  Cvis *bvis = pvis + pol->pb;  /* LL */
/*
 * If either visibility is deleted, then the combined visibility is
 * also deleted and its value is inconsequential.
 */
  if(avis->wt == 0.0f || bvis->wt == 0.0f) {
    out->re = out->im = out->wt = 0.0f;
  } else {
/*
 * Get polarization I = (RR+LL)/2.
 */
    out->re = 0.5 * (avis->re + bvis->re);
    out->im = 0.5 * (avis->im + bvis->im);
/*
 * Determine the combined weight.
 */
    out->wt = 4.0f/(1.0f/fabs(avis->wt) + 1.0f/fabs(bvis->wt));
    if(avis->wt < 0.0f || bvis->wt < 0.0f)
      out->wt = -out->wt;
  };
  return;
}

/*.......................................................................
 * The function used to combine recorded RL and LR visibilities in an
 * ob->dp->ifs[].chan[].base[].pol array to extract a visibility of
 * stokes Q = (RL+LR)/2.
 *
 * Input:
 *  pol   Obpol *  The polarization descriptor containing the index of
 *                 the polarization to be extracted.
 *  pvis   Cvis *  A ob->dp->ifs[].chan[].base[].pol array of ob->npol
 *                 polarized visibilities to extract the required
 *                 visibility from.
 * Input/Output:
 *  out    Cvis *  The container to hold the extracted visibility of
 *                 the selected polarization.
 */
static GETPOL_FN(get_qpol)
{
  Cvis *avis = pvis + pol->pa;  /* RL */
  Cvis *bvis = pvis + pol->pb;  /* LR */
/*
 * If either visibility is deleted, then the combined visibility is
 * also deleted and its value is inconsequential.
 */
  if(avis->wt == 0.0f || bvis->wt == 0.0f) {
    out->re = out->im = out->wt = 0.0f;
  } else {
/*
 * Get polarization Q = (RL+LR)/2.
 */
    out->re = 0.5 * (avis->re + bvis->re);
    out->im = 0.5 * (avis->im + bvis->im);
/*
 * Determine the combined weight.
 */
    out->wt = 4.0f/(1.0f/fabs(avis->wt) + 1.0f/fabs(bvis->wt));
    if(avis->wt < 0.0f || bvis->wt < 0.0f)
      out->wt = -out->wt;
  };
  return;
}

/*.......................................................................
 * The function used to combine recorded LR and RL visibilities in an
 * ob->dp->ifs[].chan[].base[].pol array to extract a visibility of
 * stokes U = i(LR-RL)/2.
 *
 * Input:
 *  pol   Obpol *  The polarization descriptor containing the index of
 *                 the polarization to be extracted.
 *  pvis   Cvis *  A ob->dp->ifs[].chan[].base[].pol array of ob->npol
 *                 polarized visibilities to extract the required
 *                 visibility from.
 * Input/Output:
 *  out    Cvis *  The container to hold the extracted visibility of
 *                 the selected polarization.
 */
static GETPOL_FN(get_upol)
{
  Cvis *avis = pvis + pol->pa;  /* LR */
  Cvis *bvis = pvis + pol->pb;  /* RL */
/*
 * If either visibility is deleted, then the combined visibility is
 * also deleted and its value is inconsequential.
 */
  if(avis->wt == 0.0f || bvis->wt == 0.0f) {
    out->re = out->im = out->wt = 0.0f;
  } else {
/*
 * Get polarization U = i(LR-RL)/2.
 */
    out->re = -0.5 * (avis->im - bvis->im);
    out->im =  0.5 * (avis->re - bvis->re);
/*
 * Determine the combined weight.
 */
    out->wt = 4.0f/(1.0f/fabs(avis->wt) + 1.0f/fabs(bvis->wt));
    if(avis->wt < 0.0f || bvis->wt < 0.0f)
      out->wt = -out->wt;
  };
  return;
}

/*.......................................................................
 * The function used to combine recorded RR and LL visibilities in an
 * ob->dp->ifs[].chan[].base[].pol array to extract a visibility of
 * stokes V = (RR-LL)/2.
 *
 * Input:
 *  pol   Obpol *  The polarization descriptor containing the index of
 *                 the polarization to be extracted.
 *  pvis   Cvis *  A ob->dp->ifs[].chan[].base[].pol array of ob->npol
 *                 polarized visibilities to extract the required
 *                 visibility from.
 * Input/Output:
 *  out    Cvis *  The container to hold the extracted visibility of
 *                 the selected polarization.
 */
static GETPOL_FN(get_vpol)
{
  Cvis *avis = pvis + pol->pa;  /* RR */
  Cvis *bvis = pvis + pol->pb;  /* LL */
/*
 * If either visibility is deleted, then the combined visibility is
 * also deleted and its value is inconsequential.
 */
  if(avis->wt == 0.0f || bvis->wt == 0.0f) {
    out->re = out->im = out->wt = 0.0f;
  } else {
/*
 * Get polarization V = (RR-LL)/2.
 */
    out->re = 0.5 * (avis->re - bvis->re);
    out->im = 0.5 * (avis->im - bvis->im);
/*
 * Determine the combined weight.
 */
    out->wt = 4.0f/(1.0f/fabs(avis->wt) + 1.0f/fabs(bvis->wt));
    if(avis->wt < 0.0f || bvis->wt < 0.0f)
      out->wt = -out->wt;
  };
  return;
}

/*.......................................................................
 * The function used to combine recorded RR and LL visibilities in an
 * ob->dp->ifs[].chan[].base[].pol array to extract a visibility of
 * pseudo-I polarization, assuming that RR and LL are both samples of I.
 *
 * Input:
 *  pol   Obpol *  The polarization descriptor containing the index of
 *                 the polarization to be extracted.
 *  pvis   Cvis *  A ob->dp->ifs[].chan[].base[].pol array of ob->npol
 *                 polarized visibilities to extract the required
 *                 visibility from.
 * Input/Output:
 *  out    Cvis *  The container to hold the extracted visibility of
 *                 the selected polarization.
 */
static GETPOL_FN(get_pi_pol)
{
/*
 * If only one of RR and LL was found, the single visibility to be
 * returned will be pointed to by pa, and pb will be -1.
 */
  if(pol->pb < 0) {
    *out = pvis[pol->pa];
  } else {
    Cvis *avis = pvis + pol->pa;   /* RR */
    Cvis *bvis = pvis + pol->pb;   /* LL */
/*
 * If the RR and LL visibilities are either both flagged or both unflagged,
 * merge them to produce a weighted sum visibility with the corresponding
 * flag status.
 */
    if((avis->wt > 0.0 && bvis->wt > 0.0) ||
       (avis->wt < 0.0 && bvis->wt < 0.0)) {
      float aw = fabs(avis->wt);
      float bw = fabs(bvis->wt);
      out->re = (avis->re * aw + bvis->re * bw) / (aw + bw);
      out->im = (avis->im * aw + bvis->im * bw) / (aw + bw);
      out->wt = avis->wt + bvis->wt;
/*
 * If only one of the visibilities is unflagged, assign it as the output
 * visibility.
 */
    } else if(avis->wt > 0.0) {  /* avis->wt > 0.0 && bvis->wt <= 0.0 */
      *out = *avis;
    } else if(bvis->wt > 0.0) {  /* avis->wt <= 0.0 && bvis->wt > 0.0 */
      *out = *bvis;
/*
 * If both visibilities are deleted, then the output visibility is also
 * deleted?
 */
    } else {                     /* avis->wt == 0.0 && bvis->wt == 0.0 */
      out->re = out->im = out->wt = 0.0;
    };
  };
  return;
}

/*.......................................................................
 * Locate a given STOKES parameter in the observed data. Return its
 * index.
 *
 * Input:
 *  ob  Observation *  The descriptor of the Observation.
 *  pol      Stokes    The stokes parameter to find.
 * Output:
 *  return      int    The index of the parameter in ob->pols[], or -1 if
 *                     not found.
 */
static int find_stokes(Observation *ob, Stokes pol)
{
  int i;
/*
 * Search for the requested polarization in the list of observed
 * polarizations.
 */
  for(i=0; i<ob->npol; i++) {
    if(ob->pols[i] == pol)
      return i;
  };
/*
 * Polarization not found.
 */
  return -1;
}

#include <stdio.h>
#include <stdlib.h>

#include "obs.h"
#include "logio.h"

/*.......................................................................
 * Iterate through the available indexes of IFs in a given observation.
 *
 * For example to iterate through the indexes of all sampled IFs between
 * bif and eif (inclusive), one would type:
 *
 * for(cif=bif; (cif=nextIF(ob, cif, 1, 1)) >= 0 && cif<=eif; cif++) {
 *   if(getIF(ob, cif))
 *     return 1;
 *   ...
 * };
 *
 * Input:
 *  ob    Observation *   The descriptor of the observation.
 *  cif           int     The start IF index - if 0>cif>=ob->nif then
 *                        -1 will be returned.
 *  skip_empty    int     If true don't return the indexes of
 *                        unsampled IFs (those with no selected
 *                        channels).
 *  step          int     The iteration step size, or 0 for no iteration.
 *                        If cif is not suitable, then it will be repeatedly
 *                        incremented by 'step' until a valid index is
 *                        found, or until the range of IF indexes has
 *                        been exhausted. Usual values are:
 *                         0 - Only cif is to be checked.
 *                         1 - Starting with cif, search forward for
 *                             the next suitable index.
 *                        -1 - Starting with cif, search backward for
 *                             the next suitable index.
 * Output:
 *  return        int     The index of the next IF to visit, or -1 if
 *                        all indexes have been exhausted.
 */
int nextIF(Observation *ob, int cif, int skip_empty, int step)
{
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "nextIF"))
    return -1;
/*
 * Out of range?
 */
  if(cif < 0 || cif >= ob->nif)
    return -1;
/*
 * Locate the next suitable IF.
 */
  if(step < 0) {
    while(cif>=0 && (skip_empty && ob->ifs[cif].cl == NULL))
      cif--;
    if(cif < 0)
      cif = -1;
  } else if(step > 0) {
    while(cif<ob->nif && (skip_empty && ob->ifs[cif].cl == NULL))
      cif++;
    if(cif >= ob->nif)
      cif = -1;
  } else {
    if(skip_empty && ob->ifs[cif].cl == NULL)
      cif = -1;
  };
/*
 * Return the index of the new IF.
 */
  return cif;
}

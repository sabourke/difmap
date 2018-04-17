#include <stdio.h>
#include <stdlib.h>

#include "logio.h"
#include "obs.h"

/*.......................................................................
 * Create or resize the array of IF descriptors in an Observation.
 *
 * Input:
 *  ob  Observation * The descriptor of the observation.
 *                    ob->ifs must be NULL unless previously allocated.
 *  nif         int   The number of IFs in the observation.
 * Output:
 *  return       If * The pointer to the array of IFs, or NULL on error.
 */
If *new_If(Observation *ob, int nif)
{
  int cif;  /* Index of the IF descriptor being processed */
  If *ifp;  /* The pointer into the array of IF descriptors */
/*
 * Valid Observation descriptor?
 */
  if(!ob_ready(ob, OB_ALLOC, "new_If"))
    return NULL;
/*
 * Is this the first time that IF descriptors have been allocated?
 */
  if(ob->ifs==NULL || ob->nif==0) {
    ob->nif = 0;
    ob->ifs = (If *) calloc((size_t) nif, sizeof(If));
    if(ob->ifs==NULL) {
      lprintf(stderr, "new_If: Insufficient memory for new IF array.\n");
      return NULL;
    };
/*
 * Resize the existing array of IF descriptors.
 */
  } else {
/*
 * Resize the array of IF descriptors if necessary.
 */
    if(nif != ob->nif) {
      ifp = realloc(ob->ifs, nif * sizeof(If));
      if(ifp)
	ob->ifs = ifp;
      else if(nif > ob->nif) {
	lprintf(stderr, "new_If: Insufficient memory to re-size IF array.\n");
	return NULL;
      };
    };
  };
/*
 * Default initialize any the new IF descriptors.
 */
  if(nif > ob->nif) {
    ifp = &ob->ifs[ob->nif];
    for(cif=ob->nif;  cif<nif;  cif++,ifp++) {
      ifp->freq = ifp->df = ifp->bw = 0.0;
      ifp->coff = 0;
      ifp->cl = NULL;
      ifp->wtsum_bad = 1;
    };
  };
/*
 * Record the new count of IFs.
 */
  ob->nif = nif;
  return ob->ifs;
}

/*.......................................................................
 * Delete the array of IF descriptors in an Observation.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation.
 * Output:
 *  return        If *  Allways NULL.
 */
If *del_If(Observation *ob)
{
/*
 * Anything to delete?
 */
  if(ob && ob->ifs) {
    If *ifp;
/*
 * Delete stream channel lists.
 */
    for(ifp=ob->ifs; ifp<ob->ifs + ob->nif; ifp++)
      ifp->cl = del_Chlist(ifp->cl);
/*
 * Free the IF container.
 */
    free(ob->ifs);
    ob->ifs = NULL;
  };
/*
 * Return the deleted IF array.
 */
  return NULL;
}

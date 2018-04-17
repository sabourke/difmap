#include <stdio.h>
#include <stdlib.h>
#include "logio.h"
#include "vlbconst.h"
#include "mapwin.h"
#include "model.h"
#include "winmod.h"

/*.......................................................................
 * Split a model up into two models containing, respectively, the
 * components that lie inside at least one window and the components
 * that lie outside of all windows. A new model containing the components
 * that lay within the windows will be returned, while the input model,
 * 'mod' will be left with the components that lay outside the windows.
 * If you only require the components within or those outside the windows
 * then don't forget to delete the unwanted model.
 *
 * Input:
 *   mod    Model *  The Model to be split. On return this contains the
 *                   model components that lie within the windows.
 *   wins  Mapwin *  The container of the window list.
 *   docomp   int    If true, merge appended delta components where possible
 *                   with existing delta components at the same positions,
 *                   in the returned model.
 * Output:
 *   return Model *  A new model containing the components that lay
 *                   within at least one window. NULL if memory allocation
 *                   failed.
 */
Model *win_mod(Model *mod, Mapwin *wins, int docomp)
{
  Modcmp *cmp;   /* The model component being considered. */
  Modcmp *next;  /* The next model component in the list after 'cmp' */
  Modcmp *last;  /* The last undeleted component in the list before 'cmp' */
  Model *retmod; /* The return model */
/*
 * Allocate a model for the components within the windows.
 */
  retmod = new_Model();
  if(retmod==NULL) {
    lprintf(stderr, "win_mod: Insufficient memory for temporary model\n");
    return retmod;
  };
/*
 * Loop through the input model components, remove components that lie
 * inside the windows and append them to 'retmod'. 'last' is the last 
 * undeleted component before the current one, 'cmp' and 'next' is the next
 * one after it in the list. 
 */
  last = NULL;
  for(cmp=mod->head; cmp != NULL; cmp=next) {
    next = cmp->next;
    if(!inmapwin(wins, cmp->x, cmp->y))
      last=cmp;
    else      /* Move components that lie outside all windows to 'retmod' */
      cmp=add_cmp(rem_cmp(mod,last,cmp), retmod, docomp);
  };
  return retmod;
}

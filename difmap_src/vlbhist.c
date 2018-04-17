#include <stdlib.h>
#include <stdio.h>

#include "obs.h"

#include "sphere.h"
#include "pager.h"
#include "logio.h"

/*.......................................................................
 * Display the history information of the current observation through
 * a pager.
 *
 * Input:
 *  ob   Observation *  The observation descriptor.
 *  dopage       int    If true, display through a pager. If 0, write
 *                      directly to stdout.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
int showhist(Observation *ob, int dopage)
{
  char history[81];  /* Buffer to read one history line into */
  Pager *page;       /* The descriptor of the file pager */
  int ierr=0;        /* Error status */
  int i;
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "showhist"))
    return 1;
/*
 * No history to be displayed?
 */
  if(ob->nhist<1) {
    lprintf(stdout, "There are no history lines in this observation.\n");
    return 0;
  };
/*
 * Start a file pager.
 */
  page = new_Pager();
  if(page) {
/*
 * Rewind the history file.
 */
    rec_rewind(ob->his);
/*
 * Send each history line to the pager.
 */
    for(i=0; !ierr && i<ob->nhist; i++) {
      ierr = rec_read(ob->his, 80, 1, history) < 0;
      history[80]='\0';
      ierr = ierr || pprintf(page, "%s\n", history) < 0;
    };
/*
 * Have the pager show the result.
 */
    page = end_Pager(page, !ierr, pause_output, dopage ? PAGE_INT : PAGE_OFF);
  } else {
    ierr = 1;
  };
  return ierr;
}

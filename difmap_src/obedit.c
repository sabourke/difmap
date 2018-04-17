#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "logio.h"
#include "obs.h"
#include "obedit.h"
#include "telspec.h"
#include "vlbconst.h"

/*
 * The maximum number of edit blocks to allow before flushing edits
 * See obedit.h to see the number of edit nodes in an edit block (EDBLK).
 */

enum {MAX_EDBLK=10};

static Edint *new_Edint(Observation *ob);
static Edint *del_Edint(Obedit *obed, Edint *ed);

typedef struct {
  int cmin,cmax; /* Indexes of first and last spectral-line channels. */
  int bmin,bmax; /* Indexes of first and last baselines. */
  int nedit;     /* Number of contributing edits counted */
} Edrange;

static Edrange *ed_range(Observation *ob, Integration *integ, int cif);
static int ed_uvdata(Observation *ob);
static int ed_ifdata(Observation *ob);
static int dp_edit(Observation *ob, Integration *integ, int cif);
static int ip_edit(Observation *ob, Integration *integ, int cif);

/*.......................................................................
 * Allocate a container for a list of Edint blocks to be aportioned
 * into a free-list for allocation of Edint nodes. Install the container
 * in its associated Observation descriptor.
 *
 * Input:
 *  ob  Observation *  THe descriptor of the observation for which the
 *                     edit descriptor is being allocated.
 * Output:
 *  return   Obedit *  The pointer to the initialized container, or NULL
 *                     on error.
 */
Obedit *new_Obedit(Observation *ob)
{
  Obedit *obed;  /* The pointert to the new container */
/*
 * Allocate the container.
 */
  ob->obed = obed = malloc(sizeof(Obedit));
  if(obed==NULL) {
    lprintf(stderr, "new_Obedit: Insuficient memory for edit list.\n");
    return NULL;
  };
/*
 * Initialize at least to the point at which del_Obedit() can be called.
 */
  obed->nused = 0;
  obed->free = NULL;
  obed->blocks.next = NULL;
/*
 * Initialize the free-list etc..
 */
  if(clr_Obedit(ob))
    return del_Obedit(ob);
/*
 * Return the initialized container.
 */
  return obed;
}

/*.......................................................................
 * Release the memory assigned to an Obedit container and its members.
 *
 * Input:
 *  ob Observation *  The descriptor of the Observation containing the
 *                    Obedit container to be deleted.
 * Output:
 *  return  Obedit *  Allways NULL. 
 */
Obedit *del_Obedit(Observation *ob)
{
  if(ob && ob->obed) {
    Obedit *obed = ob->obed;
/*
 * Delete any extra Edint blocks.
 */
    if(obed->blocks.next) {
      struct edblock *next;   /* The next block to be deleted */
      struct edblock *last;   /* The block being deleted */
      next = obed->blocks.next;
      while(next) {
	last = next;
	next = last->next;
	free(last);
      };
    };
/*
 * Delete the container.
 */
    free(obed);
  };
  return NULL;
}

/*.......................................................................
 * Return the Obedit structure associated with a given Observation, to
 * its original state and delete the list of edits associated with
 * each integration, ready for a new edit session.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation containing
 *                     the list of edits.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int clr_Obedit(Observation *ob)
{
  Obedit *obed; /* Local copy of ob->obed */
  int i;
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_ALLOC, "clr_Obedit"))
    return 1;
/*
 * Get the pointer to the Obedit container.
 */
  obed = ob->obed;
/*
 * Is there an Obedit container to be cleared?
 */
  if(obed) {
    Intrec *rec; /* Pointer into ob->rec[] */
    Edint *head;       /* Pointer to head of free-list */
    int ut;            /* The index of an integration */
/*
 * Delete all edit lists.
 */
    rec = ob->rec;
    for(ut=0; ut<ob->nrec; ut++,rec++)
      rec->integ->edlist = NULL;
/*
 * Delete all but the first block of the free-list.
 */
    if(obed->blocks.next) {
      struct edblock *next;   /* The next block to be deleted */
      struct edblock *last;   /* The block being deleted */
      next = obed->blocks.next;
      while(next) {
	last = next;
	next = last->next;
	free(last);
      };
    };
/*
 * Initialize the block list.
 */
    obed->nused = 0;
    obed->free = NULL;
    obed->blocks.next = NULL;
/*
 * Re-link the free-list.
 */    
    head = &obed->blocks.block[0];
    for(i=0; i<EDBLK-1; i++)
      head[i].next = &head[i+1];
    head[EDBLK-1].next = NULL;
/*
 * Make the start of the block the head of the free-list.
 */
    obed->free = head;
  };
  return 0;
}

/*.......................................................................
 * Get a new Edint node from the free-list.
 *
 * Input:
 *  ob  Observation * The observation in which the edit list resides.
 * Output:
 *  return    Edint * The pointer to the initialized Edint descriptor, or
 *                    NULL on error.
 */
static Edint *new_Edint(Observation *ob)
{
  Obedit *obed; /* The container of the free-list */
  Edint *ed;    /* The new node to be returned */
  int i;
/*
 * Get the edit-list container.
 */
  obed = ob->obed;
/*
 * Do we need to allocate a new block of Edint nodes?
 */
  if(obed->free==NULL) {
    struct edblock *tail = &obed->blocks;
/*
 * Flush pending edits if the limit on the number of cached edits
 * has been reached.
 */
    if(obed->nused >= EDBLK * MAX_EDBLK && ed_flush(ob))
      return NULL;
/*
 * Locate the tail of the list of Edint blocks.
 */
    while(tail->next)
      tail = tail->next;
/*
 * Allocate a new block.
 */
    tail->next = malloc(sizeof(struct edblock));
/*
 * No more memory?
 */
    if(tail->next==NULL)
      return NULL;
/*
 * Terminate the list of edblock's.
 */
    tail->next->next = NULL;
/*
 * The head of the free list is the first element of the new block.
 */
    obed->free = &tail->next->block[0];
/*
 * Link the free-list.
 */
    ed = obed->free;
    for(i=0; i<EDBLK-1; i++)
      ed[i].next = &ed[i+1];
    ed[EDBLK-1].next = NULL;
  };
/*
 * Extract the new node from the head of the free-list.
 */
  ed = obed->free;
  obed->free = ed->next;
  ed->next = NULL;
  obed->nused++;
/*
 * Return the new un-initialized descriptor.
 */
  return ed;
}

/*.......................................................................
 * Return a given edit node to the free-list.
 *
 * Input:
 *  obed    Obedit *  The container of the free-list.
 *  ed       Edint *  The edit node to be discarded.
 * Output:
 *  return   Edint *  Allways NULL. Use like  ed=del_Edint(ed);
 */
static Edint *del_Edint(Obedit *obed, Edint *ed)
{
/*
 * Nothing to be deleted?
 */
  if(ed==NULL)
    return ed;
/*
 * No free-list to return to?
 */
  if(obed==NULL) {
    lprintf(stderr, "del_Edint: NULL Obedit descriptor intercepted.\n");
    return NULL;
  };
/*
 * Insert the discarded node at the head of the free-list.
 */
  ed->next = obed->free;
  obed->free = ed;
  obed->nused--;
  return NULL;
}

/*.......................................................................
 * Apply a given edit operation to a given integration of a specified
 * sub-array, and if there are data scratch files, defer applying the
 * edits to them by allocating, initializing and appending a new Edint
 * node to the list of deferred edits within a given integration.
 *
 * When an editing session is completed ed_flush() must be called to
 * have the deferred edit operations applied to the data in the scratch
 * files. 
 *
 * Input:
 *  ob Observation * The descriptor of the observation ebing edited.
 *  sub   Subarray * The descriptor of the sub-array to be edited.
 *  ut         int   The index of the integration in the sub-array
 *                   to be edited.
 *  cif        int   The index of the IF to be edited if 'selif' is true.
 *                   If 'selif' is false, the value of cif is irrelevant.
 *  doflag     int   The type of operation to be recorded. If true flag
 *                   data; else restore data.
 *  selbase    int   If true record that only visibilities of the baseline
 *                   'index' be edited.
 *  selstat    int   If true record that only visibilities of the station
 *                   'index' be edited (ignored if selbase is true).
 *  selchan    int   If true record that only visibilities of the
 *                   spectral-line channels from which the current stream
 *                   was formed be edited. Otherwise the record will
 *                   specify that all channels be edited.
 *  selif      int   If true record that only visibilities of the IF refered
 *                   to by the index stored in 'cif' be edited. Otherwise
 *                   the record will specify that all IFs be edited.
 *  index      int   The index of a baseline if selbase is true, or
 *                   the index of a telescope if selstat is true.
 *                   Otherwise it is ignored.
 * Output:
 *  return     int   0 - OK.
 *                   1 - Error.
 */
int ed_integ(Observation *ob, Subarray *sub, int ut, int cif, int doflag,
	     int selbase, int selstat, int selchan, int selif, int index)
{
  Integration *integ;  /* The descriptor of the integration being edited */
/*
 * Check args - note that edits refer to the current channel range
 * selection, so we need at least OB_SELECT state.
 */
  if(!ob_ready(ob, OB_SELECT, "ed_integ") || sub_bad(sub, "ed_integ"))
    return 1;
  if(ut < 0 || ut >= sub->ntime) {
    lprintf(stderr, "ed_integ: Integration index out of range.\n");
    return 1;
  };
/*
 * Check the IF index if relevant.
 */
  if(selif && (cif < 0 || cif >= ob->nif)) {
    lprintf(stderr, "ed_integ: IF index out of range.\n");
    return 1;
  };
/*
 * Check the range of index where relevant.
 * Baseline selection.
 */
  if(selbase && (index < 0 || index >= sub->nbase)) {
    lprintf(stderr, "ed_integ: Baseline index out of range.\n");
    return 1;
  };
/*
 * Station selection.
 */
  if(!selbase && selstat && (index < 0 || index >= sub->nstat)) {
    lprintf(stderr, "ed_integ: Station index out of range.\n");
    return 1;
  };
/*
 * Mark the per-baseline weight sums of affected IFs as stale.
 */
  flag_baseline_weights(ob, selif ? cif : -1);
/*
 * Get a pointer to the descriptor of the integration being edited.
 */
  integ = &sub->integ[ut];
/*
 * If the edit refers to the IF that is currently in memory, apply it
 * immediately to the visibilities in the Observation structure.
 */
  if(ob_ready(ob, OB_GETIF, NULL) && (!selif || cif == ob->stream.cif)) {
/*
 * Loop through the baselines, editing only those selected.
 */
    Baseline *bptr = sub->base;
    Visibility *vis = integ->vis;
    int base;
    for(base=0; base<sub->nbase; base++,bptr++,vis++) {
/*
 * Edit the visibility only if selected.
 */
      if((selbase && base==index) ||
	 (selstat && (bptr->tel_a==index || bptr->tel_b==index)) ||
	 (!selbase && !selstat)) {
/*
 * Apply the requested edit operation.
 */
	if(doflag)
	  vis->bad |= FLAG_BAD;
	else
	  vis->bad &= ~FLAG_BAD;
      };
    };
  };
/*
 * Is deferred editing required?
 */
  if(ob->obed) {
    Edint *ed;          /* Pointer to the new edit node */
/*
 * Allocate the new Edint node.
 */
    ed = new_Edint(ob);
/*
 * No memory? If so apply the existing edits and flush the edit
 * list so that the free-list can be re-used. Then try to allocate the
 * new node again.
 */
    if(ed==NULL && (ed_flush(ob) || (ed=new_Edint(ob))==NULL))
      return 1;
/*
 * Initialize the new edit node.
 */
    ed->next = NULL;
    ed->cif = cif;
    ed->index = index;
    ed->doflag = doflag;
    ed->selbase = selbase;
    ed->selstat = selstat;
    ed->selchan = selchan;
    ed->selif = selif;
/*
 * Start a new edit list?
 */
    if(integ->edlist==NULL)
      integ->edlist = ed;
/*
 * Append to an existing list.
 */
    else {
      Edint *tail = integ->edlist;
/*
 * Find the tail of the list.
 */
      while(tail->next)
	tail = tail->next;
/*
 * Append to the tail of the list.
 */
      tail->next = ed;
    };
  };
  return 0;
}

/*.......................................................................
 * Apply pending edits to the data in the uvdata and ifdata scratch files,
 * and then clear the list of pending edits.
 *
 * Input:
 *  ob    Observation *  The descriptor of the Observation.
 * Output:
 *  return        int    0 - OK.
 *                       1 - Error.
 */
int ed_flush(Observation *ob)
{
/*
 * Note that it is ok to call this function before making a selection
 * if there are no edits to be flushed.
 */
  if(!ob_ready(ob, OB_ALLOC, "ed_flush"))
    return 1;
/*
 * Are there any deferred edits to be applied?
 */
  if(ob->obed && ob->obed->nused>0) {
/*
 * Note that edits are specific to a given channel-range selection, so
 * we need OB_SELECT state if there are any edits to be flushed.
 */
    if(!ob_ready(ob, OB_SELECT, NULL)) {
      lprintf(stderr,
	      "ed_flush: Can't flush edits without selection context.\n");
      clr_Obedit(ob);
      return 1;
    };
/*
 * Inform user of reason for delay.
 */
    lprintf(stdout, "Applying %d buffered edits.\n", ob->obed->nused);
/*
 * Apply edits to the uvdata scratch file.
 */
    if(ed_uvdata(ob))
      return 1;
/*
 * Apply edits to the ifdata scratch file.
 */
    if(ed_ifdata(ob))
      return 1;
/*
 * Clear the list of edits and re-create the free list.
 */
    if(clr_Obedit(ob))
      return 1;
  };
/*
 * Edits completed OK.
 */
  return 0;
}

/*.......................................................................
 * Apply deferred edits to the uvdata scratch file.
 * Private function of ed_flush().
 *
 * Input:
 *  ob   Observation *  The descriptor of the Observation being
 *                      edited.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
static int ed_uvdata(Observation *ob)
{
/*
 * Are there any edits to be applied?
 */
  if(ob->obed && ob->obed->nused>0) {
    Intrec *rec;         /* File integration record pointer */
    int ut;              /* The index of the integration being edited */
/*
 * Apply edits to each integration in turn.
 */
    rec = ob->rec;
    for(ut=0; ut<ob->nrec; ut++,rec++) {
      Integration *integ = rec->integ;
/*
 * Is there a list of edits to be applied to this integration?
 */
      if(integ->edlist) {
	int cif;        /* The index of the IF being processed */
/*
 * Edit one IF at a time.
 */
	for(cif=0; cif<ob->nif; cif++) {
/*
 * Find the range and number of edits.
 */
	  Edrange *er = ed_range(ob, integ, cif);
/*
 * Any edits to be applied?
 */
	  if(er->nedit>0) {
/*
 * Set the read/write range of the required portion of the uvdata
 * scratch file for the current IF.
 */
	    if(dp_irange(ob->dp, cif, cif) ||
	       dp_srange(ob->dp, 0, ob->npol-1) ||
	       dp_crange(ob->dp, er->cmin, er->cmax) ||
	       dp_brange(ob->dp, er->bmin, er->bmax))
	      return 1;
/*
 * Read the required portion of the IF from the uvdata scratch file.
 */
	    if(dp_read(ob->dp, ut))
	      return 1;
/*
 * Apply edits to the buffered visibilities.
 */
	    if(dp_edit(ob, integ, cif))
	      return 1;
/*
 * Write the edited data back to the uvdata scratch file.
 */
	    if(dp_write(ob->dp, ut))
	      return 1;
	  };
	};
      };
    };
/*
 * Make sure that the file is up to date.
 */
    if(dp_flush(ob->dp))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Apply deferred edits to the ifdata scratch file.
 * Private function of ed_flush().
 *
 * Input:
 *  ob   Observation *  The descriptor of the Observation being
 *                      edited.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
static int ed_ifdata(Observation *ob)
{
/*
 * Are there edits to be applied and an ifdata scratch file to apply them
 * to?
 */
  if(ob->obed && ob->obed->nused>0 && ob->ip) {
    int cif;        /* The index of the IF being processed */
/*
 * Edit each sampled IF, one IF at a time.
 */
    for(cif=0; (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
      Intrec *rec;         /* File integration record pointer */
      int ut;              /* The index of the integration being edited */
/*
 * Apply edits to each integration in turn.
 */
      rec = ob->rec;
      for(ut=0; ut<ob->nrec; ut++,rec++) {
	Integration *integ = rec->integ;
/*
 * Is there a list of edits to be applied to this integration?
 */
	if(integ->edlist) {
/*
 * Find the range and number of edits.
 */
	  Edrange *er = ed_range(ob, integ, cif);
/*
 * Any edits to be applied?
 */
	  if(er->nedit>0) {
/*
 * Set the read/write range of the required portion of the ifdata
 * scratch file for the current IF.
 */
	    if(ip_range(ob->ip, cif, er->bmin, er->bmax))
	      return 1;
/*
 * Read the required portion of the IF from the ifdata scratch file.
 */
	    if(ip_read(ob->ip, ut))
	      return 1;
/*
 * Apply edits to the buffered visibilities.
 */
	    if(ip_edit(ob, integ, cif))
	      return 1;
/*
 * Write the edited data back to the ifdata scratch file.
 */
	    if(ip_write(ob->ip, ut))
	      return 1;
	  };
	};
      };
    };
/*
 * Make sure that the file is up to date.
 */
    if(ip_flush(ob->ip))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Determine the range within the given IF to which edits are to be
 * applied.
 *
 * Input:
 *  ob    Observation *  The descriptor of the observation.
 *  integ Integration *  The descriptor of the integration being edited.
 *  cif           int    The IF being edited.
 * Output:
 *  return    Edrange *  Pointer to static internal structure used to record
 *                       the range and number of edits to be applied.
 */
static Edrange *ed_range(Observation *ob, Integration *integ, int cif)
{
  static Edrange er;  /* The container of the info to be returned */
  Chlist *cl;         /* IF spectral-line channel list */
  Subarray *sub;      /* The descriptor of the sub-array of the integration */
  Edint *ed;     /* Pointer into integ->edlist */ 
  int ca,cb;     /* Indexes of first and last required spectral-line channels.*/
  int ba,bb;     /* Indexes of first and last required baselines. */
  int first;     /* True for first iteration of the range finding loop */
/*
 * No edits counted yet.
 */
  er.nedit = 0;
/*
 * Get the descriptor of the sub-array from which the integration
 * is taken, and the list of selected spectral-line channels.
 */
  sub = integ->sub;
  cl = ob->ifs[cif].cl;
/*
 * Loop through the list of edits and determine the overall range of
 * spectral-line channels and baseline indexes.
 */
  first = 1;
  for(ed=integ->edlist; ed; ed=ed->next) {
/*
 * Can we honor any IF and channel selections?
 */
    if((ed->selif && ed->cif != cif) || (ed->selchan && cl==NULL))
      continue;
/*
 * Count the number of applicable edits.
 */
    er.nedit++;
/*
 * Find the range of baselines in the current edit.
 *
 * Single baseline selection?
 */
    if(ed->selbase) {
      ba = bb = ed->index;
/*
 * Telescope selection?
 */
    } else if(ed->selstat) {
/*
 * Find the index of the first baseline that involves the telescope.
 */
      for(ba=0; ba<sub->nbase; ba++) {
	if(sub->base[ba].tel_a == ed->index || sub->base[ba].tel_b == ed->index)
	  break;
      };
/*
 * Find the index of the last baseline that involves the telescope.
 */ 
      for(bb=sub->nbase-1; bb>ba; bb--) {
	if(sub->base[bb].tel_a == ed->index || sub->base[bb].tel_b == ed->index)
	  break;
      };
/*
 * Valid basline range?
 */
      if(bb >= sub->nbase) {
	lprintf(stderr, "ed_range: Failed to locate baseline.\n");
	er.nedit--;
	continue;
      };
    } else {
      ba = 0;
      bb = sub->nbase-1;
    };
/*
 * Find the range of spectral-line channels in the current edit.
 */
    if(ed->selchan) {      /* Source channels of current stream */
      ca = cl->ca;
      cb = cl->cb;
    } else {               /* All spectral-line channels */
      ca = 0;
      cb = ob->nchan - 1;
    };
/*
 * Expand the recorded minimum and maximum channel and baseline index
 * ranges to include the latest edit.
 */
    if(first) {     /* Initialize ranges on first iteration */
      first = 0;
      er.cmin = ca;
      er.cmax = cb;
      er.bmin = ba;
      er.bmax = bb;
    } else {        /* Expand ranges where necessary */
      if(ca < er.cmin)
	er.cmin = ca;
      if(cb > er.cmax)
	er.cmax = cb;
      if(ba < er.bmin)
	er.bmin = ba;
      if(bb > er.bmax)
	er.bmax = bb;
    };
  };
  return &er;
}

/*.......................................................................
 * Apply edits to the IF in the uvdata scratch file I/O buffer.
 *
 * Input:
 *  ob    Observation *  The descriptor of the observation.
 *  integ Integration *  The index of the integration being edited.
 *  cif           int    The index of the IF being edited.
 * Output:
 *  return        int    0 - OK.
 *                       1 - Error.
 */
static int dp_edit(Observation *ob, Integration *integ, int cif)
{
  Chlist *cl;       /* The list of selected channel ranges */
  Subarray *sub;    /* The sub-array to which the integration belongs */
  Edint *ed;        /* Pointer into edit operation list */
  Dif *ifp;         /* Pointer to the tree containing IF cif */
/*
 * Get the pointer to the branch of the buffer tree that holds the IF
 * to be edited.
 */
  ifp = &ob->dp->ifs[cif];
/*
 * Get the descriptor of the sub-array to which the integration belongs.
 */
  sub = integ->sub;
/*
 * Get the list of selected channel ranges.
 */
  cl = ob->ifs[cif].cl;
/*
 * Apply each edit in the order in the edit list.
 */
  for(ed=integ->edlist; ed; ed = ed->next) {
    int ca,cb;      /* Range of spectral-line channels */
    int fc;         /* Index of spectral-line channel being edited */
    int base;       /* Index of the baseline being edited */
    int nrange;     /* The number of channel ranges */
    int ir;         /* The index of the channel range being edited */
    Baseline *bptr; /* Pointer to baseline descriptor for index 'base' */
/*
 * Can we honor any IF and channel selections?
 */
    if((ed->selif && ed->cif != cif) || (ed->selchan && cl==NULL))
      continue;
/*
 * Finding baselines to edit is less efficient than finding spectral-line
 * channels so loop over baselines in the outer loop.
 */
    bptr = sub->base;
    for(base=0; base<sub->nbase; base++,bptr++) {
/*
 * Don't edit this baseline?
 */
      if((ed->selbase && base!=ed->index) ||
	 (ed->selstat && bptr->tel_a != ed->index && bptr->tel_b != ed->index))
	continue;
/*
 * Apply edits to each spectral-line channel in turn.
 */
      nrange = ed->selchan ? cl->nrange : 1;
      for(ir=0; ir<nrange; ir++) {
/*
 * Get the limits of the next channel range to be edited.
 */
	if(ed->selchan) {
	  ca = cl->range[ir].ca;
	  cb = cl->range[ir].cb;
	} else {
	  ca = 0;
	  cb = ob->nchan - 1;
	};
/*
 * Edit the latest channel range.
 */
	for(fc=ca; fc<=cb; fc++) {
	  Cvis *pol = ifp->chan[fc].base[base].pol;
	  int p;
	  for(p=0; p<ob->npol; p++,pol++)
	    pol->wt = ed->doflag ? -fabs(pol->wt) : fabs(pol->wt);
	};
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Apply edits to the IF in the ifdata scratch file I/O buffer.
 *
 * Input:
 *  ob     Observation *  The descriptor of the observation.
 *  integ  Integration *  The descriptor of the integration being edited.
 *  cif            int    The index of the IF being edited.
 * Output:
 *  return         int    0 - OK.
 *                        1 - Error.
 */
static int ip_edit(Observation *ob, Integration *integ, int cif)
{
  Subarray *sub; /* The sub-array to which the descriptor belongs */
  Edint *ed;     /* Pointer into edit operation list */
/*
 * Nothing to edit.
 */
  if(!ob->ifs[cif].cl)
    return 0;
/*
 * Get the sub-array to which the integratin belongs.
 */
  sub = integ->sub;
/*
 * Apply each edit in the order in the edit list.
 */
  for(ed=integ->edlist; ed; ed = ed->next) {
    int base;       /* Index of the baseline being edited */
    Baseline *bptr; /* Pointer to baseline descriptor for index 'base' */
/*
 * Can we honor any IF selections?
 */
    if(ed->selif && ed->cif != cif)
      continue;
/*
 * Locate the baselines to edit.
 */
    bptr = sub->base;
    for(base=0; base<sub->nbase; base++,bptr++) {
      Dvis *dvis;
/*
 * Don't edit this baseline?
 */
      if((ed->selbase && base!=ed->index) ||
	 (ed->selstat && bptr->tel_a != ed->index && bptr->tel_b != ed->index))
	continue;
/*
 * Locate the corresponding visibility in the ifdata buffer.
 */
      dvis = &ob->ip->dvis[base];
/*
 * Apply the requested edit operation.
 */
      dvis->wt = ed->doflag ? -fabs(dvis->wt) : fabs(dvis->wt);
    };
  };
  return 0;
}

/*.......................................................................
 * Apply defered edits to visibilities that have just been read from
 * the IF scratch file. This is meant to be called only from iniIF().
 *
 * Input:
 *  ob    Observation *  The descriptor of the observation containing
 *                       the un-edited visibilities.
 *  cif           int    The IF of the un-edited visibilities.
 * Output:
 *  return        int    0 - OK.
 *                       1 - Error.
 */
int app_Obedit(Observation *ob, int cif)
{
/*
 * Make sure that the edits have not already been applied.
 */
  if(!ob || ob->state != OB_RAWIF) {
    lprintf(stderr, "app_Obedit: Inappropriate call.\n");
    return 1;
  };
/*
 * Check the given IF index.
 */
  if(cif < 0 || cif >= ob->nif) {
    lprintf(stderr, "app_Obedit: IF index out of range.\n");
    return 1;
  };
/*
 * Is the given IF un-sampled by any of the selected channels?
 */
  if(!ob->ifs[cif].cl)
    return 0;
/*
 * Are there any edits to be applied?
 */
  if(ob->obed && ob->obed->nused>0) {
    Intrec *rec;       /* Pointer into ob->rec[] */
/*
 * Each integration has its own list of edits.
 */
    for(rec=ob->rec; rec<ob->rec + ob->nrec; rec++) {
      Integration *integ = rec->integ;
      Subarray *sub = integ->sub;
      Edint *ed;
      for(ed=integ->edlist; ed; ed=ed->next) {
/*
 * Can we honor any IF selections? Note that channel selections are
 * irrelevant at this point because if the current IF isn't sampled
 * by the current channel selection, we won't have got to this point
 * because of the above if(!ob->ifs[cif].cl) statement.
 */
	if(!ed->selif || ed->cif == cif) {
/*
 * Loop through the baselines of the current integration, editing only
 * those selected.  
 */
	  Baseline *bptr = sub->base;
	  Visibility *vis = integ->vis;
	  int base;
	  for(base=0; base<sub->nbase; base++,bptr++,vis++) {
/*
 * Edit the visibility only if selected.
 */
	    if((ed->selbase && base==ed->index) ||
	       (ed->selstat &&
		(bptr->tel_a==ed->index || bptr->tel_b==ed->index)) ||
	       (!ed->selbase && !ed->selstat)) {
/*
 * Apply the requested edit operation.
 */
	      if(ed->doflag)
		vis->bad |= FLAG_BAD;
	      else
		vis->bad &= ~FLAG_BAD;
	    };
	  };
	};
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Edit a given set of baselines over a given range of times.
 *
 * Input:
 *  ob    Observation *  The observation to be edited.
 *  doflag        int    If true flag the specified visibilities.
 *                       If false unflag them.
 *  spec         char *  A baseline specification string, usable by
 *                       read_Basespec().
 *  doall         int    If true, edit all channels and IFs. Otherwise
 *                       edit just the currently selected channels in the
 *                       currently selected IFs.
 *  mjd1, mjd2 double    The time limits of the editing, expressed in
 *                       UTC as Modified Julian Dates. If mjd1 is <=0.0,
 *                       the start time of the observation will be
 *                       substituted. If mjd2 is <=0.0, the end time
 *                       of the observation will be substituted.
 * Output:
 *  return        int    0 - OK.
 *                       1 - Error.
 */
int edit_baselines(Observation *ob, int doflag, char *spec, int doall,
		   double mjd1, double mjd2)
{
  double uta,utb;  /* The start and end times of the range to be edited, */
                   /* expressed as seconds into the year of observation. */
  int sa,sb;       /* The range of subarrays to edit */
  int irec;        /* A integration record index */
  Basespec *bs;    /* The parsed baseline specification */
/*
 * Check the arguments.
 */
  if(!ob || !spec) {
    lprintf(stderr, "edit_baselines: NULL argument(s).\n");
    return 1;
  };
/*
 * Parse the telescope specification.
 */
  bs = read_Basespec(ob, spec, NULL, 0);
  if(!bs)
    return 1;
/*
 * Inform the user about what is about to happen.
 */
  lprintf(stdout, "%s ", doflag ? "Flagging":"Unflagging");
  switch(bs->nfix) {
  case 0:
    lprintf(stdout, "all baselines");
    break;
  case 1:
    lprintf(stdout, "baselines of subarray %d", bs->isub+1);
    break;
  case 2:
    lprintf(stdout, "baselines of antenna %d:%s",
	    bs->isub+1, ob->sub[bs->isub].tel[bs->ta].name);
    break;
  case 3:
    lprintf(stdout, "baseline %d:%s-%s", bs->isub+1,
	    ob->sub[bs->isub].tel[bs->ta].name,
	    ob->sub[bs->isub].tel[bs->tb].name);
    break;
  default:
    lprintf(stderr, "\nedit_baselines: read_Basespec() error.\n");
    return 1;
  };
  if(doall)
    lprintf(stdout, " in all channels.\n");
  else
    lprintf(stdout, " in the currently selected channels.\n");
/*
 * Convert the MJD UTC's into the number of seconds into the year of
 * observation, and substitute the apropriate end time for times that
 * haven't been specified.
 */
  if(mjd1==0.0)
    uta = ob->rec[0].integ->ut;
  else
    uta = (mjd1 - ob->date.utc_ref) * daysec;
  if(mjd2==0.0)
    utb = ob->rec[ob->nrec-1].integ->ut;
  else
    utb = (mjd2 - ob->date.utc_ref) * daysec;
/*
 * Get the range of subarrays to be edited.
 */
  if(bs->nfix < 1) {
    sa = 0;
    sb = ob->nsub;
  } else {
    sa = sb = bs->isub;
  };
/*
 * Go through the data in order of integration record number.
 */
  for(irec=0; irec<ob->nrec; irec++) {
    Integration *integ = ob->rec[irec].integ;
/*
 * Get the sub-array index of the integration.
 */
    int isub = integ->sub - ob->sub;
/*
 * Is the integration inside the specified time range and sub-array
 * range?
 */
    if(integ->ut >= uta && integ->ut <= utb && isub >= sa && isub <= sb) {
      int base;   /* A baseline index */
/*
 * Get the subarray object.
 */
      Subarray *sub = integ->sub;
/*
 * Get the integration index within the sub-array.
 */
      int ut = integ - sub->integ;
/*
 * Handle the different types of baseline specification.
 */
      switch(bs->nfix) {
      case 0:
      case 1:              /* Edit all baselines */
	for(base=0; base<sub->nbase; base++) {
	  if(ed_integ(ob, integ->sub, ut, -1, doflag, 1, 0, !doall, 0, base))
	    return 1;
	};
	break;
      case 2:              /* Edit all baselines of a given antenna */
	if(ed_integ(ob, integ->sub, ut, -1, doflag, 0, 1, !doall, 0, bs->ta))
	  return 1;
	break;
      case 3:              /* Edit a specific single baseline */
	if(ed_integ(ob, integ->sub, ut, -1, doflag, 1, 0, !doall, 0,
		    loc_base(integ->sub, bs->ta, bs->tb)))
	  return 1;
	break;
      default:
	lprintf(stderr,
		"edit_baselines: Error in nfix returned by read_Baseline()\n");
	return 1;
      };
    };
  };
  return 0;
}

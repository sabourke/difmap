#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "obs.h"
#include "telspec.h"
#include "logio.h"
#include "baselist.h"

/*
 * Set the number of elements that are allocated when Basesel free-list
 * needs to be expanded.
 */

enum {ALLOC_STEP=20};

typedef struct Baseselmem Baseselmem;
struct Baseselmem {
  Basesel bsel[ALLOC_STEP];
  Baseselmem *next;
};

static struct {
  Basesel *free;
  int nused;
  Baseselmem *head;
} bsel_memory = {0,0,0};

static Baseselmem *new_Baseselmem(void);
static Basesel *new_Basesel(void);

/*.......................................................................
 * Expand the free list of Basespec structures held in bsel_memory.
 *
 * Output:
 *   return  Baseselmem *  The newly allocated block, or NULL on error.
 */
static Baseselmem *new_Baseselmem(void)
{
  Baseselmem *newmem;   /* The new memory block */
  int i;
/*
 * Attempt to allocate the new block.
 */
  newmem = (Baseselmem *) malloc(sizeof(Baseselmem));
  if(newmem == NULL) {
    lprintf(stderr, "new_Baseselmem: Insufficient memory.\n");
    return NULL;
  };
/*
 * Link all the elements of the chunk together and prepend the result to
 * the Basepec free list.
 */
  for(i=0; i<ALLOC_STEP-1; i++)
    newmem->bsel[i].next = &newmem->bsel[i+1];
  newmem->bsel[i].next = bsel_memory.free;
  bsel_memory.free = newmem->bsel;
/*
 * Prepend the new block to the Baseselmem list in 'bsel_memory'.
 */
  newmem->next = bsel_memory.head;
  bsel_memory.head = newmem;
/*
 * Return the new block.
 */
  return newmem;
}

/*.......................................................................
 * Get a Basesel node from the free list in bsel_memory.
 *
 * Output:
 *  return   Basesel *  The new node, or NULL on error.
 */
static Basesel *new_Basesel(void)
{
  Basesel *bsel;  /* The new node */
  if(!bsel_memory.free && new_Baseselmem()==NULL)
    return NULL;
  bsel = bsel_memory.free;
  bsel_memory.free = bsel->next;
  bsel_memory.nused++;
  bsel->next = NULL;
  return bsel;
}

/*.......................................................................
 * Return a baseline selection node to the free list in bsel_memory.
 *
 * Input:
 *  bgrp   Basegrp *   If the node is currently part of a list, send the
 *                     pointer to the list here to have it removed before
 *                     deletion. Otherwise send NULL.
 *  bsel  Basesel *   The node to be returned.
 * Output:
 *  return  Basesel *   Allways NULL.
 */
Basesel *del_Basesel(Basegrp *bgrp, Basesel *bsel)
{
  if(bsel) {
/*
 * Remove the node from a list?
 */
    if(bgrp) {
      Basesel *prev=NULL;
      Basesel *node;
      for(node=bgrp->bsel; node && node != bsel; prev=node,node=node->next)
	;
      if(node) {
	if(prev==NULL)
	  bgrp->bsel = node->next;
	else
	  prev->next = node->next;
	if(bgrp->tail == node)
	  bgrp->tail = prev;
	bgrp->nnode--;
      };
    };
/*
 * Prepend the node to the free list.
 */
    bsel->next = bsel_memory.free;
    bsel_memory.free = bsel;
    bsel_memory.nused--;
  };
  return NULL;
}

/*.......................................................................
 * Create an empty baseline group. This can then be used to record
 * which baselines of an observation to include and which to exclude
 * from an operation.
 *
 * Output:
 *  return    Basegrp *  The empty baseline group container.
 */
Basegrp *new_Basegrp(void)
{
  Basegrp *bgrp;   /* The return container */
/*
 * Allocate the container.
 */
  bgrp = (Basegrp *) malloc(sizeof(Basegrp));
  if(bgrp==NULL) {
    lprintf(stderr, "new_Basegrp: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize the container, at least to the point at which it can safely
 * be presented to del_Basegrp(bgrp).
 */
  bgrp->bsel = bgrp->tail = NULL;
  bgrp->nnode = 0;
  bgrp->next = NULL;
/*
 * Return the container for use.
 */
  return bgrp;
}

/*.......................................................................
 * Delete a baseline group container previously allocated with
 * new_Basegrp().
 *
 * Input:
 *  bgrp   Basegrp *  The container to be deleted.
 *  bgl   Bgrplist *  If bgrp is a member of a list of baseline groups
 *                    send the list here so that bgrp can be removed before
 *                    being deleted. Otherwise send NULL.
 * Output:
 *  return Basegrp *  The deleted container, ie. NULL.
 */
Basegrp *del_Basegrp(Basegrp *bgrp, Bgrplist *bgl)
{
  if(bgrp) {
/*
 * Extract the group from a list of groups?
 */
    if(bgl) {
      Basegrp *prev=NULL;
      Basegrp *node;
      for(node=bgl->bgrp; node && node != bgrp; prev=node,node=node->next)
	;
      if(node) {
	if(prev==NULL)
	  bgl->bgrp = node->next;
	else
	  prev->next = node->next;
	if(bgl->tail == node)
	  bgl->tail = prev;
	bgl->nsel--;
      };
    };
    clr_Basegrp(bgrp);
    free(bgrp);
  };
  return NULL;
}

/*.......................................................................
 * Clear the contents of a baseline group.
 *
 * Input:
 *  bgrp    Basegrp *  The group to be cleared.
 * Output:
 *  return  Basegrp *  The cleared group, or NULL on error.
 */
Basegrp *clr_Basegrp(Basegrp *bgrp)
{
  Basesel *node;  /* The node being deleted */
  Basesel *next;  /* The next node to be deleted */
/*
 * Delete all nodes in the list.
 */
  for(node=bgrp->bsel; node; node = next) {
    next = node->next;
    node = del_Basesel(NULL, node);
  };
  bgrp->bsel = NULL;
  bgrp->nnode = 0;
  bgrp->tail = NULL;
  return bgrp;
}

/*.......................................................................
 * Add a new entry to a baseline group.
 *
 * Input:
 *  ob  Observation *  The observation used to create bgrp.
 *  bgrp    Basegrp *  The baseline group container.
 *  bs     Basespec *  The baseline specification to be appended.
 *  include     int    If true, baselines refered to by *bs are to be
 *                     included, otherwise they are baselines to be
 *                     excluded.
 * Output:
 *  return  Basesel *  The new baseline selection node in the group, or NULL
 *                     on error.
 */
Basesel *add_Basesel(Observation *ob, Basegrp *bgrp, Basespec *bs, int include)
{
  Basesel *bsel;  /* The baseline selection node to be returned */
/*
 * Initialize the baseline specification and check that it is valid.
 */
  if(next_base(ob, FIND_FIRST, 1, bs->nfix, 1, 0, 1, bs))
    return NULL;
/*
 * Get a new baseline selection structure from the free-list in bsel_memory.
 */
  if(!(bsel = new_Basesel()))
    return NULL;
/*
 * Append the new entry at the end of the baseline group.
 */
  bsel->include = include;
  bsel->bs = *bs;
  bsel->next = NULL;
  if(bgrp->bsel==NULL) {
    bgrp->bsel = bgrp->tail = bsel;
  } else {
    bgrp->tail->next = bsel;
    bgrp->tail = bsel;
  };
  bgrp->nnode++;
/*
 * Return the new node.
 */
  return bsel;
}

/*.......................................................................
 * Create a new list of baseline groups.
 *
 * Output:
 *  return  Bgrplist *  The new container, or NULL on error.
 */
Bgrplist *new_Bgrplist(void)
{
  Bgrplist *bgl;
/*
 * Allocate the new container.
 */
  bgl = (Bgrplist *) malloc(sizeof(Bgrplist));
  if(!bgl) {
    lprintf(stderr, "new_Bgrplist: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize the list.
 */
  bgl->bgrp = bgl->tail = NULL;
  bgl->nsel = 0;
  return bgl;
}

/*.......................................................................
 * Delete a list of baseline groups.
 *
 * Input:
 *  bgl    Bgrplist *   The list to be deleted.
 * Output:
 *  return Bgrplist *   Allways NULL.
 */
Bgrplist *del_Bgrplist(Bgrplist *bgl)
{
  if(bgl) {
    Basegrp *node;  /* The node being deleted */
    Basegrp *next;  /* The next node to be deleted */
/*
 * Delete all nodes in the list.
 */
    for(node=bgl->bgrp; node; node = next) {
      next = node->next;
      node = del_Basegrp(node, NULL);
    };
    free(bgl);
  };
  return NULL;
}

/*.......................................................................
 * Append a baseline group to a list of baseline groups.
 * The new list can either be presented directly, or as a baseline
 * group string of the form expected by read_Basegrp().
 *
 * Input:
 *  ob Observation *   The observation associated with all arguments.
 *  bgl   Bgrplist *   The list to append to.
 *  bgrp   Basegrp *   The baseline group to append, or NULL if
 *                     bgrp_str is to be used. This becomes the property
 *                     of the bgl list and will be deleted when the list
 *                     is deleted. 
 *  bgrp_str  char *   If bgrp==NULL then bgrp_str must not be NULL.
 *                     It must be of the form expected by read_Basegrp(),
 *                     and will be used to initialize a new Basegrp node.
 * Output:
 *  return Basegrp *   The baseline group node added, or NULL on error.
 *                     On error 'bgrp' is deleted.
 */
Basegrp *add_Basegrp(Observation *ob, Bgrplist *bgl, Basegrp *bgrp,
		     char *bgrp_str)
{
/*
 * Nothing to add?
 */
  if((bgrp && bgrp->nnode<1) || (!bgrp && !bgrp_str)) {
    lprintf(stderr, "add_Basegrp: Empty baseline group.\n");
    return del_Basegrp(bgrp, NULL);
  };
/*
 * Allocate a baseline group.
 */
  if(!bgrp) {
    if(!(bgrp=new_Basegrp()))
      return NULL;
/*
 * Read a baseline group string.
 */
    if(read_Basegrp(ob, bgrp, bgrp_str, NULL))
      return del_Basegrp(bgrp, NULL);
  };
/*
 * Append the new entry at the end of the baseline group list.
 */
  if(bgl->bgrp==NULL) {
    bgl->bgrp = bgl->tail = bgrp;
  } else {
    bgl->tail->next = bgrp;
    bgl->tail = bgrp;
  };
  bgl->nsel++;
/*
 * Return the new node.
 */
  return bgrp;
}

/*.......................................................................
 * Parse a string of inclusive and exclusive baseline specifications
 * and append the results to a given baseline group.
 *
 * Input:
 *  ob   Observation *   The descriptor of the observation to which the
 *                       baseline selection string refers.
 *  bgrp     Basegrp *   The baseline group container to add to.
 *  string      char *   The baseline group string to parse. This
 *                       consists of a string of baseline
 *                       specifications separated by + or !.
 *                       Specifications preceded by + are inclusive,
 *                       while those preceded by !, are exclusive.
 *                       The first specification is implicitly
 *                       inclusive if not otherwise specified. If the
 *                       first specification is preceded by a !, all
 *                       baselines are included before interpretting
 *                       the first specification.
 *  endp        char **  If endp!=NULL, then a pointer to the next
 *                       unprocessed character in string[] will be
 *                       assigned to *endp. If endp==NULL an error
 *                       message will be issued if there are any input
 *                       characters remaining after the last valid
 *                       specification.
 * Output:
 *  return       int     0 - OK.
 *                       1 - Error.
 */
int read_Basegrp(Observation *ob, Basegrp *bgrp, char *string, char **endp)
{
  char *sptr = string;  /* Pointer to the next unprocessed character */
  int include;          /* True if the next selection is inclusive */
  int first=1;          /* True for the first specification */
/*
 * Check arguments.
 */
  if(!bgrp || !string) {
    lprintf(stderr, "read_Basegrp: NULL baseline group.\n");
    return 1;
  };
/*
 * Initialize the returned end pointer.
 */
  if(endp)
    *endp = string;
/*
 * Read one baseline specification from the string at a time.
 */
  while(1) {
    Basespec *bs;  /* A baseline specification descriptor */
/*
 * Skip white space to the next baseline specification separator.
 */
    while(*sptr && isspace((int)*sptr))
      sptr++;
    switch(*sptr) {
    case '+':
      include = 1;    /* The next specification is inclusive */
      sptr++;
      break;
    case '!':
      include = 0;    /* The next specification is exclusive */
      sptr++;
      break;
    case '\0':        /* Empty string, or end of input */
      if(first) {
	include = 1;
      } else {
	if(endp) *endp = sptr;
	return 0;
      };
      break;
    default:
      if(first) {
	include = 1;
      } else {
	if(!endp) {
	  lprintf(stderr, "read_Basegrp: Unexpected separator: %s.\n", sptr);
	  return 1;
	} else {
	  *endp = sptr;
	  return 0;
	};
      };
      break;
    };
/*
 * If the first specification is exclusive, start by including all
 * baselines.
 */
    if(first && !include) {
      bs = read_Basespec(ob, "", NULL, 0);
      if(!bs || !add_Basesel(ob, bgrp, bs, 1))
	return 1;
    };
/*
 * Read the next specification from the input string.
 */
    bs = read_Basespec(ob, sptr, &sptr, 0);
/*
 * Append the new specification to the current group.
 */
    if(!bs || !add_Basesel(ob, bgrp, bs, include))
      return 1;
    first = 0;
  };
  return 0;
}

/*.......................................................................
 * Write a baseline group string in the form read by read_Basegrp().
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 *  bgrp    Basegrp *  The baseline group.
 *  n           int    The max number of characters (including '\0') to
 *                     place in s[].
 * Input/Output:
 *  s          char *  An array of at least 'n' characters for the return
 *                     string.
 * Output:
 *  return      int    The number of characters used from s[], or -1
 *                     on error, or -2 on truncation.
 */
int write_Basegrp(Observation *ob, Basegrp *bgrp, int n, char *s)
{
  Basesel *bsel; /* The node being displayed */
  int nused;     /* The number of characters used in s[] */
/*
 * Check arguments.
 */
  if(!bgrp) {
    lprintf(stderr, "write_Basegrp: NULL baseline group.\n");
    return 1;
  };
  if(!s) {
    lprintf(stderr, "write_Basegrp: NULL output string.\n");
    return 1;
  };
/*
 * Terminate 's' in case there are no specifications to be written.
 */
  if(n>0)
    s[0] = '\0';
/*
 * Write one specification at a time.
 */
  nused = 0;
  for(bsel=bgrp->bsel; bsel; bsel=bsel->next) {
    if(bsel!=bgrp->bsel) {
      if(n-nused < 4)
	return -2;
      strcat(s + nused, bsel->include ? " + " : " ! ");
      nused += 3;
    };
    {
      int nnew = write_Basespec(ob, &bsel->bs, 0, 0, n-nused, s + nused);
      if(nnew < 0)
	return nnew;
      nused += nnew;
    };
  };
  return nused;
}

/*.......................................................................
 * Delete a baseline list container and its contents.
 *
 * Input:
 *  blist   Baselist *  The baseline list container to be deleted.
 * Output:
 *  return  Baselist *  The deleted container, ie. NULL.
 */
Baselist *del_Baselist(Baselist *blist)
{
  if(blist) {
    if(blist->bsub) {
      if(blist->bsub->baselines)
	free(blist->bsub->baselines);
      free(blist->bsub);
    };
    free(blist);
  };
  return NULL;
}

/*.......................................................................
 * Create and fill a container of baseline index lists (one list per
 * sub-array), based on the contents of a baseline group.
 *
 * Input:
 *   ob  Observation * The observation used to construct 'bgrp'.
 *   bgrp    Basegrp * A baseline group.
 * Output:
 *   return Baselist * A container of baseline index lists for each
 *                     sub-array.
 */
Baselist *new_Baselist(Observation *ob, Basegrp *bgrp)
{
  Baselist *blist;   /* The pointer to the container to be returned */
  int isub;          /* The index of the sub-array being processed */
  int i;
/*
 * Sanity checks.
 */
  if(bgrp==NULL) {
    lprintf(stderr, "new_Baselist: Invalid baseline group.\n");
    return NULL;
  };
/*
 * Allocate the container.
 */
  blist = (Baselist *) malloc(sizeof(Baselist));
  if(!blist) {
    lprintf(stderr, "new_Baselist: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize the container at least to the point at which it can
 * be safely sent to del_Baselist().
 */
  blist->ntotal = 0;
  blist->nsub = ob->nsub;
  blist->bsub = NULL;
/*
 * Allocate the array of baseline lists.
 */
  blist->bsub = (Bsublist *) malloc(sizeof(Bsublist) * blist->nsub);
  if(!blist->bsub) {
    lprintf(stderr, "new_Baselist: Insufficient memory.\n");
    return del_Baselist(blist);
  };
/*
 * Initialize the elements of the sub-array baseline list containers.
 */
  for(i=0; i < blist->nsub; i++) {
    blist->bsub[i].nbase = 0;
    blist->bsub[i].baselines = NULL;
  };
/*
 * Count the number of baselines included in each sub-array.
 */
  for(isub=0; isub<ob->nsub; isub++) {
    Bsublist *bsub = blist->bsub + isub;
    bsub->nbase = size_Basegrp(ob, bgrp, isub);
    blist->ntotal += bsub->nbase;
  };
/*
 * Allocate all of the sub-array baseline arrays as a single array.
 */
  blist->bsub->baselines = (int *) malloc(sizeof(int) * blist->ntotal);
  if(!blist->bsub->baselines) {
    lprintf(stderr, "new_Baselist: Insufficient memory.\n");
    return del_Baselist(blist);
  };
/*
 * Distribute the baseline array between sub-arrays.
 */
  for(isub=0; isub < ob->nsub; isub++) {
    Bsublist *bsub = blist->bsub + isub;
    bsub[1].baselines = bsub->baselines + bsub->nbase;
  };
/*
 * Check for inclusion of each baseline, and record those that are included.
 */
  for(isub=0; isub<ob->nsub; isub++) {
    Bsublist *bsub = blist->bsub + isub;
    int nused = 0;
    if(bsub->nbase > 0) {
      int base;      /* The index of the baseline being checked */
      int nbase=ob->sub[isub].nbase;
      for(base=0; base<nbase; base++) {
	if(in_Basegrp(ob, isub, base, bgrp))
	  bsub->baselines[nused++] = base;
      };
    };
  };
  return blist;
}

/*.......................................................................
 * Report whether a given baseline is included in a given baseline group.
 *
 * Input:
 *  ob  Observation *   The descriptor of the observation with which the
 *  isub        int     The index of the sub-array of the baseline.
 *  base        int     The index of the baseline in sub-array 'isub'.
 *                      baseline selection list was created.
 *  bgrp    Basegrp *   The baseline group.
 * Output:
 *  return      int     0 - The baseline is not included.
 *                      1 - The baseline is included.
 */
int in_Basegrp(Observation *ob, int isub, int base, Basegrp *bgrp)
{
  Baseline *b;     /* Baseline descriptor of baseline 'base' */
  Basesel *bsel;   /* A baseline selection list node */
  int include = 0; /* True if the baseline is included */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "in_Basegrp") || !bgrp)
    return 0;
  if(isub < 0 || isub >= ob->nsub) {
    lprintf(stderr, "in_Basegrp: Sub-array index out of range.\n");
    return 0;
  };
  if(base < 0 || base >= ob->sub[isub].nbase) {
    lprintf(stderr, "in_Basegrp: Baseline index out of range.\n");
    return 0;
  };
/*
 * Get the descriptor of the baseline.
 */
  b = &ob->sub[isub].base[base];
/*
 * See if the baseline group results in inclusion of baseline 'base'.
 */
  for(bsel=bgrp->bsel; bsel; bsel = bsel->next) {
/*
 * See if the current baseline specification cites baseline 'base'.
 */
    Basespec *bs = &bsel->bs;
    int cited = 0;
    switch(bs->nfix) {
    case 0:
      cited = 1;
      break;
    case 1:
      cited = isub==bs->isub;
      break;
    case 2:
      cited = isub==bs->isub && (bs->ta==b->tel_a || bs->ta==b->tel_b);
      break;
    case 3:
    default:
      cited = isub==bs->isub &&
	((bs->ta==b->tel_a && bs->tb==b->tel_b) ||
	 (bs->ta==b->tel_b && bs->tb==b->tel_a));
      break;
    };
/*
 * If the baseline is cited in the baseline specification, see if this
 * results in inclusion or exclusion of the baseline.
 */
    if(cited)
      include = bsel->include;
  };
  return include;
}

/*.......................................................................
 * Count the total number of baselines in a baseline group.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation used to create
 *                     'bgrp'.
 *  bgrp    Basegrp *  The baseline group to process.
 *  isub        int    The sub-array to count baselines in, or -1 for
 *                     all baselines in the observation.
 * Output:
 *  return      int    The number of baselines included.
 */
int size_Basegrp(Observation *ob, Basegrp *bgrp, int isub)
{
  int ssub;    /* The index of the first sub-array to be searched */
  int esub;    /* The index of the last sub-array to be searched */
  int nused;   /* The total number of included baselines */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "size_Basegrp") || !bgrp)
    return 0;
/*
 * Determine the range of sub-arrays to be searched.
 */
  if(isub >= 0 && isub < ob->nsub)
    ssub = esub = isub;
  else if(isub == -1) {
    ssub = 0;
    esub = ob->nsub - 1;
  } else {
    lprintf(stderr, "size_Basegrp: Sub-array index out of range.\n");
    return 0;
  };
/*
 * Count baselines in the selected sub-arrays.
 */
  nused = 0;
  for(isub=ssub; isub<=esub; isub++) {
    int nbase = ob->sub[isub].nbase;
    int base;    /* The index of the baseline being checked */
    for(base=0; base<nbase; base++) {
      if(in_Basegrp(ob, isub, base, bgrp))
	nused++;
    };
  };
  return nused;
}

/*.......................................................................
 * Search a baseline group for the next included baseline in a given
 * ordinal direction.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation used to create
 *                     'bgrp'.
 *  bgrp    Basegrp *  The baseline group.
 *  forward     int    Search forward wrt *s_isub and *s_base.
 * Input/Output:
 *  s_isub      int *  On input, the sub-array index to start searching from.
 *                     On output, the located sub-array if found.
 *  s_base      int *  On input, the baseline index in sub-array *s_isub
 *                     to start searching from. 
 *                     On output, the located baseline if found.
 * Output:
 *  return      int    0 - Baseline not found (*s_isub and *s_base are
 *                         unchanged). 
 *                     1 - A baseline was located and recorded in
 *                         *s_isub and *s_base.
 */
int srch_Basegrp(Observation *ob, Basegrp *bgrp, int forward,
		 int *s_isub, int *s_base)
{
  int base;    /* The index of the baseline being checked */
  int isub;    /* The index of the sub-array being checked */
  int first=1; /* True on the first iteration of the search loop */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "srch_Basegrp"))
    return 0;
  if(bgrp==NULL) {
    lprintf(stderr, "srch_Basegrp: Invalid baseline group.\n");
    return 0;
  };
  if(!s_isub || !s_base) {
    lprintf(stderr, "srch_Basegrp: NULL return arguments.\n");
    return 0;
  };
/*
 * Get the input indexes.
 */
  isub = *s_isub;
  base = *s_base;
/*
 * Search forward?
 */
  if(forward) {
/*
 * Start from the start of the baseline range currently below it.
 */
    if(isub < 0 || base < 0)
      isub = base = 0;
    else
      base++;
/*
 * Search forwards through the sub-arrays.
 */
    for( ; isub<ob->nsub; isub++) {
      int nbase = ob->sub[isub].nbase;
      for(base = first ? base : 0; base < nbase; base++) {
	if(in_Basegrp(ob, isub, base, bgrp)) {
	  *s_isub = isub;
	  *s_base = base;
	  return 1;
	};
      };
      first = 0;
    };
/*
 * Search backwards?
 */
  } else {
/*
 * Start from the end of the baseline range if currently above it.
 */
    if(isub > ob->nsub) {
      isub = ob->nsub - 1;
      base = ob->sub[isub].nbase - 1;
    } else {
      base--;
    };
/*
 * Search backwards through the sub-arrays.
 */
    for( ; isub >= 0; isub--) {
      for(base = first ? base : ob->sub[isub].nbase-1; base >= 0; base--) {
	if(in_Basegrp(ob, isub, base, bgrp)) {
	  *s_isub = isub;
	  *s_base = base;
	  return 1;
	};
      };
      first = 0;
    };
  };
/*
 * No included baselines were found in the specified direction.
 */
  return 0;
}

#include <stdio.h>
#include <stdlib.h>

#include "obs.h"
#include "logio.h"

/*
 * When indexing the sub-array integrations a list of ob->nsub containers,
 * one per sub-array will be maintained in order of increasing record number
 * of the next unused integration in the associated sub-arrays. Define the
 * structure of one node of such a list.
 */
typedef struct Subnode {
  long irec;           /* Record number of next un-indexed integration */
  struct Subnode *next;/* Pointer to sub-array with next higher record number */
  Integration *integ;  /* Pointer to next un-indexed integration */
  int ntime;           /* Number of integrations still available */
} Subnode;

static Subnode *sub_insert(Subnode *head, Subnode *node);

/*.......................................................................
 * (Re-)allocate and initialize an ob->rec[] array. This array contains
 * a sorted list of integrations in all sub-array integration arrays.
 * This function leaves each ob->rec[].integ NULL. Once all sub-array
 * integrations have been filled, call ini_Intrec() to initialize this
 * array.
 *
 * Input:
 *  ob      Observation *  The descriptor of the observation.
 *  nrec            int    The required number of records.
 * Output:
 *  return       Intrec *  The revised array.
 */
Intrec *new_Intrec(Observation *ob, int nrec)
{
  Intrec *rec;   /* Pointer into ob->rec */
  int irec;      /* The index of the record being processed */
/*
 * Valid Observation descriptor?
 */
  if(ob==NULL) {
    lprintf(stderr, "new_Intrec: NULL Observation descriptor intercepted.\n");
    return NULL;
  };
/*
 * Allocate a new array?
 */
  if(ob->rec==NULL) {
    ob->rec = malloc(sizeof(Intrec) * nrec);
    if(ob->rec==NULL) {
      lprintf(stderr, "new_Intrec: Insufficient memory.\n");
      return NULL;
    };
  }
/*
 * Re-size an existing array?
 */
  else if(ob->nrec != nrec) {
    rec = realloc(ob->rec, sizeof(Intrec) * nrec);
    if(rec)
      ob->rec = rec;
    else if(nrec > ob->nrec) {
      lprintf(stderr, "new_Intrec: Insufficient memory.\n");
      return NULL;
    };
  };
/*
 * Record the new number of records.
 */
  ob->nrec = nrec;
/*
 * Initialize all elements.
 */
  rec = ob->rec;
  for(irec=0; irec<nrec; irec++,rec++)
    rec->integ = NULL;
  return ob->rec;
}

/*.......................................................................
 * Delete the Intrec array of an observation.
 *
 * Input:
 *  ob    Observation *  The descriptor of the containing observation.
 * Output:
 *  return     Intrec *  The deleted array - ie. Allways NULL.
 */
Intrec *del_Intrec(Observation *ob)
{
  if(ob && ob->rec) {
    free(ob->rec);
    ob->rec = NULL;
  };
  return NULL;
}

/*.......................................................................
 * Index the integrations of all sub-arrays in increasing time-order
 * in ob->rec[]. Before calling this function all sub-arrays integrations
 * must have been initialized, and the sum of ob->sub[*].ntime must
 * equal ob->nrec.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int ini_Intrec(Observation *ob)
{
  Intrec *rec;      /* The record being initialized */
  Subarray *sub;    /* Descriptor of a sub-array */
  Subnode *node;    /* Next sub-node element */
  Subnode *head;    /* Head of sub-node list */
  Subnode *submem;  /* Pointer to start of array of sub-nodes */
  int isub;         /* Index of a sub-array */
  int irec;         /* The index of the record being intialized */
  int nrec=0;       /* Sum of sub-array ntime's */
/*
 * Do we have any data to index?
 */
  if(!ob_ready(ob, OB_DATA, "ini_Intrec"))
    return 1;
/*
 * Revert the recorded observation state to an unindexed state.
 */
  ob->state = OB_DATA;
/*
 * Count the total number of integrations.
 */
  sub = ob->sub;
  for(isub=0; isub<ob->nsub; isub++,sub++)
    nrec += sub->ntime;
/*
 * This must equal the allocated size of ob->rec.
 */
  if(nrec != ob->nrec) {
    lprintf(stderr, "ini_Intrec: Inconsistent integration count.\n");
    return 1;
  };
/*
 * Allocate ob->nsub sorting nodes.
 */
  submem = malloc(sizeof(Subnode) * ob->nsub);
  if(submem==NULL) {
    lprintf(stderr, "ini_Intrec: Insufficient memory.\n");
    return 1;
  };
/*
 * Initialize the array nodes.
 */
  node = submem;
  sub = ob->sub;
  for(isub=0; isub<ob->nsub; isub++,node++,sub++) {
    node->integ = sub->integ;
    node->irec = node->integ->irec;
    node->ntime = sub->ntime;
    node->next = NULL;
  };
/*
 * Sort the array into a linked list in order of record number.
 * A simple inefficient insertion sort should be sufficient, unless it ever
 * becomes common to have large numbers of sub-arrays in an observation.
 */
  node = submem;
  head = NULL;
  for(isub=0; isub<ob->nsub; isub++,node++)
    head = sub_insert(head, node);
/*
 * Extract pointers to integrations from each sub-array in record order
 * and record them in ob->rec[].
 */
  rec = ob->rec;
  for(irec=0; irec<ob->nrec; irec++,rec++) {
/*
 * If there are no more integrations in the current sub-array node, remove
 * the node from the list.
 */
    if(head->ntime<=0)
      head = head->next;
/*
 * If the head of the list is no longer the Subnode with the earliest
 * record number, move it to its new position in the sorted list.
 */
    if(head->next && head->irec > head->next->irec)
      head = sub_insert(head, NULL);
/*
 * Check that the record number of the new integration agrees with
 * irec.
 */
    if(irec != head->irec) {
      lprintf(stderr, "ini_Intrec: Out of order record number encountered.\n");
      free(submem);
      return 1;
    };
/*
 * Extract the next integration from the head Subnode of the list.
 */
    rec->integ = head->integ++;
    if(--head->ntime > 0)
      head->irec = head->integ->irec;
  };
/*
 * Release the temporary Subnode array.
 */
  free(submem);
/*
 * Mark the data as successfully indexed.
 */
  ob->state = OB_INDEX;
  return 0;
}

/*.......................................................................
 * Insert a Subnode node in a Subnode list such that the list is kept
 * in order of increasing time. As a common case, if node==NULL the
 * head of the list will be unlinked and re-placed in its correct
 * position. Both head and node must not be NULL at the same time.
 *
 * Input:
 *  head   Subnode *  The current head of the sorted list.
 *                    If head==NULL make 'node' the head of the list.
 *  node   Subnode *  The node to insert. If node==NULL, unlink the
 *                    head of the list and use it as the node to be
 *                    inserted.
 * Output:
 *  return Subnode *  The new head of the list.
 */
static Subnode *sub_insert(Subnode *head, Subnode *node)
{
/*
 * Remove the head of the list?
 */
  if(node==NULL) {
    node = head;
    head = head->next;
  };
/*
 * New list?
 */
  if(head==NULL) {
    head = node;
    node->next = NULL;
  } else {
    Subnode *last = NULL;    /* Pointer to node preceding the slot for node */
    Subnode *next;           /* Pointer to node following the slot for node */
    long irec = node->irec;  /* The record index of the first integration */
                             /* in the new node. */
/*
 * Find the appropriate position in the linked list for the new node.
 */
    for(next=head; next && next->irec < irec; next = next->next)
      last = next;
/*
 * Insert at head of list?
 */
    if(last==NULL) {
      node->next = head;
      head = node;
/*
 * Insert inside list?
 */
    } else if(next) {
      node->next = next;
      last->next = node;
/*
 * Append to tail of list.
 */
    } else {
      last->next = node;
      node->next = NULL;
    };
  };
/*
 * Return the new head of the list.
 */
  return head;
}

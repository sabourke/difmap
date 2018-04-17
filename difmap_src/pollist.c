#include <stdio.h>
#include <stdlib.h>

#include "obs.h"
#include "pollist.h"
#include "logio.h"

/*
 * Set the number of elements that are allocated when the polnode free-list
 * needs to be expanded.
 */
enum {ALLOC_STEP=20};

/*
 * Declare the typs of a free-list block of polarization nodes.
 */
typedef struct Polnodemem Polnodemem;
struct Polnodemem {
  Polnode polnode[ALLOC_STEP];
  Polnodemem *next;
};

static struct {
  Polnode *free;
  int nused;
  Polnodemem *head;
} polnode_memory = {0,0,0};

static Polnodemem *new_Polnodemem(void);
static Polnode *new_Polnode(void);

/*.......................................................................
 * Expand the free list of Polnode structures held in polnode_memory.
 *
 * Output:
 *   return  Polnodemem *  The newly allocated block, or NULL on error.
 */
static Polnodemem *new_Polnodemem(void)
{
  Polnodemem *newmem;   /* The new memory block */
  int i;
/*
 * Attempt to allocate the new block.
 */
  newmem = (Polnodemem *) malloc(sizeof(Polnodemem));
  if(newmem == NULL) {
    lprintf(stderr, "new_Polnodemem: Insufficient memory.\n");
    return NULL;
  };
/*
 * Link all the elements of the chunk together and prepend the result to
 * the Polnode free list.
 */
  for(i=0; i<ALLOC_STEP-1; i++)
    newmem->polnode[i].next = &newmem->polnode[i+1];
  newmem->polnode[i].next = polnode_memory.free;
  polnode_memory.free = newmem->polnode;
/*
 * Prepend the new block to the Polnodemem list in 'polnode_memory'.
 */
  newmem->next = polnode_memory.head;
  polnode_memory.head = newmem;
/*
 * Return the new block.
 */
  return newmem;
}

/*.......................................................................
 * Get a Polnode node from the free list in polnode_memory.
 *
 * Output:
 *  return   Polnode *  The new node, or NULL on error.
 */
static Polnode *new_Polnode(void)
{
  Polnode *polnode;  /* The new node */
  if(!polnode_memory.free && new_Polnodemem()==NULL)
    return NULL;
  polnode = polnode_memory.free;
  polnode_memory.free = polnode->next;
  polnode_memory.nused++;
  polnode->next = NULL;
  return polnode;
}

/*.......................................................................
 * Return a polarization node to the free list in polnode_memory.
 *
 * Input:
 *  pl       Pollist * If the node is currently part of a list, send the
 *                     pointer to the list here to have it removed before
 *                     deletion. Otherwise send NULL.
 *  polnode  Polnode * The node to be returned.
 * Output:
 *  return   Polnode * Allways NULL.
 */
Polnode *del_Polnode(Pollist *pl, Polnode *polnode)
{
  if(polnode) {
/*
 * Remove the node from a list?
 */
    if(pl) {
      Polnode *prev=NULL;
      Polnode *node;
      for(node=pl->head; node && node != polnode; prev=node,node=node->next)
	;
      if(node) {
	if(prev==NULL)
	  pl->head = node->next;
	else
	  prev->next = node->next;
	if(pl->tail == node)
	  pl->tail = prev;
	pl->npol--;
      };
    };
/*
 * Prepend the node to the free list.
 */
    polnode->next = polnode_memory.free;
    polnode_memory.free = polnode;
    polnode_memory.nused--;
  };
  return NULL;
}

/*.......................................................................
 * Create an empty list of polarizations.
 *
 * Output:
 *  return    Pollist *  The empty polarization selection list
 *                       container.
 */
Pollist *new_Pollist(void)
{
  Pollist *pl;   /* The return container */
/*
 * Allocate the container.
 */
  pl = (Pollist *) malloc(sizeof(Pollist));
  if(!pl) {
    lprintf(stderr, "new_Pollist: Insufficient memory.\n");
    return NULL;
  };
/*
 * Initialize the container, at least to the point at which it can safely
 * be presented to del_Pollist(pl).
 */
  pl->head = pl->tail = NULL;
  pl->npol = 0;
/*
 * Return the container for use.
 */
  return pl;
}

/*.......................................................................
 * Delete a polarization selection list container previously allocated with
 * new_Pollist().
 *
 * Input:
 *  pl    Pollist *  The container to be deleted.
 * Output:
 *  return Pollist *  The deleted container, ie. NULL.
 */
Pollist *del_Pollist(Pollist *pl)
{
  if(pl) {
    clr_Pollist(pl);
    free(pl);
  };
  return NULL;
}

/*.......................................................................
 * Clear the contents of a polarization selection list.
 *
 * Input:
 *  pl     Pollist *  The list to be cleared.
 * Output:
 *  return Pollist *  The cleared list, or NULL on error.
 */
Pollist *clr_Pollist(Pollist *pl)
{
  Polnode *node;  /* The node being deleted */
  Polnode *next;  /* The next node to be deleted */
/*
 * Delete all nodes in the list.
 */
  for(node=pl->head; node; node = next) {
    next = node->next;
    node = del_Polnode(NULL, node);
  };
  pl->head = NULL;
  pl->tail = NULL;
  pl->npol = 0;
  return pl;
}

/*.......................................................................
 * Add a new node to a polarization list.
 *
 * Input:
 *  ob  Observation *  The observation in which the polarization must
 *                     be found, or NULL if irrelevant.
 *  pl      Pollist *  The polarization list container.
 *  pol      Stokes    The polarization to be appended to the list.
 * Output:
 *  return   Polnode *  The new polarization selection node in the list, or NULL
 *                     on error.
 */
Polnode *add_Polnode(Observation *ob, Pollist *pl, Stokes pol)
{
  Polnode *polnode;  /* The polarization node to be returned */
/*
 * Check that the polarization can be derived from the polarizations
 * recorded in the observation.
 */
  if(ob && get_Obpol(ob, pol, 1, NULL))
    return NULL;
/*
 * Get a new polarization node from the free-list in polnode_memory.
 */
  if(!(polnode = new_Polnode()))
    return NULL;
/*
 * Append the new entry at the end of the polarization selection list.
 */
  polnode->pol = pol;
  polnode->next = NULL;
  if(pl->head==NULL) {
    pl->head = pl->tail = polnode;
  } else {
    pl->tail->next = polnode;
    pl->tail = polnode;
  };
  pl->npol++;
/*
 * Return the new node.
 */
  return polnode;
}

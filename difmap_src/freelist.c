#include <stdio.h>
#include <stdlib.h>

#include "freelist.h"
#include "logio.h"

typedef struct FreeListBlock FreeListBlock;
struct FreeListBlock {
  FreeListBlock *next;   /* The next block in the list */
  char *nodes;           /* The array of free-list nodes */
};

struct FreeList {
  size_t node_size;         /* The size of a free-list node */
  unsigned blocking_factor; /* The number of nodes per block */
  long nbusy;               /* The number of nodes that are in use */
  FreeListBlock *block;     /* The head of the list of free-list blocks */
  char *free_list;          /* The free-list of nodes */
};

static FreeListBlock *new_FreeListBlock(const char *caller, FreeList *fl);
static FreeListBlock *del_FreeListBlock(FreeListBlock *fl);
static void thread_FreeListBlock(FreeList *fl, FreeListBlock *block);

/*.......................................................................
 * Allocate a new free-list from blocks of 'blocking_factor' objects of size
 * node_size.
 *
 * Input:
 *  caller      const char *  The name of the calling function, for use in
 *                            error messages.
 *  node_size       size_t    The size of the free-list nodes to be returned
 *                            by new_FreeListNode(). Use sizeof() to
 *                            determine this.
 *  blocking_factor unsigned  The number of objects of size 'object_size'
 *                            to allocate per block.
 */
FreeList *new_FreeList(const char *caller, size_t node_size,
		       unsigned blocking_factor)
{
  FreeList *fl;  /* The new free-list container */
/*
 * When a free-list node is on the free-list, it is used as a (void *)
 * link field. Roundup node_size to a mulitple of the size of a void
 * pointer. This, plus the fact that the array of nodes is obtained via
 * malloc, which returns memory suitably aligned for any object, will
 * ensure that the first sizeof(void *) bytes of each node will be
 * suitably aligned to use as a (void *) link pointer.
 */
  node_size = sizeof(void *) *
    ((node_size + sizeof(void *) - 1) / sizeof(void *));
/*
 * Enfore a minimum block size.
 */
  if(blocking_factor < 1)
    blocking_factor = 1;
/*
 * Allocate the container of the free list.
 */
  fl = (FreeList *) malloc(sizeof(FreeList));
  if(!fl) {
    lprintf(stderr, "new_FreeList (%s): Insufficient memory.\n",
	    caller  ? caller : "unknown caller");
    return NULL;
  };
/*
 * Before attempting any operation that might fail, initialize the
 * container at least up to the point at which it can safely be passed
 * to del_FreeList().
 */
  fl->node_size = node_size;
  fl->blocking_factor = blocking_factor;
  fl->nbusy = 0;
  fl->block = NULL;
  fl->free_list = NULL;
/*
 * Allocate the first block of memory.
 */
  fl->block = new_FreeListBlock(caller, fl);
  if(!fl->block)
    return del_FreeList(caller, fl, 1);
/*
 * Add the new list of nodes to the free-list.
 */
  fl->free_list = fl->block->nodes;
/*
 * Return the free-list for use.
 */
  return fl;
}

/*.......................................................................
 * Re-thread a freelist to reclaim all allocated nodes.
 * This function should not be called unless if it is known that none
 * of the currently allocated nodes are still being used.
 *
 * Input:
 *  fl          FreeList *  The free-list to be reset, or NULL.
 */
void rst_FreeList(FreeList *fl)
{
  if(fl) {
    FreeListBlock *block;
/*
 * Re-thread the nodes of each block into individual free-lists.
 */
    for(block=fl->block; block; block=block->next)
      thread_FreeListBlock(fl, block);
/*
 * Link all of the block freelists into one large freelist.
 */
    fl->free_list = NULL;
    for(block=fl->block; block; block=block->next) {
/*
 * Locate the last node of the current block.
 */
      char *last_node = block->nodes + fl->node_size *
	(fl->blocking_factor - 1);
/*
 * Make the link-field of the last node point to the first
 * node of the current freelist, then make the first node of the
 * new block the start of the freelist. 
 */
      *(void **)last_node = fl->free_list;
      fl->free_list = block->nodes;
    };
/*
 * All allocated nodes have now been returned to the freelist.
 */
    fl->nbusy = 0;
  };
}

/*.......................................................................
 * Delete a free-list.
 *
 * Input:
 *  caller    const char *  The name of the calling function, for use in
 *                          error messages.
 *  fl          FreeList *  The free-list to be deleted, or NULL.
 *  force            int    If force==0 then del_FreeList() will complain
 *                           and refuse to delete the free-list if any
 *                           of nodes have not been returned to the free-list.
 *                          If force!=0 then del_FreeList() will not check
 *                           whether any nodes are still in use and will
 *                           always delete the list.
 * Output:
 *  return      FreeList *  Always NULL (even if the list couldn't be
 *                          deleted).
 */
FreeList *del_FreeList(const char *caller, FreeList *fl, int force)
{
  if(fl) {
/*
 * Check whether any nodes are in use.
 */
    if(!force && busy_FreeListNodes(fl) != 0) {
      lprintf(stderr, "del_FreeList (%s): %ld nodes are still in use.\n",
	      caller ? caller : "unknown caller", busy_FreeListNodes(fl));
      return NULL;
    };
/*
 * Delete the list blocks.
 */
    {
      FreeListBlock *next = fl->block;
      while(next) {
	FreeListBlock *block = next;
	next = block->next;
	block = del_FreeListBlock(block);
      };
    };
    fl->block = NULL;
    fl->free_list = NULL;
/*
 * Discard the container.
 */
    free(fl);
  };
  return NULL;
}

/*.......................................................................
 * Allocate a new object from a free-list.
 *
 * Input:
 *  caller  const char *  The name of the calling function, for use in
 *                        error messages.
 *  fl        FreeList *  The free-list to return an object from.
 * Output:
 *  return        void *  A new object of the size that was specified via
 *                        the node_size argument of new_FreeList() when
 *                        the free-list was created.
 */
void *new_FreeListNode(const char *caller, FreeList *fl)
{
  void *node;  /* The node to be returned */
/*
 * Check arguments.
 */
  if(!fl) {
    lprintf(stderr, "new_FreeListNode (%s): NULL free-list.\n",
	    caller ? caller : "unknown caller");
    return NULL;
  };
/*
 * If the free-list has been exhausted extend it by allocating
 * another block of nodes.
 */
  if(!fl->free_list) {
    FreeListBlock *block = new_FreeListBlock(caller, fl);
    if(!block)
      return NULL;
/*
 * Prepend the new block to the list of free-list blocks.
 */
    block->next = fl->block;
    fl->block = block;
/*
 * Add the new list of nodes to the free-list.
 */
    fl->free_list = fl->block->nodes;
  };
/*
 * Remove and return a node from the front of the free list.
 */
  node = fl->free_list;
  fl->free_list = *(void **)node;
/*
 * Record the loss of a node from the free-list.
 */
  fl->nbusy++;
/*
 * Return the node.
 */
  return node;
}

/*.......................................................................
 * Return an object to the free-list that it was allocated from.
 *
 * Input:
 *  caller  const char *  The name of the calling function, for use in
 *                        error messages.
 *  fl        FreeList *  The free-list from which the object was taken.
 *  object        void *  The node to be returned.
 * Output:
 *  return        void *  Always NULL.
 */
void *del_FreeListNode(const char *caller, FreeList *fl, void *object)
{
/*
 * Check arguments.
 */
  if(!fl) {
    lprintf(stderr, "del_FreeListNode (%s): NULL free-list.\n",
	    caller ? caller : "unknown caller");
    return NULL;
  };
/*
 * Return the node to the head of the free list.
 */
  if(object) {
    *(void **)object = fl->free_list;
    fl->free_list = object;
/*
 * Record the return of the node to the free-list.
 */
    fl->nbusy--;
  };
  return NULL;
}

/*.......................................................................
 * Return a count of the number of nodes that are currently allocated.
 *
 * Input:
 *  fl      FreeList *  The list to count wrt, or NULL.
 * Output:
 *  return      long    The number of nodes (or 0 if fl==NULL).
 */
long busy_FreeListNodes(FreeList *fl)
{
  return fl ? fl->nbusy : 0;
}

/*.......................................................................
 * Allocate a new list of free-list nodes. On return the nodes will
 * be linked together as a list starting with the node at the lowest
 * address and ending with a NULL next pointer.
 *
 * Input:
 *  caller    const char *  The name of the external calling function,
 *                          for use in error messages.
 *  fl          FreeList *  The free-list to allocate the list for.
 * Output:
 *  return FreeListBlock *  The new linked block of free-list nodes,
 *                          or NULL on error.
 */
static FreeListBlock *new_FreeListBlock(const char *caller, FreeList *fl)
{
  FreeListBlock *block;  /* The new block to be returned */
/*
 * Allocate the container.
 */
  block = (FreeListBlock *) malloc(sizeof(FreeListBlock));
  if(!block) {
    lprintf(stderr, "new_FreeListBlock (%s): Insufficient memory.\n",
	    caller ? caller : "unknown caller");
    return NULL;
  };
/*
 * Before attempting any operation that might fail, initialize the
 * container at least up to the point at which it can safely be passed
 * to del_FreeListBlock().
 */
  block->next = NULL;
  block->nodes = NULL;
/*
 * Allocate the block of nodes.
 */
  block->nodes = (char *) malloc(fl->node_size * fl->blocking_factor);
  if(!block->nodes) {
    lprintf(stderr, "new_FreeListBlock (%s): Insufficient memory.\n",
	    caller ? caller : "unknown caller");
    return del_FreeListBlock(block);
  };
/*
 * Initialize the block as a linked list of FreeListNode's.
 */
  thread_FreeListBlock(fl, block);
  return block;
}

/*.......................................................................
 * Link each node of a freelist block to the node that follows it.
 *
 * Input:
 *  fl         FreeList *   The freelist that contains the block.
 *  block FreeListBlock *   The block to be threaded.
 */
static void thread_FreeListBlock(FreeList *fl, FreeListBlock *block)
{
  char *mem = block->nodes;
  int i;
  for(i=0; i<fl->blocking_factor - 1; i++, mem += fl->node_size)
    *(void **)mem = mem + fl->node_size;  /* Link to the next node */
  *(void **)mem = NULL;                   /* Terminate the list */
}

/*.......................................................................
 * Delete a free-list block.
 *
 * Input:
 *  fl      FreeListBlock *  The block to be deleted, or NULL.
 * Output:
 *  return  FreeListBlock *  Always NULL.
 */
static FreeListBlock *del_FreeListBlock(FreeListBlock *fl)
{
  if(fl) {
    fl->next = NULL;
    if(fl->nodes)
      free(fl->nodes);
    fl->nodes = NULL;
    free(fl);
  };
  return NULL;
}

/*.......................................................................
 * Return non-zero if the specified free-list was created for the
 * specified node size.
 *
 * Input:
 *  fl      FreeList *   The free-list to test.
 *  node_size size_t     The node size to compare with that of the freelist.
 *                       Use sizeof() to determine the size that you need.
 */
int compatible_FreeList(FreeList *fl, size_t node_size)
{
/*
 * new_FreeList() modifies the size to make sure that it is at least
 * as big as a (void *) pointer. Do the same here before comparing
 * sizes.
 */
  if(node_size < sizeof(void *))
    node_size = sizeof(void *);
  return fl && fl->node_size == node_size;
}


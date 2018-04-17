#ifndef freelist_h
#define freelist_h

/*
 * This module provides a memory allocation scheme that helps to
 * prevent memory fragmentation by allocating large blocks of
 * fixed sized objects and forming them into a free-list for
 * subsequent allocations. The free-list is expanded as needed.
 */
typedef struct FreeList FreeList;

/*
 * Allocate a new free-list from blocks of 'blocking_factor' objects of size
 * node_size. The node_size argument should be determined by applying
 * the sizeof() operator to the object type that you intend to allocate from
 * the freelist.
 */
FreeList *new_FreeList(const char *caller, size_t node_size,
		       unsigned blocking_factor);

/*
 * If it is known that none of the nodes currently allocated from
 * a freelist are still in use, the following function can be called
 * to return all nodes to the freelist without the overhead of
 * having to call del_FreeListNode() for every allocated node. The
 * nodes of the freelist can then be reused by future callers to
 * new_FreeListNode().
 */
void rst_FreeList(FreeList *fl);

/*
 * Delete a free-list.
 */
FreeList *del_FreeList(const char *caller, FreeList *fl, int force);

/*
 * Determine the number of nodes that are currently allocated.
 */
long busy_FreeListNodes(FreeList *fl);

/*
 * Allocate a new object from a free-list.
 */
void *new_FreeListNode(const char *caller, FreeList *fl);

/*
 * Return an object to the free-list that it was allocated from.
 */
void *del_FreeListNode(const char *caller, FreeList *fl, void *object);

/*
 * Return non-zero if the specified freelist was created for the specified
 * node size. The node_size argument should be determined by applying
 * the sizeof() operator to the object type that you intend to allocate from
 * the freelist.
 */
int compatible_FreeList(FreeList *fl, size_t node_size);

#endif

#include <stdio.h>
#include <stdlib.h>

#include "logio.h"
#include "vlbconst.h"
#include "intlist.h"
#include "utbin.h"

static Intbin *add_Intbin(Intlist *ilist, double ut, long group, int isub);

enum {BINBLK=256};      /* The number of integration bins allocated at once */
typedef struct Intblk { /* A member of a list of blocks of Intbin's */
  Intbin iblk[BINBLK];  /* A block of integration bins for the free-list */
  struct Intblk *next;  /* Pointer to next block of integration bins */
} Intblk;

typedef struct Subbin { /* Sub-array integration list container */
  double beg_ut;        /* Start UT of the latest integration bin */
  double end_ut;        /* End UT of the latest integration bin */
  Intbin *head;         /* Head of the list of integration bins */
  Intbin *tail;         /* Pointer to the last integration bin appended */
  int ntime;            /* The number of bins in the sub-array */
  struct Subbin *next;  /* Pointer to a container with later integrations */
} Subbin;

struct Intlist {
  double origin;  /* The origin of the regular grid of bins (seconds) */
  double binwid;  /* The width of each bin (seconds) */
  int nsub;       /* Number of sub-arrays */
  Subbin *sbin;   /* Array of 'nsub' sub-array containers */
  Subbin *head;   /* Head of time-sorted list of sub-array containers */
  Intblk *ilist;  /* Head of list of Intbin blocks */
  Intbin *ifree;  /* Head of free list of Intbin's taken from ilist */
};

static void add_Subbin(Intlist *ilist, Subbin *sbin);
static int bad_Intlist(Intlist *ilist, char *client);

/*.......................................................................
 * Construct an integration-bin list container.
 *
 * The returned descriptor should be primed by calling add_group() for
 * each record that is to be iterated through, before use as an iterator.
 *
 * Input:
 *  nsub       int     The max expected number of sub-arrays.
 *  origin  double     The time origin of the regular bin grid (seconds).
 *  binwid  double     The width of each bin (seconds). The minimum bin
 *                     size allowed is 1 second. If binwid is less than
 *                     this then the bins will have zero width and every
 *                     different time-stamp will result in a new bin.
 * Output:
 *  ilist  Intlist *   The integration list container, or NULL on error.
 */
Intlist *new_Intlist(int nsub, double origin, double binwid)
{
  Intlist *ilist;  /* The container to be returned */
  int isub;        /* The index of a sub-array */
/*
 * Allocate the iterator container.
 */
  ilist = (Intlist *) malloc(sizeof(Intlist));
  if(ilist==NULL) {
    lprintf(stderr, "Insufficient memory for integration bin iterator.\n");
    return NULL;
  };
/*
 * Initialize the iterator at least to the point at which it will be
 * safe to call del_Intlist(ilist).
 */
  ilist->origin = origin;
  ilist->binwid = binwid < 1.0 ? 0.0 : binwid;
  ilist->nsub = nsub;
  ilist->sbin = NULL;
  ilist->head = NULL;
  ilist->ilist = NULL;
  ilist->ifree = NULL;
/*
 * Allocate the array of sub-array integration containers.
 */
  ilist->sbin = (Subbin *) malloc(sizeof(Subbin) * ilist->nsub);
  if(ilist->sbin==NULL) {
    lprintf(stderr, "Insufficient memory for integration bin iterator.\n");
    return del_Intlist(ilist);
  };
/*
 * Initialize the new sub-array containers.
 */
  for(isub=0; isub<ilist->nsub; isub++) {
    Subbin *sbin = &ilist->sbin[isub];
    sbin->beg_ut = 0.0;
    sbin->end_ut = 0.0;
    sbin->ntime = 0;
    sbin->head = NULL;
    sbin->tail = NULL;
    sbin->next = NULL;
  };
/*
 * Return the un-primed iterator.
 */
  return ilist;
}

/*.......................................................................
 * Delete an integration-bin list container and its contents.
 *
 * Input:
 *  ilist    Intlist *  The container to be deleted.
 * Output:
 *  return   Intlist *  The deleted container (ie. NULL). Use like:
 *                      ilist = del_Intlist(ilist);
 */
Intlist *del_Intlist(Intlist *ilist)
{
  if(ilist) {
/*
 * Free the array of sub-array containers.
 */
    if(ilist->sbin)
      free(ilist->sbin);
/*
 * Free the list of integration bin blocks.
 */
    if(ilist->ilist) {
      Intblk *next = ilist->ilist;
      while(next) {
	Intblk *prev = next;
	next = next->next;
	free(prev);
      };
    };
/*
 * Delete the container.
 */
    free(ilist);
  };
  return NULL;
}

/*.......................................................................
 * Return the descriptor of the next un-processed integration bin from
 * the list contained in an integration grid iterator.
 *
 * Input:
 *  ilist  Intlist *  The integration list iterator.
 * Output:
 *  return  Intbin *  The next integration bin, or NULL if there are
 *                    no further bins.
 */
Intbin *nxt_Intbin(Intlist *ilist)
{
  Intbin *ibin; /* The pointer to the integration bin to be returned */
  Subbin *sbin; /* The sub-array container with the earliest integration bin */
/*
 * Sanity check the iterator.
 */
  if(bad_Intlist(ilist, "nxt_Intbin"))
    return NULL;
/*
 * If the head of the sub-array integration-bin container list is NULL,
 * attempt to initialize it.
 */
  if(ilist->head==NULL) {
    int isub;
    sbin = ilist->sbin;
    for(isub=0; isub<ilist->nsub; isub++,sbin++)
      add_Subbin(ilist, sbin);
  };
/*
 * Get the sub-array container that is at the head of the list in ilist.
 * At the head of its list of integrations is the earliest integration bin
 * available from any sub-array.
 */
  sbin = ilist->head;
/*
 * If the list is empty then there are no remaining sub-array
 * integration containers.
 */
  if(sbin==NULL)
    return NULL;
/*
 * Remove the sub-array integration container from the sorted list.
 */
  ilist->head = sbin->next;
/*
 * Remove the head of the list of integration bins in the sub-array.
 * This is the earliest available integration bin and is the one to be
 * returned.
 */
  ibin = sbin->head;
  sbin->head = ibin->next;
/*
 * If there are any remaining integration bins in the sub-array, then
 * re-insert the sub-array container in the list of sub-array integration
 * containers, such that the list is maintained in order of the time-stamp
 * of the earliest remaining integration of each sub-array.
 */
  add_Subbin(ilist, sbin);
/*
 * Return the extracted integration bin for use.
 */
  return ibin;
}

/*.......................................................................
 * Insert a sub-array integration list container into the position in
 * the ilist->head list at which the time stamp of its first integration
 * is greater than that of the sub-array container that precedes it
 * and less than that of the container that follows it in the list.
 *
 * Input:
 *  ilist  Intlist *  The integration list iterator.
 *  sbin    Subbin *  The sub-array integration bin list container.
 */
static void add_Subbin(Intlist *ilist, Subbin *sbin)
{
/*
 * Only insert the sub-array container if its list of integration bins
 * has not been exhausted.
 */
  if(sbin->head) {
    double newut = sbin->head->ut;  /* The UT of the sub-array to be inserted */
    Subbin *prev = NULL;            /* The node preceding the insertion point */
    Subbin *next = ilist->head;     /* The node following the insertion point */
/*
 * Find the position in the list where prev and next straddle 'newut'.
 */
    for( ; next && newut > next->head->ut; next=next->next)
      prev = next;
/*
 * Re-insert the sub-array container between nodes 'prev' and 'next'.
 */
    sbin->next = next;
    if(prev)
      prev->next = sbin;
    else
      ilist->head = sbin;
  };
  return;
}

/*.......................................................................
 * Return the index of the next un-processed random-group of a FITS
 * file, from within a given integration bin.
 *
 * Input:
 *  ibin   Intbin *  The integration bin being processed.
 * Output:
 *  return   long    The index of the next group in the integration bin,
 *                   or -1, when all groups have been seen.
 */
long nxt_group(Intbin *ibin)
{
  long igroup = ibin->first <= ibin->last ? ibin->first : -1;
  ibin->first++;
  return igroup;
}

/*.......................................................................
 * Allocate a new integration bin descriptor, append it to the list
 * of integration bins maintained in the appropriate sub-array container,
 * and return the bin.
 *
 * Input:
 *  ilist   Intlist *  The integration bin list container.
 *  ut       double    The time-stamp to give the integration bin.
 *  group      long    The index of the first group to include in the
 *                     integration bin.
 *  isub        int    The index of the sub-array to which the integration
 *                     belongs.
 * Output:
 *  return   Intbin *  The new integration bin, or NULL on error.
 */
static Intbin *add_Intbin(Intlist *ilist, double ut, long group, int isub)
{
  Intbin *ibin;  /* The pointer to the integration bin to be returned */
  Subbin *sbin;  /* The container of sub-array 'isub' */
  int i;
/*
 * If the integration bin free-list is empty, allocate a new block.
 */
  if(ilist->ifree==NULL) {
/*
 * Allocate a block of NBIN integration bins.
 */
    Intblk *iblk = (Intblk *) malloc(sizeof(Intblk));
    if(iblk==NULL) {
      lprintf(stderr, "add_Intbin: Insufficient memory to bin data.\n");
      return NULL;
    };
/*
 * Add the new block to the head of the list in ilist.
 */
    iblk->next = ilist->ilist;
    ilist->ilist = iblk;
/*
 * Link all the integration bins together to form the new free-list.
 */
    ibin = ilist->ifree = &iblk->iblk[0];
    for(i=0; i<BINBLK-1; i++,ibin++)
      ibin->next = ibin + 1;
    ibin->next = NULL;
  };
/*
 * Remove an integration bin from the head of the free-list.
 */
  ibin = ilist->ifree;
  ilist->ifree = ibin->next;
/*
 * Initialize the new bin.
 */
  ibin->ut = ut;
  ibin->first = ibin->last = group;
  ibin->isub = isub;
/*
 * Append it to the list in the container of the given sub-array.
 */
  sbin = &ilist->sbin[isub];
  if(sbin->head)
    sbin->tail->next = ibin;
  else 
    sbin->head = ibin;
  sbin->tail = ibin;
  ibin->next = NULL;
/*
 * Increment the count of the number of integrations in the sub-array in
 * which the integration bin was added.
 */
  sbin->ntime++;
/*
 * Return the new bin.
 */
  return ibin;
}

/*.......................................................................
 * Append the record of new group to in the list of integration bins
 * for the given sub-array.
 *
 * Input:
 *  ilist   Intlist *  The integration list container.
 *  ut       double    The time-stamp of the group (seconds).
 *  group      long    The index of the group in the FITS file.
 *  isub        int    The index of the sub-array to which the group belongs.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int add_group(Intlist *ilist, double ut, long group, int isub)
{
  Subbin *sbin;  /* The integration bin container of sub-array 'isub' */
/*
 * Check the descriptor.
 */
  if(bad_Intlist(ilist, "add_group"))
    return 1;
/*
 * Get the appropriate sub-array integration container for the new group.
 */
  sbin = &ilist->sbin[isub];
/*
 * Does the group time-stamp lie within the integration bin last added
 * to the current sub-array?
 */
  if(sbin->tail && ut <= sbin->end_ut) {
/*
 * Ensure that the group is not out of time order.
 */
    if(ut < sbin->beg_ut) {
      lprintf(stderr,"add_group: The visibilities must be in time order.\n");
      return 1;
    };
/*
 * Expand the recorded group index range of the integration bin to encompass the
 * new group.
 */
    sbin->tail->last = group;
/*
 * Start a new integration bin for the new group in the associated sub-array.
 */
  } else {
/*
 * Get the time limits of the integration bin into which 'ut' falls.
 */
    UTbin *utbin = bintime(ilist->origin, ut, ilist->binwid);
/*
 * Add the new integration bin to the list in the current sub-array.
 */
    if(add_Intbin(ilist, utbin->mid_ut, group, isub)==NULL)
      return 1;
/*
 * Reset the sub-array container for the new integration bin.
 */
    sbin->beg_ut = utbin->beg_ut;
    sbin->end_ut = utbin->end_ut;
  };
  return 0;
}

/*.......................................................................
 * Return the number of integrations bins in a given sub-array.
 *
 * Input:
 *  ilist Intlist *  The integration list iterator.
 *  isub      int    The sub-array whose statistics are required.
 * Output:
 *  return    int    The number of integration bins in sub-array 'isub'.
 */
int ibin_count(Intlist *ilist, int isub)
{
/*
 * Sanity checks.
 */
  if(bad_Intlist(ilist, "ibin_count"))
    return 0;
  if(isub<0 || isub>=ilist->nsub) {
    lprintf(stderr, "Bad sub-array index.\n");
    return 0;
  };
/*
 * Return the number of integration bins recorded for the given sub-array.
 */
  return ilist->sbin[isub].ntime;
}

/*.......................................................................
 * Report error and return 1 if the given Intlist is invalid.
 *
 * Input:
 *  ilist  Intlist *   The descriptor to be checked.
 *  client    char *   The name of the calling function.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int bad_Intlist(Intlist *ilist, char *client)
{
  if(ilist==NULL) {
    lprintf(stderr, "%s: NULL Intlist descriptor intercepted.\n", client);
    return 1;
  };
  return 0;
}

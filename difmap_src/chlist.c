#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "logio.h"
#include "chlist.h"

/*
 * Note that channel ranges are recorded in a contiguous array of containers,
 * and the array is re-sized as ranges are added. An array is used
 * rather than a linked list for efficiency reasons. The frequency
 * channel ranges are normally accessed in inner loops, where
 * efficiency is a major concern. Array access is the easiest for
 * optimizers to optimize, and the use of a contiguous array avoids
 * paging. As ranges are added the array is re-sized in increments of
 * RANGE_INC range containers. Thus for a given value of cl->nrange>0,
 * there are actually   RANGE_INC + RANGE_INC * [(cl->nrange-1) /
 * RANGE_INC]; range elements allocated. The most frequently used
 * number of ranges is likely to be 1, whereas I doubt that people
 * will ever want more than 20 ranges [at least not very often].
 */
#define RANGE_INC 5

/*.......................................................................
 * Allocate a new [empty] channel range container.
 *
 * Output:
 *  return   Chlist *  The allocated and initialized container.
 */
Chlist *new_Chlist(void)
{
  Chlist *cl;   /* The new container */
/*
 * Allocate the container.
 */
  cl = (Chlist *) malloc(sizeof(Chlist));
  if(cl==NULL) {
    lprintf(stderr, "Insufficient memory for a new channel range list.\n");
    return cl;
  };
/*
 * Initialize it.
 */
  cl->range = NULL;
  cl->nrange = 0;
  cl->ca = cl->cb = 0;
/*
 * Return the empty container.
 */
  return cl;
}

/*.......................................................................
 * Delete a channel range list container and its contents.
 *
 * Input:
 *  cl      Chlist *   The container to be deleted.
 * Output:
 *  return  Chlist *   The deleted container (ie. NULL).
 */
Chlist *del_Chlist(Chlist *cl)
{
  if(cl) {
    if(cl->range)
      free(cl->range);
    free(cl);
  };
  return NULL;
}

/*.......................................................................
 * Add a range of channels to a channel list container.
 * If the channel range overlaps an existing one, then the two ranges
 * will be merged. New ranges are placed so as to keep the ranges sorted
 * in increasing channel order.
 *
 * Input:
 *  cl      Chlist *  The channel list container.
 *  ca         int    The (0-relative) index of one end of the new channel
 *                    range.
 *  cb         int    The (0-relative) index of the opposite end of the new
 *                    channel range to ca.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error. 
 */
int add_crange(Chlist *cl, int ca, int cb)
{
  int irange;   /* Index into cl->range at which to install the new range */
  int i;
/*
 * Sort ca and cb into ascending order.
 */
  if(ca > cb) {int new_ca=cb; cb=ca; ca=new_ca;};
/*
 * Check sanity.
 */
  if(ca < 0) {
    lprintf(stderr, "add_crange: Illegal channel index: %d.\n", ca+1);
    return 1;
  };
  if(cl==NULL) {
    lprintf(stderr, "add_crange: NULL channel-list container.\n");
    return 1;
  };
/*
 * The channel ranges are stored in ascending order. Search for the position
 * at which to insert the new range, or add to an existing range.
 * There are unlikely to be enough ranges to warrant a binary search.
 */
  for(irange=0; irange<cl->nrange && ca > cl->range[irange].cb+1; irange++);
/*
 * Extend an existing range?
 */
  if(irange<cl->nrange &&
    (ca >= cl->range[irange].ca-1 || cb >= cl->range[irange].ca-1)) {
    Chans *range = &cl->range[irange];
/*
 * Extend to lower channels?
 */
    if(ca < range->ca) {
      range->ca = ca;
    };
/*
 * Extend to higher channels?
 */
    if(cb > range->cb) {
      int ir;
/*
 * Search for the last range that is overlapped by the extension.
 */
      for(ir=cl->nrange-1; ir>irange && cb < cl->range[ir].ca-1; ir--);
/*
 * Get the upper limit of the encompassing range.
 */
      range->cb = cb > cl->range[ir].cb ? cb : cl->range[ir].cb;
/*
 * Remove superfluous ranges.
 */
      if(ir > irange) {
	int ioff = ir - irange;
	for(ir = irange+1; ir < cl->nrange; ir++)
	  cl->range[ir] = cl->range[ir+ioff];
	cl->nrange -= ioff;
      };
    };
/*
 * Add a new range at cl->range[irange].
 */
  } else {
/*
 * Extend the range list?
 */
    if(cl->nrange % RANGE_INC == 0) {
      Chans *range;
      if(cl->range) {
	range = (Chans *) realloc(cl->range, sizeof(Chans) *
				  (cl->nrange + RANGE_INC));
      } else {
	range = (Chans *) malloc(sizeof(Chans) * RANGE_INC);
      };
/*
 * Assign the re-sized array of channel ranges.
 */
      if(range) {
	cl->range = range;
      } else {
	lprintf(stderr, "Insufficient memory to store new channel range.\n");
	return 1;
      };
    };
/*
 * Make room for the new range.
 */
    for(i=irange; i<cl->nrange; i++) {
      cl->range[i+1] = cl->range[i];
    };
/*
 * Record the new channel range in the vacated slot.
 */
    cl->range[irange].ca = ca;
    cl->range[irange].cb = cb;
/*
 * Increment the count of ranges.
 */
    cl->nrange++;
  };
/*
 * Update the recorded min and max channel indexes.
 */
  cl->ca = cl->range[0].ca;
  cl->cb = cl->range[cl->nrange-1].cb;
  return 0;
}

/*.......................................................................
 * Truncate the given channel ranges to only refer to channels up to
 * channel nchan-1. Note that this function may determine that *all* of
 * the chosen channel ranges be dropped. In this case lim_Chlist will
 * return 0, and it is up to the caller to handle this case appropriately.
 *
 * Input:
 *  cl     Chlist *    The list to be truncated.
 * Output:
 *  return    int      The number of remaining channel ranges, or -1
 *                     on error.
 */
int lim_Chlist(Chlist *cl, int nchan)
{
  int ir;  /* The channel range being checked */
/*
 * Sanity check.
 */
  if(cl==NULL) {
    lprintf(stderr, "lim_Chlist: NULL channel range list container.\n");
    return -1;
  };
/*
 * Search for the first range whose upper channel bound exceeds nchan-1.
 */
  if(cl->nrange > 0) {
    for(ir=0; ir<cl->nrange && cl->range[ir].cb < nchan; ir++);
/*
 * The channel ranges are in ascending channel order, so if a bad
 * range was located truncate the range there.
 */
    if(ir < cl->nrange) {
      lprintf(stderr, "Restricting channel ranges to the available %d channels.\n", nchan);
/*
 * Remove the entire range found?
 */
      if(cl->range[ir].ca >= nchan) {
	cl->nrange = ir;
/*
 * Is just the upper bound of the range bad?
 */
      } else {
	cl->range[ir].cb = nchan - 1;
	cl->nrange = ir+1;
      };
    };
  };
/*
 * Record the encompassing channel range.
 */
  if(cl->nrange >= 1) {
    cl->ca = cl->range[0].ca;
    cl->cb = cl->range[cl->nrange-1].cb;
  } else {
    lprintf(stderr, "None of the chosen ranges of channels exist.\n");
    cl->ca = cl->cb = 0;
  };
  return cl->nrange;
}

/*.......................................................................
 * Contruct a new channel list from a sub-set of an existing list.
 *
 * Input:
 *  cl     Chlist *   The list to copy ranges from.
 *  coff      int     The offset of channel 0 in the new list wrt the
 *                    channels in 'cl'. coff can be negative.
 *  nchan     int     The potential range of channels in the new channel
 *                    list.
 * Output:
 *  return Chlist *   The new list, or NULL on error.
 */
Chlist *sub_Chlist(Chlist *cl, int coff, int nchan)
{
  Chlist *scl;  /* The new channel list */
  Chans *range; /* A single channel range from the input list */
/*
 * Check arguments.
 */
  if(!cl) {
    lprintf(stderr, "sub_Chlist: NULL channel list received.\n");
    return NULL;
  };
  if(nchan < 0) {
    lprintf(stderr, "sub_Chlist: nchan < 0\n");
    return NULL;
  };
/*
 * Create the new list.
 */
  if(!(scl = new_Chlist()))
    return NULL;
/*
 * Locate channel ranges that fall within the domain of the output list.
 */
  for(range=cl->range; range<cl->range + cl->nrange; range++) {
/*
 * Convert the input channel range into output list channel numbers.
 */
    int ca = range->ca - coff;
    int cb = range->cb - coff;
/*
 * If the input range overlaps the domain of the output list,
 * truncate the list at each end to keep it within 0..nchan-1 and
 * add the new range to the output list.
 */
    if(ca < nchan && cb >= 0) {
      if(ca < 0) ca = 0;
      if(cb >= nchan) cb = nchan - 1;
      if(add_crange(scl, ca, cb))
	return del_Chlist(scl);
    };
  };
  return scl;
}

/*.......................................................................
 * Allocate a new copy of an existing channel list.
 *
 * Input:
 *  cl      Chlist *  The channel list to be copied.
 * Output:
 *  return  Chlist *  The newly allocated copy of *cl.
 */
Chlist *cpy_Chlist(Chlist *cl)
{
  Chlist *newcl;   /* The new copy of *cl */
  int size;        /* The size of the cl->range[] array */
  int i;
/*
 * Copy a NULL list?
 */
  if(!cl)
    return NULL;
/*
 * Create the new channel list container.
 */
  newcl = new_Chlist();
  if(!newcl)
    return NULL;
/*
 * Copy non-pointer members.
 */
  *newcl = *cl;
  newcl->range = NULL;
/*
 * Compute the allocated size of the cl->range[] array. This will
 * be cl->nrange rounded up to the next multiple of RANGE_INC.
 */
  size = RANGE_INC * ((cl->nrange + RANGE_INC - 1) / RANGE_INC);
/*
 * Allocate the range array.
 */
  if(size > 0) {
    newcl->range = malloc(sizeof(Chans) * size);
    if(!newcl->range)
      return del_Chlist(newcl);
/*
 * Copy the ranges.
 */
    for(i=0; i<cl->nrange; i++)
      newcl->range[i] = cl->range[i];
  };
  return newcl;
}

/*.......................................................................
 * Return 1 if two channel lists match.
 *
 * Input:
 *  cl1    Chlist *  The first of the channel lists to compare.
 *  cl2    Chlist *  The second of the channel lists to compare.
 * Output:
 *  return    int    0 - They don't match.
 *                   1 - They do match.
 */
int eq_Chlist(Chlist *cl1, Chlist *cl2)
{
  int range;  /* The index of the range being compared */
/*
 * If the lists don't have the same number of channel ranges, then
 * they clearly don't match.
 */
  if(cl1->nrange != cl2->nrange)
    return 0;
/*
 * Compare the channel ranges.
 */
  for(range=0; range < cl1->nrange; range++) {
    Chans *c1 = cl1->range + range;
    Chans *c2 = cl2->range + range;
    if(c1->ca != c2->ca || c1->cb != c2->cb)
      return 0;
  };
/*
 * They match.
 */
  return 1;
}

/*.......................................................................
 * Write a list of channels to a text stream.
 *
 * Input:
 *  cl     Chlist *  The list of channel ranges.
 *  fp       FILE *  The text output stream to write to.
 *  filename char *  The name of the filename being written to, or NULL
 *                   to not report I/O errors.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error (unreported if filename==NULL).
 */
int write_Chlist(Chlist *cl, FILE *fp, char *filename)
{
  int i;
/*
 * Check the arguments.
 */
  if(!cl || !fp) {
    lprintf(stderr, "write_Chlist: NULL argument(s).\n");
  };
/*
 * Write the channel ranges, separated by commas.
 */
  for(i=0; i<cl->nrange; i++) {
    if(lprintf(fp, "%s%d, %d", i ? ", " : "", cl->range[i].ca+1,
	       cl->range[i].cb+1) < 0) {
      if(filename)
	lprintf(stderr, "Error writing channel ranges to file: %s\n", filename);
      return 1;
    };
  };
  return 0;
}

/*.......................................................................
 * Read a range of channels from a text input stream, as previously
 * written by write_Chlist().
 *
 * Input:
 *  fp       FILE *  The stream to read from.
 *  filename char *  The filename of the input stream, or NULL if
 *                   no error reporting is to be performed.
 *  nline     int    The current line number within the input file,
 *                   for use in error reporting when filename!=NULL.
 *                   Note that this function does not cross line
 *                   boundaries.
 * Output:
 *  cl     Chlist *  The range of channels, or NULL on error.
 */
Chlist *read_Chlist(FILE *fp, char *filename, int nline)
{
  Chlist *cl;   /* The newly allocated channel list */
  int chan[2];  /* Up to two channel numbers of a range */
  int ic=0;     /* The number of channel numbers in chan[] */
  int c = '\0'; /* The latest character read */       
  int first=1;  /* True until one channel number has been read */
/*
 * Check the arguments.
 */
  if(!fp) {
    lprintf(stderr, "read_Chlist: NULL argument(s).\n");
    return NULL;
  };
/*
 * Allocate the channel list.
 */
  cl = new_Chlist();
  if(!cl)
    return NULL;
/*
 * Read the channel ranges and add them to the channel-range list.
 */
  while(c==',' || first) {
/*
 * Skip leading spaces.
 */
    do {c = getc(fp);} while(c==' ' || c=='\t');
/*
 * Read the channel number.
 */
    if(!isdigit(c) || ungetc(c, fp)==EOF || fscanf(fp, "%d", &chan[ic++])!=1) {
      if(filename)
	lprintf(stderr, "Bad channel number on line %d or %s.\n",
		nline, filename);
      return del_Chlist(cl);
    };
/*
 * If a pair of channels has been read, add it to the channel list, and
 * prepare for the next pair.
 */
    if(ic==2) {
      if(add_crange(cl, chan[0]-1, chan[1]-1))
	return del_Chlist(cl);
      ic = 0;
    };
/*
 * Skip trailing spaces.
 */
    do {c = getc(fp);} while(c==' ' || c=='\t');
/*
 * At least one channel has been read now.
 */
    first = 0;
  };
/*
 * Attempt to put-back the character that we just read, so that
 * it can be interpretted by the caller.
 */
  if(ungetc(c, fp) == EOF) {
    lprintf(stderr, "read_Chlist: ungetc() error.\n");
    return del_Chlist(cl);
  };
/*
 * The last channel range is allowed to be a single number, which
 * is interpretted as a single-channel channel-range.
 */
  if(ic==1 && add_crange(cl, chan[0]-1, chan[0]-1))
    return del_Chlist(cl);
/*
 * Return the completed list.
 */
  return cl;
}

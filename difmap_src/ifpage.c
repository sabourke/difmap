#include <stdlib.h>
#include <stdio.h>

#include "logio.h"
#include "recio.h"
#include "dpage.h"
#include "ifpage.h"

/*
 * Define a Dvis structure used to initialize visibilities in the
 * integration buffer.
 */
static Dvis zero_vis = {0.0f,0.0f,0.0f};

static IFpage *ipmemerr(IFpage *ip);

/*.......................................................................
 * Allocate and intialize a IF-paging file descriptor. This is used
 * for buffered reading and writing of up to one integration of a given
 * IF at a time from the IF paging file.
 *
 * Input:
 *  nif       int   The number of IFs per integration.
 *  nbase     int   The number of baselines per integration.
 *  ntime     int   The number of integrations per IF.
 * Output:
 *  return IFpage * The allocated and initialized descriptor, or NULL
 *                  on error.
 */
IFpage *new_IFpage(int nif, int nbase, int ntime)
{
  IFpage *ip;   /* Pointer to the new descriptor */
/*
 * Attempt to allocate the new descriptor.
 */
  ip = malloc(sizeof(*ip));
  if(ip==NULL)
    return ipmemerr(ip);
/*
 * Intialize ip at least up to the point at which it can safely be passed to
 * del_IFpage().
 */
  ip->rio = NULL;
  ip->dvis = NULL;
  ip->first = 0;
  ip->nread = 0;
  ip->ioerr = 0;
/*
 * Install input parameters.
 */
  ip->nif = nif;
  ip->nbase = nbase;
  ip->ntime = ntime;
/*
 * Open the paging scratch file. Assign a logical record size equal to
 * one IF of integrations.
 */
  ip->rio = new_Recio("ifdata.scr", IS_SCR, 0,
		      ip->nbase * ip->ntime * sizeof(Dvis));
  if(ip->rio==NULL)
    return del_IFpage(ip);
/*
 * Allocate a visibility buffer sufficient to contain one integration of
 * a single IF.
 */
  ip->dvis = malloc(ip->nbase * sizeof(Dvis));
  if(ip->dvis==NULL)
    return ipmemerr(ip);
/*
 * Initialize the buffer.
 */
  if(ip_clear(ip))
    return del_IFpage(ip);
/*
 * Return the intialized descriptor.
 */
  return ip;
}

/*.......................................................................
 * Private cleanup function of new_IFpage() for memory allocation failures.
 *
 * Input:
 *  ip      IFpage *   The partially initialized IFpage descriptor.
 * Output:
 *  return  IFpage *   Allways NULL.
 */
static IFpage *ipmemerr(IFpage *ip)
{
  lprintf(stderr, "new_IFpage: Insufficient memory.\n");
  return del_IFpage(ip);
}

/*.......................................................................
 * Delete an IFpage descriptor and its contents.
 *
 * Input:
 *  ip      IFpage *   A descriptor originally returned by new_IFpage().
 * Output:
 *  return  IFpage *   Always NULL. Use like ip=del_IFpage(ip);
 */
IFpage *del_IFpage(IFpage *ip)
{
  if(ip) {
/*
 * Close and delete the scratch file.
 */
    if(ip->rio)
      ip->rio = del_Recio(ip->rio);
/*
 * Release memory allocated to the visibility I/O buffer.
 */
    if(ip->dvis)
      free(ip->dvis);
/*
 * Delete the container.
 */
    free(ip);
  };
  return NULL;
}

/*.......................................................................
 * Write the previously read portion of the visibility buffer ip->dvis[]
 * to the scratch file.
 *
 * Input:
 *  ip    IFpage *  The IF paging descriptor.
 *  ut      long    The (0-relative) index of the integration to be
 *                  written.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int ip_write(IFpage *ip, long ut)
{
  if(ip_error(ip, "ip_write"))
    return 1;
/*
 * Check the integration index.
 */
  if(ut < 0 || ut > ip->ntime) {
    lprintf(stderr, "ip_write: Integration index out of range.\n");
    return 1;
  };
/*
 * Position the file if necessary.
 */
  if(rec_seek(ip->rio, ip->cif, (ut * ip->nbase + ip->first) * sizeof(Dvis))) {
    ip->ioerr = 1;
    return 1;
  };
/*
 * Write the buffer to the scratch file.
 */
  if(rec_write(ip->rio, ip->nread, sizeof(Dvis), &ip->dvis[ip->first]) < ip->nread) {
    lprintf(stderr, "ip_write: Error writing to scratch file.\n");
    ip->ioerr = 1;
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Read a given portion of an integration into the corresponding
 * portion of the visibility buffer ip->dvis[].
 * The range of visibilities read will be that set by ip->first and
 * ip->nread.
 *
 * Input:
 *  ip     IFpage *  The IF paging descriptor.
 *  ut      long    The (0-relative) index of the integration to be read.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int ip_read(IFpage *ip, long ut)
{
/*
 * Check validity of ip.
 */
  if(ip_error(ip, "ip_read"))
    return 1;
/*
 * Check integration index.
 */
  if(ut < 0 || ut >= ip->ntime) {
    lprintf(stderr, "ip_read: Integration index out of range.\n");
    return 1;
  };
/*
 * Position the file if necessary.
 */
  if(rec_seek(ip->rio, ip->cif, (ut * ip->nbase + ip->first) * sizeof(Dvis))) {
    ip->ioerr = 1;
    return 1;
  };
/*
 * Read from the scratch file into the buffer.
 */
  if(rec_read(ip->rio, ip->nread, sizeof(Dvis), &ip->dvis[ip->first]) < ip->nread) {
    lprintf(stderr, "ip_read: Error reading from scratch file.\n");
    ip->ioerr = 1;
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Set the range of visibilities to be read and written in subsequent
 * reads and writes.
 *
 * Input:
 *  ip  IFpage *   The paging file descriptor.
 *  ifa    int     The index of the required IF.
 *  ba     int     The index of the first baseline.
 *  bb     int     The index of the last baseline.
 * Output:
 *  return int     0 - OK.
 *                 1 - Error.
 */
int ip_range(IFpage *ip, int ifa, int ba, int bb)
{
/*
 * Quietly ignore range setting if ip is NULL.
 */
  if(ip==NULL)
    return 0;
/*
 * Check the validity of ip.
 */
  if(ip_error(ip, "ip_range"))
    return 1;
/*
 * Ensure that ba <= bb.
 */
  if(ba > bb) {int tmp=ba; ba=bb; bb=tmp;};
/*
 * Sanity check the baseline range.
 */
  if(ba<0 || bb>ip->nbase) {
    lprintf(stderr, "ip_range: Baseline indexes out of range.\n");
    return 1;
  };
/*
 * Check the IF selection.
 */
  if(ifa < 0 || ifa >= ip->nif) {
    lprintf(stderr, "ip_range: IF index out of range.\n");
    return 1;
  };
/*
 * Set the new range.
 */
  ip->cif = ifa;
  ip->first = ba;
  ip->nread = bb - ba + 1;
  return 0;
}

/*.......................................................................
 * Check the validity of the IFpage descriptor.
 *
 * Input:
 *  ip     IFpage *   The IF paging file descriptor.
 *  fname   char *   The name of the calling function.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int ip_error(IFpage *ip, const char *fname)
{
/*
 * No descriptor?
 */
  if(ip==NULL) {
    lprintf(stderr, "%s: Intercepted NULL IFpage descriptor.\n", fname);
    return 1;
  };
/*
 * Previous I/O error?
 */
  if(ip->ioerr)
    return 1;
  return 0;
}

/*.......................................................................
 * Clear the whole output buffer.
 *
 * Input:
 *  ip    IFpage *  The IF paging file descriptor.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int ip_clear(IFpage *ip)
{
  Dvis *dvis;   /* Pointer to array of ip->nbase visibilities */
  long nvis;    /* Number of visibilities left to be cleared */
/*
 * Quietly ignore this setup call if ip is NULL.
 */
  if(ip==NULL)
    return 0;
/*
 * Check the validity of ip.
 */
  if(ip_error(ip, "ip_clear"))
    return 1;
/*
 * Clear the buffer.
 */
  nvis = ip->nbase;
  dvis = ip->dvis;
  while(nvis--)
    *dvis++ = zero_vis;
  return 0;
}

/*.......................................................................
 * Make sure that a IF paging file is up to date by flushing all
 * I/O.
 *
 * Input:
 *  ip   IFpage *  The IF paging file descriptor.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int ip_flush(IFpage *ip)
{
  return ip ? rec_flush(ip->rio) : 0;
}
